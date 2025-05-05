# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright (c) 2015 - 2025  Institut Pasteur, Paris and CNRS.                     #
# See the COPYRIGHT file for details                                               #
#                                                                                  #
# integron_finder is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU General Public License as published by             #
# the Free Software Foundation, either version 3 of the License, or                #
# (at your option) any later version.                                              #
#                                                                                  #
# integron_finder is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                   #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    #
# GNU General Public License for more details.                                     #
#                                                                                  #
# You should have received a copy of the GNU General Public License                #
# along with this program (COPYING file).                                          #
# If not, see <http://www.gnu.org/licenses/>.                                      #
####################################################################################

from abc import ABC, abstractmethod
import os
import subprocess
import shlex
from collections import namedtuple
import re
import importlib.util
from enum import Enum


import colorlog
import pandas as pd
from Bio import SeqIO
from integron_finder import IntegronError

_log = colorlog.getLogger(__name__)


"""Sequence description with fields: id strand start stop"""
SeqDesc = namedtuple('SeqDesc', ('id', 'strand', 'start', 'stop'))


class ProteinDB(ABC):
    """
    AbstractClass defining the interface for ProteinDB.
    ProteinDB provide an abstraction and a way to access to proteins corresponding
    to the replicon/contig CDS.
    """

    def __init__(self, replicon, cfg, prot_file=None):
        """

        :param replicon:
        :param cfg:
        :param prot_file:
        """
        self.cfg = cfg
        self.replicon = replicon
        self._prot_file = self._make_protfile(path=prot_file)
        self._prot_db = self._make_db()
        self._pseudo_genes = set()


    def __del__(self):
        ##### CAUTION #########
        # self._make_db() use biopython SeqIO.index
        # and index open a fileIO but not close it
        self.close()


    def close(self):
        if hasattr(self, '_prot_db'):
            self._prot_db.close()


    @abstractmethod
    def __getitem__(self, prot_seq_id):
        """

        :param str prot_seq_id: the id of a protein sequence
        :return: The Sequence corresponding to the prot_seq_id.
        :rtype: :class:`Bio.SeqRecord` object
        :raise KeyError: when seq_id does not match any sequence in DB
        """
        pass

    @abstractmethod
    def __iter__(self):
        """
        :return: a generator which iterate on the protein seq_id which constitute the contig.
        :rtype: generator
        """
        pass

    @abstractmethod
    def _make_protfile(self, path=None):
        """
        Create fasta file with protein corresponding to the nucleic sequence (replicon)

        :return: the path of the created protein file
        :rtype: str
        """
        pass


    def _make_db(self):
        """
        :return: an index of the sequence contains in protfile corresponding to the replicon
        """
        idx = SeqIO.index(self._prot_file, "fasta")
        return idx


    @abstractmethod
    def get_description(self, gene_id):
        """
        :param str gene_id: a protein/gene identifier
        :return: The description of the protein corresponding to the gene_id
        :rtype: :class:`SeqDesc` namedtuple object
        :raise IntegronError: when gene_id is not a valid Gembase gene identifier
        :raise KeyError: if gene_id is not found in GembaseDB instance
        """
        pass

    @property
    def protfile(self):
        """
        :return: The absolute path to the protein file corresponding to contig id
        :rtype: str
        """
        return self._prot_file

    @abstractmethod
    def coding_prot_ids(self):
        """
        :return: a generator which iterate on coding genes seq_ids (the non coding genes are discarded)
        :rtype: generartor
        """
        pass

    def is_pseudo_gene(self, seq_id):
        """

        :param seq_id: The sequence id to test
        :return: True if the seq_id correspond to gene which is a pseudo gene, False otherwise
        """
        return seq_id in self._pseudo_genes


class GembaseType(Enum):
    """
    Modelize the Gembase type en version
    """

    COMPLETE_2plus = 1
    DRAFT_2plus = 2
    COMPLETE_2 = 3
    DRAFT_2 = 4
    COMPLETE_1 = 5
    DRAFT_1 = 6

    def __str__(self):
        _type, vers = self.name.split('_')

        _type = _type.capitalize()
        return f"vers: {vers} {_type}"

    @property
    def complete(self) -> bool:
        """

        :return: True if the GembaseType is Complete genome or False if it's a Draft
        :rtype: bool
        """
        return self.name.startswith('COMPLETE')

    @property
    def version(self) -> str:
        """
        :return: the gembase version
        :rtype: str
        """
        return self.name.split('_')[-1]


class RepliconType(Enum):
    """
    Modelize Replicon type in Gembase
    """

    CHROMOSOME = 1
    PLASMID = 2
    PHAGE = 3
    OTHER = 4
    DRAFT = 5


    def __str__(self):
        return self.name.capitalize()

    def topology(self):
        """
        :return: the default topology of this replicon type 'circ' | 'lin'
        :rtype:  str
        """
        return {1: 'circ',
                2: 'circ',
                3: 'circ',
                4: 'circ',
                5: 'lin'
        }[self.value]


class GembaseDB(ProteinDB):
    """
    Implements :class:`ProteinDB` from a Gembase.
    Managed proteins from Proteins directory corresponding to a replicon/contig
    """

    _gene_patterns = {
        GembaseType.COMPLETE_2plus: r'(\w+).(\d{4})\.(\d{5})\.(\d{3})(?P<g_type>[CPVO])_(\w+)_(\d{5})',
        GembaseType.DRAFT_2plus: r'(\w+).(\d{4})\.(\d{5})\.(\d{3,4})(?P<g_type>D)_([ib])_(\d{5})',
        GembaseType.COMPLETE_2: r'(\w{7})\.(\d{4})\.(\d{5})\.(\d{3})(?P<g_type>[CPVO])_(\d{5})',
        GembaseType.DRAFT_2: r'(\w{4})\.(\d{4})\.(\d{5})\.(\d{4})(?P<g_type>[ib])_(\d{5})',
        GembaseType.COMPLETE_1: r'(\w{7})\.(\w|\d{4})\.(\d{5})\.(?P<g_type>[CPVO])(\d{3})_(\d{5})',
        GembaseType.DRAFT_1: r'(\w{4})\.(\d{4})\.(\d{5})\.(?P<g_type>[ib])(\d{4})_(\d{5})'
    }

    _rep_patterns = [
        r'(\w+)\.(\d{4})\.(\d{5})\.(\d{3,4})(?P<g_type>[CPVOD])',  # Complete + Draft V2_plus
        r'(\w{7})\.(\d{4})\.(\d{5})\.(\d{3})(?P<g_type>[CPVO])',  # Complete V2
        r'(\w{7})\.(\w|\d{4})\.(\d{5})\.(?P<g_type>[CPVO])(\d{3})',  # Complete V1
        r'(\w{4})\.(\d{4})\.(\d{5})\.(\d{4})',  # Draft V1 et V2
    ]

    def __init__(self, replicon, cfg, gembase_path=None, prot_file=None):
        """

        :param replicon: The replicon used to create ProteinDB (protein files and extra information)
        :type replicon: :class:`Bio.SeqRecord` object with a extra attribute path
        :param cfg: The integron_finder configuration
        :type cfg: :class:`integron_finder.config.Config` object
        :param prot_file: The path to a protein file in fasta format
                          which is the translation of the replicon

        .. warning::
            The replicon is a modified Bio.SeqRecord object.
            The attribute *path* must be injected in the object
            This attribute represent the path to a fasta file representing this replicon
        """
        _log.debug(f"call GembaseDB with gembase_path= {gembase_path}")
        self.cfg = cfg
        if gembase_path is None:
            self._gembase_path = os.path.dirname(os.path.dirname(os.path.realpath(self.cfg.input_seq_path)))
        else:
            self._gembase_path = os.path.realpath(gembase_path)
        self.replicon = replicon
        self._pseudo_genes = set()
        # in GemBase Draft the files ar based on replicon id
        # but one file can contains several contig
        # the sequence id contains the contig number
        # for instance ACBA.0917.00019 vs ACBA.0917.00019.0001
        # the filenames are based on replicon id
        # but in code the replicon id contains the contig number
        # for gembase complete both filename and seq_id contains contig number
        self._lst_dir = self.get_lst_dir(self._gembase_path)
        self._gembase_file_basename = self.find_gembase_file_basename(self._gembase_path, self.cfg.input_seq_path)
        self._lst_path = self._get_lst_path()

        self._gembase_type = self.gembase_sniffer(self._lst_path)
        self._replicon_name = os.path.splitext(os.path.basename(self.cfg.input_seq_path))[0]
        self._info = self._parse_lst(self._lst_path)
        if self._info.empty:
            msg = f"No CDS reported in {self._lst_path} for the replicon {replicon.id} ."
            _log.warning(msg)
        self.replicon_tye = self.get_replicon_type(seq_id=self.replicon.id)
        if prot_file is None:
            self._prot_file = self._make_protfile()
        else:
            self._prot_file = prot_file
        self._prot_db = self._make_db()


    @staticmethod
    def get_lst_dir( gembase_path):
        """
        :return: The path to the Gembase LST directory
        :rtype: str
        """
        found = False
        for lst_dirname in ('LST', 'LSTINF', 'LSTINFO'):
            lst_dir_path = os.path.join(gembase_path, lst_dirname)
            if os.path.exists(lst_dir_path):
                found = True
                break
        if not found:
            raise IntegronError(f"Neither 'LST' nor 'LSTINF' nor 'LSTINFO' directory found in '{gembase_path}' .")

        return lst_dir_path


    def _get_lst_path(self):
        lst_dir_path = self._lst_dir

        lst_path = os.path.join(lst_dir_path, self._gembase_file_basename + '.lst')
        if not os.path.exists(lst_path):
            raise IntegronError(f"Do not find {lst_path} in {lst_dir_path}")
        return lst_path


    def _parse_lst(self, lst_path):
        """
        Parse the LSTINFO file and extract information specific to the replicon
        :return: `class`:pandas.DataFrame` object
        :raise IntegronError: when LST dir is not found
        """
        if self.gembase_type == GembaseType.DRAFT_1:
            prots_info = self.gembase1_draft_parser(lst_path, self.replicon.id)
        elif self.gembase_type == GembaseType.COMPLETE_1:
            prots_info = self.gembase1_complete_parser(lst_path, self.replicon.id)
        elif self.gembase_type.version == '2':
            prots_info = self.gembase2_parser(lst_path, self.replicon.id)
        elif self.gembase_type.version == '2plus':
            prots_info = self.gembase2_parser(lst_path, self.replicon.id)
        else:
            msg = f"Unknow gembase format: {self.gembase_type}"
            _log.critical(msg)
            raise IntegronError(msg)
        return prots_info


    @classmethod
    def find_gembase_file_basename(cls, gembase_path, input_seq_path):
        """
        from the input file name, try to retrieve the basename which is used in gembase
        This specially useful when IF is run in parallel. The input sequence is split
        in chunks and treat in parallel. But in this case the name of the chunk does not match
        neither the lstinfo file nor the protein file.
        So this method try retrieve the original basename without extension
        for instance: ::

            ACBA.0917.00019.fna               => ACBA.0917.00019
            ACBA.0917.00019.0001.fna          => ACBA.0917.00019
            ESCO001.C.00001.C001.fst          => ESCO001.C.00001.C001
            ESCO001.C.00001.C001_chunk_1.fst  => ESCO001.C.00001.C001

        :return: the gembase basename corresponding to the input file
        :rtype: string
        """
        gembase_file_basename = os.path.splitext(os.path.basename(input_seq_path))[0]
        # when IF is run through nextflow & parallel_integron_finder
        # the input data is split the name of the chunks can vary
        # it can be
        # the name of input file with suffix _chunk_<id>
        # it can be the name of the contig
        # so wee need to find the original name to find the
        # - lstinfo file
        # - protein file
        match = re.search(r"_chunk_\d+$", gembase_file_basename)
        if match:
            gembase_file_basename = gembase_file_basename[:match.start()]

        lst_dir = cls.get_lst_dir(gembase_path)
        lst_path = os.path.join(lst_dir, gembase_file_basename + '.lst')
        if os.path.exists(lst_path):
            # it is a complete genome
            return gembase_file_basename
        else:
            # it is a contig
            # let's find the draft genome lst file
            gembase_file_basename = os.path.splitext(os.path.basename(gembase_file_basename))[0]
            lst_path = os.path.join(gembase_path, gembase_file_basename + '.lst')
            if os.path.exists(lst_path):
                return gembase_file_basename
            else:
                raise FileNotFoundError(f"cannot find lst file matching {input_seq_path} sequence")


    def _make_protfile(self, path=None):
        """
        Create fasta file with protein corresponding to this sequence, from the corresponding Gembase protfile
        This step is necessary because in Gembase1&2 Draft
        One nucleic file can contain several contigs, but all proteins are in the same file.
        and in Gembase2 replicons file can contain several chromosomes and plasmids
        :return: the path of the created protein file
        :rtype: str
        """
        if path:
            prot_file_path = path
        else:
            all_prot_path = os.path.join(self._gembase_path, 'Proteins', self._gembase_file_basename + '.prt')
            try:
                all_prots = SeqIO.index(all_prot_path, "fasta")
                if not os.path.exists(self.cfg.tmp_dir(self.replicon.id)):
                    os.makedirs(self.cfg.tmp_dir(self.replicon.id))
                prot_file_path = os.path.join(self.cfg.tmp_dir(self.replicon.id), self.replicon.id + '.prt')
                with open(prot_file_path, 'w') as prot_file:
                    for seq_id in self._info[4]:
                        try:
                            seq = all_prots[seq_id]
                            SeqIO.write(seq, prot_file, 'fasta')
                        except KeyError:
                            self._pseudo_genes.add(seq_id)
                            _log.warning(f'Sequence describe in LSTINF file {seq_id} is not present in {all_prot_path}')
            finally:
                all_prots.close()

        return prot_file_path


    @classmethod
    def gembase_sniffer(cls, lst_path):
        """
        Detect the type of gembase
        :param str lst_path: the path to the LSTINFO file corresponding to the nucleic sequence
        :returns: either the type of replicon ('Complet' or 'Draft', 1 or 2)
        :rtype: tuple
        """
        with open(lst_path) as lst_file:
            first_line = next(lst_file)

        gene_id = first_line.split()[4]
        guess_gb_type = False
        for gb_type, pattern in cls._gene_patterns.items():
            match = re.match(pattern, gene_id)
            if match:
                guess_gb_type = True
                break

        if not guess_gb_type:
            start, end , seq_id= first_line.split()[:3]
            if start == end == '0':
                msg = f"The genome {seq_id} seems empty: see {lst_path}"
                _log.critical(msg)
                raise IntegronError(msg) from None
            else:

                msg = f"Cannot detect GemBase version, check lst file '{lst_path}'."
                raise IntegronError(msg) from None
        _log.debug(f"GembaseDB sniff GemBase version:{gb_type}")
        return gb_type


    @property
    def gembase_type(self):
        return self._gembase_type

    @classmethod
    def get_replicon_type(cls, seq_id='', rep_id=''):
        """

        :param seq_id: the sequence id to parse
        :type seq_id: str
        :param rep_id: the replicon identifier
        :type rep_id: str
        :return: the kind of genome, it can be either:

            * Chromosome
            * Plasmid
            * Phage
            * Other
            * Draft

        :rtype: :class:`RepliconType` object
        """
        if not any((seq_id, rep_id)):
            raise IntegronError(f'{cls.__name__}.get_replicon_type you must provide either a seqid or a rep_id')
        elif all((seq_id, rep_id)):
            raise IntegronError(f'{cls.__name__}.get_replicon_type you must provide either a seqid or a rep_id')
        guess_gb_type = False
        if seq_id:
            patterns = cls._rep_patterns
            _id = seq_id
        elif rep_id:
            patterns = cls._gene_patterns.values()
            _id = rep_id
        for pattern in patterns:
            match = re.match(pattern, _id)
            if match:
                guess_gb_type = True
                try:
                    genome_type_letter = match.group('g_type')
                except IndexError:
                    # there is no group
                    # so a seq_id has been provided
                    # and it match a Draft
                    genome_type_letter = 'i'
                break
        if not guess_gb_type:
            msg = f"Cannot detect GemBase version, from {'seq_id' if seq_id else 'rep_id'}: '{_id}'."
            _log.error(msg)
            raise IntegronError(msg) from None
        else:
            genome_type = {
                'C': RepliconType.CHROMOSOME,
                'P': RepliconType.PLASMID,
                'V': RepliconType.PHAGE,
                'O': RepliconType.OTHER,
                'D': RepliconType.DRAFT, # V2_plus
                'i': RepliconType.DRAFT, # V2
                'b': RepliconType.DRAFT  # V2
            }[genome_type_letter]
        return genome_type


    @staticmethod
    def gembase1_complete_parser(lst_path, sequence_id):
        """
        :param str lst_path: the path of the LSTINFO file Gembase Complet
        :param str sequence_id: the id of the genomic sequence to analyse
        :return: the information related to the 'valid' CDS corresponding to the sequence_id
        :rtype: `class`:pandas.DataFrame` object
        """
        dtype = {0: 'int', # start
                 1: 'int', # end
                 2: 'str', # strand C/D
                 3: 'str', # type (CDS, ncRNA, tRNA,...)
                 4: 'str', # seq id
                 5: 'str', # Valid
                 6: 'str', # gene name
                 7: 'str'} # description

        with open(lst_path) as lst_file:
            lst_data = []
            for line in lst_file:
                start, end, strand, gene_type, seq_id, valid, gene_name, *description = line.strip().split()
                row = [start, end, strand, gene_type, seq_id, valid, gene_name, ' '.join(description)]
                lst_data.append(row)

            lst = pd.DataFrame(lst_data)
            lst = lst.astype(dtype)
            genome_info = lst.loc[lst[4].str.contains(sequence_id, regex=True)]
            prots_info = genome_info.loc[(genome_info[3] == 'CDS') & (genome_info[5] == 'Valid')]
            return prots_info


    @staticmethod
    def gembase1_draft_parser(lst_path, replicon_id):
        """
        :param str lst_path: the path of the LSTINFO file from a Gembase Draft
        :param str sequence_id: the id of the genomic sequence to analyse
        :return: the information related to the 'valid' CDS corresponding to the sequence_id
        :rtype: `class`:pandas.DataFrame` object
        """
        try:
            lst = pd.read_csv(lst_path,
                              header=None,
                              sep="\t")
        except Exception as err:
            msg = f"Error while parsing {lst_path} file: {err}"
            _log.error(msg)
            raise IntegronError(msg)
        specie, date, strain, contig = replicon_id.split('.')
        pattern = fr'{specie}\.{date}\.{strain}\.[bi]{contig}'
        try:
            genome_info = lst.loc[lst[4].str.contains(pattern, regex=True)]
        except KeyError:
            msg = f"The LST file '{lst_path}' seems not to be in gembase V1 draft format."
            _log.error(msg)
            raise IntegronError(msg) from None
        prots_info = genome_info.loc[genome_info[3] == 'CDS']
        return prots_info


    @staticmethod
    def gembase2_parser(lst_path, replicon_id):
        """
        :param str lst_path: the path of the LSTINFO file from a Gembase Draft
        :param str sequence_id: the id of the genomic sequence to analyse
        :return: the information related to the 'valid' CDS corresponding to the sequence_id
        :rtype: `class`:pandas.DataFrame` object
        """
        try:
            lst = pd.read_csv(lst_path,
                              header=None,
                              sep="\t")
        except Exception as err:
            msg = f"Error while parsing {lst_path} file: {err}"
            _log.error(msg)
            raise IntegronError(msg)

        try:
            lst = lst.loc[:, :7]
            dtype = {i: 'str' for i in range(lst.shape[1])}
            dtype[0] = dtype[1] = 'int'
            lst = lst.astype(dtype)
            genome_info = lst.loc[lst[4].str.contains(replicon_id, regex=True)]
            prots_info = genome_info.loc[genome_info[3] == 'CDS']
        except Exception:
            msg = f"The LST file '{lst_path}' seems not to be in gembase V2 draft format."
            _log.error(msg)
            raise IntegronError(msg) from None
        return prots_info


    def __getitem__(self, prot_seq_id):
        """
        :param str prot_seq_id: the id of a protein sequence
        :return: The Sequence corresponding to the prot_seq_id.
        :rtype: :class:`Bio.SeqRecord` object
        """
        return self._prot_db[prot_seq_id]


    def __iter__(self):
        """
        :return: a generator which iterate on the gene seq_ids which constitute the contig (genes and pseudogenes).
        :rtype: generator
        """
        return (seq_id for seq_id in self._info[4])


    def coding_prot_ids(self):
        """

        :return: a generator which iterate on coding genes seq_ids (the non coding genes are discarded)
        :rtype: generator
        """
        all_seq_ids = self._info[4]
        coding_seqid = all_seq_ids[~all_seq_ids.isin(self._pseudo_genes)]
        return (seq_id for seq_id in coding_seqid)


    def get_description(self, gene_id):
        """
        :param str gene_id: a protein/gene identifier
        :return: The description of the protein corresponding to the gene_id
        :rtype: :class:`SeqDesc` namedtuple object
        :raise IntegronError: when gene_id is not a valid Gembase gene identifier
        :raise KeyError: if gene_id is not found in GembaseDB instance
        """
        try:
            specie, date, strain, contig_gene = gene_id.split('.')
            contig_gene = contig_gene[1:]  # remove the first letter b/i
        except ValueError:
            raise IntegronError(f"'{gene_id}' is not a valid Gembase protein identifier.")
        pattern = fr'{specie}\.{date}\.{strain}\.\w?{contig_gene}'
        seq_info = self._info.loc[self._info[4].str.contains(pattern, regex=True)]
        if not seq_info.empty:
            return SeqDesc(seq_info[4].values[0],
                           1 if seq_info[2].values[0] == "D" else -1,
                           seq_info[0].values[0],
                           seq_info[1].values[0],
                           )
        else:
            raise KeyError(gene_id)





class ProdigalDB(ProteinDB):
    """
    Creates proteins from Replicon/contig using prodigal and provide facilities to access them.
    """


    def _make_protfile(self, path=None):
        """
        Use `prodigal` to generate proteins corresponding to the replicon

        :return: the path of the created protfile
        :rtype: str
        """
        assert self.cfg.prodigal, "'prodigal' not found."
        if path:
            prot_file_path = path
        else:
            if not os.path.exists(self.cfg.tmp_dir(self.replicon.id)):
                os.makedirs(self.cfg.tmp_dir(self.replicon.id))
            prot_file_path = os.path.join(self.cfg.tmp_dir(self.replicon.id), self.replicon.id + ".prt")
            if not os.path.exists(prot_file_path):
                prodigal_cmd = '{prodigal} {meta} -i {replicon} -a {prot} -o {out} -q '.format(
                    prodigal=self.cfg.prodigal.replace(' ', '\\ '),
                    meta='' if len(self.replicon) > 200000 else '-p meta',
                    replicon=self.replicon.path.replace(' ', '\\ '),
                    prot=prot_file_path.replace(' ', '\\ '),
                    out=os.devnull,
                )
                try:
                    _log.debug("run prodigal: {}".format(prodigal_cmd))
                    completed_process = subprocess.run(shlex.split(prodigal_cmd))
                except Exception as err:
                    raise RuntimeError(f"{prodigal_cmd} : failed : {err}")
                if completed_process.returncode != 0:
                    raise RuntimeError(f"{prodigal_cmd} : failed : prodigal returncode = {completed_process.returncode}")

        return prot_file_path


    def __getitem__(self, prot_seq_id):
        """
        :param str prot_seq_id: the id of a protein sequence
        :return: The Sequence corresponding to the prot_seq_id.
        :rtype: :class:`Bio.SeqRecord` object
        """
        try:
            return self._prot_db[prot_seq_id]
        except KeyError:
            raise IntegronError(f"protein file does not contains '{prot_seq_id}' id. "
                                f"Try again with removing previous results dir {self.cfg.result_dir}")

    def __iter__(self):
        """
        :return: a generator which iterate on the protein seq_id which constitute the contig.
        :rtype: generator
        """
        return (seq_id for seq_id in self._prot_db)


    def coding_prot_ids(self):
        """
        :return: a generator which iterate on coding genes seq_ids (the non coding genes are discarded)
        :rtype: generartor
        """
        return self.__iter__()


    def get_description(self, gene_id):
        """
        :param str gene_id: a protein/gene identifier
        :returns: The description of the protein corresponding to the gene_id
        :rtype: :class:`SeqDesc` namedtuple object
        :raise IntegronError: when gene_id is not a valid Gembase gene identifier
        :raise KeyError: if gene_id is not found in ProdigalDB instance
        """
        seq = self[gene_id]
        try:
            id_, start, stop, strand, *_ = seq.description.split(" # ")
        except ValueError:
            raise IntegronError(f"'{gene_id}' is not a valid Prodigal protein identifier.")
        start = int(start)
        stop = int(stop)
        strand = int(strand)
        return SeqDesc(id_, strand, start, stop)


class CustomDB(ProteinDB):
    """
    Creates proteins from Replicon/contig using prodigal and provide facilities to access them.
    """


    def __init__(self, replicon, cfg, prot_file):
        super().__init__(replicon, cfg, prot_file=prot_file)
        try:
            parser_path = self.cfg.annot_parser
            spec = importlib.util.spec_from_file_location('custom_module', parser_path)
            custom_module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(custom_module)
            self._parser = custom_module.description_parser
        except Exception as err:
            raise RuntimeError(f"Cannot import custom --annot-parser '{parser_path}': {err}")


    def _make_protfile(self, path=None):
        if path is None:
            raise IntegronError("If use CustomDB prot_file must be specified")
        return path


    def __getitem__(self, prot_seq_id):
        """
        :param str prot_seq_id: the id of a protein sequence
        :return: The Sequence corresponding to the prot_seq_id.
        :rtype: :class:`Bio.SeqRecord` object
        """
        try:
            return self._prot_db[prot_seq_id]
        except KeyError:
            raise IntegronError(f"protein file does not contains '{prot_seq_id}' id. "
                                f"Check if it's the right proteins file {self._prot_file} "
                                f"or remove previous results dir {self.cfg.result_dir}")


    def __iter__(self):
        """
        :return: a generator which iterate on the protein seq_id which constitute the contig.
        :rtype: generator
                """
        return (seq_id for seq_id in self._prot_db)


    def coding_prot_ids(self):
        """
        :return: a generator which iterate on coding genes seq_ids (the non coding genes are discarded)
        :rtype: generartor
        """
        return self.__iter__()


    def get_description(self, gene_id):
        """
        :param str gene_id: a protein/gene identifier
        :returns: The description of the protein corresponding to the gene_id
        :rtype: :class:`SeqDesc` namedtuple object
        :raise IntegronError: when gene_id is not a valid Gembase gene identifier
        :raise KeyError: if gene_id is not found in ProdigalDB instance
        """
        def check_id(x):
            return isinstance(x, str)

        def check_strand(s):
            return s == 1 or s == -1

        def check_start(s):
            return isinstance(s, int) and s >= 0

        check_stop = check_start

        seq = self[gene_id]

        try:
            id_, start, stop, strand = self._parser(seq.description)
        except ValueError:
            msg = f"'{gene_id}' protein is not compliant with custom --annot-parser '{self.cfg.annot_parser}'."
            _log.critical(msg)
            raise IntegronError(msg)
        except Exception as err:
            msg = f"Cannot parse protein file '{self._prot_file}' with annot-parser '{self.cfg.annot_parser}': {err}"
            _log.critical(msg)
            raise IntegronError(msg) from None

        if not all((check_id(id_), check_start(start), check_stop(stop), check_strand(strand))):
            msg = "Error during protein file parsing: expected seq_id: str, start: positive int, stop: positive int, " \
                  f"strand 1/-1. got: {id_}, {start}, {stop}, {strand}"
            _log.critical(msg)
            raise IntegronError(msg)
        return SeqDesc(id_, strand, start, stop)
