# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright (c) 2015 - 2021  Institut Pasteur, Paris and CNRS.                     #
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
from subprocess import call
from collections import namedtuple
import re
import importlib.util

import colorlog
import pandas as pd
from Bio import SeqIO, Seq
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
        try:
            # for biopython < 1.78
            idx = SeqIO.index(self._prot_file, "fasta", alphabet=Seq.IUPAC.extended_protein)
        except AttributeError:
            # for biopython > 1.76
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


class GembaseDB(ProteinDB):
    """
    Implements :class:`ProteinDB` from a Gembase.
    Managed proteins from Proteins directory corresponding to a replicon/contig
    """


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
        # in GemBase Draft the files ar based on replicon id
        # but one file can contains several contig
        # the sequence id contains the contig number
        # for instance ACBA.0917.00019 vs ACBA.0917.00019.0001
        # the filenames are based on replicon id
        # but in code the replicon id contains the contig number
        # for gembase complete both filename and seq_id contains contig number
        self._replicon_name = os.path.splitext(os.path.basename(self.cfg.input_seq_path))[0]
        self._gembase_file_basename = self._find_gembase_file_basename(self._gembase_path, self.cfg.input_seq_path)
        self._info = self._parse_lst()
        if prot_file is None:
            self._prot_file = self._make_protfile()
        else:
            self._prot_file = prot_file
        self._prot_db = self._make_db()


    def _find_gembase_file_basename(self, gembase_path, input_seq_path):
        """
        from the input file name, try to retrieve the basename which is used in gembase
        This specially useful when IF is run in parallel. The input sequence is split
        in chunks and treat in parallel. But in this case the name of the chunk does not math
        neither the lstinfo file nor the protein file.
        So this method try retrieve the original basename without extension
        for instance: ::

            ACBA.0917.00019.fna               => ACBA.0917.00019
            ACBA.0917.00019.0001.fna          => ACBA.0917.00019
            ESCO001.C.00001.C001.fst          => ESCO001.C.00001.C001
            ESCO001.C.00001.C001_chunk_1.fst  => ESCO001.C.00001.C001

        :return: the gemabse basename corresponding to the input file
        :rtype: string
        """
        lst_dir_path = os.path.join(gembase_path, 'LSTINFO')
        if not os.path.exists(lst_dir_path):
            _log.warning(f"{lst_dir_path} directory nout found. check LSTINF")
            lst_dir_path = os.path.join(gembase_path, 'LSTINF')
            if not os.path.exists(lst_dir_path):
                _log.warning(f"{lst_dir_path} directory nout found.")
                raise IntegronError(f"{lst_dir_path}O directory nout found.")
        gembase_file_basename = os.path.splitext(os.path.basename(input_seq_path))[0]
        # when IF is run through nextflow & parallel_integron_finder
        # the input data is split the name of the chunks can vary
        # it can be
        # the name of input file with suffix _chunk_<id>
        # it can be the name of the contig
        # so wee need to find the original name to find the
        # - lstinfo file
        # - protein file
        match = re.search("_chunk_\d+$", gembase_file_basename)
        if match:
            gembase_file_basename = gembase_file_basename[:match.start()]

        lst_path = os.path.join(lst_dir_path, gembase_file_basename + '.lst')
        if os.path.exists(lst_path):
            # it is a complete genome
            return gembase_file_basename
        else:
            # it is a contig
            # let's find the draft genome lst file
            gembase_file_basename = os.path.splitext(os.path.basename(gembase_file_basename))[0]
            lst_path = os.path.join(lst_dir_path, gembase_file_basename + '.lst')
            if os.path.exists(lst_path):
                return gembase_file_basename
            else:
                raise FileNotFoundError("cannot find lst file matching {} sequence".format(self.cfg.input_seq_path))


    def _make_protfile(self, path=None):
        """
        Create fasta file with protein corresponding to this sequence, from the corresponding Gembase protfile
        This step is necessary because in Gembase Draft
        One nucleic file can contains several contigs, but all proteins are in the same file.

        :return: the path of the created protein file
        :rtype: str
        """
        if path:
            prot_file_path = path
        else:
            all_prot_path = os.path.join(self._gembase_path, 'Proteins', self._gembase_file_basename + '.prt')
            try:
                all_prots = SeqIO.index(all_prot_path, "fasta", alphabet=Seq.IUPAC.extended_protein)
            except AttributeError:
                all_prots = SeqIO.index(all_prot_path, "fasta")
            if not os.path.exists(self.cfg.tmp_dir(self.replicon.id)):
                os.makedirs(self.cfg.tmp_dir(self.replicon.id))
            prot_file_path = os.path.join(self.cfg.tmp_dir(self.replicon.id), self.replicon.id + '.prt')
            with open(prot_file_path, 'w') as prot_file:
                for seq_id in self._info['seq_id']:
                    try:
                        seq = all_prots[seq_id]
                        SeqIO.write(seq, prot_file, 'fasta')
                    except KeyError:
                        _log.warning('Sequence describe in LSTINF file {} is not present in {}'.format(seq_id, all_prot_path))
        return prot_file_path


    @staticmethod
    def gembase_sniffer(lst_path):
        """
        Detect the type of gembase
        :param str lst_path: the path to the LSTINFO file corresponding to the nucleic sequence
        :returns: either 'Complet' or 'Draft'
        """
        with open(lst_path) as lst_file:
            line = lst_file.readline()
        fields = line.split()
        if fields[5] in ('Valid', 'Invalid_Size', 'Pseudo', 'Partial'):
            return 'Complet'
        else:
            if fields[0] == '0' and fields[1] == '0':
                msg = f"the genome {fields[2]} seems empty: see {lst_path}"
                _log.critical(msg)
                raise IntegronError(msg) from None
            else:
                return 'Draft'


    @staticmethod
    def gembase_complete_parser(lst_path, sequence_id):
        """
        :param str lst_path: the path of of the LSTINFO file Gembase Complet
        :param str sequence_id: the id of the genomic sequence to analyse
        :return: the information related to the 'valid' CDS corresponding to the sequence_id
        :rtype: `class`:pandas.DataFrame` object
        """
        dtype = {'start': 'int',
                 'end': 'int',
                 'strand': 'str',
                 'type': 'str',
                 'seq_id': 'str',
                 'valid': 'str',
                 'gene_name': 'str',
                 'description': 'str'}

        with open(lst_path) as lst_file:
            lst_data = []
            for line in lst_file:
                start, end, strand, gene_type, seq_id, valid, gene_name, *description = line.strip().split()
                row = [start, end, strand, gene_type, seq_id, valid, gene_name, ' '.join(description)]
                lst_data.append(row)

            lst = pd.DataFrame(lst_data,
                               columns=['start', 'end', 'strand', 'type', 'seq_id', 'valid', 'gene_name', 'description']
                               )
            lst = lst.astype(dtype)
            specie, date, strain, contig = sequence_id.split('.')
            pattern = '{}\.{}\.{}\.\w?{}'.format(specie, date, strain, contig)
            genome_info = lst.loc[lst['seq_id'].str.contains(pattern, regex=True)]
            prots_info = genome_info.loc[(genome_info['type'] == 'CDS') & (genome_info['valid'] == 'Valid')]
            return prots_info


    @staticmethod
    def gembase_draft_parser(lst_path, replicon_id):
        """
        :param str lst_path: the path of of the LSTINFO file from a Gembase Draft
        :param str sequence_id: the id of the genomic sequence to analyse
        :return: the information related to the 'valid' CDS corresponding to the sequence_id
        :rtype: `class`:pandas.DataFrame` object
        """
        try:
            lst = pd.read_csv(lst_path,
                              header=None,
                              names=['start', 'end', 'strand', 'type', 'seq_id', 'gene_name', 'description'],
                              dtype={'start': 'int',
                                     'end': 'int',
                                     'strand': 'str',
                                     'type': 'str',
                                     'seq_id': 'str',
                                     'gene_name': 'str',
                                     'description': 'str'},
                              sep="\t")
        except Exception as err:
            _log.error(f"Error while parsing {lst_path} file")
            raise err
        specie, date, strain, contig_gene = replicon_id.split('.')
        pattern = f'{specie}\.{date}\.{strain}\.\w?{contig_gene}'
        genome_info = lst.loc[lst['seq_id'].str.contains(pattern, regex=True)]
        prots_info = genome_info.loc[genome_info['type'] == 'CDS']
        return prots_info


    def _parse_lst(self):
        """
        Parse the LSTINFO file and extract information specific to the replicon
        :return: `class`:pandas.DataFrame` object
        """
        lst_path = os.path.join(self._gembase_path, 'LSTINFO', self._gembase_file_basename + '.lst')
        if not os.path.exists(lst_path):
            lst_path = os.path.join(self._gembase_path, 'LSTINF', self._gembase_file_basename + '.lst')
            if not os.path.exists(lst_path):
                raise IntegronError(f"{lst_path}O directory not found.")
        gembase_type = self.gembase_sniffer(lst_path)
        if gembase_type == 'Draft':
            prots_info = self.gembase_draft_parser(lst_path, self.replicon.id)
        else:
            prots_info = self.gembase_complete_parser(lst_path, self.replicon.id)
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
        :return: a generator which iterate on the protein seq_id which constitute the contig.
        :rtype: generator
        """
        return (seq_id for seq_id in self._info['seq_id'])


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
        pattern = f'{specie}\.{date}\.{strain}\.\w?{contig_gene}'
        seq_info = self._info.loc[self._info['seq_id'].str.contains(pattern, regex=True)]
        if not seq_info.empty:
            return SeqDesc(seq_info.seq_id.values[0],
                           1 if seq_info.strand.values[0] == "D" else -1,
                           seq_info.start.values[0],
                           seq_info.end.values[0],
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
        if path:
            prot_file_path = path
        else:
            if not os.path.exists(self.cfg.tmp_dir(self.replicon.id)):
                os.makedirs(self.cfg.tmp_dir(self.replicon.id))
            prot_file_path = os.path.join(self.cfg.tmp_dir(self.replicon.id), self.replicon.id + ".prt")
            if not os.path.exists(prot_file_path):
                prodigal_cmd = "{prodigal} {meta} -i {replicon} -a {prot} -o {out} -q ".format(
                    prodigal=self.cfg.prodigal,
                    meta='' if len(self.replicon) > 200000 else '-p meta',
                    replicon=self.replicon.path,
                    prot=prot_file_path,
                    out=os.devnull,
                )
                try:
                    _log.debug("run prodigal: {}".format(prodigal_cmd))
                    returncode = call(prodigal_cmd.split())
                except Exception as err:
                    raise RuntimeError(f"{prodigal_cmd} : failed : {err}")
                if returncode != 0:
                    raise RuntimeError(f"{prodigal_cmd} : failed : prodigal returncode = {returncode}")

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

