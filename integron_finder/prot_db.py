# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright (c) 2015 - 2018  Institut Pasteur, Paris and CNRS.                     #
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

import os
from subprocess import call
import colorlog
from collections import namedtuple
import pandas as pd

from Bio import SeqIO, Alphabet

_log = colorlog.getLogger(__name__)

"""Sequence description with fields: id strand start stop"""
SeqDesc = namedtuple('SeqDesc', ('id', 'strand', 'start', 'stop'))

class ProteinDB():

    def __init__(self):
        pass

    def __getitem__(self, seq_id):
        pass

    def __iter__(self):
        pass

    def _make_protfile(self):
        pass

    def _make_db(self):
        pass

    def get_description(self, seq_id):
        pass

    def protfile(self):
        pass


class Gembase(ProteinDB):

    def __init__(self, replicon, cfg):
        self.cfg = cfg
        self.replicon = replicon
        self._info = self.parse_lst()
        self._prot_db = self._make_db()
        self._prot_file = self._make_protfile()


    def _make_protfile(self):
        """
        Create fasta file with protein corresponding to this sequence, from the corresponding Gembase protfile
        This step is necessary because in Gembase Draft
        One nucleic file can contains several contigs, but all proteins are in the same file.

        :return: the path of the created protfile
        :rtype: str
        """
        all_prots = SeqIO.index(self._prot_file, "fasta", alphabet=Alphabet.ProteinAlphabet)
        protfile_path = os.path.join(self.cfg.tmp_res, self.replicon.id + '.prt')
        with open(protfile_path, 'w') as protfile:
            for seq_id in self._info['seq_id']:
                seq = all_prots[seq_id]
                SeqIO.write(seq, protfile, 'fasta')
        return protfile_path


    def _parse_lst(self):
        """

        :return:
        """
        lst_file = os.path.join(self.cfg.gembase, 'LSTINFO', '{}.lst'.format(self.replicon.id))
        lst = pd.read_table(lst_file,
                            header=None,
                            names=['start', 'end', 'strand', 'type', 'seq_id'],
                            usecols=(0, 1, 2, 3, 4),
                            dtype={'start': 'int',
                                   'end': 'int',
                                   'strand': 'str',
                                   'type': 'str',
                                   'seq_id': 'str'}
                            )
        specie, date, strain, contig = self.replicon.id.split('.')
        pattern = '{}\.{}\.{}\.\w?{}'.format(specie, date, strain, contig)
        genome_info = lst[lst['seq_id'].str.contains(pattern, regex=True)]
        prots_info = genome_info[genome_info['type'] == 'CDS']
        return prots_info


    def _make_db(self):
        """
        :return: an index of the sequence contains in protfile corresponding to the replicon
        """
        return SeqIO.index(self._prot_file, "fasta", alphabet=Alphabet.ProteinAlphabet)


    def __getitem__(self, prot_seq_id):
        """

        :param str protseq_id: the id of a protein sequence
        :return: The Sequence corresponding to the prot_seq_id.
        :rtype: :class:`Bio.SeqRecord` object
        """
        return self._prot_db[prot_seq_id]


    def __iter__(self):
        """
        :return: a generator which iterate on the protein seq_id which constitute
                 the contig
        :rtype: generator
        """
        return (_id for _id in self._info['seq_id'])


    @property
    def protfile(self):
        """
        :return: The absolute path to the protein file corresponding to contig id
        :rtype: str
        """
        return self._prot_file


    def get_description(self, gene_id):
        """

        :param str seq_id: the sequence identifier
        :return: SeqDesc object
        """
        specie, date, strain, contig, gene = gene_id.split('.')
        pattern = '{}\.{}\.{}\.\w?{}.{}'.format(specie, date, strain, contig, gene)
        seq_info = self._info.loc[self._info['seq_id'].str.contains(pattern, regex=True)]
        return SeqDesc(seq_info.seq_id.values[0],
                       1 if seq_info.strand.values[0] == "D" else -1,
                       seq_info.start.values[0],
                       seq_info.stop.values[0],
                       )


class ProdigalDB(ProteinDB):

    def __init__(self, replicon, cfg):
        self.cfg = cfg
        self.replicon = replicon
        self._db = self._make_db()
        self._prot_file = self._make_protfile()


    def _make_protfile(self):
        """
        Use `prodigal` to generate proteins corresponding to the replicon

        :return: the path of the created protfile
        :rtype: str
        """
        prot_file_path = os.path.join(self.cfg.tmp_dir, self.replicon.id + ".prt")
        if not os.path.isfile(prot_file_path):
            prodigal_cmd = "{prodigal} {meta} -i {replicon} -a {prot} -o {out} -q ".format(
                prodigal=self.cfg.prodigal,
                meta='' if len(self.replicon) > 200000 else '-p meta',
                replicon=self.replicon_path,
                prot=prot_file_path,
                out=os.devnull,
            )
            try:
                _log.debug("run prodigal: {}".format(prodigal_cmd))
                returncode = call(prodigal_cmd.split())
            except Exception as err:
                raise RuntimeError("{0} failed : {1}".format(prodigal_cmd, err))
            if returncode != 0:
                raise RuntimeError("{0} failed returncode = {1}".format(prodigal_cmd, returncode))

        return prot_file_path


    def _make_db(self):
        return SeqIO.index(self._prot_file, "fasta", alphabet=Alphabet.ProteinAlphabet)


    def __getitem__(self, seq_id):
        """

        :param seq_id:
        :return: :class:`Bio.SeqRecord` object
        """
        return self._db[seq_id]


    def __iter__(self):
        """
        :return: a generator which iterate on the protein seq_id which constitute
                 the contig
        :rtype: generator
                """
        return (seq_id for seq_id in self._db)

    @property
    def protfile(self):
        return self._prot_file


    def get_description(self, seq_id):
        seq = self[seq_id]
        id_, start, stop, strand, *_ = seq.description.split(" # ")
        start = int(start)
        stop = int(stop)
        strand = int(strand)
        return SeqDesc(id_, strand, start, stop)

