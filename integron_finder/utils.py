# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright © 2015 - 2018  Institut Pasteur, Paris.                                #
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
from collections import namedtuple
from Bio import Seq
from Bio import SeqIO


def make_single_fasta_reader(alphabet):

    def read_fasta(path):
        """
        :param path:the path to the fasta file
        :return: the sequence parsed
        :rtype: :class:`Bio.SeqRecord.SeqRecord` object
        """
        seq = SeqIO.read(path, "fasta",  alphabet=alphabet)
        seq.name = get_name_from_path(path)
        return seq
    return read_fasta


read_single_dna_fasta = make_single_fasta_reader(Seq.IUPAC.unambiguous_dna)

read_single_prot_fasta = make_single_fasta_reader(Seq.IUPAC.protein)


def make_multi_fasta_reader(alphabet):
    """
    fasta generator maker

    :param alphabet: the alphabet store in the fasta generator closure
    :return: generator to iterate on the fasta file in the same order as in fasta file
    """

    def fasta_iterator(path):
        """
        :param path:the path to the fasta file
        :return: the sequence parsed
        :rtype: :class:`Bio.SeqRecord.SeqRecord` object
        """
        name = get_name_from_path(path)
        seq_it = SeqIO.parse(path, "fasta",  alphabet=alphabet)
        for seq in seq_it:
            seq.name = name
            yield seq

    return fasta_iterator


read_multi_prot_fasta = make_multi_fasta_reader(Seq.IUPAC.protein)


class FastaIterator(object):
    """
    Allow to parse over a multi fasta file, and iterate over it
    to use it  create a reader with an alphabet, then iterate over it::

        read_multi_dna_fasta = FastaIterator(Seq.IUPAC.unambiguous_dna)
        sequences_db = read_multi_dna_fasta(fasta_file)
        print len(sequences_db)
        for seq in sequence_db:
            print seq.id

    .. warning:

        The order of sequences is not guarantee

    """

    def __init__(self, alphabet):
        """
        :param alphabet: The authorized alphabet
        :type alpahbet: Bio.SeqIUPAC member
        """
        self.alphabet = alphabet
        self.replicon_name = None
        self.topologies = None
        self.seq_index = None
        self.dist_threshold = 4000

    def __iter__(self):
        """
        allow to iterate over the sequences
        :return: a generator at each iteration returns a :class:`Bio.Seq.SeqRecord` object
        """
        for id_ in self.seq_index.keys():
            seq = self.seq_index[id_]
            if self.replicon_name is not None:
                seq.name = self.replicon_name
            if self.topologies:
                topology = self.topologies[seq.id]
                # If sequence is too small, it can be problematic when using circularity
                if topology == 'circ' and len(seq) <= 4 * self.distance_threshold:
                    topology = 'lin'
                seq.topology = topology
            else:
                seq.topology = 'circ'

            yield seq

    def __len__(self):
        """:returns: The nuber of sequence in the file"""
        return len(self.seq_index)

    def __call__(self, path, replicon_name=None, topologies=None, dist_threshold=4000):
        """
        :param path:the path to the fasta file
        :param replicon_name: the name to inject in each SeqRecord.name (by default is the same a SeqRecord.id)
        :return: the sequence parsed
        :rtype: :class:`Bio.SeqRecord.SeqRecord` object
        """

        if replicon_name is not None:
            self.replicon_name = replicon_name
        if topologies is not None:
            self.topologies = topologies
            self.dist_threshold = dist_threshold
        self.seq_index = SeqIO.index(path, "fasta", alphabet=self.alphabet)
        return self


read_multi_dna_fasta = FastaIterator(Seq.IUPAC.unambiguous_dna)


def model_len(path):
    """

    :param str path: the path to the covariance model file
    :return: the length of the model
    :rtype: int
    """
    if not os.path.exists(path):
        raise IOError("Path to model_attc '{}' does not exists".format(path))
    with open(path) as model_file:
        for line in model_file:
            if line.startswith('CLEN'):
                model_length = int(line.split()[1])
                return model_length
        raise RuntimeError("CLEN not found in '{}', maybe it's not infernal model file".format(path))


def get_name_from_path(path):
    """
    :param path: The path to extract name for instance the fasta file to the replicon
    :return: the name of replicon for instance
             if path = /path/to/replicon.fasta name = repliocn
    """
    return os.path.splitext(os.path.split(path)[1])[0]

"""Sequence description with fields: id strand start stop"""
SeqDesc = namedtuple('SeqDesc', ('id', 'strand', 'start', 'stop'))


def gembase_parser(description):
    """"
    :param description: description (1rst line without id) of sequence from gembase fasta file
    :return: SeqDesc
    """
    desc = description.split(" ")
    id_, strand, start, stop = desc[:2] + desc[4:6]
    strand = 1 if desc[1] == "D" else -1
    start = int(start)
    stop = int(stop)
    return SeqDesc(id_, strand, start, stop)


def non_gembase_parser(description):
    """
    :param description: description (1rst line without id) of sequence from fasta file not coming a gemabse
    :return: SeqDesc
    """
    id_, start, stop, strand, _ = description.split(" # ")
    start = int(start)
    stop = int(stop)
    strand = int(strand)
    return SeqDesc(id_, strand, start, stop)
