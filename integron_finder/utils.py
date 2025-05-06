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

import os

import colorlog
from Bio import Seq
from Bio import SeqIO

_log = colorlog.getLogger(__name__)


class MultiFastaReader:
    """
    Iterate over sequence of a fasta file
    """
    def __init__(self, path):
        """
        :param str path: the path to the fasta file
        """
        self.path = path
        self.name = get_name_from_path(path)
        self._seq_it = SeqIO.parse(path, "fasta")

    def __iter__(self):
        return self

    def __next__(self):
        """

        :return: The next sequence in the file
        :rtype: :class:`Bio.Seq.SeqRecord` object
        """
        seq = next(self._seq_it)
        seq.name = self.name
        return seq


    def __enter__(self):
        return self


    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


    def __del__(self):
        self.close()


    def close(self):
        """
        Close the file handle on the fasta file
        If it iterate until the end of the file, the file handle is automatically closed
        But if it stop before the end it is the responsability of the caller to close it.
        :return: None
        """
        # close the stream if the object is deleted before to iterate until the end
        if not self._seq_it.stream.closed:
            self._seq_it.stream.close()


class FastaIterator:
    """
    Allow to parse over a multi fasta file, and iterate over it

    .. warning::

        **The sequences order is not guarantee.**

    """

    def __init__(self, path, replicon_name=None, dist_threshold=4000, alphabet=None):
        #def __init__(self, path, alphabet=Seq.IUPAC.ambiguous_dna, replicon_name=None, dist_threshold=4000):
        """

        :param str path: The path to the file containing the sequences.
        :param alphabet: The authorized alphabet
        :type alphabet: Bio.SeqIUPAC member
        :param str replicon_name: The name of the replicon, if this specify all sequence.name will have this value
        :param int dist_threshold: The minimum length for a replicon to be considered as circular.
                                   Under this threshold even the provided topology is 'circular'
                                   the computation will be done with a 'linear' topology.
        """
        try:
            self.alphabet = Seq.IUPAC.ambiguous_dna
        except AttributeError:
            self.alphabet = None
        if self.alphabet:
            self.seq_index = SeqIO.index(path, "fasta",  alphabet=self.alphabet)
        else:
            self.seq_index = SeqIO.index(path, "fasta")
        self.seq_gen = (self.seq_index[id_] for id_ in self.seq_index.keys())
        self._topologies = None
        self.replicon_name = replicon_name
        self.dist_threshold = dist_threshold


    def _set_topologies(self, topologies):
        """

        :param topologies:
        :type topologies: :class:`integron_finder.Topology` onject
        :return:
        """
        self._topologies = topologies

    topologies = property(fset=_set_topologies)


    def _check_seq_alphabet_compliance(self, seq):
        """
        :param seq: the sequence to check
        :type seq: :class:`Bio.Seq.Seq` instance
        :return: True if sequence letters are a subset of the alphabet, False otherwise.
        """
        seq_letters = set(str(seq).upper())
        # the Bio.Alphabet has been removed from Biopython. from v1.78
        # I hard coded the IUPAC.ambiguous_dna
        # for compatibilty reasons
        alphabet = set('GATCRYWSMKHBVDN')
        return seq_letters.issubset(alphabet)

    def __next__(self):
        """
        :return: The next sequence (the order of sequences is not guaranteed).
        :rtype: a :class:`Bio.SeqRecord` object or None if the sequence is not compliant with the alphabet.
        """
        try:
            seq = next(self.seq_gen)
        except StopIteration as err:
            self.close()
            raise err from None
        if not self._check_seq_alphabet_compliance(seq.seq):
            _log.warning("sequence {} contains invalid characters, the sequence is skipped.".format(seq.id))
            return None
        if len(seq) < 50:
            _log.warning("sequence {} is too short ({} bp), the sequence is skipped (must be > 50bp).".format(seq.id,
                                                                                                              len(seq)))
            return None
        if self.replicon_name is not None:
            seq.name = self.replicon_name

        topology = None
        if self._topologies:
            topology = self._topologies[seq.id]
            # If sequence is too small, it can be problematic when using circularity
            if topology == 'circ' and len(seq) <= 4 * self.dist_threshold:
                _log.info(f"replicon '{seq.id}' length is too short ({len(seq)}<{4 * self.dist_threshold}) "
                          f"to be 'circular' set topology to 'linear'")
                topology = 'lin'

        seq.topology = topology

        if self.alphabet is None:
            seq.annotations["molecule_type"] = 'DNA'
        return seq


    def __iter__(self):
        return self


    def __len__(self):
        """:returns: The number of sequence in the file"""
        return len(self.seq_index)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self.seq_index.close()


def model_len(path):
    """

    :param str path: the path to the covariance model file
    :return: the length of the model
    :rtype: int
    """
    if not os.path.exists(path):
        msg = "Path to model_attc '{}' does not exists".format(path)
        _log.critical(msg)
        raise IOError(msg)
    with open(path) as model_file:
        for line in model_file:
            if line.startswith('CLEN'):
                model_length = int(line.split()[1])
                return model_length
        msg = "CLEN not found in '{}', maybe it's not infernal model file".format(path)
        _log.critical(msg)
        raise RuntimeError(msg)


def get_name_from_path(path):
    """
    :param path: The path to extract name for instance the fasta file to the replicon
    :return: the name of replicon for instance
             if path = /path/to/replicon.fasta name = replicon
    """
    return os.path.splitext(os.path.split(path)[1])[0]


def log_level(verbose, quiet):
    """
    :return: the level to apply to loggers. 0 <= level <=50
    :rtype: int
    """
    default = 20  # info
    level = default - (10 * verbose) + (10 * quiet)
    level = max(10, level)
    level = min(50, level)
    return level
