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
import tempfile

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.topology import Topology
from integron_finder import utils


class TestUtils(IntegronTest):

    def test_read_multi_prot_fasta(self):
        replicon_id = 'ACBA.007.P01_13'
        replicon_path = self.find_data(os.path.join('Proteins', replicon_id + '.prt'))
        replicon = utils.MultiFastaReader(replicon_path)
        expected_seq_id = ['{}_{}'.format(replicon_id, i) for i in range(1, 24)]
        received_seq_id = [seq.id for seq in replicon]
        self.assertListEqual(expected_seq_id, received_seq_id)

    def test_FastaIterator_test_topologies(self):
        file_name = 'multi_fasta'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        topologies = Topology(1, 'lin')
        with utils.FastaIterator(replicon_path) as seq_db:
            seq_db.topologies = topologies

        with utils.FastaIterator(replicon_path) as seq_db:
            seq_db.topologies = topologies
            received_seq_top = [seq.topology for seq in seq_db]
        expected_seq_top = ['lin', 'lin', 'lin']
        self.assertListEqual(expected_seq_top, received_seq_top)

        # test FastaIterator no topology is provided
        with utils.FastaIterator(replicon_path) as seq_db:
            received_seq_top = [seq.topology for seq in seq_db]
        expected_seq_top = [None, None, None]
        self.assertListEqual(expected_seq_top, received_seq_top)

        topologies_data = {'ACBA.007.P01_13': 'lin',
                           'LIAN.001.C02_10': 'circ',
                           'PSSU.001.C01_13': 'lin',
                           }
        with tempfile.NamedTemporaryFile(mode='w') as topology_file:
            for rep, topo in topologies_data.items():
                topology_file.write("{} {}\n".format(rep, topo))
            topology_file.flush()
            topologies = Topology(1, 'lin', topology_file=topology_file.name)
            with utils.FastaIterator(replicon_path) as seq_db:
                seq_db.topologies = topologies
                received_seq_top = {seq.id: seq.topology for seq in seq_db}
            self.assertDictEqual(topologies_data, received_seq_top)

        file_name = 'acba_short'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        topologies = Topology(1, 'circ')
        with utils.FastaIterator(replicon_path) as seq_db:
            seq_db.topologies = topologies
            received_seq_top = [seq.topology for seq in seq_db]
        expected_seq_top = ['lin']
        self.assertListEqual(expected_seq_top, received_seq_top)


    def test_FastaIterator(self):
        file_name = 'multi_fasta'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        topologies = Topology(1, 'lin')
        with utils.FastaIterator(replicon_path) as seq_db:
            seq_db.topologies = topologies
            received_seq_id = sorted([seq.id for seq in seq_db])

        expected_seq_id = sorted(['ACBA.007.P01_13', 'LIAN.001.C02_10', 'PSSU.001.C01_13'])
        self.assertListEqual(expected_seq_id, received_seq_id)
        self.assertEqual(len(seq_db), 3)

        expected_seq_name = expected_seq_id
        with utils.FastaIterator(replicon_path) as seq_db:
            seq_db.topologies = topologies
            received_seq_name = sorted([seq.name for seq in seq_db])
        self.assertListEqual(expected_seq_name, received_seq_name)

        replicon_name = 'foo'
        with utils.FastaIterator(replicon_path, replicon_name=replicon_name) as seq_db:
            seq_db.topologies = topologies
            received_seq_name_id = sorted([(seq.name, seq.id) for seq in seq_db])

        expected_seq_name_id = sorted(
            [(replicon_name, _id) for _id in ['ACBA.007.P01_13', 'LIAN.001.C02_10', 'PSSU.001.C01_13']]
        )
        self.assertEqual(expected_seq_name_id, received_seq_name_id)

        file_name = 'replicon_ambiguous_char'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        with utils.FastaIterator(replicon_path) as seq_db:
            received_seq_id = sorted([seq.id for seq in seq_db if seq])
        expected_seq_id = sorted(['seq_1', 'seq_2', 'seq_3', 'seq_4'])
        self.assertListEqual(expected_seq_id, received_seq_id)

        file_name = 'replicon_bad_char'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        expected_warning = """sequence seq_(3|4) contains invalid characters, the sequence is skipped.
sequence seq_(3|4) contains invalid characters, the sequence is skipped."""
        with utils.FastaIterator(replicon_path) as seq_db:
            # 2 sequences are rejected so 2 message is produced (for seq 3 and seq 4)
            with self.catch_log() as log:
                received_seq_id = sorted([seq.id for seq in seq_db if seq])
                got_warning = log.get_value().strip()
        self.assertRegex(got_warning, expected_warning)
        expected_seq_id = sorted(['seq_1', 'seq_2'])
        self.assertListEqual(expected_seq_id, received_seq_id)

        file_name = 'replicon_too_short'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        expected_warning = r"""sequence seq_(4|2) is too short \(32 bp\), the sequence is skipped \(must be > 50bp\).
sequence seq_(4|2) is too short \(32 bp\), the sequence is skipped \(must be > 50bp\)."""
        with utils.FastaIterator(replicon_path) as seq_db:
            # 2 sequences are rejected so 2 messages are produced (for seq 2 & 4)
            with self.catch_log() as log:
                received_seq_id = sorted([seq.id for seq in seq_db if seq])
                got_warning = log.get_value().strip()

        self.assertRegex(got_warning, expected_warning)
        expected_seq_id = sorted(['seq_1', 'seq_3'])
        self.assertListEqual(expected_seq_id, received_seq_id)


    def test_model_len(self):
        model_path = self.find_data(os.path.join('Models', 'attc_4.cm'))
        self.assertEqual(utils.model_len(model_path), 47)
        bad_path = 'nimportnaoik'
        with self.assertRaises(IOError) as ctx:
            with self.catch_log():
                utils.model_len(bad_path)
        self.assertEqual(str(ctx.exception),
                         "Path to model_attc '{}' does not exists".format(bad_path))
        bad_path = self.find_data(os.path.join('Models', 'phage-int.hmm'))
        with self.assertRaises(RuntimeError) as ctx:
            with self.catch_log():
                utils.model_len(bad_path)
        self.assertEqual(str(ctx.exception),
                         "CLEN not found in '{}', maybe it's not infernal model file".format(bad_path))


    def test_get_name_from_path(self):
        self.assertEqual(utils.get_name_from_path('/foo/bar.baz'), 'bar')
        self.assertEqual(utils.get_name_from_path('bar.baz'), 'bar')
        self.assertEqual(utils.get_name_from_path('../foo/bar.baz'), 'bar')
        self.assertEqual(utils.get_name_from_path('../foo/bar'), 'bar')


    def test_log_level(self):
        for verbose, quiet, expected_level in [(0, 0, 20), (0, 2, 40), (0, 5, 50), (1, 0, 10), (3, 0, 10), (2, 2, 20)]:
            self.assertEqual(utils.log_level(verbose, quiet), expected_level)
