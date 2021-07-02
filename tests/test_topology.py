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

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.topology import Topology


class TestTopology(IntegronTest):

    def test_parse_topology(self):
        topo = Topology('circ')
        for t in ('circ', 'circular', 'CIRC', 'CIRCULAR'):
            self.assertEqual(topo._parse_topology(t), 'circ')
        for t in ('lin', 'linear', 'LIN', 'LINEAR'):
            self.assertEqual(topo._parse_topology(t), 'lin')
        with self.assertRaises(RuntimeError) as ctx:
            topo._parse_topology('foo')
        self.assertEqual(str(ctx.exception), "'foo' is not allowed for topology")


    def test_parse(self):
        topo = Topology('circ')
        topo._parse(self.find_data('topology.txt'))
        expected = {
                    'seq1': 'circ',
                    'seq2': 'circ',
                    'seq3': 'lin',
                    'seq4': 'lin',
                    'seq5': 'circ',
                    'seq6': 'circ',
                    'seq7': 'lin',
                    'seq8': 'lin',
                    }
        self.assertDictEqual(expected, topo._topology)


    def test_getitem(self):
        topo = Topology('circ', topology_file=self.find_data('topology.txt'))
        expected = {
                    'seq1': 'circ',
                    'seq2': 'circ',
                    'seq3': 'lin',
                    'seq4': 'lin',
                    'seq5': 'circ',
                    'seq6': 'circ',
                    'seq7': 'lin',
                    'seq8': 'lin',
                    }
        for seqid, topology in expected.items():
            self.assertEqual(topology, topo[seqid])
        self.assertEqual('circ', topo['foo'])

