import unittest
from tests import IntegronTest
from integron.topology import Topology

class TestFunctions(IntegronTest):


    def test_parse_topology(self):
        topo = Topology('circ')
        for t in ('circ', 'circular', 'CIRC', 'CIRCULAR'):
            self.assertEqual(topo._parse_topology(t), 'circ')
        for t in ('lin', 'linear', 'LIN', 'LINEAR'):
            self.assertEqual(topo._parse_topology(t), 'lin')


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

