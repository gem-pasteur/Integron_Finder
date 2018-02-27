# -*- coding: utf-8 -*-

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

