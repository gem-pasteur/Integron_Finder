from collections import defaultdict

class Topology(object):


    def __init__(self, default, topology_file=None):
        self._default = self._parse_topology(default)
        self._topology = defaultdict(lambda: self._default)
        if topology_file:
            self._parse(topology_file)


    def _parse_topology(self, s):
        s = s.lower()
        if s in 'circular':
            return 'circ'
        elif s in 'linear':
            return 'lin'
        else:
            raise RuntimeError("'' is not allowed for topology".format(s))


    def _parse(self, topology_file):
        with open(topology_file) as topo_f:
            for entry in topo_f:
                if entry.startswith('#'):
                    continue
                seq_id, topology = entry.split()
                self._topology[seq_id] = self._parse_topology(topology)


    def __getitem__(self, item):
        return self._topology[item]

