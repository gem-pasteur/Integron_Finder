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
            raise RuntimeError("'{}' is not allowed for topology".format(s))


    def _parse(self, topology_file):
        with open(topology_file) as topo_f:
            for entry in topo_f:
                if entry.startswith('#'):
                    continue
                seq_id, topology = entry.split()
                self._topology[seq_id] = self._parse_topology(topology)


    def __getitem__(self, item):
        return self._topology[item]

