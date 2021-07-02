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

from collections import defaultdict


class Topology:
    """Class to parse and handle replicons topologies"""

    def __init__(self, default, topology_file=None):
        """

        :param str default: the default topology
        :param topology_file: the path to the file where topology for replicon are specified
        """
        self._default = self._parse_topology(default)
        self._topology = defaultdict(lambda: self._default)
        if topology_file:
            self._parse(topology_file)


    def _parse_topology(self, topo):
        """
        Parse a field topology in topology file
        the authorized values are circular, linear or circ, lin, or in uppercase

        :param topo: the field corresponding to topology in topology file
        :return: the topology in "normed" format 'circ' or 'lin'
        :rtype: str
        """
        topo = topo.lower()
        if topo in 'circular':
            return 'circ'
        elif topo in 'linear':
            return 'lin'
        else:
            raise RuntimeError("'{}' is not allowed for topology".format(topo))


    def _parse(self, topology_file):
        """
        Parse a topology file where topology is specified for replicons
        on each line a topology is specified for a replicon
        the syntax of each line is ::

            replicon_id topology

        the allowed value for toplogy are  'circ', 'circular', 'lin', 'linear'

        :param str topology_file: The path to the topology file
        """
        with open(topology_file) as topo_f:
            for entry in topo_f:
                if entry.startswith('#'):
                    continue
                seq_id, topology = entry.split()
                self._topology[seq_id] = self._parse_topology(topology)


    def __getitem__(self, replicon_id):
        """
        :param str replicon_id: The id of the replicon.
        :returns: the topology for the replicon corresponding to the replicon_id
        """
        return self._topology[replicon_id]
