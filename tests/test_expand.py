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
import tempfile
import shutil

import pandas as pd
import pandas.util.testing as pdt

from Bio import BiopythonExperimentalWarning
import warnings
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import infernal
from integron_finder.utils import FastaIterator
from integron_finder.topology import Topology

_local_max_ori = infernal.local_max


def local_max_mock():
    """expand call several time local_max which run cmserach and need lot of global variable
    to avoid to call cmsearch (we don't want to test local max), I have been cached the results
    of local_max and replace the function by this mock
    """
    _cache = {('LIAN.001.C02_10', 942899, 947099, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('local_max_right_circular_1.csv')),
              ('LIAN.001.C02_10', 946899, 951099, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('local_max_right_circular_2.csv')),
              ('LIAN.001.C02_10', 930689, 934889, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('local_max_left_circular.csv')),
              ('LIAN.001.C02_10', 930689, 930889, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('local_max_left_linear.csv'))
              }

    def fake_local_max(*args, **kwargs):
        replicon, window_beg, window_end, model_attc_path = args
        fake_args = (replicon.name,
                     window_beg, window_end,
                     os.path.splitext(os.path.split(model_attc_path)[1])[0],
                     kwargs['strand_search'])
        args = tuple(fake_args)
        return _cache[args]
    return fake_local_max


class TestExpand(IntegronTest):

    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        if os.path.exists(self.tmp_dir) and os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)

        self.model_attc_path = self.find_data(os.path.join('Models', 'attc_4.cm'))
        infernal.local_max = local_max_mock()


    def tearDown(self):
        infernal.local_max = _local_max_ori
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass


    def test_expand_circular_right(self):
        circular = True
        dist_threshold = 4000
        replicon_name = 'lian.001.c02.10'
        max_attc_size = 200

        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        sequences_db = FastaIterator(replicon_path)
        topologies = Topology('lin')
        sequences_db.topologies = topologies
        replicon = sequences_db.next()

        max_elt_input = pd.read_csv(os.path.join(self._data_dir, 'max_elt_input_1.csv'))
        df_max_input = pd.read_csv(os.path.join(self._data_dir, 'df_max_input_1.csv'))
        max_elt_expected = pd.read_csv(os.path.join(self._data_dir, 'max_elt_output_lian_right.csv'))
        max_eat_received = infernal.expand(replicon,
                                           934689, 943099, max_elt_input, df_max_input,
                                           circular, dist_threshold, max_attc_size,
                                           self.model_attc_path,
                                           search_left=False, search_right=True)
        pdt.assert_frame_equal(max_elt_expected, max_eat_received)


    def test_expand_circular_left(self):
        circular = True
        dist_threshold = 4000
        replicon_name = 'lian.001.c02.10'
        max_attc_size = 200

        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        sequences_db = FastaIterator(replicon_path)
        topologies = Topology('lin')
        sequences_db.topologies = topologies
        replicon = sequences_db.next()

        max_elt_input = pd.read_csv(os.path.join(self._data_dir, 'max_elt_input_1.csv'))
        df_max_input = pd.read_csv(os.path.join(self._data_dir, 'df_max_input_1.csv'))
        max_elt_expected = pd.read_csv(os.path.join(self._data_dir, 'max_elt_output_lian_left.csv'))
        max_eat_received = infernal.expand(replicon,
                                           934689, 943099, max_elt_input, df_max_input,
                                           circular, dist_threshold, max_attc_size,
                                           self.model_attc_path,
                                           search_left=True, search_right=False)
        pdt.assert_frame_equal(max_elt_expected, max_eat_received)


    def test_expand_linear_right(self):
        circular = False
        dist_threshold = 4000
        replicon_name = 'lian.001.c02.10'
        max_attc_size = 200

        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        sequences_db = FastaIterator(replicon_path)
        topologies = Topology('lin')
        sequences_db.topologies = topologies
        replicon = sequences_db.next()

        max_elt_input = pd.read_csv(os.path.join(self._data_dir, 'max_elt_input_1.csv'))
        df_max_input = pd.read_csv(os.path.join(self._data_dir, 'df_max_input_1.csv'))
        max_elt_expected = pd.read_csv(os.path.join(self._data_dir, 'max_elt_output_lian_right.csv'))
        max_eat_received = infernal.expand(replicon,
                                           934689, 943099, max_elt_input, df_max_input,
                                           circular, dist_threshold, max_attc_size,
                                           self.model_attc_path,
                                           search_left=False, search_right=True)
        pdt.assert_frame_equal(max_elt_expected, max_eat_received)


    def test_expand_linear_left(self):
        circular = False
        dist_threshold = 4000
        replicon_name = 'lian.001.c02.10'
        max_attc_size = 200

        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        sequences_db = FastaIterator(replicon_path)
        topologies = Topology('lin')
        sequences_db.topologies = topologies
        replicon = sequences_db.next()

        max_elt_input = pd.read_csv(os.path.join(self._data_dir, 'max_elt_input_1.csv'))
        df_max_input = pd.read_csv(os.path.join(self._data_dir, 'df_max_input_1.csv'))
        max_elt_expected = pd.read_csv(os.path.join(self._data_dir, 'max_elt_output_lian_left.csv'))
        max_eat_received = infernal.expand(replicon,
                                           934689, 943099, max_elt_input, df_max_input,
                                           circular, dist_threshold, max_attc_size,
                                           self.model_attc_path,
                                           search_left=True, search_right=False)
        pdt.assert_frame_equal(max_elt_expected, max_eat_received)
