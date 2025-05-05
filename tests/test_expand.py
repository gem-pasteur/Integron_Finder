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
import shutil

import pandas as pd
import pandas.testing as pdt

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
    _cache = {
              # test_expand_circular_no_more_attc
              ('LIAN.001.C02_10', 954528, 958728, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_no_more_attc_954528_958728_top.csv')),
              ('LIAN.001.C02_10', 958528, 962728, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_no_more_attc_958528_962728_top.csv')),
              ('LIAN.001.C02_10', 962528, 966728, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_no_more_attc_962528_966728_top.csv')),
              ('LIAN.001.C02_10', 966528, 970728, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_no_more_attc_966528_970728_top.csv')),
              ('LIAN.001.C02_10', 970528, 974728, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_no_more_attc_970528_974728_top.csv')),
              ('LIAN.001.C02_10', 974528, 978728, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_no_more_attc_974528_978728_top.csv')),
              ('LIAN.001.C02_10', 978528, 982728, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_no_more_attc_978528_982728_top.csv')),
              ('LIAN.001.C02_10', 982528, 986728, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_no_more_attc_982528_986728_top.csv')),

              # test_expand_circular_right_window_in_replicon
              ('NZ_CP016323.1', 23108, 27308, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_23108_27308_top.csv')),
              ('NZ_CP016323.1', 27108, 31308, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_27108_31308_top.csv')),
              ('NZ_CP016323.1', 31108, 35308, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_31108_35308_top.csv')),
              ('NZ_CP016323.1', 35108, 39308, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_35108_39308_top.csv')),
              ('NZ_CP016323.1', 39108, 2458, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_39108_2458_top.csv')),
              ('NZ_CP016323.1', 2258, 6458, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_2258_6458_top.csv')),
              ('NZ_CP016323.1', 6258, 10458, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_6258_10458_top.csv')),
              ('NZ_CP016323.1', 10258, 14458, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_10258_14458_top.csv')),
              ('NZ_CP016323.1', 14258, 15466, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_14258_15466_top.csv')),

              # test_expand_circular_left_window_in_replicon
              ('NZ_CP016323.1', 11266, 15466, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_in_rep_11266_15466_bottom.csv')),
              ('NZ_CP016323.1', 7266, 11466, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_in_rep_7266_11466_bottom.csv')),
              ('NZ_CP016323.1', 3266, 7466, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_in_rep_3266_7466_bottom.csv')),
              ('NZ_CP016323.1', 40116, 3466, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_in_rep_40116_3466_bottom.csv')),
              ('NZ_CP016323.1', 36116, 40316, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_in_rep_36116_40316_bottom.csv')),
              ('NZ_CP016323.1', 32116, 36316, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_in_rep_32116_36316_bottom.csv')),
              ('NZ_CP016323.1', 28116, 32316, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_in_rep_28116_32316_bottom.csv')),
              ('NZ_CP016323.1', 24116, 28316, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_in_rep_24116_28316_bottom.csv')),
              ('NZ_CP016323.1', 23108, 24316, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_in_rep_23108_24316_bottom.csv')),

              # test_expand_circular_left_window_over_ori
              ('NZ_CP016323.1', 35181, 39381, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_over_ori_35181_39381_bottom.csv')),
              ('NZ_CP016323.1', 31181, 35381, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_over_ori_31181_35381_bottom.csv')),
              ('NZ_CP016323.1', 27181, 31381, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_over_ori_27181_31381_bottom.csv')),
              ('NZ_CP016323.1', 23181, 27381, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_over_ori_23181_27381_bottom.csv')),
              ('NZ_CP016323.1', 19181, 23381, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_over_ori_19181_23381_bottom.csv')),
              ('NZ_CP016323.1', 15181, 19381, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_over_ori_15181_19381_bottom.csv')),
              ('NZ_CP016323.1', 11181, 15381, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_over_ori_11181_15381_bottom.csv')),
              ('NZ_CP016323.1', 7181, 11381, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_over_ori_7181_11381_bottom.csv')),
              ('NZ_CP016323.1', 3271, 7381, 'attc_4', 'bottom'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_left_win_over_ori_3271_7381_bottom.csv')),

              # test_expand_circular_right_window_over_ori
              ('NZ_CP016323.1', 3271, 7471, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_over_ori_3271_7471_top.csv')),
              ('NZ_CP016323.1', 7271, 11471, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_over_ori_7271_11471_top.csv')),
              ('NZ_CP016323.1', 11271, 15471, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_over_ori_11271_15471_top.csv')),
              ('NZ_CP016323.1', 15271, 19471, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_over_ori_15271_19471_top.csv')),
              ('NZ_CP016323.1', 19271, 23471, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_over_ori_19271_23471_top.csv')),
              ('NZ_CP016323.1', 23271, 27471, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_over_ori_23271_27471_top.csv')),
              ('NZ_CP016323.1', 27271, 31471, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_over_ori_27271_31471_top.csv')),
              ('NZ_CP016323.1', 31271, 35471, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_over_ori_31271_35471_top.csv')),
              ('NZ_CP016323.1', 35271, 39381, 'attc_4', 'top'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_over_ori_35271_39381_top.csv')),

              # test_expand_circular_both_window_in_rep
              # same as local_max_linear_23108_27308_both.csv
              # ('NZ_CP016323.1', 23108, 27308, 'attc_4', 'both'):
              #     pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_23108_27308_both.csv')),
              # same as local_max_linear_27108_31308_both.csv
              # ('NZ_CP016323.1', 27108, 31308, 'attc_4', 'both'):
              #     pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_27108_31308_both.csv')),
              # same as local_max_linear_31108_35308_both.csv
              # ('NZ_CP016323.1', 31108, 35308, 'attc_4', 'both'):
              #     pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_31108_35308_both.csv')),
              # same as local_max_linear_31108_35308_both.csv
              #('NZ_CP016323.1', 35108, 39308, 'attc_4', 'both'):
              #    pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_35108_39308_both.csv')),
              ('NZ_CP016323.1', 39108, 2458, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_39108_2458_both.csv')),
              ('NZ_CP016323.1', 2258, 6458, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_2258_6458_both.csv')),
              ('NZ_CP016323.1', 6258, 10458, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_6258_10458_both.csv')),
              ('NZ_CP016323.1', 10258, 14458, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_10258_14458_both.csv')),
              ('NZ_CP016323.1', 14258, 15466, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_circ_right_win_in_rep_14258_15466_both.csv')),

              # test_expand_linear
              ('NZ_CP016323.1', 23108, 27308, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_linear_23108_27308_both.csv')),
              ('NZ_CP016323.1', 27108, 31308, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_linear_27108_31308_both.csv')),
              ('NZ_CP016323.1', 31108, 35308, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_linear_31108_35308_both.csv')),
              ('NZ_CP016323.1', 35108, 39308, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_linear_35108_39308_both.csv')),
              ('NZ_CP016323.1', 39108, 40850, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_linear_39108_40850_both.csv')),
              ('NZ_CP016323.1', 11266, 15466, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_linear_11266_15466_both.csv')),
              ('NZ_CP016323.1', 7266, 11466, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_linear_7266_11466_both.csv')),
              ('NZ_CP016323.1', 3266, 7466, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_linear_3266_7466_both.csv')),
              ('NZ_CP016323.1', 0, 3466, 'attc_4', 'both'):
                  pd.read_csv(TestExpand.find_data('expand_mock', 'local_max_linear_0_3466_both.csv')),

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

        self._tmp_dir = tempfile.TemporaryDirectory(prefix='tmp_test_integron_finder')
        self.tmp_dir = self._tmp_dir.name
        if os.path.exists(self.tmp_dir) and os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)

        self.model_attc_path = self.find_data('Models', 'attc_4.cm')
        infernal.local_max = local_max_mock()


    def tearDown(self):
        infernal.local_max = _local_max_ori
        self._tmp_dir.cleanup()


    def test_expand_circular_no_more_attc(self):
        circular = True
        dist_threshold = 4000
        replicon_name = 'lian.001.c02.10'
        max_attc_size = 200
        min_attc_size = 15
        evalue_attc = 10
        wb = 941053
        we = 954728
        search_left = False
        search_right = True

        replicon_path = self.find_data('Replicons', replicon_name + '.fst')
        topologies = Topology(1, 'circ')

        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        max_elt_input = pd.read_csv(self.find_data('expand_mock', 'expand_circ_no_more_attc_input.csv'))
        max_elt_expected = pd.read_csv(self.find_data('expand_mock', 'expand_circ_no_more_attc_output.csv'))

        max_elt_received = infernal.expand(replicon,
                                           wb, we,
                                           max_elt_input,
                                           circular,
                                           dist_threshold,
                                           self.model_attc_path,
                                           max_attc_size, min_attc_size,
                                           evalue_attc=evalue_attc,
                                           search_left=search_left,
                                           search_right=search_right
                                           )

        data_type = {'Accession_number': 'str',
                     'cm_attC': 'str',
                     'cm_debut': 'int',
                     'cm_fin': 'int',
                     'pos_beg': 'int',
                     'pos_end': 'int',
                     'sens': 'str',
                     'evalue': 'float'
                    }
        max_elt_received = max_elt_received.astype(dtype=data_type)
        pdt.assert_frame_equal(max_elt_expected, max_elt_received)


    def test_expand_circular_right_window_in_replicon(self):
        # NZ_CP016323 contains 2 integrons at first pass
        # one over the origin 15266, 23308
        # in this test with this integrons and search on the right of this integron
        # left search is skip as right allow to search entire replicon
        # NZ_CP016323 contains only CALIN spread over all replicon
        # test that expand stop when a complete turn is done.
        circular = True
        dist_threshold = 4000
        replicon_name = 'NZ_CP016323'
        max_attc_size = 200
        min_attc_size = 15
        evalue_attc = 10
        wb = 15266
        we = 23308
        search_left = False
        search_right = True

        replicon_path = self.find_data('Replicons', replicon_name + '.fna')
        topologies = Topology(1, 'circ')

        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        max_elt_input = pd.read_csv(self.find_data('expand_mock', 'expand_circ_right_in_rep_input.csv'))
        max_elt_expected = pd.read_csv(self.find_data('expand_mock', 'expand_circ_right_in_rep_output.csv'))

        max_elt_received = infernal.expand(replicon,
                                           wb, we,
                                           max_elt_input,
                                           circular,
                                           dist_threshold,
                                           self.model_attc_path,
                                           max_attc_size, min_attc_size,
                                           evalue_attc=evalue_attc,
                                           search_left=search_left,
                                           search_right=search_right)
        pdt.assert_frame_equal(max_elt_expected, max_elt_received)


    def test_expand_circular_left_window_in_replicon(self):
        # NZ_CP016323 contains 2 integrons at first pass
        # one over the origin the other in middle of replicon 15266, 23308
        # in this test with this integrons and search on the left of this integron
        # left search is skip as right allow to search entire replicon
        # NZ_CP016323 contains only CALIN spread over all replicon
        # test that expand stop when a complete turn is done.
        circular = True
        dist_threshold = 4000
        replicon_name = 'NZ_CP016323'
        max_attc_size = 200
        min_attc_size = 15
        evalue_attc = 10
        wb = 15266
        we = 23308
        search_left = True
        search_right = False

        replicon_path = self.find_data('Replicons', replicon_name + '.fna')
        topologies = Topology(1, 'circ')

        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        max_elt_input = pd.read_csv(self.find_data('expand_mock', 'expand_circ_left_in_rep_input.csv'))
        max_elt_expected = pd.read_csv(self.find_data('expand_mock', 'expand_circ_left_in_rep_output.csv'))

        max_elt_received = infernal.expand(replicon,
                                           wb, we,
                                           max_elt_input,
                                           circular,
                                           dist_threshold,
                                           self.model_attc_path,
                                           max_attc_size, min_attc_size,
                                           evalue_attc=evalue_attc,
                                           search_left=search_left,
                                           search_right=search_right)
        pdt.assert_frame_equal(max_elt_expected, max_elt_received)


    def test_expand_circular_left_window_over_ori(self):
        # NZ_CP016323 contains 2 integrons at first pass
        # in this test we use one over the origin 39181, 3471 and search on the left of this integron
        # NZ_CP016323 contains only CALIN spread over all replicon
        # test that expand stop when a complete turn is done.
        circular = True
        dist_threshold = 4000
        replicon_name = 'NZ_CP016323'
        max_attc_size = 200
        min_attc_size = 15
        evalue_attc = 10
        wb = 39181
        we = 3471
        search_left = True
        search_right = False

        replicon_path = self.find_data('Replicons', replicon_name + '.fna')
        topologies = Topology(1, 'circ')

        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        max_elt_input = pd.read_csv(self.find_data('expand_mock', 'expand_circ_left_over_ori_input.csv'))
        max_elt_expected = pd.read_csv(self.find_data('expand_mock', 'expand_circ_left_over_ori_output.csv'))

        max_elt_received = infernal.expand(replicon,
                                           wb, we,
                                           max_elt_input,
                                           circular,
                                           dist_threshold,
                                           self.model_attc_path,
                                           max_attc_size, min_attc_size,
                                           evalue_attc=evalue_attc,
                                           search_left=search_left,
                                           search_right=search_right)
        pdt.assert_frame_equal(max_elt_expected, max_elt_received)


    def test_expand_circular_right_window_over_ori(self):
        # NZ_CP016323 contains 2 integrons at first pass
        # one over the origin 39181, 3471
        # in this test with we use one over the origin 39181, 3471  and search on the right of this integron
        # NZ_CP016323 contains only CALIN spread over all replicon
        # test that expand stop when a complete turn is done.
        circular = True
        dist_threshold = 4000
        replicon_name = 'NZ_CP016323'
        max_attc_size = 200
        min_attc_size = 15
        evalue_attc = 10
        wb = 39181
        we = 3471
        search_left = False
        search_right = True

        replicon_path = self.find_data('Replicons', replicon_name + '.fna')
        topologies = Topology(1, 'circ')

        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        max_elt_input = pd.read_csv(self.find_data('expand_mock', 'expand_circ_right_over_ori_input.csv'))
        max_elt_expected = pd.read_csv(self.find_data('expand_mock', 'expand_circ_right_over_ori_output.csv'))

        max_elt_received = infernal.expand(replicon,
                                           wb, we,
                                           max_elt_input,
                                           circular,
                                           dist_threshold,
                                           self.model_attc_path,
                                           max_attc_size, min_attc_size,
                                           evalue_attc=evalue_attc,
                                           search_left=search_left,
                                           search_right=search_right)
        pdt.assert_frame_equal(max_elt_expected, max_elt_received)


    def test_expand_circular_both_window_in_rep(self):
        # NZ_CP016323 contains 2 integrons at first pass
        # one over the origin the other in middle of replicon 15266, 23308
        # in this test we search from the 2nd integron and search on the right + left
        # left search is skip as right allow to search entire replicon
        # NZ_CP016323 contains only CALIN spread over all replicon
        # test that expand stop when a complete turn is done.
        circular = True
        dist_threshold = 4000
        replicon_name = 'NZ_CP016323'
        max_attc_size = 200
        min_attc_size = 15
        evalue_attc = 10
        wb = 15266
        we = 23308
        search_left = True
        search_right = True

        replicon_path = self.find_data('Replicons', replicon_name + '.fna')
        topologies = Topology(1, 'circ')

        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        max_elt_input = pd.read_csv(self.find_data('expand_mock', 'expand_circ_both_win_in_rep_input.csv'))
        max_elt_expected = pd.read_csv(self.find_data('expand_mock', 'expand_circ_both_win_in_rep_output.csv'))

        max_elt_received = infernal.expand(replicon,
                                           wb, we,
                                           max_elt_input,
                                           circular,
                                           dist_threshold,
                                           self.model_attc_path,
                                           max_attc_size, min_attc_size,
                                           evalue_attc=evalue_attc,
                                           search_left=search_left,
                                           search_right=search_right)
        pdt.assert_frame_equal(max_elt_expected, max_elt_received)


    def test_expand_linear(self):
        # NZ_CP016323 contains 2 integrons at first pass
        # one over the origin 39181, 3471 the other between 15266, 23308
        # in this test we loose the first one as the topology is set to linear
        # NZ_CP016323 contains only CALIN spread over all replicon
        # test that expand stop when reach the end or the origin of replicon.
        circular = False
        dist_threshold = 4000
        replicon_name = 'NZ_CP016323'
        max_attc_size = 200
        min_attc_size = 15
        evalue_attc = 10
        wb = 15266
        we = 23308
        search_left = True
        search_right = True

        replicon_path = self.find_data('Replicons', replicon_name + '.fna')
        topologies = Topology(1, 'lin')

        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        max_elt_input = pd.read_csv(self.find_data('expand_mock', 'expand_linear_input.csv'))
        max_elt_expected = pd.read_csv(self.find_data('expand_mock', 'expand_linear_output.csv'))

        max_elt_received = infernal.expand(replicon,
                                           wb, we,
                                           max_elt_input,
                                           circular,
                                           dist_threshold,
                                           self.model_attc_path,
                                           max_attc_size, min_attc_size,
                                           evalue_attc=evalue_attc,
                                           search_left=search_left,
                                           search_right=search_right)
        pdt.assert_frame_equal(max_elt_expected, max_elt_received)
