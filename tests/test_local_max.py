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

import os
import tempfile
import shutil
import re

import pandas as pd
import pandas.testing as pdt

# # display warning only for non installed integron_finder
# from Bio import BiopythonExperimentalWarning
# import warnings
# warnings.simplefilter('ignore', BiopythonExperimentalWarning)

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.utils import FastaIterator
from integron_finder.topology import Topology
from integron_finder import infernal

_call_ori = infernal.call
_read_infernal_ori = infernal.read_infernal

from tests import which


def call_wrapper():
    """
    hmmsearch or prodigal write lot of things on stderr or stdout 
    which noise the unit test output
    So I replace the `call` function in module integron_finder
    by a wrapper which call the original function but add redirect stderr and stdout
    in dev_null
    :return: wrapper around integron_finder.call
    :rtype: function
    """
    def wrapper(*args, **kwargs):
        with open(os.devnull, 'w') as f:
            kwargs['stderr'] = f
            kwargs['stdout'] = f
            res = _call_ori(*args, **kwargs)
        return res
    return wrapper


def read_infernal_mock(tmp_dir):
    """expand call several time local_max which run cmserach and need lot of global variable
    to avoid to call cmsearch (we don't want to test local max), I have been ccahed the results
    of local_max and replace the function by this mock
    """
    _cache = {
        (os.path.join(tmp_dir, 'LIAN.001.C02_10_942899_947099_subseq_attc_table.res'),
         'LIAN.001.C02_10', 47, 1.0, 200, 40):
            pd.DataFrame([['LIAN.001.C02_10', 'attC_4', 1, 47, 371, 496, '+', 0.130000],
                          ['LIAN.001.C02_10', 'attC_4', 1, 47, 1109, 1234, '+', 0.049000],
                          ['LIAN.001.C02_10', 'attC_4', 1, 47, 1573, 1699, '+', 0.000005]],
                         columns=['Accession_number', 'cm_attC', 'cm_debut',
                                  'cm_fin', 'pos_beg', 'pos_end', 'sens', 'evalue']),
        (os.path.join(tmp_dir, 'LIAN.001.C02_10_946899_951099_subseq_attc_table.res'),
         'LIAN.001.C02_10', 47, 1.0, 200, 40):
            pd.DataFrame(columns=['Accession_number', 'cm_attC', 'cm_debut',
                         'cm_fin', 'pos_beg', 'pos_end', 'sens', 'evalue']),
        (os.path.join(tmp_dir, 'LIAN.001.C02_10_930689_934889_subseq_attc_table.res'),
         'LIAN.001.C02_10', 47, 1.0, 200, 40):
            pd.DataFrame(columns=['Accession_number', 'cm_attC', 'cm_debut',
                         'cm_fin', 'pos_beg', 'pos_end', 'sens', 'evalue']),
              }

    def fake_read_infernal(tblout_path, replicon_name, len_model_attc,
                           evalue=None, size_max_attc=None, size_min_attc=None):
        args = (tblout_path, replicon_name, len_model_attc, evalue, size_max_attc, size_min_attc)
        return _cache[args]
    return fake_read_infernal


class TestLocalMax(IntegronTest):

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

        self.cmsearch = which('cmsearch')
        self.out_dir = self.tmp_dir
        self.model_attc_path = self.find_data(os.path.join('Models', 'attc_4.cm'))
        self.cpu = 1
        replicon_name = 'lian.001.c02.10'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        topologies = Topology('lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            self.replicon = next(sequences_db)
        self.evalue_attc = 1.
        self.max_attc_size = 200
        self.min_attc_size = 40
        self.length_cm = 47  # length in 'CLEN' (value for model attc_4.cm)
        self.call = call_wrapper()
        infernal.read_infernal = read_infernal_mock(self.tmp_dir)

    def tearDown(self):
        infernal.call = _call_ori
        infernal.read_infernal = _read_infernal_ori
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass


    def test_local_max_top(self):
        win_beg = 942899
        win_end = 947099
        strand_search = 'top'
        local_max_received = infernal.local_max(self.replicon,
                                                win_beg, win_end,
                                                self.model_attc_path,
                                                strand_search=strand_search,
                                                evalue_attc=self.evalue_attc,
                                                max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                                cmsearch_bin=self.cmsearch, out_dir=self.out_dir, cpu=self.cpu
                                                )
        local_max_expected = pd.DataFrame([['LIAN.001.C02_10', 'attC_4', 1, 47, 943270, 943395, '+', 0.13],
                                           ['LIAN.001.C02_10', 'attC_4', 1, 47, 944008, 944133, '+', 0.049],
                                           ['LIAN.001.C02_10', 'attC_4', 1, 47, 944472, 944598, '+', 5e-06]],
                                          columns=['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg',
                                                   'pos_end', 'sens', 'evalue'])
        pdt.assert_frame_equal(local_max_expected, local_max_received)

    def test_no_local_max_top(self):
        win_beg = 946899
        win_end = 951099
        strand_search = 'top'
        local_max_received = infernal.local_max(self.replicon,
                                                win_beg, win_end,
                                                self.model_attc_path,
                                                strand_search=strand_search,
                                                evalue_attc=self.evalue_attc,
                                                max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                                cmsearch_bin=self.cmsearch, out_dir=self.out_dir, cpu=self.cpu
                                                )
        local_max_expected = pd.DataFrame(columns=['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg',
                                                   'pos_end', 'sens', 'evalue'])
        pdt.assert_frame_equal(local_max_expected, local_max_received)


    def test_no_local_max_bottom(self):
        win_beg = 930689
        win_end = 934889
        strand_search = 'bottom'
        local_max_received = infernal.local_max(self.replicon,
                                                win_beg, win_end,
                                                self.model_attc_path,
                                                strand_search=strand_search,
                                                evalue_attc=self.evalue_attc,
                                                max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                                cmsearch_bin=self.cmsearch, out_dir=self.out_dir, cpu=self.cpu
                                                )
        local_max_expected = pd.DataFrame(columns=['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg',
                                                   'pos_end', 'sens', 'evalue'])
        pdt.assert_frame_equal(local_max_expected, local_max_received)


    def test_local_max_both(self):
        win_beg = 942899
        win_end = 947099
        strand_search = 'both'

        local_max_received = infernal.local_max(self.replicon,
                                                win_beg, win_end,
                                                self.model_attc_path,
                                                strand_search=strand_search,
                                                evalue_attc=self.evalue_attc,
                                                max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                                cmsearch_bin=self.cmsearch, out_dir=self.out_dir, cpu=self.cpu
                                                )
        local_max_expected = pd.DataFrame([['LIAN.001.C02_10', 'attC_4', 1, 47, 943270, 943395, '+', 0.13],
                                           ['LIAN.001.C02_10', 'attC_4', 1, 47, 944008, 944133, '+', 0.049],
                                           ['LIAN.001.C02_10', 'attC_4', 1, 47, 944472, 944598, '+', 5e-06]],
                                          columns=['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg',
                                                   'pos_end', 'sens', 'evalue'])
        pdt.assert_frame_equal(local_max_expected, local_max_received)


    def local_max_end_beg(self):
        win_beg = 942899
        win_end = 947099
        strand_search = 'top'
        local_max_received = infernal.local_max(self.replicon,
                                                win_beg, win_end,
                                                self.model_attc_path,
                                                strand_search=strand_search,
                                                evalue_attc=self.evalue_attc,
                                                max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                                cmsearch_bin=self.cmsearch, out_dir=self.out_dir, cpu=self.cpu
                                                )
        local_max_expected = pd.DataFrame([['LIAN.001.C02_10', 'attC_4', 1, 47, 943270, 943395, '+', 0.13],
                                           ['LIAN.001.C02_10', 'attC_4', 1, 47, 944008, 944133, '+', 0.049],
                                           ['LIAN.001.C02_10', 'attC_4', 1, 47, 944472, 944598, '+', 5e-06]],
                                          columns=['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg',
                                                   'pos_end', 'sens', 'evalue'])
        pdt.assert_frame_equal(local_max_expected, local_max_received)


    def test_local_max_cmsearch_failed(self):
        win_beg = 942899
        win_end = 947099
        strand_search = 'top'

        cmsearch_bin = 'failed_cmsearch'
        with self.assertRaises(RuntimeError) as ctx:
            _ = infernal.local_max(self.replicon,
                                   win_beg, win_end,
                                   self.model_attc_path,
                                   strand_search=strand_search,
                                   evalue_attc=self.evalue_attc,
                                   max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                   cmsearch_bin=cmsearch_bin, out_dir=self.out_dir, cpu=self.cpu
                                   )
        self.assertTrue(re.search("failed : \[Errno 2\] No such file or directory: '{}'".format(cmsearch_bin),
                        str(ctx.exception)))


    def test_local_max_cmsearch_bad_rc(self):
        win_beg = 942899
        win_end = 947099
        strand_search = 'top'
        infernal.call = lambda x, stdout=None: 1
        with self.assertRaises(RuntimeError) as ctx:
            _ = infernal.local_max(self.replicon,
                                   win_beg, win_end,
                                   self.model_attc_path,
                                   strand_search=strand_search,
                                   evalue_attc=self.evalue_attc,
                                   max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                   out_dir=self.out_dir, cpu=self.cpu
                                   )
        self.assertTrue(str(ctx.exception).endswith("failed returncode = {}".format(infernal.call(None))))
