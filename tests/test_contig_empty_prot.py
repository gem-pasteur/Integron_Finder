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
import shutil
import tempfile

import pandas as pd
import pandas.testing as pdt

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.scripts.finder import main
import integron_finder.scripts.finder as finder



class TestEmptyProt(IntegronTest):

    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

        self.tmp_dir = tempfile.gettempdir()
        self.out_dir = os.path.join(self.tmp_dir, 'integron_acba_test')
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
        os.makedirs(self.out_dir)
        self.which_ori = finder.shutil.which

    def tearDown(self):
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
        finder.shutil.which = self.which_ori


    def test_contig_with_empty_prot(self):
        replicon_filename = 'empty_prot_file'
        output_filename = 'Results_Integron_Finder_{}'.format(replicon_filename)
        test_result_dir = os.path.join(self.out_dir, output_filename)
        command = "integron_finder --outdir {out_dir} {replicon}".format(out_dir=self.out_dir,
                                                                         replicon=self.find_data(
                                                                             os.path.join('Replicons',
                                                                                          replicon_filename + '.fna'
                                                                                          )
                                                                               )
                                                                         )

        with self.catch_io(out=True, err=True):
            main(command.split()[1:], loglevel='WARNING')

        output_filename = '{}.integrons'.format(replicon_filename)
        expected_result_dir = self.find_data(os.path.join('Results_Integron_Finder_{}'.format(replicon_filename)))
        expected_result_path = os.path.join(expected_result_dir, output_filename)
        test_result_path = os.path.join(test_result_dir, output_filename)
        self.assertIntegronResultEqual(expected_result_path, test_result_path)

        summary_file_name = '{}.summary'.format(replicon_filename)
        test_summary_path = os.path.join(test_result_dir, summary_file_name)
        test_summary = pd.read_csv(test_summary_path, sep='\t', comment="#")
        expected_summary_path = os.path.join(expected_result_dir, summary_file_name)
        expected_summary = pd.read_csv(expected_summary_path, sep='\t', comment="#")
        pdt.assert_frame_equal(test_summary, expected_summary)


