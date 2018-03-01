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
import shutil
import tempfile

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import integrase
from integron_finder.scripts.finder import main

_prodigal_call = integrase.call

class TestAcba(IntegronTest):

    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

        self.tmp_dir = tempfile.gettempdir()
        self.out_dir = os.path.join(self.tmp_dir, 'integron_acba_test')
        os.makedirs(self.out_dir)
        integrase.call = self.call_wrapper(_prodigal_call)

    def tearDown(self):
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
            pass
        integrase.call = _prodigal_call


    def test_acba_simple(self):
        output_filename = 'Results_Integron_Finder_acba.007.p01.13'
        test_result_dir = os.path.join(self.out_dir, output_filename)
        command = "integron_finder --outdir {out_dir} {replicon}".format(out_dir=self.out_dir,
                                                                         replicon=self.find_data(
                                                                             os.path.join('Replicons',
                                                                                          'acba.007.p01.13.fst')
                                                                         )
                                                                         )
        with self.catch_io(out=True, err=True):
            main(command.split()[1:])
        results_file_to_test = ('acba.007.p01.13.gbk', 'acba.007.p01.13.integrons')
        for output_filename in results_file_to_test:
            expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_acba.007.p01.13',
                                                               output_filename))
            test_result_path = os.path.join(test_result_dir, output_filename)
            self.assertFileEqual(expected_result_path, test_result_path)


    def test_acba_annot(self):
        replicon_name = 'acba.007.p01.13'
        command = "integron_finder --outdir {out_dir} --func_annot --path_func_annot {annot_bank} {replicon}".format(
            out_dir=self.out_dir,
            annot_bank=os.path.normpath(self.find_data('Functional_annotation')),
            replicon=self.find_data(os.path.join('Replicons', '{}.fst'.format(replicon_name)))
        )
        with self.catch_io(out=True, err=True):
            main(command.split()[1:])
        results_file_to_test = ('{}.gbk'.format(replicon_name), '{}.integrons'.format(replicon_name))
        result_dir = os.path.join(self.out_dir, 'Results_Integron_Finder_{}'.format(replicon_name))
        for output_filename in results_file_to_test:
            expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_{}.annot'.format(replicon_name),
                                                               output_filename))
            test_result_path = os.path.join(result_dir, output_filename)
            self.assertFileEqual(expected_result_path, test_result_path)

        output_filename = os.path.join('other', replicon_name + '_Resfams_fa_table.res')
        expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_{}.annot'.format(replicon_name),
                                                           output_filename))
        test_result_path = os.path.join(result_dir, output_filename)
        with open(expected_result_path) as expected_result_file, open(test_result_path) as test_result_file :
            for expected_line, result_line in zip(expected_result_file, test_result_file):
                if result_line.startswith('# Program:         hmmsearch'):
                    break
                self.assertEqual(expected_line, result_line)

