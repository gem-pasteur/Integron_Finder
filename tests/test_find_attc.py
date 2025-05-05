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
import re

# # display warning only for non installed integron_finder
# from Bio import BiopythonExperimentalWarning
# import warnings
# warnings.simplefilter('ignore', BiopythonExperimentalWarning)

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import infernal

_run_ori = infernal.subprocess.run


class TestFindAttc(IntegronTest):

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
        self.cmsearch_path = shutil.which("cmsearch")
        self.cpu = 1
        self.model_attc = self.find_data(os.path.join('Models', 'attc_4.cm'))
        self.replicon_name = 'acba.007.p01.13'
        self.replicon_path = self.find_data(os.path.join('Replicons', self.replicon_name + '.fst'))
        infernal.subprocess.run = self.mute_call(_run_ori)


    def tearDown(self):
        self._tmp_dir.cleanup()
        infernal.subprocess.run = _run_ori


    def test_find_attc(self):
        infernal.find_attc(self.replicon_path, self.replicon_name, self.cmsearch_path, self.tmp_dir, self.model_attc)
        for suffix in ('_attc.res', '_attc_table.res'):
            res = os.path.join(self.tmp_dir, self.replicon_name + suffix)
            self.assertTrue(os.path.exists(res))

    def test_find_attc_no_infernal(self):
        cmsearch_bin = 'foo'
        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, 'Replicons', replicon_name + '.fst')
        with self.assertRaises(RuntimeError) as ctx:
            infernal.find_attc(replicon_path, self.replicon_name, cmsearch_bin, self.tmp_dir, self.model_attc)
        self.assertTrue(re.search(r"failed : \[Errno 2\] No such file or directory: 'foo'", str(ctx.exception)),
                        msg=str(ctx.exception))

    def test_find_attc_no_model(self):
        model_attc = 'foo'
        with self.assertRaises(RuntimeError) as ctx:
            infernal.find_attc(self.replicon_path, self.replicon_name, self.cmsearch_path, self.tmp_dir, model_attc)
        self.assertTrue(str(ctx.exception).endswith('failed returncode = 1'))
