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

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import logger_set_level
from integron_finder.hmm import scan_hmm_bank


class TestScanHmmBank(IntegronTest):

    def setUp(self):
        logger_set_level('INFO')
        self._tmp_dir = tempfile.TemporaryDirectory(prefix='tmp_test_integron_finder')
        self.tmp_dir = self._tmp_dir.name
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)

    def tearDown(self):
        self._tmp_dir.cleanup()

    def test_wrong_path(self):
        """
        Test that when the given path does not exist, the function raises an exception
        """
        with self.assertRaises(IOError) as ctx:
            scan_hmm_bank("wrong_path")
        self.assertEqual(str(ctx.exception), "wrong_path no such file or directory")

    def test_scan_dir(self):
        """
        Test that when the given path is a directory, it lists all hmm files from this
        directory (and not other files)
        """
        mypath = self.find_data("hmm_files")
        files = scan_hmm_bank(mypath)
        exp_files = ["integrase.hmm", "phage_int.hmm"]
        exp_files = [os.path.abspath(os.path.join(mypath, myfile)) for myfile in exp_files]
        self.assertEqual(set(exp_files), set(files))

    def test_scan_file_wrong(self):
        """
        Test that when the given path is a file, containing a directory name which does not
        exist, it returns a warning in stderr
        """
        mypath = os.path.join("tests", "data", "hmm_files", "list_hmm.txt")
        exp_msg = "func_annot '/my_hmms' does not match any files."
        with self.catch_log() as log:
            files = scan_hmm_bank(mypath)
            catch_msg = log.get_value().strip()

        self.assertEqual(catch_msg, exp_msg)
        self.assertEqual(files, [])


    def test_scan_file_wrong_content(self):
        """
        Test that when the given path is a file, containing something else than path to hmm file
        it raise a ValueError
        """
        mypath = self.find_data(os.path.join("Replicons", 'multi_fasta.fst'))
        with self.catch_log():
            with self.assertRaises(ValueError) as ctx:
                _ = scan_hmm_bank(mypath)
            self.assertTrue(str(ctx.exception).startswith("Too many lines with no hmm file in {}.".format(mypath)))


    def test_scan_file_names(self):
        """
        Test that when the given path is a file, containing:
        - directory names + *.hmm
        - a filename
        it returns the names of all hmm files contained in the directories + given
        """
        path1 = self.find_data("hmm_files")
        path2 = self.find_data(os.path.join("Models"))
        hmm_paths = [os.path.join(path1, "*.hmm"),
                     os.path.join(path2, "*.hmm"),
                     ]

        annot_fam = os.path.join(self.tmp_dir, 'annot_fam.HMM')
        open(annot_fam, 'w').close()
        hmm_paths.append(annot_fam)

        # Get all paths to include into the hmm_bank file
        abs_hmm = [os.path.abspath(path) for path in hmm_paths]
        # Write the hmm_bank file
        bank_file = os.path.join(self.tmp_dir, "list_hmm2.txt")
        with open(bank_file, "w") as lhmm:
            lhmm.write("# Directories containing hmm files:\n")
            lhmm.write(abs_hmm[0] + "\n")
            lhmm.write(abs_hmm[1] + "\n")
            lhmm.write("# An absolute path hmm file:\n")
            lhmm.write(abs_hmm[2] + "\n")
            lhmm.write("# " + abs_hmm[2] + ".hmm\n")
            lhmm.write("# A relative path to hmm file:\n")
            lhmm.write("# " + hmm_paths[2] + ".hmm\n")
        # Read hmm_bank file and get list of all files found
        with self.catch_log() as log:
            files = scan_hmm_bank(bank_file)
            catch_msg = log.get_value().strip()
        # Write expected list of hmm files
        exp_files1 = ["integrase.hmm", "phage_int.hmm"]
        exp_files1 = [os.path.abspath(os.path.join(path1, myfile)) for myfile in exp_files1]
        exp_files2 = ["integron_integrase.hmm", "phage-int.hmm"]
        exp_files2 = [os.path.abspath(os.path.join(path2, myfile)) for myfile in exp_files2]
        exp_files = exp_files1 + exp_files2 + [abs_hmm[2]]
        # Compare expected and output
        self.assertEqual(set(exp_files), set(files))
        out_stderr = ["the hmm {} will be used for functional annotation".format(path)
                      for path in exp_files]
        self.assertEqual(set(catch_msg.split("\n")), set(out_stderr))
