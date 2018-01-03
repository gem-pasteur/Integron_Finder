#!/usr/bin/env python
# coding: utf-8

"""
Unit tests scan_hmm_bank function of integron_finder
"""

import integron_finder
import unittest
import os
import sys
from contextlib import contextmanager
from StringIO import StringIO


class TestFunctions(unittest.TestCase):

    @contextmanager
    def catch_output(self):
        """
        Catch stderr and stdout of the code running with this function.
        They can, then, be compared to expected outputs.
        """
        new_out, new_err = StringIO(), StringIO()
        old_out, old_err = sys.stdout, sys.stderr
        try:
            sys.stdout, sys.stderr = new_out, new_err
            yield sys.stdout, sys.stderr
        finally:
            sys.stdout, sys.stderr = old_out, old_err

    def test_wrong_path(self):
        """
        Test that when the given path does not exist, the function raises an exception
        """
        with self.assertRaises(IOError) as exp:
            integron_finder.scan_hmm_bank("wrong_path")
        raised = exp.exception
        self.assertEqual(raised.message, "wrong_path no such file or directory")

    def test_scan_dir(self):
        """
        Test that when the given path is a directory, it lists all hmm files from this
        directory (and not other files)
        """
        mypath = os.path.join("tests", "data", "hmm_files")
        files = integron_finder.scan_hmm_bank(mypath)
        exp_files = ["Resfam.hmm", "integrase.hmm", "phage_int.hmm"]
        exp_files = [os.path.abspath(os.path.join(mypath, myfile)) for myfile in exp_files]
        self.assertEqual(set(exp_files), set(files))

    def test_scan_file_wrong(self):
        """
        Test that when the given path is a file, containing a directory name which does not
        exist, it returns a warning in stderr
        """
        mypath = os.path.join("tests", "data", "hmm_files", "list_hmm.txt")
        with self.catch_output() as (out, err):
            files = integron_finder.scan_hmm_bank(mypath)
        self.assertEqual(out.getvalue().strip(), "")
        self.assertEqual(err.getvalue().strip(), "WARNING func_annot '/my_hmms' does not " +
                                                  "match any files.")
        self.assertEqual(files, [])

    def test_scan_file_names(self):
        """
        Test that when the given path is a file, containing:
        - directory names + *.hmm
        - a filename
        it returns the names of all hmm files contained in the directories + given
        """
        path1 = os.path.join("tests", "data", "hmm_files")
        path2 = os.path.join("data", "Models")
        hmm_paths = [os.path.join(path1, "*.hmm"),
                     os.path.join(path2, "*.hmm"),
                     os.path.join("data", "Functional_annotation", "Resfams.hmm")]
        # Get all paths to include into the hmm_bank file
        abs_hmm = [os.path.abspath(path) for path in hmm_paths]
        # Write the hmm_bank file
        with open("list_hmm2.txt", "w") as lhmm:
            lhmm.write("# Directories containing hmm files:\n")
            lhmm.write(abs_hmm[0] + "\n")
            lhmm.write(abs_hmm[1] + "\n")
            lhmm.write("# A hmm file:\n")
            lhmm.write(abs_hmm[2] + "\n")
            lhmm.write("# " + abs_hmm[2] + ".hmm\n")
        # Read hmm_bank file and get list of all files found
        with self.catch_output() as (out, err):
            files = integron_finder.scan_hmm_bank("list_hmm2.txt")
        # Write expected list of hmm files
        exp_files1 = ["Resfam.hmm", "integrase.hmm", "phage_int.hmm"]
        exp_files1 = [os.path.abspath(os.path.join(path1, myfile)) for myfile in exp_files1]
        exp_files2 = ["integron_integrase.hmm", "phage-int.hmm"]
        exp_files2 = [os.path.abspath(os.path.join(path2, myfile)) for myfile in exp_files2]
        exp_files = exp_files1 + exp_files2 + [abs_hmm[2]]
        # Compare expected and output
        self.assertEqual(set(exp_files), set(files))
        out_stderr = ["the hmm {} will be used for functional annotation".format(path)
                      for path in exp_files]
        self.assertEqual(out.getvalue().strip(), "")
        self.assertEqual(set(err.getvalue().strip().split("\n")), set(out_stderr))
        # Remove hmm_bank file
        os.remove("list_hmm2.txt")
