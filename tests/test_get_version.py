#!/usr/bin/env python
# coding: utf-8


"""
Unit tests get_version_message function of integron_finder
"""
import unittest
import sys
import os
import subprocess

class TestFunctions(unittest.TestCase):

    def check_installed():
        """
        Check that integron_finder is installed. If not, will not run test for version
        message while installed
        """
        try:
            subprocess.call(["integron_finder"])
        except OSError:
            return False
        return True

    @unittest.skipIf(check_installed(), "integron_finder package not installed")
    def test_get_version_not_packaged(self):
        """
        test on having the version message when integron_finder is not installed
        """
        pathname = os.path.dirname(sys.argv[0])
        os.environ['INTEGRON_HOME'] = os.path.join(os.path.abspath(pathname), "..")
        p = subprocess.Popen(["./integron_finder", "-V"], stderr=subprocess.PIPE)
        version = p.communicate()[1]

        exp_version = """integron_finder version NOT packaged, it should be development version
Python {0}

 - open-source GPLv3,
 - Jean Cury, Bertrand Neron, Eduardo Rocha,
 - citation:

 Identification and analysis of integrons and cassette arrays in bacterial genomes
 Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha
 Nucleic Acids Research 2016; doi: 10.1093/nar/gkw319
 """.format(sys.version).split("\n")
        vers = version.split("\n")
        # remove the 2 first lines, corresponding to the warning for biopython, matplotlib etc.
        myversion = []
        for line in vers:
            if "warning" not in line.lower():
                myversion.append(line)
        for exp_part, part in zip(exp_version, myversion):
            self.assertEqual(exp_part, part)

    @unittest.skipIf(not check_installed(), "integron_finder package not installed")
    def test_get_version_installed(self):
        """
        test on having the version message when integron_finder is installed
        """
        int_vers = subprocess.check_output(["python", "setup.py", "--version"])
        p = subprocess.Popen("integron_finder -V",
                             stderr=subprocess.PIPE,
                             shell=True)
        version = p.communicate()[1]
        exp_version = """integron_finder version {1}
Python {0}

 - open-source GPLv3,
 - Jean Cury, Bertrand Neron, Eduardo Rocha,
 - citation:

 Identification and analysis of integrons and cassette arrays in bacterial genomes
 Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha
 Nucleic Acids Research 2016; doi: 10.1093/nar/gkw319
 """.format(sys.version, int_vers.strip()).split("\n")
        vers = version.split("\n")
        # remove the 2 first lines, corresponding to the warning for biopython
        myversion = []
        for line in vers:
            if "warning" not in line.lower():
                myversion.append(line)
        for exp_part, part in zip(exp_version, myversion):
            self.assertEqual(exp_part, part)
