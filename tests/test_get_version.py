#!/usr/bin/env python
# coding: utf-8


"""
Unit tests get_version_message function of integron_finder
"""
import sys
import subprocess
import unittest
try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

import integron_finder
from integron_finder.scripts import finder


def check_installed():
    """
    Check that integron_finder is installed. If not, will not run test for version
    message while installed
    """
    installed = False if integron_finder.__version__.endswith('VERSION') else True
    return installed


class TestFunctions(IntegronTest):

    @unittest.skipIf(check_installed(), "integron_finder package not installed")
    def test_get_version_not_packaged(self):
        """
        test on having the version message when integron_finder is not installed
        """
        real_exit = sys.exit

        sys.exit = self.fake_exit
        with self.catch_io(out=True, err=True):
            try:
                finder.main(['integron_finder', '-V'])
            except TypeError as err:
                version = sys.stderr.getvalue()
                # program exit with returncode = 0
                self.assertEqual(str(err), '0')
            finally:
                sys.exit = real_exit

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
        real_exit = sys.exit

        sys.exit = self.fake_exit
        with self.catch_io(out=True, err=True):
            try:
                finder.main(['integron_finder', '-V'])
            except TypeError as err:
                version = sys.stderr.getvalue()
                # program exit with returncode = 0
                self.assertEqual(str(err), '0')
            finally:
                sys.exit = real_exit

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
