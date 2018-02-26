#!/usr/bin/env python
# coding: utf-8


"""
Unit tests get_version_message function of integron_finder
"""
import sys

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


class TestGetVersion(IntegronTest):

    def test_get_version_not_packaged(self):
        """
        test on having the version message when integron_finder is not installed
        """
        real_exit = sys.exit
        sys.exit = self.fake_exit

        real_version = integron_finder.__version__
        forced_vers = '$' + 'VERSION'
        integron_finder.__version__ = forced_vers

        with self.catch_io(out=True, err=True):
            try:
                finder.main(['integron_finder', '-V'])
            except TypeError as err:
                version = sys.stderr.getvalue()
                # program exit with returncode = 0
                self.assertEqual(str(err), '0')
            finally:
                sys.exit = real_exit
                integron_finder.__version__ = real_version

        exp_version = """integron_finder version NOT packaged, it should be development version
Python {0}

 - open-source GPLv3,
 - Jean Cury, Bertrand Neron, Eduardo Rocha,
 - citation:

 Identification and analysis of integrons and cassette arrays in bacterial genomes
 Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha
 Nucleic Acids Research 2016; doi: 10.1093/nar/gkw319
 """.format(sys.version, forced_vers)

        self.assertEqual(exp_version.strip(), version.strip())


    def test_get_version_installed(self):
        """
        test on having the version message when integron_finder is installed
        """
        real_exit = sys.exit
        sys.exit = self.fake_exit

        real_version = integron_finder.__version__
        forced_vers = '1.5.2'
        integron_finder.__version__ = forced_vers

        with self.catch_io(out=True, err=True):
            try:
                finder.main(['integron_finder', '-V'])
            except TypeError as err:
                version = sys.stderr.getvalue()
                # program exit with returncode = 0
                self.assertEqual(str(err), '0')
            finally:
                sys.exit = real_exit
                integron_finder.__version__ = real_version

        exp_version = """integron_finder version {1}
Python {0}

 - open-source GPLv3,
 - Jean Cury, Bertrand Neron, Eduardo Rocha,
 - citation:

 Identification and analysis of integrons and cassette arrays in bacterial genomes
 Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha
 Nucleic Acids Research 2016; doi: 10.1093/nar/gkw319
 """.format(sys.version, forced_vers)
        self.assertEqual(exp_version.strip(), version.strip())