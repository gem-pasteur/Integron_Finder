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

import sys
import shutil

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)


from integron_finder.scripts import finder
import integron_finder

class TestGetVersion(IntegronTest):


    def test_get_version(self):
        """
        test on having the version message when integron_finder is installed
        """
        real_exit = sys.exit
        sys.exit = self.fake_exit

        from numpy import __version__ as np_vers
        from pandas import __version__ as pd_vers
        from matplotlib import __version__ as mplt_vers
        from Bio import __version__ as bio_vers

        with self.catch_io(out=True, err=False):
            try:
                finder.main(['integron_finder', '--version'])
            except TypeError as err:
                version = sys.stdout.getvalue()
                # program exit with returncode = 0
                self.assertEqual(str(err), '0')
            finally:
                sys.exit = real_exit

        exp_version = """integron_finder version {i_f} {commit}
Using:
 - Python {py}
 - numpy {np}
 - pandas {pd}
 - matplolib {mplt}
 - biopython {bio}

 - {prodigal}
 - {cmsearch}
 - {hmmsearch}

Authors:
 - Jean Cury, Bertrand Neron, Eduardo Rocha,

Citation:

 Néron, B.; Littner, E.; Haudiquet, M.; Perrin, A.; Cury, J.; Rocha, E.P.C.
 IntegronFinder 2.0: Identification and Analysis of Integrons across Bacteria, with a Focus on Antibiotic Resistance in Klebsiella.
 Microorganisms 2022, 10, 700. https://doi.org/10.3390/microorganisms10040700

 If you use --func-annot in conjunction with file NCBIfam-AMRFinder.hmm please also cite

 Haft, DH et al., Nucleic Acids Res. 2018 Jan 4;46(D1):D851-D860
 PMID: 29112715
""".format(i_f=integron_finder.__version__ ,
           commit=integron_finder.__commit__ if 'dev' in integron_finder.__version__ else '',
            py=sys.version.replace('\n', ' '),
            np=np_vers,
            pd=pd_vers,
            mplt=mplt_vers,
            bio=bio_vers,
            prodigal=integron_finder._prodigal_version(shutil.which("prodigal")),
            cmsearch=integron_finder._eddy_version(shutil.which("cmsearch")),
            hmmsearch=integron_finder._eddy_version(shutil.which("hmmsearch"))
            )
        self.assertEqual(exp_version, version)


    def test_get_version_no_binary(self):
        # use directly get_version_message
        # to simulate parser do not find binaries
        # shuti.which return None

        ##################################
        # case there is no binary at all #
        ##################################
        version = integron_finder.get_version_message(None,  # hmmsearch
                                                      None,  # cmsearch
                                                      None)  # prodigal

        exp_version = """hmmsearch not found in Path please use --hmmsearch option to specify it.
cmsearch not found in Path please use --cmsearch option to specify it.
prodigal not found in Path please use --prodigal option to specify it.
"""
        self.assertEqual(exp_version, version)

        #########################################
        # case there cmsearch has been  found   #
        # but not hmmsearch nor prodigal        #
        #########################################
        version = integron_finder.get_version_message(None,  # hmmsearch
                                                      shutil.which('cmsearch'),
                                                      None)  # prodigal

        exp_version = """hmmsearch not found in Path please use --hmmsearch option to specify it.
prodigal not found in Path please use --prodigal option to specify it.
"""
        self.assertEqual(exp_version, version)
