# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright (c) 2015 - 2018  Institut Pasteur, Paris and CNRS.                     #
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

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.scripts import finder


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
        from integron_finder import __version__ as if_vers

        with self.catch_io(out=True, err=False):
            try:
                finder.main(['integron_finder', '-V'])
            except TypeError as err:
                version = sys.stdout.getvalue()
                # program exit with returncode = 0
                self.assertEqual(str(err), '0')
            finally:
                sys.exit = real_exit

        exp_version = """integron_finder version {i_f}
using:    
 - Python {py}
 - numpy {np}
 - pandas {pd}
 - matplolib {mplt}
 - biopython {bio}
 
integron_finder is released under open-source GPLv3,
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions; see COPYING file for details. 

Authors:
 - Jean Cury, Bertrand Neron, Eduardo Rocha,

citation:

 Identification and analysis of integrons and cassette arrays in bacterial genomes
 Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha
 Nucleic Acids Research 2016; doi: 10.1093/nar/gkw319
 """.format(i_f=if_vers,
            py=sys.version.replace('\n', ' '),
            np=np_vers,
            pd=pd_vers,
            mplt=mplt_vers,
            bio=bio_vers
            )
        self.assertEqual(exp_version.strip(), version.strip())