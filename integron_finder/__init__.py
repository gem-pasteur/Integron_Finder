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

import sys
import colorlog

__version__ = '$VERSION'

__INTEGRON_DATA__ = '$INTEGRONDATA'


class IntegronError(Exception):
    pass


def get_version_message():
    # if week keep '$VERSION' as is
    # the setup.py will replace it by the value set in setup
    # so the test become True even if integron_finder is installed using setup.py
    if __version__ == '$' + 'VERSION':
        version = "NOT packaged, it should be development version"
    else:
        version = __version__
    version_text = """integron_finder version {0}
Python {1}

 - open-source GPLv3,
 - Jean Cury, Bertrand Neron, Eduardo Rocha,
 - citation:

 Identification and analysis of integrons and cassette arrays in bacterial genomes
 Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha
 Nucleic Acids Research 2016; doi: 10.1093/nar/gkw319
 """.format(version, sys.version)
    return version_text


def init_logger(out=sys.stderr):
    handler = colorlog.StreamHandler(out)
    formatter = colorlog.ColoredFormatter("%(log_color)s%(levelname)-8s : %(reset)s %(message)s",
                                          datefmt=None,
                                          reset=True,
                                          log_colors={
                                              'DEBUG':    'cyan',
                                              'INFO':     'green',
                                              'WARNING':  'yellow',
                                              'ERROR':    'red',
                                              'CRITICAL': 'bold_red',
                                          },
                                          secondary_log_colors={},
                                          style='%'
                                          )
    handler.setFormatter(formatter)
    logger = colorlog.getLogger('integron_finder')
    logger.addHandler(handler)
    logger.setLevel(colorlog.logging.logging.WARNING)

init_logger()


def logger_set_level(level=colorlog.logging.logging.WARNING):
    logger = colorlog.getLogger('integron_finder')
    if level <= colorlog.logging.logging.DEBUG:
        formatter = colorlog.ColoredFormatter(
            "%(log_color)s%(levelname)-8s : %(module)s: L %(lineno)d :%(reset)s %(message)s",
            datefmt=None,
            reset=True,
            log_colors={
                'DEBUG': 'cyan',
                'INFO': 'green',
                'WARNING': 'yellow',
                'ERROR': 'red',
                'CRITICAL': 'bold_red',
            },
            secondary_log_colors={},
            style='%'
            )
        handler = logger.handlers[0]
        handler.setFormatter(formatter)

    logger.setLevel(level)
