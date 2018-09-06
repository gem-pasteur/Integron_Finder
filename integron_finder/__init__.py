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

__version__ = '2.0rc3'
__INTEGRON_DATA__ = '$INTEGRONDATA'


class IntegronError(Exception):
    pass


def get_version_message():
    from numpy import __version__ as np_vers
    from pandas import __version__ as pd_vers
    from matplotlib import __version__ as mplt_vers
    from Bio import __version__ as bio_vers
    version_text = """integron_finder version {i_f}
Using:    
 - Python {py}
 - numpy {np}
 - pandas {pd}
 - matplolib {mplt}
 - biopython {bio}

Authors:
 - Jean Cury, Bertrand Neron, Eduardo Rocha,

Citation:

 Identification and analysis of integrons and cassette arrays in bacterial genomes
 Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha
 Nucleic Acids Research 2016; doi: 10.1093/nar/gkw319
 """.format(i_f=__version__,
            py=sys.version.replace('\n', ' '),
            np=np_vers,
            pd=pd_vers,
            mplt=mplt_vers,
            bio=bio_vers
    )
    return version_text


def init_logger(log_file=None, out=True):
    logger = colorlog.getLogger('integron_finder')
    logging = colorlog.logging.logging
    if out:
        stdout_handler = colorlog.StreamHandler(sys.stdout)
        stdout_formatter = colorlog.ColoredFormatter("%(log_color)s%(levelname)-8s : %(reset)s %(message)s",
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
        stdout_handler.setFormatter(stdout_formatter)
        logger.addHandler(stdout_handler)
    else:
        null_handler = logging.NullHandler()
        logger.addHandler(null_handler)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_formatter = logging.Formatter("%(levelname)-8s : %(message)s")
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    logger.setLevel(logging.WARNING)


def logger_set_level(level=colorlog.logging.logging.WARNING):

    levels = {'NOTSET': colorlog.logging.logging.NOTSET,
              'DEBUG': colorlog.logging.logging.DEBUG,
              'INFO': colorlog.logging.logging.INFO,
              'WARNING': colorlog.logging.logging.WARNING,
              'ERROR': colorlog.logging.logging.ERROR,
              'CRITICAL': colorlog.logging.logging.CRITICAL,
              }
    if level in levels:
        level = levels[level]
    elif not isinstance(level, int):
        raise IntegronError("Level must be {} or a positive integer")
    elif level < 0:
        raise IntegronError("Level must be {} or a positive integer")

    logger = colorlog.getLogger('integron_finder')
    if level <= colorlog.logging.logging.DEBUG:
        stdout_formatter = colorlog.ColoredFormatter(
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
        stdout_handler = logger.handlers[0]
        stdout_handler.setFormatter(stdout_formatter)

        logging = colorlog.logging.logging
        file_formatter = logging.Formatter("%(levelname)-8s : %(module)s: L %(lineno)d : %(message)s")
        file_handler = logger.handlers[1]
        file_handler.setFormatter(file_formatter)

    logger.setLevel(level)
