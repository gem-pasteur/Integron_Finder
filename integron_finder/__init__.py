# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright (c) 2015 - 2021  Institut Pasteur, Paris and CNRS.                     #
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
from subprocess import run
import distutils.spawn

from time import localtime, strftime

__version__ = f'2.dev{strftime("%Y%m%d", localtime())}'


class IntegronError(Exception):
    pass


class EmptyFileError(IntegronError):
    pass


def _eddy_version(path):
    process = run([path, "-h"], capture_output=True, text=True)
    vers = process.stdout.split('\n')[1].strip()[2:]
    return vers

def _hmmsearch_version(path):
    """

    :return:
    """
    return _eddy_version(path)


def _cmsearch_version(path):
    """

    :return:
    """
    return _eddy_version(path)


def _prodigal_version(path):
    """

    :return:
    """
    process = run([path, "-v"], capture_output=True, text=True)
    vers = process.stderr.strip()
    return vers



def get_version_message(hmmsearch, cmsearch, prodigal):
    from numpy import __version__ as np_vers
    from pandas import __version__ as pd_vers
    from matplotlib import __version__ as mplt_vers
    from Bio import __version__ as bio_vers

    third_parties = {'hmmsearch': hmmsearch, 'cmsearch':cmsearch, 'prodigal':prodigal}
    if not all(third_parties.values()):
        version_text = ''
        for prog_name, value in third_parties.items():
            if not value:
                version_text += f"{prog_name} not found in Path please use --{prog_name} option to specify it.\n"

    else:
        version_text = """integron_finder version {i_f}
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

 Identification and analysis of integrons and cassette arrays in bacterial genomes
 Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha
 Nucleic Acids Research 2016; doi: 10.1093/nar/gkw319

 If you use --func-annot in conjunction with file NCBIfam-AMRFinder.hmm please also cite

 Haft, DH et al., Nucleic Acids Res. 2018 Jan 4;46(D1):D851-D860
 PMID: 29112715
""".format(i_f=__version__,
           py=sys.version.replace('\n', ' '),
           np=np_vers,
           pd=pd_vers,
           mplt=mplt_vers,
           bio=bio_vers,
           prodigal=_prodigal_version(prodigal),
           cmsearch=_cmsearch_version(cmsearch),
           hmmsearch=_hmmsearch_version(hmmsearch)
          )

    return version_text


def init_logger(log_file=None, out=True):
    import colorlog
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


def logger_set_level(level='WARNING'):
    # default value must be a string
    # cannot be colorlog.logging.logging.WARNING for instance
    # because setup import __init__ to get __version__
    # so logger_set_level is defined
    # if level is colorlog.logging.logging.WARNING
    # that mean that colorlog must be already installed
    # otherwise an error occured during pip install
    #  NameError: name 'colorlog' is not defined
    import colorlog

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
