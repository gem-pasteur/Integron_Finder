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
from subprocess import run

__version__ = '2.0.3'


class IntegronError(Exception):
    pass


class EmptyFileError(IntegronError):
    pass


def get_logging_module():
    import colorlog
    try:
       logging = colorlog.logging.logging
    except AttributeError:
        logging = colorlog.wrappers.logging
    return logging


def _eddy_version(path):
    """
    execute binary from eddy's lab (hmmsearch or cmsearch) and parse version

    :param str path: the path to the binary to execute
    :return: the version
    :rtype: str
    """
    process = run([path, "-h"], capture_output=True, text=True)
    vers = process.stdout.split('\n')[1].strip()[2:]
    return vers


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
        version_text = """

                                                                                        
                                                                                        
                                                                                  
                                                                                        
                                            ██                                          
                                          ██░░██                                        
                                        ██░░░░░░██                               
                                      ██░░░░░░░░░░██                                    
                                    ██░░░░░░░░░░░░░░██           This version of Integron finder is buggy  
                                  ██░░░░░░██████░░░░░░██          https://integronfinder.readthedocs.io/en/latest/user_guide/changes.html
                                ██░░░░░░░░██████░░░░░░░░██         Except to compare results produced with this version               
                              ██░░░░░░░░░░██████░░░░░░░░░░██        with a other version.    
                            ██░░░░░░░░░░░░██████░░░░░░░░░░░░██          YOU SHOULD NOT USE IT.               
                          ██░░░░░░░░░░░░░░██████░░░░░░░░░░░░░░██                       
                        ██░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░██     Use the 2.0.5 or later instead                  
                        ██░░░░░░░░░░░░░░░░██████░░░░░░░░░░░░░░░░██                      
                      ██░░░░░░░░░░░░░░░░░░██████░░░░░░░░░░░░░░░░░░██                    
                      ██░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░██                    
                        ██████████████████████████████████████████                      
                                                                                                                                      

        
integron_finder version {i_f}

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
""".format(i_f=__version__,
           py=sys.version.replace('\n', ' '),
           np=np_vers,
           pd=pd_vers,
           mplt=mplt_vers,
           bio=bio_vers,
           prodigal=_prodigal_version(prodigal),
           cmsearch=_eddy_version(cmsearch),
           hmmsearch=_eddy_version(hmmsearch)
          )

    return version_text


def init_logger(log_file=None, out=True):
    import colorlog
    logger = colorlog.getLogger('integron_finder')
    logging = get_logging_module()

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
    logging = get_logging_module()
    levels = {'NOTSET': logging.NOTSET,
              'DEBUG': logging.DEBUG,
              'INFO': logging.INFO,
              'WARNING': logging.WARNING,
              'ERROR': logging.ERROR,
              'CRITICAL': logging.CRITICAL,
              }
    if level in levels:
        level = levels[level]
    elif not isinstance(level, int):
        raise IntegronError("Level must be {} or a positive integer")
    elif level < 0:
        raise IntegronError("Level must be {} or a positive integer")

    logger = colorlog.getLogger('integron_finder')
    if level <= logging.DEBUG:
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

        file_formatter = logging.Formatter("%(levelname)-8s : %(module)s: L %(lineno)d : %(message)s")
        file_handler = logger.handlers[1]
        file_handler.setFormatter(file_formatter)

    logger.setLevel(level)
