#!/usr/bin/env python3
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
import os
import argparse
import glob
import shutil

import integron_finder
if not integron_finder.__version__.endswith('VERSION'):
    # display warning only for non installed integron_finder
    from Bio import BiopythonExperimentalWarning
    import warnings
    warnings.simplefilter('ignore', FutureWarning)
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)

# must be done after import 'integron_finder'
import colorlog

from integron_finder import IntegronError, logger_set_level
from integron_finder import utils
from integron_finder import  results


def merge_integrons(out_file, *in_dirs):
    """

    :param in_dirs: The path of the source directories
    :type in_dirs: list of str
    :param str out_file: The path to the merged file
    :return: The The path to the merged file
    """
    integrons_files = []
    for _dir in in_dirs:
        in_files = glob.glob(os.path.join(_dir, '*' + '.integrons'))
        integrons_files.extend(in_files)
    if integrons_files:
        agg_file = results.merge_results(*integrons_files)
        agg_file.to_csv(out_file, index=False, sep="\t", na_rep="NA")
        return out_file
    else:
        msg = "No integrons file to merge"
        _log.critical(msg)
        raise IntegronError(msg)


def merge_summary(out_file, *in_dirs):
    """

    :param in_dirs: The path of the source directories
    :type in_dirs: list of str
    :param str out_file: The path to the merged file
    :return: The The path to the merged file
    """
    summaries_files = []
    for _dir in in_dirs:
        in_files = glob.glob(os.path.join(_dir, '*' + '.summary'))
        summaries_files.extend(in_files)
    if summaries_files:
        agg_file = results.merge_results(*summaries_files)
        agg_file.to_csv(out_file, sep="\t", na_rep="NA",
                        columns=['ID_replicon', 'ID_integron', 'complete', 'In0', 'CALIN'],
                        index=False)
        return out_file


def copy_file(out_dir, ext, *in_dirs):
    """
    copy files in *in_dirs* and finishing with *ext* to the *out_dir* directory

    :param in_dirs: The path of the source directories
    :type in_dirs: list of str
    :param str out_dir: The path to the destination directory
    :param str ext: the extension of files to copy
    """
    for _dir in in_dirs:
        in_files = glob.glob(os.path.join(_dir, '*' + ext))
        for one_file in in_files:
            shutil.copy(one_file, out_dir)


def parse_args(args=None):
    """

    :param args: The arguments passed on the command line (without the name of the program)
                 Typically sys.argv[1:]
    :type args: list of string.
    :return: the arguments parsed.
    :rtype: a :class:`argparse.Namespace` object.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('outdir',
                        help='path of the merged result directory')
    parser.add_argument('outfile',
                        help='the basename of the result file without extension, eg : file1')
    parser.add_argument('results',
                        nargs='+',
                        help='Path to the results dir to merge eg : path/to/ Results_Integron_Finder_acba.007.p01.13 '
                             'path/to/Results_Integron_Finder_lian.001.c02.10')

    verbosity_grp = parser.add_argument_group()
    verbosity_grp.add_argument('-v', '--verbose',
                               action='count',
                               default=0,
                               help='Increase verbosity of output (can be cumulative : -vv)')
    verbosity_grp.add_argument('-q', '--quiet',
                               action='count',
                               default=0,
                               help='Decrease verbosity of output (can be cumulative : -qq)'
                               )
    parsed_args = parser.parse_args(args)
    return parsed_args


def main(args=None, log_level=None):
    """
     main entry point to integron_merge

    :param str args: the arguments passed on the command line
    :param log_level: the output verbosity
    :type log_level: a positive int or a string among 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    """
    global _log

    args = sys.argv[1:] if args is None else args
    parsed_args = parse_args(args)

    integron_finder.init_logger()
    _log = colorlog.getLogger('integron_finder.merge')

    if not log_level:
        # logs are specify from args options
        if not log_level:
            # logs are specify from args options
            logger_set_level(utils.log_level(parsed_args.verbose, parsed_args.quiet))
    else:
        # used by unit tests to mute or unmute logs
        logger_set_level(log_level)

    outdir = os.path.realpath(parsed_args.outdir)
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            msg = "'{}' already exists and is not a directory".format(outdir)
            _log.critical(msg)
            raise IOError(msg)
    else:
        os.makedirs(parsed_args.outdir)

    integron_file_out = os.path.join(outdir, parsed_args.outfile + ".integrons")
    merge_integrons(integron_file_out, *parsed_args.results)
    summary_file_out = os.path.join(outdir, parsed_args.outfile + ".summary")
    merge_summary(summary_file_out, *parsed_args.results)
    copy_file(outdir, '.gbk', *parsed_args.results)
    copy_file(outdir, '.pdf', *parsed_args.results)


if __name__ == '__main__':
    main()