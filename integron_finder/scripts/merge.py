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
import os
import argparse
import glob
import shutil

import integron_finder

# must be done after import 'integron_finder'
import colorlog

from integron_finder import IntegronError, logger_set_level
from integron_finder import utils
from integron_finder import results

_log = None


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
        if _log:
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
        agg_file.to_csv(out_file, sep="\t")
        return out_file


def copy_file(out_dir, ext, *from_dirs):
    """
    copy files from *from_dirs* and finishing with *ext* to the *out_dir* directory

    :param str out_dir: The path to the destination directory
    :param str ext: the extension of files to copy
    :param from_dirs: The path of the source directories
    :type from_dirs: list of str
    """
    for _dir in from_dirs:
        from_files = glob.glob(os.path.join(_dir, '*' + ext))
        for one_file in from_files:
            shutil.copy(one_file, out_dir)


def copy_dir(out_dir, pattern, *from_dirs):
    """
    Look inside directories _from_dir if some dir match the pattern (glob)
    and copy the last element of each matched path to out_dir

    :param out_dir: The path to the destination directory
    :param str pattern: pattern to match
    :param from_dirs: The path of the source directories
    :type from_dirs: list of str
    """
    for _dir in from_dirs:
        dirs = glob.glob(os.path.join(_dir, pattern))
        for one_dir in dirs:
            dest_dir = os.path.basename(one_dir)
            shutil.copytree(one_dir, os.path.join(out_dir, dest_dir))


def parse_args(args=None):
    """

    :param args: The arguments passed on the command line (without the name of the program)
                 Typically sys.argv[1:]
    :type args: list of string.
    :return: the arguments parsed.
    :rtype: a :class:`argparse.Namespace` object.
    """
    description = """Merge different integron_finder results in one
     
 - merge the '.integrons' files
 - merge the '.summary' files
 - copy the *.pdf files if they exist
 - copy the *.gbk file if they exist
 - copy the temporary directory if they exist

for instance to merge the results from 3 analysis
     
  - Results_Integron_Finder_psa_1
  - Results_Integron_Finder_psa_2
  - Results_Integron_Finder_psa_3
     
into new directory Results_Integron_Finder_foo_bar

integron_merge Results_Integron_Finder_foo_bar baz  Results_Integron_Finder_psa*
    
will give :

- one directory: Results_Integron_Finder_foo_bar

containing

    - baz.integrons
    - baz.summary
    - psa_1.gbk
    - psa_2.gbk
    - psa_3.gbk
"""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=description)
    parser.add_argument('outdir',
                        help='path of directory where the merged result must be stored')
    parser.add_argument('outfile',
                        help="the basename of the result files ('.integrons' and '.summary') "
                             "without extension, eg : data_set_1")
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

    if parsed_args.outdir in parsed_args.results:
        raise ValueError("'outdir' and 'results' cannot have the same value.")
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
    copy_dir(outdir, 'tmp_*', *parsed_args.results)


if __name__ == '__main__':
    main()
