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
import math
import os
import argparse
from itertools import zip_longest
import re

from Bio import SeqIO

import integron_finder
# must be done after import 'integron_finder'
import colorlog

from integron_finder import logger_set_level
from integron_finder import utils


def split(replicon_path, chunk=None, outdir='.'):
    """
    Split the replicon_file in *chunk* chunks and write them in files.
    the name of the chunk is the input filename with suffix '_chunk_i'
    (where i is the chunk number) if the chunk contains several sequences
    or the id of the sequence if there is only one sequence in the chunk.
    There also a system that prevent to over write an existing file by appending (number)
    to the file name for instance ESCO001.B.00018.P002_(1).fst

    :param str replicon_path: The path to the replicon file.
    :param int chunk: The number of chunk desire (chunk > 0).
    :param str outdir: The path of a directory where to write chunk files.
                       The directory must exists.
    :return: The name of all chunks created.
    :rtype: List of strings.
    """
    def grouper(sequences_db, chunk_size):
        """

        :param sequences_db: The sequences to group
        :type sequences_db: A :class:`integron_finder.utils.FastaIterator` object.
        :param int chunk_size: The number of sequence by Chunk file.
        :return: a chunk of sequences.
        :rtype: An iterator of tuples.
        """
        args = [iter(sequences_db)] * chunk_size
        return zip_longest(*args)

    with utils.FastaIterator(replicon_path) as sequences_db:
        sequences_db_len = len(sequences_db)
        if not chunk:
            chunk_size = 1
        else:
            chunk_size = math.ceil(sequences_db_len / chunk)

        chunks = grouper(sequences_db, chunk_size)
        all_chunk_name = []
        for chunk_no, chunk_in in enumerate(chunks, 1):
            # if replicon contains illegal characters
            # or replicon is too short < 50 bp
            # then replicon is None
            chunk_out = []
            for rep_no, replicon in enumerate(chunk_in, 1):
                if replicon is not None:
                    replicon_name = replicon.id
                    chunk_out.append(replicon)
                else:
                    rep_no_in_db = (chunk_no - 1) * chunk_size + rep_no
                    if rep_no_in_db <= sequences_db_len:
                        _log.warning("Skipping replicon {}/{} in chunk {}".format(rep_no_in_db,
                                                                                  sequences_db_len,
                                                                                  chunk_no))
            if chunk_out:
                if chunk_size == 1:
                    chunk_name = "{}.fst".format(replicon_name)
                else:
                    replicon_name = utils.get_name_from_path(replicon_path)
                    chunk_name = "{}_chunk_{}.fst".format(replicon_name, chunk_no)
                chunk_name = os.path.join(outdir, chunk_name)
                i = 0
                while os.path.exists(chunk_name):
                    root, ext = os.path.splitext(chunk_name)
                    i += 1
                    match = re.search(r"_chunk_\d+$", root)
                    if match:
                        root = root[:match.start()]
                    chunk_name = "{}_chunk_{}{}".format(root, i, ext)

                _log.info("writing chunk '{}'".format(chunk_name))
                SeqIO.write(chunk_out, chunk_name, "fasta")
                all_chunk_name.append(chunk_name)
    return all_chunk_name


def parse_args(args):
    """

    :param args: The arguments passed on the command line (without the name of the program)
                 Typically sys.argv[1:]
    :type args: list of string.
    :return: the arguments parsed.
    :rtype: a :class:`argparse.Namespace` object.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('replicon',
                        help='Path to the replicon file (in fasta format), eg : path/to/file.fst or file.fst')

    parser.add_argument('--chunk',
                        type=int,
                        help='The number of files generate. '
                             'Each chunks will contains n replicon where '
                             'n = number of replicon in the input file / chunk.'
                             'The n may vary in some chunks because some replicon can be skip '
                             'if they contains illegal characters or are too short (<50bp)')

    parser.add_argument('-o', '--outdir',
                        default='.',
                        help='The path to the directory where to write the chunks.\n'
                             'If outdir doe not exists, it will be created.')
    parser.add_argument("--mute",
                        action='store_true',
                        default=False,
                        help="mute the log on stdout."
                             "(continue to log on integron_split.out)")

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
     main entry point to integron_split

    :param str args: the arguments passed on the command line
    :param log_level: the output verbosity
    :type log_level: a positive int or a string among 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    """
    global _log

    args = sys.argv[1:] if args is None else args
    parsed_args = parse_args(args)
    if not os.path.exists(parsed_args.outdir):
        os.makedirs(parsed_args.outdir)
    elif not os.path.isdir(parsed_args.outdir):
        raise RuntimeError(f"The outdir '{parsed_args.outdir}' already exist and is not a directory.")
    integron_finder.init_logger(log_file=os.path.join(parsed_args.outdir, 'integron_split.out'),
                                out=not parsed_args.mute)
    _log = colorlog.getLogger('integron_finder.split')

    if not log_level:
        # logs are specify from args options
        logger_set_level(utils.log_level(parsed_args.verbose, parsed_args.quiet))
    else:
        # used by unit tests to mute or unmute logs
        logger_set_level(log_level)

    chunk_names = split(parsed_args.replicon, chunk=parsed_args.chunk, outdir=parsed_args.outdir)
    print(' '.join(chunk_names))


if __name__ == '__main__':
    main()
