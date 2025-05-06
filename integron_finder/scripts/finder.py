#!/usr/bin/env python
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


"""
integron_finder is a program that looks for integron in DNA sequences.
"""
import os
import sys
import argparse
import shutil

import pandas as pd
pd.options.mode.chained_assignment = 'raise'
import integron_finder

# must be done after import 'integron_finder'
import colorlog
_log = colorlog.getLogger('integron_finder')

from Bio import SeqIO

from integron_finder import IntegronError, logger_set_level
from integron_finder import utils
from integron_finder import results
from integron_finder import __commit__ as _if_commit , __version__ as _if_version
from integron_finder.topology import Topology
from integron_finder.config import Config
from integron_finder.hmm import scan_hmm_bank
from integron_finder.integrase import find_integrase
from integron_finder.attc import find_attc_max
from integron_finder.infernal import find_attc
from integron_finder.integron import find_integron
from integron_finder.annotation import func_annot, add_feature
from integron_finder.prot_db import GembaseDB, ProdigalDB, CustomDB
from integron_finder import argparse_utils


def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("replicon",
                        type=argparse_utils.path,
                        help="Path to the replicon file (in fasta format), eg : path/to/file.fst or file.fst")

    parser.add_argument("--local-max",
                        default=False,
                        help="Allows thorough local detection "
                             "(slower but more sensitive and do not increase false positive rate).",
                        action="store_true")

    parser.add_argument("--func-annot",
                        help="Functional annotation of CDS associated "
                             "with integrons HMM files are needed in Func_annot folder.",
                        default=False,
                        action="store_true")

    parser.add_argument('--cpu',
                        default='1',
                        type=int,
                        help='Number of CPUs used by INFERNAL and HMMER. Increasing too much (usually above 4) may decrease performance.')

    parser.add_argument('-dt', '--distance-thresh',
                        dest='distance_threshold',
                        default=4000,
                        action='store',
                        type=int,
                        help='Two elements are aggregated if they are distant of DISTANCE_THRESH [4000]bp or less')

    parser.add_argument('--outdir',
                        default=".",
                        help='Set the output directory (default: current)')

    parser.add_argument("--union-integrases",
                        default=False,
                        help="Instead of taking intersection of hits from Phage_int profile (Tyr recombinases)"
                             " and integron_integrase profile, use the union of the hits",
                        action="store_true")

    parser.add_argument('--cmsearch',
                        default=shutil.which("cmsearch"),
                        type=argparse_utils.path,
                        help='Complete path to cmsearch if not in PATH. eg: /usr/local/bin/cmsearch')

    parser.add_argument('--hmmsearch',
                        default=shutil.which("hmmsearch"),
                        type=argparse_utils.path,
                        help='Complete path to hmmsearch if not in PATH. eg: /usr/local/bin/hmmsearch')

    parser.add_argument('--prodigal',
                        default=shutil.which("prodigal"),
                        type=argparse_utils.path,
                        help='Complete path to prodigal if not in PATH. eg: /usr/local/bin/prodigal')

    parser.add_argument('--path-func-annot',
                        type=argparse_utils.path,
                        help='Path to file containing all hmm bank paths (one per line)')

    parser.add_argument("--gembase",
                        default=False,
                        help="Use gembase formatted protein file instead of Prodigal."
                             " Folder structure must be preserved",
                        action="store_true")
    parser.add_argument("--gembase-path",
                        type=argparse_utils.path,
                        help="path to the gembase root directory (needed only if the replicon file is not located"
                             "in gembase-path)")
    parser.add_argument("--annot-parser",
                        type=argparse_utils.path,
                        help="the path to the parser to use to get information from protein file.")
    parser.add_argument("--prot-file",
                        type=argparse_utils.path,
                        help="The path to the proteins file used for annotations")
    parser.add_argument('--attc-model',
                        default='attc_4.cm',
                        help='Path or file to the attc model (Covariance Matrix).')

    parser.add_argument('--evalue-attc',
                        default=1.,
                        type=float,
                        help='Set evalue threshold to filter out hits above it (default: 1)')

    parser.add_argument('--calin-threshold',
                        default=2,
                        type=int,
                        help="keep 'CALIN' only if attC sites number >= calin-threshold (default: 2)")

    parser.add_argument("--keep-palindromes",
                        default=False,
                        help="For a given hit, if the palindromic version is found,"
                             " don't remove the one with highest evalue.",
                        action="store_true")

    parser.add_argument("--no-proteins",
                        help="Don't annotate CDS and don't find integrase, just look for attC sites.",
                        default=False,
                        action="store_true")

    parser.add_argument("--promoter-attI",
                        help="Search also for promoter and attI sites. (default False)",
                        default=False,
                        action="store_true")

    parser.add_argument('--max-attc-size',
                        default=200,
                        type=int,
                        help='Set maximum value fot the attC size (default: 200bp)')

    parser.add_argument('--min-attc-size',
                        default=40,
                        type=int,
                        help='set minimum value fot the attC size (default: 40bp)')

    parser.add_argument("--eagle-eyes",
                        default=False,
                        help="Synonym of --local-max. Like a soaring eagle in the sky,"
                             " catching rabbits (or attC sites) by surprise.",
                        action="store_true")
    output_options = parser.add_argument_group("Output options")
    output_options.add_argument('--pdf',
                                action='store_true',
                                default=False,
                                help='For each complete integron, a simple graphic of the region is depicted '
                                     '(in pdf format)'),
    output_options.add_argument('--gbk',
                                action='store_true',
                                default=False,
                                help='generate a GenBank file with the sequence annotated with the same annotations '
                                     'than .integrons file.')
    output_options.add_argument('--keep-tmp',
                                action='store_true',
                                default=False,
                                help='keep intermediate results. '
                                     'This results are stored in directory named tmp_<replicon id>')
    output_options.add_argument('--split-results',
                                action='store_true',
                                default=False,
                                help='Instead of merging integron results from all replicon in one file, '
                                     'keep them in separated files.')

    topology_grp = parser.add_mutually_exclusive_group()
    topology_grp.add_argument("--circ",
                              dest='circular',
                              default=False,
                              help="Set the default topology for replicons to 'circular'",
                              action="store_true")
    topology_grp.add_argument("--linear",
                              default=False,
                              help="Set the default topology for replicons to 'linear'",
                              action="store_true")
    parser.add_argument("--topology-file",
                        type=argparse_utils.path,
                        help="The path to a file where the topology for each replicon is specified.")

    parser.add_argument("--version",
                        action=argparse_utils.MyVersionAction)

    parser.add_argument("--mute",
                        action='store_true',
                        default=False,
                        help="mute the log on stdout."
                             "(continue to log on integron_finder.out)")

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

    # to ensure that --version is at the end of cmd line
    # to allow to parse cmsearch hmmsearch and prodigal options
    # before to display version
    # otherwise the version displayed used the default path for these third parties programs
    vers = False
    for i, item in enumerate(args):
        if item.startswith('--vers'):
            vers = True
            break
    if vers:
        vers = args.pop(i)
        args.append(vers)

    parsed_args = parser.parse_args(args)

    # eagle_eyes is just an alias to local_max in whole program use local_max
    parsed_args.local_max = parsed_args.local_max or parsed_args.eagle_eyes
    return Config(parsed_args)


def find_integron_in_one_replicon(replicon, config):
    """
    scan replicon for integron.

      * presence of integrase
      * presence of attC sites
      * presence of promoters and attI sites

    depending on the configuration

     * perform functional annotation

    produce a file containing presence of putative integrons

    depending on configuration

        * produce genbank file with replicon and annotations with integrons
        * produce schema of replicon with integrons (in pdf)

    :param replicon: the replicon to analyse.
    :type replicon: a :class:`Bio.SeqRecord` object.
    :param config: The configuration
    :type config: a :class:`integron_finder.config.Config` object.
    :returns: the path to the integron file (<replicon_id>.integrons)
              and the summary file (<replicon_id.summary>).
              if there is no integron the summary file is None
    :rtype: tuple (str integron_file, str summary_file) or (str integron_file, None)
    """
    result_tmp_dir = config.tmp_dir(replicon.id)
    try:
        os.mkdir(result_tmp_dir)
    except OSError:
        pass
    tmp_replicon_path = os.path.join(result_tmp_dir, replicon.id + '.fst')
    SeqIO.write(replicon, tmp_replicon_path, "fasta")
    # create attr path
    # used to generate protein file with prodigal
    replicon.path = tmp_replicon_path

    # func_annot_path is the canonical path for Functional_annotation
    # path_func_annot is the path provide on the command line
    if config.func_annot and not config.no_proteins and not config.path_func_annot:
        if os.path.exists('bank_hmm'):
            fa_hmm = scan_hmm_bank('bank_hmm')
        elif os.path.exists(config.func_annot_path):
            fa_hmm = scan_hmm_bank(config.func_annot_path)
        else:
            raise IntegronError("Neither the dir '{}' nor 'bank_hmm' exists, specify the location of hmm "
                                "profiles with --path-func-annot option".format(config.func_annot_path))
        is_func_annot = True

    elif config.path_func_annot and config.no_proteins is False:
        fa_hmm = scan_hmm_bank(config.path_func_annot)
        is_func_annot = True
    else:
        is_func_annot = False

    if is_func_annot and not fa_hmm:
        _log.warning("No hmm profiles for functional annotation detected, skip functional annotation step.")

    custom_annot_files = (config.prot_file, config.annot_parser)

    if any(custom_annot_files) and config.gembase:
        raise IntegronError("The --prot-file or --annot-parser are not compatible with --gembase option.")
    elif any(custom_annot_files) and not all(custom_annot_files):
        msg = "If you provide your own proteins file for annotation (--prot-file) " \
              "you have to provide also the parser (--annot-parser)"
        colorlog.critical(msg)
        raise IntegronError(msg)

    if config.gembase_path:
        protein_db = GembaseDB(replicon, config, gembase_path=config.gembase_path)
    elif config.gembase:
        protein_db = GembaseDB(replicon, config)
    elif config.prot_file:
        protein_db = CustomDB(replicon, config, prot_file=config.prot_file)
    else:
        protein_db = ProdigalDB(replicon, config)

    ##################
    # Default search #
    ##################
    intI_file = os.path.join(result_tmp_dir, replicon.id + "_intI.res")
    phageI_file = os.path.join(result_tmp_dir, replicon.id + "_phage_int.res")
    attC_default_file = os.path.join(result_tmp_dir, replicon.id + "_attc_table.res")

    try:
        if not config.no_proteins:
            if not os.path.isfile(intI_file) or not os.path.isfile(phageI_file):
                find_integrase(replicon.id, protein_db.protfile, result_tmp_dir, config)
        _log.info("Starting Default search ... :")
        if not os.path.isfile(attC_default_file):
            # find attc with cmsearch
            find_attc(tmp_replicon_path, replicon.name, config.cmsearch, result_tmp_dir, config.model_attc_path,
                      incE=config.evalue_attc,
                      cpu=config.cpu)

        _log.info("Default search done... : ")
        integrons = find_integron(replicon, protein_db, intI_file, phageI_file, config, attc_file=attC_default_file)

        #########################
        # Search with local_max #
        #########################
        if config.local_max:
            _log.info("Starting search with local_max...:")
            if not os.path.isfile(os.path.join(result_tmp_dir, "integron_max.pickle")):
                circular = True if replicon.topology == 'circ' else False
                integron_max = find_attc_max(integrons, replicon, config.distance_threshold,
                                             config.model_attc_path,
                                             max_attc_size=config.max_attc_size,
                                             min_attc_size=config.min_attc_size,
                                             circular=circular, out_dir=result_tmp_dir,
                                             cpu=config.cpu,
                                             evalue_attc=config.evalue_attc,
                                             cmsearch_bin=config.cmsearch)
                integron_max.to_pickle(os.path.join(result_tmp_dir, "integron_max.pickle"))
                _log.info("Search with local_max done... :")

            else:
                integron_max = pd.read_pickle(os.path.join(result_tmp_dir, "integron_max.pickle"))
                integron_max = integron_max[(integron_max.evalue < config.evalue_attc) &
                                            (abs(integron_max.pos_end - integron_max.pos_beg) < config.max_attc_size) &
                                            (config.min_attc_size < abs(integron_max.pos_end - integron_max.pos_beg))]
                _log.info("Search with local_max was already done, continue... :")

            integrons = find_integron(replicon, protein_db, intI_file, phageI_file, config, attc=integron_max)

        ##########################
        # Add promoters and attI #
        ##########################
        for integron in integrons:
            integron_type = integron.type()
            if integron_type != "In0":  # complete & CALIN
                if not config.no_proteins:
                    _log.info("Adding proteins ... :")
                    integron.add_proteins(protein_db)

            if config.promoter_attI:
                _log.info("Adding promoters and attI ... :")
                if integron_type == "complete":
                    integron.add_promoter()
                    integron.add_attI()
                elif integron_type == "In0":
                    integron.add_attI()
                    integron.add_promoter()
        #########################
        # Functional annotation #
        #########################
        if is_func_annot and fa_hmm:
            _log.info("Starting functional annotation ...:")
            func_annot(integrons, replicon, protein_db, fa_hmm, config, result_tmp_dir)

        #######################
        # Writing out results #
        #######################
        _log.info("Writing out results for replicon {}".format(replicon.id))

        if config.pdf:
            for j, integron in enumerate(integrons, 1):
                if integron.type() == "complete":
                    integron.draw_integron(file=os.path.join(config.result_dir, "{}_{}.pdf".format(replicon.id, j)))

        base_outfile = os.path.join(config.result_dir, replicon.id)
        integron_file = base_outfile + ".integrons"
        _log.debug("Writing integron_file {}".format(integron_file))
        summary_file = base_outfile + ".summary"
        if integrons:
            integrons_report = results.integrons_report(integrons)
            integrons_report.to_csv(integron_file, sep="\t", index=False, na_rep="NA")
            summary = results.summary(integrons_report)
            summary['topology'] = [replicon.topology]
            summary['size'] = [len(replicon)]
            if config.gbk:
                add_feature(replicon, integrons_report, protein_db, config.distance_threshold)
                SeqIO.write(replicon, os.path.join(config.result_dir, replicon.id + ".gbk"), "genbank")
        else:
            with open(integron_file, "w") as out_f:
                out_f.write("# No Integron found\n")
            summary = pd.DataFrame([[replicon.id, 0, 0, 0, replicon.topology, len(replicon)]],
                                   columns=['ID_replicon', 'CALIN', 'complete', 'In0', 'topology', 'size'])
            summary = summary.set_index(['ID_replicon'])
        # summary_file is a temporary file
        # all summaries are merged in results.merge_results function
        summary.to_csv(summary_file, sep="\t", na_rep="NA")

    except integron_finder.EmptyFileError:
        _log.warning('############ Skip replicon {} ############'.format(replicon.name))
        integron_file = ''
        summary_file = ''
    #########################
    # clean temporary files #
    #########################

    if not config.keep_tmp:
        try:
            shutil.rmtree(result_tmp_dir)
        except Exception as err:
            _log.warning("Cannot remove temporary results : '{} : {}'".format(result_tmp_dir, str(err)))
    protein_db.close()
    return integron_file, summary_file


def header(args, hmmsearch, cmsearch, prodigal):
    """

    :param args: the arguments passed on command line. for instance: ['--pdf' '/path/to/replicon']
    :type args: list of strings
    :return: a header containing the name of the program, information about version and licensing.
    :rtype: string
    """
    tpl = r"""
**************************************************************************
 ___       _                               _____ _           _
|_ _|_ __ | |_ ___  __ _ _ __ ___  _ __   |  ___(_)_ __   __| | ___ _ __
 | || '_ \| __/ _ \/ _` | '__/ _ \| '_ \  | |_  | | '_ \ / _` |/ _ \ '__|
 | || | | | ||  __/ (_| | | | (_) | | | | |  _| | | | | | (_| |  __/ |
|___|_| |_|\__\___|\__, |_|  \___/|_| |_| |_|   |_|_| |_|\__,_|\___|_|
                   |___/

**************************************************************************

{version} {commit}

                     =======================

integron_finder is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

integron_finder is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (COPYING file).
If not, see <http://www.gnu.org/licenses/>.

                     =======================

command used: integron_finder {cmd}

                     =======================

"""
    header = tpl.format(version=integron_finder.get_version_message(hmmsearch=hmmsearch,
                                                                    cmsearch=cmsearch,
                                                                    prodigal=prodigal),
                        commit=_if_commit,
                        cmd=' '.join(args)
                        )
    return header


def main(args=None, loglevel=None):
    """
    main entry point to integron_finder

    :param args: the arguments passed on the command line
    :type args: list of str
    :param loglevel: the output verbosity
    :type loglevel: a positive int or a string among 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    """
    global _log

    args = sys.argv[1:] if args is None else args
    config = parse_args(args)

    ###################################
    # Prepare directories for results #
    ###################################

    # need to create directory before to init logger
    # as we write log in integron_finder.out in this dir
    if not os.path.exists(config.outdir):
        os.mkdir(config.outdir)
    else:
        if not os.path.isdir(config.outdir):
            msg = "outdir '{}' already exists and is not a directory".format(config.outdir)
            # _log.critical(msg)
            # we can not log it because logger are not initialized yet.
            raise IsADirectoryError(msg)

    if not os.path.exists(config.result_dir):
        os.mkdir(config.result_dir)
    else:
        if not os.path.isdir(config.result_dir):
            msg = "result dir '{}' already exists and is not a directory".format(config.outdir)
            # _log.critical(msg)
            # we can not log it because logger are not initialized yet.
            raise IsADirectoryError(msg)
        elif not os.access(config.result_dir, os.W_OK):
            msg = "result dir '{}' already exists and is not writable".format(config.outdir)
            # _log.critical(msg)
            # we can not log it because logger are not initialized yet.
            raise PermissionError(msg)

    ####################
    # init the loggers #
    ####################
    log_file = os.path.join(config.result_dir, 'integron_finder.out')
    integron_finder.init_logger(log_file=log_file,
                                out=not config.mute)

    _log = colorlog.getLogger('integron_finder')

    if not loglevel:
        # logs are specify from args options
        logger_set_level(config.log_level)
    else:
        # used by unit tests to mute or unmute logs
        logger_set_level(loglevel)

    #######################################
    # do last config check before running #
    #######################################
    if config.cmsearch is None:
        msg = """cannot find 'cmsearch' in PATH.
Please install infernal package or setup 'cmsearch' binary path with --cmsearch option"""
        _log.critical(msg)
        raise RuntimeError(msg)

    if config.hmmsearch is None:
        msg = """cannot find 'hmmsearch' in PATH.
Please install hmmer package or setup 'hmmsearch' binary path with --hmmsearch option"""
        _log.critical(msg)
        raise RuntimeError(msg)

    if config.prodigal is None:
        msg = """cannot find 'prodigal' in PATH.
Please install prodigal package or setup 'prodigal' binary path with --prodigal option"""
        _log.critical(msg)
        raise RuntimeError(msg)

    ################
    # print Header #
    ################
    log_header = colorlog.getLogger('integron_finder.header')
    logging = integron_finder.get_logging_module()

    handlers = []
    header_log_file = logging.FileHandler(log_file)
    handlers.append(header_log_file)
    if not config.mute:
        header_stream = colorlog.StreamHandler(sys.stdout)
        handlers.append(header_stream)
    formatter = colorlog.ColoredFormatter("%(message)s")
    for h in handlers:
        h.setFormatter(formatter)
        log_header.addHandler(h)
    log_header.setLevel(logging.INFO)
    log_header.propagate = False
    log_header.info(header(args,
                           hmmsearch=config.hmmsearch,
                           cmsearch=config.cmsearch,
                           prodigal=config.prodigal)
                    )

    with utils.FastaIterator(config.input_seq_path, dist_threshold=config.distance_threshold) as sequences_db:
        ################
        # set topology #
        ################


        # the both options are mutually exclusive
        cmd_topo = None
        if config.linear:
            cmd_topo = 'lin'
        elif config.circular:
            cmd_topo = 'circ'
        topologies = Topology(len(sequences_db), cmd_topo,
                              gembase=config.gembase,
                              topology_file=config.topology_file)

        # allow sequences_db to inject topology information
        # in seq.topology attribute
        sequences_db.topologies = topologies

        ##############
        # do the job #
        ##############
        sequences_db_len = len(sequences_db)
        all_integrons = []
        all_summaries = []
        for rep_no, replicon in enumerate(sequences_db, 1):
            # if replicon contains illegal characters
            # or replicon is too short < 50 bp
            # then replicon is None
            if replicon is not None:
                _log.info("############ Processing replicon {} ({}/{}) ############\n".format(replicon.id,
                                                                                              rep_no,
                                                                                              sequences_db_len))
                integron_res, summary = find_integron_in_one_replicon(replicon, config)
                if integron_res:
                    all_integrons.append(integron_res)
                if summary:
                    all_summaries.append(summary)
            else:
                _log.warning("############ Skipping replicon {}/{} ############".format(rep_no,
                                                                                        sequences_db_len))
    if not config.split_results:
        _log.info("Merging integrons results.\n")
        agg_integrons = results.merge_results(*all_integrons)
        agg_summary = results.merge_results(*all_summaries)
        outfile_base_name = os.path.join(config.result_dir, utils.get_name_from_path(config.input_seq_path))
        merged_integron_path = outfile_base_name + ".integrons"
        if not agg_integrons.empty:
            with open(merged_integron_path, 'w') as merged_integron_file:
                merged_integron_file.write(f"# integron_finder {_if_version} {_if_commit}\n")
                merged_integron_file.write(f"# cmd: integron_finder {' '.join(args)}\n")
                agg_integrons.to_csv(merged_integron_file, sep="\t", index=False, na_rep="NA")
        else:
            with open(merged_integron_path, "w") as merged_integron_file:
                merged_integron_file.write(f"# integron_finder {_if_version} {_if_commit}\n")
                merged_integron_file.write(f"# cmd: integron_finder {' '.join(args)}\n")
                merged_integron_file.write("# No Integron found\n")

        merged_summary_path = outfile_base_name + ".summary"
        with open(merged_summary_path, 'w') as merged_summary_file:
            merged_summary_file.write(f"# integron_finder {_if_version} {_if_commit}\n")
            merged_summary_file.write(f"# cmd: integron_finder {' '.join(args)}\n")
            agg_summary.to_csv(merged_summary_file, sep="\t")

        for _file in all_integrons + all_summaries:
            if _file != merged_integron_path and _file != merged_summary_path:
                # in special case where the merged file has the same name that a replicon result file
                os.unlink(_file)


if __name__ == "__main__":
    main()
