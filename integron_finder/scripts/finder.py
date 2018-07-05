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
# Copyright (c) 2015 - 2018  Institut Pasteur, Paris.                              #
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

from __future__ import print_function

"""
integron_finder is a program that looks for integron in DNA sequences.
"""
import os
import sys
import argparse
import distutils.spawn
import shutil

import pandas as pd

import integron_finder
if not integron_finder.__version__.endswith('VERSION'):
    # display warning only for non installed integron_finder
    from Bio import BiopythonExperimentalWarning
    import warnings
    warnings.simplefilter('ignore', FutureWarning)
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)

# must be done after import 'integron_finder'
import colorlog
_log = colorlog.getLogger('integron_finder')

from Bio import SeqIO

from integron_finder import IntegronError, logger_set_level
from integron_finder import utils
from integron_finder import results
from integron_finder.topology import Topology
from integron_finder.config import Config
from integron_finder.hmm import scan_hmm_bank
from integron_finder.integrase import find_integrase
from integron_finder.attc import find_attc, find_attc_max
from integron_finder.integron import find_integron
from integron_finder.annotation import func_annot, add_feature


def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("replicon",
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
                        help='Number of CPUs used by INFERNAL and HMMER')

    parser.add_argument('-dt', '--distance-thresh',
                        dest='distance_threshold',
                        default=4000,
                        action='store',
                        type=int,
                        help='Two elements are aggregated if they are distant of DISTANCE_THRESH [4kb] or less')

    parser.add_argument('--outdir',
                        default=".",
                        help='Set the output directory (default: current)')

    parser.add_argument("--union-integrases",
                        default=False,
                        help="Instead of taking intersection of hits from Phage_int profile (Tyr recombinases)"
                             " and integron_integrase profile, use the union of the hits",
                        action="store_true")

    parser.add_argument('--cmsearch',
                        default=distutils.spawn.find_executable("cmsearch"),
                        type=str,
                        help='Complete path to cmsearch if not in PATH. eg: /usr/local/bin/cmsearch')

    parser.add_argument('--hmmsearch',
                        default=distutils.spawn.find_executable("hmmsearch"),
                        help='Complete path to hmmsearch if not in PATH. eg: /usr/local/bin/hmmsearch')

    parser.add_argument('--prodigal',
                        default=distutils.spawn.find_executable("prodigal"),
                        help='Complete path to prodigal if not in PATH. eg: /usr/local/bin/prodigal')

    parser.add_argument('--path-func-annot',
                        help='Path to file containing all hmm bank paths (one per line)')

    parser.add_argument("--gembase",
                        default=False,
                        help="Use gembase formatted protein file instead of Prodigal."
                             " Folder structure must be preserved",
                        action="store_true")

    parser.add_argument('--attc-model',
                        default='attc_4.cm',
                        help='path or file to the attc model (Covariance Matrix)')

    parser.add_argument('--evalue-attc',
                        default=1.,
                        type=float,
                        help='set evalue threshold to filter out hits above it (default: 1)')

    parser.add_argument('--calin-threshold',
                        default=2,
                        type=int,
                        help="keep 'CALIN' only if attC sites nuber >= calin-threshold (default: 2")

    parser.add_argument("--keep-palindromes",
                        default=False,
                        help="for a given hit, if the palindromic version is found,"
                             " don't remove the one with highest evalue ",
                        action="store_true")

    parser.add_argument("--no-proteins",
                        help="Don't annotate CDS and don't find integrase, just look for attC sites.",
                        default=False,
                        action="store_true")

    parser.add_argument('--max-attc-size',
                        default=200,
                        type=int,
                        help='set maximum value fot the attC size (default: 200bp)')

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
                                     'This results are stored in directory named other_<replicon id>')
    output_options.add_argument('--split-results',
                                action='store_true',
                                default=False,
                                help='Instead of merging integron results from all replicon in one file, '
                                     'keep them in separated files.')

    topology_grp = parser.add_mutually_exclusive_group()
    topology_grp.add_argument("--circ",
                              dest='circular',
                              default=False,
                              help="Set the default topology for replicons to 'cirular'",
                              action="store_true")
    topology_grp.add_argument("--linear",
                              default=False,
                              help="Set the default topology for replicons to 'linear'",
                              action="store_true")
    parser.add_argument("--topology-file",
                        help="The path to a file where the topology for each replicon is specified")

    parser.add_argument("-V", "--version",
                        action="version",
                        version=integron_finder.get_version_message())

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
    in_dir, sequence_file = os.path.split(config.replicon_path)
    result_dir_other = os.path.join(config.result_dir, "other_{}".format(replicon.id))
    try:
        os.mkdir(result_dir_other)
    except OSError:
        pass
    tmp_replicon_path = os.path.join(result_dir_other, replicon.id + '.fst')
    SeqIO.write(replicon, tmp_replicon_path, "fasta")

    if config.func_annot and not config.no_proteins and not config.path_func_annot:
        if os.path.exists('bank_hmm'):
            fa_hmm = scan_hmm_bank('bank_hmm')
        elif os.path.exists(config.func_annot_path):
            fa_hmm = scan_hmm_bank(config.func_annot_path)
        else:
            raise IntegronError("the dir '{}' neither 'bank_hmm' exists, specify the location of hmm "
                                "profile with --path-func-annot option".format(config.func_annot_path))
        is_func_annot = True

    elif config.path_func_annot and config.no_proteins is False:
        fa_hmm = scan_hmm_bank(config.path_func_annot)
        is_func_annot = True
    else:
        is_func_annot = False

    if is_func_annot and not fa_hmm:
        _log.warning("No hmm profiles for functional annotation detected, skip functional annotation step.")

    if config.gembase:
        replicon_name = utils.get_name_from_path(config.replicon_path)
        prot_file = os.path.join(in_dir, "..", "Proteins", replicon_name + ".prt")
    else:
        prot_file = os.path.join(result_dir_other, replicon.id + ".prt")

    ##################
    # Default search #
    ##################
    intI_file = os.path.join(result_dir_other, replicon.id + "_intI.res")
    phageI_file = os.path.join(result_dir_other, replicon.id + "_phage_int.res")
    attC_default_file = os.path.join(result_dir_other, replicon.id + "_attc_table.res")

    if not config.no_proteins:
        if not os.path.isfile(intI_file) or not os.path.isfile(phageI_file):
            find_integrase(tmp_replicon_path, replicon, prot_file, result_dir_other, config)

    _log.info("Starting Default search ... :")
    if not os.path.isfile(attC_default_file):
        find_attc(tmp_replicon_path, replicon.name, config.cmsearch, result_dir_other, config.model_attc_path,
                  cpu=config.cpu)

    _log.info("Default search done... : ")
    integrons = find_integron(replicon, attC_default_file, intI_file, phageI_file, config)

    #########################
    # Search with local_max #
    #########################
    if config.local_max:
        _log.info("Starting search with local_max...:")
        if not os.path.isfile(os.path.join(result_dir_other, "integron_max.pickle")):
            circular = True if replicon.topology == 'circ' else False
            integron_max = find_attc_max(integrons, replicon, config.distance_threshold,
                                         config.model_attc_path, config.max_attc_size,
                                         circular=circular, out_dir=result_dir_other,
                                         cpu=config.cpu)
            integron_max.to_pickle(os.path.join(result_dir_other, "integron_max.pickle"))
            _log.info("Search with local_max done... :")

        else:
            integron_max = pd.read_pickle(os.path.join(result_dir_other, "integron_max.pickle"))
            _log.info("Search with local_max was already done, continue... :")

        integrons = find_integron(replicon, integron_max, intI_file, phageI_file, config)

    ##########################
    # Add promoters and attI #
    ##########################
    _log.info("Adding promoters and attI ... :")
    for integron in integrons:
        integron_type = integron.type()
        if integron_type != "In0":  # complete & CALIN
            if not config.no_proteins:
                integron.add_proteins(prot_file)

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
        func_annot(integrons, replicon, prot_file, fa_hmm, config, result_dir_other)

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
    if integrons:
        integrons_report = results.integrons_report(integrons)
        integrons_report.to_csv(integron_file, sep="\t", index=False, na_rep="NA")

        summary = results.summary(integrons_report)
        summary_file = base_outfile + ".summary"
        summary.to_csv(summary_file, sep="\t", na_rep="NA", index=False,
                       columns=['ID_replicon', 'ID_integron', 'complete', 'In0', 'CALIN'])
        if config.gbk:
            add_feature(replicon, integrons_report, prot_file, config.distance_threshold)
            SeqIO.write(replicon, os.path.join(config.result_dir, replicon.id + ".gbk"), "genbank")
    else:
        with open(integron_file, "w") as out_f:
            out_f.write("# No Integron found\n")
        summary_file = None

    #########################
    # clean temporary files #
    #########################

    if not config.keep_tmp:
        try:
            shutil.rmtree(result_dir_other)
        except Exception as err:
            _log.warning("Cannot remove temporary results : '{} : {}'".format(result_dir_other, str(err)))

    return integron_file, summary_file



def header(args):
    """

    :param args: the arguments passed on command line. for instance: ['--pdf' '/path/to/replicon']
    :type args: list of strings
    :return: a header containing the name of the program, information about version and licensing.
    :rtype: string
    """
    tpl = """
**************************************************************************    
 ___       _                               _____ _           _           
|_ _|_ __ | |_ ___  __ _ _ __ ___  _ __   |  ___(_)_ __   __| | ___ _ __ 
 | || '_ \| __/ _ \/ _` | '__/ _ \| '_ \  | |_  | | '_ \ / _` |/ _ \ '__|
 | || | | | ||  __/ (_| | | | (_) | | | | |  _| | | | | | (_| |  __/ |   
|___|_| |_|\__\___|\__, |_|  \___/|_| |_| |_|   |_|_| |_|\__,_|\___|_|   
                   |___/                                                 
   
**************************************************************************    

{version}

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
    header = tpl.format(version=integron_finder.get_version_message(),
                        cmd=' '.join(args)
                        )
    return header


def main(args=None, loglevel=None):
    """
    main entry point to integron_finder

    :param str args: the arguments passed on the command line
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
    logging = colorlog.logging.logging
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
    log_header.setLevel(colorlog.logging.logging.INFO)
    log_header.propagate = False
    log_header.info(header(args))

    with utils.FastaIterator(config.replicon_path, dist_threshold=config.distance_threshold) as sequences_db:
        ################
        # set topology #
        ################
        default_topology = 'circ' if len(sequences_db) == 1 else 'lin'
        if config.linear:
            default_topology = 'lin'
        elif config.circular:
            default_topology = 'circ'
        # the both options are mutually exclusive
        topologies = Topology(default_topology, topology_file=config.topology_file)

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
        outfile_base_name = os.path.join(config.result_dir, utils.get_name_from_path(config.replicon_path))
        merged_integron_file = outfile_base_name + ".integrons"
        if not agg_integrons.empty:
            agg_integrons.to_csv(merged_integron_file, sep="\t", index=False, na_rep="NA")
        else:
            with open(merged_integron_file, "w") as out_f:
                out_f.write("# No Integron found\n")
        merged_summary_file = outfile_base_name + ".summary"
        if not agg_integrons.empty:
            agg_summary.to_csv(merged_summary_file, sep="\t", index=False, na_rep="NA",
                               columns=['ID_replicon', 'ID_integron', 'complete', 'In0', 'CALIN'])

        for _file in all_integrons + all_summaries:
            if _file != merged_integron_file and _file != merged_summary_file:
                # in special case where the merged file has the same name that a replicon result file
                os.unlink(_file)


if __name__ == "__main__":

    main()
