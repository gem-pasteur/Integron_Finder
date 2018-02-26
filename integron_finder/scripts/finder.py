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
# Copyright © 2015 - 2016  Institut Pasteur, Paris.                                #
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
import distutils.spawn

from matplotlib import use as m_use
m_use("Agg")
import pandas as pd

import integron_finder
if not integron_finder.__version__.endswith('VERSION'):
    # display warning only for non installed integron_finder
    from Bio import BiopythonExperimentalWarning
    import warnings
    warnings.simplefilter('ignore', FutureWarning)
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)

from Bio import SeqIO

from integron_finder import IntegronError
from integron_finder import utils
from integron_finder.topology import Topology
from integron_finder.config import Config
from integron_finder.hmm import scan_hmm_bank
from integron_finder.integrase import find_integrase
from integron_finder.attc import find_attc, find_attc_max
from integron_finder.integron import find_integron
from integron_finder.annotation import func_annot
from integron_finder.genbank import to_gbk


def parse_args(args):

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("replicon",
                        help="Path to the replicon file (in fasta format), eg : path/to/file.fst or file.fst")

    parser.add_argument("--local_max",
                        help="Allows thorough local detection "
                             "(slower but more sensitive and do not increase false positive rate).",
                        action="store_true")

    parser.add_argument("--func_annot",
                        help="Functional annotation of CDS associated "
                             "with integrons HMM files are needed in Func_annot folder.",
                        default=False,
                        action="store_true")

    parser.add_argument('--cpu',
                        default='1',
                        type=int,
                        help='Number of CPUs used by INFERNAL and HMMER')

    parser.add_argument('-dt', '--distance_thresh',
                        dest='distance_threshold',
                        default=4000,
                        action='store',
                        type=int,
                        help='Two elements are aggregated if they are distant of DISTANCE_THRESH [4kb] or less')

    parser.add_argument('--outdir',
                        default=".",
                        help='Set the output directory (default: current)')

    parser.add_argument("--union_integrases",
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

    parser.add_argument('--path_func_annot',
                        metavar='bank_hmm',
                        help='Path to file containing all hmm bank paths (one per line)')

    parser.add_argument("--gembase",
                        help="Use gembase formatted protein file instead of Prodigal."
                             " Folder structure must be preserved",
                        action="store_true")

    parser.add_argument('--attc_model',
                        default='attc_4.cm',
                        help='path or file to the attc model (Covariance Matrix)')

    parser.add_argument('--evalue_attc',
                        default=1.,
                        type=float,
                        help='set evalue threshold to filter out hits above it (default: 1)')

    parser.add_argument("--keep_palindromes",
                        default=False,
                        help="for a given hit, if the palindromic version is found,"
                             " don't remove the one with highest evalue ",
                        action="store_true")

    parser.add_argument("--no_proteins",
                        help="Don't annotate CDS and don't find integrase, just look for attC sites.",
                        default=True,
                        action="store_true")

    parser.add_argument('--max_attc_size',
                        default=200,
                        type=int,
                        help='set maximum value fot the attC size (default: 200bp)')

    parser.add_argument('--min_attc_size',
                        default=40,
                        type=int,
                        help='set minimum value fot the attC size (default: 40bp)')

    parser.add_argument("--eagle_eyes",
                        help="Synonym of --local_max. Like a soaring eagle in the sky,"
                             " catching rabbits(or attC sites) by surprise.",
                        action="store_true")

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

    # eagle_eyes is just an alias to local_max in whole program use local_max
    args.local_max = args.local_max or args.eagle_eyes

    args = parser.parse_args(args)

    return Config(args)


def find_integron_in_one_replicon(replicon, topology, config):
    in_dir, sequence_file = os.path.split(config.replicon_path)
    replicon_name = replicon.name

    ###################################
    # Prepare directories for results #
    ###################################
    if not os.path.exists(config.outdir):
        os.mkdir(config.outdir)
    else:
        if not os.path.isdir(config.outdir):
            raise IOError("outdir '{}' already exists and is not a directory".format(config.outdir))

    result_dir = os.path.join(config.outdir, "Results_Integron_Finder_" + replicon_name)
    try:
        os.mkdir(result_dir)
    except OSError:
        pass

    result_dir_other = os.path.join(config.outdir, "Results_Integron_Finder_" + replicon_name, "other")
    try:
        os.mkdir(result_dir_other)
    except OSError:
        pass

    # If sequence is too small, it can be problematic when using circularity
    if topology == 'circ' and len(replicon) > 4 * config.distance_threshold:
        circular = True
    else:
        circular = False

    if config.func_annot and not config.no_proteins:
        if os.path.exists('bank_hmm'):
            fa_hmm = scan_hmm_bank('bank_hmm')
        elif os.path.exists(config.func_annot_path):
            fa_hmm = scan_hmm_bank(config.func_annot_path)
        else:
            raise IntegronError("the dir '{}' neither 'bank_hmm' exists, specify the location of hmm "
                                "profile with --path_func_annot option".format(config.func_annot_path))
        is_func_annot = True

    elif config.path_func_annot and config.no_proteins is False:
        fa_hmm = scan_hmm_bank(config.path_func_annot)
        is_func_annot = True
    else:
        is_func_annot = False

    if is_func_annot and not fa_hmm:
        print >> sys.stderr, "WARNING: No hmm profiles for functional annotation detected, " \
                             "skip functional annotation step."

    model_attc_name = config.model_attc.name
    model_len = config.model_len

    max_attc_size = config.max_attc_size
    min_attc_size = config.min_attc_size

    if config.gembase:
        prot_file = os.path.join(in_dir, "..", "Proteins", replicon_name + ".prt")
    else:
        prot_file = os.path.join(result_dir_other, replicon_name + ".prt")

    ##################
    # Default search #
    ##################
    intI_file = os.path.join(result_dir_other, replicon_name + "_intI.res")
    phageI_file = os.path.join(result_dir_other, replicon_name + "_phage_int.res")
    attC_default_file = os.path.join(result_dir_other, replicon_name + "_attc_table.res")

    if not config.no_proteins:
        if not os.path.isfile(intI_file) or not os.path.isfile(phageI_file):
            find_integrase(config.replicon_path, replicon, prot_file, result_dir_other, config)

    print "\n>>> Starting Default search ... :"
    if not os.path.isfile(attC_default_file):
        find_attc(integrons, replicon, config.distance_threshold, config.model_attc_path, config.model_attc_size,
                  circular=circular, out_dir=result_dir_other)

    print ">>> Default search done... : \n"
    integrons = find_integron(replicon, attC_default_file, intI_file, phageI_file, config)

    #########################
    # Search with local_max #
    #########################
    if config.local_max:

        print "\n>>>>>> Starting search with local_max...:"
        if not os.path.isfile(os.path.join(result_dir_other, "integron_max.pickle")):

            integron_max = find_attc_max(integrons, circular=circular)
            integron_max.to_pickle(os.path.join(result_dir_other, "integron_max.pickle"))
            print ">>>>>> Search with local_max done... : \n"

        else:
            integron_max = pd.read_pickle(os.path.join(result_dir_other, "integron_max.pickle"))
            print ">>>>>> Search with local_max was already done, continue... : \n"

        integrons = find_integron(replicon, attC_default_file, intI_file, phageI_file, config)

    ##########################
    # Add promoters and attI #
    ##########################

    outfile = replicon_name + ".integrons"
    for integron in integrons:
        if integron.type() != "In0":  # complete & CALIN
            if not config.no_proteins:
                integron.add_proteins()
        elif integron.type() == "complete":
            integron.add_promoter()
            integron.add_attI()
        elif integron.type() == "In0":
            integron.add_attI()
            integron.add_promoter()

    #########################
    # Functional annotation #
    #########################
    if is_func_annot and len(fa_hmm) > 0:
        func_annot(integrons, replicon, prot_file, fa_hmm, result_dir_other, config)

    for j, integron in enumerate(integrons, 1):
        if integron.type() == "complete":
            integron.draw_integron(file=os.path.join(result_dir, "{}_{}.pd".format(replicon.name, j)))

    #######################
    # Writing out results #
    #######################

    integrons_describe = pd.concat([i.describe() for i in integrons])
    dic_id = {id_: "{:02}".format(j) for j, id_ in
              enumerate(integrons_describe.sort_values("pos_beg").id_integron.unique(), 1)}
    integrons_describe.id_integron = ["integron_" + dic_id[id_] for id_ in integrons_describe.id_integron]
    integrons_describe = integrons_describe[["ID_integron", "ID_replicon", "element",
                                             "pos_beg", "pos_end", "strand", "evalue",
                                             "type_elt", "annotation", "model",
                                             "type", "default", "distance_2attC"]]
    integrons_describe['evalue'] = integrons_describe.evalue.astype(float)
    integrons_describe.index = range(len(integrons_describe))

    integrons_describe.sort_values(["ID_integron", "pos_beg", "evalue"], inplace=True)

    integrons_describe.to_csv(os.path.join(result_dir, outfile), sep="\t", index=0, na_rep="NA")
    to_gbk(integrons_describe, replicon)
    SeqIO.write(replicon, os.path.join(result_dir, replicon_name + ".gbk"), "genbank")

    if not integrons:
        with open(os.path.join(result_dir, outfile), "w") as out_f:
            out_f.write("# No Integron found\n")


def main(args=None):
    args = sys.argv[1:] if args is None else args

    _prefix_share = '$PREFIXSHARE'

    # integron was not installed using the setup.py
    # it's a development version using environment variable
    if 'INTEGRON_HOME' in os.environ and os.environ['INTEGRON_HOME']:
        _prefix_share = os.environ['INTEGRON_HOME']
    _prefix_data = os.path.join(_prefix_share, 'data')

    config = parse_args(args)

    if not os.path.exists(_prefix_data):
        raise Exception("""cannot find integron_finder data check your installation
or define INTEGRON_HOME environment variable.""")

    if config.cmsearch is None:
        raise RuntimeError("""cannot find 'cmsearch' in PATH.
Please install infernal package or setup 'cmsearch' binary path with --cmsearch option""")

    if config.hmmsearch is None:
        raise RuntimeError("""cannot find 'hmmsearch' in PATH.
Please install hmmer package or setup 'hmmsearch' binary path with --hmmsearch option""")

    if config.prodigal is None:
        raise RuntimeError("""cannot find 'prodigal' in PATH.
Please install prodigal package or setup 'prodigal' binary path with --prodigal option""")

    sequences_db = utils.read_multi_dna_fasta(config.replicon_path)

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

    ##############
    # do the job #
    ##############
    for replicon in sequences_db:
        topology = topologies[replicon.id]
        find_integron_in_one_replicon(replicon, topology, config)


if __name__ == "__main__":

    main()
