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
from Bio import Seq
from integron_finder import IntegronError

from integron_finder.topology import Topology
from integron_finder.config import Config


def main(args=None):
    args = sys.argv[1:] if args is None else args

    _prefix_share = '$PREFIXSHARE'

    # integron was not installed using the setup.py
    # it's a development version using environment variable
    if 'INTEGRON_HOME' in os.environ and os.environ['INTEGRON_HOME']:
        _prefix_share = os.environ['INTEGRON_HOME']
    elif _prefix_share.endswith('PREFIXSHARE'):
        _prefix_share = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    _prefix_data = os.path.join(_prefix_share, 'data')

    if not os.path.exists(_prefix_data):
        raise Exception("""cannot find integron_finder data check your installation
    or define INTEGRON_HOME environment variable.""")

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
                        action='store',
                        type=str,
                        help='Number of CPUs used by INFERNAL and HMMER')

    parser.add_argument('-dt', '--distance_thresh',
                        dest='distance_threshold',
                        default=4000,
                        action='store',
                        type=int,
                        help='Two elements are aggregated if they are distant of DISTANCE_THRESH [4kb] or less')

    parser.add_argument('--outdir',
                        default=".",
                        action='store',
                        type=str,
                        metavar='.',
                        help='Set the output directory (default: current)')

    parser.add_argument("--union_integrases",
                        help="Instead of taking intersection of hits from Phage_int profile (Tyr recombinases)"
                             " and integron_integrase profile, use the union of the hits",
                        action="store_true")

    parser.add_argument('--cmsearch',
                        default=distutils.spawn.find_executable("cmsearch"),
                        action='store',
                        type=str,
                        help='Complete path to cmsearch if not in PATH. eg: /usr/local/bin/cmsearch')

    parser.add_argument('--hmmsearch',
                        default=distutils.spawn.find_executable("hmmsearch"),
                        action='store',
                        type=str,
                        help='Complete path to hmmsearch if not in PATH. eg: /usr/local/bin/hmmsearch')

    parser.add_argument('--prodigal',
                        default=distutils.spawn.find_executable("prodigal"),
                        action='store',
                        type=str,
                        help='Complete path to prodigal if not in PATH. eg: /usr/local/bin/prodigal')

    parser.add_argument('--path_func_annot',
                        action='store',
                        metavar='bank_hmm',
                        type=str,
                        help='Path to file containing all hmm bank paths (one per line)')

    parser.add_argument("--gembase",
                        help="Use gembase formatted protein file instead of Prodigal."
                             " Folder structure must be preserved",
                        action="store_true")

    parser.add_argument('--attc_model',
                        default='attc_4.cm',
                        action='store',
                        type=str,
                        metavar='file.cm',
                        help='path or file to the attc model (Covariance Matrix)')

    parser.add_argument('--evalue_attc',
                        default=1.,
                        action='store',
                        type=float,
                        metavar='1',
                        help='set evalue threshold to filter out hits above it (default: 1)')

    parser.add_argument("--keep_palindromes",
                        help="for a given hit, if the palindromic version is found,"
                             " don't remove the one with highest evalue ",
                        action="store_true")

    parser.add_argument("--no_proteins",
                        help="Don't annotate CDS and don't find integrase, just look for attC sites.",
                        default=True,
                        action="store_true")

    parser.add_argument('--max_attc_size',
                        default=200,
                        action='store',
                        type=int,
                        metavar='200',
                        help='set maximum value fot the attC size (default: 200bp)')

    parser.add_argument('--min_attc_size',
                        default=40,
                        action='store',
                        type=int,
                        metavar='40',
                        help='set minimum value fot the attC size (default: 40bp)')

    parser.add_argument("--eagle_eyes",
                        help="Synonym of --local_max. Like a soaring eagle in the sky,"
                             " catching rabbits(or attC sites) by surprise.",
                        action="store_true")

    topology_grp = parser.add_mutually_exclusive_group()
    topology_grp.add_argument("--circ",
                              dest='circular',
                              help="Set the default topology for replicons to 'cirular'",
                              action="store_true")
    topology_grp.add_argument("--linear",
                              help="Set the default topology for replicons to 'linear'",
                              action="store_true")
    parser.add_argument("--topology-file",
                        help="The path to a file where the topology for each replicon is specified")

    parser.add_argument("-V", "--version",
                        action="version",
                        version=integron_finder.get_version_message())

    args = parser.parse_args(args)

    config = Config(args)

    in_dir, sequence_file = os.path.split(config.replicon_path)
    replicon_name = config.replicon_name

    if not os.path.exists(config.outdir):
        os.mkdir(config.outdir)
    else:
        if not os.path.isdir(config.outdir):
            raise IOError("outdir '{}' is not a directory".format(config.outdir))
    try:
        os.mkdir(os.path.join(args.outdir, "Results_Integron_Finder_" + replicon_name))
    except OSError:
        pass

    try:
        os.mkdir(os.path.join(args.outdir,
                              "Results_Integron_Finder_" + replicon_name,
                              "other"))
    except OSError:
        pass

    result_dir_other = os.path.join(config.outdir,
                                    "Results_Integron_Finder_" + replicon_name,
                                    "other")
    result_dir = os.path.join(args.outdir,
                              "Results_Integron_Finder_" + replicon_name)

    sequences_db = SeqIO.index(config.replicon_path,
                               "fasta",
                               alphabet=Seq.IUPAC.unambiguous_dna)

    ################
    # set topology #
    ################
    default_topology = None

    if default_topology is None:
        if len(sequences_db) == 1:
            default_topology = 'circ'
        else:
            default_topology = 'lin'

    topologies = Topology(default_topology, topology_file=args.topology_file)

    ###############
    # Definitions #
    ###############
    for seq_id in sequences_db:
        sequence = sequences_db[seq_id]
        replicon_size = len(sequence)

        # If sequence is too small, it can be problematic when using circularity
        topology = topologies[seq_id]
        if topology == 'circ' and replicon_size > 4 * config.distance_threshold:
            circular = True
        else:
            circular = False

        if config.cmsearch is None:
            raise RuntimeError("""cannot find 'cmsearch' in PATH.
Please install infernal package or setup 'cmsearch' binary path with --cmsearch option""")

        if config.hmmsearch is None:
            raise RuntimeError("""cannot find 'hmmsearch' in PATH.
Please install hmmer package or setup 'hmmsearch' binary path with --hmmsearch option""")

        if config.prodigal is None:
            raise RuntimeError("""cannot find 'prodigal' in PATH.
Please install prodigal package or setup 'prodigal' binary path with --prodigal option""")

        if config.func_annot and not config.no_proteins:
            if os.path.exists('bank_hmm'):
                fa_hmm = integron_finder.hmm.scan_hmm_bank('bank_hmm')
            elif os.path.exists(config.func_annot_path):
                fa_hmm = integron_finder.hmm.scan_hmm_bank(config.func_annot_path)
            else:
                raise IntegronError("the dir '{}' neither '{}' exists, specify the location of hmm \
                    profile with --path_func_annot option".format(config.func_annot_path, 'bank_hmm'))
            is_func_annot = True

        elif args.path_func_annot and args.no_proteins is False:
            fa_hmm = integron_finder.hmm.scan_hmm_bank(args.path_func_annot)
            is_func_annot = True
        else:
            is_func_annot = False

        if is_func_annot and not fa_hmm:
            print >> sys.stderr, "WARNING: No hmm profiles for functional annotation detected, \
                skip functional annotation step."

        model_attc_name = config.model_attc.split("/")[-1].split(".cm")[0]

        max_attc_size = args.max_attc_size
        min_attc_size = args.min_attc_size

        with open(config.model_attc) as f:
            for i, line in enumerate(f):
                if "CLEN" in line:
                    length_cm = int(line.split()[1])
                    break

        if config.gembase:
            prot_dir = os.path.join(in_dir, "..", "Proteins")
            prot_file = os.path.join(prot_dir, replicon_name + ".prt")
        else:
            prot_file = os.path.join(result_dir_other, replicon_name + ".prt")

        ##################
        # Default search #
        ##################
        intI_file = os.path.join(result_dir_other, replicon_name + "_intI.res")
        phageI_file = os.path.join(result_dir_other, replicon_name + "_phage_int.res")
        attC_default_file = os.path.join(result_dir_other, replicon_name + "_attc_table.res")

        if not config.no_proteins:
            if (os.path.isfile(intI_file) == 0 or
                    os.path.isfile(phageI_file) == 0):
                integron_finder.integrase.find_integrase(config.replicon_path, replicon_name, result_dir_other)

        print "\n>>> Starting Default search ... :"
        if os.path.isfile(attC_default_file) == 0:
            integron_finder.attc.find_attc(config.replicon_path, replicon_name, result_dir_other)

        print ">>> Default search done... : \n"
        integrons = integron_finder.integron.find_integron(replicon_name,
                                                    attC_default_file,
                                                    intI_file,
                                                    phageI_file)

        #########################
        # Search with local_max #
        #########################
        if config.eagle_eyes or config.local_max:

            print "\n>>>>>> Starting search with local_max...:"
            if os.path.isfile(os.path.join(result_dir_other, "integron_max.pickle")) == 0:

                integron_max = integron_finder.attc.find_attc_max(integrons, circular=circular)
                integron_max.to_pickle(os.path.join(result_dir_other, "integron_max.pickle"))
                print ">>>>>> Search with local_max done... : \n"

            else:
                integron_max = pd.read_pickle(os.path.join(result_dir_other,
                                                           "integron_max.pickle"))
                print ">>>>>> Search with local_max was already done, continue... : \n"

            integrons = integron_finder.integrase.find_integron(replicon_name,
                                                         integron_max,
                                                         intI_file,
                                                         phageI_file)

        ##########################
        # Add promoters and attI #
        ##########################

        outfile = replicon_name + ".integrons"

        if len(integrons):
            j = 1
            for i in integrons:
                if i.type() != "In0":  # complete & CALIN
                    if not config.no_proteins:
                        i.add_proteins()

                if i.type() == "complete":
                    i.add_promoter()
                    i.add_attI()
                    j += 1
                if i.type() == "In0":
                    i.add_attI()
                    i.add_promoter()

            #########################
            # Functional annotation #
            #########################
            if is_func_annot and len(fa_hmm) > 0:
                integron_finder.annotation.func_annot(replicon_name, result_dir_other, fa_hmm)

            j = 1
            for i in integrons:
                if i.type() == "complete":
                    i.draw_integron(file=os.path.join(result_dir, "{}_{}.pd".format(replicon_name, j)))
                    j += 1

            #######################
            # Writing out results #
            #######################

            integrons_describe = pd.concat([i.describe() for i in integrons])
            dic_id = {i: "%02i" % (j + 1) for j, i in
                      enumerate(integrons_describe.sort_values("pos_beg").id_integron.unique())}
            integrons_describe.id_integron = ["integron_" + dic_id[i] for i in integrons_describe.id_integron]
            integrons_describe = integrons_describe[["ID_integron", "ID_replicon", "element",
                                                     "pos_beg", "pos_end", "strand", "evalue",
                                                     "type_elt", "annotation", "model",
                                                     "type", "default", "distance_2attC"]]
            integrons_describe['evalue'] = integrons_describe.evalue.astype(float)
            integrons_describe.index = range(len(integrons_describe))

            integrons_describe.sort_values(["ID_integron", "pos_beg", "evalue"], inplace=True)

            integrons_describe.to_csv(os.path.join(result_dir, outfile), sep="\t", index=0, na_rep="NA")
            integron_finder.genbank.to_gbk(integrons_describe, sequence)
            SeqIO.write(sequence, os.path.join(result_dir, replicon_name + ".gbk"), "genbank")
        else:
            out_f = open(os.path.join(result_dir, outfile), "w")
            out_f.write("# No Integron found\n")
            out_f.close()


if __name__ == "__main__":

    main()
