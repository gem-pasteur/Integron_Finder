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

import os
import argparse

import pandas as pd
import pandas.testing as pdt

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.config import Config
from integron_finder.hmm import read_hmm
from integron_finder.prot_db import GembaseDB, ProdigalDB
from integron_finder.utils import MultiFastaReader


class TestReadHMM(IntegronTest):

    def setUp(self) -> None:
        self.args = argparse.Namespace()
        self.args.gembase = False
        self.args.prot_file = False
        # need to provide any file to fool config binary check
        self.args.cmsearch = __file__
        self.args.hmmsearch = __file__
        self.args.prodigal = __file__


    def test_read_empty(self):
        # Test that when there are no hits in the hmm result file, it returns an empty
        # dataframe, without error.
        rep_name = "acba.007.p01.13"
        replicon_id = 'ACBA.007.P01_13'

        replicon_path = self.find_data(os.path.join('Replicons', rep_name + '.fst'))
        prot_file = self.find_data(os.path.join('Proteins', replicon_id + '.prt'))

        self.args.replicon = replicon_path
        cfg = Config(self.args)

        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)

        infile = self.find_data(os.path.join("fictive_results", "{}_intI-empty.res".format(replicon_id)))

        df = read_hmm(rep_name, prot_db, infile, cfg)
        exp = pd.DataFrame(columns=["Accession_number", "query_name", "ID_query", "ID_prot",
                                    "strand", "pos_beg", "pos_end", "evalue"])

        intcols = ["pos_beg", "pos_end", "strand"]
        floatcol = ["evalue"]
        exp[intcols] = exp[intcols].astype(int)
        exp[floatcol] = exp[floatcol].astype(float)
        pdt.assert_frame_equal(df, exp)


    def test_read_hmm(self):
        # Test that the hmm hits are well read
        rep_name = "acba.007.p01.13"
        replicon_id = 'ACBA.007.P01_13'

        replicon_path = self.find_data(os.path.join('Replicons', rep_name + '.fst'))
        prot_file = self.find_data(os.path.join('Proteins', replicon_id + '.prt'))

        self.args.replicon = replicon_path
        cfg = Config(self.args)

        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)

        infile = self.find_data(os.path.join("Results_Integron_Finder_{}".format(rep_name),
                                             "tmp_{}".format(replicon_id),
                                             "{}_intI.res".format(replicon_id)))

        df = read_hmm(rep_name, prot_db, infile, cfg)
        exp = pd.DataFrame(data={"Accession_number": rep_name, "query_name": "intI_Cterm",
                                 "ID_query": "-", "ID_prot": "ACBA.007.P01_13_1", "strand": 1,
                                 "pos_beg": 55, "pos_end": 1014, "evalue": 1.9e-25},
                           index=[0])
        exp = exp[["Accession_number", "query_name", "ID_query", "ID_prot",
                   "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df, exp)


    def test_read_hmm_gembase(self):
        # Test that the hmm hits are well read, when the gembase format is used (.prt file is
        # provided, prodigal is not used to find the proteins).
        replicon_id = 'ACBA.0917.00019'
        contig_id = 'ACBA.0917.00019.0001'
        result_dir_expected = self.find_data("Results_Integron_Finder_{}.gembase".format(replicon_id))
        replicon_path = self.find_data('Gembase', 'Gembase1', 'Replicons', replicon_id + '.fna')
        prot_file = os.path.join(result_dir_expected, "tmp_{}".format(contig_id), contig_id + '.prt')
        infile = os.path.join(result_dir_expected, "tmp_{}".format(contig_id), "{}_intI.res".format(contig_id))

        self.args.gembase = True
        self.args.replicon = replicon_path
        cfg = Config(self.args)

        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        with self.catch_log():
            prot_db = GembaseDB(replicon, cfg, prot_file=prot_file)

        df = read_hmm(contig_id, prot_db, infile, cfg)
        exp = pd.DataFrame(data={"Accession_number": contig_id,
                                 "query_name": "intI_Cterm",
                                 "ID_query": "-",
                                 "ID_prot": "ACBA.0917.00019.i0001_00298",
                                 "strand": -1,
                                 "pos_beg": 311597,
                                 "pos_end": 312631,
                                 "evalue": 3.6e-25},
                           index=[0])
        exp = exp[["Accession_number", "query_name", "ID_query", "ID_prot",
                   "strand", "pos_beg", "pos_end", "evalue"]]

        pdt.assert_frame_equal(df, exp)


    def test_read_hmm_evalue(self):
        # Test that the hmm hits are well read, and returned only if evalue is < to the
        # given threshold.
        rep_name = "acba.007.p01.13"
        replicon_id = 'ACBA.007.P01_13'

        replicon_path = self.find_data('Replicons', rep_name + '.fst')
        prot_file = self.find_data('Proteins', replicon_id + '.prt')

        self.args.replicon = replicon_path
        cfg = Config(self.args)

        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)

        infile = self.find_data("Results_Integron_Finder_{}".format(rep_name),
                                "tmp_{}".format(replicon_id),
                                "{}_intI.res".format(replicon_id))

        df1 = read_hmm(rep_name, prot_db, infile, cfg, evalue=1.95e-25)
        exp1 = pd.DataFrame(data={"Accession_number": rep_name, "query_name": "intI_Cterm",
                                  "ID_query": "-", "ID_prot": "ACBA.007.P01_13_1", "strand": 1,
                                  "pos_beg": 55, "pos_end": 1014, "evalue": 1.9e-25},
                            index=[0])
        exp1 = exp1[["Accession_number", "query_name", "ID_query", "ID_prot",
                     "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df1, exp1)

        df2 = read_hmm(replicon_id, prot_db, infile, cfg, evalue=1.9e-25)
        exp2 = pd.DataFrame(columns=["Accession_number", "query_name", "ID_query", "ID_prot",
                                     "strand", "pos_beg", "pos_end", "evalue"])

        intcols = ["pos_beg", "pos_end", "strand"]
        floatcol = ["evalue"]
        exp2[intcols] = exp2[intcols].astype(int)
        exp2[floatcol] = exp2[floatcol].astype(float)
        pdt.assert_frame_equal(df2, exp2)


    def test_read_hmm_evalue2(self):
        # Test that the hmm hits are well read, it returns only the hits with evalue < given
        # threshold
        rep_name = "acba.007.p01.13"
        replicon_id = 'ACBA.007.P01_13'

        replicon_path = self.find_data('Replicons', rep_name + '.fst')
        prot_file = self.find_data('Proteins', replicon_id + '.prt')

        self.args.replicon = replicon_path
        cfg = Config(self.args)

        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)

        infile = self.find_data("fictive_results", "{}_intI.res".format(replicon_id))

        df1 = read_hmm(replicon_id, prot_db, infile, cfg, evalue=1e-3)
        exp1 = pd.DataFrame(data={"Accession_number": [replicon_id] * 2,
                                  "query_name": ["intI_Cterm"] * 2,
                                  "ID_query": ["-", "-"],
                                  "ID_prot": ["ACBA.007.P01_13_1", "ACBA.007.P01_13_3"],
                                  "strand": [1, -1],
                                  "pos_beg": [55, 1722], "pos_end": [1014, 2537],
                                  "evalue": [1.9e-25, 2e-25]},
                            index=[0, 1])
        exp1 = exp1[["Accession_number", "query_name", "ID_query", "ID_prot",
                     "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df1, exp1)


    def test_read_hmm_cov(self):
        # Test that the hmm hits are well read, and returned only if coverage is > to the
        # given threshold.
        rep_name = "acba.007.p01.13"
        replicon_id = 'ACBA.007.P01_13'

        replicon_path = self.find_data('Replicons', rep_name + '.fst')
        prot_file = self.find_data('Proteins', replicon_id + '.prt')

        self.args.replicon = replicon_path
        cfg = Config(self.args)

        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)

        infile = self.find_data("Results_Integron_Finder_{}".format(rep_name),
                                "tmp_{}".format(replicon_id),
                                "{}_intI.res".format(replicon_id))

        df1 = read_hmm(rep_name, prot_db, infile, cfg, coverage=0.945)
        exp1 = pd.DataFrame(data={"Accession_number": rep_name, "query_name": "intI_Cterm",
                                  "ID_query": "-", "ID_prot": "ACBA.007.P01_13_1", "strand": 1,
                                  "pos_beg": 55, "pos_end": 1014, "evalue": 1.9e-25},
                            index=[0])
        exp1 = exp1[["Accession_number", "query_name", "ID_query", "ID_prot",
                     "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df1, exp1)
        df2 = read_hmm(rep_name, prot_db, infile, cfg, coverage=0.95)
        exp2 = pd.DataFrame(columns=["Accession_number", "query_name", "ID_query", "ID_prot",
                                     "strand", "pos_beg", "pos_end", "evalue"])
        intcols = ["pos_beg", "pos_end", "strand"]
        floatcol = ["evalue"]
        exp2[intcols] = exp2[intcols].astype(int)
        exp2[floatcol] = exp2[floatcol].astype(float)
        pdt.assert_frame_equal(df2, exp2)


    def test_read_hmm_cov2(self):
        # Test that the hmm hits are well read, it returns only the hits with coverage >
        # given threshold
        rep_name = "acba.007.p01.13"
        replicon_id = 'ACBA.007.P01_13'

        replicon_path = self.find_data('Replicons', rep_name + '.fst')
        prot_file = self.find_data('Proteins', replicon_id + '.prt')

        self.args.replicon = replicon_path
        cfg = Config(self.args)

        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)

        infile = self.find_data("fictive_results", "{}_intI.res".format(replicon_id))

        df1 = read_hmm(rep_name, prot_db, infile, cfg, coverage=0.7)
        exp1 = pd.DataFrame(data={"Accession_number": [rep_name] * 2,
                                  "query_name": ["intI_Cterm"] * 2,
                                  "ID_query": ["-", "-"],
                                  "ID_prot": ["ACBA.007.P01_13_1", "ACBA.007.P01_13_2"],
                                  "strand": [1, -1],
                                  "pos_beg": [55, 905], "pos_end": [1014, 1609],
                                  "evalue": [1.9e-25, 1e-3]},
                            index=[0, 1])
        exp1 = exp1[["Accession_number", "query_name", "ID_query", "ID_prot",
                     "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df1, exp1)


    def test_read_multi(self):
        # Test reading hmm results when there are multiple hits: 2 hits on the same protein: keep
        # only the one with the best evalue. 2 hits on 2 different proteins: keep the 2 proteins.
        replicon_id = 'ACBA.0917.00019'
        contig_id = 'ACBA.0917.00019.0001'
        result_dir_expected = self.find_data("Results_Integron_Finder_{}.gembase".format(replicon_id))
        replicon_path = self.find_data('Gembase', 'Gembase1', 'Replicons', replicon_id + '.fna')
        prot_file = os.path.join(result_dir_expected, "tmp_{}".format(contig_id), contig_id + '.prt')

        self.args.gembase = True
        self.args.replicon = replicon_path
        cfg = Config(self.args)

        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        with self.catch_log():
            prot_db = GembaseDB(replicon, cfg, prot_file=prot_file)

        infile = self.find_data('fictive_results', "{}_intI_multi.res".format(contig_id))

        df = read_hmm(contig_id, prot_db, infile, cfg)
        exp = pd.DataFrame(data={"Accession_number": [contig_id] * 2,
                                 "query_name": ["Phage_integrase"] * 2,
                                 "ID_query": ["PF00589.16"] * 2,
                                 "ID_prot": ["ACBA.0917.00019.i0001_00298", "ACBA.0917.00019.i0001_00338"],
                                 "strand": [-1, -1],
                                 "pos_beg": [311597, 350328],
                                 "pos_end": [312631, 351248],
                                 "evalue": [5.5e-66, 3.4e-51]},
                           index=[0, 1])
        exp = exp[["Accession_number", "query_name", "ID_query", "ID_prot",
                   "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df, exp)


    def test_read_hmm_hit_no_hsp(self):
        # read hmm but skip hits without any hsp.
        rep_name = "NZ_AP023221"
        replicon_id = 'NZ_AP023221.2'

        replicon_path = self.find_data('Replicons', rep_name + '.fasta')
        prot_file = self.find_data('Proteins', replicon_id + '.prt')

        self.args.replicon = replicon_path
        cfg = Config(self.args)

        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)

        infile = self.find_data("fictive_results",
                                f"{replicon_id}_phage_int.res")

        df1 = read_hmm(replicon_id, prot_db, infile, cfg)

        exp1 = pd.DataFrame(data={"Accession_number": [replicon_id] * 2,
                                  "query_name": ["Phage_integrase"] * 2,
                                  "ID_query": ["PF00589.16", "PF00589.16"],
                                  "ID_prot": ["NZ_AP023221.2_109", "NZ_AP023221.2_106"],
                                  "strand": [-1, -1],
                                  "pos_beg": [93412, 91461],
                                  "pos_end": [94152, 92201],
                                  "evalue": [7.300000e-30, 7.800000e-29]},
                            index=[0, 1])
        exp1 = exp1[["Accession_number", "query_name", "ID_query", "ID_prot",
                     "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df1, exp1)
