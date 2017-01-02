#!/usr/bin/env python
# coding: utf-8

"""
Unit tests read_hmm function of integron_finder
"""

import integron_finder
import pandas as pd
import unittest
import os
import pandas.util.testing as pdt
import argparse

class TestFunctions(unittest.TestCase):

    def setUp(self):
        """
        Define variables common to all tests
        """
        self.rep_name = "acba.007.p01.13"
        # Simulate argparse to get argument
        parser = argparse.ArgumentParser(description='Process some integers.')
        parser.add_argument("--gembase", help="gembase format", action="store_true")
        args = parser.parse_args([])
        integron_finder.args = args

    def test_read_empty(self):
        """
        Test that when there are no hits in the hmm result file, it returns an empty
        dataframe, without error.
        """
        infile = os.path.join("tests", "data", "fictive_results",
                              self.rep_name + "_intI-empty.res")
        df = integron_finder.read_hmm(self.rep_name, infile)
        exp = pd.DataFrame(columns=["Accession_number", "query_name", "ID_query", "ID_prot",
                                    "strand", "pos_beg", "pos_end", "evalue"])

        intcols = ["pos_beg", "pos_end", "strand"]
        floatcol = ["evalue"]
        exp[intcols] = exp[intcols].astype(int)
        exp[floatcol] = exp[floatcol].astype(float)
        pdt.assert_frame_equal(df, exp)

    def test_read_hmm(self):
        """
        Test that the hmm hits are well read
        """
        infile = os.path.join("tests", "data", "Results_Integron_Finder_" + self.rep_name,
                                 "other", self.rep_name + "_intI.res")
        df = integron_finder.read_hmm(self.rep_name, infile)
        exp = pd.DataFrame(data={"Accession_number": self.rep_name, "query_name": "intI_Cterm",
                                 "ID_query": "-", "ID_prot": "ACBA.007.P01_13_1", "strand": 1,
                                 "pos_beg": 55, "pos_end": 1014, "evalue": 1.9e-25},
                           index=[0])
        exp = exp[["Accession_number", "query_name", "ID_query", "ID_prot",
                   "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df, exp)

    def test_read_hmm_gembase(self):
        """
        Test that the hmm hits are well read, when the gembase format is used (.prt file is
        provided, prodigal is not used to find the proteins).
        """
        parser = argparse.ArgumentParser(description='Process some integers.')
        parser.add_argument("--gembase", help="gembase format", action="store_true")
        args = parser.parse_args(["--gembase"])
        integron_finder.args = args
        infile = os.path.join("tests", "data", "fictive_results", self.rep_name +
                                 "_intI-gembase.res")
        df = integron_finder.read_hmm(self.rep_name, infile)
        exp = pd.DataFrame(data={"Accession_number": self.rep_name, "query_name": "intI_Cterm",
                                 "ID_query": "-", "ID_prot": "ACBA007p01a_000009", "strand": 1,
                                 "pos_beg": 55, "pos_end": 1014, "evalue": 1.9e-25},
                           index=[0])
        exp = exp[["Accession_number", "query_name", "ID_query", "ID_prot",
                   "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df, exp)

    def test_read_hmm_evalue(self):
        """
        Test that the hmm hits are well read, and returned only if evalue is < to the
        given threshold.
        """
        infile = os.path.join("tests", "data", "Results_Integron_Finder_" + self.rep_name,
                                 "other", self.rep_name + "_intI.res")
        df1 = integron_finder.read_hmm(self.rep_name, infile, evalue=1.95e-25)
        exp1 = pd.DataFrame(data={"Accession_number": self.rep_name, "query_name": "intI_Cterm",
                                  "ID_query": "-", "ID_prot": "ACBA.007.P01_13_1", "strand": 1,
                                  "pos_beg": 55, "pos_end": 1014, "evalue": 1.9e-25},
                            index=[0])
        exp1 = exp1[["Accession_number", "query_name", "ID_query", "ID_prot",
                     "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df1, exp1)
        df2 = integron_finder.read_hmm(self.rep_name, infile, evalue=1.9e-25)
        exp2 = pd.DataFrame(columns=["Accession_number", "query_name", "ID_query", "ID_prot",
                                     "strand", "pos_beg", "pos_end", "evalue"])

        intcols = ["pos_beg", "pos_end", "strand"]
        floatcol = ["evalue"]
        exp2[intcols] = exp2[intcols].astype(int)
        exp2[floatcol] = exp2[floatcol].astype(float)
        pdt.assert_frame_equal(df2, exp2)

    def test_read_hmm_evalue2(self):
        """
        Test that the hmm hits are well read, it returns only the hits with evalue < given
        threshold
        """
        infile = os.path.join("tests", "data", "fictive_results",
                              self.rep_name + "_intI.res")
        df1 = integron_finder.read_hmm(self.rep_name, infile, evalue=1e-3)
        exp1 = pd.DataFrame(data={"Accession_number": [self.rep_name] * 2,
                                  "query_name": ["intI_Cterm"] * 2,
                                  "ID_query": ["-", "-"],
                                  "ID_prot": ["ACBA.007.P01_13_1", "ACBA.007.P01_13_3"],
                                  "strand": [1, 1],
                                  "pos_beg": [55, 3000], "pos_end": [1014, 3500],
                                  "evalue": [1.9e-25, 2e-25]},
                            index=[0, 1])
        exp1 = exp1[["Accession_number", "query_name", "ID_query", "ID_prot",
                     "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df1, exp1)

    def test_read_hmm_cov(self):
        """
        Test that the hmm hits are well read, and returned only if coverage is > to the
        given threshold.
        """
        infile = os.path.join("tests", "data", "Results_Integron_Finder_" + self.rep_name,
                                 "other", self.rep_name + "_intI.res")
        df1 = integron_finder.read_hmm(self.rep_name, infile, coverage=0.945)
        exp1 = pd.DataFrame(data={"Accession_number": self.rep_name, "query_name": "intI_Cterm",
                                  "ID_query": "-", "ID_prot": "ACBA.007.P01_13_1", "strand": 1,
                                  "pos_beg": 55, "pos_end": 1014, "evalue": 1.9e-25},
                            index=[0])
        exp1 = exp1[["Accession_number", "query_name", "ID_query", "ID_prot",
                     "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df1, exp1)
        df2 = integron_finder.read_hmm(self.rep_name, infile, coverage=0.95)
        exp2 = pd.DataFrame(columns=["Accession_number", "query_name", "ID_query", "ID_prot",
                                     "strand", "pos_beg", "pos_end", "evalue"])
        intcols = ["pos_beg", "pos_end", "strand"]
        floatcol = ["evalue"]
        exp2[intcols] = exp2[intcols].astype(int)
        exp2[floatcol] = exp2[floatcol].astype(float)
        pdt.assert_frame_equal(df2, exp2)

    def test_read_hmm_cov2(self):
        """
        Test that the hmm hits are well read, it returns only the hits with coverage >
        given threshold
        """
        infile = os.path.join("tests", "data", "fictive_results",
                              self.rep_name + "_intI.res")
        df1 = integron_finder.read_hmm(self.rep_name, infile, coverage=0.7)
        exp1 = pd.DataFrame(data={"Accession_number": [self.rep_name] * 2,
                                  "query_name": ["intI_Cterm"] * 2,
                                  "ID_query": ["-", "-"],
                                  "ID_prot": ["ACBA.007.P01_13_1", "ACBA.007.P01_13_2"],
                                  "strand": [1, 1],
                                  "pos_beg": [55, 2000], "pos_end": [1014, 2500],
                                  "evalue": [1.9e-25, 1e-3]},
                            index=[0, 1])
        exp1 = exp1[["Accession_number", "query_name", "ID_query", "ID_prot",
                     "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df1, exp1)

    def test_read_multi(self):
        """
        Test reading hmm results when there are multiple hits: 2 hits on the same protein: keep
        only the one with the best evalue. 2 hits on 2 different proteins: keep the 2 proteins.
        """
        parser = argparse.ArgumentParser(description='Process some integers.')
        parser.add_argument("--gembase", help="gembase format", action="store_true")
        args = parser.parse_args(["--gembase"])
        integron_finder.args = args
        infile = os.path.join("tests", "data", "fictive_results", self.rep_name +
                              "_intI-multi.res")
        df = integron_finder.read_hmm(self.rep_name, infile)
        exp = pd.DataFrame(data={"Accession_number": [self.rep_name] * 2,
                                 "query_name": ["intI_Cterm"] * 2, "ID_query": ["-"] * 2,
                                 "ID_prot": ["ACBA007p01a_000009", "ACBA007p01a_000008"],
                                 "strand": [1, -1],
                                 "pos_beg": [55, 1], "pos_end": [1014, 50],
                                 "evalue": [4.5e-25, 2.3e-25]},
                           index=[0, 1])
        exp = exp[["Accession_number", "query_name", "ID_query", "ID_prot",
                   "strand", "pos_beg", "pos_end", "evalue"]]
        pdt.assert_frame_equal(df, exp)
