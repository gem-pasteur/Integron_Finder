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
        Test that when there are no attC sites detected, the attc array is empty.
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

    # test with gembase format of .prt
    # test filtering evalue
    # test filtering coverage
    # test various hits for same protein: keep best evalue


