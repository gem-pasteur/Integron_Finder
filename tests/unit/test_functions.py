#!/usr/bin/env python
# coding: utf-8


"""
Unit tests for functions in integron_finder.py
"""

import integron_finder
import pandas as pd
import unittest
import os
import pandas.util.testing as pdt


class TestFunctions(unittest.TestCase):

    def setUp(self):
        """
        Define variables common to all tests
        """
        self.rep_name = "acba.007.p01.13"
        integron_finder.replicon_name = self.rep_name  # name of replicon in given fasta file
        integron_finder.length_cm = 47  # length in 'CLEN' (value for model attc_4.cm)
        integron_finder.DISTANCE_THRESHOLD = 4000  # (4kb at least between 2 different arrays)
        integron_finder.SIZE_REPLICON = 20301  # size of acba.007.p01.13


    def test_search_attc_empty(self):
        """
        Test that when there are no attC sites detected, the attc array is empty.
        """
        attc_file = os.path.join("tests", "data", "Results_Integron_Finder_" + self.rep_name,
                                 "other", self.rep_name + "_attc_table-empty.res")
        # Construct attC dataframe (read from infernal file)
        attc_df = integron_finder.read_infernal(attc_file)
        attc_array = integron_finder.search_attc(attc_df, True)
        self.assertEqual(len(attc_array), 0)
        attc_res = []
        self.assertEqual(attc_array, attc_res)


    def test_search_attc_uniq(self):
        """
        Test that it finds a unique attc array when giving a table with 3 attC sites
        on the same strand and separated by less than 4kb each.
        """
        attc_file = os.path.join("tests", "data", "Results_Integron_Finder_" + self.rep_name,
                                 "other", self.rep_name + "_attc_table.res")
        # Construct attC dataframe (read from infernal file)
        attc_df = integron_finder.read_infernal(attc_file)
        # search attC arrays, keeping palindromes
        # 2 attc sites are in the same array if they are on the same strand, and separated by
        # a distance less than 4kb

        attc_array = integron_finder.search_attc(attc_df, True)
        self.assertEqual(len(attc_array), 1)

        # Construct expected output:
        attc_res = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut", "cm_fin", "pos_beg", "pos_end", "sens", "evalue"], dtype='int')
        attc_res = attc_res.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 17825,
                                    "pos_end": 17884, "sens": "-", "evalue": 1e-9},
                                    ignore_index=True)
        attc_res = attc_res.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 19080,
                                    "pos_end": 19149, "sens": "-", "evalue": 1e-4},
                                    ignore_index=True)
        attc_res = attc_res.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 19618,
                                    "pos_end": 19726, "sens": "-", "evalue": 1.1e-7},
                                    ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        attc_res[intcols] = attc_res[intcols].astype(int)
        pdt.assert_frame_equal(attc_res, attc_array[0])

    def test_search_attc_diff_strand(self):
        """
        Test that it finds a 2 attc arrays when giving a table with 3 attC sites
        on the same strand and 1 on the other strand, all separated by less than 4kb each.
        """
        attc_file = os.path.join("tests", "data", "Results_Integron_Finder_" + self.rep_name,
                                 "other", self.rep_name + "_attc_table.res")
        # Construct attC dataframe (read from infernal file)
        attc_df = integron_finder.read_infernal(attc_file)
        # Add another attC on the opposite strand
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                  "cm_debut": 1, "cm_fin": 47, "pos_beg": 15000,
                                  "pos_end": 16000, "sens": "+", "evalue": 0.19},
                                  ignore_index=True)
        # search attC arrays, keeping palindromes
        # 2 attc sites are in the same array if they are on the same strand, and separated by
        # a distance less than 4kb
        attc_array = integron_finder.search_attc(attc_df, True)
        self.assertEqual(len(attc_array), 2)

        # Construct expected outputs:
        attc_res = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut", "cm_fin",
                                         "pos_beg", "pos_end", "sens", "evalue"])
        attc_res = attc_res.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 17825,
                                    "pos_end": 17884, "sens": "-", "evalue": 1e-9},
                                    ignore_index=True)
        attc_res = attc_res.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 19080,
                                    "pos_end": 19149, "sens": "-", "evalue": 1e-4},
                                    ignore_index=True)
        attc_res = attc_res.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 19618,
                                    "pos_end": 19726, "sens": "-", "evalue": 1.1e-7},
                                    ignore_index=True)
        attc_res2 = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut", "cm_fin",
                                          "pos_beg", "pos_end", "sens", "evalue"])
        attc_res2 = attc_res2.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                      "cm_debut": 1, "cm_fin": 47, "pos_beg": 15000,
                                      "pos_end": 16000, "sens": "+", "evalue": 0.19},
                                     ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        attc_res[intcols] = attc_res[intcols].astype(int)
        attc_res2[intcols] = attc_res2[intcols].astype(int)
        pdt.assert_frame_equal(attc_res, attc_array[1])
        pdt.assert_frame_equal(attc_res2, attc_array[0])

