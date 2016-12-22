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
        Test that it finds a size 2 attc array when giving a table with 3 attC sites
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


    def test_search_attc_dist_same_strand(self):
        """
        Test that it finds a size 2 attc arrays when giving a table with 3 attC sites
        on the same strand and separated by less than 4 kb, and another 1 separated from the
        others by more than 4kb.
        """
        attc_file = os.path.join("tests", "data", "Results_Integron_Finder_" + self.rep_name,
                                 "other", self.rep_name + "_attc_table.res")
        # Construct attC dataframe (read from infernal file)
        attc_df = integron_finder.read_infernal(attc_file)
        # Add another attC at more than 4kb, same strand
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                  "cm_debut": 1, "cm_fin": 47, "pos_beg": 12900,
                                  "pos_end": 13800, "sens": "-", "evalue": 1e-3},
                                  ignore_index=True)
        attc_df.sort_values(["Accession_number", "pos_beg", "evalue"], inplace=True)
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
                                      "cm_debut": 1, "cm_fin": 47, "pos_beg": 12900,
                                      "pos_end": 13800, "sens": "-", "evalue": 1e-03},
                                     ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        attc_res[intcols] = attc_res[intcols].astype(int)
        attc_res2[intcols] = attc_res2[intcols].astype(int)
        pdt.assert_frame_equal(attc_res, attc_array[1])
        pdt.assert_frame_equal(attc_res2, attc_array[0])


    def test_search_attc_dist_diff_strand(self):
        """
        Test that it finds a size 3 attc array when giving a table with:
        - 3 attC sites on the same strand (-) and separated by less than 4 kb
        - 2 other attC sites separated by less than 4kb but on the other strand (+)
        - 1 other attC site , also on strand +, but separated by more than 4kb.
        """
        attc_file = os.path.join("tests", "data", "Results_Integron_Finder_" + self.rep_name,
                                 "other", self.rep_name + "_attc_table.res")
        # Construct attC dataframe (read from infernal file)
        attc_df = integron_finder.read_infernal(attc_file)
        # Add another attC at more than 4kb, same strand
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                  "cm_debut": 1, "cm_fin": 47, "pos_beg": 15800,
                                  "pos_end": 16000, "sens": "+", "evalue": 1e-3},
                                  ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                  "cm_debut": 1, "cm_fin": 47, "pos_beg": 12000,
                                  "pos_end": 12500, "sens": "+", "evalue": 1e-3},
                                  ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                  "cm_debut": 1, "cm_fin": 47, "pos_beg": 7100,
                                  "pos_end": 8200, "sens": "+", "evalue": 1e-3},
                                  ignore_index=True)
        attc_df.sort_values(["Accession_number", "pos_beg", "evalue"], inplace=True)
        # search attC arrays, keeping palindromes
        # 2 attc sites are in the same array if they are on the same strand, and separated by
        # a distance less than 4kb
        attc_array = integron_finder.search_attc(attc_df, True)
        self.assertEqual(len(attc_array), 3)

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
                                      "cm_debut": 1, "cm_fin": 47, "pos_beg": 12000,
                                      "pos_end": 12500, "sens": "+", "evalue": 1e-03},
                                     ignore_index=True)
        attc_res2 = attc_res2.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                      "cm_debut": 1, "cm_fin": 47, "pos_beg": 15800,
                                      "pos_end": 16000, "sens": "+", "evalue": 1e-03},
                                     ignore_index=True)
        attc_res3 = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut", "cm_fin",
                                          "pos_beg", "pos_end", "sens", "evalue"])
        attc_res3 = attc_res3.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                      "cm_debut": 1, "cm_fin": 47, "pos_beg": 7100,
                                      "pos_end": 8200, "sens": "+", "evalue": 1e-03},
                                     ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        attc_res[intcols] = attc_res[intcols].astype(int)
        attc_res2[intcols] = attc_res2[intcols].astype(int)
        attc_res3[intcols] = attc_res3[intcols].astype(int)
        pdt.assert_frame_equal(attc_res, attc_array[2])
        pdt.assert_frame_equal(attc_res2, attc_array[1])
        pdt.assert_frame_equal(attc_res3, attc_array[0])

    def test_search_attc_uniq_circ(self):
        """
        Test that it finds a unique attc array when giving a table with:
        - 2 attC sites at the begining of the genome, separated by less than 4kb
        - 2 attC sites at the end of the genome, separated by less than 4kb from the
        1st attC site of the genome if we take into account its circularity.
        All 3 attC sites are on the strand -
        """
        attc_df = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut", "cm_fin", "pos_beg", "pos_end", "sens", "evalue"], dtype='int')
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 1000,
                                    "pos_end": 2000, "sens": "-", "evalue": 1e-9},
                                    ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 3000,
                                    "pos_end": 4000, "sens": "-", "evalue": 1e-4},
                                    ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 16000,
                                    "pos_end": 17000, "sens": "-", "evalue": 1.1e-7},
                                    ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 19815,
                                    "pos_end": 20000, "sens": "-", "evalue": 1.1e-7},
                                    ignore_index=True)
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        attc_df[intcols] = attc_df[intcols].astype(int)

        attc_array = integron_finder.search_attc(attc_df, True)
        self.assertEqual(len(attc_array), 1)

        # Output of search_attc is ordered as the cluster is:
        # - 2 last attC of genome
        # - 2 first attC of genomes
        # Whereas input was ordered by begin position. Reorder attc_df to compare with output.
        attc_df = attc_df.reindex([2, 3, 0, 1])
        attc_df.reset_index(inplace=True, drop=True)
        pdt.assert_frame_equal(attc_df, attc_array[0])

    def test_search_attc_uniq_circ_plus(self):
        """
        Test that it finds a unique attc array when giving a table with:
        - 2 attC sites at the begining of the genome, separated by less than 4kb
        - 2 attC sites at the end of the genome, separated by less than 4kb from the
        1st attC site of the genome if we take into account its circularity.
        All 3 attC sites are on the same strand +
        """
        attc_df = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut", "cm_fin", "pos_beg", "pos_end", "sens", "evalue"], dtype='int')
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 1000,
                                    "pos_end": 2000, "sens": "+", "evalue": 1e-9},
                                    ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 3000,
                                    "pos_end": 4000, "sens": "+", "evalue": 1e-4},
                                    ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 16000,
                                    "pos_end": 17000, "sens": "+", "evalue": 1.1e-7},
                                    ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 19815,
                                    "pos_end": 20000, "sens": "+", "evalue": 1.1e-7},
                                    ignore_index=True)
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        attc_df[intcols] = attc_df[intcols].astype(int)

        attc_array = integron_finder.search_attc(attc_df, True)
        self.assertEqual(len(attc_array), 1)

        # Output of search_attc is ordered as the cluster is:
        # - 2 last attC of genome
        # - 2 first attC of genomes
        # Whereas input was ordered by begin position. Reorder attc_df to compare with output.
        attc_df = attc_df.reindex([2, 3, 0, 1])
        attc_df.reset_index(inplace=True, drop=True)
        pdt.assert_frame_equal(attc_df, attc_array[0])


    def test_search_attc_drop_pal(self):
        """
        If there is 1 palindrome attC, check that it keeps the one with the highest evalue,
        and that clusters are then found according to it.
        """
        attc_df = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut", "cm_fin", "pos_beg", "pos_end", "sens", "evalue"], dtype='int')
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 1000,
                                    "pos_end": 2000, "sens": "+", "evalue": 1e-9},
                                    ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 3000,
                                    "pos_end": 4000, "sens": "-", "evalue": 1e-4},
                                    ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 3000,
                                    "pos_end": 4000, "sens": "+", "evalue": 1e-9},
                                    ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 5500,
                                    "pos_end": 7000, "sens": "+", "evalue": 1.1e-7},
                                    ignore_index=True)
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        attc_df[intcols] = attc_df[intcols].astype(int)

        attc_array = integron_finder.search_attc(attc_df, False)
        self.assertEqual(len(attc_array), 1)

        # Construct expected outputs:
        attc_res = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut", "cm_fin",
                                         "pos_beg", "pos_end", "sens", "evalue"])
        attc_res = attc_res.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 1000,
                                    "pos_end": 2000, "sens": "+", "evalue": 1e-9},
                                    ignore_index=True)
        attc_res = attc_res.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 3000,
                                    "pos_end": 4000, "sens": "+", "evalue": 1e-9},
                                    ignore_index=True)
        attc_res = attc_res.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 5500,
                                    "pos_end": 7000, "sens": "+", "evalue": 1.1e-7},
                                    ignore_index=True)

        attc_res[intcols] = attc_res[intcols].astype(int)
        attc_array[0].reset_index(inplace=True, drop=True)
        pdt.assert_frame_equal(attc_res, attc_array[0])


    def test_search_attc_drop_pal_break(self):
        """
        If there is 1 palindrome attC, check that it keeps the one with the highest evalue,
        and that clusters are then found according to it.
        """
        attc_df = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut", "cm_fin", "pos_beg", "pos_end", "sens", "evalue"], dtype='int')
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 1000,
                                    "pos_end": 2000, "sens": "-", "evalue": 1e-9},
                                    ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 3000,
                                    "pos_end": 4000, "sens": "-", "evalue": 1e-4},
                                    ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 3000,
                                    "pos_end": 4000, "sens": "+", "evalue": 1e-9},
                                    ignore_index=True)
        attc_df = attc_df.append({"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 5500,
                                    "pos_end": 7000, "sens": "-", "evalue": 1.1e-7},
                                    ignore_index=True)
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        attc_df[intcols] = attc_df[intcols].astype(int)

        attc_array = integron_finder.search_attc(attc_df, False)
        self.assertEqual(len(attc_array), 3)

        # Construct expected outputs:
        columns = ["Accession_number", "cm_attC", "cm_debut", "cm_fin", "pos_beg",
                   "pos_end", "sens", "evalue"]
        attc_res = pd.DataFrame(data={"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 1000,
                                    "pos_end": 2000, "sens": "-", "evalue": 1e-9},
                                index=[0])
        attc_res = attc_res[columns]

        attc_res2 = pd.DataFrame(data={"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 3000,
                                    "pos_end": 4000, "sens": "+", "evalue": 1e-9},
                                index=[0])
        attc_res2 = attc_res2[columns]

        attc_res3 = pd.DataFrame(data={"Accession_number": self.rep_name, "cm_attC": "attC_4",
                                    "cm_debut": 1, "cm_fin": 47, "pos_beg": 5500,
                                    "pos_end": 7000, "sens": "-", "evalue": 1.1e-7},
                                index=[0])
        attc_res3 = attc_res3[columns]

        attc_res[intcols] = attc_res[intcols].astype(int)
        attc_res2[intcols] = attc_res2[intcols].astype(int)
        attc_res3[intcols] = attc_res3[intcols].astype(int)
        pdt.assert_frame_equal(attc_res2, attc_array[0])
        pdt.assert_frame_equal(attc_res, attc_array[1])
        pdt.assert_frame_equal(attc_res3, attc_array[2])

