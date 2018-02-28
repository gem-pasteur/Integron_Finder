#!/usr/bin/env python
# coding: utf-8


"""
Unit tests read_infernal function of integron_finder
"""

import os
import pandas as pd
import pandas.util.testing as pdt

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import infernal


class TestReadInfernal(IntegronTest):

    def setUp(self):
        """
        Define variables common to all tests
        """
        self.replicon_name = "acba.007.p01.13"
        self.length_cm = 47  # length in 'CLEN' (value for model attc_4.cm)


    def test_nofile(self):
        """
        Test that the function returns an empty dataframe if the given infernal file does not
        exist.

        """
        filename = "infernal.txt"
        df = infernal.read_infernal(filename, self.replicon_name, self.length_cm)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        pdt.assert_frame_equal(df, expect)

    def test_nohit(self):
        """
        Test that if the infernal file exists but there is no hit
        inside, it returns an empty dataframe.
        """
        filename = self.find_data(os.path.join("fictive_results", self.replicon_name + "_attc_table-empty.res"))
        df = infernal.read_infernal(filename, self.replicon_name, self.length_cm)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        pdt.assert_frame_equal(df, expect)

    def test_evalue_thres(self):
        """
        Test that if the infernal file exists and there are hits, but
        the given evalue threshold is smaller than the hits thresholds:
        no hit kept, should return an empty dataframe.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_" + self.replicon_name,
                                  "other", self.replicon_name + "_attc_table.res"))
        df = infernal.read_infernal(filename, self.replicon_name, self.length_cm, evalue=1e-10)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        pdt.assert_frame_equal(df, expect)

    def test_generate_df(self):
        """
        Test that if the infernal file exists and there are hits, it returns the
        dataframe corresponding to it.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_" + self.replicon_name,
                                  "other", self.replicon_name + "_attc_table.res"))
        df = infernal.read_infernal(filename, self.replicon_name, self.length_cm)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 17825,
                                "pos_end": 17884, "sens": "-", "evalue": 1e-9},
                               ignore_index=True)
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 19080,
                                "pos_end": 19149, "sens": "-", "evalue": 1e-4},
                               ignore_index=True)
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 19618,
                                "pos_end": 19726, "sens": "-", "evalue": 1.1e-7},
                               ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        expect[intcols] = expect[intcols].astype(int)
        pdt.assert_frame_equal(df, expect)

    def test_attcsize_minthres(self):
        """
        Test that the filter by a minimum attc size works.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_" + self.replicon_name,
                                  "other", self.replicon_name + "_attc_table.res"))
        df = infernal.read_infernal(filename, self.replicon_name, self.length_cm, size_min_attc=60)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 19080,
                                "pos_end": 19149, "sens": "-", "evalue": 1e-4},
                               ignore_index=True)
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 19618,
                                "pos_end": 19726, "sens": "-", "evalue": 1.1e-7},
                               ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        expect[intcols] = expect[intcols].astype(int)
        pdt.assert_frame_equal(df, expect)

    def test_attcsize_maxthres(self):
        """
        Test that the filter by a maximum attc size works.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_" + self.replicon_name,
                                  "other", self.replicon_name + "_attc_table.res"))
        df = infernal.read_infernal(filename, self.replicon_name, self.length_cm, size_max_attc=100)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 17825,
                                "pos_end": 17884, "sens": "-", "evalue": 1e-9},
                               ignore_index=True)
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 19080,
                                "pos_end": 19149, "sens": "-", "evalue": 1e-4},
                               ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        expect[intcols] = expect[intcols].astype(int)
        pdt.assert_frame_equal(df, expect)

    def test_filter_evalue_thres(self):
        """
        Test that the filter by a maximum attc size works.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_" + self.replicon_name,
                                  "other", self.replicon_name + "_attc_table.res"))
        df = infernal.read_infernal(filename, self.replicon_name, self.length_cm, evalue=1e-8)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 17825,
                                "pos_end": 17884, "sens": "-", "evalue": 1e-9},
                               ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        expect[intcols] = expect[intcols].astype(int)
        pdt.assert_frame_equal(df, expect)

    def test_no_total_cm_match_strandp(self):
        """
        Test that when the model did not completely match on the sequence,
        the start and end positions of hit are well recalculated. All hits are on strand +
        """
        filename = self.find_data(os.path.join("fictive_results", self.replicon_name + "_attc_table-partial.res"))
        df = infernal.read_infernal(filename, self.replicon_name, self.length_cm)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 40, "pos_beg": 17825,
                                "pos_end": 17891, "sens": "+", "evalue": 1e-9},
                               ignore_index=True)
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 19080,
                                "pos_end": 19149, "sens": "+", "evalue": 1e-4},
                               ignore_index=True)
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 10, "cm_fin": 47, "pos_beg": 19609,
                                "pos_end": 19726, "sens": "+", "evalue": 1.1e-7},
                               ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        expect[intcols] = expect[intcols].astype(int)
        pdt.assert_frame_equal(df, expect)

    def test_no_total_cm_match_strandm(self):
        """
        Test that when the model did not completely match on the sequence,
        the start and end positions of hit are well recalculated. All hits are on strand -
        """
        filename = self.find_data(os.path.join("fictive_results", self.replicon_name + "_attc_table-partialm.res"))
        df = infernal.read_infernal(filename, self.replicon_name, self.length_cm)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 40, "pos_beg": 17818,
                                "pos_end": 17884, "sens": "-", "evalue": 1e-9},
                               ignore_index=True)
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 19080,
                                "pos_end": 19149, "sens": "-", "evalue": 1e-4},
                               ignore_index=True)
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 10, "cm_fin": 47, "pos_beg": 19618,
                                "pos_end": 19735, "sens": "-", "evalue": 1.1e-7},
                               ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        expect[intcols] = expect[intcols].astype(int)
        pdt.assert_frame_equal(df, expect)

    def test_attcsize_minthres(self):
        """
        Test that the filter by a minimum attc size works.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_" + self.replicon_name,
                                  "other", self.replicon_name + "_attc_table.res"))
        df = infernal.read_infernal(filename, self.replicon_name, self.length_cm, size_min_attc=60)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 19080,
                                "pos_end": 19149, "sens": "-", "evalue": 1e-4},
                               ignore_index=True)
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 19618,
                                "pos_end": 19726, "sens": "-", "evalue": 1.1e-7},
                               ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        expect[intcols] = expect[intcols].astype(int)
        pdt.assert_frame_equal(df, expect)

    def test_attcsize_maxthres(self):
        """
        Test that the filter by a maximum attc size works.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_" + self.replicon_name,
                                  "other", self.replicon_name + "_attc_table.res"))
        df = infernal.read_infernal(filename, self.replicon_name, self.length_cm, size_max_attc=100)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 17825,
                                "pos_end": 17884, "sens": "-", "evalue": 1e-9},
                               ignore_index=True)
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 19080,
                                "pos_end": 19149, "sens": "-", "evalue": 1e-4},
                               ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        expect[intcols] = expect[intcols].astype(int)
        pdt.assert_frame_equal(df, expect)

    def test_filter_evalue_thres(self):
        """
        Test that the filter by a maximum attc size works.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_" + self.replicon_name,
                                  "other", self.replicon_name + "_attc_table.res"))
        df = infernal.read_infernal(filename, self.replicon_name, self.length_cm, evalue=1e-8)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        expect = expect.append({"Accession_number": self.replicon_name, "cm_attC": "attC_4",
                                "cm_debut": 1, "cm_fin": 47, "pos_beg": 17825,
                                "pos_end": 17884, "sens": "-", "evalue": 1e-9},
                               ignore_index=True)
        # convert positions to int
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        expect[intcols] = expect[intcols].astype(int)
        pdt.assert_frame_equal(df, expect)

