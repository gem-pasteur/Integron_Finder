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
import pandas as pd
import pandas.testing as pdt

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
        self.replicon_id = "ACBA.007.P01_13"
        self.replicon_size = 20301
        self.length_cm = 47  # length in 'CLEN' (value for model attc_4.cm)
        self.dtype = {"Accession_number": "str",
                      "cm_attC": "str",
                      "cm_debut": "int",
                      "cm_fin": "int",
                      "pos_beg": "int",
                      "pos_end": "int",
                      "evalue": "float",
           }

    def test_nofile(self):
        """
        Test that the function returns an empty dataframe if the given infernal file does not
        exist.

        """
        filename = "infernal.txt"
        df = infernal.read_infernal(filename,
                                    self.replicon_id, self.replicon_size,
                                    self.length_cm)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        expect = expect.astype(self.dtype)
        pdt.assert_frame_equal(df, expect)

    def test_nohit(self):
        """
        Test that if the infernal file exists but there is no hit
        inside, it returns an empty dataframe.
        """
        filename = self.find_data(os.path.join("fictive_results", "{}_attc_table-empty.res".format(self.replicon_id)))
        df = infernal.read_infernal(filename,
                                    self.replicon_id, self.replicon_size,
                                    self.length_cm)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        expect = expect.astype(self.dtype)
        pdt.assert_frame_equal(df, expect)

    def test_evalue_thres(self):
        """
        Test that if the infernal file exists and there are hits, but
        the given evalue threshold is smaller than the hits thresholds:
        no hit kept, should return an empty dataframe.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_{}".format(self.replicon_name),
                                               "tmp_{}".format(self.replicon_id),
                                               "{}_attc_table.res".format(self.replicon_id)))
        df = infernal.read_infernal(filename,
                                    self.replicon_id, self.replicon_size,
                                    self.length_cm, evalue=1e-10)
        expect = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                       "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        expect = expect.astype(self.dtype)
        pdt.assert_frame_equal(df, expect)


    def test_generate_df(self):
        """
        Test that if the infernal file exists and there are hits, it returns the
        dataframe corresponding to it.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_{}".format(self.replicon_name),
                                               "tmp_{}".format(self.replicon_id),
                                               "{}_attc_table.res".format(self.replicon_id)))
        df = infernal.read_infernal(filename,
                                    self.replicon_id, self.replicon_size,
                                    self.length_cm)
        expect = pd.DataFrame({"Accession_number": [self.replicon_id, self.replicon_id, self.replicon_id,],
                                "cm_attC": ["attc_4", "attc_4", "attc_4",],
                                "cm_debut": [1, 1, 1],
                                "cm_fin": [47, 47, 47],
                                "pos_beg": [17825, 19080, 19618],
                                "pos_end": [17884, 19149, 19726],
                                "sens": ["-", "-", "-"],
                                "evalue": [1e-9, 1e-4, 1.1e-7]})
        expect = expect.astype(self.dtype)
        pdt.assert_frame_equal(df, expect)

    def test_attcsize_minthres(self):
        """
        Test that the filter by a minimum attc size works.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_{}".format(self.replicon_name),
                                               "tmp_{}".format(self.replicon_id),
                                               "{}_attc_table.res".format(self.replicon_id)))
        df = infernal.read_infernal(filename,
                                    self.replicon_id, self.replicon_size,
                                    self.length_cm, size_min_attc=60)
        expect = pd.DataFrame({"Accession_number": [self.replicon_id, self.replicon_id],
                                "cm_attC": ["attc_4", "attc_4"],
                                "cm_debut": [1, 1],
                                "cm_fin": [47, 47],
                                "pos_beg": [19080, 19618],
                                "pos_end": [19149, 19726],
                                "sens": ["-", "-"],
                                "evalue": [1e-4, 1.1e-7]})

        expect = expect.astype(self.dtype)
        pdt.assert_frame_equal(df, expect)

    def test_attcsize_maxthres(self):
        """
        Test that the filter by a maximum attc size works.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_{}".format(self.replicon_name),
                                               "tmp_{}".format(self.replicon_id),
                                               "{}_attc_table.res".format(self.replicon_id)))
        df = infernal.read_infernal(filename,
                                    self.replicon_id, self.replicon_size,
                                    self.length_cm, size_max_attc=100)
        expect = pd.DataFrame({"Accession_number": [self.replicon_id, self.replicon_id],
                               "cm_attC": ["attc_4", "attc_4"],
                               "cm_debut": [1, 1],
                               "cm_fin": [47, 47],
                               "pos_beg": [17825, 19080],
                               "pos_end": [17884, 19149],
                               "sens": ["-", "-"],
                               "evalue": [1e-9, 1e-4]})
        expect = expect.astype(self.dtype)
        pdt.assert_frame_equal(df, expect)

    def test_filter_evalue_thres(self):
        """
        Test that the filter by a maximum attc size works.
        """
        filename = self.find_data(os.path.join("Results_Integron_Finder_{}".format(self.replicon_name),
                                               "tmp_{}".format(self.replicon_id),
                                               "{}_attc_table.res".format(self.replicon_id)))
        df = infernal.read_infernal(filename,
                                    self.replicon_id, self.replicon_size,
                                    self.length_cm, evalue=1e-8)
        expect = pd.DataFrame({"Accession_number": [self.replicon_id],
                               "cm_attC": ["attc_4"],
                               "cm_debut": [1],
                               "cm_fin": [47],
                               "pos_beg": [17825],
                               "pos_end": [17884],
                               "sens": ["-"],
                               "evalue": [1e-9]}
                              )
        expect = expect.astype(self.dtype)
        pdt.assert_frame_equal(df, expect)

    def test_no_total_cm_match_strandp(self):
        """
        Test that when the model did not completely match on the sequence,
        the start and end positions of hit are well recalculated. All hits are on strand +
        """
        filename = self.find_data("fictive_results",
                                  f"{self.replicon_id}_attc_table-partial.res")

        df = infernal.read_infernal(filename,
                                    self.replicon_id, self.replicon_size,
                                    self.length_cm)
        expect = pd.DataFrame({"Accession_number": [self.replicon_id, self.replicon_id, self.replicon_id],
                               "cm_attC": ["attc_4", "attc_4", "attc_4"],
                               "cm_debut": [1, 1, 10],
                               "cm_fin": [40, 47, 47],
                               "pos_beg": [17825, 19080, 19609],
                               "pos_end": [17891, 19149, 19726],
                               "sens": ["+", "+", "+"],
                               "evalue": [1e-9,  1e-4, 1.1e-7]})
        expect = expect.astype(self.dtype)
        pdt.assert_frame_equal(df, expect)

    def test_no_total_cm_match_strandm(self):
        """
        Test that when the model did not completely match on the sequence,
        the start and end positions of hit are well recalculated. All hits are on strand -
        """
        filename = self.find_data(
            os.path.join("fictive_results", "{}_attc_table-partialm.res".format(self.replicon_id)))
        df = infernal.read_infernal(filename,
                                    self.replicon_id, self.replicon_size,
                                    self.length_cm)
        expect = pd.DataFrame({"Accession_number": [self.replicon_id, self.replicon_id, self.replicon_id],
                               "cm_attC": ["attc_4", "attc_4", "attc_4"],
                               "cm_debut": [1, 1, 10],
                               "cm_fin": [40, 47, 47],
                               "pos_beg": [17818, 19080, 19618],
                               "pos_end": [17884, 19149, 19735],
                               "sens": ["-", "-", "-"],
                               "evalue": [1e-9, 1e-4, 1.1e-7]})

        expect = expect.astype(self.dtype)
        pdt.assert_frame_equal(df, expect)

    def test_attc_overflow_pos(self):
        """test when model is truncated and on very first or last replicon pos"""
        filename = self.find_data('fictive_results', '37_0_200_subseq_overflow_attc_table.res')
        replicon_id = '37'
        replicon_size = 3109
        model_len = 47
        df = infernal.read_infernal(filename,
                                    replicon_id, replicon_size,
                                    model_len)

        expect = pd.DataFrame({"Accession_number": [replicon_id, replicon_id],
                               "cm_attC": ["attc_4", "attc_4"],
                               "cm_debut": [4, 1],
                               "cm_fin": [44, 41],
                               "pos_beg": [1, 3065],
                               "pos_end": [126, 3109],
                               "sens": ["-", "+"],
                               "evalue": [0.0024,0.0023]})

        expect = expect.astype(self.dtype)
        pdt.assert_frame_equal(df, expect)