#!/usr/bin/env python
# coding: utf-8

"""
Unit tests func_annot function of integron_finder
"""

import integron_finder
import os
import unittest
import distutils.spawn
import shutil
import argparse
import glob
import pandas as pd
import pandas.util.testing as pdt
import numpy as np


class TestFunctions(unittest.TestCase):

    def setUp(self):
        """
        Define variables common to all tests
        """
        self.replicon_path = os.path.join("tests", "data", "acba.007.p01.13.fst")
        self.replicon_name = "acba.007.p01.13"
        self.out_dir = "tmpdir"
        self.hmm_files = [os.path.join("data", "Functional_annotation", "Resfams.hmm")]
        if os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
        os.mkdir(self.out_dir)
        # Define integron_finder variables
        parser = argparse.ArgumentParser(description='Process some integers.')
        parser.add_argument("--gembase", action="store_true")
        args = parser.parse_args([])
        integron_finder.args = args
        integron_finder.SIZE_REPLICON = 20301  # size of acba.007.p01.13
        integron_finder.N_CPU = "1"
        integron_finder.PRODIGAL = distutils.spawn.find_executable("prodigal")
        integron_finder.HMMSEARCH = distutils.spawn.find_executable("hmmsearch")
        integron_finder.MODEL_integrase = os.path.join("data", "Models", "integron_integrase.hmm")
        integron_finder.MODEL_phage_int = os.path.join("data", "Models", "phage-int.hmm")
        integron_finder.PROT_file = os.path.join(self.out_dir, self.replicon_name + ".prt")
        self.exp_files = ["acba.007.p01.13.prt", "acba.007.p01.13_Resfams_fa_table.res",
                          "acba.007.p01.13_intI_table.res", "acba.007.p01.13_phage_int_table.res",
                          "acba.007.p01.13_Resfams_fa.res", "acba.007.p01.13_intI.res",
                          "acba.007.p01.13_phage_int.res", "acba.007.p01.13_subseqprot.tmp"]
        self.exp_files = [os.path.join(self.out_dir, file) for file in self.exp_files]

    def tearDown(self):
        """
        To do after each test. remove output directory if it was generated
        """
        if os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)

    def test_annot_calin(self):
        """
        Test func_annot when the integron is a CALIN (attC but no integrase), with 4 proteins:
        for 3 of them resfam annotations are found, and not for the last 1.
        """
        # Run prodigal to find CDS on replicon (and run hmmsearch on integrase (2 profiles))
        integron_finder.find_integrase(self.replicon_path, self.replicon_name, self.out_dir)
        # Create integron
        integron1 = integron_finder.Integron(self.replicon_name)
        integron_finder.integrons = [integron1]
        # Add only attc sites (no integrase)
        integron1.add_attC(17825, 17884, -1, 7e-9, "attc_4")
        integron1.add_attC(19080, 19149, -1, 7e-4, "attc_4")
        integron1.add_attC(19618, 19726, -1, 7e-7, "attc_4")
        # Add proteins between attC sites
        integron1.add_proteins()
        # Check that proteins dataframe is as expected before annotation
        proteins = pd.DataFrame({"pos_beg": [17375, 17886, 19090, 19721],
                                 "pos_end": [17722, 18665, 19749, 20254],
                                 "strand": [-1] * 4,
                                 "evalue": [np.nan] * 4,
                                 "type_elt": ["protein"] * 4,
                                 "model": ["NA"] * 4,
                                 "distance_2attC": [np.nan] * 4,
                                 "annotation": ["protein"] * 4},
                                index=["ACBA.007.P01_13_20", "ACBA.007.P01_13_21",
                                       "ACBA.007.P01_13_22", "ACBA.007.P01_13_23"])
        proteins = proteins[["pos_beg", "pos_end", "strand", "evalue", "type_elt",
                             "model", "distance_2attC", "annotation"]]
        pdt.assert_frame_equal(proteins, integron1.proteins)

        # Annotate proteins
        integron_finder.func_annot(self.replicon_name, self.out_dir, self.hmm_files)

        # Check that all files generated are as expected
        files_created = glob.glob(os.path.join(self.out_dir, "*"))
        self.assertEqual(set(self.exp_files), set(files_created))

        # Check that annotated proteins are as expected
        proteins.loc["ACBA.007.P01_13_20"] = [17375, 17722, -1, 4.5e-31, "protein",
                                              "RF0066", np.nan, "emrE"]
        proteins.loc["ACBA.007.P01_13_21"] = [17886, 18665, -1, 7.4e-168, "protein",
                                              "RF0027", np.nan, "ANT3"]
        proteins.loc["ACBA.007.P01_13_23"] = [19721, 20254, -1, 6.2e-110, "protein",
                                              "RF0003", np.nan, "AAC3-I"]
        pdt.assert_frame_equal(proteins, integron1.proteins)

    # TODO:
    # annotate in0
    # annotate complete
    # create a _subseqprot.tmp file before running the function
    # several integrons
    # error with hmmsearch
    # evalue threshold (should not work...)
