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
        # Run prodigal to find CDS on replicon (and run hmmsearch on integrase (2 profiles))
        integron_finder.find_integrase(self.replicon_path, self.replicon_name, self.out_dir)

    def tearDown(self):
        """
        To do after each test. remove output directory if it was generated
        """
        if os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
        integron_finder.integrons = []

    def test_annot_calin(self):
        """
        Test func_annot when the integron is a CALIN (attC but no integrase), with 4 proteins:
        for 3 of them resfam annotations are found, and not for the last 1.
        """
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

    def test_annot_calin_empty(self):
        """
        Test func_annot when the integron is a CALIN (attC but no integrase), without any protein:
        nothing to annotate
        """
        # Create integron
        integron1 = integron_finder.Integron(self.replicon_name)
        integron_finder.integrons = [integron1]
        # Add only attc sites (no integrase)
        integron1.add_attC(17825, 17884, -1, 7e-9, "attc_4")
        integron1.add_attC(19080, 19149, -1, 7e-4, "attc_4")
        integron1.add_attC(19618, 19726, -1, 7e-7, "attc_4")

        # check proteins before annotation
        proteins = pd.DataFrame(columns=["pos_beg", "pos_end", "strand",
                                         "evalue", "type_elt", "model",
                                         "distance_2attC", "annotation"])
        proteins = proteins.astype(dtype={"pos_beg": "int", "pos_end": "int", "strand": "int",
                                          "evalue": "float", "type_elt": "str", "model": "str",
                                          "distance_2attC": "float", "annotation": "str"})
        pdt.assert_frame_equal(proteins, integron1.proteins)

        # Annotate proteins
        integron_finder.func_annot(self.replicon_name, self.out_dir, self.hmm_files)

        # Check that all files generated are as expected
        files_created = glob.glob(os.path.join(self.out_dir, "*"))
        exp_files = ["acba.007.p01.13.prt", "acba.007.p01.13_intI_table.res",
                     "acba.007.p01.13_phage_int_table.res", "acba.007.p01.13_intI.res",
                     "acba.007.p01.13_phage_int.res"]
        exp_files = [os.path.join(self.out_dir, file) for file in exp_files]
        self.assertEqual(set(exp_files), set(files_created))

        # Check proteins after annotation
        pdt.assert_frame_equal(proteins, integron1.proteins)

    def test_annot_in0(self):
        """
        Test func_annot when the integron is a in0: only an integrase. There are no proteins
        to annotate, but the _subseqprot.tmp file already exists (not deleted in a last run
        for example...)
        """
        # create empty _subseqprot.tmp file (must be deleted by func_annot)
        open(os.path.join(self.out_dir, "acba.007.p01.13_subseqprot.tmp"), "w").close()
        # Create integron
        integron1 = integron_finder.Integron(self.replicon_name)
        integron_finder.integrons = [integron1]
        # Add integrase
        integron1.add_integrase(55, 1014, "ACBA.007.P01_13_1", 1, 1.9e-25,
                                      "intersection_tyr_intI")
        # check proteins before annotation
        proteins = pd.DataFrame(columns=["pos_beg", "pos_end", "strand",
                                         "evalue", "type_elt", "model",
                                         "distance_2attC", "annotation"])
        proteins = proteins.astype(dtype={"pos_beg": "int", "pos_end": "int", "strand": "int",
                                          "evalue": "float", "type_elt": "str", "model": "str",
                                          "distance_2attC": "float", "annotation": "str"})
        pdt.assert_frame_equal(proteins, integron1.proteins)

        # Annotate proteins
        integron_finder.func_annot(self.replicon_name, self.out_dir, self.hmm_files)
        # Check that all files generated are as expected
        files_created = glob.glob(os.path.join(self.out_dir, "*"))
        exp_files = ["acba.007.p01.13.prt", "acba.007.p01.13_intI_table.res",
                     "acba.007.p01.13_phage_int_table.res", "acba.007.p01.13_intI.res",
                     "acba.007.p01.13_phage_int.res"]
        exp_files = [os.path.join(self.out_dir, file) for file in exp_files]
        self.assertEqual(set(exp_files), set(files_created))
        # check proteins after annotation
        pdt.assert_frame_equal(proteins, integron1.proteins)

    def test_annot_multi(self):
        """
        Test func_annot when there are 4 integrons:
        - 1 calin with 4 proteins, 2 having a resfam annotation
        - 1 calin with 2 proteins, none having a resfam annotation
        - 1 in0
        - 1 complete with 4 proteins, 3 having a resfam annotation
        """
        # resfam pour: 16, 13, 3, 12
        #
        # Create integron in0
        integron1 = integron_finder.Integron(self.replicon_name)
        integron1.add_integrase(56, 1014, "ACBA.007.P01_13_1", 1, 1.9e-25,
                                      "intersection_tyr_intI")
        # Create integron CALIN with resfam proteins
        integron2 = integron_finder.Integron(self.replicon_name)
        integron2.add_attC(7400, 7650, -1, 7e-9, "attc_4")
        integron2.add_attC(8600, 8650, -1, 7e-4, "attc_4")
        integron2.add_attC(10200, 10400, -1, 7e-7, "attc_4")
        integron2.add_attC(10800, 10900, -1, 7e-7, "attc_4")
        integron2.add_proteins()

        # Create integron CALIN without any resfam proteins
        integron3 = integron_finder.Integron(self.replicon_name)
        integron3.add_attC(4320, 4400, -1, 7e-9, "attc_4")
        integron3.add_proteins()

        # Create complete integron
        integron4 = integron_finder.Integron(self.replicon_name)
        integron4.add_attC(17825, 17884, -1, 7e-9, "attc_4")
        integron4.add_attC(19080, 19149, -1, 7e-4, "attc_4")
        integron4.add_attC(19618, 19726, -1, 7e-7, "attc_4")
        integron4.add_integrase(16542, 17381, "ACBA.007.P01_13_19", -1, 1.9e-25,
                                      "intersection_tyr_intI")
        integron4.add_proteins()

        integron_finder.integrons = [integron1, integron2, integron3, integron4]

        # Create dataframes for expected proteins before annotation
        proteins1 = pd.DataFrame(columns=["pos_beg", "pos_end", "strand",
                                          "evalue", "type_elt", "model",
                                          "distance_2attC", "annotation"])
        proteins1 = proteins1.astype(dtype={"pos_beg": "int", "pos_end": "int", "strand": "int",
                                          "evalue": "float", "type_elt": "str", "model": "str",
                                          "distance_2attC": "float", "annotation": "str"})
        proteins1 = proteins1[["pos_beg", "pos_end", "strand", "evalue", "type_elt",
                               "model", "distance_2attC", "annotation"]]
        proteins2 = pd.DataFrame({"pos_beg": [7088, 7710, 8650, 10524],
                                 "pos_end": [7351, 8594, 10125, 11699],
                                 "strand": [1, -1, -1, -1],
                                 "evalue": [np.nan] * 4,
                                 "type_elt": ["protein"] * 4,
                                 "model": ["NA"] * 4,
                                 "distance_2attC": [np.nan] * 4,
                                 "annotation": ["protein"] * 4},
                                index=["ACBA.007.P01_13_11", "ACBA.007.P01_13_12",
                                       "ACBA.007.P01_13_13", "ACBA.007.P01_13_14"])
        proteins2 = proteins2[["pos_beg", "pos_end", "strand", "evalue", "type_elt",
                               "model", "distance_2attC", "annotation"]]
        proteins3 = pd.DataFrame({"pos_beg": [3546, 4380],
                                 "pos_end": [4313, 4721],
                                 "strand": [1, 1],
                                 "evalue": [np.nan] * 2,
                                 "type_elt": ["protein"] * 2,
                                 "model": ["NA"] * 2,
                                 "distance_2attC": [np.nan] * 2,
                                 "annotation": ["protein"] * 2},
                                index=["ACBA.007.P01_13_6", "ACBA.007.P01_13_7"])
        proteins3 = proteins3[["pos_beg", "pos_end", "strand", "evalue", "type_elt",
                               "model", "distance_2attC", "annotation"]]
        proteins4 = pd.DataFrame({"pos_beg": [17375, 17886, 19090, 19721],
                                 "pos_end": [17722, 18665, 19749, 20254],
                                 "strand": [-1] * 4,
                                 "evalue": [np.nan] * 4,
                                 "type_elt": ["protein"] * 4,
                                 "model": ["NA"] * 4,
                                 "distance_2attC": [np.nan] * 4,
                                 "annotation": ["protein"] * 4},
                                index=["ACBA.007.P01_13_20", "ACBA.007.P01_13_21",
                                       "ACBA.007.P01_13_22", "ACBA.007.P01_13_23"])
        proteins4 = proteins4[["pos_beg", "pos_end", "strand", "evalue", "type_elt",
                               "model", "distance_2attC", "annotation"]]

        # Check proteins before annotation
        int_proteins = [proteins1, proteins2, proteins3, proteins4]
        for inte, prots in zip(integron_finder.integrons, int_proteins):
            pdt.assert_frame_equal(inte.proteins, prots)

        # Annotate proteins with evalue threshold
        integron_finder.func_annot(self.replicon_name, self.out_dir, self.hmm_files, evalue=1e-32)

        # Check that all files generated are as expected
        files_created = glob.glob(os.path.join(self.out_dir, "*"))
        self.assertEqual(set(self.exp_files), set(files_created))

        # Check that annotated proteins are as expected
        proteins2.loc["ACBA.007.P01_13_13"] = [8650, 10125, -1, 2.4e-86, "protein",
                                              "RF0007", np.nan, "ABC_efflux"]
        proteins4.loc["ACBA.007.P01_13_21"] = [17886, 18665, -1, 7.4e-168, "protein",
                                              "RF0027", np.nan, "ANT3"]
        proteins4.loc["ACBA.007.P01_13_23"] = [19721, 20254, -1, 6.2e-110, "protein",
                                              "RF0003", np.nan, "AAC3-I"]
        for inte, prots in zip(integron_finder.integrons, int_proteins):
            pdt.assert_frame_equal(inte.proteins, prots)

        # Annotate proteins with default evalue (1 more annotation)
        integron_finder.func_annot(self.replicon_name, self.out_dir, self.hmm_files)
        proteins4.loc["ACBA.007.P01_13_20"] = [17375, 17722, -1, 4.5e-31, "protein",
                                              "RF0066", np.nan, "emrE"]
        for inte, prots in zip(integron_finder.integrons, int_proteins):
            pdt.assert_frame_equal(inte.proteins, prots)

        # Annotate proteins with lower coverage threshold (1 more annotation)
        integron_finder.func_annot(self.replicon_name, self.out_dir, self.hmm_files, coverage=0.4)

        proteins2.loc["ACBA.007.P01_13_12"] = [7710, 8594, -1, 1.6e-5, "protein",
                                              "RF0033", np.nan, "APH3"]
        for inte, prots in zip(integron_finder.integrons, int_proteins):
            pdt.assert_frame_equal(inte.proteins, prots)

    def test_annot_wrong_hmm(self):
        """
        Test that when the given hmm file does not exist, it returns an error specifying that
        the hmm command ended with a non-zero return code.
        """
        hmm_files = ["myhmm.hmm"]
        # Create integron
        integron1 = integron_finder.Integron(self.replicon_name)
        integron_finder.integrons = [integron1]
        # Add only attc sites (no integrase)
        integron1.add_attC(17825, 17884, -1, 7e-9, "attc_4")
        integron1.add_attC(19080, 19149, -1, 7e-4, "attc_4")
        integron1.add_attC(19618, 19726, -1, 7e-7, "attc_4")
        # Add proteins between attC sites
        integron1.add_proteins()
        # Annotate proteins
        with self.assertRaises(RuntimeError) as exp:
            integron_finder.func_annot(self.replicon_name, self.out_dir, hmm_files)

        raised = exp.exception
        self.assertEqual(raised.message, integron_finder.HMMSEARCH + " failed return code = 1")

    def test_annot_wrong_hmmsearch(self):
        """
        Test that when the given HMMSEARCH command does not exist, it raises an exception
        specifying that the given command could not run.
        """
        integron_finder.HMMSEARCH = "hmmsearchh"
        # Create integron
        integron1 = integron_finder.Integron(self.replicon_name)
        integron_finder.integrons = [integron1]
        # Add only attc sites (no integrase)
        integron1.add_attC(17825, 17884, -1, 7e-9, "attc_4")
        integron1.add_attC(19080, 19149, -1, 7e-4, "attc_4")
        integron1.add_attC(19618, 19726, -1, 7e-7, "attc_4")
        # Add proteins between attC sites
        integron1.add_proteins()
        # Annotate proteins
        with self.assertRaises(RuntimeError) as exp:
            integron_finder.func_annot(self.replicon_name, self.out_dir, self.hmm_files)
        raised = exp.exception
        self.assertEqual(raised.message, integron_finder.HMMSEARCH +\
                         " failed : [Errno 2] No such file or directory")
