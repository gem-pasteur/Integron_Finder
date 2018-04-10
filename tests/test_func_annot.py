# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright © 2015 - 2018  Institut Pasteur, Paris.                                #
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
import distutils.spawn
import shutil
import argparse
import glob
import tempfile
import re

import numpy as np
import pandas as pd
import pandas.util.testing as pdt

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.config import Config
from integron_finder.topology import Topology
from integron_finder.utils import FastaIterator
from integron_finder.integrase import find_integrase
from integron_finder.integron import Integron
from integron_finder.annotation import func_annot
from integron_finder import annotation

_annot_call_ori = annotation.call

class TestFuncAnnot(IntegronTest):

    @classmethod
    def setUpClass(cls):
        # mute the logs
        cls.set_log_level('WARNING')

    @classmethod
    def tearDownClass(cls):
        # restore the log level for the other tests
        # loggers are singleton
        cls.set_log_level('INFO')

    def setUp(self):
        """
        Define variables common to all tests
        """
        replicon_name = "acba.007.p01.13"
        self.replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        sequences_db = FastaIterator(self.replicon_path)
        topologies = Topology('lin')
        sequences_db.topologies = topologies
        self.replicon = next(sequences_db)

        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)

        self.hmm_files = [self.find_data(os.path.join("Functional_annotation", "Resfams.hmm"))]

        # Define integron_finder variables
        args = argparse.Namespace()
        args.gembase = False
        args.hmmsearch = distutils.spawn.find_executable("hmmsearch")
        args.prodigal = distutils.spawn.find_executable("prodigal")
        args.cpu = 1
        args.out_dir = self.tmp_dir
        self.cfg = Config(args)
        self.cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

        prot_dir = os.path.join(self.tmp_dir, 'Proteins')
        os.makedirs(prot_dir)
        self.prot_file = os.path.join(prot_dir, self.replicon.name + ".prt")
        shutil.copyfile(self.find_data(os.path.join('Proteins', self.replicon.id + ".prt")), self.prot_file)

        self.exp_files = ["{}{}".format(self.replicon.id, suffix) for suffix in (".prt",
                                                                                 "_Resfams_fa_table.res",
                                                                                 "_intI_table.res",
                                                                                 "_phage_int_table.res",
                                                                                 "_Resfams_fa.res",
                                                                                 "_intI.res",
                                                                                 "_phage_int.res",
                                                                                 "_subseqprot.tmp")]
        self.exp_files = [os.path.join(self.tmp_dir, file) for file in self.exp_files]

        self.prot_dtype = {"pos_beg": 'int',
                           "pos_end": 'int',
                           "strand": 'int',
                           "evalue": 'float',
                           "type_elt": 'str',
                           "annotation": 'str',
                           "model": 'str',
                           "distance_2attC": 'float'}

        # Run prodigal to find CDS on replicon (and run hmmsearch on integrase (2 profiles))
        self.integrases = find_integrase(self.replicon_path, self.replicon, self.prot_file, self.tmp_dir, self.cfg)
        annotation.call = self.mute_call(_annot_call_ori)

    def tearDown(self):
        """
        To do after each test. remove output directory if it was generated
        """
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
            pass
        annotation.call = _annot_call_ori

    def test_annot_calin(self):
        """
        Test func_annot when the integron is a CALIN (attC but no integrase), with 4 proteins:
        for 3 of them resfam annotations are found, and not for the last 1.
        """
        # Create integron
        integron1 = Integron(self.replicon, self.cfg)
        integrons = [integron1]
        # Add only attc sites (no integrase)
        integron1.add_attC(17825, 17884, -1, 7e-9, "attc_4")
        integron1.add_attC(19080, 19149, -1, 7e-4, "attc_4")
        integron1.add_attC(19618, 19726, -1, 7e-7, "attc_4")
        # Add proteins between attC sites
        integron1.add_proteins(self.prot_file)
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

        func_annot(integrons, self.replicon, self.prot_file, self.hmm_files, self.cfg, self.tmp_dir)

        # Check that all files generated are as expected
        files_created = [f for f in glob.glob(os.path.join(self.tmp_dir, "*")) if os.path.isfile(f)]
        self.assertEqual(set(self.exp_files), set(files_created))

        # Check that annotated proteins are as expected
        proteins.loc["ACBA.007.P01_13_20"] = [17375, 17722, -1, 4.5e-31, "protein", "RF0066", np.nan, "emrE"]
        proteins.loc["ACBA.007.P01_13_21"] = [17886, 18665, -1, 7.4e-168, "protein", "RF0027", np.nan, "ANT3"]
        proteins.loc["ACBA.007.P01_13_23"] = [19721, 20254, -1, 6.2e-110, "protein", "RF0003", np.nan, "AAC3-I"]
        pdt.assert_frame_equal(proteins, integron1.proteins)


    def test_annot_calin_empty(self):
        """
        Test func_annot when the integron is a CALIN (attC but no integrase), without any protein:
        nothing to annotate
        """
        # Create integron
        integron1 = Integron(self.replicon, self.cfg)
        integrons = [integron1]
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
        func_annot(integrons, self.replicon, self.prot_file, self.hmm_files, self.cfg, self.tmp_dir)

        # Check that all files generated are as expected
        files_created = [f for f in glob.glob(os.path.join(self.tmp_dir, "*")) if os.path.isfile(f)]
        exp_files = ["{}{}".format(self.replicon.id, suffix) for suffix in (".prt",
                                                                            "_intI_table.res",
                                                                            "_phage_int_table.res",
                                                                            "_intI.res",
                                                                            "_phage_int.res")]
        exp_files = [os.path.join(self.tmp_dir, file) for file in exp_files]
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
        open(os.path.join(self.tmp_dir, "{}_subseqprot.tmp".format(self.replicon.id)), "w").close()
        # Create integron
        integron1 = Integron(self.replicon, self.cfg)
        integrons = [integron1]
        # Add integrase
        integron1.add_integrase(55, 1014, "ACBA.007.P01_13_1", 1, 1.9e-25, "intersection_tyr_intI")
        # check proteins before annotation
        proteins = pd.DataFrame(columns=["pos_beg", "pos_end", "strand",
                                         "evalue", "type_elt", "model",
                                         "distance_2attC", "annotation"])
        proteins = proteins.astype(dtype={"pos_beg": "int", "pos_end": "int", "strand": "int",
                                          "evalue": "float", "type_elt": "str", "model": "str",
                                          "distance_2attC": "float", "annotation": "str"})
        pdt.assert_frame_equal(proteins, integron1.proteins)

        # Annotate proteins
        func_annot(integrons, self.replicon, self.prot_file, self.hmm_files, self.cfg, self.tmp_dir)
        # Check that all files generated are as expected
        files_created = [f for f in glob.glob(os.path.join(self.tmp_dir, "*")) if os.path.isfile(f)]
        exp_files = ["{}{}".format(self.replicon.id, suffix) for suffix in (".prt",
                                                                            "_intI_table.res",
                                                                            "_phage_int_table.res",
                                                                            "_intI.res",
                                                                            "_phage_int.res")]
        exp_files = [os.path.join(self.tmp_dir, file) for file in exp_files]
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

        # Create integron in0
        integron1 = Integron(self.replicon.name, self.cfg)
        integron1.add_integrase(56, 1014, "ACBA.007.P01_13_1", 1, 1.9e-25, "intersection_tyr_intI")

        # Create integron CALIN with resfam proteins
        integron2 = Integron(self.replicon, self.cfg)
        integron2.add_attC(7400, 7650, -1, 7e-9, "attc_4")
        integron2.add_attC(8600, 8650, -1, 7e-4, "attc_4")
        integron2.add_attC(10200, 10400, -1, 7e-7, "attc_4")
        integron2.add_attC(10800, 10900, -1, 7e-7, "attc_4")
        integron2.add_proteins(self.prot_file)

        # Create integron CALIN without any resfam proteins
        integron3 = Integron(self.replicon, self.cfg)
        integron3.add_attC(4320, 4400, -1, 7e-9, "attc_4")
        integron3.add_proteins(self.prot_file)

        # Create complete integron
        integron4 = Integron(self.replicon, self.cfg)
        integron4.add_attC(17825, 17884, -1, 7e-9, "attc_4")
        integron4.add_attC(19080, 19149, -1, 7e-4, "attc_4")
        integron4.add_attC(19618, 19726, -1, 7e-7, "attc_4")
        integron4.add_integrase(16542, 17381, "ACBA.007.P01_13_19", -1, 1.9e-25, "intersection_tyr_intI")
        integron4.add_proteins(self.prot_file)

        integrons = [integron1, integron2, integron3, integron4]

        # Create dataframes for expected proteins before annotation
        proteins1 = pd.DataFrame(columns=["pos_beg", "pos_end", "strand",
                                          "evalue", "type_elt", "model",
                                          "distance_2attC", "annotation"])
        proteins1 = proteins1.astype(dtype={"pos_beg": "int", "pos_end": "int", "strand": "int",
                                            "evalue": "float", "type_elt": "str", "model": "str",
                                            "distance_2attC": "float", "annotation": "str"})
        proteins1 = proteins1[["pos_beg", "pos_end", "strand", "evalue", "type_elt",
                               "model", "distance_2attC", "annotation"]]
        proteins1 = proteins1.astype(dtype=self.prot_dtype)

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
        proteins2 = proteins2.astype(dtype=self.prot_dtype)

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
        proteins3 = proteins3.astype(dtype=self.prot_dtype)

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
        proteins4 = proteins4.astype(dtype=self.prot_dtype)

        # Check proteins before annotation
        expected_proteins = [proteins1, proteins2, proteins3, proteins4]

        for inte, exp_prot in zip(integrons, expected_proteins):
            pdt.assert_frame_equal(inte.proteins, exp_prot)

        # Annotate proteins with evalue threshold
        func_annot(integrons, self.replicon, self.prot_file, self.hmm_files, self.cfg, self.tmp_dir, evalue=1e-32)

        # Check that all files generated are as expected
        files_created = [f for f in glob.glob(os.path.join(self.tmp_dir, "*")) if os.path.isfile(f)]
        self.assertEqual(set(self.exp_files), set(files_created))

        # Check that annotated proteins are as expected
        proteins2.loc["ACBA.007.P01_13_13"] = [8650, 10125, -1, 2.4e-86, "protein",
                                               "RF0007", np.nan, "ABC_efflux"]
        proteins4.loc["ACBA.007.P01_13_21"] = [17886, 18665, -1, 7.4e-168, "protein",
                                               "RF0027", np.nan, "ANT3"]
        proteins4.loc["ACBA.007.P01_13_23"] = [19721, 20254, -1, 6.2e-110, "protein",
                                               "RF0003", np.nan, "AAC3-I"]
        for inte, prots in zip(integrons, expected_proteins):
            pdt.assert_frame_equal(inte.proteins, prots)

        # Annotate proteins with default evalue (1 more annotation)
        with self.catch_io(out=True):
            func_annot(integrons, self.replicon, self.prot_file, self.hmm_files, self.cfg, self.tmp_dir)
        proteins4.loc["ACBA.007.P01_13_20"] = [17375, 17722, -1, 4.5e-31, "protein",
                                               "RF0066", np.nan, "emrE"]
        for inte, prots in zip(integrons, expected_proteins):
            pdt.assert_frame_equal(inte.proteins, prots)

        # Annotate proteins with lower coverage threshold (1 more annotation)
        with self.catch_io(out=True):
            func_annot(integrons, self.replicon, self.prot_file, self.hmm_files, self.cfg, self.tmp_dir, coverage=0.4)

        proteins2.loc["ACBA.007.P01_13_12"] = [7710, 8594, -1, 1.6e-5, "protein",
                                               "RF0033", np.nan, "APH3"]
        for inte, prots in zip(integrons, expected_proteins):
            pdt.assert_frame_equal(inte.proteins, prots)


    def test_annot_wrong_hmm(self):
        """
        Test that when the given hmm file does not exist, it returns an error specifying that
        the hmm command ended with a non-zero return code.
        """
        wrong_hmm_files = ["myhmm.hmm"]
        # Create integron
        integron1 = Integron(self.replicon, self.cfg)
        integrons = [integron1]
        # Add only attc sites (no integrase)
        integron1.add_attC(17825, 17884, -1, 7e-9, "attc_4")
        integron1.add_attC(19080, 19149, -1, 7e-4, "attc_4")
        integron1.add_attC(19618, 19726, -1, 7e-7, "attc_4")
        # Add proteins between attC sites
        integron1.add_proteins(self.prot_file)
        # Annotate proteins
        with self.assertRaises(RuntimeError) as ctx:
            func_annot(integrons, self.replicon, self.prot_file, wrong_hmm_files, self.cfg, self.tmp_dir)
        self.assertTrue(str(ctx.exception).endswith(" failed return code = 1"))


    def test_annot_wrong_hmmsearch(self):
        """
        Test that when the given HMMSEARCH command does not exist, it raises an exception
        specifying that the given command could not run.
        """
        self.cfg._args.hmmsearch = "nimportnaoik"
        # Create integron
        integron1 = Integron(self.replicon.name, self.cfg)
        integrons = [integron1]
        # Add only attc sites (no integrase)
        integron1.add_attC(17825, 17884, -1, 7e-9, "attc_4")
        integron1.add_attC(19080, 19149, -1, 7e-4, "attc_4")
        integron1.add_attC(19618, 19726, -1, 7e-7, "attc_4")
        # Add proteins between attC sites
        integron1.add_proteins(self.prot_file)
        # Annotate proteins
        with self.assertRaises(RuntimeError) as ctx:
            func_annot(integrons, self.replicon, self.prot_file, self.hmm_files, self.cfg, self.tmp_dir)
        self.assertTrue(re.search("failed : \[Errno 2\] No such file or directory: 'nimportnaoik'", str(ctx.exception)))
