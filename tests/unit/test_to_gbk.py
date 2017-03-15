#!/usr/bin/env python
# coding: utf-8

"""
Unit tests func_annot function of integron_finder
"""

import integron_finder
import unittest
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
from Bio import Seq


class TestFunctions(unittest.TestCase):
    def setUp(self):
        """
        Define variables common to all tests
        """
        self.replicon_path = os.path.join("tests", "data", "acba.007.p01.13.fst")
        integron_finder.PROT_file = os.path.join("tests", "data",
                                                 "Results_Integron_Finder_acba.007.p01.13",
                                                 "acba.007.p01.13.prt")

    def tearDown(self):
        """
        To do after each test. remove output directory if it was generated
        """
        pass
    # def test_blabla(self):
    #     """
    #     Test blabla
    #     """
    #     # to_gbk(df, sequence)
    #     # sequence is     SEQUENCE = SeqIO.read(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna) -> sequence du replicon donné en entrée
    #     # df is :
    #     # pour chaque integron, tous les éléments de l'integron (integrase, attC, attI, promoteur, proteins)
    #     associate an integron number to each integron, sorted by start position of all proteins -> replace ID_integron
    #     reorder columns, and put evalue column to float type
    #     sort all elements by ID_integron, pos_beg, evalue

    #     first case : integron contains only 1 element

    def test_integron_1elem_prot(self):
        """
        Test to_gbk when the only element is an integron composed of 1 protein only.

        """
        df = pd.DataFrame({"ID_replicon": "acba.007.p01.13", "element": "ACBA.007.P01_13_1",
                           "pos_beg": 55, "pos_end": 1014, "strand": -1,
                           "evalue": 1.9e-25, "type_elt": "protein", "annotation": "protein",
                           "model": "NA", "type": "unknown", "default": "Yes",
                           "distance_2attC": np.nan, "ID_integron": "integron_01"},
                          index = [0])
        seq = SeqIO.read(self.replicon_path, "fasta",
                         alphabet=Seq.IUPAC.unambiguous_dna)

        integron_finder.to_gbk(df, seq)
        # import IPython
        # IPython.embed()
        self.assertEqual(len(seq.features), 2)
        # check each feature
        # print seq
        # print df.shape



