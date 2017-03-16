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
from Bio import SeqFeature


class TestFunctions(unittest.TestCase):
    def setUp(self):
        """
        Define variables common to all tests
        """
        self.replicon_path = os.path.join("tests", "data", "acba.007.p01.13.fst")
        self.seq = SeqIO.read(self.replicon_path, "fasta",
                              alphabet=Seq.IUPAC.unambiguous_dna)
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
        infos = {"ID_replicon": "acba.007.p01.13",
                 "ID_integron": "integron_01",
                 "element": "ACBA.007.P01_13_20",
                 "pos_beg": 17375,
                 "pos_end": 17375,
                 "strand": -1,
                 "evalue": np.nan,
                 "type_elt": "protein",
                 "annotation": "protein",
                 "model": "NA",
                 "type": "complete",
                 "default": "Yes",
                 "distance_2attC": np.nan
                }

        df = pd.DataFrame(infos, index = [0])

        start_seq = self.seq.seq
        start_id = self.seq.id

        integron_finder.to_gbk(df, self.seq)

        # Translation should be protein ACBA.007.P01_13_20 in
        # tests/data/Results_Integron_Finder_acba.007.p01.13/acba.007.p01.13.prt
        translate = ("MKGWLFLVIAIVGEVIATSALKSSEGFTKLAPSAVVIIGYGIAFYFLSLVLKSIPVGVAY"
                     "AVWSGLGVVIITAIAWLLHGQKLDAWGFVGMGLIIAAFLLARSPSWKSLRRPTPW*")

        # Check that there are 2 features (integron and protein)
        self.assertEqual(len(self.seq.features), 2)
        # Check that initial sequence and id are not modified
        self.assertEqual(self.seq.seq, start_seq)
        self.assertEqual(self.seq.id, start_id)
        # Check first feature: integron
        self.assertEqual(self.seq.features[0].location.start, infos["pos_beg"] - 1)
        self.assertEqual(self.seq.features[0].location.end, infos["pos_end"])
        self.assertEqual(self.seq.features[0].strand, 0)
        self.assertEqual(self.seq.features[0].type, "integron")
        self.assertEqual(self.seq.features[0].qualifiers["integron_id"], infos["ID_integron"])
        self.assertEqual(self.seq.features[0].qualifiers["integron_type"], infos["type"])
        # Check second feature: protein
        self.assertEqual(self.seq.features[1].location.start, infos["pos_beg"] - 1)
        self.assertEqual(self.seq.features[1].location.end, infos["pos_end"])
        self.assertEqual(self.seq.features[1].strand, infos["strand"])
        self.assertEqual(self.seq.features[1].type, "CDS")
        self.assertEqual(self.seq.features[1].qualifiers["protein_id"], infos["element"])
        self.assertEqual(self.seq.features[1].qualifiers["gene"], infos["annotation"])
        self.assertEqual(self.seq.features[1].qualifiers["model"], infos["model"])
        self.assertEqual(self.seq.features[1].qualifiers["translation"].tostring(), translate)


    def test_integron_1elem_int(self):
        """
        Test to_gbk when the only element is an integron composed of 1 integrase only.

        """
        infos = {"ID_replicon": "acba.007.p01.13",
                 "ID_integron": "integron_01",
                 "element": "ACBA.007.P01_13_1",
                 "pos_beg": 55,
                 "pos_end": 1014,
                 "strand": 1,
                 "evalue": 1.9e-25,
                 "type_elt": "protein",
                 "annotation": "intI",
                 "model": "intersection_tyr_intI",
                 "type": "complete",
                 "default": "Yes",
                 "distance_2attC": np.nan
                }

        df = pd.DataFrame(infos, index = [0])

        start_seq = self.seq.seq
        start_id = self.seq.id

        integron_finder.to_gbk(df, self.seq)

        # Translation should be protein ACBA.007.P01_13_1 in
        # tests/data/Results_Integron_Finder_acba.007.p01.13/acba.007.p01.13.prt
        translate = ("MKTATAPLPPLRSVKVLDQLRERIRYLHYSLRTEQAYVNWVRAFIRFHGVRHPATLGSSE"
                     "VEAFLSWLANERKVSVSTHRQALAALLFFYGKVLCTDLPWLQEIGRPRPSRRLPVVLTPD"
                     "EVVRILGFLEGEHRLFAQLLYGTGMRISEGLQLRVKDLDFDHGTIIVREGKGSKDRALML"
                     "PESLAPSLREQLSRARAWWLKDQAEGRSGVALPDALERKYPRAGHSWPWFWVFAQHTHST"
                     "DPRSGVVRRHHMYDQTFQRAFKRAVEGTVAKLAMRQPFVLFKGLTFQKLCLPGAFRPGDH"
                     "HNKMLRPGLCVVHASPQYL*")

        # Check that there are 2 features (integron and protein)
        self.assertEqual(len(self.seq.features), 2)
        # Check that initial sequence and id are not modified
        self.assertEqual(self.seq.seq, start_seq)
        self.assertEqual(self.seq.id, start_id)
        # Check first feature: integron
        self.assertEqual(self.seq.features[0].location.start, infos["pos_beg"] - 1)
        self.assertEqual(self.seq.features[0].location.end, infos["pos_end"])
        self.assertEqual(self.seq.features[0].strand, 0)
        self.assertEqual(self.seq.features[0].type, "integron")
        self.assertEqual(self.seq.features[0].qualifiers["integron_id"], infos["ID_integron"])
        self.assertEqual(self.seq.features[0].qualifiers["integron_type"], infos["type"])
        # Check second feature: protein
        self.assertEqual(self.seq.features[1].location.start, infos["pos_beg"] - 1)
        self.assertEqual(self.seq.features[1].location.end, infos["pos_end"])
        self.assertEqual(self.seq.features[1].strand, infos["strand"])
        self.assertEqual(self.seq.features[1].type, "integrase")
        self.assertEqual(self.seq.features[1].qualifiers["protein_id"], infos["element"])
        self.assertEqual(self.seq.features[1].qualifiers["gene"], infos["annotation"])
        self.assertEqual(self.seq.features[1].qualifiers["model"], infos["model"])
        self.assertEqual(self.seq.features[1].qualifiers["translation"].tostring(), translate)


    def test_integron_1elem_prom(self):
        """
        Test to_gbk when the only element is an integron composed of 1 promoter only.

        """
        infos = {"ID_replicon": "acba.007.p01.13",
                 "ID_integron": "integron_01",
                 "element": "Pc_int1",
                 "pos_beg": 25,
                 "pos_end": 51,
                 "strand": -1,
                 "evalue": np.nan,
                 "type_elt": "Promoter",
                 "annotation": "Pc_1",
                 "model": "NA",
                 "type": "complete",
                 "default": "Yes",
                 "distance_2attC": np.nan
                }

        df = pd.DataFrame(infos, index = [0])

        start_seq = self.seq.seq
        start_id = self.seq.id

        integron_finder.to_gbk(df, self.seq)

        # Check that there are 2 features (integron and promoter)
        self.assertEqual(len(self.seq.features), 2)
        # Check that initial sequence and id are not modified
        self.assertEqual(self.seq.seq, start_seq)
        self.assertEqual(self.seq.id, start_id)
        # Check first feature: integron
        self.assertEqual(self.seq.features[0].location.start, infos["pos_beg"] - 1)
        self.assertEqual(self.seq.features[0].location.end, infos["pos_end"])
        self.assertEqual(self.seq.features[0].strand, 0)
        self.assertEqual(self.seq.features[0].type, "integron")
        self.assertEqual(self.seq.features[0].qualifiers["integron_id"], infos["ID_integron"])
        self.assertEqual(self.seq.features[0].qualifiers["integron_type"], infos["type"])
        # Check second feature: protein
        self.assertEqual(self.seq.features[1].location.start, infos["pos_beg"] - 1)
        self.assertEqual(self.seq.features[1].location.end, infos["pos_end"])
        self.assertEqual(self.seq.features[1].strand, infos["strand"])
        self.assertEqual(self.seq.features[1].type, infos["type_elt"])
        self.assertEqual(self.seq.features[1].qualifiers["Promoter"], infos["element"])
        self.assertEqual(self.seq.features[1].qualifiers["model"], infos["model"])


