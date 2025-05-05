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
import argparse
import pandas as pd
import numpy as np

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.annotation import add_feature
from integron_finder.utils import FastaIterator
from integron_finder.topology import Topology
from integron_finder.config import Config
from integron_finder.prot_db import ProdigalDB

class TestAddFeature(IntegronTest):

    def setUp(self):
        """
        Define variables common to all tests
        """
        self.replicon_path = self.find_data(os.path.join('Replicons', "acba.007.p01.13.fst"))
        self.replicon_id = 'ACBA.007.P01_13'
        topologies = Topology(1, 'lin')
        with FastaIterator(self.replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            self.seq = next(sequences_db)

        self.prot_file = self.find_data(os.path.join("Results_Integron_Finder_acba.007.p01.13",
                                                     "tmp_{}".format(self.replicon_id),
                                                     "{}.prt".format(self.replicon_id)))
        args = argparse.Namespace()
        args.replicon = self.replicon_path
        args.gembase = False
        args.prot_file = False
        args.cmsearch = __file__
        args.hmmsearch = __file__
        args.prodigal = __file__
        cfg = Config(args)
        self.prot_db = ProdigalDB(self.seq, cfg, prot_file=self.prot_file)
        self.dist_threshold = 4000


    def test_integron_1elem_prot(self):
        # Test add_feature when the only element is an integron composed of 1 protein only.

        infos = {"ID_replicon": self.replicon_id,
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

        df = pd.DataFrame(infos, index=[0])

        start_seq = self.seq.seq
        start_id = self.seq.id

        add_feature(self.seq, df, self.prot_db, self.dist_threshold)

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
        self.assertEqual(self.seq.features[0].location.strand, 0)
        self.assertEqual(self.seq.features[0].type, "integron")
        self.assertEqual(self.seq.features[0].qualifiers["integron_id"], infos["ID_integron"])
        self.assertEqual(self.seq.features[0].qualifiers["integron_type"], infos["type"])
        # Check second feature: protein
        self.assertEqual(self.seq.features[1].location.start, infos["pos_beg"] - 1)
        self.assertEqual(self.seq.features[1].location.end, infos["pos_end"])
        self.assertEqual(self.seq.features[1].location.strand, infos["strand"])
        self.assertEqual(self.seq.features[1].type, "CDS")
        self.assertEqual(self.seq.features[1].qualifiers["protein_id"], infos["element"])
        self.assertEqual(self.seq.features[1].qualifiers["gene"], infos["annotation"])
        self.assertEqual(self.seq.features[1].qualifiers["model"], infos["model"])
        self.assertEqual(str(self.seq.features[1].qualifiers["translation"]), translate)


    def test_integron_1elem_int(self):
        # Test add_feature when the only element is an integron composed of 1 integrase only.

        infos = {"ID_replicon": self.replicon_id,
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

        df = pd.DataFrame(infos, index=[0])

        start_seq = self.seq.seq
        start_id = self.seq.id

        add_feature(self.seq, df, self.prot_db, self.dist_threshold)

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
        self.assertEqual(self.seq.features[0].location.strand, 0)
        self.assertEqual(self.seq.features[0].type, "integron")
        self.assertEqual(self.seq.features[0].qualifiers["integron_id"], infos["ID_integron"])
        self.assertEqual(self.seq.features[0].qualifiers["integron_type"], infos["type"])
        # Check second feature: protein
        self.assertEqual(self.seq.features[1].location.start, infos["pos_beg"] - 1)
        self.assertEqual(self.seq.features[1].location.end, infos["pos_end"])
        self.assertEqual(self.seq.features[1].location.strand, infos["strand"])
        self.assertEqual(self.seq.features[1].type, "integrase")
        self.assertEqual(self.seq.features[1].qualifiers["protein_id"], infos["element"])
        self.assertEqual(self.seq.features[1].qualifiers["gene"], infos["annotation"])
        self.assertEqual(self.seq.features[1].qualifiers["model"], infos["model"])
        self.assertEqual(str(self.seq.features[1].qualifiers["translation"]), translate)


    def test_integron_1elem_prom(self):
        # Test add_feature when the only element is an integron composed of 1 promoter only.
        infos = {"ID_replicon": self.replicon_id,
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

        df = pd.DataFrame(infos, index=[0])

        start_seq = self.seq.seq
        start_id = self.seq.id

        add_feature(self.seq, df, self.prot_file, self.dist_threshold)

        # Check that there are 2 features (integron and promoter)
        self.assertEqual(len(self.seq.features), 2)
        # Check that initial sequence and id are not modified
        self.assertEqual(self.seq.seq, start_seq)
        self.assertEqual(self.seq.id, start_id)
        # Check first feature: integron
        self.assertEqual(self.seq.features[0].location.start, infos["pos_beg"] - 1)
        self.assertEqual(self.seq.features[0].location.end, infos["pos_end"])
        self.assertEqual(self.seq.features[0].location.strand, 0)
        self.assertEqual(self.seq.features[0].type, "integron")
        self.assertEqual(self.seq.features[0].qualifiers["integron_id"], infos["ID_integron"])
        self.assertEqual(self.seq.features[0].qualifiers["integron_type"], infos["type"])
        # Check second feature: promotor
        self.assertEqual(self.seq.features[1].location.start, infos["pos_beg"] - 1)
        self.assertEqual(self.seq.features[1].location.end, infos["pos_end"])
        self.assertEqual(self.seq.features[1].location.strand, infos["strand"])
        self.assertEqual(self.seq.features[1].type, infos["type_elt"])
        self.assertEqual(self.seq.features[1].qualifiers["Promoter"], infos["element"])
        self.assertEqual(self.seq.features[1].qualifiers["model"], infos["model"])


    def test_integron_nelem(self):
        # Test add_feature when there are several elements in the integron:
        # protein, integrase, promotor, attC

        # Integron is composed of: 1 promoter, 1 integrase, 1protein, 1attC
        # promoter and integrase are at the begining of the genome
        # others are at the end
        int_id = "integron_01"
        int_type = "complete"
        infos_prom = {"ID_replicon": self.replicon_id,
                      "ID_integron": int_id,
                      "element": "Pc_int1",
                      "pos_beg": 25,
                      "pos_end": 51,
                      "strand": -1,
                      "evalue": np.nan,
                      "type_elt": "Promoter",
                      "annotation": "Pc_1",
                      "model": "NA",
                      "type": int_type,
                      "default": "Yes",
                      "distance_2attC": np.nan
                      }
        infos_int = {"ID_replicon": self.replicon_id,
                     "ID_integron": int_id,
                     "element": "ACBA.007.P01_13_1",
                     "pos_beg": 55,
                     "pos_end": 1014,
                     "strand": 1,
                     "evalue": 1.9e-25,
                     "type_elt": "protein",
                     "annotation": "intI",
                     "model": "NA",
                     "type": int_type,
                     "default": "Yes",
                     "distance_2attC": np.nan
                     }
        infos_prot = {"ID_replicon": self.replicon_id,
                      "ID_integron": int_id,
                      "element": "ACBA.007.P01_13_20",
                      "pos_beg": 17375,
                      "pos_end": 17375,
                      "strand": -1,
                      "evalue": np.nan,
                      "type_elt": "protein",
                      "annotation": "protein",
                      "model": "intersection_tyr_intI",
                      "type": int_type,
                      "default": "Yes",
                      "distance_2attC": np.nan
                      }
        infos_attC = {"ID_replicon": self.replicon_id,
                      "ID_integron": int_id,
                      "element": "attc_001",
                      "pos_beg": 17825,
                      "pos_end": 17884,
                      "strand": -1,
                      "evalue": 1e-9,
                      "type_elt": "attC",
                      "annotation": "attC",
                      "model": "attc_4",
                      "type": int_type,
                      "default": "Yes",
                      "distance_2attC": np.nan
                      }

        df1 = pd.DataFrame(infos_prom, index=[0])
        df2 = pd.DataFrame(infos_int, index=[0])
        df3 = pd.DataFrame(infos_prot, index=[0])
        df4 = pd.DataFrame(infos_attC, index=[0])

        df = pd.concat([df1, df2, df3, df4])

        start_seq = self.seq.seq
        start_id = self.seq.id
        lenseq = len(self.seq)

        tr_int = ("MKTATAPLPPLRSVKVLDQLRERIRYLHYSLRTEQAYVNWVRAFIRFHGVRHPATLGSSE"
                  "VEAFLSWLANERKVSVSTHRQALAALLFFYGKVLCTDLPWLQEIGRPRPSRRLPVVLTPD"
                  "EVVRILGFLEGEHRLFAQLLYGTGMRISEGLQLRVKDLDFDHGTIIVREGKGSKDRALML"
                  "PESLAPSLREQLSRARAWWLKDQAEGRSGVALPDALERKYPRAGHSWPWFWVFAQHTHST"
                  "DPRSGVVRRHHMYDQTFQRAFKRAVEGTVAKLAMRQPFVLFKGLTFQKLCLPGAFRPGDH"
                  "HNKMLRPGLCVVHASPQYL*")
        tr_prot = ("MKGWLFLVIAIVGEVIATSALKSSEGFTKLAPSAVVIIGYGIAFYFLSLVLKSIPVGVAY"
                   "AVWSGLGVVIITAIAWLLHGQKLDAWGFVGMGLIIAAFLLARSPSWKSLRRPTPW*")

        add_feature(self.seq, df, self.prot_db, self.dist_threshold)

        # Check that there are 5 features (integron, promoter, integrase, protein, attC)
        self.assertEqual(len(self.seq.features), 5)
        # Check that initial sequence and id are not modified
        self.assertEqual(self.seq.seq, start_seq)
        self.assertEqual(self.seq.id, start_id)
        # Check first feature: integron
        self.assertEqual(self.seq.features[0].location.parts[0].start, infos_prot["pos_beg"] - 1)
        self.assertEqual(self.seq.features[0].location.parts[0].end,  lenseq)
        self.assertEqual(self.seq.features[0].location.parts[1].start, 0)
        self.assertEqual(self.seq.features[0].location.parts[1].end, infos_int["pos_end"])
        self.assertEqual(self.seq.features[0].location.strand, 0)
        self.assertEqual(self.seq.features[0].type, "integron")
        self.assertEqual(self.seq.features[0].qualifiers["integron_id"], int_id)
        self.assertEqual(self.seq.features[0].qualifiers["integron_type"], int_type)
        # Check second feature: promoter
        self.assertEqual(self.seq.features[1].location.start, infos_prom["pos_beg"] - 1)
        self.assertEqual(self.seq.features[1].location.end, infos_prom["pos_end"])
        self.assertEqual(self.seq.features[1].location.strand, infos_prom["strand"])
        self.assertEqual(self.seq.features[1].type, "Promoter")
        self.assertEqual(self.seq.features[1].qualifiers["Promoter"], infos_prom["element"])
        self.assertEqual(self.seq.features[1].qualifiers["model"], infos_prom["model"])
        # Check second feature: integrase
        self.assertEqual(self.seq.features[2].location.start, infos_int["pos_beg"] - 1)
        self.assertEqual(self.seq.features[2].location.end, infos_int["pos_end"])
        self.assertEqual(self.seq.features[2].location.strand, infos_int["strand"])
        self.assertEqual(self.seq.features[2].type, "integrase")
        self.assertEqual(self.seq.features[2].qualifiers["protein_id"], infos_int["element"])
        self.assertEqual(self.seq.features[2].qualifiers["gene"], infos_int["annotation"])
        self.assertEqual(self.seq.features[2].qualifiers["model"], infos_int["model"])
        self.assertEqual(str(self.seq.features[2].qualifiers["translation"]), tr_int)
        # Check second feature: protein
        self.assertEqual(self.seq.features[3].location.start, infos_prot["pos_beg"] - 1)
        self.assertEqual(self.seq.features[3].location.end, infos_prot["pos_end"])
        self.assertEqual(self.seq.features[3].location.strand, infos_prot["strand"])
        self.assertEqual(self.seq.features[3].type, "CDS")
        self.assertEqual(self.seq.features[3].qualifiers["protein_id"], infos_prot["element"])
        self.assertEqual(self.seq.features[3].qualifiers["gene"], infos_prot["annotation"])
        self.assertEqual(self.seq.features[3].qualifiers["model"], infos_prot["model"])
        self.assertEqual(str(self.seq.features[3].qualifiers["translation"]), tr_prot)
        # Check second feature: attC
        self.assertEqual(self.seq.features[4].location.start, infos_attC["pos_beg"] - 1)
        self.assertEqual(self.seq.features[4].location.end, infos_attC["pos_end"])
        self.assertEqual(self.seq.features[4].location.strand, infos_attC["strand"])
        self.assertEqual(self.seq.features[4].type, "attC")
        self.assertEqual(self.seq.features[4].qualifiers["attC"], infos_attC["element"])
        self.assertEqual(self.seq.features[4].qualifiers["model"], infos_attC["model"])


    def test_integron_2int_nelem(self):
        # Test add_feature when there are 2 integrons:
        #     integron 1 with several elements: protein, integrase, promoter
        #     integron 2 with only 1 attC site
        # Integrons are not over the edge of sequence

        # integron 1
        int_id = "integron_01"
        int_type = "complete"
        infos_prom = {"ID_replicon": self.replicon_id,
                      "ID_integron": int_id,
                      "element": "Pc_int1",
                      "pos_beg": 25,
                      "pos_end": 51,
                      "strand": -1,
                      "evalue": np.nan,
                      "type_elt": "Promoter",
                      "annotation": "Pc_1",
                      "model": "NA",
                      "type": int_type,
                      "default": "Yes",
                      "distance_2attC": np.nan
                      }
        infos_int = {"ID_replicon": self.replicon_id,
                     "ID_integron": int_id,
                     "element": "ACBA.007.P01_13_1",
                     "pos_beg": 55,
                     "pos_end": 1014,
                     "strand": 1,
                     "evalue": 1.9e-25,
                     "type_elt": "protein",
                     "annotation": "intI",
                     "model": "NA",
                     "type": int_type,
                     "default": "Yes",
                     "distance_2attC": np.nan
                     }
        infos_prot = {"ID_replicon": self.replicon_id,
                      "ID_integron": int_id,
                      "element": "ACBA.007.P01_13_20",
                      "pos_beg": 2000,
                      "pos_end": 2056,
                      "strand": -1,
                      "evalue": np.nan,
                      "type_elt": "protein",
                      "annotation": "protein",
                      "model": "intersection_tyr_intI",
                      "type": int_type,
                      "default": "Yes",
                      "distance_2attC": np.nan
                      }
        # integron 2
        infos_attC = {"ID_replicon": self.replicon_id,
                      "ID_integron": "integron_02",
                      "element": "attc_001",
                      "pos_beg": 17825,
                      "pos_end": 17884,
                      "strand": -1,
                      "evalue": 1e-9,
                      "type_elt": "attC",
                      "annotation": "attC",
                      "model": "attc_4",
                      "type": int_type,
                      "default": "Yes",
                      "distance_2attC": np.nan
                      }

        df1 = pd.DataFrame(infos_prom, index=[0])
        df2 = pd.DataFrame(infos_int, index=[0])
        df3 = pd.DataFrame(infos_prot, index=[0])
        df4 = pd.DataFrame(infos_attC, index=[0])

        df = pd.concat([df1, df2, df3, df4])

        start_seq = self.seq.seq
        start_id = self.seq.id

        tr_int = ("MKTATAPLPPLRSVKVLDQLRERIRYLHYSLRTEQAYVNWVRAFIRFHGVRHPATLGSSE"
                  "VEAFLSWLANERKVSVSTHRQALAALLFFYGKVLCTDLPWLQEIGRPRPSRRLPVVLTPD"
                  "EVVRILGFLEGEHRLFAQLLYGTGMRISEGLQLRVKDLDFDHGTIIVREGKGSKDRALML"
                  "PESLAPSLREQLSRARAWWLKDQAEGRSGVALPDALERKYPRAGHSWPWFWVFAQHTHST"
                  "DPRSGVVRRHHMYDQTFQRAFKRAVEGTVAKLAMRQPFVLFKGLTFQKLCLPGAFRPGDH"
                  "HNKMLRPGLCVVHASPQYL*")
        tr_prot = ("MKGWLFLVIAIVGEVIATSALKSSEGFTKLAPSAVVIIGYGIAFYFLSLVLKSIPVGVAY"
                   "AVWSGLGVVIITAIAWLLHGQKLDAWGFVGMGLIIAAFLLARSPSWKSLRRPTPW*")

        add_feature(self.seq, df, self.prot_db, self.dist_threshold)

        # Check that there are 6 features (integron1, promoter, integrase, protein,
        #                                  integron2, attC)
        self.assertEqual(len(self.seq.features), 6)
        # Check that initial sequence and id are not modified
        self.assertEqual(self.seq.seq, start_seq)
        self.assertEqual(self.seq.id, start_id)
        # Check first feature: integron1
        self.assertEqual(self.seq.features[0].location.start, infos_prom["pos_beg"] - 1)
        self.assertEqual(self.seq.features[0].location.end, infos_prot["pos_end"])
        self.assertEqual(self.seq.features[0].location.strand, 0)
        self.assertEqual(self.seq.features[0].type, "integron")
        self.assertEqual(self.seq.features[0].qualifiers["integron_id"], int_id)
        self.assertEqual(self.seq.features[0].qualifiers["integron_type"], int_type)
        # Check feature 2: promoter
        self.assertEqual(self.seq.features[1].location.start, infos_prom["pos_beg"] - 1)
        self.assertEqual(self.seq.features[1].location.end, infos_prom["pos_end"])
        self.assertEqual(self.seq.features[1].location.strand, infos_prom["strand"])
        self.assertEqual(self.seq.features[1].type, "Promoter")
        self.assertEqual(self.seq.features[1].qualifiers["Promoter"], infos_prom["element"])
        self.assertEqual(self.seq.features[1].qualifiers["model"], infos_prom["model"])
        # Check feature 3: integrase
        self.assertEqual(self.seq.features[2].location.start, infos_int["pos_beg"] - 1)
        self.assertEqual(self.seq.features[2].location.end, infos_int["pos_end"])
        self.assertEqual(self.seq.features[2].location.strand, infos_int["strand"])
        self.assertEqual(self.seq.features[2].type, "integrase")
        self.assertEqual(self.seq.features[2].qualifiers["protein_id"], infos_int["element"])
        self.assertEqual(self.seq.features[2].qualifiers["gene"], infos_int["annotation"])
        self.assertEqual(self.seq.features[2].qualifiers["model"], infos_int["model"])
        self.assertEqual(str(self.seq.features[2].qualifiers["translation"]), tr_int)
        # Check feature 4: protein
        self.assertEqual(self.seq.features[3].location.start, infos_prot["pos_beg"] - 1)
        self.assertEqual(self.seq.features[3].location.end, infos_prot["pos_end"])
        self.assertEqual(self.seq.features[3].location.strand, infos_prot["strand"])
        self.assertEqual(self.seq.features[3].type, "CDS")
        self.assertEqual(self.seq.features[3].qualifiers["protein_id"], infos_prot["element"])
        self.assertEqual(self.seq.features[3].qualifiers["gene"], infos_prot["annotation"])
        self.assertEqual(self.seq.features[3].qualifiers["model"], infos_prot["model"])
        self.assertEqual(str(self.seq.features[3].qualifiers["translation"]), tr_prot)
        # Check feature 5: integron2
        self.assertEqual(self.seq.features[4].location.start, infos_attC["pos_beg"] - 1)
        self.assertEqual(self.seq.features[4].location.end, infos_attC["pos_end"])
        self.assertEqual(self.seq.features[4].location.strand, 0)
        self.assertEqual(self.seq.features[4].type, "integron")
        self.assertEqual(self.seq.features[4].qualifiers["integron_id"], "integron_02")
        self.assertEqual(self.seq.features[4].qualifiers["integron_type"], int_type)
        # Check feature 6: attC
        self.assertEqual(self.seq.features[5].location.start, infos_attC["pos_beg"] - 1)
        self.assertEqual(self.seq.features[5].location.end, infos_attC["pos_end"])
        self.assertEqual(self.seq.features[5].location.strand, infos_attC["strand"])
        self.assertEqual(self.seq.features[5].type, "attC")
        self.assertEqual(self.seq.features[5].qualifiers["attC"], infos_attC["element"])
        self.assertEqual(self.seq.features[5].qualifiers["model"], infos_attC["model"])


    def test_integron_long_seqname(self):
        # Test add_feature when the only element is an integron composed of 1 protein only.

        infos = {"ID_replicon": self.replicon_id,
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

        df = pd.DataFrame(infos, index=[0])

        start_seq = self.seq.seq
        start_id = self.seq.id
        seq_name = self.seq.name
        self.seq.name = "abcdefgh" + seq_name

        add_feature(self.seq, df, self.prot_db, self.dist_threshold)

        # Translation should be protein ACBA.007.P01_13_20 in
        # tests/data/Results_Integron_Finder_acba.007.p01.13/acba.007.p01.13.prt
        # translate = ("MKGWLFLVIAIVGEVIATSALKSSEGFTKLAPSAVVIIGYGIAFYFLSLVLKSIPVGVAY"
        #             "AVWSGLGVVIITAIAWLLHGQKLDAWGFVGMGLIIAAFLLARSPSWKSLRRPTPW*")

        # Check that there are 2 features (integron and protein)
        self.assertEqual(len(self.seq.features), 2)
        # Check that initial sequence and id are not modified
        self.assertEqual(self.seq.seq, start_seq)
        self.assertEqual(self.seq.id, start_id)
        # Check that sequence name has been shortened
        self.assertEqual(self.seq.name, "h" + seq_name)
