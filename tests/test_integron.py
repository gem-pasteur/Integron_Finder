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
import pandas.testing as pdt
import numpy as np
from Bio import Seq
from Bio.SeqRecord import SeqRecord

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.config import Config
from integron_finder.utils import FastaIterator
from integron_finder.topology import Topology
from integron_finder.integron import Integron
from integron_finder.prot_db import ProdigalDB

class TestIntegron(IntegronTest):


    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

        self.columns = ['pos_beg', 'pos_end', 'strand', 'evalue', 'type_elt', 'model', 'distance_2attC', 'annotation']
        self.dtype = {"pos_beg": 'int',
                      "pos_end": 'int',
                      "strand": 'int',
                      "evalue": 'float',
                      "type_elt": 'str',
                      "annotation": 'str',
                      "model": 'str',
                      "distance_2attC": 'float'}
        args = argparse.Namespace()
        args.gembase = False
        args.prot_file = False
        args.cmsearch = __file__
        args.hmmsearch = __file__
        args.prodigal = __file__
        self.cfg = Config(args)
        self.cfg._args.eagle_eyes = False
        self.cfg._args.local_max = False
        self._prefix_data = os.path.join(os.path.dirname(__file__), 'data')


    def test_add_integrase(self):
        replicon_name = "acba.007.p01.13"
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        topologies = Topology(2, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        data_integrase = {"pos_beg": 55,
                          "pos_end": 1014,
                          "strand": 1,
                          "evalue": 1.900000e-25,
                          "type_elt": "protein",
                          "annotation": "intI",
                          "model": "intersection_tyr_intI",
                          "distance_2attC": np.nan}
        id_int = "ACBA.007.P01_13_1"

        df = pd.DataFrame(data_integrase,
                          columns=self.columns,
                          index=[id_int])
        df = df.astype(dtype=self.dtype)

        integron = Integron(replicon, self.cfg)
        integron.add_integrase(data_integrase["pos_beg"],
                               data_integrase["pos_end"],
                               id_int,
                               data_integrase["strand"],
                               data_integrase["evalue"],
                               data_integrase["model"]
                               )
        pdt.assert_frame_equal(df, integron.integrase)

        with self.assertRaises(RuntimeError) as ctx:
            integron.add_integrase(data_integrase["pos_beg"],
                                   data_integrase["pos_end"],
                                   id_int,
                                   data_integrase["strand"],
                                   data_integrase["evalue"],
                                   data_integrase["model"]
                                   )
        self.assertEqual(str(ctx.exception), "add_integrase should be called once.")


    def test_add_attc(self):
        replicon_name = "acba.007.p01.13"
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        topologies = Topology(2, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        data_attc_1 = {"pos_beg": 10,
                       "pos_end": 100,
                       "strand": -1,
                       "evalue": 1.1e-07,
                       "type_elt": "attC",
                       "annotation": "attC",
                       "model": "attc_4",
                       "distance_2attC": np.nan}

        attc_1 = pd.DataFrame(data_attc_1,
                              columns=self.columns,
                              index=['attc_001'])
        attc_1 = attc_1.astype(dtype=self.dtype)

        integron = Integron(replicon, self.cfg)
        integron.add_attC(attc_1.loc['attc_001', 'pos_beg'],
                          attc_1.loc['attc_001', 'pos_end'],
                          attc_1.loc['attc_001', 'strand'],
                          attc_1.loc['attc_001', 'evalue'],
                          attc_1.loc['attc_001', 'model'])

        pdt.assert_frame_equal(attc_1, integron.attC)

        attc_2 = pd.DataFrame(data_attc_1,
                              columns=self.columns,
                              index=['attc_002'])
        attc_2 = attc_2.astype(dtype=self.dtype)
        attc_2['pos_beg'] = attc_2['pos_beg'] + 100
        attc_2['pos_end'] = attc_2['pos_end'] + 100
        attc_2["distance_2attC"] = (attc_2.loc['attc_002', 'pos_beg'] - attc_1.loc['attc_001', 'pos_end']) % len(replicon)

        attc = pd.concat((attc_1, attc_2))

        integron.add_attC(attc_2.loc['attc_002', 'pos_beg'],
                          attc_2.loc['attc_002', 'pos_end'],
                          attc_2.loc['attc_002', 'strand'],
                          attc_2.loc['attc_002', 'evalue'],
                          attc_2.loc['attc_002', 'model'])
        pdt.assert_frame_equal(attc, integron.attC)


    def test_type(self):
        replicon = SeqRecord(Seq.Seq(''), id='foo')
        no_integrase = Integron(replicon, self.cfg)
        self.assertIsNone(no_integrase.type())

        replicon = SeqRecord(Seq.Seq(''), id='just_one_integrase')
        just_one_integrase = Integron(replicon, self.cfg)
        just_one_integrase.add_integrase(10,
                                         100,
                                         'foo',
                                         1,
                                         1e-2,
                                         "intersection_tyr_intI")
        self.assertEqual(just_one_integrase.type(), "In0")

        replicon = SeqRecord(Seq.Seq(''), id='just_one_attC')
        just_one_attC = Integron(replicon, self.cfg)
        just_one_attC.add_attC(10,
                               100,
                               1,
                               1e-2,
                               "intersection_tyr_intI")
        self.assertEqual(just_one_attC.type(), "CALIN")

        replicon = SeqRecord(Seq.Seq(''), id='one_integrase_one_attC')
        one_integrase_one_attC = Integron(replicon, self.cfg)
        one_integrase_one_attC.add_integrase(10,
                                             100,
                                             'foo',
                                             1,
                                             1e-2,
                                             "intersection_tyr_intI")
        one_integrase_one_attC.add_attC(10,
                                        100,
                                        1,
                                        1e-2,
                                        "intersection_tyr_intI")
        self.assertEqual(one_integrase_one_attC.type(), "complete")


    def test_add_promoter(self):
        replicon_name = 'saen.040.p01.10'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        topologies = Topology(2, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        # to test promoter we need to ad attC and integrase first
        # as add_promoter use attc and integrase
        attC = pd.DataFrame({'pos_beg': [104651, 105162, 106018, 107567, 108423, 108743],
                             'pos_end': [104710, 105221, 106087, 107626, 108482, 108832],
                             'strand': [-1] * 6,
                             'evalue': [3.400000e-06, 7.500000e-09, 6.800000e-06, 2.800000e-07, 6.600000e-06, 1.800000e-04],
                             'type_elt': ['attC'] * 6,
                             'annotation': ['attC'] * 6,
                             'model': ['attc_4'] * 6,
                             'distance_2attC': [np.nan, 452.0, 797.0, 1480.0, 797.0, 261.0]
                             },
                            index=['attc_00{}'.format(i) for i in range(1, 7)],
                            columns=self.columns)
        attC = attC.astype(dtype=self.dtype)

        integrase = pd.DataFrame({'pos_beg': 109469,
                                  'pos_end': 110482,
                                  'strand': 1,
                                  'evalue': 1.600000e-24,
                                  'type_elt': 'protein',
                                  'annotation': 'intI',
                                  'model': 'intersection_tyr_intI',
                                  'distance_2attC': np.nan
                                  },
                                 index=['SAEN.040.P01_10_135'],
                                 columns=self.columns)
        integrase = integrase.astype(dtype=self.dtype)

        ##########################################
        # test promoter with attC with integrase #
        ##########################################
        integron = Integron(replicon, self.cfg)
        integron.attC = attC
        integron.integrase = integrase

        integron.add_promoter()

        exp_promoters = pd.DataFrame({'pos_beg': [109413, 109439],
                                      'pos_end': [109447, 109465],
                                      'strand': [1, -1],
                                      'evalue': [np.nan] * 2,
                                      'type_elt': ['Promoter'] * 2,
                                      'annotation': ['Pint_1', 'Pc_1'],
                                      'model': ['NA'] * 2,
                                      'distance_2attC': [np.nan] * 2
                                      },
                                     index=['P_intI1', 'Pc_int1'],
                                     columns=self.columns
                                     )
        exp_promoters = exp_promoters.astype(dtype=self.dtype)

        pdt.assert_frame_equal(exp_promoters, integron.promoter)

        #############################################
        # test promoter with attC without integrase #
        #############################################
        integron = Integron(replicon, self.cfg)
        integron.attC = attC
        integron.add_promoter()

        empty_promoter = pd.DataFrame(columns=self.columns)
        empty_promoter = empty_promoter.astype(dtype=self.dtype)

        pdt.assert_frame_equal(empty_promoter, integron.promoter)

        #############################################
        # test promoter without attC with integrase #
        #############################################
        integron = Integron(replicon, self.cfg)
        integron.integrase = integrase

        integron.add_promoter()

        pdt.assert_frame_equal(exp_promoters, integron.promoter)


    def test_attI(self):
        replicon_name = 'saen.040.p01.10'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        topologies = Topology(2, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        attC = pd.DataFrame({'pos_beg': [104651, 105162, 106018, 107567, 108423, 108743],
                             'pos_end': [104710, 105221, 106087, 107626, 108482, 108832],
                             'strand': [-1] * 6,
                             'evalue': [3.400000e-06, 7.500000e-09, 6.800000e-06, 2.800000e-07, 6.600000e-06, 1.800000e-04],
                             'type_elt': ['attC'] * 6,
                             'annotation': ['attC'] * 6,
                             'model': ['attc_4'] * 6,
                             'distance_2attC': [np.nan, 452.0, 797.0, 1480.0, 797.0, 261.0]
                             },
                            index=['attc_00{}'.format(i) for i in range(1, 7)],
                            columns=self.columns)
        attC = attC.astype(dtype=self.dtype)


        integrase = pd.DataFrame({'pos_beg': 109469,
                                  'pos_end': 110482,
                                  'strand': 1,
                                  'evalue': 1.600000e-24,
                                  'type_elt': 'protein',
                                  'annotation': 'intI',
                                  'model': 'intersection_tyr_intI',
                                  'distance_2attC': np.nan
                                  },
                                 index=['SAEN.040.P01_10_135'],
                                 columns=self.columns)
        integrase = integrase.astype(dtype=self.dtype)

        ##########################################
        # test promoter with attC with integrase #
        ##########################################
        integron = Integron(replicon, self.cfg)
        integron.attC = attC
        integron.integrase = integrase

        exp_attI = pd.DataFrame({'pos_beg': [109330],
                                 'pos_end': [109388],
                                 'strand': [-1],
                                 'evalue': [np.nan],
                                 'type_elt': 'attI',
                                 'annotation': 'attI_1',
                                 'model': 'NA',
                                 'distance_2attC': [np.nan]
                                 },
                                index=['attI1'],
                                columns=self.columns)
        exp_attI = exp_attI.astype(dtype=self.dtype)

        integron.add_attI()

        pdt.assert_frame_equal(exp_attI, integron.attI)

        #############################################
        # test promoter with attC without integrase #
        #############################################
        integron = Integron(replicon, self.cfg)
        integron.attC = attC

        empty_attI = pd.DataFrame(columns=self.columns)
        empty_attI = empty_attI.astype(dtype=self.dtype)

        integron.add_attI()

        pdt.assert_frame_equal(empty_attI, integron.attI)

        #############################################
        # test promoter without attC with integrase #
        #############################################
        integron = Integron(replicon, self.cfg)
        integron.integrase = integrase

        integron.add_attI()

        pdt.assert_frame_equal(exp_attI, integron.attI)


    def test_add_proteins(self):
        replicon_name = 'pssu.001.c01.13'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        topologies = Topology(2, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        prot_file = os.path.join(self._data_dir,
                                 '{}.prt.short'.format(replicon_name))

        args = argparse.Namespace()
        args.gembase = False
        args.prot_file = False
        args.cmsearch = __file__
        args.hmmsearch = __file__
        args.prodigal = __file__
        args.annot_parser_name = None
        cfg = Config(args)
        integron = Integron(replicon, cfg)

        data_attc = {"pos_beg": [3072863, 3073496, 3074121, 3075059, 3075593, 3076281, 3076659],
                     "pos_end": [3072931, 3073555, 3074232, 3075118, 3075652, 3076340, 3076718],
                     "strand": [-1] * 7,
                     "evalue": [2.5e-06, 7e-08, 6.5e-08, 3.2e-06, 4.1e-07, 1.4e-08, 4e-08],
                     "type_elt": ['attC'] * 7,
                     "annotation": ['attC'] * 7,
                     "model": ['attc_4'] * 7,
                     "distance_2attC": [np.nan, 565.0, 566.0, 827.0, 475.0, 629.0, 319.0]}

        attC = pd.DataFrame(data_attc,
                            columns=self.columns,
                            index=['attc_00{}'.format(i) for i in range(len(data_attc['pos_beg']))])
        attC = attC.astype(dtype=self.dtype)

        integron.attC = attC
        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)
        integron.add_proteins(prot_db)

        exp_proteins = pd.DataFrame({'pos_beg': [3071974, 3072950, 3074243, 3076720],
                                     'pos_end': [3072855, 3073468, 3075055, 3077511],
                                     'strand': [-1] * 4,
                                     'evalue': [np.nan] * 4,
                                     'type_elt': ['protein'] * 4,
                                     'annotation': ['protein'] * 4,
                                     'model': ['NA'] * 4,
                                     'distance_2attC': [np.nan] *4
                                     },
                                    index=['PSSU.001.C01_13_281{}'.format(i) for i in range(5, 9)],
                                    columns=self.columns
                                    )
        exp_proteins = exp_proteins.astype(dtype=self.dtype)
        pdt.assert_frame_equal(exp_proteins.sort_index(), integron.proteins.sort_index())


    def test_describe(self):
        replicon_name = "acba.007.p01.13"
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        topologies = Topology(2, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        args = argparse.Namespace()
        args.gembase = False
        args.prot_file = False
        args.cmsearch = __file__
        args.hmmsearch = __file__
        args.prodigal = __file__
        args.eagle_eyes = False
        args.local_max = False
        cfg = Config(args)

        integron = Integron(replicon, cfg)

        data_integrase = {"pos_beg": 55,
                          "pos_end": 1014,
                          "strand": 1,
                          "evalue": 1.900000e-25,
                          "type_elt": "protein",
                          "annotation": "intI",
                          "model": "intersection_tyr_intI",
                          "distance_2attC": np.nan}

        id_int = "ACBA.007.P01_13_1"
        integrase = pd.DataFrame(data_integrase, columns=self.columns, index=[id_int])
        integrase = integrase.astype(dtype=self.dtype)

        data_attc = {"pos_beg": 10,
                     "pos_end": 100,
                     "strand": -1,
                     "evalue": 1.1e-07,
                     "type_elt": "attC",
                     "annotation": "attC",
                     "model": "attc_4",
                     "distance_2attC": np.nan}

        attC = pd.DataFrame(data_attc, columns=self.columns, index=['attc_001'])
        attC = attC.astype(dtype=self.dtype)
        promoter = pd.DataFrame(data_attc, columns=self.columns, index=['prom_001'])
        promoter = promoter.astype(dtype=self.dtype)
        attI = pd.DataFrame(data_attc, columns=self.columns, index=['attI_001'])
        attI = attI.astype(dtype=self.dtype)
        proteins = pd.DataFrame(data_attc, columns=self.columns, index=['prot_001'])
        proteins = proteins.astype(dtype=self.dtype)

        excp_description = pd.concat([integrase, attC, promoter, attI, proteins], ignore_index=False)
        excp_description = excp_description.reset_index()
        excp_description.columns = ["element"] + list(excp_description.columns[1:])
        excp_description["type"] = "complete"
        excp_description["ID_replicon"] = replicon.id
        excp_description["ID_integron"] = id(integron)  # uniq identifier of a given Integron
        excp_description["default"] = "Yes"
        excp_description["considered_topology"] = replicon.topology
        excp_description.drop_duplicates(subset=["element"], inplace=True)

        self.cfg._args.eagle_eyes = False
        self.cfg._args.eagle_eyes = False
        integron.integrase = integrase
        integron.attC = attC
        integron.promoter = promoter
        integron.attI = attI
        integron.proteins = proteins

        recieved_description = integron.describe()
        pdt.assert_frame_equal(recieved_description, excp_description)


    # def test_draw_integron(self):
    #     pass


    def test_has_integrase(self):
        replicon = SeqRecord(Seq.Seq(''), id='foo', name='bar')
        integron = Integron(replicon, self.cfg)
        self.assertFalse(integron.has_integrase())

        replicon = SeqRecord(Seq.Seq(''), id='just_one_integrase', name='bar')
        just_one_integrase = Integron(replicon, self.cfg)
        just_one_integrase.add_integrase(10,
                                         100,
                                         'foo',
                                         1,
                                         1e-2,
                                         "intersection_tyr_intI")
        self.assertTrue(just_one_integrase.has_integrase())


    def test_has_attC(self):
        replicon = SeqRecord(Seq.Seq(''), id='foo')
        integron = Integron(replicon, self.cfg)
        self.assertFalse(integron.has_attC())

        replicon = SeqRecord(Seq.Seq(''), id='just_one_attC', name='bar')
        just_one_attC = Integron(replicon, self.cfg)
        just_one_attC.add_attC(10,
                               100,
                               1,
                               1e-2,
                               "intersection_tyr_intI")
        self.assertTrue(just_one_attC.has_attC())
