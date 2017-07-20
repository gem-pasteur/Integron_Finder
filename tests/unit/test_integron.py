__author__ = 'bneron'

import os
import unittest
import argparse

import pandas as pd
import pandas.util.testing as pdt
import numpy as np
from Bio import SeqIO, Seq

from integron_finder import Integron

### TODO ###
# to remove when global variables will be replaced
# still N_CPU, SIZE_REPLICON, DISTANCE_THRESHOLD, ...

import integron_finder
##########

class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.join(os.path.dirname(__file__), '..', "data"))

    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..' '..')))

        self.columns = ['pos_beg', 'pos_end', 'strand', 'evalue', 'type_elt', 'model', 'distance_2attC', 'annotation']
        self.dtype = {"pos_beg": 'int',
                      "pos_end": 'int',
                      "strand": 'int',
                      "evalue": 'float',
                      "type_elt": 'str',
                      "annotation": 'str',
                      "model": 'str',
                      "distance_2attC": 'float'}


    def test_add_integrase(self):
        id_replicon = "acba.007.p01.13"
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

        integron = Integron(id_replicon)
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
        id_replicon = "acba.007.p01.13"

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

        integron = Integron(id_replicon)
        integron.add_attC(attc_1.loc['attc_001', 'pos_beg'],
                          attc_1.loc['attc_001', 'pos_end'],
                          attc_1.loc['attc_001', 'strand'],
                          attc_1.loc['attc_001', 'evalue'],
                          attc_1.loc['attc_001', 'model'])

        pdt.assert_frame_equal(attc_1, integron.attC)

        integron_finder.SIZE_REPLICON = 3

        attc_2 = pd.DataFrame(data_attc_1,
                              columns=self.columns,
                              index=['attc_002'])
        attc_2 = attc_2.astype(dtype=self.dtype)
        attc_2['pos_beg'] = attc_2['pos_beg'] + 100
        attc_2['pos_end'] = attc_2['pos_end'] + 100
        attc_2["distance_2attC"] = (attc_2.loc['attc_002', 'pos_beg'] - attc_1.loc['attc_001', 'pos_end']) % \
                                   integron_finder.SIZE_REPLICON
        attc = attc_1.append(attc_2)

        integron.add_attC(attc_2.loc['attc_002', 'pos_beg'],
                          attc_2.loc['attc_002', 'pos_end'],
                          attc_2.loc['attc_002', 'strand'],
                          attc_2.loc['attc_002', 'evalue'],
                          attc_2.loc['attc_002', 'model'])
        pdt.assert_frame_equal(attc, integron.attC)


    def test_type(self):
        no_integrase = Integron("foo")
        self.assertIsNone(no_integrase.type())

        just_one_integrase = Integron("just_one_integrase")
        just_one_integrase.add_integrase(10,
                                         100,
                                         'foo',
                                         1,
                                         1e-2,
                                         "intersection_tyr_intI")
        self.assertEqual(just_one_integrase.type(), "In0")

        just_one_attC = Integron("just_one_attC")
        just_one_attC.add_attC(10,
                               100,
                               1,
                               1e-2,
                               "intersection_tyr_intI")
        self.assertEqual(just_one_attC.type(), "CALIN")

        one_integrase_one_attC = Integron("one_integrase_one_attC")
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
        replicon_path = os.path.join(self._data_dir, "Replicons", replicon_name + '.fst')
        integron_finder.SIZE_REPLICON = 148711
        integron_finder.PROT_file = os.path.join(self._data_dir,
                                                 'Proteins',
                                                 '{}.prt'.format(replicon_name))
        integron_finder.MODEL_DIR = os.path.join(self.integron_home, "data", "Models")
        integron_finder.SEQUENCE = SeqIO.read(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)

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
        integron = Integron(replicon_name)
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
        integron = Integron(replicon_name)
        integron.attC = attC
        integron.add_promoter()

        empty_promoter = pd.DataFrame(columns=self.columns)
        empty_promoter = empty_promoter.astype(dtype=self.dtype)

        pdt.assert_frame_equal(empty_promoter, integron.promoter)

        #############################################
        # test promoter without attC with integrase #
        #############################################
        integron = Integron(replicon_name)
        integron.integrase = integrase

        integron.add_promoter()

        pdt.assert_frame_equal(exp_promoters, integron.promoter)


    def test_attI(self):
        replicon_name = 'saen.040.p01.10'
        replicon_path = os.path.join(self._data_dir, "Replicons", replicon_name + '.fst')
        integron_finder.SIZE_REPLICON = 148711
        integron_finder.PROT_file = os.path.join(self._data_dir,
                                                 'Proteins',
                                                 '{}.prt'.format(replicon_name))
        integron_finder.MODEL_DIR = os.path.join(self.integron_home, "data", "Models")
        integron_finder.SEQUENCE = SeqIO.read(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)

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
        integron = Integron(replicon_name)
        integron.attC = attC
        integron.integrase = integrase

        empty_proteins = pd.DataFrame(columns=self.columns)
        empty_proteins = empty_proteins.astype(dtype=self.dtype)
        integron.proteins = empty_proteins

        empty_promoters = pd.DataFrame(columns=self.columns)
        empty_promoters = empty_promoters.astype(dtype=self.dtype)
        integron.promoter = empty_promoters

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
        integron = Integron(replicon_name)
        integron.attC = attC

        empty_attI = pd.DataFrame(columns=self.columns)
        empty_attI = empty_attI.astype(dtype=self.dtype)

        integron.add_attI()

        pdt.assert_frame_equal(empty_attI, integron.attI)

        #############################################
        # test promoter without attC with integrase #
        #############################################
        integron = Integron(replicon_name)
        integron.integrase = integrase

        integron.add_attI()

        pdt.assert_frame_equal(exp_attI, integron.attI)


    def test_add_proteins(self):
        replicon_name = 'pssu.001.c01.13'
        integron_finder.SIZE_REPLICON = 3419049
        integron_finder.PROT_file = os.path.join(self._data_dir,
                                                 '{}.prt.short'.format(replicon_name))
        args = argparse.Namespace()
        args.gembase = False
        integron_finder.args = args

        integron = Integron(replicon_name)

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

        integron.add_proteins()

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
        pdt.assert_frame_equal(exp_proteins, integron.proteins)


    def test_describe(self):
        id_replicon = "acba.007.p01.13"
        integron = Integron(id_replicon)

        data_integrase = {"pos_beg": 55,
                          "pos_end": 1014,
                          "strand": 1,
                          "evalue": 1.900000e-25,
                          "type_elt": "protein",
                          "annotation": "intI",
                          "model": "intersection_tyr_intI",
                          "distance_2attC": np.nan}

        id_int = "ACBA.007.P01_13_1"
        integrase = pd.DataFrame(data_integrase,
                                 columns=self.columns,
                                 index=[id_int])
        integrase = integrase.astype(dtype=self.dtype)

        data_attc = {"pos_beg": 10,
                     "pos_end": 100,
                     "strand": -1,
                     "evalue": 1.1e-07,
                     "type_elt": "attC",
                     "annotation": "attC",
                     "model": "attc_4",
                     "distance_2attC": np.nan}

        attC = pd.DataFrame(data_attc,
                            columns=self.columns,
                            index=['attc_001'])
        attC = attC.astype(dtype=self.dtype)
        promoter = pd.DataFrame(data_attc,
                                columns=self.columns,
                                index=['prom_001'])
        promoter = promoter.astype(dtype=self.dtype)
        attI = pd.DataFrame(data_attc,
                            columns=self.columns,
                            index=['attI_001'])
        attI = attI.astype(dtype=self.dtype)
        proteins = pd.DataFrame(data_attc,
                                columns=self.columns,
                                index=['prot_001'])
        proteins = proteins.astype(dtype=self.dtype)

        excp_description = pd.concat([integrase, attC, promoter, attI, proteins], ignore_index=False)
        excp_description = excp_description.reset_index()
        excp_description.columns = ["element"] + list(excp_description.columns[1:])
        excp_description["type"] = "complete"
        excp_description["ID_replicon"] = id_replicon
        excp_description["ID_integron"] = id(integron)  # uniq identifier of a given Integron
        excp_description["default"] = "Yes"
        excp_description.drop_duplicates(subset=["element"], inplace=True)

        args = argparse.Namespace
        args.eagle_eyes = False
        args.local_max = False
        integron_finder.args = args

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
        integron = Integron("foo")
        self.assertFalse(integron.has_integrase())

        just_one_integrase = Integron("just_one_integrase")
        just_one_integrase.add_integrase(10,
                                         100,
                                         'foo',
                                         1,
                                         1e-2,
                                         "intersection_tyr_intI")
        self.assertTrue(just_one_integrase.has_integrase())


    def test_has_attC(self):
        integron = Integron("foo")
        self.assertFalse(integron.has_attC())

        just_one_attC = Integron("just_one_attC")
        just_one_attC.add_attC(10,
                               100,
                               1,
                               1e-2,
                               "intersection_tyr_intI")
        self.assertTrue(just_one_attC.has_attC())
