__author__ = 'bneron'

import os
import unittest
import argparse

import pandas as pd
import pandas.util.testing as pdt
import numpy as np

from integron_finder import Integron

### TODO ###
# to remove when global variables will be replaced
# still N_CPU, SIZE_REPLICON, DISTANCE_THRESHOLD, ...

import integron_finder
##########

class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "data")


    def test_add_integrase(self):
        #                       pos_beg, pos_end, strand,   evalue,   type_elt, annotation, model,           distance_2attC
        # ACBA.007.P01_13_1,    55     ,    1014,    1,   1.900000e-25, protein, intI,    intersection_tyr_intI, NaN
        # LIAN.001.C02_10_825,  934165 ,  934689,   -1,   1.900000e-26, protein, intI,    intersection_tyr_intI, NaN
        # PSSU.001.C01_13_1919, 2131128, 2132135,    1,   5.100000e-27, protein, intI,    intersection_tyr_intI, NaN

        id_replicon = "acba.007.p01.13"
        data_integrase = {"pos_beg": 55,
                          "pos_end": 1014,
                          "strand": 1,
                          "evalue": 1.900000e-25,
                          "type_elt": "protein",
                          "annotation": "intI",
                          "model": "intersection_tyr_intI",
                          "distance_2attC": np.nan}
        integrase_dtype = {"pos_beg": 'int',
                           "pos_end": 'int',
                           "strand": 'int',
                           "evalue": 'float',
                           "type_elt": 'str',
                           "annotation": 'str',
                           "model": 'str',
                           "distance_2attC": 'float'}
        id_int = "ACBA.007.P01_13_1"

        df = pd.DataFrame(data_integrase, index=[id_int])
        df = df.astype(dtype=integrase_dtype)

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
        attc_dtype = {"pos_beg": 'int',
                      "pos_end": 'int',
                      "strand": 'int',
                      "evalue": 'float',
                      "type_elt": 'str',
                      "annotation": 'str',
                      "model": 'str',
                      "distance_2attC": 'float'}

        attc_1 = pd.DataFrame(data_attc_1, index=['attc_001'])
        attc_1 = attc_1.astype(dtype=attc_dtype)
        integron = Integron(id_replicon)
        integron.add_attC(attc_1.loc['attc_001', 'pos_beg'],
                          attc_1.loc['attc_001', 'pos_end'],
                          attc_1.loc['attc_001', 'strand'],
                          attc_1.loc['attc_001', 'evalue'],
                          attc_1.loc['attc_001', 'model'])

        pdt.assert_frame_equal(attc_1, integron.attC)

        integron_finder.SIZE_REPLICON = 3

        attc_2 = pd.DataFrame(data_attc_1, index=['attc_002'])
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


    # def test_add_promoter(self):
    #     pass
    #
    # def test_attI(self):
    #     pass
    #
    # def add_proteins(self):
    #     pass
    #
    def test_describe(self):
        id_replicon = "acba.007.p01.13"
        integron = Integron(id_replicon)

        dtype = {"pos_beg": 'int',
                 "pos_end": 'int',
                 "strand": 'int',
                 "evalue": 'float',
                 "type_elt": 'str',
                 "annotation": 'str',
                 "model": 'str',
                 "distance_2attC": 'float'}

        data_integrase = {"pos_beg": 55,
                          "pos_end": 1014,
                          "strand": 1,
                          "evalue": 1.900000e-25,
                          "type_elt": "protein",
                          "annotation": "intI",
                          "model": "intersection_tyr_intI",
                          "distance_2attC": np.nan}

        id_int = "ACBA.007.P01_13_1"
        integrase = pd.DataFrame(data_integrase, index=[id_int])
        integrase = integrase.astype(dtype=dtype)

        data_attc = {"pos_beg": 10,
                     "pos_end": 100,
                     "strand": -1,
                     "evalue": 1.1e-07,
                     "type_elt": "attC",
                     "annotation": "attC",
                     "model": "attc_4",
                     "distance_2attC": np.nan}

        attC = pd.DataFrame(data_attc, index=['attc_001'])
        attC = attC.astype(dtype=dtype)
        promoter = pd.DataFrame(data_attc, index=['prom_001'])
        promoter = promoter.astype(dtype=dtype)
        attI = pd.DataFrame(data_attc, index=['attI_001'])
        attI = attI.astype(dtype=dtype)
        proteins = pd.DataFrame(data_attc, index=['prot_001'])
        proteins = proteins.astype(dtype=dtype)

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
