__author__ = 'bneron'

import os
import pandas as pd
import numpy as np
from integron_finder import Integron

class Test(pd.util.testing.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "data")


    def test_add_integrase(self):
        #                       pos_beg, pos_end, strand,   evalue,   type_elt, annotation, model,           distance_2attC
        # ACBA.007.P01_13_1,    55     ,    1014,    1,   1.900000e-25, protein, intI,    intersection_tyr_intI, NaN
        # LIAN.001.C02_10_825,  934165 ,  934689,   -1,   1.900000e-26, protein, intI,    intersection_tyr_intI, NaN
        # PSSU.001.C01_13_1919, 2131128, 2132135,    1,   5.100000e-27, protein, intI,    intersection_tyr_intI, NaN

        id_replicon = "acba.007.p01.13"
        data_integrase = {"pos_beg" : 55,
                          "pos_end" : 1014,
                          "strand" : 1,
                          "evalue" : 1.900000e-25,
                          "type_elt" : "protein",
                          "annotation" : "intI",
                          "model" : "intersection_tyr_intI",
                          "distance_2attC" : np.nan}
        id_int = "ACBA.007.P01_13_1"
        df = pd.DataFrame(columns=["pos_beg", "pos_end", "strand", "evalue", "type_elt",
                                   "model", "distance_2attC", "annotation"])
        tmp_df = pd.DataFrame()
        tmp_df["pos_beg"] = [data_integrase["pos_beg"]]
        tmp_df["pos_end"] = [data_integrase["pos_end"]]
        tmp_df["strand"] = [data_integrase["strand"]]
        tmp_df["evalue"] = [data_integrase["evalue"]]
        tmp_df["type_elt"] = data_integrase["type_elt"]
        tmp_df["annotation"] = data_integrase["annotation"]
        tmp_df["model"] = [data_integrase["model"]]
        tmp_df.index = [id_int]
        tmp_df["distance_2attC"] = [np.nan]
        df = df.append(tmp_df)
        integron = Integron(id_replicon)
        integron.add_integrase(data_integrase["pos_beg"],
                               data_integrase["pos_end"],
                               [id_int],
                               data_integrase["strand"],
                               data_integrase["evalue"],
                               data_integrase["model"],
                               )
        self.assert_frame_equal(df, integron.integrase, check_dtype=False)

    # def test_add_attc(self):
    #     pass
    #
    # def test_type(self):
    #     pass
    #
    # def test_add_promoter(self):
    #     pass
    #
    # def test_attI(self):
    #     pass
    #
    # def add_proteins(self):
    #     pass
    #
    # def test_describe(self):
    #     pass
    #
    # def test_draw_integron(self):
    #     pass

