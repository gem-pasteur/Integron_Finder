# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright (c) 2015 - 2018  Institut Pasteur, Paris.                              #
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
import tempfile
import shutil
import argparse

import numpy as np
import pandas as pd
import pandas.util.testing as pdt

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.integron import Integron
from integron_finder.topology import Topology
from integron_finder.utils import FastaIterator
from integron_finder.config import Config
from integron_finder import results


class TestUtils(IntegronTest):

    def setUp(self):
        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        if os.path.exists(self.tmp_dir) and os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)
        
    def tearDown(self):
        if os.path.exists(self.tmp_dir):
            try:
                shutil.rmtree(self.tmp_dir)
            except Exception:
                pass

    def test_integrons_report(self):
        replicon_name = "acba.007.p01.13"
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        topologies = Topology('circ')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        args = argparse.Namespace()
        cfg = Config(args)
        cfg._args.eagle_eyes = False
        cfg._args.eagle_eyes = False
        cfg._args.local_max = False

        integron = Integron(replicon, cfg)
        columns = ['pos_beg', 'pos_end', 'strand', 'evalue', 'type_elt', 'model', 'distance_2attC', 'annotation']
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

        integrase = pd.DataFrame(data_integrase, columns=columns, index=[id_int])
        integrase = integrase.astype(dtype=dtype)

        data_attc = {"pos_beg": [17825, 19080, 19618],
                     "pos_end": [17884, 19149, 19726],
                     "strand": [-1] * 3,
                     "evalue": [1.000000e-09, 1.000000e-04, 1.100000e-07],
                     "type_elt": ["attC"] * 3,
                     "annotation": ["attC"] * 3,
                     "model": ["attc_4"] * 3,
                     "distance_2attC": [np.nan, 1196.0, 469.0]}

        attC = pd.DataFrame(data_attc,
                            columns=columns,
                            index=['attc_00{}'.format(i) for i in range(1, 4)])
        attC = attC.astype(dtype=dtype)

        promoter = pd.DataFrame({'pos_beg': 25,
                                 'pos_end': 51,
                                 'strand': -1,
                                 'evalue': np.nan,
                                 'type_elt': 'Promoter',
                                 'annotation': 'Pc_1',
                                 'model': np.nan,
                                 'distance_2attC': np.nan
                                 },
                                index=['Pc_int1'],
                                columns=columns
                                )
        promoter = promoter.astype(dtype=dtype)

        proteins = pd.DataFrame({'pos_beg': [17375, 17886, 19090, 19721],
                                 'pos_end': [17722, 18665, 19749, 20254],
                                 'strand': [-1] * 4,
                                 'evalue': [np.nan] * 4,
                                 'type_elt': ['protein'] * 4,
                                 'annotation': ['protein'] * 4,
                                 'model': [np.nan] * 4,
                                 'distance_2attC': [np.nan] * 4
                                 },
                                index=['ACBA.007.P01_13_2{}'.format(i) for i in range(0, 4)],
                                columns=columns
                                )
        proteins = proteins.astype(dtype=dtype)

        integron.integrase = integrase
        integron.attC = attC
        integron.promoter = promoter
        integron.proteins = proteins
        report = results.integrons_report([integron])
        exp_report = pd.read_table(
            self.find_data(os.path.join('Results_Integron_Finder_{}'.format(replicon_name),
                                        '{}.integrons'.format(replicon_name)
                                        ))
        )
        exp_report = exp_report.astype(dtype=dtype)
        pdt.assert_frame_equal(exp_report, report)


    def test_merge_integrons(self):
        f1 = os.path.join(self.tmp_dir, 'f1')
        dtype = {"ID_integron": 'str',
                 "ID_replicon": 'str',
                 "element": 'str',
                 "pos_beg": 'int',
                 "pos_end": 'int',
                 "strand": 'int',
                 "evalue": 'float',
                 "type_elt": 'str',
                 "annotation": 'str',
                 "model": 'str',
                 "type": 'str',
                 "default": 'str',
                 "distance_2attC": 'float',
                 "considered_topology": 'str'
                 }
        res1 = pd.DataFrame({"ID_integron": ['integron_01'],
                             "ID_replicon": ['ACBA.007.P01_13'],
                             "element": ['Pc_int1'],
                             "pos_beg": [25],
                             "pos_end": [51],
                             "strand": [-1],
                             "evalue": [None],
                             "type_elt": ['Promoter'],
                             "annotation": ['Pc_1'],
                             "model": [None],
                             "type": ['complete'],
                             "default": ['No'],
                             "distance_2attC": [None],
                             "considered_topology": ['lin'],
                             })
        res1 = res1.astype(dtype=dtype)
        res1.to_csv(f1, sep="\t", index=False, na_rep="NA")

        f2 = os.path.join(self.tmp_dir, 'f2')
        res2 = pd.DataFrame({"ID_integron": ['integron_02'],
                             "ID_replicon": ['LIAN.001.C02_10'],
                             "element": ['LIAN.001.C02_10_825'],
                             "pos_beg": [934165],
                             "pos_end": [934689],
                             "strand": [-1],
                             "evalue": [1.9e-26],
                             "type_elt": ['protein intI'],
                             "annotation": ['intersection_tyr_intI'],
                             "model": [None],
                             "type": ['complete'],
                             "default": ['No'],
                             "distance_2attC": [None],
                             "considered_topology": ['lin']
                             })
        res2 = res2.astype(dtype=dtype)
        res2.to_csv(f2, sep="\t", index=False, na_rep="NA")

        expected_res = pd.concat([res1, res2])
        res = results.merge_results(f1, f2)
        pdt.assert_frame_equal(expected_res, res)

        f3 = os.path.join(self.tmp_dir, 'f3')
        with open(f3, 'w') as f3_fo:
            f3_fo.write("# No Integron found \n")

        res = results.merge_results(f1, f3, f2)
        pdt.assert_frame_equal(expected_res, res)

        res = results.merge_results()
        expected_res = pd.DataFrame(columns=['ID_integron', 'ID_replicon', 'element',
                                             'pos_beg', 'pos_end', 'strand', 'evalue',
                                             'type_elt annotation', 'model', 'type', 'default',
                                             'distance_2attC', 'considered_topology'])
        pdt.assert_frame_equal(expected_res, res)


    def test_merge_summary(self):
        f1 = os.path.join(self.tmp_dir, 'f1')
        dtype = {'ID_replicon': 'str',
                 'ID_integron': 'str',
                 'complete': 'int',
                 'In0': 'int',
                 'CALIN': 'int'
                 }
        res1 = pd.DataFrame({'ID_replicon': ['ACBA.007.P01_13'],
                             'ID_integron': ['integron_01'],
                             'complete': [1],
                             'In0': [0],
                             'CALIN': [0],
                             }, columns=['ID_replicon', 'ID_integron', 'complete', 'In0', 'CALIN'])
        res1 = res1.astype(dtype=dtype)
        res1.to_csv(f1, sep="\t", na_rep="NA", index=False,
                    columns=['ID_replicon', 'ID_integron', 'complete', 'In0', 'CALIN'])

        f2 = os.path.join(self.tmp_dir, 'f2')
        res2 = pd.DataFrame({'ID_replicon': ['LIAN.001.C02_10'] * 6,
                             'ID_integron': ['integron_0{}'.format(i) for i in range(1,7)],
                             'complete': [0, 1, 0, 0, 0, 0],
                             'In0': [0] * 6,
                             'CALIN': [1, 0, 1, 1, 1, 1],
                             }, columns=['ID_replicon', 'ID_integron', 'complete', 'In0', 'CALIN'])
        res2 = res2.astype(dtype=dtype)
        res2.to_csv(f2, sep="\t", na_rep="NA", index=False,
                    columns=['ID_replicon', 'ID_integron', 'complete', 'In0', 'CALIN'])

        expected_res = pd.concat([res1, res2])
        res = results.merge_results(f1, f2)

        pdt.assert_frame_equal(expected_res, res)


    def test_summary(self):
        acba_res = self.find_data('Results_Integron_Finder_acba.007.p01.13/acba.007.p01.13.integrons')
        acba_df = pd.read_table(acba_res)
        summary = results.summary(acba_df)
        dtype = {'ID_replicon': 'str',
                 'ID_integron': 'str',
                 'complete': 'int',
                 'In0': 'int',
                 'CALIN': 'int'
                 }
        exp = pd.DataFrame({'ID_replicon': ['ACBA.007.P01_13'],
                             'ID_integron': ['integron_01'],
                             'complete': [1],
                             'In0': [0],
                             'CALIN': [0],
                            }, columns=['ID_replicon', 'ID_integron', 'complete', 'In0', 'CALIN'])
        exp = exp.astype(dtype=dtype)

        lian_res = self.find_data('lian.001.c02.10_simple.integrons')
        lian_df = pd.read_table(lian_res)
        summary = results.summary(lian_df)
        exp = pd.DataFrame({'ID_replicon': ['LIAN.001.C02_10'] * 6,
                             'ID_integron': ['integron_0{}'.format(i) for i in range(1, 7)],
                             'complete': [0, 1, 0, 0, 0, 0],
                             'In0': [0] * 6,
                             'CALIN': [1, 0, 1, 1, 1, 1],
                            }, columns=['ID_replicon', 'ID_integron', 'complete', 'In0', 'CALIN'])
        exp = exp.astype(dtype=dtype)
        pdt.assert_frame_equal(exp, summary)


    def test_filter_calin(self):
        lian_res = self.find_data('lian.001.c02.10_simple.integrons')
        lian_df = pd.read_table(lian_res)
        filtered = results.filter_calin(lian_df, 4)
        exp = lian_df[lian_df.ID_integron != 'integron_03']
        pdt.assert_frame_equal(exp, filtered)