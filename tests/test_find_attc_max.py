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
import tempfile
import shutil
import argparse

import numpy as np
import pandas as pd
import pandas.testing as pdt


try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.integron import Integron
from integron_finder.config import Config
from integron_finder.utils import FastaIterator
from integron_finder.topology import Topology
from integron_finder.attc import find_attc_max


class TestFindAttCMax(IntegronTest):

    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
        self._tmp_dir = tempfile.TemporaryDirectory(prefix='tmp_test_integron_finder')
        self.tmp_dir = self._tmp_dir.name
        if os.path.exists(self.tmp_dir) and os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)

        args = argparse.Namespace()
        args.gembase = False
        args.prot_file = False
        args.cmsearch = __file__
        args.hmmsearch = __file__
        args.prodigal = __file__
        args.attc_model = 'attc_4.cm'
        args.max_attc_size = 200
        args.min_attc_size = 40
        args.distance_threshold = 4000  # (4kb at least between 2 different arrays)
        args.eagle_eyes = False
        args.local_max = False
        self.cfg = Config(args)
        self.cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

        replicon_name = 'OBAL001.B.00005.C001'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))

        topologies = Topology(1, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            self.replicon = next(sequences_db)

        self.integron = Integron(self.replicon, self.cfg)

        self.columns = ['pos_beg', 'pos_end', 'strand', 'evalue', 'type_elt', 'model', 'distance_2attC', 'annotation']
        self.dtype = {"pos_beg": 'int',
                      "pos_end": 'int',
                      "strand": 'int',
                      "evalue": 'float',
                      "type_elt": 'str',
                      "annotation": 'str',
                      "model": 'str',
                      "distance_2attC": 'float'}

        self.max_dtype = {'Accession_number': 'str',
                          'cm_attC': 'str',
                          'cm_debut': 'int',
                          'cm_fin': 'int',
                          'pos_beg': 'int',
                          'pos_end': 'int',
                          'sens': 'str',
                          'evalue': 'float'}
        self.max_cols = ['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg', 'pos_end', 'sens', 'evalue']

    def tearDown(self):
        self._tmp_dir.cleanup()


    def test_find_attc_max_linear(self):

        integrase = pd.DataFrame({'pos_beg': 1545830,
                                  'pos_end': 1546807,
                                  'strand': -1,
                                  'evalue': 1.100000e-21,
                                  'type_elt': 'protein',
                                  'annotation': 'intI',
                                  'model': 'intersection_tyr_intI',
                                  'distance_2attC': np.nan
                                  },
                                 index=['OBAL001.B.00005.C001_141'],
                                 columns=self.columns)
        integrase = integrase.astype(dtype=self.dtype)
        self.integron.integrase = integrase

        attC = pd.DataFrame({'pos_beg': [1547800, 1548775],
                             'pos_end': [1547859, 1548834],
                             'strand': [1, 1],
                             'evalue': [0.00049, 0.00009],
                             'type_elt': ['attC', 'attC'],
                             'annotation': ['attC', 'attC'],
                             'model': ['attc_4', 'attc_4'],
                             'distance_2attC': [np.nan, 916.0]
                             },
                            index=['attc_001', 'attc_002'],
                            columns=self.columns)
        attC = attC.astype(dtype=self.dtype)
        self.integron.attC = attC
        integrons = [self.integron]

        max_final = find_attc_max(integrons, self.replicon,
                                  self.cfg.distance_threshold, self.cfg.model_attc_path,
                                  self.cfg.max_attc_size, self.cfg.min_attc_size,
                                  circular=False,
                                  out_dir=self.tmp_dir)

        exp = pd.DataFrame({'Accession_number': ['OBAL001.B.00005.C001', 'OBAL001.B.00005.C001'],
                            'cm_attC': ['attc_4', 'attc_4'],
                            'cm_debut': [4, 1],
                            'cm_fin': [44, 47],
                            'pos_beg': [1547800, 1548775],
                            'pos_end': [1547859, 1548834],
                            'sens': ['+', '+'],
                            'evalue': [0.000240,  0.000045]
                            },
                           index=[0, 1],
                           columns=self.max_cols
                           )
        exp = exp.astype(dtype=self.max_dtype)
        pdt.assert_frame_equal(max_final, exp)


    def test_find_attc_max_complete(self):
        integrase = pd.DataFrame({'pos_beg': 1545830,
                                  'pos_end': 1546807,
                                  'strand': -1,
                                  'evalue': 1.100000e-21,
                                  'type_elt': 'protein',
                                  'annotation': 'intI',
                                  'model': 'intersection_tyr_intI',
                                  'distance_2attC': np.nan
                                  },
                                 index=['OBAL001.B.00005.C001_141'],
                                 columns=self.columns)
        integrase = integrase.astype(dtype=self.dtype)
        self.integron.integrase = integrase

        attC = pd.DataFrame({'pos_beg': [1547800, 1548775],
                             'pos_end': [1547859, 1548834],
                             'strand': [1, 1],
                             'evalue': [0.00049, 0.00009],
                             'type_elt': ['attC', 'attC'],
                             'annotation': ['attC', 'attC'],
                             'model': ['attc_4', 'attc_4'],
                             'distance_2attC': [np.nan, 916.0]
                             },
                            index=['attc_001', 'attc_002'],
                            columns=self.columns)
        attC = attC.astype(dtype=self.dtype)
        self.integron.attC = attC
        integrons = [self.integron]

        max_final = find_attc_max(integrons, self.replicon,
                                  self.cfg.distance_threshold, self.cfg.model_attc_path,
                                  self.cfg.max_attc_size, self.cfg.min_attc_size,
                                  circular=True,
                                  out_dir=self.tmp_dir)

        exp = pd.DataFrame({'Accession_number': ['OBAL001.B.00005.C001', 'OBAL001.B.00005.C001'],
                            'cm_attC': ['attc_4', 'attc_4'],
                            'cm_debut': [4, 1],
                            'cm_fin': [44, 47],
                            'pos_beg': [1547800, 1548775],
                            'pos_end': [1547859, 1548834],
                            'sens': ['+', '+'],
                            'evalue': [0.000240,  0.000045]
                            },
                           index=[0, 1],
                           columns=self.max_cols
                           )
        exp = exp.astype(dtype=self.max_dtype)
        pdt.assert_frame_equal(max_final, exp)


    def test_find_attc_max_calin(self):
        attC = pd.DataFrame({'pos_beg': [421689],
                             'pos_end': [421764],
                             'strand': [1],
                             'evalue': [0.13],
                             'type_elt': ['attC'],
                             'annotation': ['attC'],
                             'model': ['attc_4'],
                             'distance_2attC': [np.nan]
                             },
                            index=['attc_001'],
                            columns=self.columns)
        attC = attC.astype(dtype=self.dtype)
        self.integron.attC = attC
        integrons = [self.integron]

        for is_circ in True, False:
            with self.subTest(is_circ=is_circ):
                max_final = find_attc_max(integrons, self.replicon,
                                          self.cfg.distance_threshold, self.cfg.model_attc_path,
                                          self.cfg.max_attc_size, self.cfg.min_attc_size,
                                          circular=is_circ,
                                          out_dir=self.tmp_dir)

                exp = pd.DataFrame({'Accession_number': ['OBAL001.B.00005.C001'],
                                    'cm_attC': ['attc_4'],
                                    'cm_debut': [1],
                                    'cm_fin': [47],
                                    'pos_beg': [421689],
                                    'pos_end': [421764],
                                    'sens': ['+'],
                                    'evalue': [0.062]
                                    },
                                   index=[0],
                                   columns=self.max_cols
                                   )
                exp = exp.astype(dtype=self.max_dtype)

                cols_2_compare = ['Accession_number', 'cm_attC', 'pos_beg', 'pos_end', 'sens', 'evalue']
                # cm_debut and cm_fin can vary depending if the data are generated with
                # local_max or get from previous regular search (in this last case cm_debut & cm_fin are dummy)
                pdt.assert_frame_equal(max_final[cols_2_compare], exp[cols_2_compare])


    def test_find_attc_max_In0(self):
        replicon_name = 'ESCO001.B.00018.P002'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))

        topologies = Topology(1, 'circ')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        integron = Integron(replicon, self.cfg)

        integrase = pd.DataFrame({'pos_beg': [90229],
                                  'pos_end': [91242],
                                  'strand': -1,
                                  'evalue': 1.400000e-24,
                                  'type_elt': 'protein',
                                  'annotation': 'intI',
                                  'model': 'intersection_tyr_intI',
                                  'distance_2attC': np.nan
                                  },
                                 index=['ESCO001.B.00018.P002_106'],
                                 columns=self.columns)
        integrase = integrase.astype(dtype=self.dtype)
        integron.integrase = integrase
        integrons = [integron]
        for is_circ in True, False:
            with self.subTest(is_circ=is_circ):
                max_final = find_attc_max(integrons, replicon,
                                          self.cfg.distance_threshold, self.cfg.model_attc_path,
                                          self.cfg.max_attc_size, self.cfg.min_attc_size,
                                          circular=is_circ,
                                          out_dir=self.tmp_dir)

                exp = pd.DataFrame(columns=self.max_cols)
                exp = exp.astype(dtype=self.max_dtype)
                pdt.assert_frame_equal(max_final, exp)
