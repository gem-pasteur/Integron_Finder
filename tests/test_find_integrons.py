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

# # display warning only for non installed integron_finder
# from Bio import BiopythonExperimentalWarning
# import warnings
# warnings.simplefilter('ignore', BiopythonExperimentalWarning)

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.integron import find_integron
from integron_finder.config import Config
from integron_finder.utils import FastaIterator
from integron_finder.topology import Topology
from integron_finder.infernal import read_infernal
from integron_finder.prot_db import ProdigalDB


class TestFindIntegons(IntegronTest):


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
        # by default use empty Namespace just to match with api
        # if need some values create a local config with these values
        # this rule should avoid side effects
        args = argparse.Namespace()
        args.gembase = False
        args.prot_file = False
        args.cmsearch = __file__
        args.hmmsearch = __file__
        args.prodigal = __file__
        self.cfg = Config(args)
        self.cfg._args.eagle_eyes = False
        self.cfg._args.local_max = False
        self.cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

        self.columns = ['pos_beg', 'pos_end', 'strand', 'evalue', 'type_elt', 'model', 'distance_2attC', 'annotation']
        self.dtype = {"pos_beg": 'int',
                      "pos_end": 'int',
                      "strand": 'int',
                      "evalue": 'float',
                      "type_elt": 'str',
                      "annotation": 'str',
                      "model": 'str',
                      "distance_2attC": 'float'}
        self.set_log_level('INFO')


    def tearDown(self):
        self.set_log_level('WARNING')
        self._tmp_dir.cleanup()


    def test_find_integron(self):
        replicon_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        prot_file = self.find_data(os.path.join('Proteins', prot_name + '.prt'))
        topologies = Topology(1, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)
        replicon_results_path = self.find_data(os.path.join('Results_Integron_Finder_{}'.format(replicon_name),
                                                            'tmp_{}'.format(replicon.id)))
        attc_file = os.path.join(replicon_results_path, '{}_attc_table.res'.format(replicon.id))
        intI_file = os.path.join(replicon_results_path, '{}_intI.res'.format(replicon.id))
        phageI_file = os.path.join(replicon_results_path, '{}_phage_int.res'.format(replicon.id))

        args = argparse.Namespace()
        args.gembase = False
        args.prot_file = False
        args.cmsearch = __file__
        args.hmmsearch = __file__
        args.prodigal = __file__
        args.no_proteins = True
        args.keep_palindromes = True
        args.distance_threshold = 4000
        args.attc_model = 'attc_4.cm'
        args.evalue_attc = 1.0
        args.max_attc_size = 200
        args.min_attc_size = 40
        args.calin_threshold = 2
        args.local_max = False
        cfg = Config(args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')
        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)

        exp_msg = """In replicon {}, there are:
- 0 complete integron(s) found with a total 0 attC site(s)
- 1 CALIN element(s) found with a total of 3 attC site(s)
- 0 In0 element(s) found with a total of 0 attC site""".format(replicon.id)
        with self.catch_log() as log:
            integrons = find_integron(replicon, prot_db, intI_file, phageI_file, cfg, attc_file=attc_file)
            catch_msg = log.get_value().strip()

        self.assertEqual(catch_msg, exp_msg)
        self.assertEqual(len(integrons), 1)
        integron = integrons[0]
        self.assertEqual(integron.replicon.id, replicon.id)

        exp = pd.DataFrame({'annotation': ['attC'] * 3,
                            'distance_2attC': [np.nan, 1196.0, 469.0],
                            'evalue': [1.000000e-09, 1.000000e-04, 1.100000e-07],
                            'model': ['attc_4'] * 3,
                            'pos_beg': [17825, 19080, 19618],
                            'pos_end': [17884, 19149, 19726],
                            'strand': [-1, -1, -1],
                            'type_elt': 'attC'},
                           columns=self.columns,
                           index=['attc_001', 'attc_002', 'attc_003'])
        pdt.assert_frame_equal(integron.attC, exp)

        exp = pd.DataFrame(columns=self.columns,)
        exp = exp.astype(dtype=self.dtype)

        pdt.assert_frame_equal(integron.integrase, exp)
        pdt.assert_frame_equal(integron.promoter, exp)
        pdt.assert_frame_equal(integron.attI, exp)
        pdt.assert_frame_equal(integron.proteins, exp)


    def test_find_integron_calin_threshold(self):
        replicon_name = 'ESCO001.B.00018.P002'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        prot_file = self.find_data(os.path.join('Proteins', replicon_name + '.prt'))
        topologies = Topology(1, 'circ')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)
        replicon_results_path = self.find_data(os.path.join('Results_Integron_Finder_{}'.format(replicon_name),
                                                            'tmp_{}'.format(replicon.id)))
        attc_file = os.path.join(replicon_results_path, '{}_attc_table.res'.format(replicon.id))
        intI_file = os.path.join(replicon_results_path, '{}_intI.res'.format(replicon.id))
        phageI_file = os.path.join(replicon_results_path, '{}_phage_int.res'.format(replicon.id))

        args = argparse.Namespace()
        args.no_proteins = False
        args.keep_palindromes = True
        args.distance_threshold = 4000
        args.attc_model = 'attc_4.cm'
        args.evalue_attc = 1.0
        args.max_attc_size = 200
        args.min_attc_size = 40
        args.local_max = False
        args.gembase = False
        args.prot_file = False
        args.cmsearch = __file__
        args.hmmsearch = __file__
        args.prodigal = __file__
        args.union_integrases = False
        args.calin_threshold = 2

        cfg = Config(args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')
        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)

        with self.catch_log():
            integrons = find_integron(replicon, prot_db, intI_file, phageI_file, cfg, attc_file=attc_file)

        self.assertEqual(len(integrons), 2)

        args.calin_threshold = 3
        cfg = Config(args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

        with self.catch_log():
            integrons = find_integron(replicon, prot_db, intI_file, phageI_file, cfg, attc_file=attc_file)
        self.assertEqual(len(integrons), 1)


    def test_find_integron_attC_is_df(self):
        replicon_name = 'acba.007.p01.13'
        replicon_id = 'ACBA.007.P01_13'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        prot_file = self.find_data(os.path.join('Proteins', replicon_id + '.prt'))
        topologies = Topology(1, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)
        attc_file = self.find_data(os.path.join('Results_Integron_Finder_{}'.format(replicon_name),
                                                'tmp_{}'.format(replicon.id),
                                                '{}_attc_table.res'.format(replicon.id)))

        intI_file = self.find_data(os.path.join('Results_Integron_Finder_{}'.format(replicon_name),
                                                'tmp_{}'.format(replicon.id),
                                                '{}_intI.res'.format(replicon.id)))
        phageI_file = self.find_data(os.path.join('Results_Integron_Finder_{}'.format(replicon_name),
                                                  'tmp_{}'.format(replicon.id),
                                                  '{}_phage_int.res'.format(replicon.id)))

        args = argparse.Namespace()
        args.gembase = False
        args.prot_file = False
        args.cmsearch = __file__
        args.hmmsearch = __file__
        args.prodigal = __file__
        args.no_proteins = True
        args.keep_palindromes = True
        args.attc_model = 'attc_4.cm'
        args.evalue_attc = 1.0
        args.max_attc_size = 200
        args.min_attc_size = 40
        args.distance_threshold = 4000
        args.calin_threshold = 2
        args.local_max = False
        cfg = Config(args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

        len_model_attc = 47  # length in 'CLEN' (value for model attc_4.cm)

        attc = read_infernal(attc_file,
                             replicon_name,
                             len(replicon),
                             len_model_attc,
                             evalue=cfg.evalue_attc,
                             size_max_attc=cfg.max_attc_size,
                             size_min_attc=cfg.min_attc_size)
        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)

        exp_msg = """In replicon {}, there are:
- 0 complete integron(s) found with a total 0 attC site(s)
- 1 CALIN element(s) found with a total of 3 attC site(s)
- 0 In0 element(s) found with a total of 0 attC site""".format(replicon.id)
        with self.catch_log() as log:
            integrons = find_integron(replicon,
                                      prot_db,
                                      intI_file,
                                      phageI_file,
                                      cfg,
                                      attc=attc)
            catch_msg = log.get_value().strip()
        self.assertEqual(catch_msg, exp_msg)

        self.assertEqual(len(integrons), 1)
        integron = integrons[0]
        self.assertEqual(integron.replicon.name, replicon_id)

        exp = pd.DataFrame({'annotation': ['attC'] * 3,
                            'distance_2attC': [np.nan, 1196.0, 469.0],
                            'evalue': [1.000000e-09, 1.000000e-04, 1.100000e-07],
                            'model': ['attc_4'] * 3,
                            'pos_beg': [17825, 19080, 19618],
                            'pos_end': [17884, 19149, 19726],
                            'strand': [-1, -1, -1],
                            'type_elt': 'attC'},
        columns=self.columns,
        index=['attc_001', 'attc_002', 'attc_003'])
        pdt.assert_frame_equal(integron.attC, exp)

        exp = pd.DataFrame(columns=self.columns)
        exp = exp.astype(dtype=self.dtype)

        pdt.assert_frame_equal(integron.integrase, exp)
        pdt.assert_frame_equal(integron.promoter, exp)
        pdt.assert_frame_equal(integron.attI, exp)
        pdt.assert_frame_equal(integron.proteins, exp)


    def test_find_integron_proteins_lin_replicon(self):
        replicon_name = 'acba.007.p01.13'
        replicon_id = 'ACBA.007.P01_13'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        prot_file = self.find_data(os.path.join('Proteins', replicon_id + '.prt'))
        topologies = Topology(1, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)
        exp_result_dir = 'Results_Integron_Finder_acba.007.p01.13.linear'
        attc_file = self.find_data(os.path.join(exp_result_dir,
                                                'tmp_{}'.format(replicon.id),
                                                '{}_attc_table.res'.format(replicon.id)))
        intI_file = self.find_data(os.path.join(exp_result_dir,
                                                'tmp_{}'.format(replicon.id),
                                                '{}_intI.res'.format(replicon.id)))
        phageI_file = self.find_data(os.path.join(exp_result_dir,
                                                  'tmp_{}'.format(replicon.id),
                                                  '{}_phage_int.res'.format(replicon.id)))
        args = argparse.Namespace()
        args.no_proteins = False
        args.keep_palindromes = True
        args.union_integrases = False
        args.gembase = False  # needed by read_hmm which is called when no_proteins == False

        args = argparse.Namespace()
        args.evalue_attc = 1.
        args.max_attc_size = 200
        args.min_attc_size = 40
        args.distance_threshold = 4000  # (4kb at least between 2 different arrays)
        args.attc_model = 'attc_4.cm'
        args.no_proteins = False
        args.gembase = False  # needed by read_hmm which is called when no_proteins == False
        args.prot_file = False
        args.cmsearch = __file__
        args.hmmsearch = __file__
        args.prodigal = __file__
        args.union_integrases = False
        args.keep_palindromes = True
        args.calin_threshold = 2
        args.local_max = False

        cfg = Config(args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')
        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)

        exp_msg = """In replicon {}, there are:
- 0 complete integron(s) found with a total 0 attC site(s)
- 1 CALIN element(s) found with a total of 3 attC site(s)
- 1 In0 element(s) found with a total of 0 attC site""".format(replicon.id)
        with self.catch_log() as log:
            integrons = find_integron(replicon,
                                      prot_db,
                                      intI_file,
                                      phageI_file,
                                      cfg,
                                      attc_file=attc_file)
            catch_msg = log.get_value().strip()
        self.assertEqual(catch_msg, exp_msg)
        self.assertEqual(len(integrons), 2)

        exp_int = []
        exp = pd.DataFrame({'annotation': 'intI',
                            'distance_2attC': np.nan,
                            'evalue':  1.900000e-25,
                            'model': 'intersection_tyr_intI',
                            'pos_beg': 55,
                            'pos_end': 1014,
                            'strand': 1,
                            'type_elt': 'protein'},
                           columns=self.columns,
                           index=['ACBA.007.P01_13_1'])
        exp = exp.astype(dtype=self.dtype)
        exp_int.append(exp)
        exp_int.append(pd.DataFrame(columns=self.columns).astype(dtype=self.dtype))

        exp_attC = [pd.DataFrame(columns=self.columns).astype(dtype=self.dtype)]
        exp = pd.DataFrame({'annotation': ['attC'] * 3,
                            'distance_2attC': [np.nan, 1196.0,  469.0],
                            'evalue':  [1.000000e-09, 1.000000e-04, 1.100000e-07],
                            'model': ['attc_4'] * 3,
                            'pos_beg': [17825, 19080, 19618],
                            'pos_end': [17884, 19149, 19726],
                            'strand': [-1, -1, -1],
                            'type_elt': 'attC'},
                           columns=self.columns,
                           index=['attc_001', 'attc_002', 'attc_003'])
        exp = exp.astype(dtype=self.dtype)
        exp_attC.append(exp)
        empty = pd.DataFrame(columns=self.columns).astype(dtype=self.dtype)

        for i, integron in enumerate(integrons):
            self.assertEqual(integron.replicon.name, replicon_id)
            pdt.assert_frame_equal(integron.integrase, exp_int[i])
            pdt.assert_frame_equal(integron.attC, exp_attC[i])

            pdt.assert_frame_equal(integron.promoter, empty)
            pdt.assert_frame_equal(integron.attI, empty)
            pdt.assert_frame_equal(integron.proteins, empty)


    def test_find_integron_proteins_circ_replicon(self):
        replicon_name = 'acba.007.p01.13'
        replicon_id = 'ACBA.007.P01_13'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        prot_file = self.find_data(os.path.join('Proteins', replicon_id + '.prt'))
        topologies = Topology(1, 'circ')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)
        exp_result_dir = 'Results_Integron_Finder_acba.007.p01.13.circular'
        attc_file = self.find_data(os.path.join(exp_result_dir,
                                                'tmp_{}'.format(replicon.id),
                                                '{}_attc_table.res'.format(replicon.id)))
        intI_file = self.find_data(os.path.join(exp_result_dir,
                                                'tmp_{}'.format(replicon.id),
                                                '{}_intI.res'.format(replicon.id)))
        phageI_file = self.find_data(os.path.join(exp_result_dir,
                                                  'tmp_{}'.format(replicon.id),
                                                  '{}_phage_int.res'.format(replicon.id)))
        args = argparse.Namespace()
        args.evalue_attc = 1.
        args.max_attc_size = 200
        args.min_attc_size = 40
        args.distance_threshold = 4000  # (4kb at least between 2 different arrays)
        args.attc_model = 'attc_4.cm'
        args.no_proteins = False
        args.gembase = False  # needed by read_hmm which is called when no_proteins == False
        args.prot_file = False
        args.cmsearch = __file__
        args.hmmsearch = __file__
        args.prodigal = __file__

        args.union_integrases = False
        args.keep_palindromes = True
        args.calin_threshold = 2
        args.local_max = False

        cfg = Config(args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)
        exp_msg = """In replicon {}, there are:
- 1 complete integron(s) found with a total 3 attC site(s)
- 0 CALIN element(s) found with a total of 0 attC site(s)
- 0 In0 element(s) found with a total of 0 attC site""".format(replicon.id)
        with self.catch_log() as log:
            integrons = find_integron(replicon,
                                      prot_db,
                                      intI_file,
                                      phageI_file,
                                      cfg,
                                      attc_file=attc_file)
            catch_msg = log.get_value().strip()
        self.assertEqual(catch_msg, exp_msg)
        self.assertEqual(len(integrons), 1)
        integron = integrons[0]
        self.assertEqual(integron.replicon.name, replicon_id)

        exp = pd.DataFrame({'annotation': 'intI',
                            'distance_2attC': np.nan,
                            'evalue':  1.900000e-25,
                            'model': 'intersection_tyr_intI',
                            'pos_beg': 55,
                            'pos_end': 1014,
                            'strand': 1,
                            'type_elt': 'protein'},
                           columns=self.columns,
                           index=['ACBA.007.P01_13_1'])
        exp = exp.astype(dtype=self.dtype)
        pdt.assert_frame_equal(integron.integrase, exp)

        exp = pd.DataFrame({'annotation': ['attC'] * 3,
                            'distance_2attC': [np.nan, 1196.0,  469.0],
                            'evalue':  [1.000000e-09, 1.000000e-04, 1.100000e-07],
                            'model': ['attc_4'] * 3,
                            'pos_beg': [17825, 19080, 19618],
                            'pos_end': [17884, 19149, 19726],
                            'strand': [-1, -1, -1],
                            'type_elt': 'attC'},
                           columns=self.columns,
                           index=['attc_001', 'attc_002', 'attc_003'])
        exp = exp.astype(dtype=self.dtype)
        pdt.assert_frame_equal(integron.attC, exp)

        exp = pd.DataFrame(columns=self.columns)
        exp = exp.astype(dtype=self.dtype)

        pdt.assert_frame_equal(integron.promoter, exp)
        pdt.assert_frame_equal(integron.attI, exp)
        pdt.assert_frame_equal(integron.proteins, exp)


    def test_find_integron_proteins_n_union_integrase(self):
        replicon_name = 'OBAL001.B.00005.C001'
        replicon_id = 'OBAL001.B.00005.C001'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        prot_file = self.find_data(os.path.join('Proteins', replicon_name + '.prt'))
        topologies = Topology(1, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)
        result_dir = 'Results_Integron_Finder_{}.union'.format(replicon_name)
        attc_file = self.find_data(os.path.join(result_dir,
                                                'tmp_{}'.format(replicon.id),
                                                '{}_attc_table.res'.format(replicon.id)))
        intI_file = self.find_data(os.path.join(result_dir,
                                                'tmp_{}'.format(replicon.id),
                                                '{}_intI.res'.format(replicon.id)))
        phageI_file = self.find_data(os.path.join(result_dir,
                                                  'tmp_{}'.format(replicon.id),
                                                  '{}_phage_int.res'.format(replicon.id)))
        args = argparse.Namespace()
        args.evalue_attc = 1.
        args.max_attc_size = 200
        args.min_attc_size = 40
        args.distance_threshold = 4000  # (4kb at least between 2 different arrays)
        args.calin_threshold = 2
        args.attc_model = 'attc_4.cm'
        args.no_proteins = False
        args.keep_palindromes = True
        args.union_integrases = True
        args.gembase = False  # needed by read_hmm which is called when no_proteins == False
        args.prot_file = False
        args.cmsearch = __file__
        args.hmmsearch = __file__
        args.prodigal = __file__
        args.local_max = False
        cfg = Config(args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

        prot_db = ProdigalDB(replicon, cfg, prot_file=prot_file)
        exp_msg = """In replicon {}, there are:
- 3 complete integron(s) found with a total 4 attC site(s)
- 0 CALIN element(s) found with a total of 0 attC site(s)
- 2 In0 element(s) found with a total of 0 attC site""".format(replicon.id)
        with self.catch_log() as log:
            integrons = find_integron(replicon,
                                      prot_db,
                                      intI_file,
                                      phageI_file,
                                      cfg,
                                      attc_file=attc_file)
            catch_msg = log.get_value().strip()
        self.assertEqual(catch_msg, exp_msg)
        self.assertEqual(len(integrons), 5)
        integron = integrons[0]
        self.assertEqual(integron.replicon.name, replicon_id)

        empty = pd.DataFrame(columns=self.columns).astype(dtype=self.dtype)

        exp_int = []
        exp_int.append(pd.DataFrame(
            [[418072, 419283, 1, 5.400000e-25, 'protein', 'Phage_integrase', np.nan, 'intI']],
            columns=self.columns,
            index=['OBAL001.B.00005.C001_388']).astype(dtype=self.dtype))
        exp_int.append(pd.DataFrame(
            [[434671, 440118, -1, 0.085, 'protein', 'Phage_integrase', np.nan, 'intI']],
            columns=self.columns,
            index=['OBAL001.B.00005.C001_399']).astype(dtype=self.dtype))
        exp_int.append(pd.DataFrame(
            [[516941, 517834, -1, 1.200000e-54, 'protein', 'Phage_integrase', np.nan, 'intI']],
            columns=self.columns,
            index=['OBAL001.B.00005.C001_472']).astype(dtype=self.dtype))
        exp_int.append(pd.DataFrame(
            [[1545830, 1546807, -1, 1.100000e-21, 'protein', 'intersection_tyr_intI', np.nan, 'intI']],
            columns=self.columns,
            index=['OBAL001.B.00005.C001_1416']).astype(dtype=self.dtype))
        exp_int.append(pd.DataFrame(
            [[1940269, 1941171, 1, 4.200000e-43, 'protein', 'Phage_integrase', np.nan, 'intI']],
            columns=self.columns,
            index=['OBAL001.B.00005.C001_1793']).astype(dtype=self.dtype))

        exp_attC = []
        exp_attC.append(pd.DataFrame(
            [[421689, 421764, 1, 0.13, 'attC', 'attc_4', np.nan, 'attC']],
            columns=self.columns,
            index=['attc_001']).astype(dtype=self.dtype))
        exp_attC.append(pd.DataFrame(
            [[442458, 442514, -1, 7.000000e-07, 'attC', 'attc_4', np.nan, 'attC']],
            columns=self.columns,
            index=['attc_001']).astype(dtype=self.dtype))
        exp_attC.append(empty)
        exp_attC.append(pd.DataFrame(
            [[1547800, 1547859, 1, 0.00049, 'attC', 'attc_4', np.nan, 'attC'],
             [1548775, 1548834, 1, 0.00009, 'attC', 'attc_4', 916.0, 'attC']],
            columns=self.columns,
            index=['attc_001', 'attc_002']).astype(dtype=self.dtype))
        exp_attC.append(empty)

        for i, integron in enumerate(integrons):
            self.assertEqual(integron.replicon.name, replicon_id)
            pdt.assert_frame_equal(integron.integrase, exp_int[i])
            pdt.assert_frame_equal(integron.attC, exp_attC[i])
            pdt.assert_frame_equal(integron.promoter, empty)
            pdt.assert_frame_equal(integron.attI, empty)
            pdt.assert_frame_equal(integron.proteins, empty)
