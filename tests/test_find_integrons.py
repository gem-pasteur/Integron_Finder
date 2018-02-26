import os
import sys
import tempfile
import shutil
import argparse
from StringIO import StringIO
from contextlib import contextmanager

import numpy as np
import pandas as pd
import pandas.util.testing as pdt

# display warning only for non installed integron_finder
from Bio import BiopythonExperimentalWarning
import warnings
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.integron import Integron, find_integron
from integron_finder.config import Config
from integron_finder.utils import read_single_dna_fasta
from integron_finder.infernal import read_infernal


class TestFindIntegons(IntegronTest):


    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        os.makedirs(self.tmp_dir)

        # by default use empty Namespace just to match with api
        # if need some values create a local config with these values
        # this rule should avoid side effects
        args = argparse.Namespace()
        self.cfg = Config(args)
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

    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_dir)
            pass
        except:
            pass


    @contextmanager
    def catch_output(self):
        """
        Catch stderr and stdout of the code running with this function.
        They can, then, be compared to expected outputs.
        """
        new_out, new_err = StringIO(), StringIO()
        old_out, old_err = sys.stdout, sys.stderr
        try:
            sys.stdout, sys.stderr = new_out, new_err
            yield sys.stdout, sys.stderr
        finally:
            sys.stdout, sys.stderr = old_out, old_err


    def test_find_integron(self):
        replicon_name = 'acba.007.p01.13'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        replicon_results_path = self.find_data(os.path.join('Results_Integron_Finder_{}'.format(replicon_name), 'other'))
        attc_file = os.path.join(replicon_results_path, '{}_attc_table.res'.format(replicon_name))
        intI_file = os.path.join(replicon_results_path, '{}_intI.res'.format(replicon_name))
        phageI_file = os.path.join(replicon_results_path, '{}_phage_int.res'.format(replicon_name))

        replicon = read_single_dna_fasta(replicon_path)

        args = argparse.Namespace()
        args.no_proteins = True
        args.keep_palindromes = True
        args.distance_threshold = 4000
        args.attc_model = 'attc_4.cm'
        args.evalue_attc = 1.0
        args.max_attc_size = 200
        args.min_attc_size = 40
        cfg = Config(args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')
        with self.catch_output() as (out, err):
            integrons = find_integron(replicon, attc_file, intI_file, phageI_file, cfg)

        self.assertEqual(err.getvalue().strip(), "")
        self.assertEqual(out.getvalue().strip(), """In replicon acba.007.p01.13, there are:
- 0 complete integron(s) found with a total 0 attC site(s)
- 1 CALIN element(s) found with a total of 3 attC site(s)
- 0 In0 element(s) found with a total of 0 attC site""")

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


    def test_find_integron_attC_is_df(self):
        replicon_name = 'acba.007.p01.13'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        replicon = read_single_dna_fasta(replicon_path)
        attc_file = self.find_data(os.path.join('Results_Integron_Finder_{}'.format(replicon_name),
                                                'other', '{}_attc_table.res'.format(replicon_name)))

        intI_file = self.find_data(os.path.join('Results_Integron_Finder_{}'.format(replicon_name),
                                                'other', '{}_intI.res'.format(replicon_name)))
        phageI_file = self.find_data(os.path.join('Results_Integron_Finder_{}'.format(replicon_name),
                                                  'other', '{}_phage_int.res'.format(replicon_name)))

        # integron_finder.replicon_name = replicon_name
        # integron_finder.SEQUENCE = SeqIO.read(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)
        # integron_finder.SIZE_REPLICON = len(integron_finder.SEQUENCE)
        # integron_finder.args = args
        # integron_finder.DISTANCE_THRESHOLD = 4000  # (4kb at least between 2 different arrays)
        # integron_finder.model_attc_name = integron_finder.MODEL_attc.split("/")[-1].split(".cm")[0]

        args = argparse.Namespace()
        args.no_proteins = True
        args.keep_palindromes = True
        args.attc_model = 'attc_4.cm'
        args.evalue_attc = 1.0
        args.max_attc_size = 200
        args.min_attc_size = 40
        args.distance_threshold = 4000
        cfg = Config(args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')
        len_model_attc = 47  # length in 'CLEN' (value for model attc_4.cm)

        attc_file = read_infernal(attc_file, replicon_name,
                                  len_model_attc,
                                  evalue=cfg.evalue_attc,
                                  size_max_attc=cfg.max_attc_size,
                                  size_min_attc=cfg.min_attc_size)

        with self.catch_output() as (out, err):
            integrons = find_integron(replicon,
                                      attc_file,
                                      intI_file,
                                      phageI_file,
                                      cfg)

        self.assertEqual(err.getvalue().strip(), "")
        self.assertEqual(out.getvalue().strip(), """In replicon acba.007.p01.13, there are:
- 0 complete integron(s) found with a total 0 attC site(s)
- 1 CALIN element(s) found with a total of 3 attC site(s)
- 0 In0 element(s) found with a total of 0 attC site""")

        self.assertEqual(len(integrons), 1)
        integron = integrons[0]
        self.assertEqual(integron.replicon.name, replicon_name)

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


    def test_find_integron_proteins(self):
        replicon_name = 'acba.007.p01.13'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        replicon = read_single_dna_fasta(replicon_path)
        attc_file = os.path.join(self._data_dir,
                                 'Results_Integron_Finder_acba.007.p01.13', 'other', 'acba.007.p01.13_attc_table.res')
        intI_file = os.path.join(self._data_dir,
                                 'Results_Integron_Finder_acba.007.p01.13', 'other', 'acba.007.p01.13_intI.res')
        phageI_file = os.path.join(self._data_dir,
                                   'Results_Integron_Finder_acba.007.p01.13', 'other', 'acba.007.p01.13_phage_int.res')
        args = argparse.Namespace()
        args.no_proteins = False
        args.keep_palindromes = True
        args.union_integrases = False
        args.gembase = False  # needed by read_hmm which is called when no_proteins == False

        replicon = read_single_dna_fasta(replicon_path)
        args = argparse.Namespace()
        args.evalue_attc = 1.
        args.max_attc_size = 200
        args.min_attc_size = 40
        args.distance_threshold = 4000  # (4kb at least between 2 different arrays)
        args.attc_model = 'attc_4.cm'
        args.no_proteins = False
        args.gembase = False  # needed by read_hmm which is called when no_proteins == False
        args.union_integrases = False
        args.keep_palindromes = True
        cfg = Config(args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

        with self.catch_output() as (out, err):
            integrons = find_integron(replicon,
                                      attc_file,
                                      intI_file,
                                      phageI_file,
                                      cfg)

        self.assertEqual(err.getvalue().strip(), "")
        self.assertEqual(out.getvalue().strip(), """In replicon acba.007.p01.13, there are:
- 1 complete integron(s) found with a total 3 attC site(s)
- 0 CALIN element(s) found with a total of 0 attC site(s)
- 0 In0 element(s) found with a total of 0 attC site""")

        self.assertEqual(len(integrons), 1)
        integron = integrons[0]
        self.assertEqual(integron.replicon.name, replicon_name)

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
        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, 'Replicons', replicon_name + '.fst')
        attc_file = os.path.join(self._data_dir,
                                 'Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_attc_table.res')
        intI_file = os.path.join(self._data_dir,
                                 'Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_intI.res')
        phageI_file = os.path.join(self._data_dir,
                                   'Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_phage_int.res')
        args = argparse.Namespace()
        args.no_proteins = False
        args.keep_palindromes = True
        args.union_integrases = True
        args.gembase = False  # needed by read_hmm which is called when no_proteins == False

        replicon = read_single_dna_fasta(replicon_path)

        args = argparse.Namespace()
        args.evalue_attc = 1.
        args.max_attc_size = 200
        args.min_attc_size = 40
        args.distance_threshold = 4000  # (4kb at least between 2 different arrays)
        args.attc_model = 'attc_4.cm'
        args.no_proteins = False
        args.keep_palindromes = True
        args.union_integrases = True
        args.gembase = False  # needed by read_hmm which is called when no_proteins == False
        cfg = Config(args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

        with self.catch_output() as (out, err):
            integrons = find_integron(replicon,
                                      attc_file,
                                      intI_file,
                                      phageI_file,
                                      cfg)

        self.assertEqual(out.getvalue().strip(), """In replicon acba.007.p01.13, there are:
- 1 complete integron(s) found with a total 3 attC site(s)
- 0 CALIN element(s) found with a total of 0 attC site(s)
- 0 In0 element(s) found with a total of 0 attC site""")

        self.assertEqual(len(integrons), 1)
        integron = integrons[0]
        self.assertEqual(integron.replicon.name, replicon_name)

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
