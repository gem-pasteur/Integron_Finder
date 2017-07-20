import os
import tempfile
import shutil
import unittest
import argparse

import numpy as np

import pandas as pd
import pandas.util.testing as pdt

# display warning only for non installed integron_finder
from Bio import BiopythonExperimentalWarning
from Bio import Seq, SeqIO
import warnings
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

import integron_finder
from integron_finder import Integron

from tests import which


class TestFindFindAttCMax(unittest.TestCase):

    _data_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', "data"))


    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..' '..')))

        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        os.makedirs(self.tmp_dir)

        integron_finder.PRODIGAL = which('prodigal')
        integron_finder.HMMSEARCH = which('hmmsearch')
        integron_finder.N_CPU = '1'
        integron_finder.MODEL_DIR = os.path.join(self.integron_home, "data", "Models")
        integron_finder.MODEL_integrase = os.path.join(integron_finder.MODEL_DIR, "integron_integrase.hmm")
        integron_finder.MODEL_phage_int = os.path.join(integron_finder.MODEL_DIR, "phage-int.hmm")
        integron_finder.MODEL_attc = os.path.join(self.integron_home, 'data', 'Models', 'attc_4.cm')

        integron_finder.circular = True
        integron_finder.out_dir = self.tmp_dir
        integron_finder.CMSEARCH = which('cmsearch')
        integron_finder.evalue_attc = 1.
        integron_finder.max_attc_size = 200
        integron_finder.min_attc_size = 40
        integron_finder.length_cm = 47  # length in 'CLEN' (value for model attc_4.cm)
        integron_finder.DISTANCE_THRESHOLD = 4000  # (4kb at least between 2 different arrays)

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
        try:
            shutil.rmtree(self.tmp_dir)
            pass
        except:
            pass


    def test_find_attc_max_linear(self):
        replicon_name = 'OBAL001.B.00005.C001'
        replicon_path = os.path.join(self._data_dir, 'Replicons', replicon_name + '.fst')

        args = argparse.Namespace
        args.eagle_eyes = False
        args.local_max = False

        integron_finder.replicon_name = replicon_name
        integron_finder.SEQUENCE = SeqIO.read(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)
        integron_finder.SIZE_REPLICON = len(integron_finder.SEQUENCE)
        integron_finder.args = args
        integron_finder.circular = False

        integron = Integron(replicon_name)
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
        integron.integrase = integrase

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
        integron.attC = attC
        integrons = [integron]

        max_final = integron_finder.find_attc_max(integrons)

        exp = pd.DataFrame({'Accession_number': ['OBAL001.B.00005.C001', 'OBAL001.B.00005.C001'],
                            'cm_attC': ['attC_4', 'attC_4'],
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
        replicon_name = 'OBAL001.B.00005.C001'
        replicon_path = os.path.join(self._data_dir, 'Replicons', replicon_name + '.fst')

        args = argparse.Namespace
        args.eagle_eyes = False
        args.local_max = False

        integron_finder.replicon_name = replicon_name
        integron_finder.SEQUENCE = SeqIO.read(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)
        integron_finder.SIZE_REPLICON = len(integron_finder.SEQUENCE)
        integron_finder.args = args

        integron = Integron(replicon_name)
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
        integron.integrase = integrase

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
        integron.attC = attC
        integrons = [integron]

        max_final = integron_finder.find_attc_max(integrons)

        exp = pd.DataFrame({'Accession_number': ['OBAL001.B.00005.C001', 'OBAL001.B.00005.C001'],
                            'cm_attC': ['attC_4', 'attC_4'],
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
        replicon_name = 'OBAL001.B.00005.C001'
        replicon_path = os.path.join(self._data_dir, 'Replicons', replicon_name + '.fst')

        args = argparse.Namespace
        args.eagle_eyes = False
        args.local_max = False

        integron_finder.replicon_name = replicon_name
        integron_finder.SEQUENCE = SeqIO.read(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)
        integron_finder.SIZE_REPLICON = len(integron_finder.SEQUENCE)
        integron_finder.args = args

        integron = Integron(replicon_name)

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
        integron.attC = attC
        integrons = [integron]

        max_final = integron_finder.find_attc_max(integrons)

        exp = pd.DataFrame({'Accession_number': ['OBAL001.B.00005.C001'],
                            'cm_attC': ['attC_4'],
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

        pdt.assert_frame_equal(max_final, exp)