import os
import tempfile
import shutil
import unittest
import argparse

import pandas as pd
import pandas.util.testing as pdt

# display warning only for non installed integron_finder
from Bio import BiopythonExperimentalWarning
from Bio import Seq, SeqIO
import warnings
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

import integron_finder
from tests import which





class TestFindIntegons(unittest.TestCase):

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


    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_dir)
            pass
        except:
            pass


    def test_find_integron(self):
        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, replicon_name + '.fst')
        attc_file = os.path.join(self._data_dir,
                                 'Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_attc_table.res')
        intI_file = os.path.join(self._data_dir,
                                 'Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_intI.res')
        phageI_file = os.path.join(self._data_dir,
                                   'Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_phage_int.res')
        args = argparse.Namespace
        args.no_proteins = True
        args.keep_palindromes = True

        integron_finder.replicon_name = replicon_name
        integron_finder.SEQUENCE = SeqIO.read(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)
        integron_finder.SIZE_REPLICON = len(integron_finder.SEQUENCE)
        integron_finder.args = args
        integron_finder.evalue_attc = 1.
        integron_finder.max_attc_size = 200
        integron_finder.min_attc_size = 40
        integron_finder.length_cm = 47  # length in 'CLEN' (value for model attc_4.cm)
        integron_finder.DISTANCE_THRESHOLD = 4000  # (4kb at least between 2 different arrays)
        integron_finder.model_attc_name = integron_finder.MODEL_attc.split("/")[-1].split(".cm")[0]

        integrons = integron_finder.find_integron(replicon_name,
                                                  attc_file,
                                                  intI_file,
                                                  phageI_file)
        self.assertEqual(len(integrons), 1)
        integron = integrons[0]
        self.assertEqual(integron.ID_replicon, replicon_name)
        pdt.assert_frame_equal(integron.integrase,
                               pd.read_csv(os.path.join(self._data_dir,'{}_integrase.csv'.format(replicon_name))))
        pdt.assert_frame_equal(integron.attC,
                               pd.read_csv(os.path.join(self._data_dir,'{}_attC.csv'.format(replicon_name))))
        pdt.assert_frame_equal(integron.promoter,
                               pd.read_csv(os.path.join(self._data_dir,'{}_promoter.csv'.format(replicon_name))))
        pdt.assert_frame_equal(integron.attI,
                               pd.read_csv(os.path.join(self._data_dir,'{}_attI.csv'.format(replicon_name))))
        pdt.assert_frame_equal(integron.proteins,
                               pd.read_csv(os.path.join(self._data_dir,'{}_proteins.csv'.format(replicon_name))))
