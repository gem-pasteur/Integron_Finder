import os
import unittest
import tempfile
import shutil

from Bio import SeqIO, Seq
import pandas as pd
import pandas.util.testing as pdt


from Bio import BiopythonExperimentalWarning
import warnings
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

import integron_finder
_local_max_ori = integron_finder.local_max


def local_max_mock():
    """expand call several time local_max which run cmserach and need lot of global variable
    to avoid to call cmsearch (we don't want to test local max), I have been ccahed the results
    of local_max and replace the function by this mock
    """
    _cache = {('lian.001.c02.10', 942899, 947099, 'top'): pd.read_csv(os.path.join(TestExpand._data_dir,
                                                                                   'local_max_right_circular_1.csv')),
              ('lian.001.c02.10', 946899, 951099, 'top'): pd.read_csv(os.path.join(TestExpand._data_dir,
                                                                                   'local_max_right_circular_2.csv')),
              ('lian.001.c02.10', 930689, 934889, 'bottom'): pd.read_csv(os.path.join(TestExpand._data_dir,
                                                                                      'local_max_left_circular.csv')),
              ('lian.001.c02.10', 930689, 930889, 'bottom'): pd.read_csv(os.path.join(TestExpand._data_dir,
                                                                                      'local_max_left_linear.csv'))
              }
    def fake_local_max(*args):
        return _cache[args]
    return fake_local_max


class TestExpand(unittest.TestCase):

    _data_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', "data"))

    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..' '..')))

        integron_finder.max_attc_size = 200
        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        os.makedirs(self.tmp_dir)

        integron_finder.local_max = local_max_mock()


    def tearDown(self):
        integron_finder.local_max = _local_max_ori
        try:
            shutil.rmtree(self.tmp_dir)
            pass
        except:
            pass


    def test_expand_circular_right(self):
        integron_finder.circular = True
        integron_finder.DISTANCE_THRESHOLD = 4000
        integron_finder.replicon_name = 'lian.001.c02.10'

        replicon_path = os.path.join(self._data_dir, integron_finder.replicon_name +'.fst')
        sequence = SeqIO.read(replicon_path, "fasta",
                                              alphabet=Seq.IUPAC.unambiguous_dna)
        integron_finder.SIZE_REPLICON = len(sequence)

        max_elt_input = pd.read_csv(os.path.join(self._data_dir, 'max_elt_input_1.csv'))
        df_max_input = pd.read_csv(os.path.join(self._data_dir, 'df_max_input_1.csv'))
        max_elt_expected = pd.read_csv(os.path.join(self._data_dir, 'max_elt_output_lian_right.csv'))
        max_eat_received = integron_finder.expand(934689, 943099, max_elt_input, df_max_input,
                                                  search_left=False, search_right=True)
        pdt.assert_frame_equal(max_elt_expected, max_eat_received)


    def test_expand_circular_left(self):
        integron_finder.circular = True
        integron_finder.DISTANCE_THRESHOLD = 4000
        integron_finder.replicon_name = 'lian.001.c02.10'

        replicon_path = os.path.join(self._data_dir, integron_finder.replicon_name + '.fst')
        sequence = SeqIO.read(replicon_path, "fasta",
                              alphabet=Seq.IUPAC.unambiguous_dna)
        integron_finder.SIZE_REPLICON = len(sequence)

        max_elt_input = pd.read_csv(os.path.join(self._data_dir, 'max_elt_input_1.csv'))
        df_max_input = pd.read_csv(os.path.join(self._data_dir, 'df_max_input_1.csv'))
        max_elt_expected = pd.read_csv(os.path.join(self._data_dir, 'max_elt_output_lian_left.csv'))
        max_eat_received = integron_finder.expand(934689, 943099, max_elt_input, df_max_input,
                                                  search_left=True, search_right=False)
        pdt.assert_frame_equal(max_elt_expected, max_eat_received)


    def test_expand_linear_right(self):
        integron_finder.circular = False
        integron_finder.DISTANCE_THRESHOLD = 4000
        integron_finder.replicon_name = 'lian.001.c02.10'

        replicon_path = os.path.join(self._data_dir, integron_finder.replicon_name +'.fst')
        sequence = SeqIO.read(replicon_path, "fasta",
                                              alphabet=Seq.IUPAC.unambiguous_dna)
        integron_finder.SIZE_REPLICON = len(sequence)

        max_elt_input = pd.read_csv(os.path.join(self._data_dir, 'max_elt_input_1.csv'))
        df_max_input = pd.read_csv(os.path.join(self._data_dir, 'df_max_input_1.csv'))
        max_elt_expected = pd.read_csv(os.path.join(self._data_dir, 'max_elt_output_lian_right.csv'))
        max_eat_received = integron_finder.expand(934689, 943099, max_elt_input, df_max_input,
                                                  search_left=False, search_right=True)
        pdt.assert_frame_equal(max_elt_expected, max_eat_received)


    def test_expand_linear_left(self):
        integron_finder.circular = False
        integron_finder.DISTANCE_THRESHOLD = 4000
        integron_finder.replicon_name = 'lian.001.c02.10'

        replicon_path = os.path.join(self._data_dir, integron_finder.replicon_name + '.fst')
        sequence = SeqIO.read(replicon_path, "fasta",
                              alphabet=Seq.IUPAC.unambiguous_dna)
        integron_finder.SIZE_REPLICON = len(sequence)

        max_elt_input = pd.read_csv(os.path.join(self._data_dir, 'max_elt_input_1.csv'))
        df_max_input = pd.read_csv(os.path.join(self._data_dir, 'df_max_input_1.csv'))
        max_elt_expected = pd.read_csv(os.path.join(self._data_dir, 'max_elt_output_lian_left.csv'))
        max_eat_received = integron_finder.expand(934689, 943099, max_elt_input, df_max_input,
                                                  search_left=True, search_right=False)
        pdt.assert_frame_equal(max_elt_expected, max_eat_received)
