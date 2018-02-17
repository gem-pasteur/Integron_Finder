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

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import infernal


_local_max_ori = infernal.local_max


def local_max_mock():
    """expand call several time local_max which run cmserach and need lot of global variable
    to avoid to call cmsearch (we don't want to test local max), I have been ccahed the results
    of local_max and replace the function by this mock
    """
    _cache = {('lian.001.c02.10', 942899, 947099, 'top'):
                  pd.read_csv(TestExpand.find_data('local_max_right_circular_1.csv')),
              ('lian.001.c02.10', 946899, 951099, 'top'):
                  pd.read_csv(TestExpand.find_data('local_max_right_circular_2.csv')),
              ('lian.001.c02.10', 930689, 934889, 'bottom'):
                  pd.read_csv(TestExpand.find_data('local_max_left_circular.csv')),
              ('lian.001.c02.10', 930689, 930889, 'bottom'):
                  pd.read_csv(TestExpand.find_data('local_max_left_linear.csv'))
              }

    def fake_local_max(*args):
        return _cache[args]
    return fake_local_max


class TestExpand(IntegronTest):

    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        os.makedirs(self.tmp_dir)

        infernal.local_max = local_max_mock()


    def tearDown(self):
        infernal.local_max = _local_max_ori
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass


    def test_expand_circular_right(self):
        circular = True
        dist_threshold = 4000
        replicon_name = 'lian.001.c02.10'
        max_attc_size = 200

        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        sequence = SeqIO.read(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)

        max_elt_input = pd.read_csv(os.path.join(self._data_dir, 'max_elt_input_1.csv'))
        df_max_input = pd.read_csv(os.path.join(self._data_dir, 'df_max_input_1.csv'))
        max_elt_expected = pd.read_csv(os.path.join(self._data_dir, 'max_elt_output_lian_right.csv'))
        max_eat_received = infernal.expand(replicon_name, len(sequence),
                                           934689, 943099, max_elt_input, df_max_input,
                                           circular, dist_threshold, max_attc_size,
                                           search_left=False, search_right=True)
        pdt.assert_frame_equal(max_elt_expected, max_eat_received)


    def test_expand_circular_left(self):
        circular = True
        dist_threshold = 4000
        replicon_name = 'lian.001.c02.10'
        max_attc_size = 200

        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        sequence = SeqIO.read(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)

        max_elt_input = pd.read_csv(os.path.join(self._data_dir, 'max_elt_input_1.csv'))
        df_max_input = pd.read_csv(os.path.join(self._data_dir, 'df_max_input_1.csv'))
        max_elt_expected = pd.read_csv(os.path.join(self._data_dir, 'max_elt_output_lian_left.csv'))
        max_eat_received = infernal.expand(replicon_name, len(sequence),
                                           934689, 943099, max_elt_input, df_max_input,
                                           circular, dist_threshold, max_attc_size,
                                           search_left=True, search_right=False)
        pdt.assert_frame_equal(max_elt_expected, max_eat_received)


    def test_expand_linear_right(self):
        circular = False
        dist_threshold = 4000
        replicon_name = 'lian.001.c02.10'
        max_attc_size = 200

        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        sequence = SeqIO.read(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)

        max_elt_input = pd.read_csv(os.path.join(self._data_dir, 'max_elt_input_1.csv'))
        df_max_input = pd.read_csv(os.path.join(self._data_dir, 'df_max_input_1.csv'))
        max_elt_expected = pd.read_csv(os.path.join(self._data_dir, 'max_elt_output_lian_right.csv'))
        max_eat_received = infernal.expand(replicon_name, len(sequence),
                                           934689, 943099, max_elt_input, df_max_input,
                                           circular, dist_threshold, max_attc_size,
                                           search_left=False, search_right=True)
        pdt.assert_frame_equal(max_elt_expected, max_eat_received)


    def test_expand_linear_left(self):
        circular = False
        dist_threshold = 4000
        replicon_name = 'lian.001.c02.10'
        max_attc_size = 200

        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
        sequence = SeqIO.read(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)

        max_elt_input = pd.read_csv(os.path.join(self._data_dir, 'max_elt_input_1.csv'))
        df_max_input = pd.read_csv(os.path.join(self._data_dir, 'df_max_input_1.csv'))
        max_elt_expected = pd.read_csv(os.path.join(self._data_dir, 'max_elt_output_lian_left.csv'))
        max_eat_received = infernal.expand(replicon_name, len(sequence),
                                           934689, 943099, max_elt_input, df_max_input,
                                           circular, dist_threshold, max_attc_size,
                                           search_left=True, search_right=False)
        pdt.assert_frame_equal(max_elt_expected, max_eat_received)
