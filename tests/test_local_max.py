import os
import tempfile
import shutil

import pandas as pd
import pandas.util.testing as pdt

# display warning only for non installed integron_finder
from Bio import BiopythonExperimentalWarning
from Bio import Seq, SeqIO
import warnings
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import infernal

_call_ori = infernal.call
_read_infernal_ori = infernal.read_infernal

from tests import which


def call_wrapper():
    """
    hmmsearch or prodigal write lot of things on stderr or stdout 
    which noise the unit test output
    So I replace the `call` function in module integron_finder
    by a wrapper which call the original function but add redirect stderr and stdout
    in dev_null
    :return: wrapper around integron_finder.call
    :rtype: function
    """
    def wrapper(*args, **kwargs):
        with open(os.devnull, 'w') as f:
            kwargs['stderr'] = f
            kwargs['stdout'] = f
            res = _call_ori(*args, **kwargs)
        return res
    return wrapper


def read_infernal_mock(tmp_dir):
    """expand call several time local_max which run cmserach and need lot of global variable
    to avoid to call cmsearch (we don't want to test local max), I have been ccahed the results
    of local_max and replace the function by this mock
    """
    _cache = {
        (os.path.join(tmp_dir, 'lian.001.c02.10_942899_947099_subseq_attc_table.res'),
         'lian.001.c02.10', 47, 1.0, 200, 40):
            pd.DataFrame([['lian.001.c02.10', 'attC_4', 1, 47, 371, 496, '+', 0.130000],
                          ['lian.001.c02.10', 'attC_4', 1, 47, 1109, 1234, '+', 0.049000],
                          ['lian.001.c02.10', 'attC_4', 1, 47, 1573, 1699, '+', 0.000005]],
                         columns=['Accession_number', 'cm_attC', 'cm_debut',
                                  'cm_fin', 'pos_beg', 'pos_end', 'sens', 'evalue']),
        (os.path.join(tmp_dir, 'lian.001.c02.10_946899_951099_subseq_attc_table.res'),
         'lian.001.c02.10', 47, 1.0, 200, 40):
            pd.DataFrame(columns=['Accession_number', 'cm_attC', 'cm_debut',
                         'cm_fin', 'pos_beg', 'pos_end', 'sens', 'evalue']),
        (os.path.join(tmp_dir, 'lian.001.c02.10_930689_934889_subseq_attc_table.res'),
         'lian.001.c02.10', 47, 1.0, 200, 40):
            pd.DataFrame(columns=['Accession_number', 'cm_attC', 'cm_debut',
                         'cm_fin', 'pos_beg', 'pos_end', 'sens', 'evalue']),
              }

    def fake_read_infernal(tblout_path, replicon_name, len_model_attc,
                           evalue=None, size_max_attc=None, size_min_attc=None):
        args = (tblout_path, replicon_name, len_model_attc, evalue, size_max_attc, size_min_attc)
        return _cache[args]
    return fake_read_infernal


class TestLocalMax(IntegronTest):

    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        os.makedirs(self.tmp_dir)

        self.PRODIGAL = which('prodigal')
        self.HMMSEARCH = which('hmmsearch')
        self.cmsearch = which('cmsearch')
        self.out_dir = self.tmp_dir
        self.model_attc = self.find_data(os.path.join('Models', 'attc_4.cm'))
        self.cpu_nb = 1
        self.replicon_name = 'lian.001.c02.10'
        self.replicon_path = self.find_data(os.path.join('Replicons', self.replicon_name + '.fst'))
        self.sequence = SeqIO.read(self.replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)
        self.evalue_attc = 1.
        self.max_attc_size = 200
        self.min_attc_size = 40
        self.length_cm = 47  # length in 'CLEN' (value for model attc_4.cm)
        self.call = call_wrapper()
        infernal.read_infernal = read_infernal_mock(self.tmp_dir)

    def tearDown(self):
        infernal.call = _call_ori
        infernal.read_infernal = _read_infernal_ori
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass


    def test_local_max_top(self):
        win_beg = 942899
        win_end = 947099
        strand_search = 'top'
        local_max_received = infernal.local_max(self.replicon_name, self.sequence,
                                                win_beg, win_end, self.length_cm,
                                                model_attc=self.model_attc, strand_search=strand_search,
                                                evalue_attc=self.evalue_attc,
                                                max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                                cmsearch_bin=self.cmsearch, out_dir=self.out_dir, cpu_nb=self.cpu_nb
                                                )
        local_max_expected = pd.DataFrame([['lian.001.c02.10', 'attC_4', 1, 47, 943270, 943395, '+', 0.13],
                                           ['lian.001.c02.10', 'attC_4', 1, 47, 944008, 944133, '+', 0.049],
                                           ['lian.001.c02.10', 'attC_4', 1, 47, 944472, 944598, '+', 4.7e-06]],
                                          columns=['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg',
                                                   'pos_end', 'sens', 'evalue'])
        pdt.assert_frame_equal(local_max_expected, local_max_received)


    def test_no_local_max_top(self):
        win_beg = 946899
        win_end = 951099
        strand_search = 'top'
        local_max_received = infernal.local_max(self.replicon_name, self.sequence,
                                                win_beg, win_end, self.length_cm,
                                                model_attc=self.model_attc, strand_search=strand_search,
                                                evalue_attc=self.evalue_attc,
                                                max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                                cmsearch_bin=self.cmsearch, out_dir=self.out_dir, cpu_nb=self.cpu_nb
                                                )
        local_max_expected = pd.DataFrame(columns=['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg',
                                                   'pos_end', 'sens', 'evalue'])
        pdt.assert_frame_equal(local_max_expected, local_max_received)


    def test_no_local_max_bottom(self):
        win_beg = 946899
        win_end = 951099
        strand_search = 'bottom'
        local_max_received = infernal.local_max(self.replicon_name, self.sequence,
                                                win_beg, win_end, self.length_cm,
                                                model_attc=self.model_attc, strand_search=strand_search,
                                                evalue_attc=self.evalue_attc,
                                                max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                                cmsearch_bin=self.cmsearch, out_dir=self.out_dir, cpu_nb=self.cpu_nb
                                                )
        local_max_expected = pd.DataFrame(columns=['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg',
                                                   'pos_end', 'sens', 'evalue'])
        pdt.assert_frame_equal(local_max_expected, local_max_received)


    def test_local_max_both(self):
        win_beg = 946899
        win_end = 951099
        strand_search = 'both'

        local_max_received = infernal.local_max(self.replicon_name, self.sequence,
                                                win_beg, win_end, self.length_cm,
                                                model_attc=self.model_attc, strand_search=strand_search,
                                                evalue_attc=self.evalue_attc,
                                                max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                                cmsearch_bin=self.cmsearch, out_dir=self.out_dir, cpu_nb=self.cpu_nb
                                                )
        local_max_expected = pd.DataFrame([['lian.001.c02.10', 'attC_4', 1, 47, 943270, 943395, '+', 0.13],
                                           ['lian.001.c02.10', 'attC_4', 1, 47, 944008, 944133, '+', 0.049],
                                           ['lian.001.c02.10', 'attC_4', 1, 47, 944472, 944598, '+', 4.7e-06]],
                                          columns=['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg',
                                                   'pos_end', 'sens', 'evalue'])
        pdt.assert_frame_equal(local_max_expected, local_max_received)


    def local_max_end_beg(self):
        win_beg = 946899
        win_end = 951099
        strand_search = 'top'
        local_max_received = infernal.local_max(self.replicon_name, self.sequence,
                                                win_beg, win_end, self.length_cm,
                                                model_attc=self.model_attc, strand_search=strand_search,
                                                evalue_attc=self.evalue_attc,
                                                max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                                cmsearch_bin=self.cmsearch, out_dir=self.out_dir, cpu_nb=self.cpu_nb
                                                )
        local_max_expected = pd.DataFrame([['lian.001.c02.10', 'attC_4', 1, 47, 943270, 943395, '+', 0.13],
                                           ['lian.001.c02.10', 'attC_4', 1, 47, 944008, 944133, '+', 0.049],
                                           ['lian.001.c02.10', 'attC_4', 1, 47, 944472, 944598, '+', 4.7e-06]],
                                          columns=['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg',
                                                   'pos_end', 'sens', 'evalue'])
        pdt.assert_frame_equal(local_max_expected, local_max_received)

    def test_local_max_cmsearch_failed(self):
        win_beg = 946899
        win_end = 951099
        strand_search = 'top'

        cmsearch_bin = 'failed_cmsearch'
        with self.assertRaises(RuntimeError) as ctx:
            _ = infernal.local_max(self.replicon_name, self.sequence,
                                   win_beg, win_end, self.length_cm,
                                   model_attc=self.model_attc, strand_search=strand_search,
                                   evalue_attc=self.evalue_attc,
                                   max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                   cmsearch_bin=cmsearch_bin, out_dir=self.out_dir, cpu_nb=self.cpu_nb
                                   )
        self.assertEqual(ctx.exception.message,
                         "{} failed : [Errno 2] No such file or directory".format(cmsearch_bin))

        infernal.call = lambda x: 1
        with self.assertRaises(RuntimeError) as ctx:
            _ = infernal.local_max(self.replicon_name, self.sequence,
                                   win_beg, win_end, self.length_cm,
                                   model_attc=self.model_attc, strand_search=strand_search,
                                   evalue_attc=self.evalue_attc,
                                   max_attc_size=self.max_attc_size, min_attc_size=self.min_attc_size,
                                   cmsearch_bin=cmsearch_bin, out_dir=self.out_dir, cpu_nb=self.cpu_nb
                                   )
        self.assertEqual(ctx.exception.message, "{} failed returncode = {}".format(cmsearch_bin,
                                                                                   infernal.call(None)))
