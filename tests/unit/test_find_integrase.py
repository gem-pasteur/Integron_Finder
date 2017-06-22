import os
import tempfile
import shutil
import unittest
import platform
from collections import namedtuple

# display warning only for non installed integron_finder
from Bio import BiopythonExperimentalWarning
import warnings
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

import integron_finder
_call_ori = integron_finder.call

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


class TestFindIntegrase(unittest.TestCase):

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

        integron_finder.call = call_wrapper()

    def tearDown(self):
        integron_finder.call = _call_ori
        try:
            shutil.rmtree(self.tmp_dir)
            pass
        except:
            pass


    def test_find_integrase_gembase(self):
        FakeArgs = namedtuple('FakeArgs', 'gembase')
        integron_finder.args = FakeArgs(True)

        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, replicon_name + '.fst')

        integron_finder.PROT_file = os.path.join(self.tmp_dir, replicon_name + ".prt")
        shutil.copyfile(os.path.join(self._data_dir, replicon_name + ".prt"),
                        integron_finder.PROT_file)

        integron_finder.find_integrase(replicon_path, replicon_name, self.tmp_dir)
        for suffix in ('_intI.res', '_intI_table.res','_phage_int.res', '_phage_int_table.res'):
            res = os.path.join(self.tmp_dir, replicon_name + suffix)
            self.assertTrue(os.path.exists(res))


    def test_find_integrase_no_gembase_with_protfile(self):
        FakeArgs = namedtuple('FakeArgs', 'gembase')
        integron_finder.args = FakeArgs(False)
        integron_finder.SIZE_REPLICON = 200

        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, replicon_name + '.fst')

        integron_finder.PROT_file = os.path.join(self.tmp_dir, replicon_name + ".prt")
        shutil.copyfile(os.path.join(self._data_dir, replicon_name + ".prt"),
                        integron_finder.PROT_file)

        integron_finder.find_integrase(replicon_path, replicon_name, self.tmp_dir)
        for suffix in ('_intI.res', '_intI_table.res','_phage_int.res', '_phage_int_table.res'):
            res = os.path.join(self.tmp_dir, replicon_name + suffix)
            self.assertTrue(os.path.exists(res))


    def test_find_integrase_no_gembase_no_protfile_short_seq(self):
        FakeArgs = namedtuple('FakeArgs', 'gembase')
        integron_finder.args = FakeArgs(False)
        integron_finder.SIZE_REPLICON = 200

        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, replicon_name + '.fst')

        integron_finder.PROT_file = os.path.join(self.tmp_dir, replicon_name + ".prt")

        integron_finder.find_integrase(replicon_path, replicon_name, self.tmp_dir)
        for suffix in ('_intI.res', '_intI_table.res','_phage_int.res', '_phage_int_table.res'):
            res = os.path.join(self.tmp_dir, replicon_name + suffix)
            self.assertTrue(os.path.exists(res))


    def test_find_integrase_no_gembase_no_protfile_long_seq(self):
        FakeArgs = namedtuple('FakeArgs', 'gembase')
        integron_finder.args = FakeArgs(False)
        integron_finder.SIZE_REPLICON = 500000

        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, replicon_name + '.fst')

        integron_finder.PROT_file = os.path.join(self.tmp_dir, replicon_name + ".prt")
        integron_finder.find_integrase(replicon_path, replicon_name, self.tmp_dir)
        for suffix in ('_intI.res', '_intI_table.res','_phage_int.res', '_phage_int_table.res'):
            res = os.path.join(self.tmp_dir, replicon_name + suffix)
            self.assertTrue(os.path.exists(res))


    def test_find_integrase_no_gembase_no_protfile_no_prodigal(self):
        FakeArgs = namedtuple('FakeArgs', 'gembase')
        integron_finder.args = FakeArgs(False)
        integron_finder.SIZE_REPLICON = 500000

        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, replicon_name + '.fst')

        integron_finder.PRODIGAL = 'foo'
        with self.assertRaises(RuntimeError) as ctx:
            integron_finder.find_integrase(replicon_path, replicon_name, self.tmp_dir)
        self.assertEqual(str(ctx.exception), 'foo failed : [Errno 2] No such file or directory')


    def test_find_integrase_no_gembase_no_protfile_no_replicon(self):
        FakeArgs = namedtuple('FakeArgs', 'gembase')
        integron_finder.args = FakeArgs(False)
        integron_finder.SIZE_REPLICON = 500000

        replicon_name = 'acba.007.p01.13'

        with self.assertRaises(RuntimeError) as ctx:
            integron_finder.find_integrase('foo', replicon_name, self.tmp_dir)
        self.assertEqual(str(ctx.exception), '{} failed returncode = 5'.format(integron_finder.PRODIGAL))


    def test_find_integrase_gembase_no_hmmer(self):
        FakeArgs = namedtuple('FakeArgs', 'gembase')
        integron_finder.args = FakeArgs(True)
        integron_finder.HMMSEARCH = 'foo'

        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, replicon_name + '.fst')

        integron_finder.PROT_file = os.path.join(self.tmp_dir, replicon_name + ".prt")
        shutil.copyfile(os.path.join(self._data_dir, replicon_name + ".prt"),
                        integron_finder.PROT_file)
        with self.assertRaises(RuntimeError) as ctx:
            integron_finder.find_integrase(replicon_path, replicon_name, self.tmp_dir)

        self.assertEqual(str(ctx.exception), 'foo failed : [Errno 2] No such file or directory')


    def test_find_integrase_gembase_no_hmmer(self):
        FakeArgs = namedtuple('FakeArgs', 'gembase')
        integron_finder.args = FakeArgs(True)
        integron_finder.HMMSEARCH = 'foo'

        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, replicon_name + '.fst')

        integron_finder.PROT_file = os.path.join(self.tmp_dir, replicon_name + ".prt")
        shutil.copyfile(os.path.join(self._data_dir, replicon_name + ".prt"),
                        integron_finder.PROT_file)
        with self.assertRaises(RuntimeError) as ctx:
            integron_finder.find_integrase(replicon_path, replicon_name, self.tmp_dir)
        self.assertEqual(str(ctx.exception), 'foo failed : [Errno 2] No such file or directory')


    def test_find_integrase_gembase_hmmer_error(self):
        FakeArgs = namedtuple('FakeArgs', 'gembase')
        integron_finder.args = FakeArgs(True)
        integron_finder.N_CPU = 'foo'

        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, replicon_name + '.fst')

        integron_finder.PROT_file = os.path.join(self.tmp_dir, replicon_name + ".prt")
        shutil.copyfile(os.path.join(self._data_dir, replicon_name + ".prt"),
                        integron_finder.PROT_file)
        with self.assertRaises(RuntimeError) as ctx:
            integron_finder.find_integrase(replicon_path, replicon_name, self.tmp_dir)
        self.assertEqual(str(ctx.exception),  '{} failed return code = 1'.format(integron_finder.HMMSEARCH))


