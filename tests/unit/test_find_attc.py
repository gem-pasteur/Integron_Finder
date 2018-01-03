import os
import tempfile
import shutil

# display warning only for non installed integron_finder
from Bio import BiopythonExperimentalWarning
import warnings
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)


import unittest
import integron_finder
from tests import which

class TestFindAttc(unittest.TestCase):

    _data_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', "data"))

    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        os.makedirs(self.tmp_dir)
        integron_finder.CMSEARCH = which('cmsearch')
        integron_finder.N_CPU = '1'
        integron_finder.MODEL_attc = os.path.join(self.integron_home, 'data', 'Models', 'attc_4.cm')


    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass

    def test_find_attc(self):
        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, 'Replicons', replicon_name + '.fst')
        integron_finder.find_attc(replicon_path, replicon_name, self.tmp_dir)
        for suffix in ('_attc.res', '_attc_table.res'):
            res = os.path.join(self.tmp_dir, replicon_name + suffix)
            self.assertTrue(os.path.exists(res))

    def test_find_attc_no_infernal(self):
        integron_finder.CMSEARCH = 'foo'
        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, 'Replicons', replicon_name + '.fst')
        with self.assertRaises(RuntimeError) as ctx:
            integron_finder.find_attc(replicon_path, replicon_name, self.tmp_dir)
        self.assertEqual(str(ctx.exception), 'foo failed : [Errno 2] No such file or directory')

    def test_find_attc_no_model(self):
        integron_finder.MODEL_attc = 'foo'
        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, 'Replicons', replicon_name + '.fst')
        with self.assertRaises(RuntimeError) as ctx:
            integron_finder.find_attc(replicon_path, replicon_name, self.tmp_dir)
        self.assertEqual(str(ctx.exception), '{} failed returncode = 1'.format(integron_finder.CMSEARCH))
