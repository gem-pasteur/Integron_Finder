import os
import tempfile
import shutil
import distutils

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

from integron_finder import attc

_call_ori = attc.call


class TestFindAttc(IntegronTest):

    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        os.makedirs(self.tmp_dir)
        self.cmsearch_path = distutils.spawn.find_executable("cmsearch")
        self.cpu_nb = 1
        self.model_attc = self.find_data(os.path.join('Models', 'attc_4.cm'))
        self.replicon_name = 'acba.007.p01.13'
        self.replicon_path = self.find_data(os.path.join('Replicons', self.replicon_name + '.fst'))
        attc.call = self.call_wrapper(_call_ori)

    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass
        attc.call = _call_ori

    def test_find_attc(self):
        attc.find_attc(self.replicon_path, self.replicon_name, self.cmsearch_path, self.tmp_dir, self.model_attc)
        for suffix in ('_attc.res', '_attc_table.res'):
            res = os.path.join(self.tmp_dir, self.replicon_name + suffix)
            self.assertTrue(os.path.exists(res))

    def test_find_attc_no_infernal(self):
        cmsearch_bin = 'foo'
        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, 'Replicons', replicon_name + '.fst')
        with self.assertRaises(RuntimeError) as ctx:
            attc.find_attc(replicon_path, self.replicon_name, cmsearch_bin, self.tmp_dir, self.model_attc)
        self.assertTrue(str(ctx.exception).endswith('failed : [Errno 2] No such file or directory'),
                        msg=str(ctx.exception))

    def test_find_attc_no_model(self):
        model_attc = 'foo'
        with self.assertRaises(RuntimeError) as ctx:
            attc.find_attc(self.replicon_path, self.replicon_name, self.cmsearch_path, self.tmp_dir, model_attc)
        self.assertTrue(str(ctx.exception).endswith('failed returncode = 1'))
