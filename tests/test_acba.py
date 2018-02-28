import os
import shutil
import tempfile

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import integrase
from integron_finder.scripts.finder import main

_prodigal_call = integrase.call

class TestAcba(IntegronTest):

    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

        self.tmp_dir = tempfile.gettempdir()
        self.out_dir = os.path.join(self.tmp_dir, 'integron_acba_test')
        os.makedirs(self.out_dir)
        integrase.call = self.call_wrapper(_prodigal_call)

    def tearDown(self):
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
        integrase.call = _prodigal_call


    def test_acba_simple(self):
        output_filename = 'Results_Integron_Finder_acba.007.p01.13'
        test_result_dir = os.path.join(self.out_dir, output_filename)
        command = "integron_finder --outdir {out_dir} {replicon}".format(out_dir=self.out_dir,
                                                                         replicon=self.find_data(
                                                                             os.path.join('Replicons',
                                                                                          'acba.007.p01.13.fst')
                                                                         )
                                                                         )
        with self.catch_io(out=True, err=True):
            main(command.split()[1:])
        results_file_to_test = ('acba.007.p01.13.gbk', 'acba.007.p01.13.integrons')
        for output_filename in results_file_to_test:
            expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_acba.007.p01.13',
                                                               output_filename))
            test_result_path = os.path.join(test_result_dir, output_filename)
            self.assertFileEqual(expected_result_path, test_result_path)
