# -*- coding: utf-8 -*-

from __future__ import print_function
import shutil
import tempfile
import os
import sys
from subprocess import Popen, PIPE
from itertools import izip

from tests import IntegronTest, which


class Test(IntegronTest):


    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..' '..')))
        self.tmp_dir = tempfile.gettempdir()
        self.bin = os.path.join(self.integron_home, 'integron_finder') if self.local_install else which('integron_finder')
        self.replicon_name = 'acba.007.p01.13'


    def tearDown(self):
        try:
            shutil.rmtree(self.out_dir)
            pass
        except:
            pass


    def test_acba_annot(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'integron_acba_annot_test')
        os.makedirs(self.out_dir)
        output_filename = 'Results_Integron_Finder_' + self.replicon_name
        test_result_dir = os.path.join(self.out_dir, output_filename)
        command = "{python_exe} {bin} --outdir {out_dir} --func_annot --path_func_annot {annot_bank} {replicon}".format(
            python_exe=sys.executable,
            bin=self.bin,
            out_dir=self.out_dir,
            annot_bank=os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'data', 'Functional_annot')),
            replicon=os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'data', self.replicon_name + '.fst'))
        )
        print("\n", command)
        if not self.bin:
            raise RuntimeError('coverage not found, INTEGRON_HOME must be either in your path or INTEGRON_HOME'
                               ' must be defined command launch: \n{}'.format(command))

        try:
            integron_process = Popen(command,
                                     shell=True,
                                     stdin=None,
                                     stderr=PIPE,
                                     close_fds=False
                                     )
        except Exception as err:
            msg = "integron_finder execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err

        integron_process.wait()
        self.assertEqual(integron_process.returncode, 0,
                         "integron_finder finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                             integron_process.returncode,
                             command,
                             ''.join([l.decode('utf-8') for l in integron_process.stderr.readlines()]),
                                                                                                                   )
                         )
        results_file_to_test = (self.replicon_name + '.gbk',
                                self.replicon_name + '.integrons',
                                )
        for output_filename in results_file_to_test:
            expected_result_path = os.path.join(self._data_dir, 'Results_Integron_Finder_' + self.replicon_name +'.annot', output_filename)
            test_result_path = os.path.join(test_result_dir, output_filename)
            with open(expected_result_path) as expected_result_file, open(test_result_path) as test_result_file :
                # test equality line by line
                for expected_line, result_line in izip(expected_result_file, test_result_file):
                    self.assertEqual(expected_line, result_line)

        output_filename = os.path.join('other', self.replicon_name + '_Resfams_fa_table.res')
        expected_result_path = os.path.join(self._data_dir, 'Results_Integron_Finder_' + self.replicon_name + '.annot', output_filename)
        test_result_path = os.path.join(test_result_dir, output_filename)
        with open(expected_result_path) as expected_result_file, open(test_result_path) as test_result_file :
            for expected_line, result_line in izip(expected_result_file, test_result_file):
                if result_line.startswith('# Program:         hmmsearch'):
                    break
                self.assertEqual(expected_line, result_line)