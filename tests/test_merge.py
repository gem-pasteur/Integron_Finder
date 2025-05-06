# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright (c) 2015 - 2025  Institut Pasteur, Paris and CNRS.                     #
# See the COPYRIGHT file for details                                               #
#                                                                                  #
# integron_finder is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU General Public License as published by             #
# the Free Software Foundation, either version 3 of the License, or                #
# (at your option) any later version.                                              #
#                                                                                  #
# integron_finder is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                   #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    #
# GNU General Public License for more details.                                     #
#                                                                                  #
# You should have received a copy of the GNU General Public License                #
# along with this program (COPYING file).                                          #
# If not, see <http://www.gnu.org/licenses/>.                                      #
####################################################################################

import tempfile
import os
import shutil

import pandas as pd
import pandas.testing as pdt

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import IntegronError, logger_set_level
import integron_finder.scripts.merge as merge


class TestMerge(IntegronTest):

    def setUp(self):
        tmp_dir = tempfile.gettempdir()
        self.out_dir = os.path.join(tmp_dir, 'test_integron_merge')
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
        os.makedirs(self.out_dir)
        self.replicons = ('acba.007.p01.13', 'lian.001.c02.10', 'pssu.001.c01.13')
        self.ids = ('ACBA.007.P01_13', 'LIAN.001.C02_10', 'PSSU.001.C01_13')
        self.res_dirs = []
        # copy *.intgrons in tmp
        for rep in self.replicons:
            res_dir = os.path.join(self.out_dir, f'Result_{rep}')
            os.makedirs(res_dir)
            self.res_dirs.append(res_dir)
            res_file = self.find_data(f"{rep}_local_max_lin.integrons")
            shutil.copyfile(res_file, os.path.join(res_dir, f"{rep}.integrons"))
        # copy *.summary *.gbk in tmp
        for _id in self.ids:
            for ext in ('gbk', 'summary'):
                src_file = self.find_data("{}_local_max_lin.{}".format(_id, ext))
                shutil.copyfile(src_file, os.path.join(res_dir, "{}.{}".format(_id, ext)))

            pdf_file = self.find_data("{}_1_local_max_lin.pdf".format(_id))
            shutil.copyfile(pdf_file, os.path.join(res_dir, "{}_1.pdf".format(_id)))


    def tearDown(self):
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
            pass
        logger_set_level('INFO')

    def test_merge_integrons(self):
        outfile = os.path.join(self.out_dir, 'merged.integrons')
        merge.merge_integrons(outfile, *self.res_dirs)
        agg_results = pd.read_csv(outfile, sep="\t", comment="#")

        expe_result_file = self.find_data('acba_lian_pssu_merged.integrons')
        expected_results = pd.read_csv(expe_result_file, sep="\t", comment="#")
        pdt.assert_frame_equal(agg_results, expected_results)

    def test_merge_summary(self):
        outfile = os.path.join(self.out_dir, 'merged.summary')
        merge.merge_summary(outfile, *self.res_dirs)
        agg_results = pd.read_csv(outfile, sep="\t", comment="#")
        agg_results.sort_values(by=['ID_replicon'], inplace=True)
        agg_results.reset_index(inplace=True, drop=True)

        expe_result_file = self.find_data('acba_lian_pssu_merged.summary')
        expected_results = pd.read_csv(expe_result_file, sep="\t", comment="#")
        expected_results.sort_values(by=['ID_replicon'], inplace=True)
        expected_results.reset_index(inplace=True, drop=True)
        pdt.assert_frame_equal(agg_results, expected_results)


    def test_merge_no_files(self):
        outfile = os.path.join(self.out_dir, 'merged.integrons')
        for rep, res_dir in zip(self.replicons, self.res_dirs):
            os.unlink(os.path.join(res_dir, "{}.integrons".format(rep)))
        exp_msg = 'No integrons file to merge'
        with self.assertRaises(IntegronError) as ctx:
            # 100 is to disable CRITICAL log message
            logger_set_level(100)
            merge.merge_integrons(outfile, *self.res_dirs)
        self.assertEqual(str(ctx.exception), exp_msg)

    def test_copy_file(self):
        merge.copy_file(self.out_dir, '.gbk', *self.res_dirs)
        merge.copy_file(self.out_dir, '.pdf', *self.res_dirs)
        for _id in self.ids:
            self.assertTrue(os.path.exists(os.path.join(self.out_dir, '{}.gbk'.format(_id))))
            self.assertTrue(os.path.exists(os.path.join(self.out_dir, '{}_1.pdf'.format(_id))))

class TestParseArgs(IntegronTest):

    def test_parse_one_result(self):
        parsed_args = merge.parse_args(['outdir', 'outfile', 'result_1'])
        self.assertEqual(parsed_args.outdir, 'outdir')
        self.assertEqual(parsed_args.outfile, 'outfile')
        self.assertListEqual(parsed_args.results, ['result_1'])
        self.assertEqual(parsed_args.quiet, 0)
        self.assertEqual(parsed_args.verbose, 0)

    def test_parse_2_results(self):
        parsed_args = merge.parse_args(['outdir', 'outfile', 'result_1', 'result_2'])
        self.assertEqual(parsed_args.outdir, 'outdir')
        self.assertEqual(parsed_args.outfile, 'outfile')
        self.assertListEqual(parsed_args.results, ['result_1', 'result_2'])
        self.assertEqual(parsed_args.quiet, 0)
        self.assertEqual(parsed_args.verbose, 0)

    def test_verbose(self):
        parsed_args = merge.parse_args(['outdir', 'outfile', 'result_1', 'result_2'])
        self.assertEqual(parsed_args.verbose, 0)
        parsed_args = merge.parse_args(['--verbose', 'outdir', 'outfile', 'result'])
        self.assertEqual(parsed_args.verbose, 1)
        parsed_args = merge.parse_args(['-vv', 'outdir', 'outfile', 'result'])
        self.assertEqual(parsed_args.verbose, 2)

    def test_quiet(self):
        parsed_args = merge.parse_args(['outdir', 'outfile', 'result'])
        self.assertEqual(parsed_args.quiet, 0)
        parsed_args = merge.parse_args(['--quiet', 'outdir', 'outfile', 'result'])
        self.assertEqual(parsed_args.quiet, 1)
        parsed_args = merge.parse_args(['-qq', 'outdir', 'outfile', 'result'])
        self.assertEqual(parsed_args.quiet, 2)

    def test_outdir_results_are_equals(self):
        with self.assertRaises(ValueError) as ctx:
            _ = merge.parse_args(['outdir', 'outfile', 'outdir'])
        self.assertEqual(str(ctx.exception),
                         "'outdir' and 'results' cannot have the same value.")
        
class TestMain(IntegronTest):

    def setUp(self):
        tmp_dir = tempfile.gettempdir()
        self.out_dir = os.path.join(tmp_dir, 'test_integron_merge')
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
        os.makedirs(self.out_dir)
        self.replicons = ('acba.007.p01.13', 'lian.001.c02.10', 'pssu.001.c01.13')
        self.ids = ('ACBA.007.P01_13', 'LIAN.001.C02_10', 'PSSU.001.C01_13')
        self.res_dirs = []
        for rep in self.replicons:
            res_dir = os.path.join(self.out_dir, 'Result_{}'.format(rep))
            os.makedirs(res_dir)
            self.res_dirs.append(res_dir)

            res_file = self.find_data("{}_local_max_lin.integrons".format(rep))
            shutil.copyfile(res_file, os.path.join(res_dir, "{}.integrons".format(rep)))

        for _id in self.ids:
            gbk_file = self.find_data("{}_local_max_lin.gbk".format(_id))
            shutil.copyfile(gbk_file, os.path.join(res_dir, "{}.gbk".format(_id)))
            pdf_file = self.find_data("{}_1_local_max_lin.pdf".format(_id))
            shutil.copyfile(pdf_file, os.path.join(res_dir, "{}_1.pdf".format(_id)))


    def tearDown(self):
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)


    def test_main(self):
        outfile = 'merged.integrons'
        command = 'integron_merge {} {} {}'.format(self.out_dir, outfile, ' '.join(self.res_dirs))
        with self.catch_io(out=True, err=True):
            merge.main(command.split()[1:])
            # out = sys.stdout.getvalue()

        integron_file = os.path.join(self.out_dir, outfile)
        merge.merge_integrons(integron_file, *self.res_dirs)
        agg_results = pd.read_csv(integron_file, sep="\t", comment="#")

        expe_result_file = self.find_data('acba_lian_pssu_merged.integrons')
        expected_results = pd.read_csv(expe_result_file, sep="\t", comment="#")

        pdt.assert_frame_equal(agg_results, expected_results)

        merge.copy_file(self.out_dir, '.gbk', *self.res_dirs)
        merge.copy_file(self.out_dir, '.pdf', *self.res_dirs)
        for _id in self.ids:
            self.assertTrue(os.path.exists(os.path.join(self.out_dir, '{}.gbk'.format(_id))))
            self.assertTrue(os.path.exists(os.path.join(self.out_dir, '{}_1.pdf'.format(_id))))

    def test_main_outdir_exists(self):
        bad_outdir = os.path.join(self.out_dir, 'bad_out_dir')
        open(bad_outdir, 'w').close()
        exp_msg = "'{}' already exists and is not a directory".format(bad_outdir)
        with self.assertRaises(IOError) as ctx:
            command = 'integron_merge {} {} {}'.format(bad_outdir, 'whatever', ' '.join(self.res_dirs))
            # 100 is to disable CRITICAL log message
            merge.main(command.split()[1:], log_level=100)
        self.assertEqual(str(ctx.exception), exp_msg)

    def test_main_no_outdir(self):
        outfile = 'merged.integrons'
        out_dir = os.path.join(self.out_dir, 'real_outdir')
        command = 'integron_merge {} {} {}'.format(out_dir, outfile, ' '.join(self.res_dirs))
        with self.catch_io(out=True, err=True):
            merge.main(command.split()[1:])
            # out = sys.stdout.getvalue()

        integron_file = os.path.join(out_dir, outfile)
        merge.merge_integrons(integron_file, *self.res_dirs)
        agg_results = pd.read_csv(integron_file, sep="\t", comment="#")

        expe_result_file = self.find_data('acba_lian_pssu_merged.integrons')
        expected_results = pd.read_csv(expe_result_file, sep="\t", comment="#")

        pdt.assert_frame_equal(agg_results, expected_results)

        merge.copy_file(self.out_dir, '.gbk', *self.res_dirs)
        merge.copy_file(self.out_dir, '.pdf', *self.res_dirs)
        for _id in self.ids:
            self.assertTrue(os.path.exists(os.path.join(out_dir, '{}.gbk'.format(_id))))
            self.assertTrue(os.path.exists(os.path.join(out_dir, '{}_1.pdf'.format(_id))))
