# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright (c) 2015 - 2021  Institut Pasteur, Paris and CNRS.                     #
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

import os
import shutil
import tempfile
import argparse
import unittest
import pkg_resources

import pandas as pd
import pandas.testing as pdt

# from Bio import BiopythonExperimentalWarning
# import warnings
# warnings.simplefilter('ignore', BiopythonExperimentalWarning)
from Bio import SeqIO

try:
    from tests import IntegronTest
    from tests import hide_executable
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import integrase
from integron_finder import config
from integron_finder import IntegronError
from integron_finder.scripts.finder import main
import integron_finder.scripts.finder as finder

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
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
        os.makedirs(self.out_dir)
        integrase.call = self.mute_call(_prodigal_call)
        self.find_executable_ori = finder.distutils.spawn.find_executable
        self._prefix_data = pkg_resources.resource_filename('integron_finder', "data")
        self.func_annot_dir = os.path.join(self._prefix_data, "Functional_annotation")

    def tearDown(self):
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
        integrase.call = _prodigal_call
        finder.distutils.spawn.find_executable = self.find_executable_ori


    def test_acba_simple_linear(self):
        replicon_filename = 'acba.007.p01.13'
        output_filename = 'Results_Integron_Finder_{}'.format(replicon_filename)
        test_result_dir = os.path.join(self.out_dir, output_filename)
        command = "integron_finder --outdir {out_dir} --linear {replicon}".format(out_dir=self.out_dir,
                                                                                  replicon=self.find_data(
                                                                                      os.path.join('Replicons',
                                                                                         replicon_filename + '.fst'
                                                                                         )
                                                                                       )
                                                                                  )

        with self.catch_io(out=True, err=True):
            main(command.split()[1:], loglevel='WARNING')

        output_filename = '{}.integrons'.format(replicon_filename)
        expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_acba.007.p01.13.linear',
                                                           output_filename))
        test_result_path = os.path.join(test_result_dir, output_filename)
        self.assertIntegronResultEqual(expected_result_path, test_result_path)

        summary_file_name = '{}.summary'.format(replicon_filename)
        exp_summary_path = self.find_data(
            os.path.join('Results_Integron_Finder_acba.007.p01.13.linear', summary_file_name))
        exp_summary = pd.read_csv(exp_summary_path, sep="\t", comment="#")
        test_summary_path = os.path.join(test_result_dir, summary_file_name)
        test_summary = pd.read_csv(test_summary_path, sep="\t", comment="#")
        pdt.assert_frame_equal(exp_summary, test_summary)


    def test_acba_custom_parser(self):
        replicon_filename = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        output_filename = 'Results_Integron_Finder_{}'.format(replicon_filename)
        test_result_dir = os.path.join(self.out_dir, output_filename)
        command = f"integron_finder --outdir {self.out_dir} --linear " \
                  f"--prot-file {self.find_data('Proteins', prot_name)} " \
                  f"--annot-parser {self.find_data('prodigal_annot_parser.py')} " \
                  f"{self.find_data('Replicons', replicon_filename + '.fst')}"

        with self.catch_io(out=True, err=True):
            main(command.split()[1:], loglevel='WARNING')

        output_filename = '{}.integrons'.format(replicon_filename)
        expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_acba.007.p01.13.linear',
                                                           output_filename))
        test_result_path = os.path.join(test_result_dir, output_filename)
        self.assertIntegronResultEqual(expected_result_path, test_result_path)

        summary_file_name = '{}.summary'.format(replicon_filename)
        exp_summary_path = self.find_data(
            os.path.join('Results_Integron_Finder_acba.007.p01.13.linear', summary_file_name))
        exp_summary = pd.read_csv(exp_summary_path, sep="\t", comment="#")
        test_summary_path = os.path.join(test_result_dir, summary_file_name)
        test_summary = pd.read_csv(test_summary_path, sep="\t", comment="#")
        pdt.assert_frame_equal(exp_summary, test_summary)

    def test_acba_custom_parser_no_parser(self):
        replicon_filename = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        output_filename = 'Results_Integron_Finder_{}'.format(replicon_filename)
        test_result_dir = os.path.join(self.out_dir, output_filename)
        command = f"integron_finder --outdir {self.out_dir} --linear " \
                  f"--prot-file {self.find_data('Proteins', prot_name)} " \
                  f"{self.find_data('Replicons', replicon_filename + '.fst')}"

        with self.catch_io(out=True, err=True):
            with self.assertRaises(IntegronError) as ctx:
                main(command.split()[1:], loglevel='WARNING')
        self.assertEqual(str(ctx.exception),
                         "If you provide your own proteins file for annotation (--prot-file) "
                         "you have to provide also the parser (--annot-parser)"
                         )

    def test_acba_custom_parser_no_protfile(self):
        replicon_filename = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        output_filename = 'Results_Integron_Finder_{}'.format(replicon_filename)
        test_result_dir = os.path.join(self.out_dir, output_filename)
        command = f"integron_finder --outdir {self.out_dir} --linear " \
                  f"--annot-parser {self.find_data('prodigal_annot_parser.py')} " \
                  f"{self.find_data('Replicons', replicon_filename + '.fst')}"

        with self.catch_io(out=True, err=True):
            with self.assertRaises(IntegronError) as ctx:
                main(command.split()[1:], loglevel='WARNING')
        self.assertEqual(str(ctx.exception),
                         "If you provide your own proteins file for annotation (--prot-file) "
                         "you have to provide also the parser (--annot-parser)"
                         )



    def test_acba_sequential_eq_isolate(self):
        """
        test if we find the same results if we run IF in sequential as isolated on each seq
        ACBA.0917.00019 contains 2 contigs 0001 and 0002.
        0002 does not contains integrons

        .integrons file should be identical
        .summary file should contains the 001 contig and 0 calin, 0 complete and 0 in0
        :return:
        """
        seq_replicon_filename = 'ACBA.0917.00019'
        seq_output_dir = f'Results_Integron_Finder_{seq_replicon_filename}'
        seq_test_result_dir = os.path.join(self.out_dir, seq_output_dir)
        seq_cmd = "integron_finder --outdir {out_dir} " \
                  "--keep-tmp {replicon}".format(out_dir=self.out_dir,
                                                 replicon=self.find_data(
                                                          os.path.join('Gembase',
                                                                       'Replicons',
                                                                       seq_replicon_filename + '.fna')
                                                          )
                                                 )
        with self.catch_io(out=True, err=True):
            main(seq_cmd.split()[1:], loglevel='WARNING')

        contig_name = 'ACBA.0917.00019.0001'
        iso_output_dir = 'Results_Integron_Finder_{}'.format(contig_name)
        iso_test_result_dir = os.path.join(self.out_dir, iso_output_dir)
        iso_cmd = "integron_finder --outdir {out_dir} " \
                  "--keep-tmp --lin {replicon}".format(out_dir=self.out_dir,
                                                 replicon=self.find_data(
                                                          os.path.join('Replicons',
                                                                        contig_name + '.fst')
                                                          )
                                                 )
        with self.catch_io(out=True, err=True):
            main(iso_cmd.split()[1:], loglevel='WARNING')

        seq_integron_result_path = os.path.join(seq_test_result_dir, seq_replicon_filename + '.integrons')
        iso_integron_result_path = os.path.join(iso_test_result_dir, contig_name + '.integrons')
        self.assertIntegronResultEqual(seq_integron_result_path, iso_integron_result_path)

        seq_summary_result_path = os.path.join(seq_test_result_dir, seq_replicon_filename + '.summary')
        iso_summary_result_path = os.path.join(iso_test_result_dir, contig_name + '.summary')
        seq_summary = pd.read_csv(seq_summary_result_path, sep="\t", comment='#')
        iso_summary = pd.read_csv(iso_summary_result_path, sep="\t", comment='#')
        summary_first_contig = seq_summary.loc[seq_summary['ID_replicon'] == 'ACBA.0917.00019.0001']
        pdt.assert_frame_equal(summary_first_contig, iso_summary)

        summary_2nd_contig = seq_summary.loc[seq_summary['ID_replicon'] == 'ACBA.0917.00019.0002']
        # the index are different 0/1 so I fix this by indexing by 'ID_replicon'
        summary_2nd_contig.set_index(['ID_replicon'], inplace=True)
        iso_summary = pd.DataFrame([['ACBA.0917.00019.0002', 0, 0, 0, 'lin', 8729]],
                                   columns=['ID_replicon', 'CALIN', 'complete', 'In0', 'topology', 'size'])
        iso_summary.set_index(['ID_replicon'], inplace=True)
        pdt.assert_frame_equal(summary_2nd_contig, iso_summary)


    def test_acba_simple_gembase(self):
        """
        ACBA.0917.00019 contains 2 contigs 0001 and 0002.
        0002 does not contains integrons
        """
        replicon_filename = 'ACBA.0917.00019'
        contig_id = 'ACBA.0917.00019.0001'
        output_filename = 'Results_Integron_Finder_{}'.format(replicon_filename)
        test_result_dir = os.path.join(self.out_dir, output_filename)
        command = "integron_finder --outdir {out_dir} --keep-tmp " \
                  "--promoter-attI --gembase {replicon}".format(out_dir=self.out_dir,
                                                                replicon=self.find_data(
                                                                   os.path.join('Gembase',
                                                                                'Replicons',
                                                                                replicon_filename + '.fna')
                                                                )
                                                                )
        with self.catch_io(out=True, err=True):
            main(command.split()[1:], loglevel='WARNING')

        output_filename = '{}.integrons'.format(replicon_filename)
        expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_{}.gembase'.format(replicon_filename),
                                                           output_filename))
        test_result_path = os.path.join(test_result_dir, output_filename)
        self.assertIntegronResultEqual(expected_result_path, test_result_path)

        summary_file_name = '{}.summary'.format(replicon_filename)
        exp_summary_path = self.find_data(
            os.path.join('Results_Integron_Finder_{}.gembase'.format(replicon_filename), summary_file_name))
        exp_summary = pd.read_csv(exp_summary_path, sep="\t", comment='#')
        test_summary_path = os.path.join(test_result_dir, summary_file_name)
        test_summary = pd.read_csv(test_summary_path, sep="\t", comment='#')
        pdt.assert_frame_equal(exp_summary, test_summary)


    def test_acba_simple_no_gbk_no_pdf(self):
        replicon_filename = 'acba.007.p01.13'
        replicon_id = 'ACBA.007.P01_13'
        output_filename = 'Results_Integron_Finder_{}'.format(replicon_filename)
        test_result_dir = os.path.join(self.out_dir, output_filename)
        command = "integron_finder --outdir {out_dir} --promoter-attI {replicon}".format(out_dir=self.out_dir,
                                                                         replicon=self.find_data(
                                                                             os.path.join('Replicons',
                                                                                          '{}.fst'.format(replicon_filename))
                                                                         )
                                                                         )

        with self.catch_io(out=True, err=True):
            main(command.split()[1:], loglevel='WARNING')

        output_filename = '{}.integrons'.format(replicon_filename)
        expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_acba.007.p01.13',
                                                           output_filename))
        test_result_path = os.path.join(test_result_dir, output_filename)
        self.assertIntegronResultEqual(expected_result_path, test_result_path)

        summary_file_name = '{}.summary'.format(replicon_filename)
        exp_summary_path = self.find_data(os.path.join('Results_Integron_Finder_acba.007.p01.13', summary_file_name))
        exp_summary = pd.read_csv(exp_summary_path, sep="\t", comment='#')
        test_summary_path = os.path.join(test_result_dir, summary_file_name)
        test_summary = pd.read_csv(test_summary_path, sep="\t", comment='#')
        pdt.assert_frame_equal(exp_summary, test_summary)

        gbk = '{}.gbk'.format(replicon_filename)
        gbk_test = os.path.join(test_result_dir, gbk)
        self.assertFalse(os.path.exists(gbk_test))

        pdf = '{}_1.pdf'.format(replicon_filename)
        pdf_test = os.path.join(test_result_dir, pdf)
        self.assertFalse(os.path.exists(pdf_test))


    def test_acba_simple_with_gbk(self):
        replicon_filename = 'acba.007.p01.13'
        replicon_id = 'ACBA.007.P01_13'
        output_dirname = 'Results_Integron_Finder_{}'.format(replicon_filename)
        test_result_dir = os.path.join(self.out_dir, output_dirname)
        command = "integron_finder --outdir {out_dir} --gbk --promoter-attI {replicon}".format(out_dir=self.out_dir,
                                                                               replicon=self.find_data(
                                                                                   os.path.join('Replicons',
                                                                                                '{}.fst'.format(replicon_filename))
                                                                                )
                                                                               )

        with self.catch_io(out=True, err=True):
            main(command.split()[1:], loglevel='WARNING')

        gbk = '{}.gbk'.format(replicon_id)
        expected_gbk = self.find_data(os.path.join('Results_Integron_Finder_{}'.format(replicon_filename), gbk))
        gbk_test = os.path.join(test_result_dir, gbk)
        expected_gbk = SeqIO.read(expected_gbk, 'gb')
        gbk_test = SeqIO.read(gbk_test, 'gb')
        self.assertSeqRecordEqual(expected_gbk, gbk_test)

        output_filename = '{}.integrons'.format(replicon_filename)
        expected_result_path = self.find_data(os.path.join(output_dirname, output_filename))
        test_result_path = os.path.join(test_result_dir, output_filename)
        self.assertIntegronResultEqual(expected_result_path, test_result_path)

    def test_acba_simple_with_gbk_without_promoter(self):
        replicon_filename = 'acba.007.p01.13'
        replicon_id = 'ACBA.007.P01_13'
        command = "integron_finder --outdir {out_dir} --gbk {replicon}".format(out_dir=self.out_dir,
                                                                               replicon=self.find_data(
                                                                                   os.path.join('Replicons',
                                                                                                '{}.fst'.format(replicon_filename))
                                                                                )
                                                                               )

        with self.catch_io(out=True, err=True):
            main(command.split()[1:], loglevel='WARNING')

        output_dirname = 'Results_Integron_Finder_{}'.format(replicon_filename)
        test_result_dir = os.path.join(self.out_dir, output_dirname)
        gbk = '{}.gbk'.format(replicon_id)
        expected_gbk = self.find_data(os.path.join(output_dirname + ".wo_promoter", gbk))
        gbk_test = os.path.join(test_result_dir, gbk)
        expected_gbk = SeqIO.read(expected_gbk, 'gb')
        gbk_test = SeqIO.read(gbk_test, 'gb')
        self.assertSeqRecordEqual(expected_gbk, gbk_test)

        output_filename = '{}.integrons'.format(replicon_filename)
        expected_result_path = self.find_data(os.path.join(output_dirname + ".wo_promoter", output_filename))
        test_result_path = os.path.join(test_result_dir, output_filename)
        self.assertIntegronResultEqual(expected_result_path, test_result_path)


    def test_acba_simple_with_pdf(self):
        replicon_filename = 'acba.007.p01.13'
        replicon_id = 'ACBA.007.P01_13'
        output_dirname = 'Results_Integron_Finder_{}'.format(replicon_filename)
        test_result_dir = os.path.join(self.out_dir, output_dirname)
        command = "integron_finder --outdir {out_dir} --pdf " \
                  "--promoter-attI {replicon}".format(out_dir=self.out_dir,
                                                      replicon=self.find_data(
                                                          os.path.join('Replicons', '{}.fst'.format(replicon_filename))
                                                        )
                                                      )

        with self.catch_io(out=True, err=True):
            main(command.split()[1:], loglevel='WARNING')
        pdf = '{}_1.pdf'.format(replicon_id)
        pdf_test = os.path.join(test_result_dir, pdf)
        self.assertTrue(os.path.exists(pdf_test))


    def test_acba_annot(self):
        replicon_filename = 'acba.007.p01.13'
        replicon_id = 'ACBA.007.P01_13'
        command = "integron_finder --outdir {out_dir} --func-annot --path-func-annot {annot_bank} --promoter-attI " \
                  "--gbk --keep-tmp " \
                  "{replicon}".format(out_dir=self.out_dir,
                                      annot_bank=self.func_annot_dir,
                                      replicon=self.find_data(os.path.join('Replicons', '{}.fst'.format(replicon_filename)))
                                      )

        with self.catch_io(out=True, err=False):
            main(command.split()[1:], loglevel='WARNING')

        result_dir = os.path.join(self.out_dir, 'Results_Integron_Finder_{}'.format(replicon_filename))

        gbk = '{}.gbk'.format(replicon_id)
        expected_gbk = self.find_data(os.path.join('Results_Integron_Finder_{}.annot'.format(replicon_filename), gbk))
        gbk_test = os.path.join(result_dir, gbk)
        expected_gbk = SeqIO.read(expected_gbk, 'gb')
        gbk_test = SeqIO.read(gbk_test, 'gb')
        self.assertSeqRecordEqual(expected_gbk, gbk_test)

        output_filename = '{}.integrons'.format(replicon_filename)
        expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_{}.annot'.format(replicon_filename),
                                                           output_filename))
        test_result_path = os.path.join(result_dir, output_filename)
        self.assertIntegronResultEqual(expected_result_path, test_result_path)

        output_filename = os.path.join('tmp_{}'.format(replicon_id), replicon_id + '_NCBIfam-AMRFinder_fa_table.res')
        expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_{}.annot'.format(replicon_filename),
                                                           output_filename))
        test_result_path = os.path.join(result_dir, output_filename)
        self.assertHmmEqual(expected_result_path, test_result_path)


    def test_acba_local_max(self):
        replicon_filename = 'acba.007.p01.13'
        replicon_id = 'ACBA.007.P01_13'
        command = "integron_finder --outdir {out_dir} --func-annot --path-func-annot {annot_bank} --local-max --gbk " \
                  "--keep-tmp --promoter-attI {replicon}".format(
                                    out_dir=self.out_dir,
                                    annot_bank=self.func_annot_dir,
                                    replicon=self.find_data(os.path.join('Replicons', '{}.fst'.format(replicon_filename)))
                                )
        with self.catch_io(out=True, err=True):
            main(command.split()[1:], loglevel='WARNING')

        result_dir = os.path.join(self.out_dir, 'Results_Integron_Finder_{}'.format(replicon_filename))

        gbk = '{}.gbk'.format(replicon_id)
        expected_gbk = self.find_data(os.path.join('Results_Integron_Finder_{}.local_max'.format(replicon_filename), gbk))
        gbk_test = os.path.join(result_dir, gbk)
        expected_gbk = SeqIO.read(expected_gbk, 'gb')
        gbk_test = SeqIO.read(gbk_test, 'gb')
        self.assertSeqRecordEqual(expected_gbk, gbk_test)

        output_filename = '{}.integrons'.format(replicon_filename)
        expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_{}.local_max'.format(replicon_filename),
                                                           output_filename))
        test_result_path = os.path.join(result_dir, output_filename)
        self.assertIntegronResultEqual(expected_result_path, test_result_path)

        output_filename = os.path.join('tmp_{}'.format(replicon_id), '{}_NCBIfam-AMRFinder_fa_table.res'.format(replicon_id))
        expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_{}.local_max'.format(replicon_filename),
                                                           output_filename))

        test_result_path = os.path.join(result_dir, output_filename)
        self.assertHmmEqual(expected_result_path, test_result_path)

        output_filename = os.path.join('tmp_{}'.format(replicon_id), '{}_13825_1014_subseq_attc_table.res'.format(replicon_id))
        expected_result_path = self.find_data(os.path.join('Results_Integron_Finder_{}.local_max'.format(replicon_filename),
                                                           output_filename))
        test_result_path = os.path.join(result_dir, output_filename)
        with open(expected_result_path) as expected_result_file, open(test_result_path) as test_result_file:
            for expected_line, result_line in zip(expected_result_file, test_result_file):
                if result_line.startswith('# Program: '):
                    break
                self.assertEqual(expected_line, result_line)


    def test_no_integron(self):
        replicon_filename = 'fake_seq'
        replicon_id = 'fake_seq'
        output_filename = 'Results_Integron_Finder_{}'.format(replicon_filename)
        test_result_dir = os.path.join(self.out_dir, output_filename)
        command = "integron_finder --outdir {out_dir} {replicon}".format(out_dir=self.out_dir,
                                                                         replicon=self.find_data(
                                                                             os.path.join('Replicons',
                                                                                          '{}.fst'.format(replicon_filename))
                                                                         )
                                                                         )
        with self.catch_io(out=True, err=True):
            main(command.split()[1:], loglevel='WARNING')

        output_filename = 'fake_seq.integrons'
        test_result_path = os.path.join(test_result_dir, output_filename)
        with open(test_result_path) as tested_file:
            test_line = next(tested_file)
            self.assertTrue(test_line.startswith('# cmd: integron_finder'))
            test_line = next(tested_file)
            self.assertEqual(test_line.strip(), '# No Integron found')


    def test_acba_no_hmmer(self):
        replicon_filename = 'acba.007.p01.13'
        decorator = hide_executable('hmmsearch')
        finder.distutils.spawn.find_executable = decorator(finder.distutils.spawn.find_executable)
        command = "integron_finder --outdir {out_dir} {replicon}".format(out_dir=self.out_dir,
                                                                         replicon=self.find_data(
                                                                             os.path.join('Replicons',
                                                                                          '{}.fst'.format(replicon_filename))
                                                                         )
                                                                         )
        with self.assertRaises(RuntimeError) as ctx:
            with self.catch_io(out=True):
                # in case the error is not raised
                # anyway do not want to mess up the test output
                # I cannot catch log because loggers are reinitialized in main
                # I need to catch stdout as log are write on
                main(command.split()[1:])

        err_msg = "cannot find 'hmmsearch' in PATH.\n" \
                  "Please install hmmer package or setup 'hmmsearch' binary path with --hmmsearch option"
        self.assertEqual(err_msg, str(ctx.exception))


    def test_acba_no_prodigal(self):
        replicon_filename = 'acba.007.p01.13'
        decorator = hide_executable('prodigal')
        finder.distutils.spawn.find_executable = decorator(finder.distutils.spawn.find_executable)
        command = "integron_finder --outdir {out_dir} {replicon}".format(out_dir=self.out_dir,
                                                                         replicon=self.find_data(
                                                                             os.path.join('Replicons',
                                                                                          '{}.fst'.format(replicon_filename))
                                                                         )
                                                                         )
        with self.assertRaises(RuntimeError) as ctx:
            with self.catch_io(out=True):
                # in case the error is not raised
                # anyway do not want to mess up the test output
                # I cannot catch log because loggers are reinitialized in main
                # I need to catch stdout as log are write on
                main(command.split()[1:])

        err_msg = "cannot find 'prodigal' in PATH.\n" \
                  "Please install prodigal package or setup 'prodigal' binary path with --prodigal option"
        self.assertEqual(err_msg, str(ctx.exception))


    def test_acba_no_cmsearch(self):
        replicon_filename = 'acba.007.p01.13'
        decorator = hide_executable('cmsearch')
        finder.distutils.spawn.find_executable = decorator(finder.distutils.spawn.find_executable)
        command = "integron_finder --outdir {out_dir} {replicon}".format(out_dir=self.out_dir,
                                                                         replicon=self.find_data(
                                                                             os.path.join('Replicons',
                                                                                          '{}.fst'.format(replicon_filename))
                                                                         )
                                                                         )
        with self.assertRaises(RuntimeError) as ctx:
            with self.catch_io(out=True):
                # in case the error is not raised
                # anyway do not want to mess up the test output
                # I cannot catch log because loggers are reinitialized in main
                # I need to catch stdout as log are write on
                main(command.split()[1:])
        err_msg = "cannot find 'cmsearch' in PATH.\n" \
                  "Please install infernal package or setup 'cmsearch' binary path with --cmsearch option"
        self.assertEqual(err_msg, str(ctx.exception))


    def test_outdir_is_file(self):
        replicon_filename = 'acba.007.p01.13'
        bad_out_dir = os.path.join(self.out_dir, 'bad_out_dir')
        open(bad_out_dir, 'w').close()
        command = "integron_finder --outdir {out_dir} {replicon}".format(out_dir=bad_out_dir,
                                                                         replicon=self.find_data(
                                                                             os.path.join('Replicons',
                                                                                          '{}.fst'.format(replicon_filename))
                                                                         )
                                                                         )
        with self.assertRaises(IsADirectoryError) as ctx:
            with self.catch_io(out=True):
                # in case the error is not raised
                # anyway do not want to mess up the test output
                # I cannot catch log because loggers are reinitialized in main
                # I need to catch stdout as log are write on
                main(command.split()[1:])
        err_msg = "outdir '{}' already exists and is not a directory".format(bad_out_dir)
        self.assertEqual(err_msg, str(ctx.exception))


    def test_resultdir_is_file(self):
        replicon_filename = 'acba.007.p01.13'
        args = argparse.Namespace()
        args.replicon = self.find_data(os.path.join('Replicons', '{}.fst'.format(replicon_filename)))
        args.outdir = self.out_dir
        cf = config.Config(args)
        result_dir_is_file = cf.result_dir
        open(result_dir_is_file, 'w').close()
        command = "integron_finder --outdir {out_dir} {replicon}".format(out_dir=self.out_dir,
                                                                         replicon=self.find_data(
                                                                             os.path.join('Replicons',
                                                                                          '{}.fst'.format(replicon_filename))
                                                                         )
                                                                         )
        with self.assertRaises(IsADirectoryError) as ctx:
            with self.catch_io(out=True):
                # in case the error is not raised
                # anyway do not want to mess up the test output
                # I cannot catch log because loggers are reinitialized in main
                # I need to catch stdout as log are write on
                main(command.split()[1:])
        err_msg = "result dir '{}' already exists and is not a directory".format(self.out_dir)
        self.assertEqual(err_msg, str(ctx.exception))


    @unittest.skipIf(os.getuid() == 0, "root have always permission to write")
    def test_resultdir_not_writable(self):
        replicon_filename = 'acba.007.p01.13'
        args = argparse.Namespace()
        args.replicon = self.find_data(os.path.join('Replicons', '{}.fst'.format(replicon_filename)))
        args.outdir = self.out_dir
        cf = config.Config(args)
        os.mkdir(cf.result_dir, mode=0o500)
        command = "integron_finder --outdir {out_dir} {replicon}".format(out_dir=self.out_dir,
                                                                         replicon=self.find_data(
                                                                             os.path.join('Replicons',
                                                                                          '{}.fst'.format(replicon_filename))
                                                                         )
                                                                         )
        with self.assertRaises(PermissionError) as ctx:
            with self.catch_io(out=True):
                # in case the error is not raised
                # anyway do not want to mess up the test output
                # I cannot catch log because loggers are reinitialized in main
                # I need to catch stdout as log are write on
                main(command.split()[1:])
        err_msg = "result dir '{}' already exists and is not writable".format(self.out_dir)
        self.assertEqual(err_msg, str(ctx.exception))

