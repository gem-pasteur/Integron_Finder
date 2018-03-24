# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright © 2015 - 2018  Institut Pasteur, Paris.                                #
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
import distutils
import sys

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.scripts.finder import parse_args


class TestParseArgs(IntegronTest):

    def test_replicon(self):
        replicon = 'foo'
        cfg = parse_args([replicon])
        self.assertEqual(cfg.replicon_path, os.path.abspath(replicon))

    def test_wo_replicon(self):
        real_exit = sys.exit

        sys.exit = self.fake_exit
        with self.catch_io(err=True):
            try:
                _ = parse_args([])
            except TypeError as err:
                msg = sys.stderr.getvalue()
                msg_end = 'error: too few arguments\n'
                self.assertTrue(msg.endswith(msg_end), "{} != {}".format(msg[len(msg) - len(msg_end):], msg_end))
                # program exit with returncode = 2
                self.assertEqual(str(err), '2')
            finally:
                sys.exit = real_exit

    def test_local_max(self):
        cfg = parse_args(['replicon'])
        self.assertFalse(cfg.local_max)
        cfg = parse_args(['--local_max', 'replicon'])
        self.assertTrue(cfg.local_max)

    def test_func_annot(self):
        cfg = parse_args(['replicon'])
        self.assertFalse(cfg.func_annot)
        cfg = parse_args(['--func_annot', 'replicon'])
        self.assertTrue(cfg.func_annot)

    def test_cpu(self):
        cfg = parse_args(['replicon'])
        self.assertEqual(cfg.cpu, 1)
        cpu = 10
        cfg = parse_args(['--cpu', str(cpu), 'replicon'])
        self.assertEqual(cfg.cpu, cpu)

    def test_distance_threshold(self):
        cfg = parse_args(['replicon'])
        self.assertEqual(cfg.distance_threshold, 4000)
        dt = 50
        cfg = parse_args(['--distance_thresh', str(dt), 'replicon'])
        self.assertEqual(cfg.distance_threshold, dt)

    def test_outdir(self):
        cfg = parse_args(['replicon'])
        self.assertEqual(cfg.outdir, os.path.abspath('.'))
        outdir = 'foo'
        cfg = parse_args(['--outdir', outdir, 'replicon'])
        self.assertEqual(cfg.outdir, os.path.abspath(outdir))

    def test_union_integrase(self):
        cfg = parse_args(['replicon'])
        self.assertFalse(cfg.union_integrases)
        cfg = parse_args(['--union_integrases', 'replicon'])
        self.assertTrue(cfg.union_integrases)

    def test_cmsearch(self):
        cfg = parse_args(['replicon'])
        self.assertEqual(cfg.cmsearch, distutils.spawn.find_executable("cmsearch"))
        cmsearch = 'foo'
        cfg = parse_args(['--cmsearch', cmsearch, 'replicon'])
        self.assertEqual(cfg.cmsearch, cmsearch)

    def test_hmmsearch(self):
        cfg = parse_args(['replicon'])
        self.assertEqual(cfg.hmmsearch, distutils.spawn.find_executable("hmmsearch"))
        hmmsearch = 'foo'
        cfg = parse_args(['--hmmsearch', hmmsearch, 'replicon'])
        self.assertEqual(cfg.hmmsearch, hmmsearch)

    def test_prodigal(self):
        cfg = parse_args(['replicon'])
        prodigal = 'foo'
        self.assertEqual(cfg.prodigal, distutils.spawn.find_executable('prodigal'))
        cfg = parse_args(['--prodigal', prodigal, 'replicon'])
        self.assertEqual(cfg.prodigal, prodigal)

    def test_path_func_annot(self):
        cfg = parse_args(['replicon'])
        self.assertIsNone(cfg.path_func_annot)
        func_annot = 'foo'
        cfg = parse_args(['--path_func_annot', func_annot, 'replicon'])
        self.assertEqual(cfg.path_func_annot, func_annot)

    def test_gembase(self):
        cfg = parse_args(['replicon'])
        self.assertFalse(cfg.gembase)
        cfg = parse_args(['--gembase', 'replicon'])
        self.assertTrue(cfg.gembase)

    def test_attc_model(self):
        cfg = parse_args(['replicon'])
        self.assertEqual(cfg.model_attc_name, 'attc_4')
        model = 'foo'
        cfg = parse_args(['--attc_model', model, 'replicon'])
        self.assertEqual(cfg.model_attc_name, model)

    def test_evalues_attc(self):
        cfg = parse_args(['replicon'])
        self.assertEqual(cfg.evalue_attc, 1.0)
        evalue = 5.0
        cfg = parse_args(['--evalue_attc', str(evalue), 'replicon'])
        self.assertEqual(cfg.evalue_attc, evalue)

    def test_keep_palindromes(self):
        cfg = parse_args(['replicon'])
        self.assertFalse(cfg.keep_palindromes)
        cfg = parse_args(['--keep_palindromes', 'replicon'])
        self.assertTrue(cfg.keep_palindromes)

    def test_no_proteins(self):
        cfg = parse_args(['replicon'])
        self.assertFalse(cfg.no_proteins)
        cfg = parse_args(['--no_proteins', 'replicon'])
        self.assertTrue(cfg.no_proteins)

    def test_max_attc_size(self):
        cfg = parse_args(['replicon'])
        self.assertEqual(cfg.max_attc_size, 200)
        max_attc_size = 50
        cfg = parse_args(['--max_attc_size', str(max_attc_size), 'replicon'])
        self.assertEqual(cfg.max_attc_size, max_attc_size)


    def test_min_attc_size(self):
        cfg = parse_args(['replicon'])
        self.assertEqual(cfg.min_attc_size, 40)
        min_attc_size = 50
        cfg = parse_args(['--min_attc_size', str(min_attc_size), 'replicon'])
        self.assertEqual(cfg.min_attc_size, min_attc_size)

    def test_eagle_eyes(self):
        cfg = parse_args(['replicon'])
        self.assertFalse(cfg.local_max)
        # eagle option is an "alias" for local_max
        cfg = parse_args(['--eagle_eyes', 'replicon'])
        self.assertTrue(cfg.local_max)

    def test_circular(self):
        cfg = parse_args(['replicon'])
        self.assertFalse(cfg.circular)
        cfg = parse_args(['--circ', 'replicon'])
        self.assertTrue(cfg.circular)

    def test_linear(self):
        cfg = parse_args(['replicon'])
        self.assertFalse(cfg.linear)
        cfg = parse_args(['--linear', 'replicon'])
        self.assertTrue(cfg.linear)

    def test_cirular_n_linear(self):
        real_exit = sys.exit

        sys.exit = self.fake_exit
        with self.catch_io(err=True):
            try:
                _ = parse_args(['--circ', '--linear', 'replicon'])
            except TypeError as err:
                msg = sys.stderr.getvalue()
                msg_end = 'error: argument --linear: not allowed with argument --circ\n'
                self.assertTrue(msg.endswith(msg_end), "{} != {}".format(msg[len(msg) - len(msg_end):], msg_end))
                # program exit with returncode = 2
                self.assertEqual(str(err), '2')
            finally:
                sys.exit = real_exit

    def test_topology_file(self):
        cfg = parse_args(['replicon'])
        self.assertIsNone(cfg.topology_file)
        topo = 'foo'
        cfg = parse_args(['--topology_file', topo, 'replicon'])
        self.assertEqual(cfg.topology_file, topo)

    def test_verbose(self):
        cfg = parse_args(['replicon'])
        self.assertEqual(cfg.verbose, 0)
        cfg = parse_args(['--verbose', 'replicon'])
        self.assertEqual(cfg.verbose, 1)
        cfg = parse_args(['-vv', 'replicon'])
        self.assertEqual(cfg.verbose, 2)

    def test_quiet(self):
        cfg = parse_args(['replicon'])
        self.assertEqual(cfg.quiet, 0)
        cfg = parse_args(['--quiet', 'replicon'])
        self.assertEqual(cfg.quiet, 1)
        cfg = parse_args(['-qq', 'replicon'])
        self.assertEqual(cfg.quiet, 2)

    def test_pdf(self):
        cfg = parse_args(['replicon'])
        self.assertFalse(cfg.pdf)
        cfg = parse_args(['--pdf', 'replicon'])
        self.assertTrue(cfg.pdf)

    def test_gbk(self):
        cfg = parse_args(['replicon'])
        self.assertFalse(cfg.gbk)
        cfg = parse_args(['--gbk', 'replicon'])
        self.assertTrue(cfg.gbk)

    def test_keep_tmp(self):
        cfg = parse_args(['replicon'])
        self.assertFalse(cfg.keep_tmp)
        cfg = parse_args(['--keep-tmp', 'replicon'])
        self.assertTrue(cfg.keep_tmp)

    def test_version(self):
        real_exit = sys.exit

        sys.exit = self.fake_exit
        with self.catch_io(err=True):
            try:
                _ = cfg = parse_args(['--version'])
            except TypeError as err:
                msg = sys.stderr.getvalue()
                msg_expected = """integron_finder version NOT packaged, it should be development version
Python {}

 - open-source GPLv3,
 - Jean Cury, Bertrand Neron, Eduardo Rocha,
 - citation:

 Identification and analysis of integrons and cassette arrays in bacterial genomes
 Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha
 Nucleic Acids Research 2016; doi: 10.1093/nar/gkw319
 """.format(sys.version)

                self.assertEqual(msg.strip(), msg_expected.strip())
                # program exit with returncode = 2
                self.assertEqual(str(err), '0')
            finally:
                sys.exit = real_exit

