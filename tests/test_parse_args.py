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

import os
import shutil
import sys

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.scripts.finder import parse_args


class TestParseArgs(IntegronTest):

    def setUp(self):
        self.replicon = self.find_data('Replicons', 'acba_short.fst')

    def test_replicon(self):
        # whatever the replicon the path must exists
        # otherwise the parser raise an error as it check the path
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.input_seq_path, self.replicon)

    def test_wo_replicon(self):
        real_exit = sys.exit

        sys.exit = self.fake_exit
        with self.catch_io(err=True):
            try:
                _ = parse_args([])
            except TypeError as err:
                msg = sys.stderr.getvalue()
                msg_end = 'error: the following arguments are required: replicon\n'
                self.assertTrue(msg.endswith(msg_end), "{} != {}".format(msg[len(msg) - len(msg_end):], msg_end))
                # program exit with returncode = 2
                self.assertEqual(str(err), '2')
            finally:
                sys.exit = real_exit

    def test_local_max(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.local_max)
        cfg = parse_args(['--local-max', self.replicon])
        self.assertTrue(cfg.local_max)

    def test_func_annot(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.func_annot)
        cfg = parse_args(['--func-annot', self.replicon])
        self.assertTrue(cfg.func_annot)

    def test_cpu(self):
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.cpu, 1)
        cpu = 10
        cfg = parse_args(['--cpu', str(cpu), self.replicon])
        self.assertEqual(cfg.cpu, cpu)

    def test_distance_threshold(self):
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.distance_threshold, 4000)
        dt = 50
        cfg = parse_args(['--distance-thresh', str(dt), self.replicon])
        self.assertEqual(cfg.distance_threshold, dt)

    def test_outdir(self):
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.outdir, os.path.abspath('.'))
        outdir = 'foo'
        cfg = parse_args(['--outdir', outdir, self.replicon])
        self.assertEqual(cfg.outdir, os.path.abspath(outdir))

    def test_union_integrase(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.union_integrases)
        cfg = parse_args(['--union-integrases', self.replicon])
        self.assertTrue(cfg.union_integrases)

    def test_cmsearch(self):
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.cmsearch, shutil.which("cmsearch"))
        cmsearch = self.find_data('Replicons', 'fake_seq.fst')  # whatever the value the file must exists
        cfg = parse_args(['--cmsearch', cmsearch, self.replicon])
        self.assertEqual(cfg.cmsearch, cmsearch)

    def test_hmmsearch(self):
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.hmmsearch, shutil.which("hmmsearch"))
        hmmsearch = self.find_data('Replicons', 'fake_seq.fst')  # whatever the value the file must exists
        cfg = parse_args(['--hmmsearch', hmmsearch, self.replicon])
        self.assertEqual(cfg.hmmsearch, hmmsearch)

    def test_prodigal(self):
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.prodigal, shutil.which('prodigal'))
        prodigal = self.find_data('Replicons', 'fake_seq.fst')  # whatever the value the file must exists
        cfg = parse_args(['--prodigal', prodigal, self.replicon])
        self.assertEqual(cfg.prodigal, prodigal)

    def test_path_func_annot(self):
        cfg = parse_args([self.replicon])
        self.assertIsNone(cfg.path_func_annot)
        func_annot = self.find_data('hmm_files')  # whatever the value the file must exists
        cfg = parse_args(['--path-func-annot', func_annot, self.replicon])
        self.assertEqual(cfg.path_func_annot, func_annot)

    def test_gembase(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.gembase)
        cfg = parse_args(['--gembase', self.replicon])
        self.assertTrue(cfg.gembase)

    def test_attc_model(self):
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.model_attc_name, 'attc_4')
        model = 'foo'
        cfg = parse_args(['--attc-model', model, self.replicon])
        self.assertEqual(cfg.model_attc_name, model)

    def test_evalues_attc(self):
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.evalue_attc, 1.0)
        evalue = 5.0
        cfg = parse_args(['--evalue-attc', str(evalue), self.replicon])
        self.assertEqual(cfg.evalue_attc, evalue)

    def test_keep_palindromes(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.keep_palindromes)
        cfg = parse_args(['--keep-palindromes', self.replicon])
        self.assertTrue(cfg.keep_palindromes)

    def test_no_proteins(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.no_proteins)
        cfg = parse_args(['--no-proteins', self.replicon])
        self.assertTrue(cfg.no_proteins)

    def test_max_attc_size(self):
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.max_attc_size, 200)
        max_attc_size = 50
        cfg = parse_args(['--max-attc-size', str(max_attc_size), self.replicon])
        self.assertEqual(cfg.max_attc_size, max_attc_size)


    def test_min_attc_size(self):
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.min_attc_size, 40)
        min_attc_size = 50
        cfg = parse_args(['--min-attc-size', str(min_attc_size), self.replicon])
        self.assertEqual(cfg.min_attc_size, min_attc_size)

    def test_eagle_eyes(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.local_max)
        # eagle option is an "alias" for local_max
        cfg = parse_args(['--eagle-eyes', self.replicon])
        self.assertTrue(cfg.local_max)

    def test_circular(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.circular)
        cfg = parse_args(['--circ', self.replicon])
        self.assertTrue(cfg.circular)

    def test_linear(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.linear)
        cfg = parse_args(['--linear', self.replicon])
        self.assertTrue(cfg.linear)

    def test_cirular_n_linear(self):
        real_exit = sys.exit

        sys.exit = self.fake_exit
        with self.catch_io(err=True):
            try:
                _ = parse_args(['--circ', '--linear', self.replicon])
            except TypeError as err:
                msg = sys.stderr.getvalue()
                msg_end = 'error: argument --linear: not allowed with argument --circ\n'
                self.assertTrue(msg.endswith(msg_end), "{} != {}".format(msg[len(msg) - len(msg_end):], msg_end))
                # program exit with returncode = 2
                self.assertEqual(str(err), '2')
            finally:
                sys.exit = real_exit

    def test_topology_file(self):
        cfg = parse_args([self.replicon])
        self.assertIsNone(cfg.topology_file)
        topo = self.find_data('topology.txt')  # whatever the value the file must exists
        cfg = parse_args(['--topology-file', topo, self.replicon])
        self.assertEqual(cfg.topology_file, topo)

    def test_mute(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.mute)
        cfg = parse_args(['--mute', self.replicon])
        self.assertTrue(cfg.mute)

    def test_verbose(self):
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.verbose, 0)
        cfg = parse_args(['--verbose', self.replicon])
        self.assertEqual(cfg.verbose, 1)
        cfg = parse_args(['-vv', self.replicon])
        self.assertEqual(cfg.verbose, 2)

    def test_quiet(self):
        cfg = parse_args([self.replicon])
        self.assertEqual(cfg.quiet, 0)
        cfg = parse_args(['--quiet', self.replicon])
        self.assertEqual(cfg.quiet, 1)
        cfg = parse_args(['-qq', self.replicon])
        self.assertEqual(cfg.quiet, 2)

    def test_pdf(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.pdf)
        cfg = parse_args(['--pdf', self.replicon])
        self.assertTrue(cfg.pdf)

    def test_gbk(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.gbk)
        cfg = parse_args(['--gbk', self.replicon])
        self.assertTrue(cfg.gbk)

    def test_keep_tmp(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.keep_tmp)
        cfg = parse_args(['--keep-tmp', self.replicon])
        self.assertTrue(cfg.keep_tmp)

    def test_split_results(self):
        cfg = parse_args([self.replicon])
        self.assertFalse(cfg.split_results)
        cfg = parse_args(['--split-results', self.replicon])
        self.assertTrue(cfg.split_results)


    def test_version(self):
        real_exit = sys.exit
        sys.exit = self.fake_exit

        from numpy import __version__ as np_vers
        from pandas import __version__ as pd_vers
        from matplotlib import __version__ as mplt_vers
        from Bio import __version__ as bio_vers
        import integron_finder

        with self.catch_io(out=True):
            try:
                _ = parse_args(['--version'])
            except TypeError:
                msg = sys.stdout.getvalue()
                msg_expected = """integron_finder version {i_f} {commit}
Using:
 - Python {py}
 - numpy {np}
 - pandas {pd}
 - matplolib {mplt}
 - biopython {bio}

 - {prodigal}
 - {cmsearch}
 - {hmmsearch}

Authors:
 - Jean Cury, Bertrand Neron, Eduardo Rocha,

Citation:

 Néron, B.; Littner, E.; Haudiquet, M.; Perrin, A.; Cury, J.; Rocha, E.P.C.
 IntegronFinder 2.0: Identification and Analysis of Integrons across Bacteria, with a Focus on Antibiotic Resistance in Klebsiella.
 Microorganisms 2022, 10, 700. https://doi.org/10.3390/microorganisms10040700

 If you use --func-annot in conjunction with file NCBIfam-AMRFinder.hmm please also cite

 Haft, DH et al., Nucleic Acids Res. 2018 Jan 4;46(D1):D851-D860
 PMID: 29112715
""".format(i_f=integron_finder.__version__,
           commit=integron_finder.__commit__ if 'dev' in integron_finder.__version__ else '',
           py=sys.version.replace('\n', ' '),
           np=np_vers,
           pd=pd_vers,
           mplt=mplt_vers,
           bio=bio_vers,
           prodigal=integron_finder._prodigal_version(shutil.which("prodigal")),
           cmsearch=integron_finder._eddy_version(shutil.which("cmsearch")),
           hmmsearch=integron_finder._eddy_version(shutil.which("hmmsearch"))
           )
            finally:
                sys.exit = real_exit
            self.assertEqual(msg.strip(), msg_expected.strip())
