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
import glob

from Bio import SeqIO, Seq

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

import integron_finder.scripts.split as split


class TestSplit(IntegronTest):

    def setUpClass():
        split._log= split.colorlog.getLogger('integron_finder.split')


    def setUp(self):
        self._tmp_dir = tempfile.TemporaryDirectory(prefix='tmp_test_integron_finder')
        self.tmp_dir = self._tmp_dir.name
        self.out_dir = os.path.join(self.tmp_dir, 'test_integron_split')
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
        os.makedirs(self.out_dir)


    def tearDown(self):
        self._tmp_dir.cleanup()
        

    def test_split_wo_chunk(self):
        replicon_path = self.find_data(os.path.join('Replicons', 'multi_fasta.fst'))
        chunk_names = split.split(replicon_path, outdir=self.out_dir)

        try:
            seq_index = SeqIO.index(replicon_path, "fasta")
            files_expected = [os.path.join(self.out_dir, r + '.fst') for r in seq_index]
            self.assertListEqual(files_expected, chunk_names)
            for one_chunk in chunk_names:
                with open(one_chunk) as f:
                    seq_it = SeqIO.parse(f, 'fasta')
                    for s in seq_it:
                        ref_seq = seq_index[s.id]
                        self.assertEqual(s.id, ref_seq.id)
                        self.assertEqual(s.description, ref_seq.description)
                        self.assertEqual(s.seq, ref_seq.seq)
        finally:
            seq_index.close()


    def test_split_avoid_overwriting(self):
        replicon_path = self.find_data(os.path.join('Replicons', 'ESCO001.B.00018.P002.fst'))
        chunk_names = split.split(replicon_path, outdir=self.out_dir)
        files_expected = [os.path.join(self.out_dir, 'ESCO001.B.00018.P002.fst')]
        self.assertListEqual(files_expected, chunk_names)

        chunk_names = split.split(replicon_path, outdir=self.out_dir)
        files_expected = [os.path.join(self.out_dir, 'ESCO001.B.00018.P002_chunk_1.fst')]
        self.assertListEqual(files_expected, chunk_names)

    def test_split_w_chunk(self):
        replicon_path = self.find_data(os.path.join('Replicons', 'multi_fasta.fst'))
        chunk = 2
        chunk_names = split.split(replicon_path, outdir=self.out_dir, chunk=chunk)

        try:
            seq_index = SeqIO.index(replicon_path, "fasta")
            files_expected = [os.path.join(self.out_dir, "multi_fasta_chunk_{}.fst".format(i)) for i in range(1, chunk + 1)]
            self.assertListEqual(files_expected, chunk_names)
            for one_chunk in chunk_names:
                with open(one_chunk) as f:
                    seq_it = SeqIO.parse(f, 'fasta')
                    for s in seq_it:
                        ref_seq = seq_index[s.id]
                        self.assertEqual(s.id, ref_seq.id)
                        self.assertEqual(s.description, ref_seq.description)
                        self.assertEqual(s.seq, ref_seq.seq)
        finally:
            seq_index.close()


class TestParseArgs(IntegronTest):

    def test_parse_replicon(self):
        parsed_args = split.parse_args(['replicon'])
        self.assertEqual(parsed_args.chunk, None)
        self.assertEqual(parsed_args.outdir, '.')
        self.assertEqual(parsed_args.quiet, 0)
        self.assertEqual(parsed_args.verbose, 0)
        self.assertEqual(parsed_args.replicon, 'replicon')

    def test_parse_outdir(self):
        parsed_args = split.parse_args(['--outdir', 'foo', 'replicon'])
        self.assertEqual(parsed_args.chunk, None)
        self.assertEqual(parsed_args.outdir, 'foo')
        self.assertEqual(parsed_args.quiet, 0)
        self.assertEqual(parsed_args.verbose, 0)
        self.assertEqual(parsed_args.replicon, 'replicon')

    def test_parse_chunk(self):
        parsed_args = split.parse_args(['--outdir', 'foo', '--chunk', '10', 'replicon'])
        self.assertEqual(parsed_args.chunk, 10)
        self.assertEqual(parsed_args.outdir, 'foo')
        self.assertEqual(parsed_args.quiet, 0)
        self.assertEqual(parsed_args.verbose, 0)
        self.assertEqual(parsed_args.replicon, 'replicon')

    def test_mute(self):
        parsed_args = split.parse_args(['replicon'])
        self.assertFalse(parsed_args.mute)
        parsed_args = split.parse_args(['--mute', 'replicon'])
        self.assertTrue(parsed_args.mute)

    def test_verbose(self):
        parsed_args = split.parse_args(['--outdir', 'foo', '--chunk', '10', 'replicon'])
        self.assertEqual(parsed_args.verbose, 0)
        parsed_args = split.parse_args(['--outdir', 'foo', '--chunk', '10', '--verbose', 'replicon'])
        self.assertEqual(parsed_args.verbose, 1)
        parsed_args = split.parse_args(['--outdir', 'foo', '--chunk', '10', '-vv', 'replicon'])
        self.assertEqual(parsed_args.verbose, 2)

    def test_quiet(self):
        parsed_args = split.parse_args(['--outdir', 'foo', '--chunk', '10', 'replicon'])
        self.assertEqual(parsed_args.quiet, 0)
        parsed_args = split.parse_args(['--outdir', 'foo', '--chunk', '10', '--quiet', 'replicon'])
        self.assertEqual(parsed_args.quiet, 1)
        parsed_args = split.parse_args(['--outdir', 'foo', '--chunk', '10', '-qq', 'replicon'])
        self.assertEqual(parsed_args.quiet, 2)


class TestMain(IntegronTest):

    def setUp(self):
        tmp_dir = tempfile.gettempdir()
        self.out_dir = os.path.join(tmp_dir, 'test_integron_split')
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)
        os.makedirs(self.out_dir)

    def tearDown(self):
        if os.path.exists(self.out_dir) and os.path.isdir(self.out_dir):
            shutil.rmtree(self.out_dir)

    def test_out_dir_isnot_dir(self):
        out_dir_not_dir  = os.path.join(self.out_dir, "foo_file.txt")
        open(out_dir_not_dir, 'w').close()
        replicon_path = self.find_data(os.path.join('Replicons', 'multi_fasta.fst'))
        command = 'integron_split --outdir {} {}'.format(out_dir_not_dir, replicon_path)
        with self.assertRaises(RuntimeError) as ctx:
            split.main(command.split()[1:], log_level="WARNING")
        self.assertEqual(str(ctx.exception),
                         "The outdir '/tmp/test_integron_split/foo_file.txt' already exist and is not a directory.")

    def test_wo_chunk(self):
        replicon_path = self.find_data(os.path.join('Replicons', 'multi_fasta.fst'))
        command = 'integron_split --outdir {} {}'.format(self.out_dir, replicon_path)
        with self.catch_io(out=True, err=True):
            split.main(command.split()[1:], log_level="WARNING")
        try:
            seq_index = SeqIO.index(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)
        except AttributeError:
            seq_index = SeqIO.index(replicon_path, "fasta")
        files_expected = sorted([os.path.join(self.out_dir, r + '.fst') for r in seq_index])
        chunk_names = sorted(glob.glob(os.path.join(self.out_dir, '*.fst')))
        self.assertListEqual(files_expected, chunk_names)
        for f in chunk_names:
            try:
                seq_it = SeqIO.parse(f, 'fasta')
                for s in seq_it:
                    ref_seq = seq_index[s.id]
                    self.assertEqual(s.id, ref_seq.id)
                    self.assertEqual(s.description, ref_seq.description)
                    self.assertEqual(s.seq, ref_seq.seq)
            finally:
                seq_it.stream.close()
        seq_index.close()


    def test_w_chunk(self):
        replicon_path = self.find_data(os.path.join('Replicons', 'multi_fasta.fst'))
        chunk = 2
        command = 'integron_split --outdir {} --chunk {} {}'.format(self.out_dir, chunk, replicon_path)
        with self.catch_io(out=True, err=True):
            split.main(command.split()[1:], log_level="WARNING")
        try:
            seq_index = SeqIO.index(replicon_path, "fasta", alphabet=Seq.IUPAC.unambiguous_dna)
        except AttributeError:
            seq_index = SeqIO.index(replicon_path, "fasta")

        files_expected = sorted([os.path.join(self.out_dir, "multi_fasta_chunk_{}.fst".format(i))
                                 for i in range(1, chunk + 1)])
        chunk_names = sorted(glob.glob(os.path.join(self.out_dir, '*.fst')))
        self.assertListEqual(files_expected, chunk_names)
        for f in chunk_names:
            try:
                seq_it = SeqIO.parse(f, 'fasta')
                for s in seq_it:
                    ref_seq = seq_index[s.id]
                    self.assertEqual(s.id, ref_seq.id)
                    self.assertEqual(s.description, ref_seq.description)
                    self.assertEqual(s.seq, ref_seq.seq)
            finally:
                seq_it.stream.close()
        seq_index.close()
