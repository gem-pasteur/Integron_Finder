# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright (c) 2015 - 2018  Institut Pasteur, Paris and CNRS.                     #
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
import argparse
import tempfile
import distutils.spawn
import shutil

from Bio import SeqIO, Seq

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.config import Config
from integron_finder.utils import FastaIterator, read_multi_prot_fasta
from integron_finder.prot_db import GembaseDB, ProdigalDB, SeqDesc


class TestGemBase(IntegronTest):

    def setUp(self):
        """
        Define variables common to all tests
        """
        # Simulate argparse to get argument
        self.args = argparse.Namespace()
        self.args.gembase = True
        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        self.args.outdir = self.tmp_dir

        if os.path.exists(self.tmp_dir) and os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)


    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_dir)
            pass
        except:
            pass


    def test_ProteinDB(self):
        # From Gembase Draft , Gembase Complete
        file_names = ('ACBA.0917.00019.fna', 'ESCO001.C.00001.C001.fst')
        for file_name in file_names:
            replicon_path = self.find_data(os.path.join('Gembase', 'Replicons', file_name))
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            os.makedirs(cfg.tmp_dir)

            seq_db = FastaIterator(replicon_path)
            replicon = next(seq_db)
            db = GembaseDB(replicon, cfg)
            self.assertTrue(db.replicon.id, replicon.id)



    # def test_make_protfile(self):
    #     file_name = 'ACBA.0917.00019'
    #     replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
    #     self.args.replicon = replicon_path
    #     cfg = Config(self.args)
    #     os.makedirs(cfg.tmp_dir)
    #
    #     seq_db = FastaIterator(replicon_path)
    #     replicon = next(seq_db)
    #
    #     db = ProdigalDB(replicon, cfg)
    #     for seq_nb, seqs in enumerate(zip(
    #             read_multi_prot_fasta(self.find_data(os.path.join('Proteins', file_name +'.prt'))),
    #             read_multi_prot_fasta(db.protfile)), 1):
    #         expected, test = seqs
    #         self.assertEqual(expected.id, test.id)
    #     self.assertEqual(seq_nb, 23)


    # def test_protfile(self):
    #     file_name = 'acba.007.p01.13'
    #     prot_name = 'ACBA.007.P01_13.prt'
    #     replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
    #     self.args.replicon = replicon_path
    #     cfg = Config(self.args)
    #     os.makedirs(cfg.tmp_dir)
    #
    #     seq_db = FastaIterator(replicon_path)
    #     replicon = next(seq_db)
    #
    #     db = ProdigalDB(replicon, cfg)
    #     self.assertEqual(os.path.join(cfg.tmp_dir, prot_name), db.protfile)
    #
    #
    # def test_getitem(self):
    #     file_name = 'acba.007.p01.13'
    #     prot_name = 'ACBA.007.P01_13.prt'
    #     replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
    #     self.args.replicon = replicon_path
    #     cfg = Config(self.args)
    #     os.makedirs(cfg.tmp_dir)
    #
    #     seq_db = FastaIterator(replicon_path)
    #     replicon = next(seq_db)
    #
    #     db = ProdigalDB(replicon, cfg)
    #     exp = read_multi_prot_fasta(self.find_data(os.path.join('Proteins', prot_name)))
    #     for prot_expected in exp:
    #         prot_received = db[prot_expected.id]
    #         self.assertEqual(prot_received.id,
    #                          prot_expected.id)
    #         self.assertEqual(prot_received.seq,
    #                          prot_expected.seq)
    #     with self.assertRaises(KeyError) as ctx:
    #         db['nimport_naoik']
    #     self.assertEqual(str(ctx.exception), "'nimport_naoik'")
    #
    #
    # def test_iter(self):
    #     file_name = 'acba.007.p01.13'
    #     prot_name = 'ACBA.007.P01_13.prt'
    #     replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
    #     self.args.replicon = replicon_path
    #     cfg = Config(self.args)
    #     os.makedirs(cfg.tmp_dir)
    #
    #     seq_db = FastaIterator(replicon_path)
    #     replicon = next(seq_db)
    #
    #     db = ProdigalDB(replicon, cfg)
    #     idx = SeqIO.index(self.find_data(os.path.join('Proteins', prot_name)), 'fasta',
    #                       alphabet=Seq.IUPAC.extended_protein)
    #     for exp_seq_id, get_seq_id in zip(idx, db):
    #         self.assertEqual(exp_seq_id, get_seq_id)
    #
    #
    # def test_get_description(self):
    #     # SeqDesc(id, strand, strat, stop)
    #     file_name = 'acba.007.p01.13'
    #     replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
    #     self.args.replicon = replicon_path
    #     cfg = Config(self.args)
    #     os.makedirs(cfg.tmp_dir)
    #
    #     seq_db = FastaIterator(replicon_path)
    #     replicon = next(seq_db)
    #
    #     db = ProdigalDB(replicon, cfg)
    #
    #     descriptions = {'ACBA.007.P01_13_23': SeqDesc('ACBA.007.P01_13_23', -1, 19721, 20254),
    #                     'ACBA.007.P01_13_1': SeqDesc('ACBA.007.P01_13_1', 1, 55, 1014)}
    #     for seq_id, desc in descriptions.items():
    #         self.assertEqual(desc, db.get_description(seq_id))
    #     with self.assertRaises(KeyError) as ctx:
    #         db.get_description('nimport_naoik')
    #     self.assertEqual(str(ctx.exception), "'nimport_naoik'")


class TestProdigalDB(IntegronTest):

    def setUp(self):
        """
        Define variables common to all tests
        """
        # Simulate argparse to get argument
        self.args = argparse.Namespace()
        self.args.gembase = False
        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        self.args.outdir = self.tmp_dir
        self.args.prodigal = distutils.spawn.find_executable("prodigal")

        if os.path.exists(self.tmp_dir) and os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)


    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_dir)
            pass
        except:
            pass

    def test_ProteinDB(self):
        file_name = 'acba.007.p01.13'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        os.makedirs(cfg.tmp_dir)

        seq_db = FastaIterator(replicon_path)
        replicon = next(seq_db)
        db = ProdigalDB(replicon, cfg)
        self.assertTrue(db.replicon.id, replicon.id)


    def test_make_protfile(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        os.makedirs(cfg.tmp_dir)

        seq_db = FastaIterator(replicon_path)
        replicon = next(seq_db)

        db = ProdigalDB(replicon, cfg)
        for seq_nb, seqs in enumerate(zip(
                read_multi_prot_fasta(self.find_data(os.path.join('Proteins', prot_name))),
                read_multi_prot_fasta(db.protfile)), 1):
            expected, test = seqs
            self.assertEqual(expected.id, test.id)
        self.assertEqual(seq_nb, 23)

    def test_protfile(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        os.makedirs(cfg.tmp_dir)

        seq_db = FastaIterator(replicon_path)
        replicon = next(seq_db)

        db = ProdigalDB(replicon, cfg)
        self.assertEqual(os.path.join(cfg.tmp_dir, prot_name), db.protfile)


    def test_getitem(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        os.makedirs(cfg.tmp_dir)

        seq_db = FastaIterator(replicon_path)
        replicon = next(seq_db)

        db = ProdigalDB(replicon, cfg)
        exp = read_multi_prot_fasta(self.find_data(os.path.join('Proteins', prot_name)))
        for prot_expected in exp:
            prot_received = db[prot_expected.id]
            self.assertEqual(prot_received.id,
                             prot_expected.id)
            self.assertEqual(prot_received.seq,
                             prot_expected.seq)
        with self.assertRaises(KeyError) as ctx:
            db['nimport_naoik']
        self.assertEqual(str(ctx.exception), "'nimport_naoik'")


    def test_iter(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        os.makedirs(cfg.tmp_dir)

        seq_db = FastaIterator(replicon_path)
        replicon = next(seq_db)

        db = ProdigalDB(replicon, cfg)
        idx = SeqIO.index(self.find_data(os.path.join('Proteins', prot_name)), 'fasta',
                          alphabet=Seq.IUPAC.extended_protein)
        for exp_seq_id, get_seq_id in zip(idx, db):
            self.assertEqual(exp_seq_id, get_seq_id)


    def test_get_description(self):
        # SeqDesc(id, strand, strat, stop)
        file_name = 'acba.007.p01.13'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        os.makedirs(cfg.tmp_dir)

        seq_db = FastaIterator(replicon_path)
        replicon = next(seq_db)

        db = ProdigalDB(replicon, cfg)

        descriptions = {'ACBA.007.P01_13_23': SeqDesc('ACBA.007.P01_13_23', -1, 19721, 20254),
                        'ACBA.007.P01_13_1':  SeqDesc('ACBA.007.P01_13_1', 1, 55, 1014)}
        for seq_id, desc in descriptions.items():
            self.assertEqual(desc, db.get_description(seq_id))
        with self.assertRaises(KeyError) as ctx:
            db.get_description('nimport_naoik')
        self.assertEqual(str(ctx.exception), "'nimport_naoik'")