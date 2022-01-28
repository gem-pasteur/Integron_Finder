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
import argparse
import tempfile
import distutils.spawn
import shutil
import re

from Bio import SeqIO, Seq

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import IntegronError
from integron_finder.config import Config
from integron_finder.utils import read_multi_prot_fasta
from integron_finder.prot_db import GembaseDB, ProdigalDB, SeqDesc, CustomDB


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
        except:
            pass


    def test_ProteinDB(self):
        # From Gembase Draft , Gembase Complete
        file_names = ('ACBA.0917.00019.fna', 'ESCO001.C.00001.C001.fst')
        for file_name in file_names:
            replicon_path = self.find_data(os.path.join('Gembase', 'Replicons', file_name))
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            seq_db = read_multi_prot_fasta(replicon_path)
            replicon = next(seq_db)
            replicon.path = replicon_path
            os.makedirs(cfg.tmp_dir(replicon.id))

            with self.catch_log():
                db = GembaseDB(replicon, cfg)
            self.assertTrue(db.replicon.id, replicon.id)

    def test_find_gembase_file_basename(self):
        """
        test if find_gembase_file_basename get the the right basename
        for files in gembase
        """
        gembase_path = self.find_data('Gembase')
        file_names = ('ACBA.0917.00019.fna', 'ESCO001.C.00001.C001.fst')
        print()
        for file_name in file_names:
            replicon_path = self.find_data(os.path.join('Gembase', 'Replicons', file_name))
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            seq_db = read_multi_prot_fasta(replicon_path)
            replicon = next(seq_db)
            replicon.path = replicon_path
            os.makedirs(cfg.tmp_dir(replicon.id))

            with self.subTest(file_name=file_name):
                with self.catch_log() as log:
                    db = GembaseDB(replicon, cfg)
                    catch_msg = log.get_value().strip()
                self.assertTrue(db._find_gembase_file_basename,
                                os.path.splitext(file_name)[0])


    def test_find_gembase_file_basename_file_not_in_gembase(self):
        """
        test if find_gembase_file_basename get the the right basename
        for files not located in gembase and file name is the output of split operation
        a file containing one contig
        a file representing a chunk
        """
        gembase_path = self.find_data('Gembase')

        file_names = {'ACBA.0917.00019': self.find_data(os.path.join('Replicons', 'ACBA.0917.00019.0001.fst')),
                      'ESCO001.C.00001.C001.fst': os.path.join(self.tmp_dir, 'ESCO001.C.00001.C001_chunk_1.fst')
                      }

        shutil.copyfile(os.path.join(gembase_path, 'Replicons', 'ESCO001.C.00001.C001.fst'),
                        file_names['ESCO001.C.00001.C001.fst'])

        for base_file_name, replicon_path in file_names.items():
            self.args.replicon = replicon_path
            self.args.gembase_path = gembase_path
            cfg = Config(self.args)
            seq_db = read_multi_prot_fasta(replicon_path)
            replicon = next(seq_db)
            replicon.path = replicon_path
            os.makedirs(cfg.tmp_dir(replicon.id))
            with self.subTest(base_file_name=base_file_name):
                with self.catch_log():
                    db = GembaseDB(replicon, cfg, gembase_path=gembase_path)
                self.assertTrue(db._find_gembase_file_basename,
                                base_file_name)

        replicon_path = self.find_data(os.path.join('Replicons', 'acba.007.p01.13.fst'))
        self.args.replicon = replicon_path
        self.args.gembase_path = gembase_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        with self.assertRaises(FileNotFoundError) as ctx:
            with self.catch_log():
                GembaseDB(replicon, cfg, gembase_path=gembase_path)
        self.assertEqual(str(ctx.exception),
                         'cannot find lst file matching {} sequence'.format(replicon_path))


    def test_gembase_sniffer(self):
        file_names = (('ACBA.0917.00019', 'Draft'),
                      ('ESCO001.C.00001.C001', 'Complet'),
                      ('ACJO001.0321.00008.P008', 'Complet'))
        for file_name, gem_type in file_names:
            lst_path = self.find_data(os.path.join('Gembase', 'LSTINF', file_name + '.lst'))
            with self.subTest(lst_path=lst_path):
                type_recieved = GembaseDB.gembase_sniffer(lst_path)
                self.assertEqual(type_recieved, gem_type)

        lst_path = self.find_data(os.path.join('Gembase', 'LSTINF', 'SAEN001.0321.00753.P003.lst'))
        with self.assertRaises(IntegronError) as ctx:
            with self.catch_log():
                GembaseDB.gembase_sniffer(lst_path)
        self.assertEqual(str(ctx.exception),
                         f'the genome SAEN001.0321.00753.P003_00001 seems empty: see {lst_path}')

    def test_gembase_complete_parser(self):
        replicon_id = 'ESCO001.C.00001.C001'
        lst_path = self.find_data(os.path.join('Gembase', 'LSTINF', replicon_id + '.lst'))
        prots_info = GembaseDB.gembase_complete_parser(lst_path, replicon_id)
        columns = ['start', 'end', 'strand', 'type', 'seq_id', 'valid', 'gene_name', 'description']
        self.assertListEqual(list(prots_info.columns), columns)
        self.assertEqual(prots_info.shape, (4139, len(columns)))
        first_row = [190, 255, 'D', 'CDS', 'ESCO001.C.00001.C001_00001', 'Valid', 'thrL',
                     '@b0001@NP_414542.1@ b0001 1 190 255 | leader; Amino acid biosynthesis:'
                     ' Threonine thr operon leader peptide | ..']

        recieved_first_row = prots_info.iloc[0].values.tolist()
        self.assertListEqual(first_row, recieved_first_row)

        last_row = [4640942, 4641628, 'D', 'CDS', 'ESCO001.C.00001.C001_04495', 'Valid', 'yjtD',
                    '@b4403@NP_418820.1@ b4403 1 4640942 4641628 | putative methyltransferase | ..']
        recieved_last_row = prots_info.iloc[len(prots_info) - 1].values.tolist()
        self.assertListEqual(last_row, recieved_last_row)


    def test_gembase_draft_parser(self):
        replicon_name = 'ACBA.0917.00019'
        replicon_id = 'ACBA.0917.00019.0001'
        lst_path = self.find_data(os.path.join('Gembase', 'LSTINF', replicon_name + '.lst'))
        prots_info = GembaseDB.gembase_draft_parser(lst_path, replicon_id)
        columns = ['start', 'end', 'strand', 'type', 'seq_id', 'gene_name', 'description']
        self.assertListEqual(list(prots_info.columns), columns)
        self.assertEqual(prots_info.shape, (3870, len(columns)))
        first_row = [266, 1480, 'C', 'CDS', 'ACBA.0917.00019.b0001_00001', 'tyrS',
                     '| Tyrosine--tRNA ligase | 6.1.1.1 | similar to AA sequence:UniProtKB:P41256']

        recieved_first_row = prots_info.iloc[0].values.tolist()
        self.assertListEqual(first_row, recieved_first_row)

        last_row = [4043755, 4044354, 'C', 'CDS', 'ACBA.0917.00019.i0001_03957', 'yfcG_3',
                    '| Disulfide-bond oxidoreductase YfcG | 1.8.4.- | similar to AA sequence:UniProtKB:P77526']
        recieved_last_row = prots_info.iloc[len(prots_info) - 1].values.tolist()
        self.assertListEqual(last_row, recieved_last_row)


    def test_make_protfile(self):
        file_name = (('ACBA.0917.00019', '.fna', 3870), ('ESCO001.C.00001.C001', '.fst', 3870))
        for seq_name, ext, seq_nb in file_name:
            replicon_path = self.find_data(os.path.join('Gembase', 'Replicons', seq_name + ext))
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            seq_db = read_multi_prot_fasta(replicon_path)
            replicon = next(seq_db)
            replicon.path = replicon_path
            os.makedirs(cfg.tmp_dir(replicon.id))

            with self.catch_log():
                db = GembaseDB(replicon, cfg)
            for seq_nb, seqs in enumerate(zip(
                    read_multi_prot_fasta(self.find_data(os.path.join('Gembase', 'Proteins', seq_name + '.prt'))),
                    read_multi_prot_fasta(db.protfile)), 1):
                expected, test = seqs
                self.assertEqual(expected.id, test.id)
            self.assertEqual(seq_nb, seq_nb)


    def test_protfile(self):
        file_name = (('ACBA.0917.00019', '.fna'), ('ESCO001.C.00001.C001', '.fst'))
        for seq_name, ext in file_name:
            replicon_path = self.find_data(os.path.join('Gembase', 'Replicons', seq_name + ext))
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            seq_db = read_multi_prot_fasta(replicon_path)
            replicon = next(seq_db)
            replicon.path = replicon_path
            os.makedirs(cfg.tmp_dir(replicon.id))

            with self.catch_log():
                db = GembaseDB(replicon, cfg)
            self.assertEqual(os.path.join(cfg.tmp_dir(replicon.id), replicon.id + '.prt'), db.protfile)


    def test_getitem(self):
        file_name = (('ACBA.0917.00019', '.fna'), ('ESCO001.C.00001.C001', '.fst'))
        for seq_name, ext in file_name:
            replicon_path = self.find_data(os.path.join('Gembase', 'Replicons', seq_name + ext))
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            seq_db = read_multi_prot_fasta(replicon_path)
            replicon = next(seq_db)
            os.makedirs(cfg.tmp_dir(replicon.id))

            with self.catch_log():
                db = GembaseDB(replicon, cfg)
            exp = read_multi_prot_fasta(self.find_data(os.path.join('Gembase', 'Proteins', seq_name + '.prt')))

            specie, date, strain, contig = replicon.id.split('.')
            pattern = '{}\.{}\.{}\.\w?{}'.format(specie, date, strain, contig)

            for prot_expected in exp:
                if re.match(pattern, prot_expected.id):
                    prot_received = db[prot_expected.id]
                    self.assertEqual(prot_received.id,
                                     prot_expected.id)
                    self.assertEqual(prot_received.seq,
                                     prot_expected.seq)
        with self.assertRaises(KeyError) as ctx:
            db['nimport_naoik']
        self.assertEqual(str(ctx.exception), "'nimport_naoik'")


    def test_iter(self):
        # test Gembase Draft
        seq_name = 'ACBA.0917.00019'
        ext = '.fna'
        replicon_path = self.find_data(os.path.join('Gembase', 'Replicons', seq_name + ext))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        with self.catch_log():
            db = GembaseDB(replicon, cfg)

        try:
            idx = SeqIO.index(self.find_data(os.path.join('Gembase', 'Proteins', seq_name + '.prt')), 'fasta',
                              alphabet=Seq.IUPAC.extended_protein)
        except AttributeError:
            idx = SeqIO.index(self.find_data(os.path.join('Gembase', 'Proteins', seq_name + '.prt')), 'fasta')
        specie, date, strain, contig = replicon.id.split('.')
        pattern = '{}\.{}\.{}\.\w?{}'.format(specie, date, strain, contig)
        self.assertListEqual(sorted([i for i in idx if re.match(pattern, i)]), sorted([i for i in db]))

        # test Gembase Complet
        seq_name = 'ESCO001.C.00001.C001'
        ext = '.fst'
        replicon_path = self.find_data(os.path.join('Gembase', 'Replicons', seq_name + ext))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        with self.catch_log():
            db = GembaseDB(replicon, cfg)
        try:
            idx = SeqIO.index(self.find_data(os.path.join('Gembase', 'Proteins', seq_name + '.prt')), 'fasta',
                              alphabet=Seq.IUPAC.extended_protein)
        except AttributeError:
            idx = SeqIO.index(self.find_data(os.path.join('Gembase', 'Proteins', seq_name + '.prt')), 'fasta')

        specie, date, strain, contig = replicon.id.split('.')
        pattern = '{}\.{}\.{}\.\w?{}'.format(specie, date, strain, contig)
        seqid_from_gembase_protfile = set([i for i in idx if re.match(pattern, i)])
        seqid_from_if = set([i for i in db])
        non_common_seq = seqid_from_gembase_protfile ^ seqid_from_if
        # in Gembase complete the annotation from lstinfo provided from genbank
        # it appear some times that some CDS are not translate in proteins
        # So in data I have 3 genes from LSTINF are not in .prt file
        diff = {'ESCO001.C.00001.C001_03974', 'ESCO001.C.00001.C001_01509', 'ESCO001.C.00001.C001_04162'}
        self.assertSetEqual(non_common_seq, diff)


    def test_get_description(self):
        # SeqDesc(id, strand, strat, stop)
        file_name = {('ACBA.0917.00019', '.fna'):
                         {'ACBA.0917.00019.b0001_00001': SeqDesc('ACBA.0917.00019.b0001_00001', -1, 266, 1480),
                          'ACBA.0917.00019.i0001_03957': SeqDesc('ACBA.0917.00019.i0001_03957', -1, 4043755, 4044354)},
                     }

        for seq_name, ext in file_name:
            replicon_path = self.find_data(os.path.join('Gembase', 'Replicons', seq_name + ext))
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            seq_db = read_multi_prot_fasta(replicon_path)
            replicon = next(seq_db)
            replicon.path = replicon_path
            os.makedirs(cfg.tmp_dir(replicon.id))
            with self.catch_log():
                db = GembaseDB(replicon, cfg)

            descriptions = file_name[(seq_name, ext)]
            for seq_id, desc in descriptions.items():
                with self.subTest(seq_id=seq_id, desc=desc):
                    self.assertEqual(desc, db.get_description(seq_id))

        with self.assertRaises(IntegronError) as ctx:
            db.get_description('nimport_naoik')
        self.assertEqual(str(ctx.exception), "'nimport_naoik' is not a valid Gembase protein identifier.")

        with self.assertRaises(KeyError) as ctx:
            db.get_description('FOO.BAR.00019.i0001_03924')
        self.assertEqual(str(ctx.exception), "'FOO.BAR.00019.i0001_03924'")


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
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = ProdigalDB(replicon, cfg)
        self.assertTrue(db.replicon.id, replicon.id)


    def test_ProteinDB_no_prodigal(self):
        file_name = 'acba.007.p01.13'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        self.args.prodigal = None
        with self.assertRaises(RuntimeError) as ctx:
            ProdigalDB(replicon, cfg)


    def test_make_protfile(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'

        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = ProdigalDB(replicon, cfg)
        for seq_nb, seqs in enumerate(zip(
                read_multi_prot_fasta(self.find_data(os.path.join('Proteins', prot_name))),
                read_multi_prot_fasta(db.protfile)), 1):
            expected, test = seqs
            self.assertEqual(expected.id, test.id)
        self.assertEqual(seq_nb, 23)


    def test_make_protfile_no_dir(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path

        db = ProdigalDB(replicon, cfg)
        for seq_nb, seqs in enumerate(zip(
                read_multi_prot_fasta(self.find_data(os.path.join('Proteins', prot_name))),
                read_multi_prot_fasta(db.protfile)), 1):
            expected, test = seqs
            self.assertEqual(expected.id, test.id)
        self.assertEqual(seq_nb, 23)


    def test_make_protfile_no_prodigal(self):
        file_name = 'acba.007.p01.13'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        self.args.prodigal = 'foo_bar'
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path

        with self.assertRaises(RuntimeError) as ctx:
            ProdigalDB(replicon, cfg)

    def test_make_protfile_prodigal_failed(self):
        file_name = 'acba.007.p01.13'
        replicon_path = self.find_data('Replicons', file_name + '.fst')
        self.args.replicon = replicon_path
        self.args.prodigal = self.find_data('fake_prodigal')
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path

        with self.assertRaises(RuntimeError) as ctx:
            ProdigalDB(replicon, cfg)

        self.assertTrue(str(ctx.exception).strip().endswith(": failed : prodigal returncode = 50"))

    def test_protfile(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = ProdigalDB(replicon, cfg)
        self.assertEqual(os.path.join(cfg.tmp_dir(replicon.id), prot_name), db.protfile)


    def test_getitem(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = ProdigalDB(replicon, cfg)
        exp = read_multi_prot_fasta(self.find_data(os.path.join('Proteins', prot_name)))
        for prot_expected in exp:
            prot_received = db[prot_expected.id]
            self.assertEqual(prot_received.id,
                             prot_expected.id)
            self.assertEqual(prot_received.seq,
                             prot_expected.seq)
        with self.assertRaises(IntegronError) as ctx:
            gene_id = 'nimport_naoik'
            db[gene_id]
        self.assertEqual(str(ctx.exception),
                         f"protein file does not contains '{gene_id}' id. "
                         f"Try again with removing previous results dir {cfg.result_dir}"
                         )


    def test_iter(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = ProdigalDB(replicon, cfg)
        try:
            idx = SeqIO.index(self.find_data(os.path.join('Proteins', prot_name)), 'fasta',
                               alphabet=Seq.IUPAC.extended_protein)
        except AttributeError:
            idx = SeqIO.index(self.find_data(os.path.join('Proteins', prot_name)), 'fasta')
        for exp_seq_id, get_seq_id in zip(idx, db):
            self.assertEqual(exp_seq_id, get_seq_id)


    def test_get_description(self):
        # SeqDesc(id, strand, strat, stop)
        file_name = 'acba.007.p01.13'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = ProdigalDB(replicon, cfg)

        descriptions = {'ACBA.007.P01_13_23': SeqDesc('ACBA.007.P01_13_23', -1, 19721, 20254),
                        'ACBA.007.P01_13_1':  SeqDesc('ACBA.007.P01_13_1', 1, 55, 1014)}
        for seq_id, desc in descriptions.items():
            self.assertEqual(desc, db.get_description(seq_id))


class TestCustomDB(IntegronTest):

    def setUp(self):
        """
        Define variables common to all tests
        """
        # Simulate argparse to get argument
        self.args = argparse.Namespace()
        self.args.gembase = False
        self.tmp_dir = os.path.join(tempfile.gettempdir(), 'tmp_test_integron_finder')
        self.args.outdir = self.tmp_dir
        self.args.prodigal = None
        self.args.annot_parser = self.find_data('prodigal_annot_parser.py')

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
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data('Replicons', file_name + '.fst')
        protein_path = self.find_data('Proteins', prot_name)
        self.args.replicon = replicon_path
        self.args.prot_file = protein_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = CustomDB(replicon, cfg, protein_path)
        self.assertTrue(db.replicon.id, replicon.id)


    def test_ProteinDB_bad_parser(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        protein_path = self.find_data('Proteins', prot_name)
        self.args.replicon = replicon_path
        self.args.prot_file = protein_path
        self.args.annot_parser = self.find_data('df_max_input_1.csv')
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        with self.assertRaises(RuntimeError) as ctx:
            with self.catch_log():
                db = CustomDB(replicon, cfg, protein_path)
        self.assertEqual(str(ctx.exception),
                         f"Cannot import custom --annot-parser '{self.args.annot_parser}': "
                         f"'NoneType' object has no attribute 'loader'")


    def test_protfile(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data('Replicons', file_name + '.fst')
        protein_path = self.find_data('Proteins', prot_name)
        self.args.replicon = replicon_path
        self.args.prot_file = protein_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = CustomDB(replicon, cfg, protein_path)
        self.assertEqual(protein_path,
                         db.protfile)


    def test_getitem(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        protein_path = self.find_data('Proteins', prot_name)
        self.args.replicon = replicon_path
        self.args.prot_file = protein_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = CustomDB(replicon, cfg, protein_path)
        exp = read_multi_prot_fasta(self.find_data(os.path.join('Proteins', prot_name)))
        for prot_expected in exp:
            prot_received = db[prot_expected.id]
            self.assertEqual(prot_received.id,
                             prot_expected.id)
            self.assertEqual(prot_received.seq,
                             prot_expected.seq)
        with self.assertRaises(IntegronError) as ctx:
            gene_id = 'nimport_naoik'
            db[gene_id]
        self.assertEqual(str(ctx.exception),
                         f"protein file does not contains '{gene_id}' id. "
                         f"Check if it's the right proteins file {protein_path} "
                         f"or remove previous results dir {cfg.result_dir}"
                         )

    def test_iter(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        protein_path = self.find_data('Proteins', prot_name)
        self.args.replicon = replicon_path
        self.args.prot_file = protein_path
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = CustomDB(replicon, cfg, protein_path)
        try:
            idx = SeqIO.index(self.find_data(os.path.join('Proteins', prot_name)), 'fasta',
                               alphabet=Seq.IUPAC.extended_protein)
        except AttributeError:
            idx = SeqIO.index(self.find_data(os.path.join('Proteins', prot_name)), 'fasta')
        for exp_seq_id, get_seq_id in zip(idx, db):
            self.assertEqual(exp_seq_id, get_seq_id)


    def test_get_description(self):
        # SeqDesc(id, strand, strat, stop)
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        protein_path = self.find_data('Proteins', prot_name)
        self.args.replicon = replicon_path
        self.args.prot_file = protein_path
        self.args.annot_parser = self.find_data('prodigal_annot_parser.py')
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = CustomDB(replicon, cfg, protein_path)

        descriptions = {'ACBA.007.P01_13_23': SeqDesc('ACBA.007.P01_13_23', -1, 19721, 20254),
                        'ACBA.007.P01_13_1':  SeqDesc('ACBA.007.P01_13_1', 1, 55, 1014)}
        for seq_id, desc in descriptions.items():
            self.assertEqual(desc, db.get_description(seq_id))


    def test_get_description_stupid_parser(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        protein_path = self.find_data('Proteins', prot_name)
        self.args.replicon = replicon_path
        self.args.prot_file = protein_path
        self.args.annot_parser = self.find_data('stupid_annot_parser.py')
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        seq_id = 'ACBA.007.P01_13_23'
        db = CustomDB(replicon, cfg, protein_path)
        with self.assertRaises(IntegronError) as ctx:
            with self.catch_log():
                db.get_description(seq_id)
        self.assertEqual(str(ctx.exception),
                         f"Cannot parse protein file '{self.args.prot_file}' with annot-parser '{self.args.annot_parser}': "
                          "cannot unpack non-iterable NoneType object")


    def test_get_description_stupid_parser2(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        protein_path = self.find_data('Proteins', prot_name)
        self.args.replicon = replicon_path
        self.args.prot_file = protein_path
        self.args.annot_parser = self.find_data('stupid_annot_parser2.py')
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        seq_id = 'ACBA.007.P01_13_23'
        db = CustomDB(replicon, cfg, protein_path)
        with self.assertRaises(IntegronError) as ctx:
            with self.catch_log():
                db.get_description(seq_id)
        self.assertEqual(str(ctx.exception),
                         f"'{seq_id}' protein is not compliant with custom --annot-parser '{self.args.annot_parser}'."
                         )


    def test_get_description_lazy_parser(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        protein_path = self.find_data('Proteins', prot_name)
        self.args.replicon = replicon_path
        self.args.prot_file = protein_path
        self.args.annot_parser = self.find_data('lazy_annot_parser.py')
        cfg = Config(self.args)
        seq_db = read_multi_prot_fasta(replicon_path)
        replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        seq_id = 'ACBA.007.P01_13_23'
        db = CustomDB(replicon, cfg, protein_path)
        with self.assertRaises(IntegronError) as ctx:
            with self.catch_log():
                db.get_description(seq_id)
        self.assertEqual(str(ctx.exception),
                         f"Error during protein file parsing: "
                         f"expected seq_id: str, start: positive int, stop: positive int, strand 1/-1. got: "
                         f"ACBA.007.P01_13_23, 19721, 20254, -1")

