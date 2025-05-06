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
import argparse
import tempfile
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
from integron_finder.utils import MultiFastaReader
from integron_finder.prot_db import GembaseDB, ProdigalDB, SeqDesc, CustomDB, GembaseType, RepliconType

class TestGemBase(IntegronTest):

    def setUp(self):
        """
        Define variables common to all tests
        """
        # Simulate argparse to get argument
        self.args = argparse.Namespace()
        self.args.gembase = True
        self.args.prot_file = False
        self._tmp_dir = tempfile.TemporaryDirectory(prefix='tmp_test_integron_finder')
        self.tmp_dir = self._tmp_dir.name
        self.args.outdir = self.tmp_dir
        self.args.prodigal = None
        self.args.cmsearch = __file__
        self.args.hmmsearch = __file__
        if os.path.exists(self.tmp_dir) and os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)
        self.maxDiff = 30


    def tearDown(self):
        self._tmp_dir.cleanup()


    def test_ProteinDB(self):
        # From Gembase1 Draft , Gembase1 Complete
        file_names = (
            (1, 'ACBA.0917.00019.fna'),
            (1, 'ESCO001.C.00001.C001.fst'),
            (2, 'VIBR.0322.11443.fna'),
            (2, 'VICH001.0523.00090.fna')
        )
        for gb_v, file_name in file_names:
            replicon_path = self.find_data(os.path.join('Gembase', f'Gembase{gb_v}', 'Replicons', file_name))
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            seq_db = MultiFastaReader(replicon_path)
            replicon = next(seq_db)
            seq_db.close()
            replicon.path = replicon_path
            os.makedirs(cfg.tmp_dir(replicon.id))
            try:
                with self.catch_log():
                    db = GembaseDB(replicon, cfg)
                self.assertTrue(db.replicon.id, replicon.id)
            finally:
                db.close()

    def test_no_LST_dir(self):
        dest_replicons_path = os.path.join(self.tmp_dir, 'Gembase', 'Replicons')
        os.makedirs(dest_replicons_path)
        file_name = 'VICH001.0523.00090.fna'
        replicon_path = self.find_data(os.path.join('Gembase', 'Gembase2', 'Replicons', file_name))
        replicon_copy = os.path.join(dest_replicons_path, file_name)
        shutil.copy(replicon_path, replicon_copy)
        self.args.replicon = replicon_copy
        cfg = Config(self.args)
        with MultiFastaReader(replicon_copy) as seq_db:
            replicon = next(seq_db)
        with self.assertRaises(IntegronError) as ctx:
            GembaseDB(replicon, cfg)
        self.assertEqual(str(ctx.exception),
                         f"Neither 'LST' nor 'LSTINF' nor 'LSTINFO' directory found in '{os.path.dirname(dest_replicons_path)}' .")


    def test_find_gembase_file_basename(self):
        # test if *find_gembase_file_basename* get the right basename
        # for files in gembase
        file_names = (
            ('1', 'ACBA.0917.00019.fna', 'ACBA.0917.00019'),                    # Draft v1
            ('1', 'ESCO001.C.00001.C001.fst', 'ESCO001.C.00001.C001'),          # Complete v1
            ('2', 'VICH001.0523.00090.fna', 'VICH001.0523.00090'),              # Complete V2
            ('2', 'VIBR.0322.11443.fna', 'VIBR.0322.11443'),                    # Draft V2
            ('2plus', 'VICH001.0523.00090.fna', 'VICH001.0523.00090'),          # Complete V2plus
            ('2plus', 'VIBR001.0322.11443.fna', 'VIBR001.0322.11443'),          # Draft V2plus
            ('1', 'ESCO001.C.00001.C001_chunk_1.fst' , 'ESCO001.C.00001.C001'), # Complete v1 chunk
        )
        for gbv, rep_file_name, basename in file_names:
            gembase_path = self.find_data('Gembase', f'Gembase{gbv}')
            replicon_path = self.find_data('Gembase', f'Gembase{gbv}', 'Replicons', rep_file_name)
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            with MultiFastaReader(replicon_path) as seq_db:
                replicon = next(seq_db)
            replicon.path = replicon_path

            with self.subTest(file_name=rep_file_name):
                with self.catch_log():
                    # if sequence present in lst not present in prt file
                    # a warning is raised
                    db = GembaseDB(replicon, cfg)
                self.assertEqual(db.find_gembase_file_basename(gembase_path, replicon_path),
                                 basename)


    def test_find_gembase_file_basename_file_not_in_gembase(self):
        # test if find_gembase_file_basename get the right basename
        # for files not located in gembase and file name is the output of split operation
        # a file containing one contig
        # a file representing a chunk
        gembase_path_ori = self.find_data('Gembase', 'Gembase2')
        gembase_path_dest = os.path.join(self.tmp_dir, 'Gembase')
        shutil.copytree(gembase_path_ori, gembase_path_dest)

        extra_replicon_ori = self.find_data('Gembase', 'Gembase1', 'Replicons', 'ACBA.0917.00019.fna')
        extra_replicon_path = os.path.join(gembase_path_dest, 'Replicons', os.path.basename(extra_replicon_ori))
        shutil.copy(extra_replicon_ori, extra_replicon_path)

        self.args.replicon = extra_replicon_path
        self.args.gembase_path = gembase_path_dest
        cfg = Config(self.args)
        with MultiFastaReader(extra_replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = extra_replicon_path

        with self.assertRaises(FileNotFoundError) as ctx:
            with self.catch_log():
                GembaseDB(replicon, cfg, gembase_path=gembase_path_dest)
        self.assertEqual(str(ctx.exception),
                         f'cannot find lst file matching {extra_replicon_path} sequence')


    def test_gembase_sniffer(self):
        file_names = (
                      (1, 'ACBA.0917.00019', GembaseType.DRAFT_1),
                      (1, 'ESCO001.C.00001.C001', GembaseType.COMPLETE_1),
                      (1, 'ACJO001.0321.00008.P008', GembaseType.COMPLETE_1),
                      (2, 'VICH001.0523.00090', GembaseType.COMPLETE_2),
                      (2, 'VIBR.0322.11443', GembaseType.DRAFT_2),
                      )
        for gbv, file_name, gem_type in file_names:
            lst_dir = 'LSTINF' if gbv == 1 else 'LST'
            lst_path = self.find_data(os.path.join('Gembase', f'Gembase{gbv}', lst_dir, file_name + '.lst'))
            with self.subTest(lst_path=lst_path):
                type_recieved = GembaseDB.gembase_sniffer(lst_path)
                self.assertEqual(type_recieved, gem_type)

        lst_path = self.find_data(os.path.join('Gembase', 'Gembase1', 'LSTINF', 'SAEN001.0321.00753.P003.lst'))
        with self.assertRaises(IntegronError) as ctx:
            with self.catch_log():
                GembaseDB.gembase_sniffer(lst_path)
        self.maxDiff = None
        self.assertEqual(str(ctx.exception),
                         f"The genome SAEN001.0321.00753.P003_00001 seems empty: see {lst_path}")

        lst_path = self.find_data(os.path.join('Gembase', 'Gembase1', 'LSTINF', 'Bad_Fmt.0321.00008.P008.lst'))
        with self.assertRaises(IntegronError) as ctx:
            with self.catch_log():
                GembaseDB.gembase_sniffer(lst_path)
        self.maxDiff = None
        self.assertEqual(str(ctx.exception),
                         f"Cannot detect GemBase version, check lst file '{lst_path}'.")

    def test_gembase1_complete_parser(self):
        replicon_id = 'ESCO001.C.00001.C001'
        lst_path = self.find_data(os.path.join('Gembase', 'Gembase1', 'LSTINF', replicon_id + '.lst'))
        prots_info = GembaseDB.gembase1_complete_parser(lst_path, replicon_id)
        columns = list(range(8))
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


    def test_gembase1_draft_parser(self):
        replicon_name = 'ACBA.0917.00019'
        replicon_id = 'ACBA.0917.00019.0001'
        lst_path = self.find_data('Gembase', 'Gembase1', 'LSTINF', replicon_name + '.lst')
        prots_info = GembaseDB.gembase1_draft_parser(lst_path, replicon_id)
        columns = list(range(7))
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

        # check thet if it's not right format it raise an error
        replicon_id = 'ESCO001.C.00001.C001'
        lst_path = self.find_data('Gembase', 'Gembase1', 'LSTINF', f'{replicon_id}.lst')
        with self.assertRaises(IntegronError) as ctx:
            with self.catch_log() as log:
                GembaseDB.gembase1_draft_parser(lst_path, replicon_id)
                log_msg = log.get_value().strip()
            self.assertEqual(log_msg,
                             f"The LST file '{lst_path}' seems not to be in gembase V1 draft format.")
        self.assertEqual(str(ctx.exception),
                         f"The LST file '{lst_path}' seems not to be in gembase V1 draft format.")

        replicon_id = 'SAEN001.0321.00753.P003'
        lst_path = os.path.join(self.tmp_dir, f'{replicon_id}.fst')
        open(lst_path, 'w').close()
        with self.assertRaises(IntegronError) as ctx:
            with self.catch_log() as log:
                GembaseDB.gembase1_draft_parser(lst_path, replicon_id)
                log_msg = log.get_value().strip()
            self.assertEqual(log_msg,
                             f" Error while parsing {lst_path} file: No columns to parse from file")


    def test_gembase2_complete_parser(self):
        replicon_name = 'VICH001.0523.00090'
        replicon_id = 'VICH001.0523.00090.001C'
        lst_path = self.find_data('Gembase', 'Gembase2', 'LST', replicon_name + '.lst')
        prots_info = GembaseDB.gembase2_parser(lst_path, replicon_id)
        self.assertEqual(prots_info.shape, (2570, 8))
        first_row = [1, 2961008, 'C', 'CDS', 'VICH001.0523.00090.001C_00001',
                     'FKV26_RS00005', 'mltC', 'WP_001230176.1',
                     ]
        recieved_first_row = prots_info.iloc[0].values.tolist()
        self.assertListEqual(first_row, recieved_first_row)

        last_row = [2958304, 2959929, 'C', 'CDS', 'VICH001.0523.00090.001C_02700',
                    'FKV26_RS13450', 'nan', 'WP_000769916.1' ]
        recieved_last_row = prots_info.iloc[len(prots_info) - 1].values.tolist()
        self.assertListEqual(last_row, recieved_last_row)

        replicon_id = 'SAEN001.0321.00753.P003'
        lst_path = os.path.join(self.tmp_dir, f'{replicon_id}.fst')
        open(lst_path, 'w').close()
        with self.assertRaises(IntegronError) as ctx:
            with self.catch_log() as log:
                GembaseDB.gembase2_parser(lst_path, replicon_id)
                log_msg = log.get_value().strip()
            self.assertEqual(log_msg,
                             f" Error while parsing {lst_path} file: No columns to parse from file")
        self.assertEqual(str(ctx.exception),
                         f"Error while parsing {lst_path} file: No columns to parse from file")

        replicon_id = 'SAEN001.0321.00753.P003'
        lst_path = self.find_data('Gembase', 'Gembase1', 'LSTINF', f'{replicon_id}.lst')
        with self.assertRaises(IntegronError) as ctx:
            with self.catch_log() as log:
                GembaseDB.gembase2_parser(lst_path, replicon_id)
                log_msg = log.get_value().strip()
            self.assertEqual(log_msg,
                             f" The LST file '{lst_path}' seems not to be in gembase V2 draft format.")
        self.assertEqual(str(ctx.exception), f"The LST file '{lst_path}' seems not to be in gembase V2 draft format.")


    def test_gembase2_draft_parser(self):
        replicon_name = 'VIBR.0322.11443'
        replicon_id = 'VIBR.0322.11443.0001'
        lst_path = self.find_data(os.path.join('Gembase', 'Gembase2', 'LST', replicon_name + '.lst'))
        prots_info = GembaseDB.gembase2_parser(lst_path, replicon_id)
        self.assertEqual(prots_info.shape, (114, 7))
        recieved_first_row = prots_info.iloc[0].values.tolist()[:-1]
        first_row = [2, 655, 'C', 'CDS', 'VIBR.0322.11443.0001b_00001', 'nan']
        recieved_last_row = prots_info.iloc[-1].values.tolist()[:-1]
        last_row = [138318, 140534, 'C', 'CDS', 'VIBR.0322.11443.0001b_00114', 'nan']
        self.assertListEqual(first_row, recieved_first_row)
        self.assertListEqual(last_row, recieved_last_row)


    def test_get_replicon_type(self):
        expected_types= (
            ('VIBR.0322.11443.0002', RepliconType.DRAFT),
            ('VICH001.0523.00090.001C', RepliconType.CHROMOSOME),
            ('ACBA.0917.00019.0001', RepliconType.DRAFT),
            ('ACJO001.0321.00008.P008', RepliconType.PLASMID),
            ('ESCO001.C.00001.C001', RepliconType.CHROMOSOME)
        )

        for seq_id, exp_typ in expected_types:
            with self.subTest(seq_id=seq_id):
                get_typ = GembaseDB.get_replicon_type(seq_id=seq_id)
                self.assertEqual(get_typ, exp_typ)

        with self.assertRaises(IntegronError) as ctx:
            seq_id = "foo"
            with self.catch_log() as log:
                GembaseDB.get_replicon_type(seq_id=seq_id)
                log_msg = log.get_value().strip()
            exp_msg = f"Cannot detect GemBase version, from seq_id: {seq_id}."
            self.assertEqual(log_msg,
                             exp_msg
                             )
            self.assertEqual(str(ctx.exception),
                             exp_msg
                             )

        expected_types =(
            ('VIBR.0322.11443.0149b_03432', RepliconType.DRAFT),
            ('VICH001.0523.00090.001C_00024', RepliconType.CHROMOSOME),
            ('ACBA.0917.00019.i0001_00026', RepliconType.DRAFT),
            ('ACJO001.0321.00008.P008_00009', RepliconType.PLASMID),
            ('ESCO001.C.00001.C001_00021', RepliconType.CHROMOSOME)
        )

        for gene_id, exp_typ in expected_types:
            with self.subTest(gene_id=gene_id):
                get_typ = GembaseDB.get_replicon_type(rep_id=gene_id)
                self.assertEqual(get_typ, exp_typ)

        with self.assertRaises(IntegronError) as ctx:
            GembaseDB.get_replicon_type()
            self.assertEqual(str(ctx.exception),
                             'GembaseDB.get_replicon_type you must provide either a seqid or a rep_id'
                             )

        with self.assertRaises(IntegronError) as ctx:
            GembaseDB.get_replicon_type(seq_id='foo', rep_id='bar')
            self.assertEqual(str(ctx.exception),
                             'GembaseDB.get_replicon_type you must provide either a seqid or a rep_id'
                             )


    def test_make_protfile(self):
        file_name = (('ACBA.0917.00019', '.fna', 3870), ('ESCO001.C.00001.C001', '.fst', 3870))
        for seq_name, ext, seq_nb in file_name:
            replicon_path = self.find_data(os.path.join('Gembase', 'Gembase1', 'Replicons', seq_name + ext))
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            with MultiFastaReader(replicon_path) as seq_db:
                replicon = next(seq_db)
            replicon.path = replicon_path
            os.makedirs(cfg.tmp_dir(replicon.id))

            with self.catch_log():
                db = GembaseDB(replicon, cfg)
            for seq_nb, seqs in enumerate(zip(
                    MultiFastaReader(self.find_data(os.path.join('Gembase', 'Gembase1', 'Proteins', seq_name + '.prt'))),
                    MultiFastaReader(db.protfile)),
                    1):
                expected, test = seqs
                self.assertEqual(expected.id, test.id)
            self.assertEqual(seq_nb, seq_nb)


    def test_protfile(self):
        file_name = (('ACBA.0917.00019', '.fna'), ('ESCO001.C.00001.C001', '.fst'))
        for seq_name, ext in file_name:
            replicon_path = self.find_data(os.path.join('Gembase', 'Gembase1', 'Replicons', seq_name + ext))
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            with MultiFastaReader(replicon_path) as seq_db:
                replicon = next(seq_db)
            replicon.path = replicon_path
            os.makedirs(cfg.tmp_dir(replicon.id))

            with self.catch_log():
                db = GembaseDB(replicon, cfg)
            self.assertEqual(os.path.join(cfg.tmp_dir(replicon.id), replicon.id + '.prt'), db.protfile)


    def test_getitem(self):
        file_name = (('ACBA.0917.00019', '.fna'), ('ESCO001.C.00001.C001', '.fst'))
        for seq_name, ext in file_name:
            replicon_path = self.find_data(os.path.join('Gembase', 'Gembase1', 'Replicons', seq_name + ext))
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            with MultiFastaReader(replicon_path) as seq_db:
                replicon = next(seq_db)
            os.makedirs(cfg.tmp_dir(replicon.id))

            with self.catch_log():
                db = GembaseDB(replicon, cfg)
            with MultiFastaReader(self.find_data('Gembase', 'Gembase1', 'Proteins', seq_name + '.prt')) as exp:
                specie, date, strain, contig = replicon.id.split('.')
                pattern = rf'{specie}\.{date}\.{strain}\.\w?{contig}'

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
        replicon_path = self.find_data(os.path.join('Gembase', 'Gembase1', 'Replicons', seq_name + ext))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = replicon_path
        try:
            with self.catch_log():
                db = GembaseDB(replicon, cfg)

            idx = SeqIO.index(self.find_data(os.path.join('Gembase', 'Gembase1', 'Proteins', seq_name + '.prt')), 'fasta')
            specie, date, strain, contig = replicon.id.split('.')
            pattern = fr'{specie}\.{date}\.{strain}\.\w?{contig}'
            self.assertListEqual(sorted([i for i in idx if re.match(pattern, i)]), sorted([i for i in db]))
        finally:
            db.close()
            idx.close()
        # test Gembase Complet
        seq_name = 'ESCO001.C.00001.C001'
        ext = '.fst'
        replicon_path = self.find_data(os.path.join('Gembase', 'Gembase1', 'Replicons', seq_name + ext))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = replicon_path
        try:
            with self.catch_log():
                db = GembaseDB(replicon, cfg)
            idx = SeqIO.index(self.find_data(os.path.join('Gembase', 'Gembase1', 'Proteins', seq_name + '.prt')), 'fasta')
            specie, date, strain, contig = replicon.id.split('.')
            pattern = fr'{specie}\.{date}\.{strain}\.\w?{contig}'
            seqid_from_gembase_protfile = set([i for i in idx if re.match(pattern, i)])
            seqid_from_if = set([i for i in db])
            non_common_seq = seqid_from_gembase_protfile ^ seqid_from_if
            # in Gembase complete the annotation from lstinfo provided from genbank
            # it appear some times that some CDS are not translate in proteins
            # So in data I have 3 genes from LSTINF are not in .prt file
            diff = {'ESCO001.C.00001.C001_03974', 'ESCO001.C.00001.C001_01509', 'ESCO001.C.00001.C001_04162'}
            self.assertSetEqual(non_common_seq, diff)
        finally:
            db.close()
            idx.close()

    def test_codig_prot_ids(self):
        seq_name = 'ESCO001.C.00001.C001'
        ext = '.fst'
        replicon_path = self.find_data(os.path.join('Gembase', 'Gembase1', 'Replicons', seq_name + ext))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        seq_db.close()
        replicon.path = replicon_path
        try:
            with self.catch_log():
                db = GembaseDB(replicon, cfg)
            idx = SeqIO.index(self.find_data(os.path.join('Gembase', 'Gembase1', 'Proteins', seq_name + '.prt')), 'fasta')

            specie, date, strain, contig = replicon.id.split('.')
            pattern = fr'{specie}\.{date}\.{strain}\.\w?{contig}'
            seqid_from_gembase_protfile = set([i for i in idx if re.match(pattern, i)])
            seqid_from_if = set([i for i in db.coding_prot_ids()])
            self.assertSetEqual(seqid_from_if, seqid_from_gembase_protfile)
        finally:
            db.close()
            idx.close()


    def test_is_pseudo_gene(self):
        seq_name = 'ESCO001.C.00001.C001'
        ext = '.fst'
        replicon_path = self.find_data(os.path.join('Gembase', 'Gembase1', 'Replicons', seq_name + ext))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = replicon_path
        with self.catch_log():
            db = GembaseDB(replicon, cfg)
            self.assertTrue(db.is_pseudo_gene('ESCO001.C.00001.C001_03974'))
            self.assertFalse(db.is_pseudo_gene('ESCO001.C.00001.C001_03975'))


    def test_get_description(self):
        # SeqDesc(id, strand, statt, stop)
        file_name = {('ACBA.0917.00019', '.fna'):
                         {'ACBA.0917.00019.b0001_00001': SeqDesc('ACBA.0917.00019.b0001_00001', -1, 266, 1480),
                          'ACBA.0917.00019.i0001_03957': SeqDesc('ACBA.0917.00019.i0001_03957', -1, 4043755, 4044354)},
                     }

        for seq_name, ext in file_name:
            replicon_path = self.find_data(os.path.join('Gembase', 'Gembase1', 'Replicons', seq_name + ext))
            self.args.replicon = replicon_path
            cfg = Config(self.args)
            with MultiFastaReader(replicon_path) as seq_db:
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
        self.args.prot_file = False
        self._tmp_dir = tempfile.TemporaryDirectory(prefix='tmp_test_integron_finder')
        self.tmp_dir = self._tmp_dir.name
        self.args.outdir = self.tmp_dir
        self.args.prodigal = shutil.which("prodigal")
        self.args.cmsearch = __file__
        self.args.hmmsearch = __file__
        if os.path.exists(self.tmp_dir) and os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)


    def tearDown(self):
        self._tmp_dir.cleanup()


    def test_ProteinDB(self):
        file_name = 'acba.007.p01.13'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        with MultiFastaReader(replicon_path) as seq_db:
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
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        self.args.prodigal = None
        with self.assertRaises(AssertionError) as ctx:
            ProdigalDB(replicon, cfg)
        self.assertEqual(str(ctx.exception),
                         "'prodigal' not found.")

    def test_make_protfile(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'

        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        del seq_db
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        try:
            db = ProdigalDB(replicon, cfg)
            for seq_nb, seqs in enumerate(zip(
                    MultiFastaReader(self.find_data(os.path.join('Proteins', prot_name))),
                    MultiFastaReader(db.protfile)), 1):
                expected, test = seqs
                self.assertEqual(expected.id, test.id)
            self.assertEqual(seq_nb, 23)
        finally:
            db.close()

    def test_make_protfile_no_dir(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = replicon_path

        db = ProdigalDB(replicon, cfg)
        for seq_nb, seqs in enumerate(zip(
                MultiFastaReader(self.find_data(os.path.join('Proteins', prot_name))),
                MultiFastaReader(db.protfile)), 1):
            expected, test = seqs
            self.assertEqual(expected.id, test.id)
        self.assertEqual(seq_nb, 23)


    def test_make_protfile_prodigal_failed(self):
        file_name = 'acba.007.p01.13'
        replicon_path = self.find_data('Replicons', file_name + '.fst')
        self.args.replicon = replicon_path
        self.args.prodigal = self.find_data('fake_prodigal')
        cfg = Config(self.args)
        with MultiFastaReader(replicon_path) as seq_db:
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
        with MultiFastaReader(replicon_path) as seq_db:
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
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = ProdigalDB(replicon, cfg)
        with MultiFastaReader(self.find_data('Proteins', prot_name)) as exp:
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
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = ProdigalDB(replicon, cfg)

        idx = SeqIO.index(self.find_data(os.path.join('Proteins', prot_name)), 'fasta')
        try:
            for exp_seq_id, get_seq_id in zip(idx, db):
                self.assertEqual(exp_seq_id, get_seq_id)
        finally:
            idx.close()

    def test_coding_prot_ids(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = ProdigalDB(replicon, cfg)

        idx = SeqIO.index(self.find_data(os.path.join('Proteins', prot_name)), 'fasta')
        # there is no pseudo genes in prodigalDB so iterating on coding_prot_ids == iteratiing on whole DB
        try:
            for exp_seq_id, get_seq_id in zip(idx, db.coding_prot_ids()):
                self.assertEqual(exp_seq_id, get_seq_id)
        finally:
            idx.close()

    def test_get_description(self):
        # SeqDesc(id, strand, strat, stop)
        file_name = 'acba.007.p01.13'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        self.args.replicon = replicon_path
        cfg = Config(self.args)
        with MultiFastaReader(replicon_path) as seq_db:
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
        self.args.prot_file = True

        self.args.prodigal = None
        self.args.cmsearch = __file__
        self.args.hmmsearch = __file__
        self._tmp_dir = tempfile.TemporaryDirectory(prefix='tmp_test_integron_finder')
        self.tmp_dir = self._tmp_dir.name
        self.args.outdir = self.tmp_dir
        self.args.annot_parser = self.find_data('prodigal_annot_parser.py')

        if os.path.exists(self.tmp_dir) and os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)


    def tearDown(self):
        self._tmp_dir.cleanup()


    def test_ProteinDB(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data('Replicons', file_name + '.fst')
        protein_path = self.find_data('Proteins', prot_name)
        self.args.replicon = replicon_path
        self.args.prot_file = protein_path
        cfg = Config(self.args)
        with MultiFastaReader(replicon_path) as seq_db:
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
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        with self.assertRaises(RuntimeError) as ctx:
            with self.catch_log():
                CustomDB(replicon, cfg, protein_path)
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
        with MultiFastaReader(replicon_path) as seq_db:
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
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = CustomDB(replicon, cfg, protein_path)
        exp = MultiFastaReader(self.find_data(os.path.join('Proteins', prot_name)))
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
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))
        try:
            db = CustomDB(replicon, cfg, protein_path)
            idx = SeqIO.index(self.find_data(os.path.join('Proteins', prot_name)), 'fasta')
            for exp_seq_id, get_seq_id in zip(idx, db):
                self.assertEqual(exp_seq_id, get_seq_id)
        finally:
            db.close()
            idx.close()

    def test_coding_prot_ids(self):
        file_name = 'acba.007.p01.13'
        prot_name = 'ACBA.007.P01_13.prt'
        replicon_path = self.find_data(os.path.join('Replicons', file_name + '.fst'))
        protein_path = self.find_data('Proteins', prot_name)
        self.args.replicon = replicon_path
        self.args.prot_file = protein_path
        cfg = Config(self.args)
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        seq_db.close()
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        db = CustomDB(replicon, cfg, protein_path)
        try:
            idx = SeqIO.index(self.find_data(os.path.join('Proteins', prot_name)), 'fasta',
                              alphabet=Seq.IUPAC.extended_protein)
        except AttributeError:
            idx = SeqIO.index(self.find_data(os.path.join('Proteins', prot_name)), 'fasta')
        finally:
            idx.close()
        # there is no pseudgenes in customdb
        for exp_seq_id, get_seq_id in zip(idx, db.coding_prot_ids()):
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
        with MultiFastaReader(replicon_path) as seq_db:
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
        with MultiFastaReader(replicon_path) as seq_db:
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
        with MultiFastaReader(replicon_path) as seq_db:
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
        with MultiFastaReader(replicon_path) as seq_db:
            replicon = next(seq_db)
        replicon.path = replicon_path
        os.makedirs(cfg.tmp_dir(replicon.id))

        seq_id = 'ACBA.007.P01_13_23'
        db = CustomDB(replicon, cfg, protein_path)
        with self.assertRaises(IntegronError) as ctx:
            with self.catch_log():
                db.get_description(seq_id)
        self.assertEqual(str(ctx.exception),
                         "Error during protein file parsing: "
                         "expected seq_id: str, start: positive int, stop: positive int, strand 1/-1. got: "
                         "ACBA.007.P01_13_23, 19721, 20254, -1")


class TestRepliconType(IntegronTest):


    def test_str(self):
        self.assertEqual(str(RepliconType.CHROMOSOME), 'Chromosome')
        self.assertEqual(str(RepliconType.OTHER), 'Other')

    def topology(self):
        topos = (
            (RepliconType.CHROMOSOME, 'circ'),
            (RepliconType.PLASMID, 'circ'),
            (RepliconType.PHAGE, 'circ'),
            (RepliconType.OTHER, 'circ'),
            (RepliconType.DRAFT, 'lin'),
        )

        for r_type, topo in topos:
            with self.subTest(topology=topo):
                self.assertEqual(r_type.topology(), topo)


class TestGembaseType(IntegronTest):

    def test_str(self):
        self.assertEqual(str(GembaseType.COMPLETE_2), 'vers: 2 Complete')
        self.assertEqual(str(GembaseType.DRAFT_1), 'vers: 1 Draft')
        self.assertEqual(str(GembaseType.DRAFT_2plus), 'vers: 2plus Draft')


    def test_complete(self):
        for gt in (GembaseType.COMPLETE_2, GembaseType.COMPLETE_1, GembaseType.COMPLETE_2plus):
            with self.subTest(gt=str(gt)):
                self.assertTrue(gt.complete)

        for gt in (GembaseType.DRAFT_1, GembaseType.DRAFT_2, GembaseType.DRAFT_2plus):
            with self.subTest(gt=str(gt)):
                self.assertFalse(gt.complete)

    def test_version(self):

        for gt , v in (
                (GembaseType.DRAFT_1, '1'),
                (GembaseType.COMPLETE_1, '1'),
                (GembaseType.DRAFT_2, '2'),
                (GembaseType.DRAFT_2, '2'),
                (GembaseType.COMPLETE_2, '2'),
                (GembaseType.DRAFT_2plus, '2plus'),
                (GembaseType.COMPLETE_2plus, '2plus')
        ):
            with self.subTest(gt=str(gt)):
                self.assertEqual(gt.version, v)
