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
import tempfile
import shutil
import argparse
import re

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder.utils import FastaIterator
from integron_finder.topology import Topology
from integron_finder.config import Config
from integron_finder import integrase
from integron_finder import EmptyFileError

_run_ori = integrase.subprocess.run


class TestFindIntegrase(IntegronTest):

    def setUp(self):
        if 'INTEGRON_HOME' in os.environ:
            self.integron_home = os.environ['INTEGRON_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.integron_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
        self._tmp_dir = tempfile.TemporaryDirectory(prefix='tmp_test_integron_finder')
        self.tmp_dir = self._tmp_dir.name
        if os.path.exists(self.tmp_dir) and os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)

        self.args = argparse.Namespace()
        self.args.attc_model = 'attc_4.cm'
        self.args.cpu = 1
        self.args.hmmsearch = shutil.which('hmmsearch')
        self.args.prodigal = shutil.which("prodigal")
        self.args.gembase = False
        self.args.prot_file = False
        self.args.cmsearch = __file__
        integrase.subprocess.run = self.mute_call(_run_ori)

    def tearDown(self):
        integrase.subprocess.run = _run_ori
        self._tmp_dir.cleanup()


    def test_find_integrase_gembase(self):
        cfg = Config(self.args)
        self.args.gembase = True
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

        replicon_name = 'acba.007.p01.13'
        replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))

        topologies = Topology(1, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        prot_file = os.path.join(self.tmp_dir, replicon_name + ".prt")

        shutil.copyfile(self.find_data(os.path.join('Proteins', replicon.id + ".prt")), prot_file)

        integrase.find_integrase(replicon.id, prot_file, self.tmp_dir, cfg)

        for suffix in ('_intI.res', '_intI_table.res', '_phage_int.res', '_phage_int_table.res'):
            res = os.path.join(self.tmp_dir, replicon.id + suffix)
            self.assertTrue(os.path.exists(res))


    def test_find_integrase_no_gembase_with_protfile(self):
        try:
            cfg = Config(self.args)
            self.args.gembase = False
            cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

            replicon_name = 'acba.007.p01.13'
            replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
            topologies = Topology(1, 'lin')
            with FastaIterator(replicon_path) as sequences_db:
                sequences_db.topologies = topologies
                replicon = next(sequences_db)

            len_ori = replicon.__class__.__len__
            replicon.__class__.__len__ = lambda x: 200

            prot_file = os.path.join(self.tmp_dir, replicon.id + ".prt")
            shutil.copyfile(self.find_data(os.path.join('Proteins', replicon.id + ".prt")), prot_file)

            integrase.find_integrase(replicon.id, prot_file, self.tmp_dir, cfg)
            for suffix in ('_intI.res', '_intI_table.res', '_phage_int.res', '_phage_int_table.res'):
                res = os.path.join(self.tmp_dir, replicon.id + suffix)
                self.assertTrue(os.path.exists(res))
        finally:
            replicon.__class__.__len__ = len_ori


    def test_find_integrase_no_gembase_with_protfile_empty(self):
        try:
            cfg = Config(self.args)
            self.args.gembase = False
            cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

            replicon_name = 'acba.007.p01.13'
            replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
            topologies = Topology(1, 'lin')
            with FastaIterator(replicon_path) as sequences_db:
                sequences_db.topologies = topologies
                replicon = next(sequences_db)

            len_ori = replicon.__class__.__len__
            replicon.__class__.__len__ = lambda x: 200

            prot_file = os.path.join(self.tmp_dir, replicon.id + ".prt")
            open(prot_file, 'w').close()
            with self.assertRaises(EmptyFileError) as ctx:
                with self.catch_log():
                    integrase.find_integrase(replicon.id, prot_file, self.tmp_dir, cfg)
            self.assertTrue(re.match("^The protein file: '.*' is empty cannot perform hmmsearch on it.$",
                                     str(ctx.exception)))
        finally:
            replicon.__class__.__len__ = len_ori


    def test_find_integrase_no_gembase_no_protfile_short_seq(self):
        try:
            cfg = Config(self.args)
            self.args.gembase = False
            cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

            replicon_name = 'acba.007.p01.13'
            replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
            prot_name = 'ACBA.007.P01_13'
            prot_path = self.find_data(os.path.join('Proteins', prot_name + '.prt'))

            topologies = Topology(1, 'lin')
            with FastaIterator(replicon_path) as sequences_db:
                sequences_db.topologies = topologies
                replicon = next(sequences_db)

            len_ori = replicon.__class__.__len__
            replicon.__class__.__len__ = lambda x: 200

            prot_file = os.path.join(self.tmp_dir, replicon.id + ".prt")
            shutil.copyfile(prot_path, prot_file)

            integrase.find_integrase(replicon.id, prot_file, self.tmp_dir, cfg)
            for suffix in ('_intI.res', '_intI_table.res', '_phage_int.res', '_phage_int_table.res'):
                res = os.path.join(self.tmp_dir, replicon.id + suffix)
                self.assertTrue(os.path.exists(res))
        finally:
            replicon.__class__.__len__ = len_ori


    def test_find_integrase_no_gembase_no_protfile_long_seq(self):
        try:
            cfg = Config(self.args)
            self.args.gembase = False
            cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

            replicon_name = 'acba.007.p01.13'
            replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
            topologies = Topology(1, 'lin')
            with FastaIterator(replicon_path) as sequences_db:
                sequences_db.topologies = topologies
                replicon = next(sequences_db)

            len_ori = replicon.__class__.__len__
            replicon.__class__.__len__ = lambda x: 500000

            prot_file = os.path.join(self.tmp_dir, replicon.id + ".prt")

            shutil.copyfile(self.find_data(os.path.join('Proteins', replicon.id + ".prt")), prot_file)

            integrase.find_integrase(replicon.id, prot_file, self.tmp_dir, cfg)
            for suffix in ('_intI.res', '_intI_table.res', '_phage_int.res', '_phage_int_table.res'):
                res = os.path.join(self.tmp_dir, replicon.id + suffix)
                self.assertTrue(os.path.exists(res))
        finally:
            replicon.__class__.__len__ = len_ori


    def test_find_integrase_no_gembase_no_protfile_no_prodigal(self):
        len_ori = None
        try:
            self.args.hmmsearch = __file__
            self.args.gembase = False
            cfg = Config(self.args)
            cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')
            cfg._args.hmmsearch = 'foo'

            replicon_name = 'acba.007.p01.13'
            replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
            topologies = Topology(1, 'lin')
            with FastaIterator(replicon_path) as sequences_db:
                sequences_db.topologies = topologies
                replicon = next(sequences_db)

            len_ori = replicon.__class__.__len__
            replicon.__class__.__len__ = lambda x: 500000

            prot_file = os.path.join(self.tmp_dir, replicon.id + ".prt")

            shutil.copyfile(self.find_data(os.path.join('Proteins', replicon.id + ".prt")), prot_file)

            with self.assertRaises(RuntimeError) as ctx:
                integrase.find_integrase(replicon.id, prot_file, self.tmp_dir, cfg)
            self.assertTrue(re.search(r"failed : \[Errno 2\] No such file or directory: 'foo'", str(ctx.exception)))
        finally:
            if len_ori:
                replicon.__class__.__len__ = len_ori


    def test_find_integrase_no_gembase_no_protfile(self):
        try:
            cfg = Config(self.args)
            self.args.gembase = False
            cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

            replicon_name = 'acba.007.p01.13'
            replicon_path = self.find_data(os.path.join('Replicons', replicon_name + '.fst'))
            topologies = Topology(1, 'lin')
            with FastaIterator(replicon_path) as sequences_db:
                sequences_db.topologies = topologies
                replicon = next(sequences_db)

            len_ori = replicon.__class__.__len__
            replicon.__class__.__len__ = lambda x: 500000

            prot_file = os.path.join(self.tmp_dir, "foo.prt")
            open(prot_file, 'w').close()
            with self.catch_log():
                with self.assertRaises(EmptyFileError):
                    integrase.find_integrase(replicon.id,  prot_file, self.tmp_dir, cfg)
        finally:
            replicon.__class__.__len__ = len_ori


    def test_find_integrase_gembase_no_hmmer_no_replicon(self):
        self.args.gembase = True
        self.args.hmmsearch = __file__
        cfg = Config(self.args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, 'Replicons', replicon_name + '.fst')
        topologies = Topology(1, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        prot_file = os.path.join(self.tmp_dir, replicon.id + ".prt")

        with self.catch_log():
            with self.assertRaises(RuntimeError) as ctx:
                integrase.find_integrase(replicon.id, prot_file, self.tmp_dir, cfg)
            self.assertEqual("The protein file: '{}' does not exists cannot perform hmmsearch on it.".format(prot_file),
                             str(ctx.exception))


    def test_find_integrase_gembase_hmmer_error(self):
        self.args.gembase = True
        self.args.cpu = 'foo'
        cfg = Config(self.args)
        cfg._prefix_data = os.path.join(os.path.dirname(__file__), 'data')

        replicon_name = 'acba.007.p01.13'
        replicon_path = os.path.join(self._data_dir, 'Replicons', replicon_name + '.fst')
        topologies = Topology(1, 'lin')
        with FastaIterator(replicon_path) as sequences_db:
            sequences_db.topologies = topologies
            replicon = next(sequences_db)

        prot_file = os.path.join(self.tmp_dir, replicon.id + ".prt")
        shutil.copyfile(os.path.join(self._data_dir, 'Proteins', replicon.id + ".prt"),
                        prot_file)
        with self.assertRaises(RuntimeError) as ctx:
            integrase.find_integrase(replicon.id, prot_file, self.tmp_dir, cfg)
        self.assertTrue(str(ctx.exception).endswith('failed return code = 1'))
