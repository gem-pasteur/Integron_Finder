# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright (c) 2015 - 2025  Institut Pasteur, Paris and CNRS                      #
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

import argparse
import os

try:
    from tests import IntegronTest
except ImportError as err:
    msg = "Cannot import integron_finder: {0!s}".format(err)
    raise ImportError(msg)

from integron_finder import config


class TestConfig(IntegronTest):

    def setUp(self):
        self.args = argparse.Namespace()
        self.args.gembase = False
        self.args.prot_file = False
        self.args.prodigal = __file__
        self.args.cmsearch = __file__
        self.args.hmmsearch = __file__


    def test_config(self):
        # check that we can instanciate a config
        # if gembase is True even prodigal is not provided
        self.args = argparse.Namespace()
        self.args.gembase = True
        self.args.prot_file = False
        self.args.prodigal = None
        self.args.cmsearch = __file__
        self.args.hmmsearch = __file__
        cf = config.Config(self.args)

        # check that we cannot instanciate a config
        # if cmsearch is not provided
        self.args = argparse.Namespace()
        self.args.gembase = True
        self.args.prot_file = False
        self.args.prodigal = None
        self.args.cmsearch = None
        self.args.hmmsearch = __file__
        with self.assertRaises(RuntimeError) as ctx:
            with self.catch_log():
                config.Config(self.args)
        self.assertEqual(str(ctx.exception),
                         "Cannot find 'cmsearch' in PATH.\n" \
                          "Please install infernal package or setup 'cmsearch' binary path with --cmsearch option"
                         )
        # check that we cannot instanciate a config
        # if hmmsearch is not provided
        self.args = argparse.Namespace()
        self.args.gembase = True
        self.args.prot_file = False
        self.args.prodigal = None
        self.args.cmsearch = __file__
        self.args.hmmsearch = None
        with self.assertRaises(RuntimeError) as ctx:
            with self.catch_log():
                config.Config(self.args)
        self.assertEqual(str(ctx.exception),
                         "Cannot find 'hmmsearch' in PATH.\n" \
                          "Please install hmmer package or setup 'hmmsearch' binary path with --hmmsearch option"
                         )
        # check that we can instanciate a config
        # if a custom prot_file is provided even prodigal is not provided
        self.args = argparse.Namespace()
        self.args.gembase = False
        self.args.prot_file = True
        self.args.prodigal = None
        self.args.cmsearch = __file__
        self.args.hmmsearch = __file__
        cf = config.Config(self.args)
        self.assertTrue(isinstance(cf, config.Config))

        # check that we can instanciate a config
        # if a prodigal, hmmsearch and cmsearch are provided
        self.args = argparse.Namespace()
        self.args.gembase = False
        self.args.prot_file = False
        self.args.prodigal = __file__
        self.args.cmsearch = __file__
        self.args.hmmsearch = __file__
        cf = config.Config(self.args)
        self.assertTrue(isinstance(cf, config.Config))

        # check that we cannot instanciate a config
        # if prodigal is not provided and it's not a gambase nor a custom prot_file
        self.args = argparse.Namespace()
        self.args.gembase = False
        self.args.prot_file = False
        self.args.prodigal = None
        self.args.cmsearch = __file__
        self.args.hmmsearch = __file__
        with self.assertRaises(RuntimeError) as ctx:
            with self.catch_log():
                config.Config(self.args)
        self.assertEqual(str(ctx.exception),
                         "Cannot find 'prodigal' in PATH.\n" \
                          "Please install Prodigal package or setup 'prodigal' binary path with --prodigal option"
                         )

    def test_getattr(self):
        self.args.replicon = 'foo'
        self.args.bar = 'baz'
        cf = config.Config(self.args)
        self.assertEqual(cf.input_seq_path, os.path.abspath('foo'))
        self.assertEqual(cf.bar, 'baz')
        with self.assertRaises(AttributeError) as ctx:
            self.assertEqual(cf.foobar, 'foobar')
        self.assertEqual(str(ctx.exception), "config object has no attribute 'foobar'")

    def test_replicon_path(self):
        self.args.replicon = '../foo'
        cf = config.Config(self.args)
        self.assertEqual(cf.input_seq_path, os.path.abspath('../foo'))

    def test_input_dir(self):
        self.args.replicon = '../foo.fasta'
        cf = config.Config(self.args)
        self.assertEqual(cf.input_dir, os.path.split(os.path.abspath('../foo'))[0])

    def test_result_dir(self):
        replicon = '../foo.fasta'
        outdir = 'outdir'
        self.args.replicon = replicon
        self.args.outdir = outdir
        cf = config.Config(self.args)
        exp_result_dir = os.path.abspath(os.path.join(outdir,
                                                      "Results_Integron_Finder_" +
                                                      os.path.splitext(os.path.split(replicon)[1])[0]))
        self.assertEqual(cf.result_dir, exp_result_dir)

    def test_tmp_dir(self):
        replicon_id = 'foo'
        replicon_path = '../{}.fasta'.format(replicon_id)
        outdir = 'outdir'
        self.args.replicon = replicon_path
        self.args.outdir = outdir
        cf = config.Config(self.args)
        exp_result_tmp_dir = os.path.abspath(os.path.join(outdir,
                                                          "Results_Integron_Finder_{}".format(replicon_id),
                                                          "tmp_{}".format(replicon_id))
                                             )
        self.assertEqual(cf.tmp_dir(replicon_id), exp_result_tmp_dir)

    def test_default_topology(self):
        cf = config.Config(self.args)
        self.assertIsNone(cf.default_topology)
        self.args.circular = True
        self.args.linear = False
        cf = config.Config(self.args)
        self.assertEqual(cf.default_topology, 'circ')
        self.args.circular = False
        self.args.linear = True
        cf = config.Config(self.args)
        self.assertEqual(cf.default_topology, 'lin')
        self.args.circular = False
        self.args.linear = False
        cf = config.Config(self.args)
        self.assertIsNone(cf.default_topology)

    def test_model_dir(self):
        cf = config.Config(self.args)
        cf._prefix_data = 'foo'
        self.assertEqual(cf.model_dir, os.path.join('foo', 'Models'))

    def test_model_integrase(self):
        cf = config.Config(self.args)
        cf._prefix_data = 'foo'
        self.assertEqual(cf.model_integrase, os.path.join('foo', 'Models', 'integron_integrase.hmm'))

    def test_model_phage_int(self):
        cf = config.Config(self.args)
        cf._prefix_data = 'foo'
        self.assertEqual(cf.model_phage_int, os.path.join('foo', 'Models', 'phage-int.hmm'))

    def test_model_attc_path(self):
        cf = config.Config(self.args)
        with self.assertRaises(RuntimeError) as ctx:
            cf.model_attc_path
        self.assertEqual(str(ctx.exception), "'model_attc' is not define.")
        
        self.args.attc_model = 'bar'
        cf = config.Config(self.args)
        cf._prefix_data = 'foo'
        self.assertEqual(cf.model_attc_path, os.path.join('foo', 'Models', 'bar'))
        self.args.attc_model = 'bar/baz'
        self.assertEqual(cf.model_attc_path, 'bar/baz')

    def test_model_attc_name(self):
        cf = config.Config(self.args)
        with self.assertRaises(RuntimeError) as ctx:
            cf.model_attc_name
        self.assertEqual(str(ctx.exception), "'model_attc' is not define.")

        self.args.attc_model = 'bar'
        cf = config.Config(self.args)
        cf._prefix_data = 'foo'
        self.assertEqual(cf.model_attc_name, 'bar')
        self.args.attc_model = 'bar/baz'
        self.assertEqual(cf.model_attc_name, 'baz')

    def test_model_len(self):
        cf = config.Config(self.args)
        with self.assertRaises(RuntimeError) as ctx:
            with self.catch_log():
                cf.model_len
        self.assertEqual(str(ctx.exception), "'model_attc' is not define.")

        self.args.attc_model = 'foo'
        cf = config.Config(self.args)
        with self.assertRaises(IOError) as ctx:
            with self.catch_log():
                cf.model_len
        self.assertEqual(str(ctx.exception), "Path to model_attc '{}' does not exists".format(cf.model_attc_path))

        model_path = os.path.join(os.path.dirname(__file__), 'data', 'Replicons', 'acba.007.p01.13.fst')
        self.args.attc_model = model_path
        cf = config.Config(self.args)
        with self.assertRaises(RuntimeError) as ctx:
            with self.catch_log():
                cf.model_len
        self.assertEqual(str(ctx.exception), "CLEN not found in '{}', maybe it's not infernal model file".format(model_path))

        self.args.attc_model = 'attc_4.cm'
        cf = config.Config(self.args)
        cf._prefix_data = os.path.join(os.path.dirname(__file__), 'data')
        self.assertEqual(cf.model_len, 47)
        # test the model_len cache
        self.assertEqual(cf.model_len, 47)

    def test_func_annot_path(self):
        cf = config.Config(self.args)
        cf._prefix_data = 'foo'
        self.assertEqual(cf.func_annot_path, os.path.join('foo', 'Functional_annotation'))


    def test_log_level(self):
        for v, q, exp_level in [(0, 0, 20), (0, 2, 40), (0, 5, 50), (1, 0, 10), (3, 0, 10), (2, 2, 20)]:
            self.args.verbose = v
            self.args.quiet = q
            cf = config.Config(self.args)
            self.assertEqual(cf.log_level, exp_level)

    def test_prot_file(self):
        self.args.prot_file = 'nimportnaoik'
        cf = config.Config(self.args)
        self.assertEqual(cf.prot_file, self.args.prot_file)