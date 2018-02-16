# -*- coding: utf-8 -*-

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

    def test_getattr(self):
        self.args.replicon = 'foo'
        self.args.bar = 'baz'
        cf = config.Config(self.args)
        self.assertEqual(cf.replicon_path, os.path.abspath('foo'))
        self.assertEqual(cf.bar, 'baz')
        with self.assertRaises(AttributeError) as ctx:
            self.assertEqual(cf.foobar, 'foobar')
        self.assertEqual(str(ctx.exception), "config object has no attribute 'foobar'")

    def test_replicon_path(self):
        self.args.replicon = '../foo'
        cf = config.Config(self.args)
        self.assertEqual(cf.replicon_path, os.path.abspath('../foo'))

    def test_replicon_name(self):
        self.args.replicon = '../foo.fasta'
        cf = config.Config(self.args)
        self.assertEqual(cf.replicon_name, 'foo')

    def test_input_dir(self):
        self.args.replicon = '../foo.fasta'
        cf = config.Config(self.args)
        self.assertEqual(cf.input_dir, os.path.split(os.path.abspath('../foo'))[0])

    def test_mode_name(self):
        self.args.eagle_eyes = True
        self.args.local_max = False
        cf = config.Config(self.args)
        self.assertEqual(cf.mode_name, 'local_max')
        self.args.eagle_eyes = False
        self.args.local_max = True
        cf = config.Config(self.args)
        self.assertEqual(cf.mode_name, 'local_max')
        self.args.eagle_eyes = False
        self.args.local_max = False
        cf = config.Config(self.args)
        self.assertEqual(cf.mode_name, 'default')

    def test_default_topology(self):
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

    def test_model_attc(self):
        self.args.attc_model = 'bar'
        cf = config.Config(self.args)
        cf._prefix_data = 'foo'
        self.assertEqual(cf.model_attc, os.path.join('foo', 'Models', 'bar'))
        self.args.attc_model = 'bar/baz'
        self.assertEqual(cf.model_attc, 'bar/baz')

    def test_func_annot_path(self):
        cf = config.Config(self.args)
        cf._prefix_data = 'foo'
        self.assertEqual(cf.func_annot_path, os.path.join('foo', 'Functional_annotation'))
