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
import sys
import unittest
import functools
import colorlog
from io import StringIO
from contextlib import contextmanager

import pandas as pd

import pandas.testing as pdt

from integron_finder import IntegronError, logger_set_level, get_logging_module


class IntegronTest(unittest.TestCase):
    # the data for tests
    _data_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "data"))
    maxDiff = None
    
    @classmethod
    def find_data(cls, *names):
        data_path = os.path.join(cls._data_dir, *names)
        if os.path.exists(data_path):
            return data_path
        else:
            raise IOError("data '{}' does not exists".format(data_path))


    @contextmanager
    def catch_io(self, out=False, err=False):
        """
        Catch stderr and stdout of the code running within this block.
        """
        old_out = sys.stdout
        new_out = old_out
        old_err = sys.stderr
        new_err = old_err
        if out:
            new_out = StringIO()
        if err:
            new_err = StringIO()
        try:
            sys.stdout, sys.stderr = new_out, new_err
            yield sys.stdout, sys.stderr
        finally:
            sys.stdout, sys.stderr = old_out, old_err

    @contextmanager
    def catch_log(self):
        logger = colorlog.getLogger('integron_finder')
        handlers_ori = logger.handlers
        fake_handler = colorlog.StreamHandler(StringIO())
        try:
            logger.handlers = [fake_handler]
            yield LoggerWrapper(logger)
        finally:
            logger.handlers = handlers_ori

    @classmethod
    def set_log_level(cls, level):
        logging = get_logging_module()
        levels = {'NOTSET': logging.NOTSET,
                  'DEBUG': logging.DEBUG,
                  'INFO': logging.INFO,
                  'WARNING': logging.WARNING,
                  'ERROR': logging.ERROR,
                  'CRITICAL': logging.CRITICAL,
        }
        if level in levels:
            level = levels[level]
        elif not isinstance(level, int):
            raise IntegronError("Level must be {} or a positive integer")
        elif level < 0:
            raise IntegronError("Level must be {} or a positive integer")

        logger_set_level(level)


    @staticmethod
    def fake_exit(*args, **kwargs):
        returncode = args[0]
        raise TypeError(returncode)


    @staticmethod
    def mute_call(call_ori):
        """
        hmmsearch or prodigal write lot of things on stderr or stdout
        which noise the unit test output
        So I replace the `call` function in module integron_finder
        by a wrapper which call the original function but add redirect stderr and stdout
        in dev_null
        :return: wrapper around call function
        :rtype: function
        """
        def wrapper(*args, **kwargs):
            with open(os.devnull, 'w') as f:
                kwargs['stderr'] = f
                kwargs['stdout'] = f
                res = call_ori(*args, **kwargs)
            return res
        return wrapper


    def assertFileEqual(self, f1, f2, msg=None):
        #self.maxDiff = None
        with open(f1) as fh1, open(f2) as fh2:
            self.assertMultiLineEqual(fh1.read(), fh2.read(), msg=msg)

    def assertIntegronResultEqual(self, f1, f2, debug=False):
        try:
            df1 = pd.read_csv(f1, sep="\t", comment='#')
        except pd.errors.EmptyDataError:
            df1 = pd.DataFrame()
        try:
            df2 = pd.read_csv(f2, sep="\t", comment='#')
        except pd.errors.EmptyDataError:
            df2 = pd.DataFrame()
        if debug:
            print("\n############################")
            print(df1)
            print("==============================")
            print(df2)
            print("##############################")
        pdt.assert_frame_equal(df1, df2)


    def assertSeqRecordEqual(self, s1, s2):
        for attr in ('id', 'name', 'seq'):
            s1_attr = getattr(s1, attr)
            s2_attr = getattr(s2, attr)
            self.assertEqual(s1_attr, s2_attr, msg="{} are different: {} != {}".format(attr, s1_attr, s2_attr))

        # there is a bug in some biopython version
        self.assertEqual(s1.description.rstrip('.'), s2.description.rstrip('.'))
        for s1_feat, s2_feat in zip(s1.features, s2.features):
            # location cannot be directly compared
            # strand is now in location not anylonger in feature
            self.assertEqual(str(s1_feat.location), str(s2_feat.location))

            # for attr in ('strand', 'type'):
            #     f1_attr = getattr(s1_feat, attr)
            #     f2_attr = getattr(s2_feat, attr)
            self.assertEqual(s1_feat.type, s2_feat.type, msg="{} are different: {} != {}".format(attr, s1_feat, s2_feat))

            # The order of qualifers does not matter
            # ('integron_type', ['complete']), ('integron_id', ['integron_01'])
            # ('integron_id', ['integron_01']), ('integron_type', ['complete'])
            self.assertDictEqual(dict(s1_feat.qualifiers), dict(s2_feat.qualifiers))

    def assertHmmEqual(self, hmm1, hmm2):
        with open(hmm1) as hmm1_file, open(hmm2) as hmm2_file:
            for hmm1_line, hmm2_line in zip(hmm1_file, hmm2_file):
                if hmm1_line.startswith('#') and hmm2_line.startswith('#'):
                    continue
                hmm1_fields = hmm1_line.split('#')[:-1]
                hmm2_fields = hmm2_line.split('#')[:-1]
                self.assertListEqual(hmm1_fields, hmm2_fields)


class LoggerWrapper(object):

    def __init__(self, logger):
        self.logger = logger

    def __getattr__(self, item):
        return getattr(self.logger, item)

    def get_value(self):
        return self.logger.handlers[0].stream.getvalue()


def hide_executable(bin_2_hide):
    """
    This a decorator maker, it return a decorator which can be used to decorate a "which" like function
    the decorator call the "which" like function except for value of bin_2_hide in this case it return None
    to simulate that the "which" like function does not find any corresponding executable.

    :param bin_2_hide: the name of the binary to hide
    :return: a decorator
    """
    def which(func):
        @functools.wraps(func)
        def wrapper(exe):
            if exe == bin_2_hide:
                return None
            else:
                return func(exe)
        return wrapper
    return which

