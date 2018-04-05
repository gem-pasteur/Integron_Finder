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
import sys
import unittest
import platform
import colorlog
from io import StringIO
from contextlib import contextmanager

from integron_finder import IntegronError, logger_set_level

class IntegronTest(unittest.TestCase):

    _data_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "data"))

    @classmethod
    def find_data(cls, name):
        data_path = os.path.join(cls._data_dir, name)
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
        levels = {'NOTSET': colorlog.logging.logging.NOTSET,
                  'DEBUG': colorlog.logging.logging.DEBUG,
                  'INFO': colorlog.logging.logging.INFO,
                  'WARNING': colorlog.logging.logging.WARNING,
                  'ERROR': colorlog.logging.logging.ERROR,
                  'CRITICAL': colorlog.logging.logging.CRITICAL,
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
        self.maxDiff = None
        with open(f1) as fh1, open(f2) as fh2:
            self.assertMultiLineEqual(fh1.read(), fh2.read(), msg=msg)

    def assertSeqRecordEqual(self, s1, s2):
        for attr in ('id', 'name', 'seq'):
            s1_attr = getattr(s1, attr)
            s2_attr = getattr(s2, attr)
            self.assertEqual(s1_attr, s2_attr, msg="{} are different: {} != {}".format(attr, s1_attr, s2_attr))

        # there is a bug in some biopython version
        self.assertEqual(s1.description.rstrip('.'), s2.description.rstrip('.'))
        for s1_feat, s2_feat in zip(s1.features, s2.features):
            # location cannot be directly compared
            self.assertEqual(str(s1_feat.location), str(s2_feat.location))

            for attr in ('strand', 'type'):
                f1_attr = getattr(s1_feat, attr)
                f2_attr = getattr(s2_feat, attr)
                self.assertEqual(f1_attr, f2_attr, msg="{} are different: {} != {}".format(attr, f1_attr, f2_attr))

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


def which(name, flags=os.X_OK):
    """
    Search PATH for executable files with the given name.

    :param name: the name of the executable to search
    :type name: str
    :param flags: os mod the name must have, default is executable (os.X_OK).
    :type flags: os file mode R_OK|R_OK|W_OK|X_OK
    :return: the path of the executable
    :rtype: string or None
    """
    result = None
    path = os.environ.get('PATH', None)
    if path is None:
        return result
    for p in os.environ.get('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if platform.system() == 'Windows':
            p += '.exe'
        if os.access(p, flags):
            result = p
            break
    return result
