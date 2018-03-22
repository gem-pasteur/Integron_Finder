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

from . import utils
from . import __INTEGRON_DATA__

class Config(object):
    """
    Config object hold values issue from command lines
    """

    def __init__(self, args):
        self._model_len = None  # model_len cache, because it's computation is "heavy" (open file)
        self._args = args

        if __INTEGRON_DATA__ == '$' + 'INTEGRONDATA':
            self._prefix_data = os.path.join(os.path.join(os.path.dirname(__file__)), '..', 'data')
        else:
            self._prefix_data = os.path.join(__INTEGRON_DATA__, 'data')

    def __getattr__(self, item):
        try:
            attr = getattr(self._args, item)
            return attr
        except AttributeError:
            raise AttributeError("config object has no attribute '{}'".format(item))

    @property
    def replicon_path(self):
        """The absolute path to the replicon"""
        return os.path.abspath(self._args.replicon)

    @property
    def input_dir(self):
        """The absolute path to the directory where is located the replicon"""
        in_dir, sequence_file = os.path.split(self.replicon_path)
        return in_dir

    @property
    def outdir(self):
        """The absolute path where to write the results directory"""
        return os.path.abspath(self._args.outdir)

    @property
    def result_dir(self):
        """The absolute path to results directory"""
        file_name = utils.get_name_from_path(self.replicon_path)
        result_dir = os.path.join(self.outdir, "Results_Integron_Finder_" + file_name)
        return result_dir

    @property
    def default_topology(self):
        """The default topology
           available values are: 'circ' for circular or 'lin' for linear."""
        try:
            if self._args.circular:
                return 'circ'
            elif self._args.linear:
                return 'lin'
            else:
                return None
        except AttributeError:
            return None

    @property
    def model_dir(self):
        """The absolute path to the directory containing the models"""
        return os.path.join(self._prefix_data, "Models")

    @property
    def model_integrase(self):
        """The absolute path to the integrase model file"""
        return os.path.join(self.model_dir, "integron_integrase.hmm")

    @property
    def model_phage_int(self):
        return os.path.join(self.model_dir, "phage-int.hmm")

    @property
    def model_attc_path(self):
        """The absolute path to the attC model file"""
        try:
            self._args.attc_model
        except AttributeError:
            raise RuntimeError("'model_attc' is not define.")
        if len(self._args.attc_model.split(os.sep)) > 1:  # contain path
            model_attc = self._args.attc_model
        else:
            model_attc = os.path.join(self.model_dir, self._args.attc_model)
        return model_attc

    @property
    def model_attc_name(self):
        """The name of the attc model"""
        try:
            self._args.attc_model
        except AttributeError:
            raise RuntimeError("'model_attc' is not define.")
        return utils.get_name_from_path(self.model_attc_path)

    @property
    def model_len(self):
        """
        :return: The length of the attc model (corresponding to CLEN field).
        :raises: IOError if model_attc_path does match an existing file
                 RuntimeError if the file doe not content CLEN field.
        """
        try:
            self._args.attc_model
        except AttributeError:
            raise RuntimeError("'model_attc' is not define.")
        if self._model_len is None:
            model_len = utils.model_len(self.model_attc_path)
            self._model_len = model_len
            return model_len
        else:
            return self._model_len

    @property
    def func_annot_path(self):
        """The absolute path to the directory containing file needed for the functional annotation"""
        return os.path.join(self._prefix_data, "Functional_annotation")


    @property
    def log_level(self):
        """
        :return: the level to apply to loggers. 0 <= level <=50
        :rtype: int
        """
        default = 20  # info
        level = default - (10 * self._args.verbose) + (10 * self._args.quiet)
        level = max(0, level)
        level = min(50, level)
        return level
