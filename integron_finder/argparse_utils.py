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

import os.path
import argparse
import integron_finder


class MyVersionAction(argparse._VersionAction):

    def __call__(self, parser, namespace, values, option_string=None):
        """
        Override the :meth:`argparse._VersionAction.__call__`
        to allow to call get_version_message with the

        * --hmmsearch
        * --cmsearch
        * --prodigal

        value provided on the command line

        I ensure that the --version option is the last argument
        so all other options are parsed and are in namespace
        then --version is parsed then code below is executed.
        """
        version = integron_finder.get_version_message(hmmsearch=namespace.hmmsearch,
                                                      cmsearch=namespace.cmsearch,
                                                      prodigal=namespace.prodigal)
        formatter = parser._get_formatter()
        formatter.add_text(version)
        parser._print_message(formatter.format_help(), argparse._sys.stdout)
        parser.exit()


def path(value):
    if os.path.exists(value):
        return value
    else:
        raise argparse.ArgumentTypeError(f"{value} no such file or directory")