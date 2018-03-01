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
from subprocess import call


def find_integrase(replicon_path, replicon, prot_file, out_dir, cfg):
    """
    Call Prodigal for Gene annotation and hmmer to find integrase, either with phage_int
    HMM profile or with intI profile.

    :param replicon_path: the path of the replicon to analyse
    :type replicon_path: string
    :param replicon_name: the name of the replicon to analyse
    :type replicon_name: string
    :param out_dir: the relative path to the directory where prodigal outputs will be stored
    :type out_dir: str
    :returns: None, the results are written on the disk
    """
    if not cfg.gembase:
        # Test whether the protein file exist to avoid new annotation for each run on the same replicon
        prot_tr_path = os.path.join(out_dir, replicon.name + ".prt")
        if not os.path.isfile(prot_tr_path):
            prodigal_cmd = "{prodigal} {meta} -i {replicon} -a {prot} -o {out}".format(
                prodigal=cfg.prodigal,
                meta='' if len(replicon) > 200000 else '-p meta',
                replicon=replicon_path,
                prot=prot_tr_path,
                out=os.devnull)
            try:
                returncode = call(prodigal_cmd.split())
            except Exception as err:
                raise RuntimeError("{0} failed : {1}".format(prodigal_cmd, err))
            if returncode != 0:
                raise RuntimeError("{0} failed returncode = {1}".format(prodigal_cmd, returncode))

    intI_hmm_out = os.path.join(out_dir, replicon.name + "_intI.res")
    hmm_cmd = []
    if not os.path.isfile(intI_hmm_out):
        hmm_cmd.append([cfg.hmmsearch,
                        "--cpu", str(cfg.cpu),
                        "--tblout", os.path.join(out_dir, replicon.name + "_intI_table.res"),
                        "-o", intI_hmm_out,
                        cfg.model_integrase,
                        prot_file])

    phage_hmm_out = os.path.join(out_dir, replicon.name + "_phage_int.res")
    if not os.path.isfile(phage_hmm_out):
        hmm_cmd.append([cfg.hmmsearch,
                        "--cpu", str(cfg.cpu),
                        "--tblout", os.path.join(out_dir, replicon.name + "_phage_int_table.res"),
                        "-o", phage_hmm_out,
                        cfg.model_phage_int,
                        prot_file])

    for cmd in hmm_cmd:
        try:
            returncode = call(cmd)
        except Exception as err:
            raise RuntimeError("{0} failed : {1}".format(' '.join(cmd), err))
        if returncode != 0:
            raise RuntimeError("{0} failed return code = {1}".format(' '.join(cmd), returncode))

