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
import glob
import colorlog

import numpy as np
import pandas as pd

# display warning only for non installed integron_finder
from Bio import BiopythonExperimentalWarning
import warnings
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

from Bio import SearchIO


_log = colorlog.getLogger(__name__)


def scan_hmm_bank(path):
    """
    :param str path: - if the path is a dir:
                       return all files ending with .hmm in the dir
                     - if the path is a file:
                       parse the file, each line must be an expression (glob)
                       pointing to hmm files
    :return: lists of hmm files to consider for annotation
    :rtype: list of str
    :raise IOError: if the path does not exists
    """
    real_path = os.path.realpath(path)
    files = []
    if os.path.exists(real_path):
        if os.path.isdir(real_path):
            files = glob.glob(os.path.join(real_path, '*.hmm'))
        elif os.path.isfile(real_path):
            with open(real_path) as hmm_bank:
                for bank_path in hmm_bank:
                    if bank_path.startswith('#'):
                        continue
                    elif not os.path.isabs(bank_path):
                        if "_prefix_share" in globals():
                            prefix = _prefix_share
                        else:
                            prefix = os.environ['INTEGRON_HOME']
                        bank_path = os.path.normpath(os.path.join(prefix, bank_path))
                    bank_files = glob.glob(os.path.expanduser(bank_path.strip("\n").strip()))
                    if not bank_files:
                        _log.warning("func_annot '{}' does not match any files.".format(bank_path))
                    else:
                        for path in bank_files:
                            _log.warning("the hmm {} will be used for functional annotation".format(path))
                        files.extend(bank_files)

        return files
    else:
        raise IOError("{} no such file or directory".format(path))


def read_hmm(replicon_name, infile, cfg, evalue=1., coverage=0.5):
    """
    Function that parse hmmer --out output and returns a pandas DataFrame
    filter output by evalue and coverage. (Being % of the profile aligned)

    :param str replicon_name: the name of the replicon
    :param str infile: the hmm output (in tabulated format) to parse
    :param cfg: the config
    :type cfg: :class:`integron_finder.config.Config` object.
    :param float evalue: filter out hitst with evalue greater tha evalue.
    :param float coverage: filter out hits with coverage under coverage (% of the profile aligned)
    :returns: data Frame with columns:
              | "Accession_number", "query_name", "ID_query", "ID_prot", "strand", "pos_beg", "pos_end", "evalue"
              | each row correspond to a hit.
    :rtype: a :class:`pandas.DataFrame`
    """

    df = pd.DataFrame(columns=["Accession_number", "query_name", "ID_query",
                               "ID_prot", "strand", "pos_beg", "pos_end",
                               "evalue", "hmmfrom", "hmmto", "alifrom",
                               "alito", "len_profile"])
    gen = SearchIO.parse(infile, 'hmmer3-text')
    for idx, query_result in enumerate(gen):
        len_profile = query_result.seq_len
        query = query_result.id

        try:
            id_query = query_result.accession
        except AttributeError:
            id_query = "-"
        for idx2, hit in enumerate(query_result.hits):
            id_prot = hit.id
            if not cfg.gembase:
                pos_beg, pos_end, strand = [x.strip() for x in hit.description_all[0].split('#') if x][:-1]
            else:
                desc = [j for j in hit.description_all[0].split(" ")]
                strand = 1 if desc[0] == "D" else -1
                pos_beg = int(desc[3])
                pos_end = int(desc[4])

            evalue_tmp = []
            hmmfrom = []
            hmmto = []
            alifrom = []
            alito = []

            for hsp in hit.hsps:
                # strand = hsp.query_strand
                evalue_tmp.append(hsp.evalue)
                hmmfrom.append(hsp.query_start + 1)
                hmmto.append(hsp.query_end)
                alifrom.append(hsp.hit_start + 1)
                alito.append(hsp.hit_end)

            best_evalue = np.argmin(evalue_tmp)

            df.loc[idx+idx2, "ID_prot"] = id_prot
            df.loc[idx+idx2, "ID_query"] = id_query # "-"  # remnant of ancient parsing function to keep data structure
            df.loc[idx+idx2, "pos_beg"] = int(pos_beg)
            df.loc[idx+idx2, "pos_end"] = int(pos_end)
            df.loc[idx+idx2, "strand"] = int(strand)
            df.loc[idx+idx2, "evalue"] = evalue_tmp[best_evalue]   # i-evalue
            df.loc[idx+idx2, "hmmfrom"] = hmmfrom[best_evalue]   # hmmfrom
            df.loc[idx+idx2, "hmmto"] = hmmto[best_evalue]     # hmm to
            df.loc[idx+idx2, "alifrom"] = alifrom[best_evalue]   # alifrom
            df.loc[idx+idx2, "alito"] = alito [best_evalue]  # ali to
            df.loc[idx+idx2, "len_profile"] = float(len_profile)
            df.loc[idx+idx2, "Accession_number"] = replicon_name
            df.loc[idx+idx2, "query_name"] = query

    intcols = ["pos_beg", "pos_end", "strand"]
    floatcol = ["evalue", "len_profile"]
    df[intcols] = df[intcols].astype(int)
    df[floatcol] = df[floatcol].astype(float)

    df = df[(((df.hmmto - df.hmmfrom) / df.len_profile) > coverage) & (df.evalue < evalue)]
    df.index = range(len(df))

    df_out = df[["Accession_number", "query_name", "ID_query", "ID_prot",
                 "strand", "pos_beg", "pos_end", "evalue"]]

    return df_out
