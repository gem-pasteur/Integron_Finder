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
import glob
import warnings

import colorlog
import numpy as np
import pandas as pd

# display warning only for non installed integron_finder
from Bio import BiopythonExperimentalWarning
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
            files.extend(glob.glob(os.path.join(real_path, '*.HMM')))
        elif os.path.isfile(real_path):
            with open(real_path) as hmm_bank:
                wrong_lines = 0
                for bank_path in hmm_bank:
                    bank_path = bank_path.strip()
                    if bank_path.startswith('#'):
                        continue
                    elif not os.path.isabs(bank_path):
                        bank_path = os.path.normpath(bank_path)
                    bank_files = glob.glob(os.path.expanduser(bank_path))
                    if not bank_files:
                        wrong_lines += 1
                        _log.warning("func_annot '{}' does not match any files.".format(bank_path))
                        if wrong_lines > 10:
                            msg = "Too many lines with no hmm file in {}.\nIs there right file?" \
                                  "\nsee documentation tutorial functional annotation".format(real_path)
                            _log.error(msg)
                            raise ValueError(msg)
                    else:
                        for path in bank_files:
                            _log.warning("the hmm {} will be used for functional annotation".format(path))
                        files.extend(bank_files)

        return files
    else:
        raise IOError("{} no such file or directory".format(path))


def read_hmm(replicon_id, prot_db, infile, cfg, evalue=1., coverage=0.5):
    """
    Function that parse hmmer --out output and returns a pandas DataFrame
    filter output by evalue and coverage. (Being % of the profile aligned)

    :param str replicon_id: the id of the replicon
    :param prot_db: The protein database corresponding to the replicon translation
    :type prot_db: :class:`integron_finder.prot_db.ProteinDB` object.
    :param str infile: the hmm output (in tabulated format) to parse
    :param cfg: the config
    :type cfg: :class:`integron_finder.config.Config` object.
    :param float evalue: filter out hits with evalue greater tha evalue.
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
    _log.debug("Parse {}".format(infile))
    gen = SearchIO.parse(infile, 'hmmer3-text')
    for idx, query_result in enumerate(gen):
        len_profile = query_result.seq_len
        query = query_result.id

        try:
            id_query = query_result.accession
        except AttributeError:
            id_query = "-"
        for idx2, hit in enumerate(query_result.hits):
            if not hit.hsps:
                _log.debug(f"{hit.id} does not contains any hsps. skip it.")
                continue
            id_prot = hit.id

            _, strand, pos_beg, pos_end = prot_db.get_description(hit.id)

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
            df.loc[idx+idx2, "ID_query"] = id_query  # "-"  # remnant of ancient parsing function to keep data structure
            df.loc[idx+idx2, "pos_beg"] = pos_beg
            df.loc[idx+idx2, "pos_end"] = pos_end
            df.loc[idx+idx2, "strand"] = strand
            df.loc[idx+idx2, "evalue"] = evalue_tmp[best_evalue]   # i-evalue
            df.loc[idx+idx2, "hmmfrom"] = hmmfrom[best_evalue]   # hmmfrom
            df.loc[idx+idx2, "hmmto"] = hmmto[best_evalue]     # hmm to
            df.loc[idx+idx2, "alifrom"] = alifrom[best_evalue]   # alifrom
            df.loc[idx+idx2, "alito"] = alito[best_evalue]  # ali to
            df.loc[idx+idx2, "len_profile"] = float(len_profile)
            df.loc[idx+idx2, "Accession_number"] = replicon_id
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
