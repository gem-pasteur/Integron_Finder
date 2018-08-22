# -*- coding: utf-8 -*-

####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright (c) 2015 - 2018  Institut Pasteur, Paris and CNRS.                     #
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
import colorlog
import pandas as pd
from Bio import SeqIO

from .utils import model_len

_log = colorlog.getLogger(__name__)


def read_infernal(infile, replicon_id, len_model_attc,
                  evalue=1, size_max_attc=200, size_min_attc=40):
    """
    Function that parse cmsearch --tblout output and returns a pandas DataFrame

    :param str infile: the path to the output of cmsearch in tabulated format (--tblout)
    :param str replicon_id: the id of the replicon are the integrons were found.
    :param int len_model_attc: the length of the attc model
    :param float evalue: evalue threshold to filter out hits above it
    :param int size_max_attc: The maximum value fot the attC size
    :param int size_min_attc: The minimum value fot the attC size
    :return: table with columns:

            | "Accession_number", "cm_attC", "cm_debut", "cm_fin", "pos_beg", "pos_end", "sens", "evalue"
            | and each row is a hit that match the attc covariance model.
    :rtype: :class:`pandas.DataFrame` object
    """
    _log.debug("read_infernal {}, {}, {}, evalue={}, size_max_attc={}, size_min_attc={}".format(
        infile, replicon_id, len_model_attc, evalue, size_max_attc, size_min_attc
    ))
    try:
        _ = pd.read_table(infile, comment="#")
    except:
        return pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                     "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])

    df = pd.read_table(infile, sep="\s+", engine="python",  header=None,
                       skipfooter=10, skiprows=2, usecols=[2, 5, 6, 7, 8, 9, 15])
    # some line can have different number of columns due to difference in description
    # we do not use this columns so we must parse only cols we need
    # Keep only columns: query_name(2), mdl from(5), mdl to(6), seq from(7),
    # seq to(8), strand(9), E-value(15)
    df.columns = ["cm_attC", "cm_debut", "cm_fin", "pos_beg_tmp", "pos_end_tmp", "sens", "evalue"]
    df["Accession_number"] = replicon_id
    df = df[df.evalue < evalue]  # filter on evalue
    df = df[(abs(df.pos_end_tmp - df.pos_beg_tmp) < size_max_attc) & (size_min_attc < abs(df.pos_end_tmp - df.pos_beg_tmp))]
    if not df.empty:
        df.sort_values(['pos_end_tmp', 'evalue'], inplace=True)
        df.index = list(range(0, len(df)))
        idx = (df.pos_beg_tmp > df.pos_end_tmp)
        df.loc[idx, "pos_beg"] = df.loc[idx].apply(lambda x: x["pos_end_tmp"] - (len_model_attc - x["cm_fin"]), axis=1)
        df.loc[idx, "pos_end"] = df.loc[idx].apply(lambda x: x["pos_beg_tmp"] + (x["cm_debut"] - 1), axis=1)

        df.loc[~idx, "pos_end"] = df.loc[~idx].apply(lambda x: x["pos_end_tmp"] + (len_model_attc - x["cm_fin"]), axis=1)
        df.loc[~idx, "pos_beg"] = df.loc[~idx].apply(lambda x: x["pos_beg_tmp"] - (x["cm_debut"] - 1), axis=1)

        return df[["Accession_number", "cm_attC", "cm_debut", "cm_fin", "pos_beg", "pos_end", "sens", "evalue"]]
    else:
        return pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                     "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])


def local_max(replicon,
              window_beg, window_end,
              model_attc_path,
              strand_search="both",
              evalue_attc=1., max_attc_size=200, min_attc_size=40,
              cmsearch_bin='cmsearch', out_dir='.', cpu_nb=1):
    """
    :param replicon: The name of replicon (without suffix)
    :type replicon: :class:`Bio.Seq.SeqRecord` object.
    :param int window_beg: start of window to search for attc (position of protein)
    :param int window_end: end of window to search for attc (position of protein)
    :param str strand_search: The strand on which to looking for attc.
                              Available values:

                                * 'top': Only search the top (Watson) strand of target sequences.
                                * 'bottom': Only search the bottom (Crick) strand of target sequences
                                * 'both': search on both strands

    :param float evalue_attc: evalue threshold to filter out hits above it
    :param int max_attc_size: The maximum value fot the attC size
    :param int min_attc_size: The minimum value fot the attC size
    :param str cmsearch_bin: the path to cmsearch
    :param str out_dir: the path to directory where to write results
    :param int cpu_nb: The number of cpu used by cmsearch
    :return: DataFrame with same structure as the DataFrame returns by :func:`read_infernal`
             where position are converted on position on replicon and attc are filtered
             by evalue, min_attc_size, max_attc_size
             also write a file with intermediate results <replicon_id>_subseq_attc_table_end.res
             this file store the local_max results before filtering by max_attc_size and min_attc_size
    :rtype: :class:`pandas.DataFrame` object
    """
    replicon_size = len(replicon)
    if window_beg < window_end:
        subseq = replicon[window_beg:window_end]
    else:
        # the window overlap the replicon origin
        subseq1 = replicon[window_beg:]
        subseq2 = replicon[:window_end]
        subseq = subseq1 + subseq2

    infile_path = os.path.join(out_dir, replicon.id + "_subseq.fst")
    with open(infile_path, "w") as f:
        SeqIO.write(subseq, f, "fasta")

    output_path = os.path.join(out_dir, "{name}_{win_beg}_{win_end}_subseq_attc.res".format(name=replicon.id,
                                                                                            win_beg=window_beg,
                                                                                            win_end=window_end))
    tblout_path = os.path.join(out_dir, "{name}_{win_beg}_{win_end}_subseq_attc_table.res".format(name=replicon.id,
                                                                                                  win_beg=window_beg,
                                                                                                  win_end=window_end))

    cmsearch_cmd = \
        "{bin} -Z {size} {strand} --max --cpu {cpu} -o {out} --tblout {tblout} -E 10 {mod_attc_path} {infile}".format(
            bin=cmsearch_bin,
            size=replicon_size / 1000000.,
            strand={"both": "", "top": "--toponly", "bottom": "--bottomonly"}[strand_search],
            cpu=cpu_nb,
            out=output_path,
            tblout=tblout_path,
            mod_attc_path=model_attc_path,
            infile=infile_path)
    try:
        _log.debug("run cmsearch: {}".format(cmsearch_cmd))
        returncode = call(cmsearch_cmd.split())
    except Exception as err:
        raise RuntimeError("{0} failed : {1}".format(cmsearch_cmd, err))
    if returncode != 0:
        raise RuntimeError("{0} failed returncode = {1}".format(cmsearch_cmd, returncode))

    df_max = read_infernal(tblout_path,
                           replicon.id, model_len(model_attc_path),
                           evalue=evalue_attc,
                           size_max_attc=max_attc_size,
                           size_min_attc=min_attc_size)

    # if replicon is linear
    # df_max.pos_beg + window_beg is always < replicon_size
    # (df_max.pos_beg + window_beg) % replicon_size = (df_max.pos_beg + window_beg)
    # if replicon is circular and attc site overlap origin
    # df_max.pos_beg + window_beg  > replicon_size
    # (df_max.pos_beg + window_beg) % replicon_size is position on replicon
    # for instance with pos = 100 and replicon size = 90
    # 100 % 90 = 10
    df_max.pos_beg = (df_max.pos_beg + window_beg) % replicon_size
    df_max.pos_end = (df_max.pos_end + window_beg) % replicon_size
    df_max.to_csv(os.path.join(out_dir, replicon.id + "_subseq_attc_table_end.res"),
                  sep="\t", index=0, mode="a", header=0)
    # filter on size
    df_max = df_max[(abs(df_max.pos_end - df_max.pos_beg) > min_attc_size) &
                    (abs(df_max.pos_end - df_max.pos_beg) < max_attc_size)]
    return df_max


def expand(replicon,
           window_beg, window_end, max_elt, df_max,
           circular, dist_threshold, max_attc_size,
           model_attc_path,
           search_left=False, search_right=False,
           out_dir='.', cpu=1):
    """
    for a given element, we can search on the left hand side (if integrase is on the right for instance)
    or right hand side (opposite situation) or both side (only integrase or only attC sites)

    :param replicon: The Replicon to annotate
    :type replicon: a :class:`Bio.Seq.SeqRecord` object.
    :param int window_beg: start of window to search for attc (position of protein)
    :param int window_end: end of window to search for attc (position of protein)
    :param max_elt: DataFrame with columns:
        ::

            Accession_number cm_attC  cm_debut  cm_fin   pos_beg   pos_end sens   evalue

        and each row is an occurrence of attc site

    :type max_elt: :class:`pandas.DataFrame` object
    :param df_max: DataFrame with columns
        ::

            Accession_number cm_attC  cm_debut  cm_fin   pos_beg   pos_end sens   evalue

        and each row is an occurrence of attc site

    :type df_max: :class:`pandas.DataFrame` object
    :param bool circular: True if replicon topology is circular otherwise False.
    :param int dist_threshold: Two elements are aggregated if they are distant of dist_threshold [4kb] or less
    :param int max_attc_size: The maximum value fot the attC size
    :param str model_attc_path: the path to the attc model file
    :param bool search_left: need to search on right of ???
    :param bool search_right: need to search on right of ???
    :param str out_dir: The path to directory where to write results
    :param int cpu: the number of cpu use by expand
    :return: a copy of max_elt with attC hits
    :rtype: :class:`pandas.DataFrame` object

    """
    replicon_size = len(replicon)
    # for a given element, we can search on the left hand side (if integrase is on the right for instance)
    # or right hand side (opposite situation) or both side (only integrase or only attC sites)
    wb = window_beg
    we = window_end

    if search_right:

        if circular:
            window_beg = (window_end - max_attc_size) % replicon_size
            # max_attc_size (200bo by default) to allow the detection of sites that would overlap 2 consecutive windows
            window_end = (window_end + dist_threshold) % replicon_size
        else:
            window_beg = max(0, window_end - max_attc_size)
            # max_attc_size (200bo by default) to allow the detection of sites that would overlap 2 consecutive windows
            window_end = min(replicon_size, window_end + dist_threshold)

        searched_strand = "both" if search_left else "top"  # search on both strands if search in both directions

        while not df_max.empty and 0 < (window_beg and window_end) < replicon_size:
            df_max = local_max(replicon,
                               window_beg, window_end,
                               model_attc_path,
                               strand_search=searched_strand,
                               out_dir=out_dir, cpu_nb=cpu
                               )
            max_elt = pd.concat([max_elt, df_max])

            if circular:
                window_beg = (window_end - max_attc_size) % replicon_size
                window_end = (window_end + dist_threshold) % replicon_size
            else:
                window_beg = max(0, window_end - max_attc_size)
                window_end = min(replicon_size, window_end + dist_threshold)

        # re-initialize in case we enter search left too.
        df_max = max_elt.copy()
        window_beg = wb
        window_end = we

    if search_left:
        if circular:
            window_end = (window_beg + 200) % replicon_size
            window_beg = (window_beg - dist_threshold) % replicon_size
        else:
            window_beg = max(0, window_beg - dist_threshold)
            window_end = min(replicon_size, window_beg + 200)

        searched_strand = "both" if search_right else "bottom"

        while not df_max.empty and 0 < (window_beg and window_end) < replicon_size:

            df_max = local_max(replicon,
                               window_beg, window_end,
                               model_attc_path,
                               strand_search=searched_strand,
                               out_dir=out_dir,
                               cpu_nb=cpu)
            max_elt = pd.concat([max_elt, df_max])  # update of attC list of hits.

            if circular:
                window_end = (window_beg + 200) % replicon_size
                window_beg = (window_beg - dist_threshold) % replicon_size
            else:
                window_end = min(replicon_size, window_beg + 200)
                window_beg = max(0, window_beg - dist_threshold)

    max_elt.drop_duplicates(inplace=True)
    max_elt.index = list(range(len(max_elt)))
    return max_elt


