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
import subprocess
import shlex
import colorlog
import pandas as pd
from Bio import SeqIO

from .utils import model_len

_log = colorlog.getLogger(__name__)


def read_infernal(infile, replicon_id, replicon_size,
                  len_model_attc,
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
    _log.debug(f"read_infernal {infile}, {replicon_id}, {replicon_size}, {len_model_attc}, "
               f"evalue={evalue}, size_max_attc={size_max_attc}, size_min_attc={size_min_attc}")
    dtype = {"Accession_number": "str",
             "cm_attC": "str",
             "cm_debut": "int",
             "cm_fin": "int",
             "pos_beg": "int",
             "pos_end": "int",
             "evalue": "float",
             }
    try:
        _ = pd.read_csv(infile, comment="#", sep="\t")
    except Exception:
        df = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                   "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
        return df.astype(dtype)

    df = pd.read_csv(infile, sep="\\s+", engine="python",  header=None,
                     comment='#',
                     usecols=[2, 5, 6, 7, 8, 9, 15])
    # some line can have different number of columns due to difference in description
    # we do not use these columns, so we must parse only cols we need
    # Keep only columns: query_name(2), mdl from(5), mdl to(6), seq from(7),
    # seq to(8), strand(9), E-value(15)
    df.columns = ["cm_attC", "cm_debut", "cm_fin", "pos_beg_tmp", "pos_end_tmp", "sens", "evalue"]
    df["Accession_number"] = replicon_id
    _log.debug(f"Before filtering on evalue {evalue}, there were {len(df)} attC sites")
    df = df[df.evalue < evalue]  # filter on evalue
    _log.debug(f"After filtering on evalue {evalue}, there are now {len(df)} attC sites")
    df = df[(abs(df.pos_end_tmp - df.pos_beg_tmp) < size_max_attc) &
            (size_min_attc < abs(df.pos_end_tmp - df.pos_beg_tmp))]
    _log.debug(f"After filtering on size max: {size_max_attc} and size min: {size_min_attc}, "
               f"there are now {len(df)} attC sites")
    if not df.empty:
        df.sort_values(['pos_end_tmp', 'evalue'], inplace=True)
        df.index = list(range(0, len(df)))
        idx = (df.pos_beg_tmp > df.pos_end_tmp)
        df.loc[idx, "pos_beg"] = df.loc[idx].apply(lambda x: max(x["pos_end_tmp"] - (len_model_attc - x["cm_fin"]),
                                                                 1),
                                                   axis=1)
        df.loc[idx, "pos_end"] = df.loc[idx].apply(lambda x: min(x["pos_beg_tmp"] + (x["cm_debut"] - 1),
                                                                 replicon_size),
                                                   axis=1)

        df.loc[~idx, "pos_beg"] = df.loc[~idx].apply(lambda x: max(x["pos_beg_tmp"] - (x["cm_debut"] - 1),
                                                                   1),
                                                     axis=1)
        df.loc[~idx, "pos_end"] = df.loc[~idx].apply(lambda x: min(x["pos_end_tmp"] + (len_model_attc - x["cm_fin"]),
                                                                   replicon_size),
                                                     axis=1)

        df = df[["Accession_number", "cm_attC", "cm_debut", "cm_fin", "pos_beg", "pos_end", "sens", "evalue"]]
        df["cm_attC"] = df["cm_attC"].str.lower()
    else:
        df = pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                     "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])
    return df.astype(dtype)


def find_attc(replicon_path, replicon_id, cmsearch_path, out_dir, model_attc, incE=1., cpu=1):
    """
    Call cmsearch to find attC sites in a single replicon.

    :param str replicon_path: the path of the fasta file representing the replicon to analyse.
    :param str replicon_id: the id of the replicon to analyse.
    :param str cmsearch_path: the path to the cmsearch executable.
    :param str out_dir: the path to the directory where cmsearch outputs will be stored.
    :param str model_attc: path to the attc model (Covariance Matrix).
    :param float incE: consider sequences <= this E-value threshold as significant (to get the alignment with -A)
    :param int cpu: the number of cpu used by cmsearch.
    :returns: None, the results are written on the disk.
    :raises RuntimeError: when cmsearch run failed.
    """
    cmsearch_cmd = '{cmsearch} --cpu {cpu} -A {out} --tblout {tblout_path} ' \
                   '-E 10 --incE {incE} {mod_attc} {infile}'.format(cmsearch=cmsearch_path.replace(' ', '\\ '),
                                                                    cpu=cpu,
                                                                    out=os.path.join(out_dir,
                                                                                     replicon_id + "_attc.res").replace(' ', '\\ '),
                                                                    tblout_path=os.path.join(out_dir,
                                                                                             replicon_id +
                                                                                             "_attc_table.res").replace(' ', '\\ '),
                                                                    incE=incE,
                                                                    mod_attc=model_attc.replace(' ', '\\ '),
                                                                    infile=replicon_path.replace(' ', '\\ '))
    try:
        _log.debug(f"run cmsearch: {cmsearch_cmd}")
        with open(os.devnull, 'w') as dev_null:
            cmd = shlex.split(cmsearch_cmd)
            completed_process = subprocess.run(cmd, stdout=dev_null)
    except Exception as err:
        raise RuntimeError(f"{' '.join(cmd)} failed : {err}")
    if completed_process.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed returncode = {completed_process.returncode}")


def local_max(replicon,
              window_beg, window_end,
              model_attc_path,
              strand_search="both",
              evalue_attc=1., max_attc_size=200, min_attc_size=40,
              cmsearch_bin='cmsearch', out_dir='.', cpu=1):
    """
    :param replicon: The replicon to analyse
    :type replicon: :class:`Bio.Seq.SeqRecord` object.
    :param int window_beg: Start of window to search for attc (position of protein).
    :param int window_end: End of window to search for attc (position of protein).
    :param str model_attc_path: The path to the covariance model for attc (eg: attc_4.cm)
                                used by cmsearch to find attC sites
    :param str strand_search: The strand on which to looking for attc.
                              Available values:

                                * 'top': Only search the top (Watson) strand of target sequences.
                                * 'bottom': Only search the bottom (Crick) strand of target sequences
                                * 'both': search on both strands

    :param float evalue_attc: evalue threshold to filter out hits above it
    :param int max_attc_size: The maximum value fot the attC size
    :param int min_attc_size: The minimum value fot the attC size
    :param str cmsearch_bin: The path to cmsearch
    :param str out_dir: The path to directory where to write results
    :param int cpu: The number of cpu used by cmsearch
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
        '{bin} -Z {size} {strand} --max --cpu {cpu} -A {out} --tblout {tblout} -E 10 ' \
        '--incE {incE} {mod_attc_path} {infile}'.format(bin=cmsearch_bin.replace(' ', '\\ '),
                                                        size=replicon_size / 1000000.,  # search space size in *Mb*
                                                        strand={"both": "",
                                                                "top": "--toponly",
                                                                "bottom": "--bottomonly"}[strand_search],
                                                        cpu=cpu,
                                                        out=output_path.replace(' ', '\\ '),
                                                        tblout=tblout_path.replace(' ', '\\ '),
                                                        incE=evalue_attc,
                                                        mod_attc_path=model_attc_path.replace(' ', '\\ '),
                                                        infile=infile_path.replace(' ', '\\ '))
    try:
        _log.debug("run cmsearch: {}".format(cmsearch_cmd))
        with open(os.devnull, 'w') as dev_null:
            completed_process = subprocess.run(shlex.split(cmsearch_cmd), stdout=dev_null)
    except Exception as err:
        raise RuntimeError(f"{cmsearch_cmd} failed : {err}")
    if completed_process.returncode != 0:
        raise RuntimeError(f"{cmsearch_cmd} failed returncode = {completed_process.returncode}")
    df_max = read_infernal(tblout_path,
                           replicon.id,
                           replicon_size,
                           model_len(model_attc_path),
                           evalue=evalue_attc,
                           size_max_attc=max_attc_size,
                           size_min_attc=min_attc_size)
    # if replicon is linear
    # df_max.pos_beg + window_beg is always < replicon_size
    # (df_max.pos_beg + window_beg) % replicon_size = (df_max.pos_beg + window_beg)
    #
    # if replicon is circular and attc site overlap origin
    # df_max.pos_beg + window_beg  > replicon_size
    # (df_max.pos_beg + window_beg) % replicon_size is position on replicon
    # for instance with pos = 100 and replicon size = 90
    # 100 % 90 = 10
    if replicon.topology == 'circ':
        df_max.pos_beg = (df_max.pos_beg + window_beg) % replicon_size
        df_max.pos_end = (df_max.pos_end + window_beg) % replicon_size
    else:
        df_max.pos_beg = (df_max.pos_beg + window_beg).clip(lower=0, upper=replicon_size)
        df_max.pos_end = (df_max.pos_end + window_beg).clip(lower=0, upper=replicon_size)
    df_max.to_csv(os.path.join(out_dir, replicon.id + "_subseq_attc_table_end.res"),
                  sep="\t", index=0, mode="a", header=0)
    # filter on size
    df_max = df_max[(abs(df_max.pos_end - df_max.pos_beg) > min_attc_size) &
                    (abs(df_max.pos_end - df_max.pos_beg) < max_attc_size)]
    return df_max


def expand(replicon,
           window_beg, window_end, max_elt,
           circular, dist_threshold, model_attc_path,
           max_attc_size=200, min_attc_size=40, evalue_attc=1.,
           search_left=False, search_right=False,
           out_dir='.', cpu=1, cmsearch_bin='cmsearch'):
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
    :param int max_attc_size: The maximum value for the attC size
    :param int min_attc_size: The minimum value for the attC size
    :param str model_attc_path: the path to the attc model file
    :param float evalue_attc: evalue threshold to filter out hits above it
    :param bool search_left: trigger the local_max search on the left of the already detected element
    :param bool search_right: trigger the local_max search on the right of the already detected element
    :param str out_dir: The path to directory where to write results
    :param int cpu: the number of cpu use by expand
    :param str cmsearch_bin: the path to the `cmsearch` binary to use
    :return: a copy of max_elt with attC hits
    :rtype: :class:`pandas.DataFrame` object

    """
    replicon_size = len(replicon)
    # for a given element, we can search on the left hand side of it
    # (if the integrase is on the right and attC sites on the left for instance),
    # on the right hand side of it (opposite situation), or on both sides (only integrase or only attC sites)
    wb = window_beg
    we = window_end

    if search_right:
        # max_attc_size (200bo by default) to allow the detection of sites that would overlap 2 consecutive windows
        if circular:
            # to allow the detection of sites that overlap last window and the first one (given as argument)
            end_of_turn = wb + max_attc_size
            pass_through_ori = False
            window_beg = (window_end - max_attc_size) % replicon_size
            window_end = (window_end + dist_threshold) % replicon_size
        else:
            window_beg = max(0, window_end - max_attc_size)
            window_end = min(replicon_size, window_end + dist_threshold)

        searched_strand = "both" if search_left else "top"  # search on both strands if search in both directions
        while True:
            df_max = local_max(replicon,
                               window_beg, window_end,
                               model_attc_path,
                               max_attc_size=max_attc_size,
                               min_attc_size=min_attc_size,
                               strand_search=searched_strand,
                               out_dir=out_dir,
                               cpu=cpu,
                               evalue_attc=evalue_attc,
                               cmsearch_bin=cmsearch_bin
                               )
            to_concat = [df for df in [max_elt, df_max] if not df.empty]
            if to_concat:
                max_elt = pd.concat(to_concat)
            _log.info(f"\tsearched {window_beg}->{window_end} in {searched_strand} strand(s) found {len(df_max)} attc sites")

            if circular:
                if window_end == end_of_turn:
                    if searched_strand == "both":
                        # we loop all over the replicon and search in both strand
                        # searching left is useless
                        search_left = False
                    break
                elif df_max.empty:
                    break

                pass_through_ori = pass_through_ori or (window_end + dist_threshold) >= replicon_size

                window_beg = (window_end - max_attc_size) % replicon_size
                window_end = (window_end + dist_threshold) % replicon_size
                if we > wb:
                    if pass_through_ori:
                        window_end = min(end_of_turn, window_end)
                else:
                    window_end = min(end_of_turn, window_end)
            else:
                if df_max.empty or window_end == replicon_size:
                    break
                window_beg = max(0, window_end - max_attc_size)
                window_end = min(replicon_size, window_end + dist_threshold)

    if search_left:
        window_beg = wb
        if circular:
            #  to allow the detection of sites that overlap last window and the first one (given as argument)
            end_of_turn = we - max_attc_size
            pass_through_ori = False
            window_end = (window_beg + max_attc_size) % replicon_size
            window_beg = (window_beg - dist_threshold) % replicon_size
        else:
            window_end = min(replicon_size, window_beg + max_attc_size)
            window_beg = max(0, window_beg - dist_threshold)

        searched_strand = "both" if search_right else "bottom"
        while True:
            df_max = local_max(replicon,
                               window_beg, window_end,
                               model_attc_path,
                               max_attc_size=max_attc_size,
                               min_attc_size=min_attc_size,
                               strand_search=searched_strand,
                               out_dir=out_dir,
                               cpu=cpu,
                               evalue_attc=evalue_attc,
                               cmsearch_bin=cmsearch_bin)
            to_concat = [df for df in (max_elt, df_max) if not df.empty]
            if to_concat:
                max_elt = pd.concat(to_concat) # update of attC list of hits.
            _log.info(f"\tsearched {window_beg}->{window_end} in {searched_strand} strand(s) found {len(df_max)} attc sites")

            if circular:
                if df_max.empty or window_beg == end_of_turn:
                    break

                pass_through_ori = pass_through_ori or (window_end - dist_threshold) <= 0

                window_end = (window_beg + max_attc_size) % replicon_size
                window_beg = (window_beg - dist_threshold) % replicon_size

                if we > wb:
                    if pass_through_ori:
                        window_beg = max(window_beg, end_of_turn)
                else:
                    window_beg = max(window_beg, end_of_turn)
            else:
                if df_max.empty or window_beg == 0:
                    break
                window_end = min(replicon_size, window_beg + max_attc_size)
                window_beg = max(0, window_beg - dist_threshold)
    max_elt.drop_duplicates(inplace=True)
    max_elt.index = list(range(len(max_elt)))
    return max_elt
