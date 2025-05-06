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

import colorlog
import numpy as np
import pandas as pd

from .infernal import local_max, expand

_log = colorlog.getLogger(__name__)


def search_attc(attc_df, keep_palindromes, dist_threshold, replicon_size, rep_topology):
    """
    Parse the attc data set (sorted along start site) for the given replicon and return list of arrays.
    One array is composed of attC sites on the same strand and separated by a distance less than dist_threshold.

    :param attc_df:
    :type attc_df: :class:`pandas.DataFrame`
    :param bool keep_palindromes: True if the palindromes must be kept in attc result, False otherwise
    :param int dist_threshold: the maximal distance between 2 elements to aggregate them
    :param int replicon_size: the replicon number of base pair
    :param str rep_topology: the replicon topology should be 'lin' or 'circ'
    :return: a list attC sites found on replicon
    :rtype: list of :class:`pandas.DataFrame` objects
    """
    def overlap(attc1, attc2):
        """
        compare two attc sites (attc2 must start at smae position or after attc1)
        and remove attc sites which overlap more than 50%
        keep the one with lower evalue

        :param attc1: the first attc to compare
        :type attc1: pandas.Series with fields
                     "Accession_number", "cm_attC", "cm_debut", "cm_fin", "pos_beg", "pos_end", "sens", "evalue"
        :param attc2: the second attc to compare
        :type attc2: pandas.Series with fields
                     "Accession_number", "cm_attC", "cm_debut", "cm_fin", "pos_beg", "pos_end", "sens", "evalue"
        :return: a tuple of attc (pandas Series)
        """
        if attc2.pos_beg >= attc1.pos_end:
            return attc1, attc2
        else:
            attc_1_len = attc1.pos_end - attc1.pos_beg
            overlap = min(attc1.pos_end, attc2.pos_end) - attc2.pos_beg # by def attc_2.pos_beg >= attc_1.pos_beg
            ratio = overlap / attc_1_len
            if ratio > .5:
                return (attc1, ) if attc1.evalue < attc2.evalue else (attc2, )
            else:
                return attc1, attc2
    ok = False

    position_bkp_minus = []
    position_bkp_plus = []

    attc_plus = attc_df[attc_df.sens == "+"].copy()
    attc_minus = attc_df[attc_df.sens == "-"].copy()

    if not keep_palindromes:
        attc_df.sort_values(by=["pos_beg", "pos_end", "evalue"], inplace=True)
        if not attc_df.empty:
            attc_to_keep = [attc_df.iloc[0]]
            for row_idx in range(1, len(attc_df)):
                previous_attc = attc_to_keep.pop()
                attc_to_keep.extend(overlap(previous_attc, attc_df.iloc[row_idx]))
            attc_df = pd.DataFrame(attc_to_keep, index=list(range(len(attc_to_keep))))
        attc_plus = attc_df[attc_df.sens == "+"].copy()
        attc_minus = attc_df[attc_df.sens == "-"].copy()

    # can be reordered
    if (attc_plus.pos_beg.diff() > dist_threshold).any() or (attc_minus.pos_beg.diff() > dist_threshold).any():
        if not attc_plus.empty:
            bkp_plus = attc_plus[attc_plus.pos_beg.diff() > dist_threshold].index
            position_bkp_plus = [attc_plus.index.get_loc(i) for i in bkp_plus]
        if not attc_minus.empty:
            bkp_minus = attc_minus[(attc_minus.pos_beg.diff() > dist_threshold)].index
            position_bkp_minus = [attc_minus.index.get_loc(i) for i in bkp_minus]
        ok = True
    if not attc_plus.empty and not attc_minus.empty:
        ok = True

    if not ok:
        if attc_df.empty:
            attc_array = []
        else:
            attc_array = [attc_df]
    else:
        if attc_plus.empty:
            array_plus = []
        else:
            array_plus = np.split(attc_plus.values, position_bkp_plus)
            # array_plus is a list of np.array
            first_pos_beg = array_plus[0][0][4]
            last_pos_beg = array_plus[-1][-1][4]
            if len(array_plus) > 1 \
                    and rep_topology == 'circ' \
                    and (first_pos_beg - last_pos_beg) % replicon_size < dist_threshold:

                array_plus[0] = np.concatenate((array_plus[-1], array_plus[0]))
                del array_plus[-1]

        if attc_minus.empty:
            array_minus = []
        else:
            array_minus = np.split(attc_minus.values, position_bkp_minus)
            # array_minus is a list of np.array
            first_pos_beg = array_minus[0][0][4]
            last_pos_beg = array_minus[-1][-1][4]
            if len(array_minus) > 1 \
                    and rep_topology == 'circ' \
                    and (first_pos_beg - last_pos_beg) % replicon_size < dist_threshold:
                array_minus[0] = np.concatenate((array_minus[-1], array_minus[0]))
                del array_minus[-1]
            last_pos_beg = array_minus[-1][-1][4]

        tmp = array_plus + array_minus

        attc_array = [pd.DataFrame(i, columns=["Accession_number", "cm_attC", "cm_debut",
                                               "cm_fin", "pos_beg", "pos_end", "sens", "evalue"]) for i in tmp]
        # convert positions to int, and evalue to float
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        for a in attc_array:
            a[intcols] = a[intcols].astype(int)
            a["evalue"] = a["evalue"].astype(float)
    return attc_array


def find_attc_max(integrons, replicon, distance_threshold,
                  model_attc_path,
                  max_attc_size, min_attc_size,
                  evalue_attc=1.,
                  circular=True, out_dir='.',
                  cmsearch_bin='cmsearch',
                  cpu=1):
    """
    Look for attC site with cmsearch --max option which remove all heuristic filters.
    As this option make the algorithm way slower, we only run it in the region around a
    hit. We call it local_max or eagle_eyes.

    **Default hit**

    .. code-block:: text
    
                         attC
        __________________-->____-->_________-->_____________
        ______<--------______________________________________
                 intI
                      ^-------------------------------------^
                     Search-space with --local_max

    **Updated hit**

    .. code-block:: text

                         attC          ***         ***
        __________________-->____-->___-->___-->___-->_______
        ______<--------______________________________________
                 intI

    :param integrons: the integrons may contain or not attC or intI.
    :type integrons: list of :class:`Integron` objects.
    :param replicon: replicon where the integrons were found (genomic fasta file).
    :type replicon: :class:`Bio.Seq.SeqRecord` object.
    :param int distance_threshold: the maximal distance between 2 elements to aggregate them.
    :param float evalue_attc: evalue threshold to filter out hits above it.
    :param str model_attc_path: path to the attc model (Covariance Matrix).
    :param int max_attc_size: maximum value for the attC size.
    :param int min_attc_size: minimum value for the attC size.
    :param bool circular: True if replicon is circular, False otherwise.
    :param str out_dir: The directory where to write results
                        used indirectly by some called functions as :func:`infernal.local_max` or `infernal.expand`.
    :param str cmsearch_bin: The path to the `cmsearch_bin` binary to use
    :param int cpu: call local_max with the right number of cpu
    :return: a table of attC site
    :rtype: :class:`pd.DataFrame` object with monotonic indexes

    """
    size_replicon = len(replicon)
    columns = ['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg', 'pos_end', 'sens', 'evalue']
    data_type = {'Accession_number': 'str', 'cm_attC': 'str',
                 'cm_debut': 'int', 'cm_fin': 'int',
                 'pos_beg': 'int', 'pos_end': 'int',
                 'sens': 'str', 'evalue': 'float'
                 }

    max_final = pd.DataFrame(columns=columns)

    def merge_previous_attc_w_local_max(integron, local_max):
        """

        :param integron:
        :param local_max:
        :return:
        """
        previous_attc = pd.DataFrame({
            'Accession_number': [replicon.id] * len(integron.attC),
            'cm_attC': integron.attC.model,
            'cm_debut': [-1] * len(integron.attC),  # not use
            'cm_fin': [-1] * len(integron.attC),  # not use
            'pos_beg': integron.attC.pos_beg,
            'pos_end': integron.attC.pos_end,
            'sens': ['-' if strand == -1 else '+' for strand in integron.attC.strand],
            'evalue': integron.attC.evalue
            }
        )
        previous_attc = previous_attc.astype(dtype=data_type)

        all_attc = pd.concat((previous_attc, local_max))
        all_attc.sort_values(by=["pos_beg", "pos_end", "evalue"], inplace=True)
        attc_init = pd.DataFrame([all_attc.iloc[0]])
        # there is only one row in attc_init
        attc_init["pos_beg"] = attc_init["pos_beg"] - 10
        attc_init["pos_end"] = attc_init["pos_end"] - 10
        # I add a first row otherwise the first row is always discarded
        # as diff return NaN
        all_attc = pd.concat((attc_init, all_attc), ignore_index=True)
        all_attc = all_attc.loc[((all_attc["pos_beg"].diff().abs() > 3) |
                                 (all_attc["pos_end"].diff().abs() > 3)) &
                                (all_attc["sens"].eq(all_attc["sens"].shift())
                                 )].copy()
        return all_attc

    for i in integrons:
        max_elt = pd.DataFrame(columns=columns)
        max_elt = max_elt.astype(dtype=data_type)
        full_element = i.describe()  # dataframe
        # full_element structure
        # "pos_beg", "pos_end", "strand", "evalue", "type_elt", "model",
        # "distance_2attC", "annotation", "considered_topology"
        if all(full_element.type == "complete"):
            # Where is the integrase compared to the attc sites (no matter the strand) :
            integrase_is_left = ((full_element[full_element.type_elt == "attC"].pos_beg.values[0] -
                                  full_element[full_element.annotation == "intI"].pos_end.values[0]) % size_replicon <
                                 (full_element[full_element.annotation == "intI"].pos_beg.values[0] -
                                  full_element[full_element.type_elt == "attC"].pos_end.values[-1]) % size_replicon)

            if integrase_is_left:
                window_beg = full_element[full_element.annotation == "intI"].pos_end.values[0]
                distance_threshold_left = 0
                window_end = full_element[full_element.type_elt == "attC"].pos_end.values[-1]
                distance_threshold_right = distance_threshold

            else:  # is right
                window_beg = full_element[full_element.type_elt == "attC"].pos_beg.values[0]
                distance_threshold_left = distance_threshold
                window_end = full_element[full_element.annotation == "intI"].pos_end.values[-1]
                distance_threshold_right = 0

            if circular:
                window_beg = (window_beg - distance_threshold_left) % size_replicon
                window_end = (window_end + distance_threshold_right) % size_replicon
            else:
                window_beg = max(0, window_beg - distance_threshold_left)
                window_end = min(size_replicon, window_end + distance_threshold_right)

            strand = "top" if full_element[full_element.type_elt == "attC"].strand.values[0] == 1 else "bottom"
            df_max = local_max(replicon, window_beg, window_end, model_attc_path,
                               strand_search=strand,
                               evalue_attc=evalue_attc,
                               max_attc_size=max_attc_size,
                               min_attc_size=min_attc_size,
                               out_dir=out_dir,
                               cmsearch_bin=cmsearch_bin,
                               cpu=cpu)

            all_attc = merge_previous_attc_w_local_max(i, df_max)

            max_elt = pd.concat([df for df in (max_elt, all_attc) if not df.empty])

            # If we find new attC after the last found with default algo and if the integrase is on the left
            # (We don't expand over the integrase) :
            # pos_beg - pos_end so it's the same, the distance will always be > distance_threshold

            go_left = (full_element[full_element.type_elt == "attC"].pos_beg.values[0] - all_attc.pos_end.values[0]
                       ) % size_replicon < distance_threshold and not integrase_is_left
            go_right = (all_attc.pos_beg.values[-1] - full_element[full_element.type_elt == "attC"].pos_end.values[-1]
                        ) % size_replicon < distance_threshold and integrase_is_left
            max_elt = expand(replicon,
                             window_beg, window_end, max_elt,
                             circular, distance_threshold,
                             model_attc_path,
                             max_attc_size=max_attc_size,
                             min_attc_size=min_attc_size,
                             evalue_attc=evalue_attc,
                             search_left=go_left, search_right=go_right,
                             out_dir=out_dir,
                             cmsearch_bin=cmsearch_bin,
                             cpu=cpu)

        elif all(full_element.type == "CALIN"):
            if full_element[full_element.pos_beg.isin(max_final.pos_beg)].empty:
                # if cluster don't overlap already max-searched region
                window_beg = full_element[full_element.type_elt == "attC"].pos_beg.values[0]
                window_end = full_element[full_element.type_elt == "attC"].pos_end.values[-1]
                if circular:
                    window_beg = (window_beg - distance_threshold) % size_replicon
                    window_end = (window_end + distance_threshold) % size_replicon
                else:
                    window_beg = max(0, window_beg - distance_threshold)
                    window_end = min(size_replicon, window_end + distance_threshold)
                strand = "top" if full_element[full_element.type_elt == "attC"].strand.values[0] == 1 else "bottom"

                df_max = local_max(replicon, window_beg, window_end,
                                   model_attc_path,
                                   strand_search=strand,
                                   evalue_attc=evalue_attc,
                                   max_attc_size=max_attc_size,
                                   min_attc_size=min_attc_size,
                                   out_dir=out_dir,
                                   cmsearch_bin=cmsearch_bin,
                                   cpu=cpu)

                all_attc = merge_previous_attc_w_local_max(i, df_max)
                to_concat = [df for df in [max_elt, all_attc] if not df.empty]
                if to_concat:
                    max_elt = pd.concat(to_concat)

                if not df_max.empty:  # Max can sometimes find bigger attC than permitted
                    go_left = (full_element[full_element.type_elt == "attC"].pos_beg.values[0] - all_attc.pos_end.values[0]
                               ) % size_replicon < distance_threshold
                    go_right = (all_attc.pos_beg.values[-1] - full_element[full_element.type_elt == "attC"].pos_end.values[-1]
                                ) % size_replicon < distance_threshold
                    max_elt = expand(replicon,
                                     window_beg, window_end, max_elt,
                                     circular, distance_threshold,
                                     model_attc_path,
                                     max_attc_size=max_attc_size,
                                     min_attc_size=min_attc_size,
                                     evalue_attc=evalue_attc,
                                     search_left=go_left, search_right=go_right,
                                     out_dir=out_dir,
                                     cmsearch_bin=cmsearch_bin,
                                     cpu=cpu)

        elif all(full_element.type == "In0"):
            if all(full_element.model != "Phage_integrase"):
                window_beg = full_element[full_element.annotation == "intI"].pos_beg.values[0]
                window_end = full_element[full_element.annotation == "intI"].pos_end.values[-1]
                if circular:
                    window_beg = (window_beg - distance_threshold) % size_replicon
                    window_end = (window_end + distance_threshold) % size_replicon
                else:
                    window_beg = max(0, window_beg - distance_threshold)
                    window_end = min(size_replicon, window_end + distance_threshold)
                df_max = local_max(replicon,
                                   window_beg, window_end,
                                   model_attc_path,
                                   evalue_attc=evalue_attc,
                                   max_attc_size=max_attc_size,
                                   min_attc_size=min_attc_size,
                                   out_dir=out_dir,
                                   cmsearch_bin=cmsearch_bin,
                                   cpu=cpu)
                to_concat = [df for df in [max_elt, df_max] if not df.empty]
                if to_concat:
                    max_elt = pd.concat(to_concat)
                if not max_elt.empty:
                    max_elt = expand(replicon,
                                     window_beg, window_end, max_elt,
                                     circular, distance_threshold,
                                     model_attc_path,
                                     max_attc_size=max_attc_size,
                                     min_attc_size=min_attc_size,
                                     evalue_attc=evalue_attc,
                                     search_left=True, search_right=True,
                                     out_dir=out_dir,
                                     cmsearch_bin=cmsearch_bin,
                                     cpu=cpu)
        to_concat = [df for df in [max_final, max_elt] if not df.empty]
        if to_concat:
            max_final = pd.concat(to_concat, ignore_index=True)
        max_final.drop_duplicates(subset=max_final.columns[:-1], inplace=True)
        max_final.index = range(len(max_final)) # this line is important for next step the index must be monotonic
    max_final = max_final.astype(dtype=data_type)
    return max_final
