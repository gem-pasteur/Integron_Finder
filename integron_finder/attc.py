import os
from subprocess import call

import numpy as np
import pandas as pd

from .infernal import local_max, expand

def search_attc(attc_df, keep_palindromes, dist_threshold, replicon_size):
    """
    Parse the attc dataset (sorted along start site) for the given replicon and return list of arrays.
    One array is composed of attC sites on the same strand and separated by a
    distance less than 5kb
    """
    ok = False

    position_bkp_minus = []
    position_bkp_plus = []

    attc_plus = attc_df[attc_df.sens == "+"].copy()
    attc_minus = attc_df[attc_df.sens == "-"].copy()

    if not keep_palindromes:
        attc_df = attc_df.sort_values(["pos_beg", "evalue"]).drop_duplicates(subset=["pos_beg"]).copy()
        attc_plus = attc_df[attc_df.sens == "+"].copy()
        attc_minus = attc_df[attc_df.sens == "-"].copy()

    # can be reordered
    if (attc_plus.pos_beg.diff() > dist_threshold).any() or (attc_minus.pos_beg.diff() > dist_threshold).any():
        if len(attc_plus) > 0:
            bkp_plus = attc_plus[(attc_plus.pos_beg.diff() > dist_threshold)].index
            position_bkp_plus = [attc_plus.index.get_loc(i) for i in bkp_plus]
        if len(attc_minus) > 0:
            bkp_minus = attc_minus[(attc_minus.pos_beg.diff() > dist_threshold)].index
            position_bkp_minus = [attc_minus.index.get_loc(i) for i in bkp_minus]
        ok = True

    if len(attc_plus) > 0 and len(attc_minus) > 0:
        ok = True

    if not ok:
        if len(attc_df) == 0:
            return []
        else:
            return [attc_df]
    else:
        if len(attc_plus) > 0:
            array_plus = np.split(attc_plus.values, position_bkp_plus)
            if (array_plus[0][0][4] - array_plus[-1][-1][4]) % replicon_size < dist_threshold and len(array_plus) > 1:
                array_plus[0] = np.concatenate((array_plus[-1], array_plus[0]))
                del array_plus[-1]
        else:
            array_plus = np.array([])
        if len(attc_minus) > 0:
            array_minus = np.split(attc_minus.values, position_bkp_minus)
            if (array_minus[0][0][4]-array_minus[-1][-1][4]) % replicon_size < dist_threshold and len(array_minus) > 1:
                array_minus[0] = np.concatenate((array_minus[-1], array_minus[0]))
                del array_minus[-1]
        else:
            array_minus = np.array([])

        if len(array_minus) > 0 and len(array_plus) > 0:
            tmp = array_plus + array_minus
        elif len(array_minus) == 0:
            tmp = array_plus
        elif len(array_plus) == 0:
            tmp = array_minus

        attc_array = [pd.DataFrame(i, columns=["Accession_number", "cm_attC", "cm_debut",
                                               "cm_fin", "pos_beg", "pos_end", "sens", "evalue"]) for i in tmp]
        # convert positions to int, and evalue to float
        intcols = ["cm_debut", "cm_fin", "pos_beg", "pos_end"]
        for a in attc_array:
            a[intcols] = a[intcols].astype(int)
            a["evalue"] = a["evalue"].astype(float)
        return attc_array


def find_attc(replicon_path, replicon_name, cmsearch_path, out_dir, model_attc, cpu=1):
    """
    Call cmsearch to find attC sites in a single replicon.

    :param replicon_path: the path of the replicon to analyse
    :type replicon_path: string
    :param replicon_name: the name of the replicon to analyse
    :type replicon_name: string
    :param out_dir: the relative path to the directory where cmsearch outputs will be stored
    :type out_dir: str
    :returns: None, the results are written on the disk
    :raises RuntimeError: when cmsearch run failed
    """
    cmsearch_cmd = [cmsearch_path,
                    "--cpu", str(cpu),
                    "-o", os.path.join(out_dir, replicon_name + "_attc.res"),
                    "--tblout", os.path.join(out_dir, replicon_name + "_attc_table.res"),
                    "-E", "10",
                    model_attc,
                    replicon_path]
    try:
        returncode = call(cmsearch_cmd)
    except Exception as err:
        raise RuntimeError("{0} failed : {1}".format(' '.join(cmsearch_cmd), err))
    if returncode != 0:
        raise RuntimeError("{0} failed returncode = {1}".format(' '.join(cmsearch_cmd), returncode))


def find_attc_max(integrons, replicon, distance_threshold,
                  model_attc_path, max_attc_size, circular=True, outfile="attC_max_1.res", out_dir='.'):
    """
    Look for attC site with cmsearch --max option wich remove all heuristic filters.
    As this option make the algorithm way slower, we only run it in the region around a
    hit. We call it local_max or eagle_eyes.


    Default hit :
    =============
                     attC
    __________________-->____-->_________-->_____________
    ______<--------______________________________________
             intI
                  ^-------------------------------------^
                 Search-space with --local_max

    Updated hit :
    =============

                     attC          ***         ***
    __________________-->____-->___-->___-->___-->_______
    ______<--------______________________________________
             intI

    :param integrons: the integrons may contain or not attC or intI.
    :type integrons: list of :class:`Integron` objects.
    :param bool circular: True if replicon is circular, False otherwise.
    :param outfile: the name of cmsearch result file
    :type outfile: string
    :return:
    :rtype: :class:`pd.DataFrame`
    """
    size_replicon = len(replicon)
    columns = ['Accession_number', 'cm_attC', 'cm_debut', 'cm_fin', 'pos_beg', 'pos_end', 'sens', 'evalue']
    data_type = {'Accession_number': 'str', 'cm_attC': 'str',
                 'cm_debut': 'int', 'cm_fin': 'int',
                 'pos_beg': 'int', 'pos_end': 'int', }
    max_final = pd.DataFrame(columns=columns)
    max_final = max_final.astype(dtype=data_type)
    for i in integrons:
        max_elt = pd.DataFrame(columns=columns)
        max_elt = max_elt.astype(dtype=data_type)
        full_element = i.describe()

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
            df_max = local_max(replicon, window_beg, window_end, model_attc_path, strand_search=strand, out_dir=out_dir)
            max_elt = pd.concat([max_elt, df_max])

            # If we find new attC after the last found with default algo and if the integrase is on the left
            # (We don't expand over the integrase) :
            # pos_beg - pos_end so it's the same, the distance will always be > distance_threshold

            go_left = (full_element[full_element.type_elt == "attC"].pos_beg.values[0] - df_max.pos_end.values[0]
                       ) % size_replicon < distance_threshold and not integrase_is_left
            go_right = (df_max.pos_beg.values[-1] - full_element[full_element.type_elt == "attC"].pos_end.values[-1]
                        ) % size_replicon < distance_threshold and integrase_is_left
            max_elt = expand(replicon,
                             window_beg, window_end, max_elt, df_max,
                             circular, distance_threshold, max_attc_size,
                             model_attc_path,
                             search_left=go_left, search_right=go_right)

        elif all(full_element.type == "CALIN"):
            if len(full_element[full_element.pos_beg.isin(max_final.pos_beg)]) == 0:
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
                                   out_dir=out_dir)
                max_elt = pd.concat([max_elt, df_max])

                if len(df_max) > 0:  # Max can sometimes find bigger attC than permitted
                    go_left = (full_element[full_element.type_elt == "attC"].pos_beg.values[0] - df_max.pos_end.values[0]
                               ) % size_replicon < distance_threshold
                    go_right = (df_max.pos_beg.values[-1] - full_element[full_element.type_elt == "attC"].pos_end.values[-1]
                                ) % size_replicon < distance_threshold
                    max_elt = expand(replicon,
                                     window_beg, window_end, max_elt, df_max,
                                     circular, distance_threshold, max_attc_size,
                                     model_attc_path,
                                     search_left=go_left, search_right=go_right)

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
                                   strand_search=strand,
                                   out_dir=out_dir)
                max_elt = pd.concat([max_elt, df_max])
                if len(max_elt) > 0:
                    max_elt = expand(replicon,
                                     window_beg, window_end, max_elt, df_max,
                                     circular, distance_threshold, max_attc_size,
                                     model_attc_path,
                                     search_left=True, search_right=True)

        max_final = pd.concat([max_final, max_elt])
        max_final.drop_duplicates(subset=max_final.columns[:-1], inplace=True)
        max_final.index = range(len(max_final))
    return max_final
