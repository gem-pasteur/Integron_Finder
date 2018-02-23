import os
from subprocess import call
import pandas as pd
from Bio import SeqIO

from integron_finder.utils import model_len


def read_infernal(infile, replicon_name, len_model_attc,
                  evalue=1, size_max_attc=200, size_min_attc=40):
    """
    Function that parse cmsearch --tblout output and returns a pandas DataFrame
    """

    try:
        _ = pd.read_table(infile, comment="#")
    except:
        return pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                     "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])

    df = pd.read_table(infile, sep="\s*", engine="python",  header=None, skipfooter=10, skiprows=2)
    # Keep only columns: query_name(2), mdl from(5), mdl to(6), seq from(7),
    # seq to(8), strand(9), E-value(15)
    df = df[[2, 5, 6, 7, 8, 9, 15]]
    df = df[(df[15] < evalue)]  # filter on evalue
    df = df[(abs(df[8] - df[7]) < size_max_attc) & (size_min_attc < abs(df[8] - df[7]))]
    if len(df) > 0:
        df["Accession_number"] = replicon_name
        c = df.columns.tolist()
        df = df[c[-1:] + c[:-1]]
        df.sort_values([8, 15], inplace=True)
        df.index = range(0, len(df))
        df.columns = ["Accession_number", "cm_attC", "cm_debut", "cm_fin",
                      "pos_beg_tmp", "pos_end_tmp",
                      "sens", "evalue"]
        idx = (df.pos_beg_tmp > df.pos_end_tmp)
        df.loc[idx, "pos_beg"] = df.loc[idx].apply(lambda x: x["pos_end_tmp"] - (len_model_attc - x["cm_fin"]), axis=1)
        df.loc[idx, "pos_end"] = df.loc[idx].apply(lambda x: x["pos_beg_tmp"] + (x["cm_debut"] - 1), axis=1)

        df.loc[~idx, "pos_end"] = df.loc[~idx].apply(lambda x: x["pos_end_tmp"] + (len_model_attc - x["cm_fin"]), axis=1)
        df.loc[~idx, "pos_beg"] = df.loc[~idx].apply(lambda x: x["pos_beg_tmp"] - (x["cm_debut"] - 1), axis=1)

        return df[["Accession_number", "cm_attC", "cm_debut", "cm_fin",
                      "pos_beg", "pos_end",
                      "sens", "evalue"]]
    else:
        return pd.DataFrame(columns=["Accession_number", "cm_attC", "cm_debut",
                                     "cm_fin", "pos_beg", "pos_end", "sens", "evalue"])


def local_max(replicon,
              window_beg, window_end,
              model_attc_path,
              strand_search="both",
              evalue_attc=1, max_attc_size=200, min_attc_size=40,
              cmsearch_bin='cmsearch', out_dir='.', cpu_nb=1):
    """

    :param replicon_name: the name of replicon (without suffix)
    :type replicon_name: str
    :param window_beg:
    :type window_beg: int
    :param window_end:
    :type window_end: int
    :param strand_search:
    :type strand_search: str
    :return:
    :rtype: :class:`pd.DataFrame` object
    """
    replicon_size = len(replicon)
    if window_beg < window_end:
        subseq = replicon[window_beg:window_end]
    else:
        subseq1 = replicon[window_beg:replicon_size]
        subseq2 = replicon[:window_end]
        subseq = subseq1 + subseq2

    with open(os.path.join(out_dir, replicon.name + "_subseq.fst"), "w") as f:
        SeqIO.write(subseq, f, "fasta")

    output_path = os.path.join(out_dir, "{name}_{win_beg}_{win_end}_subseq_attc.res".format(name=replicon.name,
                                                                                            win_beg=window_beg,
                                                                                            win_end=window_end))
    tblout_path = os.path.join(out_dir, "{name}_{win_beg}_{win_end}_subseq_attc_table.res".format(name=replicon.name,
                                                                                                  win_beg=window_beg,
                                                                                                  win_end=window_end))

    infile_path = os.path.join(out_dir, replicon.name + "_subseq.fst")
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
        returncode = call(cmsearch_cmd.split())
    except Exception as err:
        raise RuntimeError("{0} failed : {1}".format(cmsearch_cmd, err))
    if returncode != 0:
        raise RuntimeError("{0} failed returncode = {1}".format(cmsearch_cmd, returncode))

    df_max = read_infernal(tblout_path,
                           replicon.name, model_len(model_attc_path),
                           evalue=evalue_attc,
                           size_max_attc=max_attc_size,
                           size_min_attc=min_attc_size)

    df_max.pos_beg = (df_max.pos_beg + window_beg) % replicon_size
    df_max.pos_end = (df_max.pos_end + window_beg) % replicon_size
    df_max.to_csv(os.path.join(out_dir, replicon.name + "_subseq_attc_table_end.res"),
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
           out_dir='.'):
    """
    for a given element, we can search on the left hand side (if integrase is on the right for instance)
    or right hand side (opposite situation) or both side (only integrase or only attC sites)

    :param window_beg:
    :type window_beg: int
    :param window_end:
    :type window_end: int
    :param max_elt:
    :type max_elt: int
    :param df_max:
    :type df_max: :class:`pandas.DataFrame` object
    :param search_left: need to search on right of ???
    :type search_left: bool
    :param search_right: need to search on right of ???
    :type search_right: bool
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

        while len(df_max) > 0 and 0 < (window_beg and window_end) < replicon_size:

            df_max = local_max(replicon,
                               window_beg, window_end,
                               model_attc_path,
                               strand_search=searched_strand,
                               out_dir=out_dir
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

        while len(df_max) > 0 and 0 < (window_beg and window_end) < replicon_size:

            df_max = local_max(replicon,
                               window_beg, window_end,
                               model_attc_path,
                               strand_search=searched_strand,
                               out_dir=out_dir)
            max_elt = pd.concat([max_elt, df_max])  # update of attC list of hits.

            if circular:
                window_end = (window_beg + 200) % replicon_size
                window_beg = (window_beg - dist_threshold) % replicon_size

            else:
                window_end = min(replicon_size, window_beg + 200)
                window_beg = max(0, window_beg - dist_threshold)

    max_elt.drop_duplicates(inplace=True)
    max_elt.index = range(len(max_elt))
    return max_elt


