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

import numpy as np
import pandas as pd
from Bio import SeqFeature
from Bio import SeqIO

from .utils import get_name_from_path
from .hmm import read_hmm

_log = colorlog.getLogger(__name__)


def func_annot(integrons, replicon, prot_db, hmm_files, cfg, out_dir='.', evalue=10, coverage=0.5):
    """
    | Call hmmmer to annotate CDS associated with the integron.
    | Use NCBIfam-AMRFinder.LIB available here https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/
    | check for the appropriate version

    :param integrons: integrons list to annotate
    :type integrons: list of :class:`integron_finder.integron.Integron` objects.
    :param replicon: replicon where the integrons were found (genomic fasta file)
    :type replicon: :class:`Bio.Seq.SeqRecord` object
    :param prot_db: the protein database corresponding to the replicon translation
    :type prot_db: :class:`integron.prot_db.ProteinDB` object.
    :param hmm_files: List of path of hmm profiles to use to scan the prot_file
    :type hmm_files: List[str]
    :param cfg: the configuration for this analyse
    :type cfg: :class:`integron_finder.config.Config`
    :param str out_dir: the path of the directory where to store the results
    :param float evalue:
    :param float coverage:
    :return: None.

             But several files per hmm file are produced.

             * subseqprot.tmp: fasta file containing a subset of protfile (the proteins belonging to the integron)
             * <hmm>_fa.res: an output of the hmm search.
             * <hmm>_fa_table.res: an output of the hmm search in tabulated format.

    """

    prot_tmp = f'{os.path.join(out_dir, replicon.id + "_subseqprot.tmp")}'

    for integron in integrons:
        if os.path.isfile(prot_tmp):
            os.remove(prot_tmp)

        if integron.type() != "In0" and not integron.proteins.empty:

            func_annotate_res = pd.DataFrame(columns=["Accession_number",
                                                      "query_name", "ID_query",
                                                      "ID_prot", "strand",
                                                      "pos_beg", "pos_end", "evalue"])

            prot_to_annotate = []
            for prot_nb, prot_id in enumerate(prot_db.coding_prot_ids(), 1):
                if prot_id in integron.proteins.index:
                    prot_to_annotate.append(prot_db[prot_id])

            SeqIO.write(prot_to_annotate, prot_tmp, "fasta")
            for hmm in hmm_files:
                name_wo_ext = "{}_{}".format(replicon.id, get_name_from_path(hmm))
                hmm_out = os.path.join(out_dir, "{}_fa.res".format(name_wo_ext))
                hmm_tableout = os.path.join(out_dir, "{}_fa_table.res".format(name_wo_ext))
                hmm_cmd = "{hmmsearch} --cut_ga -Z {prot_nb} --cpu {cpu} --tblout {tblout} -o {hmm_out}" \
                          " {hmm} {prot_tmp}".format(
                    hmmsearch=cfg.hmmsearch.replace(' ', '\\ '),
                    prot_nb=prot_nb,
                    cpu=cfg.cpu,
                    tblout=hmm_tableout.replace(' ', '\\ '),
                    hmm_out=hmm_out.replace(' ', '\\ '),
                    hmm=hmm.replace(' ', '\\ '),
                    prot_tmp=prot_tmp.replace(' ', '\\ ')
                )

                try:
                    _log.debug(f"run hmmsearch: {hmm_cmd}")
                    cmd = shlex.split(hmm_cmd)
                    completed_process = subprocess.run(cmd)
                except Exception as err:
                    raise RuntimeError(f"{cmd} failed : {err}")
                if completed_process.returncode != 0:
                    raise RuntimeError(f"{cmd} failed return code = {completed_process.returncode}")
                hmm_in = read_hmm(replicon.id, prot_db, hmm_out, cfg, evalue=evalue, coverage=coverage
                                  ).sort_values("evalue").drop_duplicates(subset="ID_prot")
                to_concat = [df for df in [func_annotate_res, hmm_in] if not df.empty]
                if to_concat:
                    func_annotate_res = pd.concat(to_concat)
            func_annotate_res = func_annotate_res.sort_values("evalue").drop_duplicates(subset="ID_prot")

            integron.proteins.loc[func_annotate_res.ID_prot, "evalue"] = func_annotate_res.evalue.values
            integron.proteins.loc[func_annotate_res.ID_prot, "annotation"] = func_annotate_res.query_name.values
            integron.proteins.loc[func_annotate_res.ID_prot, "model"] = func_annotate_res.ID_query.values
            integron.proteins = integron.proteins.astype(dtype=integron.dtype)


def add_feature(replicon, integron_desc, prot_db, dist_threshold):
    """
    Add integron annotation to the replicon.

    :param replicon: The Replicon to annotate
    :type replicon: a :class:`Bio.Seq.SeqRecord` object.
    :param integron_desc: integron description
    :type integron_desc: a :class:`pandas.DataFrame`
    :param prot_db: the path to the fasta file containing the translation of the replicon.
    :type prot_db: a :class:`integron_finder.prot_db.ProteinDB` object.
    :param int dist_threshold: Two elements are aggregated if they are distant of dist_threshold or less.
    """
    integron_desc = integron_desc.set_index("ID_integron").copy()
    for i in integron_desc.index.unique():

        if isinstance(integron_desc.loc[i], pd.Series):
            type_integron = integron_desc.loc[i].type
            start_integron = integron_desc.loc[i].pos_beg
            end_integron = integron_desc.loc[i].pos_end
            tmp = SeqFeature.SeqFeature(location=SeqFeature.SimpleLocation(int(start_integron) - 1,
                                                                           int(end_integron),
                                                                           strand=0),
                                        type="integron",
                                        qualifiers={"integron_id": i, "integron_type": type_integron}
                                        )
            replicon.features.append(tmp)
            if integron_desc.loc[i].type_elt == "protein":

                tmp = SeqFeature.SeqFeature(location=SeqFeature.SimpleLocation(int(integron_desc.loc[i].pos_beg) - 1,
                                                                               int(integron_desc.loc[i].pos_end),
                                                                               strand=integron_desc.loc[i].strand),
                                            type="CDS" if integron_desc.loc[i].annotation != "intI" else "integrase",
                                            qualifiers={"protein_id": integron_desc.loc[i].element,
                                                        "gene": integron_desc.loc[i].annotation,
                                                        "model": integron_desc.loc[i].model}
                                            )
                tmp.qualifiers["translation"] = [prot_db[prt_id] for prt_id in prot_db
                                                 if prt_id == integron_desc.loc[i].element][0].seq
                replicon.features.append(tmp)

            else:
                tmp = SeqFeature.SeqFeature(location=SeqFeature.SimpleLocation(int(integron_desc.loc[i].pos_beg) - 1,
                                                                               int(integron_desc.loc[i].pos_end),
                                                                               strand=integron_desc.loc[i].strand),
                                            type=integron_desc.loc[i].type_elt,
                                            qualifiers={integron_desc.loc[i].type_elt: integron_desc.loc[i].element,
                                                        "model": integron_desc.loc[i].model}
                                            )

                replicon.features.append(tmp)

        else:
            type_integron = integron_desc.loc[i].type.values[0]
            # Should only be true if integron over edge of replicon:
            diff = integron_desc.loc[i].pos_beg.diff() > dist_threshold

            if diff.any():
                pos = np.where(diff)[0][0]
                start_integron_1 = int(integron_desc.loc[i].pos_beg.values[pos])
                end_integron_1 = len(replicon)
                start_integron_2 = 1
                end_integron_2 = int(integron_desc.loc[i].pos_end.values[pos-1])

                f1 = SeqFeature.SimpleLocation(start_integron_1 - 1, end_integron_1, strand=0)
                f2 = SeqFeature.SimpleLocation(start_integron_2 - 1, end_integron_2, strand=0)
                tmp = SeqFeature.SeqFeature(location=f1 + f2,
                                            type="integron",
                                            qualifiers={"integron_id": i, "integron_type": type_integron}
                                            )
            else:
                start_integron = int(integron_desc.loc[i].pos_beg.values[0])
                end_integron = int(integron_desc.loc[i].pos_end.values[-1])

                tmp = SeqFeature.SeqFeature(location=SeqFeature.SimpleLocation(start_integron - 1,
                                                                               end_integron,
                                                                               strand=0),
                                            type="integron",
                                            qualifiers={"integron_id": i, "integron_type": type_integron}
                                            )
            replicon.features.append(tmp)
            for r in integron_desc.loc[i].iterrows():
                if r[1].type_elt == "protein":
                    tmp = SeqFeature.SeqFeature(location=SeqFeature.SimpleLocation(int(r[1].pos_beg) - 1,
                                                                                   int(r[1].pos_end),
                                                                                   strand=r[1].strand),
                                                type="CDS" if r[1].annotation != "intI" else "integrase",
                                                qualifiers={"protein_id": r[1].element,
                                                            "gene": r[1].annotation,
                                                            "model": r[1].model}
                                                )
                    tmp.qualifiers["translation"] = [prot_db[prt_id] for prt_id in prot_db
                                                     if prt_id == r[1].element][0].seq
                    replicon.features.append(tmp)
                else:
                    tmp = SeqFeature.SeqFeature(location=SeqFeature.SimpleLocation(int(r[1].pos_beg) - 1,
                                                                                   int(r[1].pos_end),
                                                                                   strand=r[1].strand),
                                                type=r[1].type_elt,
                                                qualifiers={r[1].type_elt: r[1].element, "model": r[1].model}
                                                )

                    replicon.features.append(tmp)

    # We get a ValueError otherwise, eg:
    # ValueError: Locus identifier 'gi|00000000|gb|XX123456.2|' is too long
    if len(replicon.name) > 16:
        replicon.name = replicon.name[-16:]
