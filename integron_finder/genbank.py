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

import numpy as np
import pandas as pd
from Bio import SeqFeature
from Bio import SeqIO


def add_feature(replicon, integron_desc, prot_file, dist_threshold):
    """
    Add integron annotation to the replicon.

    :param replicon: The Replicon to annotate
    :type replicon: a :class:`Bio.Seq.SeqRecord` object.
    :param integron_desc:
    :type integron_desc: a :class:`pandas.DataFrame`
    :param prot_file: the path to the fasta file containing the traduction of the replicon.
    :param int dist_threshold: Two elements are aggregated if they are distant of dist_threshold or less.
    """

    integron_desc = integron_desc.set_index("ID_integron").copy()
    for i in integron_desc.index.unique():

        if isinstance(integron_desc.loc[i], pd.Series):
            type_integron = integron_desc.loc[i].type
            start_integron = integron_desc.loc[i].pos_beg
            end_integron = integron_desc.loc[i].pos_end
            tmp = SeqFeature.SeqFeature(location=SeqFeature.FeatureLocation(start_integron - 1, end_integron),
                                        strand=0,
                                        type="integron",
                                        qualifiers={"integron_id": i, "integron_type": type_integron}
                                        )
            replicon.features.append(tmp)
            if integron_desc.loc[i].type_elt == "protein":

                tmp = SeqFeature.SeqFeature(location=SeqFeature.FeatureLocation(integron_desc.loc[i].pos_beg - 1,
                                                                                integron_desc.loc[i].pos_end),
                                            strand=integron_desc.loc[i].strand,
                                            type="CDS" if integron_desc.loc[i].annotation != "intI" else "integrase",
                                            qualifiers={"protein_id": integron_desc.loc[i].element,
                                                        "gene": integron_desc.loc[i].annotation,
                                                        "model": integron_desc.loc[i].model}
                                            )

                tmp.qualifiers["translation"] = [prt for prt in SeqIO.parse(prot_file, "fasta")
                                                 if prt.id == integron_desc.loc[i].element][0].seq
                replicon.features.append(tmp)

            else:
                tmp = SeqFeature.SeqFeature(location=SeqFeature.FeatureLocation(integron_desc.loc[i].pos_beg - 1,
                                                                                integron_desc.loc[i].pos_end),
                                            strand=integron_desc.loc[i].strand,
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
                start_integron_1 = integron_desc.loc[i].pos_beg.values[pos]
                end_integron_1 = len(replicon)
                start_integron_2 = 1
                end_integron_2 = integron_desc.loc[i].pos_end.values[pos-1]

                f1 = SeqFeature.FeatureLocation(start_integron_1 - 1, end_integron_1)
                f2 = SeqFeature.FeatureLocation(start_integron_2 - 1, end_integron_2)
                tmp = SeqFeature.SeqFeature(location=f1 + f2,
                                            strand=0,
                                            type="integron",
                                            qualifiers={"integron_id": i, "integron_type": type_integron}
                                            )
            else:
                start_integron = integron_desc.loc[i].pos_beg.values[0]
                end_integron = integron_desc.loc[i].pos_end.values[-1]

                tmp = SeqFeature.SeqFeature(location=SeqFeature.FeatureLocation(start_integron - 1, end_integron),
                                            strand=0,
                                            type="integron",
                                            qualifiers={"integron_id": i, "integron_type": type_integron}
                                            )
            replicon.features.append(tmp)
            for r in integron_desc.loc[i].iterrows():
                if r[1].type_elt == "protein":
                    tmp = SeqFeature.SeqFeature(location=SeqFeature.FeatureLocation(r[1].pos_beg - 1, r[1].pos_end),
                                                strand=r[1].strand,
                                                type="CDS" if r[1].annotation != "intI" else "integrase",
                                                qualifiers={"protein_id": r[1].element,
                                                            "gene": r[1].annotation,
                                                            "model": r[1].model}
                                                )

                    tmp.qualifiers["translation"] = [prt for prt in SeqIO.parse(prot_file, "fasta")
                                                     if prt.id == r[1].element][0].seq
                    replicon.features.append(tmp)
                else:
                    tmp = SeqFeature.SeqFeature(location=SeqFeature.FeatureLocation(r[1].pos_beg - 1, r[1].pos_end),
                                                strand=r[1].strand,
                                                type=r[1].type_elt,
                                                qualifiers={r[1].type_elt: r[1].element, "model": r[1].model}
                                                )

                    replicon.features.append(tmp)

    # We get a ValueError otherwise, eg:
    # ValueError: Locus identifier 'gi|00000000|gb|XX123456.2|' is too long
    if len(replicon.name) > 16:
        replicon.name = replicon.name[-16:]
