import numpy as np
import pandas as pd
from Bio import SeqFeature
from Bio import SeqIO


def to_gbk(integron_desc, sequence, prot_file, dist_threshold):

    """ from a dataframe like integrons_describe and a sequence, create an genbank file with integron annotation """

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
            sequence.features.append(tmp)
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
                sequence.features.append(tmp)

            else:
                tmp = SeqFeature.SeqFeature(location=SeqFeature.FeatureLocation(integron_desc.loc[i].pos_beg - 1,
                                                                                integron_desc.loc[i].pos_end),
                                            strand=integron_desc.loc[i].strand,
                                            type=integron_desc.loc[i].type_elt,
                                            qualifiers={integron_desc.loc[i].type_elt: integron_desc.loc[i].element,
                                                        "model": integron_desc.loc[i].model}
                                            )

                sequence.features.append(tmp)

        else:
            type_integron = integron_desc.loc[i].type.values[0]
            # Should only be true if integron over edge of sequence:
            diff = integron_desc.loc[i].pos_beg.diff() > dist_threshold

            if diff.any():
                pos = np.where(diff)[0][0]
                start_integron_1 = integron_desc.loc[i].pos_beg.values[pos]
                end_integron_1 = len(sequence)
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
            sequence.features.append(tmp)
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
                    sequence.features.append(tmp)
                else:
                    tmp = SeqFeature.SeqFeature(location=SeqFeature.FeatureLocation(r[1].pos_beg - 1, r[1].pos_end),
                                                strand=r[1].strand,
                                                type=r[1].type_elt,
                                                qualifiers={r[1].type_elt: r[1].element, "model": r[1].model}
                                                )

                    sequence.features.append(tmp)

    # We get a ValueError otherwise, eg:
    # ValueError: Locus identifier 'gi|00000000|gb|XX123456.2|' is too long
    if len(sequence.name) > 16:
        sequence.name = sequence.name[-16:]
