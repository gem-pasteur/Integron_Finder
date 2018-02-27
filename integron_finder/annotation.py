import os
from subprocess import call

from Bio import BiopythonExperimentalWarning
import warnings
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

from Bio import SeqIO
import pandas as pd

from .utils import get_name_from_path, read_multi_prot_fasta
from .hmm import read_hmm


def func_annot(integrons, replicon, prot_file, hmm_files, cfg, out_dir='.', evalue=10, coverage=0.5):
    """
    Call hmmmer to annotate CDS associated with the integron. Use Resfams per default (Gibson et al, ISME J.,  2014)
    """
    print "# Start Functional annotation... : "
    prot_tmp = os.path.join(out_dir, replicon.name + "_subseqprot.tmp")

    for integron in integrons:
        if os.path.isfile(prot_tmp):
            os.remove(prot_tmp)

        if integron.type() != "In0" and not integron.proteins.empty:

            func_annotate_res = pd.DataFrame(columns=["Accession_number",
                                                      "query_name", "ID_query",
                                                      "ID_prot", "strand",
                                                      "pos_beg", "pos_end", "evalue"])

            prot_to_annotate = []
            # It's protein file, fasta_reader is dedicated fr dna
            all_prots = read_multi_prot_fasta(prot_file)

            for prot_nb, prot in enumerate(all_prots, 1):
                if prot.id in integron.proteins.index:
                    prot_to_annotate.append(prot)

            SeqIO.write(prot_to_annotate, prot_tmp, "fasta")
            for hmm in hmm_files:
                name_wo_ext = "{}_{}".format(replicon.name, get_name_from_path(hmm))
                hmm_out = os.path.join(out_dir, "{}_fa.res".format(name_wo_ext))
                hmm_tableout = os.path.join(out_dir, "{}_fa_table.res".format(name_wo_ext))
                hmm_cmd = [cfg.hmmsearch,
                            "-Z", str(prot_nb),
                            "--cpu", str(cfg.cpu),
                            "--tblout", hmm_tableout,
                            "-o", hmm_out,
                            hmm,
                            prot_tmp]

                try:
                    returncode = call(hmm_cmd)
                except Exception as err:
                    raise RuntimeError("{0} failed : {1}".format(hmm_cmd[0], err))
                if returncode != 0:
                    raise RuntimeError("{0} failed return code = {1}".format(hmm_cmd[0], returncode))
                hmm_in = read_hmm(replicon.name, hmm_out, cfg, evalue=evalue, coverage=coverage
                                  ).sort_values("evalue").drop_duplicates(subset="ID_prot")
                func_annotate_res = pd.concat([func_annotate_res, hmm_in])
            func_annotate_res = func_annotate_res.sort_values("evalue").drop_duplicates(subset="ID_prot")

            integron.proteins.loc[func_annotate_res.ID_prot, "evalue"] = func_annotate_res.evalue.values
            integron.proteins.loc[func_annotate_res.ID_prot, "annotation"] = func_annotate_res.query_name.values
            integron.proteins.loc[func_annotate_res.ID_prot, "model"] = func_annotate_res.ID_query.values
            integron.proteins = integron.proteins.astype(dtype=integron.dtype)
