def func_annot(replicon_name, out_dir, hmm_files, evalue=10, coverage=0.5):
    """
    Call hmmmer to annotate CDS associated with the integron. Use Resfams per default (Gibson et al, ISME J.,  2014)
    """
    print "# Start Functional annotation... : "
    prot_tmp = os.path.join(out_dir, replicon_name + "_subseqprot.tmp")

    for integron in integrons:
        if os.path.isfile(prot_tmp):
            os.remove(prot_tmp)

        if integron.type() != "In0" and len(integron.proteins) > 0:

            func_annotate_res = pd.DataFrame(columns=["Accession_number",
                                                      "query_name", "ID_query",
                                                      "ID_prot", "strand",
                                                      "pos_beg", "pos_end", "evalue"])


            prot_to_annotate = []
            prot = SeqIO.parse(PROT_file, "fasta")
            n_prot = 0
            for p in prot:
                n_prot += 1
                if p.id in integron.proteins.index:
                    prot_to_annotate.append(p)
            SeqIO.write(prot_to_annotate, prot_tmp, "fasta")
            for hmm in hmm_files:
                hmm_out = os.path.join(out_dir, "_".join([replicon_name,
                                                         hmm.split("/")[-1].split(".")[0],
                                                         "fa.res"]))
                hmm_tableout = os.path.join(out_dir, "_".join([replicon_name,
                                                         hmm.split("/")[-1].split(".")[0],
                                                         "fa_table.res"]))
                hmm_cmd = [HMMSEARCH,
                            "-Z", str(n_prot),
                            "--cpu", N_CPU,
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
                hmm_in = read_hmm(replicon_name, hmm_out, evalue=evalue, coverage=coverage).sort_values("evalue").drop_duplicates(subset="ID_prot")
                func_annotate_res = pd.concat([func_annotate_res, hmm_in])
            func_annotate_res = func_annotate_res.sort_values("evalue").drop_duplicates(subset="ID_prot")

            integron.proteins.loc[func_annotate_res.ID_prot, "evalue"] = func_annotate_res.evalue.values
            integron.proteins.loc[func_annotate_res.ID_prot, "annotation"] = func_annotate_res.query_name.values
            integron.proteins.loc[func_annotate_res.ID_prot, "model"] = func_annotate_res.ID_query.values
            integron.proteins = integron.proteins.astype(dtype=integron.dtype)
