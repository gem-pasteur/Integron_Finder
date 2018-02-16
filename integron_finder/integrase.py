def find_integrase(replicon_path, replicon_name, out_dir):
    """
    Call Prodigal for Gene annotation and hmmer to find integrase, either with phage_int
    HMM profile or with intI profile.

    :param replicon_path: the path of the replicon to analyse
    :type replicon_path: string
    :param replicon_name: the name of the replicon to analyse
    :type replicon_name: string
    :param out_dir: the relative path to the directory where prodigal outputs will be stored
    :type out_dir: str
    :returns: None, the results are written on the disk
    """
    if not args.gembase:
        # Test whether the protein file exist to avoid new annotation for each run on the same replicon
        prot_tr_path = os.path.join(out_dir, replicon_name + ".prt")
        if not os.path.isfile(prot_tr_path):
            dev_null = os.devnull
            if SIZE_REPLICON > 200000:
                prodigal_cmd = [PRODIGAL,
                                "-i", replicon_path,
                                "-a", prot_tr_path,
                                "-o", dev_null]

            else: # if small genome, prodigal annotate it as contig.
                prodigal_cmd = [PRODIGAL,
                                "-p", "meta",
                                "-i", replicon_path,
                                "-a", prot_tr_path,
                                "-o", dev_null]
            try:
                returncode = call(prodigal_cmd)
            except Exception as err:
                raise RuntimeError("{0} failed : {1}".format(prodigal_cmd[0], err))
            if returncode != 0:
                raise RuntimeError("{0} failed returncode = {1}".format(prodigal_cmd[0], returncode))

    intI_hmm_out = os.path.join(out_dir, replicon_name + "_intI.res")
    hmm_cmd = []
    if not os.path.isfile(intI_hmm_out):
        hmm_cmd.append([HMMSEARCH,
                        "--cpu", N_CPU,
                        "--tblout", os.path.join(out_dir, replicon_name + "_intI_table.res"),
                        "-o", intI_hmm_out,
                        MODEL_integrase,
                        PROT_file])

    phage_hmm_out = os.path.join(out_dir, replicon_name + "_phage_int.res")
    if not os.path.isfile(phage_hmm_out):
        hmm_cmd.append([HMMSEARCH,
                        "--cpu", N_CPU,
                        "--tblout", os.path.join(out_dir, replicon_name + "_phage_int_table.res"),
                        "-o", phage_hmm_out,
                        MODEL_phage_int,
                        PROT_file])

    for cmd in hmm_cmd:
        try:
            returncode = call(cmd)
        except Exception as err:
            raise RuntimeError("{0} failed : {1}".format(' '.join(cmd), err))
        if returncode != 0:
            raise RuntimeError("{0} failed return code = {1}".format(' '.join(cmd), returncode))

