[![install](https://img.shields.io/github/downloads/gem-pasteur/Integron_finder/total.svg)](https://github.com/gem-pasteur/Integron_Finder/releases/latest) [![Doc] (https://readthedocs.org/projects/integronfinder/badge/?version=latest)](https://integronfinder.readthedocs.io/en/latest) [![License (GPL version 3)](https://img.shields.io/badge/license-GNU%20GPL%20version%203-blue.svg?style=flat-square)](https://opensource.org/licenses/GPL-3.0)

Integron Finder
===============

Find integrons in DNA sequences

See the [latest release](https://github.com/gem-pasteur/Integron_Finder/releases/latest) to download the program. Or, use it online on the [Mobyle Webserver](http://mobyle.pasteur.fr/cgi-bin/portal.py#forms::integron_finder)

See Documentation for how to install it and how to use it: [![Doc] (https://readthedocs.org/projects/integronfinder/badge/?version=latest)](https://integronfinder.readthedocs.io/en/latest)

# Licence:

Integron Finder is under [open source licence GPLv3](https://opensource.org/licenses/GPL-3.0)

# Dependencies :

- Python 2.7
   - Pandas 0.18.0
   - Numpy 1.9.1
   - Biopython 1.65
   - Matplotlib 1.4.3
   - psutils 2.1.3
- HMMER 3.1b1
- INFERNAL 1.1
- Prodigal V2.6.2

# Usage

```
usage: integron_finder [-h] [--local_max] [--func_annot] [--cpu CPU]
                       [-dt DISTANCE_THRESH] [--outdir .] [--linear]
                       [--union_integrases] [--cmsearch CMSEARCH]
                       [--hmmsearch HMMSEARCH] [--prodigal PRODIGAL]
                       [--path_func_annot bank_hmm] [--gembase]
                       [--attc_model file.cm] [--evalue_attc 1]
                       [--keep_palindromes] [--no_proteins]
                       [--max_attc_size 200] [--min_attc_size 40]
                       [--eagle_eyes] [-V]
                       replicon

positional arguments:
  replicon              Path to the replicon file (in fasta format), eg :
                        path/to/file.fst or file.fst

optional arguments:
  -h, --help            show this help message and exit
  --local_max           Allows thorough local detection (slower but more
                        sensitive and do not increase false positive rate).
  --func_annot          Functional annotation of CDS associated with integrons
                        HMM files are needed in Func_annot folder.
  --cpu CPU             Number of CPUs used by INFERNAL and HMMER
  -dt DISTANCE_THRESH, --distance_thresh DISTANCE_THRESH
                        Two elements are aggregated if they are distant of
                        DISTANCE_THRESH [4kb] or less
  --outdir .            Set the output directory (default: current)
  --linear              Consider replicon as linear. If replicon smaller than
                        20kb, it will be considered as linear
  --union_integrases    Instead of taking intersection of hits from Phage_int
                        profile (Tyr recombinases) and integron_integrase
                        profile, use the union of the hits
  --cmsearch CMSEARCH   Complete path to cmsearch if not in PATH. eg:
                        /usr/local/bin/cmsearch
  --hmmsearch HMMSEARCH
                        Complete path to hmmsearch if not in PATH. eg:
                        /usr/local/bin/hmmsearch
  --prodigal PRODIGAL   Complete path to prodigal if not in PATH. eg:
                        /usr/local/bin/prodigal
  --path_func_annot bank_hmm
                        Path to file containing all hmm bank paths (one per
                        line)
  --gembase             Use gembase formatted protein file instead of
                        Prodigal. Folder structure must be preserved
  --attc_model file.cm  path or file to the attc model (Covariance Matrix)
  --evalue_attc 1       set evalue threshold to filter out hits above it
                        (default: 1)
  --keep_palindromes    for a given hit, if the palindromic version is found,
                        don't remove the one with highest evalue
  --no_proteins         Don't annotate CDS and don't find integrase, just look
                        for attC sites.
  --max_attc_size 200   set maximum value fot the attC size (default: 200bp)
  --min_attc_size 40    set minimum value fot the attC size (default: 40bp)
  --eagle_eyes          Synonym of --local_max. Like a soaring eagle in the
                        sky, catching rabbits(or attC sites) by surprise.
  -V, --version         show program's version number and exit

```


### Example

    integron_finder myfastafile.fst --local_max --func_annot

## Output :

A folder name Results\_id\_genome, inside there are different files :

- *.gbk : contains the input sequence with all integrons and features found.
- *.integrons : contain list of all element detected (attc, protein near attC, integrase, Pc, attI, Pint) with position, strand, evalue, etc...
- *.pdf : representation of complete integrons detected (with integrase (redish) and at least one attc (blueish)). If a protein has a hit with an antibiotic resistance gene, it's yellow, otherwise grey.

 and one folder, `other`, containing the different outputs of the different steps of the program.

 # Mobyle

You can use this program whithout installing it, through a webserver:

http://mobyle.pasteur.fr/cgi-bin/portal.py#forms::integron_finder

# Citation

 The paper is published in Nucleic Acid Research.

Identification and analysis of integrons and cassette arrays in bacterial genomes
Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha
Nucleic Acids Research 2016; [doi: 10.1093/nar/gkw319](http://nar.oxfordjournals.org/cgi/content/full/gkw319)

 Please cite also the following articles:

 - Nawrocki, E.P. and Eddy, S.R. (2013) Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics, 29, 2933-2935.
 - Eddy, S.R. (2011) Accelerated Profile HMM Searches. PLoS Comput Biol, 7, e1002195.
 - Hyatt, D., Chen, G.L., Locascio, P.F., Land, M.L., Larimer, F.W. and Hauser, L.J. (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics, 11, 119.

 and if you use the function `--func_annot` which uses Resfams:

 - Gibson, M.K., Forsberg, K.J. and Dantas, G. (2015) Improved annotation of antibiotic resistance determinants reveals microbial resistomes cluster by ecology. ISME J, 9, 207-216.
