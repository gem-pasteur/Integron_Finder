[![install](https://img.shields.io/github/downloads/gem-pasteur/Integron_finder/total.svg)](https://github.com/gem-pasteur/Integron_Finder/releases/latest) [![Doc](https://readthedocs.org/projects/integronfinder/badge/?version=latest)](https://integronfinder.readthedocs.io/en/latest) [![License (GPL version 3)](https://img.shields.io/badge/license-GNU%20GPL%20version%203-blue.svg?style=flat-square)](https://opensource.org/licenses/GPL-3.0) [![Build Status](https://travis-ci.org/gem-pasteur/Integron_Finder.svg?branch=dev)](https://travis-ci.org/gem-pasteur/Integron_Finder)
  
Integron Finder
===============

Finds integrons in DNA sequences

You can use it in command line, see *installation* below, 
or you can use it online on the 
[Galaxy Pasteur](https://galaxy.pasteur.fr/root?tool_id=toolshed.pasteur.fr%2Frepos%2Fkhillion%2Fintegron_finder%2Fintegron_finder%2F1.5.1).

See Documentation for how to use it: 
[![Doc](https://readthedocs.org/projects/integronfinder/badge/?version=latest)](https://integronfinder.readthedocs.io/en/latest)

# Installation

for user

     pip install integron_finder

for more installation option, for developer see documentaiton

# Licence:

Integron Finder is under [open source licence GPLv3](https://opensource.org/licenses/GPL-3.0)

# Dependencies :

- Python 2.7
- Pandas 0.18.0
- Numpy 1.9.1
- Biopython 1.69
- Matplotlib 1.4.3
- psutils 2.1.3
- HMMER 3.1b1
- INFERNAL 1.1
- Prodigal V2.6.2

# Usage

```
usage: integron_finder [-h] [--local-max] [--func-annot] [--cpu CPU]
                       [-dt DISTANCE_THRESHOLD] [--outdir OUTDIR]
                       [--union-integrases] [--cmsearch CMSEARCH]
                       [--hmmsearch HMMSEARCH] [--prodigal PRODIGAL]
                       [--path-func-annot PATH_FUNC_ANNOT] [--gembase]
                       [--attc-model ATTC_MODEL] [--evalue-attc EVALUE_ATTC]
                       [--keep-palindromes] [--no-proteins]
                       [--max-attc-size MAX_ATTC_SIZE]
                       [--min-attc-size MIN_ATTC_SIZE] [--eagle-eyes] [--pdf]
                       [--gbk] [--keep-tmp] [--split-results]
                       [--circ | --linear] [--topology-file TOPOLOGY_FILE]
                       [-V] [-v] [-q]
                       replicon

positional arguments:
  replicon              Path to the replicon file (in fasta format), eg :
                        path/to/file.fst or file.fst

optional arguments:
  -h, --help            show this help message and exit
  --local-max           Allows thorough local detection (slower but more
                        sensitive and do not increase false positive rate).
  --func-annot          Functional annotation of CDS associated with integrons
                        HMM files are needed in Func_annot folder.
  --cpu CPU             Number of CPUs used by INFERNAL and HMMER
  -dt DISTANCE_THRESHOLD, --distance-thresh DISTANCE_THRESHOLD
                        Two elements are aggregated if they are distant of
                        DISTANCE_THRESH [4kb] or less
  --outdir OUTDIR       Set the output directory (default: current)
  --union-integrases    Instead of taking intersection of hits from Phage_int
                        profile (Tyr recombinases) and integron_integrase
                        profile, use the union of the hits
  --cmsearch CMSEARCH   Complete path to cmsearch if not in PATH. eg:
                        /usr/local/bin/cmsearch
  --hmmsearch HMMSEARCH
                        Complete path to hmmsearch if not in PATH. eg:
                        /usr/local/bin/hmmsearch
  --prodigal PRODIGAL   Complete path to prodigal if not in PATH. eg:
                        /usr/local/bin/prodigal
  --path-func-annot PATH_FUNC_ANNOT
                        Path to file containing all hmm bank paths (one per
                        line)
  --gembase             Use gembase formatted protein file instead of
                        Prodigal. Folder structure must be preserved
  --attc-model ATTC_MODEL
                        path or file to the attc model (Covariance Matrix)
  --evalue-attc EVALUE_ATTC
                        set evalue threshold to filter out hits above it
                        (default: 1)
  --keep-palindromes    for a given hit, if the palindromic version is found,
                        don't remove the one with highest evalue
  --no-proteins         Don't annotate CDS and don't find integrase, just look
                        for attC sites.
  --max-attc-size MAX_ATTC_SIZE
                        set maximum value fot the attC size (default: 200bp)
  --min-attc-size MIN_ATTC_SIZE
                        set minimum value fot the attC size (default: 40bp)
  --eagle-eyes          Synonym of --local-max. Like a soaring eagle in the
                        sky, catching rabbits (or attC sites) by surprise.
  --circ                Set the default topology for replicons to 'cirular'
  --linear              Set the default topology for replicons to 'linear'
  --topology-file TOPOLOGY_FILE
                        The path to a file where the topology for each
                        replicon is specified
  -V, --version         show program's version number and exit

Output options:
  --pdf                 For each complete integron, a simple graphic of the
                        region is depicted (in pdf format)
  --gbk                 generate a GenBank file with the sequence annotated
                        with the same annotations than .integrons file.
  --keep-tmp            keep intermediate results. This results are stored in
                        directory named other_<replicon id>
  --split-results       Instead of merging integron results from all replicon
                        in one file, keep them in separated files.

  -v, --verbose         Increase verbosity of output (can be cumulative : -vv)
  -q, --quiet           Decrease verbosity of output (can be cumulative : -qq)
```


### Example

    integron_finder myfastafile.fst --local_max --func_annot

## Output :

A folder name `Results_<id_genome>`, inside there are different files :

- ***.gbk** : contains the input sequence with all integrons and features found.
- ***.integrons** : contain list of all element detected (attc, protein near attC, integrase, Pc, attI, Pint) with position, 
  strand, evalue, etc...
- ***.pdf** : representation of complete integrons detected (with integrase (redish) and at least one attc (blueish)).
  If a protein has a hit with an antibiotic resistance gene, it's yellow, otherwise grey.

 and one folder, `other`, containing the different outputs of the different steps of the program.

# Galaxy

You can use this program without installing it, through the pasteur galaxy webserver instance:

* [Galaxy Pasteur](https://galaxy.pasteur.fr/root?tool_id=toolshed.pasteur.fr%2Frepos%2Fkhillion%2Fintegron_finder%2Fintegron_finder%2F1.5.1)

# Citation

The paper is published in Nucleic Acid Research.

**Identification and analysis of integrons and cassette arrays in bacterial genomes**  
Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha  
*Nucleic Acids Research 2016; [doi: 10.1093/nar/gkw319](http://nar.oxfordjournals.org/cgi/content/full/gkw319)*

 Please cite also the following articles:

 - Nawrocki, E.P. and Eddy, S.R. (2013)  
   **Infernal 1.1: 100-fold faster RNA homology searches.**  
   *Bioinformatics, 29, 2933-2935.*
 - Eddy, S.R. (2011)  
   **Accelerated Profile HMM Searches.**  
   *PLoS Comput Biol, 7, e1002195.*
 - Hyatt, D., Chen, G.L., Locascio, P.F., Land, M.L., Larimer, F.W. and Hauser, L.J. (2010)  
   **Prodigal: prokaryotic gene recognition and translation initiation site identification.**  
   *BMC Bioinformatics, 11, 119.*

and if you use the function `--func_annot` which uses Resfams:

 - Gibson, M.K., Forsberg, K.J. and Dantas, G. (2015)  
   **Improved annotation of antibiotic resistance determinants reveals microbial resistomes cluster by ecology.**  
   *ISME J, 9, 207-216.*
