[![install](https://img.shields.io/github/downloads/gem-pasteur/Integron_finder/total.svg)](https://github.com/gem-pasteur/Integron_Finder/releases/latest) [![Doc](https://readthedocs.org/projects/integronfinder/badge/?version=latest)](https://integronfinder.readthedocs.io/en/latest) [![License (GPL version 3)](https://img.shields.io/badge/license-GNU%20GPL%20version%203-blue.svg?style=flat-square)](https://opensource.org/licenses/GPL-3.0) [![Build Status](https://travis-ci.org/gem-pasteur/Integron_Finder.svg?branch=master)](https://travis-ci.org/gem-pasteur/Integron_Finder)

# Integron Finder

Finds integrons in DNA sequences

You can use it in command line, see *installation* below,
or you can use it online on the
[Galaxy Pasteur](https://galaxy.pasteur.fr/root?tool_id=toolshed.pasteur.fr%2Frepos%2Fkhillion%2Fintegron_finder%2Fintegron_finder%2F1.5.1).

See Documentation for how to use it:
[![Doc](https://readthedocs.org/projects/integronfinder/badge/?version=latest)](https://integronfinder.readthedocs.io/en/latest)

## Installation

*Nota Bene*: version 2.0 isn't fully released yet. The pip installation below is for the version 1.5 with python2.7. 
If you want the version 2.0, you can install it with pip, but there are still a few minor bugs around (see issues). 

    pip install git+https://github.com/gem-pasteur/Integron_Finder/#egg=integron_finder

### For user (version 1.5)

    pip install integron_finder

for more installation options, or for developer installation see documentation

#### Singularity container

For reproducibility and easy way to use integron_finder without installing
third party software (hmmsearch, ...) or libraries, we provide containers based on singularity

[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/1314)

So you just have to install singularity (https://github.com/singularityware/singularity)

then

    singularity run shub:gem-pasteur/Integron_Finder -h

to run the last version (from master branch)
if you prefer download a specific version (with renaming the container on the fly)

    singularity pull --name integron_finder shub:gem-pasteur/Integron_Finder:2.0

when you have the image in local you can use it as if integron_finder has been installed

    ./integron_finder -h

#### Conda installation [not yet available]

From 2.0 version, Integron_Finder is available as [conda](https://conda.io/docs/index.html) package.
Integron_finder is in [bioconda](https://bioconda.github.io/) channel.
(The advantage with this solution is that it will install prodigal, hmmer, and infernal too.)

1. install conda
2. Set up channels:

        conda config --add channels defaults
        conda config --add channels conda-forge
        conda config --add channels bioconda

3. install integron_finder:

        conda install integron_finder

### For developer

If you want to develop or submit a patch on this software you are welcome.
See [Developer installation](https://integronfinder.readthedocs.io/en/latest/developer_guide/dev_guide.html#developer-installation) in documentation.


## Licence:

* Integron Finder is developed and released under [open source licence GPLv3](https://opensource.org/licenses/GPL-3.0)
  (see COPYING file)
* *Resfams* is provided with the courtesy of Dantas, Gautam http://www.dantaslab.org/resfams
* The other data files:
    * *attc_4.cm*
    * *integron_integrase.hmm*
    * *phage-int.hmm*     
  are licensed under [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/)
  [![CC BY-NC-SA](https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png)](http://creativecommons.org/licenses/by-nc-sa/4.0/)
   
## Contributing

We encourage contributions, bug report, enhancement ...

But before to do that we encourage to read [the contributing guide](CONTRIBUTING.md).

## Dependencies

- Python >=3.4
- Pandas >=0.22.0
- Numpy >=1.14.2
- Biopython >=1.70
- Matplotlib >=2.2.2
- HMMER >=3.1b2
- INFERNAL >=1.1.2
- Prodigal >=2.6.2

## Usage

```
usage: integron_finder [-h] [--local-max] [--func-annot] [--cpu CPU]
                       [-dt DISTANCE_THRESHOLD] [--outdir OUTDIR]
                       [--union-integrases] [--cmsearch CMSEARCH]
                       [--hmmsearch HMMSEARCH] [--prodigal PRODIGAL]
                       [--path-func-annot PATH_FUNC_ANNOT] [--gembase]
                       [--annot-parser ANNOT_PARSER_NAME]
                       [--attc-model ATTC_MODEL] [--evalue-attc EVALUE_ATTC]
                       [--calin-threshold CALIN_THRESHOLD]
                       [--keep-palindromes] [--no-proteins] [--promoter-attI]
                       [--max-attc-size MAX_ATTC_SIZE]
                       [--min-attc-size MIN_ATTC_SIZE] [--eagle-eyes] [--pdf]
                       [--gbk] [--keep-tmp] [--split-results]
                       [--circ | --linear] [--topology-file TOPOLOGY_FILE]
                       [-V] [--mute] [-v] [-q]
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
                        DISTANCE_THRESH [4000]bp or less
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
  --annot-parser ANNOT_PARSER_NAME
                        the name of the parser to use to get information from
                        protein file.
  --attc-model ATTC_MODEL
                        Path or file to the attc model (Covariance Matrix).
  --evalue-attc EVALUE_ATTC
                        Set evalue threshold to filter out hits above it
                        (default: 1)
  --calin-threshold CALIN_THRESHOLD
                        keep 'CALIN' only if attC sites nuber >= calin-
                        threshold (default: 2)
  --keep-palindromes    For a given hit, if the palindromic version is found,
                        don't remove the one with highest evalue.
  --no-proteins         Don't annotate CDS and don't find integrase, just look
                        for attC sites.
  --promoter-attI       Search also for promoter and attI sites. (default
                        False)
  --max-attc-size MAX_ATTC_SIZE
                        Set maximum value fot the attC size (default: 200bp)
  --min-attc-size MIN_ATTC_SIZE
                        set minimum value fot the attC size (default: 40bp)
  --eagle-eyes          Synonym of --local-max. Like a soaring eagle in the
                        sky, catching rabbits (or attC sites) by surprise.
  --circ                Set the default topology for replicons to 'circular'
  --linear              Set the default topology for replicons to 'linear'
  --topology-file TOPOLOGY_FILE
                        The path to a file where the topology for each
                        replicon is specified.
  -V, --version         show program's version number and exit
  --mute                mute the log on stdout.(continue to log on
                        integron_finder.out)

Output options:
  --pdf                 For each complete integron, a simple graphic of the
                        region is depicted (in pdf format)
  --gbk                 generate a GenBank file with the sequence annotated
                        with the same annotations than .integrons file.
  --keep-tmp            keep intermediate results. This results are stored in
                        directory named tmp_<replicon id>
  --split-results       Instead of merging integron results from all replicon
                        in one file, keep them in separated files.

  -v, --verbose         Increase verbosity of output (can be cumulative : -vv)
  -q, --quiet           Decrease verbosity of output (can be cumulative : -qq)
```


### Example

    integron_finder --local-max --func-annot mysequences.fst

### Output :

By default, integron_finder will output 3 files under Results_Integron_Finder_mysequences:

- `mysequences.integrons` : A file with all integrons and their elements detected in all sequences in the input file.
- `mysequences.summary` : A summary file with the number and type of integrons per sequence.
- `integron_finder.out` : A copy standard output. The stdout can be silenced with the argument --mute

The amount of log in the standard output can be controlled with `--verbose` for more or `--quiet` for less, and both are cumulative arguments, eg. `-vv` or `-qq`.

Other files can be created on demand:

- `--gbk`: Creates a Genbank files with all the annotations found (present in the .integrons file)
- `--pdf`: Creates a simple pdf graphic with complete integrons
- `--keep-tmp`: Keep temporary files. See Keep intermediate files for more.

# Galaxy

You can use this program without installing it, through the pasteur galaxy webserver instance:

* [Galaxy Pasteur](https://galaxy.pasteur.fr/root?tool_id=toolshed.pasteur.fr%2Frepos%2Fkhillion%2Fintegron_finder%2Fintegron_finder%2F1.5.1)

# Citation

The paper is published in Nucleic Acid Research.

**Identification and analysis of integrons and cassette arrays in bacterial genomes**
Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha
*Nucleic Acids Research 2016; [doi: 10.1093/nar/gkw319](https://academic.oup.com/nar/article/44/10/4539/2516972)*

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
