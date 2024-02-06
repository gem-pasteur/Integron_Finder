[![testing](https://github.com/gem-pasteur/Integron_Finder/actions/workflows/testing.yml/badge.svg?branch=master)](https://github.com/gem-pasteur/Integron_Finder/actions/workflows/testing.yml)
[![codecov](https://codecov.io/gh/gem-pasteur/Integron_Finder/branch/master/graph/badge.svg?token=hZbOx1MEM7)](https://codecov.io/gh/gem-pasteur/Integron_Finder)
[![Doc](https://readthedocs.org/projects/integronfinder/badge/?version=latest)](https://integronfinder.readthedocs.io/en/latest) 
[![License (GPL version 3)](https://img.shields.io/badge/license-GNU%20GPL%20version%203-blue.svg?style=flat-square)](https://opensource.org/licenses/GPL-3.0)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/integron_finder)](https://pypi.org/project/integron_finder/)
[![PyPI](https://img.shields.io/pypi/v/integron_finder)](https://pypi.org/project/integron_finder/)
[![Downloads](https://pepy.tech/badge/integron-finder)](https://github.com/gem-pasteur/Integron_Finder/releases/latest) 
[![Docker Image Version (tag latest semver)](https://img.shields.io/docker/v/gempasteur/integron_finder?label=docker&sort=semver)](https://hub.docker.com/r/gempasteur/integron_finder/tags)
![Conda](https://img.shields.io/conda/v/bioconda/integron_finder)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/gem-pasteur/Integron_Finder/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/gem-pasteur/Integron_Finder)


# Integron Finder

Finds integrons in DNA sequences

You can use it in command line, see *installation* below,
or you can use it online on the
[Galaxy Pasteur](https://galaxy.pasteur.fr/root?tool_id=toolshed.pasteur.fr%2Frepos%2Fkhillion%2Fintegron_finder%2Fintegron_finder%2F2.0.1).

See Documentation for how to use it:
[![Doc](https://readthedocs.org/projects/integronfinder/badge/?version=latest)](https://integronfinder.readthedocs.io/en/latest)

## Installation

Although a system wide installation is possible and supported, many distribution do not allow it.
So we describe bellow some user wide installation procedures.

### For user

    pip install --user integron_finder==2.xx

for more installation options, or for developer installation see documentation

### In a virtualenv

To avoid interaction with the system libraries you can install integron_finder in a [virtualenv](https://docs.python.org/3/library/venv.html).

1. create and activate the virtualenv
   ```bash
    python -m venv Integron_Finder
    ./Integron_Finder/bin/activate
   ```
2. install integron_finder
    ```bash
    (Integron_Finder) python -m pip install integron_finder
    ```
   all libraries will be located in `Integron_Finder` directory
3. when you want to quit the virtualenv
    ```bash
    (Integron_Finder) deactivate
    ```

#### Container

For reproducibility and easy way to use integron_finder without installing
third party software (hmmsearch, prodigal, ...) or libraries, we provide containers based on docker.

https://hub.docker.com/r/gempasteur/integron_finder

##### Docker

The computation are perform under IF user in /home/IF inside the container. 
So You have to mount a directory from the host in the container to exchange data 
(inputs data, and results) from the host and the container. 

The shared directory must be writable by the IF user or overwrite the user in the container by your id (see example below)

```
mkdir shared_dir
cd shared_dir
docker run -v $PWD:/home/IF -u $(id -u ${USER}):$(id -g ${USER}) integron_finder:2.0rc9 --local-max --circ --keep-tmp NZ_CP016323.fna
```

##### Singularity

As the docker image is registered in docker hub you can also use it directly with *Singularity*. 
Unlike *docker*, you have not to worry about shared directory, your `home` and `/tmp` are automatically shared.

```
singularity run -H ${HOME} docker://gempasteur/integron_finder:2.0rc9  --local-max --circ --keep-tmp NZ_CP016323.fna
```

or use *-b* option if the data is not in your home.

```
singularity run -H ${HOME} -b <the directory containing data> docker://gempasteur/integron_finder:2.0rc9 --local-max --circ --keep-tmp NZ_CP016323.fna
```

#### Conda installation

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
See [Developer installation](https://integronfinder.readthedocs.io/en/latest/developer_guide/dev_guide.html#developer-installation) 
in documentation.


## Licence:

* Integron Finder is developed and released under [open source licence GPLv3](https://opensource.org/licenses/GPL-3.0)
  (see COPYING file)
* *NCBIfam-AMRFinder* is provided by NCBI and accessible here: https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/2021-06-01.1/NCBIfam-AMRFinder.LIB
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

- Python >=3.10
- Pandas >=2
- Numpy >=1.26
- Biopython >=1.82
- Matplotlib >=3.8
- colorlog
- HMMER >=3.1b2,<=3.3.2
- INFERNAL >=1.1.2,<=1.1.4
- Prodigal >=2.6.2,<=V2.6.3

## Usage

```
usage: integron_finder [-h] [--local-max] [--func-annot] [--cpu CPU] [-dt DISTANCE_THRESHOLD] [--outdir OUTDIR] [--union-integrases]
                       [--cmsearch CMSEARCH] [--hmmsearch HMMSEARCH] [--prodigal PRODIGAL] [--path-func-annot PATH_FUNC_ANNOT] [--gembase]
                       [--gembase-path GEMBASE_PATH] [--annot-parser ANNOT_PARSER] [--prot-file PROT_FILE] [--attc-model ATTC_MODEL]
                       [--evalue-attc EVALUE_ATTC] [--calin-threshold CALIN_THRESHOLD] [--keep-palindromes] [--no-proteins]
                       [--promoter-attI] [--max-attc-size MAX_ATTC_SIZE] [--min-attc-size MIN_ATTC_SIZE] [--eagle-eyes] [--pdf] [--gbk]
                       [--keep-tmp] [--split-results] [--circ | --linear] [--topology-file TOPOLOGY_FILE] [--version] [--mute] [-v] [-q]
                       replicon

positional arguments:
  replicon              Path to the replicon file (in fasta format), eg : path/to/file.fst or file.fst

optional arguments:
  -h, --help            show this help message and exit
  --local-max           Allows thorough local detection (slower but more sensitive and do not increase false positive rate).
  --func-annot          Functional annotation of CDS associated with integrons HMM files are needed in Func_annot folder.
  --cpu CPU             Number of CPUs used by INFERNAL and HMMER. Increasing too much (usually above 4) may decrease performance. 
  -dt DISTANCE_THRESHOLD, --distance-thresh DISTANCE_THRESHOLD
                        Two elements are aggregated if they are distant of DISTANCE_THRESH [4000]bp or less
  --outdir OUTDIR       Set the output directory (default: current)
  --union-integrases    Instead of taking intersection of hits from Phage_int profile (Tyr recombinases) and integron_integrase profile,
                        use the union of the hits
  --cmsearch CMSEARCH   Complete path to cmsearch if not in PATH. eg: /usr/local/bin/cmsearch
  --hmmsearch HMMSEARCH
                        Complete path to hmmsearch if not in PATH. eg: /usr/local/bin/hmmsearch
  --prodigal PRODIGAL   Complete path to prodigal if not in PATH. eg: /usr/local/bin/prodigal
  --path-func-annot PATH_FUNC_ANNOT
                        Path to file containing all hmm bank paths (one per line)
  --gembase             Use gembase formatted protein file instead of Prodigal. Folder structure must be preserved
  --gembase-path GEMBASE_PATH
                        path to the gembase root directory (needed only if the replicon file is not locatedin gembase-path)
  --annot-parser ANNOT_PARSER
                        the path to the parser to use to get information from protein file.
  --prot-file PROT_FILE
                        The path to the proteins file used for annotations
  --attc-model ATTC_MODEL
                        Path or file to the attc model (Covariance Matrix).
  --evalue-attc EVALUE_ATTC
                        Set evalue threshold to filter out hits above it (default: 1)
  --calin-threshold CALIN_THRESHOLD
                        keep 'CALIN' only if attC sites number >= calin-threshold (default: 2)
  --keep-palindromes    For a given hit, if the palindromic version is found, don't remove the one with highest evalue.
  --no-proteins         Don't annotate CDS and don't find integrase, just look for attC sites.
  --promoter-attI       Search also for promoter and attI sites. (default False)
  --max-attc-size MAX_ATTC_SIZE
                        Set maximum value fot the attC size (default: 200bp)
  --min-attc-size MIN_ATTC_SIZE
                        set minimum value fot the attC size (default: 40bp)
  --eagle-eyes          Synonym of --local-max. Like a soaring eagle in the sky, catching rabbits (or attC sites) by surprise.
  --circ                Set the default topology for replicons to 'circular'
  --linear              Set the default topology for replicons to 'linear'
  --topology-file TOPOLOGY_FILE
                        The path to a file where the topology for each replicon is specified.
  --version             show program's version number and exit
  --mute                mute the log on stdout.(continue to log on integron_finder.out)

Output options:
  --pdf                 For each complete integron, a simple graphic of the region is depicted (in pdf format)
  --gbk                 generate a GenBank file with the sequence annotated with the same annotations than .integrons file.
  --keep-tmp            keep intermediate results. This results are stored in directory named tmp_<replicon id>
  --split-results       Instead of merging integron results from all replicon in one file, keep them in separated files.

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

The amount of log in the standard output can be controlled with `--verbose` for more or `--quiet` for less, 
and both are cumulative arguments, eg. `-vv` or `-qq`.

Other files can be created on demand:

- `--gbk`: Creates a Genbank files with all the annotations found (present in the .integrons file)
- `--pdf`: Creates a simple pdf graphic with complete integrons
- `--keep-tmp`: Keep temporary files. See Keep intermediate files for more.

# Galaxy

You can use this program without installing it, through the pasteur galaxy webserver instance:

* [Galaxy Pasteur](https://galaxy.pasteur.fr/root?tool_id=toolshed.pasteur.fr%2Frepos%2Fkhillion%2Fintegron_finder%2Fintegron_finder%2F2.0.1)

# Citation

The paper is published in Microorganism.

Néron, Bertrand, Eloi Littner, Matthieu Haudiquet, Amandine Perrin, Jean Cury, and Eduardo P.C. Rocha. 2022. 
**IntegronFinder 2.0: Identification and Analysis of Integrons across Bacteria, with a Focus on Antibiotic Resistance in Klebsiella** 
Microorganisms 10, no. 4: 700. https://doi.org/10.3390/microorganisms10040700

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

and if you use the function `--func_annot` which uses *NCBIfam-AMRFinder* hmm profiles:

 - Haft, DH et al., Nucleic Acids Res. 2018 Jan 4;46(D1):D851-D860 PMID: 29112715 
 
