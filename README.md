Integron_Finder
===============

Find integrons in genomes

# Dependencies :

- Python 2.7
   - Pandas 0.15.1
   - Numpy 1.9.1
   - Biopython 1.63
   - Matplotlib 1.4.3
- HMMER 3.1b1 
- INFERNAL 1.1 
- Prodigal V2.6.2 

# Usage 

```
usage: Integron_Finder.py [-h] [--max] [--resfams] [--cpu CPU]
                          [-dt DISTANCE_THRESH] [--outdir .] [--linear]
                          [--union_integrases] [--cmsearch CMSEARCH]
                          [--hmmsearch HMMSEARCH] [--prodigal PRODIGAL]
                          [--gembase] [--attc_model file.cm] [--evalue_attc 1]
                          [--keep_palindromes]
                          replicon

positional arguments:
  replicon              Path and/or fasta file to the replicon, eg :
                        path/to/file.fst or file.fst

optional arguments:
  -h, --help            show this help message and exit
  --max                 Allows exact local detection (slower)
  --resfams             Detect antibiotic resistances with Resfams HMM
                        profiles
  --cpu CPU             Number of CPUs used by INFERNAL and HMMER
  -dt DISTANCE_THRESH, --distance_thresh DISTANCE_THRESH
                        Two element are aggregated if they are distant of
                        DISTANCE_THRESH [4kb] or less
  --outdir .            set the output directory (default: current)
  --linear              consider replicon as linear. If replicon smaller than
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
  --gembase             Use gembase formatted protein file instead of
                        Prodigal. Folder structure must be preserved
  --attc_model file.cm  path or file to the attc model (Covariance Matrix)
  --evalue_attc 1       set evalue threshold to filter out hits above it
                        (default: 1)
  --keep_palindromes    for a given hit, if the palindromic version is found,
                        don't remove the one with highest evalue

```


### Example

python attc_finder.py myfastafile.fst --resfams --max

## Output :

A folder name Results\_id\_genome, inside their are many files, the one you care are :

- *.gbk : contains the input sequence with all integrons and features found.
- *.integrons : contain list of all element detected (attc, protein near attC, integrase, Pc, attI, Pint) with position, strand, evalue, etc... 
- *.pdf : representation of complete integrons detected (with integrase (redish) and at least one attc (blueish)). If a protein has a hit with an antibiotic resistance gene, it's yellow, otherwise grey.




