Integron_Finder
===============

Find integrons in genomes

Dependencies :

- Python 2.7
   - Pandas 0.15.1
   - Numpy 1.9.1
   - Biopython 1.63
   - Matplotlib 1.4.3
- HMMER 3.1b1 
- INFERNAL 1.1 
- Prodigal V2.6.2 

## Usage 

```
usage: attc_finder.py [-h] [--max] [--resfams] [--gembase]
                      [--union_integrases] [--attc_model file.cm]
                      [--evalue_attc 1] [--outdir .] [--keep_palindromes]
                      [--linear]
                      replicon

positional arguments:
  replicon              Path and/or fasta file to the replicon, eg :
                        path/to/file.fst or file.fst

optional arguments:
  -h, --help            show this help message and exit
  --max                 Allows exact local detection (slower)
  --resfams             Detect antibiotic resistances with Resfams HMM
                        profiles
  --gembase             Use gembase formatted protein file. Folder structure
                        must be preserved
  --union_integrases    Instead of taking intersection of hits from Phage_int
                        profile (Tyr recombinases) and integron_integrase
                        profile, use the union of the hits
  --attc_model file.cm  path or file to the attc model (Covariance Matrix)
  --evalue_attc 1       set evalue threshold to filter out hits above it
                        (default: 1)
  --outdir .            set the output directory (default: current)
  --keep_palindromes    for a given hit, if the palindromic version is found,
                        don't remove the one with highest evalue
  --linear              consider replicon as linear. If replicon smaller than
                        20kb, it will be considered as linear

```


### Example

python attc_finder.py myfastafile.fst --resfams --max

## Output :

A folder name Results\_id\_genome, inside their are many files, the one you care are :

- *.gbk : contains the input sequence with all integrons and features found.
- *.integrons : contain list of all element detected (attc, protein near attC, integrase, Pc, attI, Pint) with position, strand, evalue, etc... 
- *.pdf : representation of complete integrons detected (with integrase (redish) and at least one attc (blueish)). If a protein has a hit with an antibiotic resistance gene, it's yellow, otherwise grey.




