# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/data/Models/integron_integrase.hmm
# target sequence database:        /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_OBAL001.B.00005.C001/other_OBAL001.B.00005.C001/OBAL001.B.00005.C001.prt
# output directed to file:         /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_OBAL001.B.00005.C001/other_OBAL001.B.00005.C001/OBAL001.B.00005.C001_intI.res
# per-seq hits tabular output:     /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_OBAL001.B.00005.C001/other_OBAL001.B.00005.C001/OBAL001.B.00005.C001_intI_table.res
# number of worker threads:        1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       intI_Cterm  [M=59]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence                  Description
    ------- ------ -----    ------- ------ -----   ---- --  --------                  -----------
    6.7e-22   73.9   0.1    1.1e-21   73.2   0.1    1.4  1  OBAL001.B.00005.C001_1416  # 1545830 # 1546807 # -1 # ;gc_con


Domain annotation for each sequence (and alignments):
>> OBAL001.B.00005.C001_1416  # 1545830 # 1546807 # -1 # ;gc_cont=0.456
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   73.2   0.1   5.1e-25   1.1e-21       2      59 .]     187     243 ..     186     243 .. 0.96

  Alignments for each domain:
  == domain 1  score: 73.2 bits;  conditional E-value: 5.1e-25
                 intI_Cterm   2 lhekDlaegyggVyLPnaLarKYPnaakelaWqylFPsaklsvdprsgelrRHHldes 59 
                                ++++D+ +g+g  +LP aL+rKYP+a ++ aW+++FPs++l+++p++g+l+RHHl++s
  OBAL001.B.00005.C001_1416 187 IQQDDNLQGVGP-SLPFALDRKYPSAYRQAAWMFVFPSSTLCNHPYNGKLCRHHLHDS 243
                                689999******.******************************************985 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (59 nodes)
Target sequences:                         2218  (714377 residues searched)
Passed MSV filter:                        48  (0.0216411); expected 44.4 (0.02)
Passed bias filter:                       38  (0.0171326); expected 44.4 (0.02)
Passed Vit filter:                         3  (0.00135257); expected 2.2 (0.001)
Passed Fwd filter:                         1  (0.000450857); expected 0.0 (1e-05)
Initial search space (Z):               2218  [actual number of targets]
Domain search space  (domZ):               1  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# Mc/sec: inf
//
[ok]
