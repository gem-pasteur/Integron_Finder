# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /Users/aperrin/Softwares/gemsrc/Integron_Finder/data/Models/integron_integrase.hmm
# target sequence database:        /Users/aperrin/Softwares/gemsrc/Integron_Finder/datatest/Replicons/../Proteins/acba.007.p01.13.prt
# output directed to file:         acba-gembase/Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_intI.res
# per-seq hits tabular output:     acba-gembase/Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_intI_table.res
# number of worker threads:        1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       intI_Cterm  [M=59]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence           Description
    ------- ------ -----    ------- ------ -----   ---- --  --------           -----------
    1.3e-25   79.7   3.3    2.2e-25   78.9   3.3    1.4  1  ACBA007p01a_000009  D GTG TAA 55 1014


Domain annotation for each sequence (and alignments):
>> ACBA007p01a_000009  D GTG TAA 55 1014
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   78.9   3.3   8.4e-27   2.2e-25       2      58 ..     198     254 ..     197     255 .. 0.96

  Alignments for each domain:
  == domain 1  score: 78.9 bits;  conditional E-value: 8.4e-27
          intI_Cterm   2 lhekDlaegyggVyLPnaLarKYPnaakelaWqylFPsaklsvdprsgelrRHHlde 58 
                         ++ kD+aeg +gV LP+aL+rKYP+a++++ W+++F++++ s+dprsg++rRHH+ +
  ACBA007p01a_000009 198 WWLKDQAEGRSGVALPDALERKYPRAGHSWPWFWVFAQHTHSTDPRSGVVRRHHMYD 254
                         6889***************************************************87 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (59 nodes)
Target sequences:                           26  (5968 residues searched)
Passed MSV filter:                         2  (0.0769231); expected 0.5 (0.02)
Passed bias filter:                        1  (0.0384615); expected 0.5 (0.02)
Passed Vit filter:                         1  (0.0384615); expected 0.0 (0.001)
Passed Fwd filter:                         1  (0.0384615); expected 0.0 (1e-05)
Initial search space (Z):                 26  [actual number of targets]
Domain search space  (domZ):               1  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# Mc/sec: inf
//
[ok]
