# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /Users/aperrin/Softwares/gemsrc/Integron_Finder/data/Models/integron_integrase.hmm
# target sequence database:        /Users/aperrin/Softwares/gemsrc/Integron_Finder/datatest/Replicons/../Proteins/acba.007.p01.13.prt
# output directed to file:         acba-gembase2/Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_intI.res
# per-seq hits tabular output:     acba-gembase2/Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_intI_table.res
# number of worker threads:        1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       intI_Cterm  [M=59]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence           Description
    ------- ------ -----    ------- ------ -----   ---- --  --------           -----------
    1.3e-42  134.1  10.1    4.5e-25   78.0   3.3    2.4  2  ACBA007p01a_000009  D GTG TAA 55 1014
    1.2e-25   79.8   3.3    2.3e-25   78.9   3.3    1.5  1  ACBA007p01a_000008  C GTG TAA 1 50


Domain annotation for each sequence (and alignments):
>> ACBA007p01a_000009  D GTG TAA 55 1014
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   60.1   0.6   1.2e-20   1.6e-19       2      49 ..     457     504 ..     456     505 .. 0.96
   2 !   78.0   3.3   3.3e-26   4.5e-25       2      58 ..     198     254 ..     197     255 .. 0.96

  Alignments for each domain:
  == domain 1  score: 60.1 bits;  conditional E-value: 1.2e-20
          intI_Cterm   2 lhekDlaegyggVyLPnaLarKYPnaakelaWqylFPsaklsvdprsg 49
                         ++ kD+aeg +gV LP+aL+rKYP+a++++ W+++F++++ s+dprsg
  ACBA007p01a_000009 457 WWLKDQAEGRSGVALPDALERKYPRAGHSWPWFWVFAQHTHSTDPRSG 504
                         6889******************************************98 PP

  == domain 2  score: 78.0 bits;  conditional E-value: 3.3e-26
          intI_Cterm   2 lhekDlaegyggVyLPnaLarKYPnaakelaWqylFPsaklsvdprsgelrRHHlde 58
                         ++ kD+aeg +gV LP+aL+rKYP+a++++ W+++F++++ s+dprsg++rRHH+ +
  ACBA007p01a_000009 198 WWLKDQAEGRSGVALPDALERKYPRAGHSWPWFWVFAQHTHSTDPRSGVVRRHHMYD 254
                         6889***************************************************87 PP

>> ACBA007p01a_000008  C GTG TAA 1 50
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   78.9   3.3   1.7e-26   2.3e-25       2      58 ..     198     254 ..     197     255 .. 0.96

  Alignments for each domain:
  == domain 1  score: 78.9 bits;  conditional E-value: 1.7e-26
          intI_Cterm   2 lhekDlaegyggVyLPnaLarKYPnaakelaWqylFPsaklsvdprsgelrRHHlde 58
                         ++ kD+aeg +gV LP+aL+rKYP+a++++ W+++F++++ s+dprsg++rRHH+ +
  ACBA007p01a_000008 198 WWLKDQAEGRSGVALPDALERKYPRAGHSWPWFWVFAQHTHSTDPRSGVVRRHHMYD 254
                         6889***************************************************87 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (59 nodes)
Target sequences:                           27  (6474 residues searched)
Passed MSV filter:                         3  (0.111111); expected 0.5 (0.02)
Passed bias filter:                        2  (0.0740741); expected 0.5 (0.02)
Passed Vit filter:                         2  (0.0740741); expected 0.0 (0.001)
Passed Fwd filter:                         2  (0.0740741); expected 0.0 (1e-05)
Initial search space (Z):                 27  [actual number of targets]
Domain search space  (domZ):               2  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# Mc/sec: inf
//
[ok]
