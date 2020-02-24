# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/data/Models/integron_integrase.hmm
# target sequence database:        /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_ACBA.0917.00019/tmp_ACBA.0917.00019.0001/ACBA.0917.00019.0001.prt
# output directed to file:         /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_ACBA.0917.00019/tmp_ACBA.0917.00019.0001/ACBA.0917.00019.0001_intI.res
# per-seq hits tabular output:     /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_ACBA.0917.00019/tmp_ACBA.0917.00019.0001/ACBA.0917.00019.0001_intI_table.res
# number of worker threads:        1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       intI_Cterm  [M=59]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence                    Description
    ------- ------ -----    ------- ------ -----   ---- --  --------                    -----------
      2e-23   79.6   3.3    3.6e-23   78.8   3.3    1.5  1  ACBA.0917.00019.i0001_00298  1035 xerD_1 | Tyrosine recombina


Domain annotation for each sequence (and alignments):
>> ACBA.0917.00019.i0001_00298  1035 xerD_1 | Tyrosine recombinase XerD | NA | similar to AA sequence:UniProtKB:P9WF33
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   78.8   3.3   9.4e-27   3.6e-23       2      58 ..     198     254 ..     197     255 .. 0.96

  Alignments for each domain:
  == domain 1  score: 78.8 bits;  conditional E-value: 9.4e-27
                   intI_Cterm   2 lhekDlaegyggVyLPnaLarKYPnaakelaWqylFPsaklsvdprsgelrRHHlde 58 
                                  ++ kD+aeg +gV LP+aL+rKYP+a++++ W+++F++++ s+dprsg++rRHH+ +
  ACBA.0917.00019.i0001_00298 198 WWLKDQAEGRSGVALPDALERKYPRAGHSWPWFWVFAQHTHSTDPRSGVVRRHHMYD 254
                                  6889***************************************************87 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (59 nodes)
Target sequences:                         3870  (1172415 residues searched)
Passed MSV filter:                        84  (0.0217054); expected 77.4 (0.02)
Passed bias filter:                       73  (0.018863); expected 77.4 (0.02)
Passed Vit filter:                         9  (0.00232558); expected 3.9 (0.001)
Passed Fwd filter:                         1  (0.000258398); expected 0.0 (1e-05)
Initial search space (Z):               3870  [actual number of targets]
Domain search space  (domZ):               1  [number of targets reported over threshold]
# CPU time: 0.01u 0.00s 00:00:00.01 Elapsed: 00:00:00.00
# Mc/sec: inf
//
[ok]
