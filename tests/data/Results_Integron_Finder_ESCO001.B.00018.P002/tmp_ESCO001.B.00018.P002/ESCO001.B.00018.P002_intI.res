# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/data/Models/integron_integrase.hmm
# target sequence database:        /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_ESCO001.B.00018.P002/other_ESCO001.B.00018.P002/ESCO001.B.00018.P002.prt
# output directed to file:         /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_ESCO001.B.00018.P002/other_ESCO001.B.00018.P002/ESCO001.B.00018.P002_intI.res
# per-seq hits tabular output:     /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_ESCO001.B.00018.P002/other_ESCO001.B.00018.P002/ESCO001.B.00018.P002_intI_table.res
# number of worker threads:        1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       intI_Cterm  [M=59]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence                 Description
    ------- ------ -----    ------- ------ -----   ---- --  --------                 -----------
    7.6e-25   79.6   3.3    1.4e-24   78.8   3.3    1.4  1  ESCO001.B.00018.P002_106  # 90229 # 91242 # -1 # ID=1_106;par


Domain annotation for each sequence (and alignments):
>> ESCO001.B.00018.P002_106  # 90229 # 91242 # -1 # ID=1_106;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   78.8   3.3   9.1e-27   1.4e-24       2      58 ..     198     254 ..     197     255 .. 0.96

  Alignments for each domain:
  == domain 1  score: 78.8 bits;  conditional E-value: 9.1e-27
                intI_Cterm   2 lhekDlaegyggVyLPnaLarKYPnaakelaWqylFPsaklsvdprsgelrRHHlde 58 
                               ++ kD+aeg +gV LP+aL+rKYP+a++++ W+++F++++ s+dprsg++rRHH+ +
  ESCO001.B.00018.P002_106 198 WWLKDQAEGRSGVALPDALERKYPRAGHSWPWFWVFAQHTHSTDPRSGVVRRHHMYD 254
                               6889***************************************************87 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (59 nodes)
Target sequences:                          151  (35057 residues searched)
Passed MSV filter:                         3  (0.0198675); expected 3.0 (0.02)
Passed bias filter:                        2  (0.013245); expected 3.0 (0.02)
Passed Vit filter:                         1  (0.00662252); expected 0.2 (0.001)
Passed Fwd filter:                         1  (0.00662252); expected 0.0 (1e-05)
Initial search space (Z):                151  [actual number of targets]
Domain search space  (domZ):               1  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.01
# Mc/sec: 206.84
//
[ok]
