# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/data/Models/phage-int.hmm
# target sequence database:        /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_acba.007.p01.13/tmp_ACBA.007.P01_13/ACBA.007.P01_13.prt
# output directed to file:         /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_acba.007.p01.13/tmp_ACBA.007.P01_13/ACBA.007.P01_13_phage_int.res
# per-seq hits tabular output:     /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_acba.007.p01.13/tmp_ACBA.007.P01_13/ACBA.007.P01_13_phage_int_table.res
# number of worker threads:        1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       Phage_integrase  [M=173]
Accession:   PF00589.16
Description: Phage integrase family
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence          Description
    ------- ------ -----    ------- ------ -----   ---- --  --------          -----------
    4.4e-39  124.7   0.0    1.2e-38  123.2   0.0    1.7  2  ACBA.007.P01_13_1  # 55 # 1014 # 1 # ;gc_cont=0.585


Domain annotation for each sequence (and alignments):
>> ACBA.007.P01_13_1  # 55 # 1014 # 1 # ;gc_cont=0.585
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?   -2.2   0.0      0.18       4.2       3      26 ..      56      83 ..      54      88 .. 0.66
   2 !  123.2   0.0   5.4e-40   1.2e-38       2     121 ..     116     267 ..     115     279 .. 0.97

  Alignments for each domain:
  == domain 1  score: -2.2 bits;  conditional E-value: 0.18
                       HHHHHHHHHHHCCCT....HHHHHHHHH CS
    Phage_integrase  3 Ltedeverllaalee....slsirdrll 26
                       L ++eve++l +l +    s s+++++l
  ACBA.007.P01_13_1 56 LGSSEVEAFLSWLANerkvSVSTHRQAL 83
                       6678999999999985555555555555 PP

  == domain 2  score: 123.2 bits;  conditional E-value: 5.4e-40
                        HHHHHHHHHHHHCCCTHHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHHHHH......H CS
    Phage_integrase   2 vLtedeverllaaleeslsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkeilsdr.....k 88 
                        vLt+dev+r+l++le+    ++rl+++ll++Tg+RisE l+lrvkdld+d+gti+v  re+K++k+r+++l+e+l+++l+e+ls++     k
  ACBA.007.P01_13_1 116 VLTPDEVVRILGFLEG----EHRLFAQLLYGTGMRISEGLQLRVKDLDFDHGTIIV--REGKGSKDRALMLPESLAPSLREQLSRArawwlK 201
                        8***************....************************************..*******************************999 PP

                        HTTSTTS......................BSSBEC...........TSSB..HHHHHHHHHHHHHH CS
    Phage_integrase  89 keaeere......................llfvsk...........rgkplsdstvnrafkravke 121
                        ++ae+r+                      ++f+++           r+++++d+t++rafkrav+ 
  ACBA.007.P01_13_1 202 DQAEGRSgvalpdalerkypraghswpwfWVFAQHthstdprsgvvRRHHMYDQTFQRAFKRAVEG 267
                        9999999*********************************************************97 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (173 nodes)
Target sequences:                           23  (5678 residues searched)
Passed MSV filter:                         1  (0.0434783); expected 0.5 (0.02)
Passed bias filter:                        1  (0.0434783); expected 0.5 (0.02)
Passed Vit filter:                         1  (0.0434783); expected 0.0 (0.001)
Passed Fwd filter:                         1  (0.0434783); expected 0.0 (1e-05)
Initial search space (Z):                 23  [actual number of targets]
Domain search space  (domZ):               1  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# Mc/sec: inf
//
[ok]
