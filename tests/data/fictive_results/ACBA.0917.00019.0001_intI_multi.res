# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  tests/data/Models/phage-int.hmm
# target sequence database:        tests/data/Results_Integron_Finder_ACBA.0917.00019.gembase/tmp_ACBA.0917.00019.0001/ACBA.0917.00019.0001.prt
# output directed to file:         hmm_result/ACBA.0917.00019.0001.res
# per-seq hits tabular output:     hmm_result/ACBA.0917.00019.0001_intI_table.res
# sequence reporting threshold:    E-value <= 1e-50
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       Phage_integrase  [M=173]
Accession:   PF00589.16
Description: Phage integrase family
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence                    Description
    ------- ------ -----    ------- ------ -----   ---- --  --------                    -----------
    2.3e-66  220.8   0.0    5.5e-66  219.5   0.0    1.6  2  ACBA.0917.00019.i0001_00298  1035 xerD_1 | Tyrosine recombina
    1.7e-51  172.4   0.1    3.4e-51  171.4   0.1    1.5  2  ACBA.0917.00019.i0001_00338  921 xerD_2 | Tyrosine recombinas


Domain annotation for each sequence (and alignments):
>> ACBA.0917.00019.i0001_00298  1035 xerD_1 | Tyrosine recombinase XerD | NA | similar to AA sequence:UniProtKB:P9WF33
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?   -2.3   0.0       0.4   7.7e+02       3      26 ..      56      83 ..      54      88 .. 0.66
   2 !  219.5   0.0   2.9e-69   5.5e-66       2     171 ..     116     317 ..     115     319 .. 0.99

  Alignments for each domain:
  == domain 1  score: -2.3 bits;  conditional E-value: 0.4
                                 HHHHHHHHHHHCCCT....HHHHHHHHH CS
              Phage_integrase  3 Ltedeverllaalee....slsirdrll 26
                                 L ++eve++l +l +    s s+++++l
  ACBA.0917.00019.i0001_00298 56 LGSSEVEAFLSWLANerkvSVSTHRQAL 83
                                 6678999999999985555555555555 PP

  == domain 2  score: 219.5 bits;  conditional E-value: 2.9e-69
                                  HHHHHHHHHHHHCCCTHHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHH CS
              Phage_integrase   2 vLtedeverllaaleeslsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkei 83 
                                  vLt+dev+r+l++le+    ++rl+++ll++Tg+RisE l+lrvkdld+d+gti+v  re+K++k+r+++l+e+l+++l+e+
  ACBA.0917.00019.i0001_00298 116 VLTPDEVVRILGFLEG----EHRLFAQLLYGTGMRISEGLQLRVKDLDFDHGTIIV--REGKGSKDRALMLPESLAPSLREQ 191
                                  8***************....************************************..************************ PP

                                  HHH......HHTTSTTS......................BSSBEC...........TSSB..HHHHHHHHHHHHHHTT--CC CS
              Phage_integrase  84 lsdr.....kkeaeere......................llfvsk...........rgkplsdstvnrafkravkeagieke 127
                                  ls++     k++ae+r+                      ++f+++           r+++++d+t++rafkrav++agi+k+
  ACBA.0917.00019.i0001_00298 192 LSRArawwlKDQAEGRSgvalpdalerkypraghswpwfWVFAQHthstdprsgvvRRHHMYDQTFQRAFKRAVEQAGITKP 273
                                  *******9999999999***************************************************************** PP

                                  -HHHHHHHHHHHHHHH----HHHHHHH----SHHHHHHHHCCSH CS
              Phage_integrase 128 ltpHtLRhsfatallesGvdlkvvqkllGHssisttkiYthvak 171
                                  +tpHtLRhsfatall+sG+d+++vq+llGHs++stt+iYthv k
  ACBA.0917.00019.i0001_00298 274 ATPHTLRHSFATALLRSGYDIRTVQDLLGHSDVSTTMIYTHVLK 317
                                  *****************************************987 PP

>> ACBA.0917.00019.i0001_00338  921 xerD_2 | Tyrosine recombinase XerD | NA | similar to AA sequence:UniProtKB:P0A8P8
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?   -3.0   0.0      0.67   1.3e+03      59      83 ..      73      99 ..      59     107 .. 0.64
   2 !  171.4   0.1   1.7e-54   3.4e-51       3     172 ..     123     293 ..     121     294 .. 0.98

  Alignments for each domain:
  == domain 1  score: -3.0 bits;  conditional E-value: 0.67
                                 ..CCSSSCCEEEEE-..HHHHHHHHHH CS
              Phage_integrase 59 aretKtkkertvpls..eellevlkei 83
                                 ++ +K+ ++++  ls  ++++++l+e+
  ACBA.0917.00019.i0001_00338 73 TKVGKSPRSIARCLSalRQFYKFLREQ 99
                                 445666666777777777777777775 PP

  == domain 2  score: 171.4 bits;  conditional E-value: 1.7e-54
                                  HHHHHHHHHHHCCCT..HHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHH CS
              Phage_integrase   3 Ltedeverllaalee..slsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlke 82 
                                  L e++ve+l++a +    l +rdr+++e+l+++glR+sEll+lr++ ++l++g+ r+    +K++ker vpl++ +++++++
  ACBA.0917.00019.i0001_00338 123 LSEEDVEALIQAPDIttALGLRDRAMFEVLYACGLRVSELLNLRLELINLKQGYLRI---TGKGNKERLVPLGQYACDWVER 201
                                  899***********99999**************************************...********************** PP

                                  HHHH...HHTTSTTSBSSBECTSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHH CS
              Phage_integrase  83 ilsdr..kkeaeerellfvskrgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHssist 162
                                  +l+++  +  ++++++lf++++g  +s+++++ a+kr++ +a+i+ el+pHtLRh fat+ll++G+dl+vvq llGHs++st
  ACBA.0917.00019.i0001_00338 202 YLNEArpQLYKSSTDYLFLTQHGGIMSRQNFWYAIKRYALQANIQAELSPHTLRHAFATHLLNHGADLRVVQMLLGHSDLST 283
                                  *****9999999********************************************************************** PP

                                  HHHHHCCSHH CS
              Phage_integrase 163 tkiYthvake 172
                                  t+iYthva+ 
  ACBA.0917.00019.i0001_00338 284 TQIYTHVAQV 293
                                  *******975 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (173 nodes)
Target sequences:                         3870  (1172415 residues searched)
Passed MSV filter:                       129  (0.0333333); expected 77.4 (0.02)
Passed bias filter:                       95  (0.0245478); expected 77.4 (0.02)
Passed Vit filter:                        18  (0.00465116); expected 3.9 (0.001)
Passed Fwd filter:                        10  (0.00258398); expected 0.0 (1e-05)
Initial search space (Z):               3870  [actual number of targets]
Domain search space  (domZ):               2  [number of targets reported over threshold]
# CPU time: 0.04u 0.00s 00:00:00.04 Elapsed: 00:00:00.02
# Mc/sec: 10141.39
//
[ok]
