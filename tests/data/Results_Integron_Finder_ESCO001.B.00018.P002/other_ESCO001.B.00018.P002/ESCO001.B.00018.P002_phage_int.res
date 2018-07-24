# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/data/Models/phage-int.hmm
# target sequence database:        /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_ESCO001.B.00018.P002/other_ESCO001.B.00018.P002/ESCO001.B.00018.P002.prt
# output directed to file:         /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_ESCO001.B.00018.P002/other_ESCO001.B.00018.P002/ESCO001.B.00018.P002_phage_int.res
# per-seq hits tabular output:     /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_ESCO001.B.00018.P002/other_ESCO001.B.00018.P002/ESCO001.B.00018.P002_phage_int_table.res
# number of worker threads:        1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       Phage_integrase  [M=173]
Accession:   PF00589.16
Description: Phage integrase family
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence                 Description
    ------- ------ -----    ------- ------ -----   ---- --  --------                 -----------
    8.1e-68  220.9   0.0    2.1e-67  219.6   0.0    1.7  2  ESCO001.B.00018.P002_106  # 90229 # 91242 # -1 # ID=1_106;par
    1.3e-37  122.6   0.0    2.1e-37  121.9   0.0    1.3  1  ESCO001.B.00018.P002_99   # 84261 # 85052 # -1 # ID=1_99;part
    7.8e-29   94.0   0.1    1.1e-28   93.6   0.1    1.2  1  ESCO001.B.00018.P002_82   # 70340 # 71080 # -1 # ID=1_82;part


Domain annotation for each sequence (and alignments):
>> ESCO001.B.00018.P002_106  # 90229 # 91242 # -1 # ID=1_106;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?   -2.3   0.0      0.59        30       3      26 ..      56      83 ..      54      88 .. 0.66
   2 !  219.6   0.0   4.1e-69   2.1e-67       2     171 ..     116     317 ..     115     319 .. 0.99

  Alignments for each domain:
  == domain 1  score: -2.3 bits;  conditional E-value: 0.59
                              HHHHHHHHHHHCCCT....HHHHHHHHH CS
           Phage_integrase  3 Ltedeverllaalee....slsirdrll 26
                              L ++eve++l +l +    s s+++++l
  ESCO001.B.00018.P002_106 56 LGSSEVEAFLSWLANerkvSVSTHRQAL 83
                              6678999999999985555555555555 PP

  == domain 2  score: 219.6 bits;  conditional E-value: 4.1e-69
                               HHHHHHHHHHHHCCCTHHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHHHHH CS
           Phage_integrase   2 vLtedeverllaaleeslsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkeilsd 86 
                               vLt+dev+r+l++le+    ++rl+++ll++Tg+RisE l+lrvkdld+d+gti+v  re+K++k+r+++l+e+l+++l+e+ls+
  ESCO001.B.00018.P002_106 116 VLTPDEVVRILGFLEG----EHRLFAQLLYGTGMRISEGLQLRVKDLDFDHGTIIV--REGKGSKDRALMLPESLAPSLREQLSR 194
                               8***************....************************************..*************************** PP

                               ......HHTTSTTS......................BSSBEC...........TSSB..HHHHHHHHHHHHHHTT--CC-HHHHH CS
           Phage_integrase  87 r.....kkeaeere......................llfvsk...........rgkplsdstvnrafkravkeagiekeltpHtL 133
                               +     k++ae+r+                      ++f+++           r+++++d+t++rafkrav++agi+k++tpHtL
  ESCO001.B.00018.P002_106 195 ArawwlKDQAEGRSgvalpdalerkypraghswpwfWVFAQHthstdprsgvvRRHHMYDQTFQRAFKRAVEQAGITKPATPHTL 279
                               ****9999999999*********************************************************************** PP

                               HHHHHHHHHH----HHHHHHH----SHHHHHHHHCCSH CS
           Phage_integrase 134 RhsfatallesGvdlkvvqkllGHssisttkiYthvak 171
                               Rhsfatall+sG+d+++vq+llGHs++stt+iYthv k
  ESCO001.B.00018.P002_106 280 RHSFATALLRSGYDIRTVQDLLGHSDVSTTMIYTHVLK 317
                               ***********************************987 PP

>> ESCO001.B.00018.P002_99  # 84261 # 85052 # -1 # ID=1_99;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  121.9   0.0   4.2e-39   2.1e-37       1     171 [.      44     236 ..      44     238 .. 0.86

  Alignments for each domain:
  == domain 1  score: 121.9 bits;  conditional E-value: 4.2e-39
                              -HHHHHHHHHHHHCCCTHHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..................CCSSSCC... CS
          Phage_integrase   1 kvLtedeverllaaleeslsirdrllvellleTglRisEllslrvkdldldngtirvp................aretKtkke... 67 
                              k+L + e+++ll+++ +   +++++l+ +l++Tg+Ri+E+l+l+++d++l  ++ +v+                            
  ESCO001.B.00018.P002_99  44 KYLLAPEISALLHYVPD---LHRKMLLATLWNTGARINEALALTRGDFSLAPPYPFVQlatlkqrtekatrtagR------VPags 120
                              68999************...**************************************99988888775554321......22444 PP

                              ...EEEEE-.HHHHHHHHHHHHH....HHTTSTTSBSSBECTSSB..HHHHHHHHHHHHHHTT--.....CC-HHHHHHHHHHHHH CS
          Phage_integrase  68 ...rtvpls.eellevlkeilsdr...kkeaeerellfvskrgkplsdstvnrafkravkeagie.....keltpHtLRhsfatal 141
                                 r vpls  +++++l+ +++      +  ++r+  + ++r + ++d+tv++++++av++a+ +      ++tpHt+Rhs+a+++
  ESCO001.B.00018.P002_99 121 qvhRLVPLSdTQYVSQLQMMVATLkipLERRNKRTGRMEKARIWGITDRTVRTWLNEAVEAAAADgvtfsVPVTPHTFRHSYAMHM 206
                              567999999888889998888888665444444444444559*******************9999999****************** PP

                              HH----HHHHHHH----SHHHHHHHHCCSH CS
          Phage_integrase 142 lesGvdlkvvqkllGHssisttkiYthvak 171
                              l +G++lkv+q+l+GH+sis+t++Yt+v++
  ESCO001.B.00018.P002_99 207 LYAGIPLKVLQSLMGHKSISSTEVYTKVFA 236
                              ***************************986 PP

>> ESCO001.B.00018.P002_82  # 70340 # 71080 # -1 # ID=1_82;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_c
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   93.6   0.1   2.1e-30   1.1e-28       2     171 ..      44     217 ..      43     219 .. 0.88

  Alignments for each domain:
  == domain 1  score: 93.6 bits;  conditional E-value: 2.1e-30
                              HHHHHHHHHHHHCCCTHHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE.........CCSSSCC..EEEEE-.HHHH CS
          Phage_integrase   2 vLtedeverllaaleeslsirdrllvellleTglRisEllslrvkdldldngtirvp.......aretKtkke..rtvpls.eell 77 
                              +L + ev++ll ++ +    r+++l+ +l++Tg+Ri+E++ l+ + +dl+  + +v+       ar++   k+  r vpl+  +++
  ESCO001.B.00018.P002_82  44 YLLAPEVSALLFYMPD---QRHHMLFATLWNTGMRIGEARMLTPESFDLNGVRPFVRilsekvrARRGRPPKDevRLVPLTdISYV 126
                              78899***********...9*************************************88888776777766667788888877888 PP

                              HHHHHHHHH.HHTTSTTSBSSBECTSSB..HHHHHHHHHHHHHHTT--.....CC-HHHHHHHHHHHHHHH----HHHHHHH---- CS
          Phage_integrase  78 evlkeilsdrkkeaeerellfvskrgkplsdstvnrafkravkeagie.....keltpHtLRhsfatallesGvdlkvvqkllGHs 158
                               +++ ++        +re      + + ++d+t+++++k+av++a+ +      ++tpHt+Rhs+++++l +  + kv+q l+GH+
  ESCO001.B.00018.P002_82 127 RQMESWMITT--RPRRRE------PLWAVTDETMRNWLKQAVRRAEADgvhfsIPVTPHTFRHSYIMHMLYHRQPRKVIQALVGHR 204
                              8888888777..333333......34456799************9999999999******************************** PP

                              SHHHHHHHHCCSH CS
          Phage_integrase 159 sisttkiYthvak 171
                              + +++++Yt+v++
  ESCO001.B.00018.P002_82 205 DPRSMEVYTRVFA 217
                              *********9986 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (173 nodes)
Target sequences:                          151  (35057 residues searched)
Passed MSV filter:                         8  (0.0529801); expected 3.0 (0.02)
Passed bias filter:                        7  (0.0463576); expected 3.0 (0.02)
Passed Vit filter:                         4  (0.0264901); expected 0.2 (0.001)
Passed Fwd filter:                         3  (0.0198675); expected 0.0 (1e-05)
Initial search space (Z):                151  [actual number of targets]
Domain search space  (domZ):               3  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# Mc/sec: inf
//
[ok]
