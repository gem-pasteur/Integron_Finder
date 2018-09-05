# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/data/Models/phage-int.hmm
# target sequence database:        /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_OBAL001.B.00005.C001/other_OBAL001.B.00005.C001/OBAL001.B.00005.C001.prt
# output directed to file:         /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_OBAL001.B.00005.C001/other_OBAL001.B.00005.C001/OBAL001.B.00005.C001_phage_int.res
# per-seq hits tabular output:     /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_OBAL001.B.00005.C001/other_OBAL001.B.00005.C001/OBAL001.B.00005.C001_phage_int_table.res
# number of worker threads:        1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       Phage_integrase  [M=173]
Accession:   PF00589.16
Description: Phage integrase family
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence                  Description
    ------- ------ -----    ------- ------ -----   ---- --  --------                  -----------
    8.5e-55  182.3   0.0    1.2e-54  181.8   0.0    1.2  1  OBAL001.B.00005.C001_472   # 516941 # 517834 # -1 # ;gc_cont=
    3.8e-45  150.9   0.0    7.3e-45  150.0   0.0    1.4  1  OBAL001.B.00005.C001_1416  # 1545830 # 1546807 # -1 # ;gc_con
    2.2e-43  145.2   0.0    4.2e-43  144.3   0.0    1.5  1  OBAL001.B.00005.C001_1793  # 1940269 # 1941171 # 1 # ;gc_cont
    2.9e-25   86.2   0.0    5.4e-25   85.3   0.0    1.4  1  OBAL001.B.00005.C001_388   # 418072 # 419283 # 1 # ;gc_cont=0
    2.5e-07   27.8   0.6      0.085    9.8   0.0    5.5  5  OBAL001.B.00005.C001_399   # 434671 # 440118 # -1 # ;gc_cont=


Domain annotation for each sequence (and alignments):
>> OBAL001.B.00005.C001_472  # 516941 # 517834 # -1 # ;gc_cont=0.496
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  181.8   0.0   2.8e-57   1.2e-54       2     172 ..     112     284 ..     111     285 .. 0.97

  Alignments for each domain:
  == domain 1  score: 181.8 bits;  conditional E-value: 2.8e-57
                               HHHHHHHHHHHHCCCT..HHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHHH CS
           Phage_integrase   2 vLtedeverllaalee..slsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkeil 84 
                                Lte++ve+ll a +   +l +rdr+++e+l++ glR+sEl++lr+++l+l++g+irv    +K++ker vpl+ee++++l+++l
  OBAL001.B.00005.C001_472 112 TLTEQDVEQLLSAPDCedTLGMRDRCMLEVLYASGLRVSELVALRLDELNLRQGVIRV---IGKGGKERLVPLGEEAIQWLERYL 193
                               69**********9877779***************************************...************************ PP

                               HH....HHTTSTTSBSSBECTSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHHHHH CS
           Phage_integrase  85 sdr...kkeaeerellfvskrgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHssisttkiY 166
                               +++       + ++++f+s rg++++++t+++++kr +++agi+k+l+pHtLRh fat+ll++G+dl+vvq llGH+++stt+iY
  OBAL001.B.00005.C001_472 194 QQArpiLLGLKLSDVVFPSLRGSQMTRQTFWHRIKRHAAAAGITKPLSPHTLRHAFATHLLNHGADLRVVQMLLGHTDLSTTQIY 278
                               ***999999999************************************************************************* PP

                               HCCSHH CS
           Phage_integrase 167 thvake 172
                               thva+ 
  OBAL001.B.00005.C001_472 279 THVANV 284
                               ***975 PP

>> OBAL001.B.00005.C001_1416  # 1545830 # 1546807 # -1 # ;gc_cont=0.456
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  150.0   0.0   1.6e-47   7.3e-45       2     172 ..     105     307 ..     104     308 .. 0.95

  Alignments for each domain:
  == domain 1  score: 150.0 bits;  conditional E-value: 1.6e-47
                                HHHHHHHHHHHHCCCTHHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHHHH CS
            Phage_integrase   2 vLtedeverllaaleeslsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkeils 85 
                                v+ ++ev+r+l+ +++    r+++++ ll++ glRi+E+l+lrvkd+d+dng+i v  +++K++k+r+  l+ +l++++k+ ++
  OBAL001.B.00005.C001_1416 105 VISANEVQRILQVMDT----RNQVIFALLYGAGLRINECLRLRVKDFDFDNGCITV--HDGKGGKSRNSLLPTRLIPAIKQLIE 182
                                67899***********....************************************..************************** PP

                                H.HHTTSTTS..........................BSSBEC...........TSSB..HHHHHHHHHHHHHHTT--.CC-HHH CS
            Phage_integrase  86 drkkeaeere..........................llfvsk...........rgkplsdstvnrafkravkeagie.keltpH 131
                                ++   +++++                          ++f+s             +++l ds +++a+k+av++ gi  k++t H
  OBAL001.B.00005.C001_1416 183 QARLIQQDDNlqgvgpslpfaldrkypsayrqaawmFVFPSStlcnhpyngklCRHHLHDSVARKALKAAVQKVGIVsKRVTCH 266
                                **555555556789999************************9****************************************** PP

                                HHHHHHHHHHHH----HHHHHHH----SHHHHHHHHCCSHH CS
            Phage_integrase 132 tLRhsfatallesGvdlkvvqkllGHssisttkiYthvake 172
                                t+ hsfat+ll++G d+++vq+llGH++++tt+iYthv  +
  OBAL001.B.00005.C001_1416 267 TFHHSFATHLLQAGRDIRTVQELLGHTDVKTTQIYTHVLGQ 307
                                *************************************9876 PP

>> OBAL001.B.00005.C001_1793  # 1940269 # 1941171 # 1 # ;gc_cont=0.505
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  144.3   0.0   9.4e-46   4.2e-43       2     172 ..     113     281 ..     112     282 .. 0.95

  Alignments for each domain:
  == domain 1  score: 144.3 bits;  conditional E-value: 9.4e-46
                                HHHHHHHHHHHHCCCT..HHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHH CS
            Phage_integrase   2 vLtedeverllaalee..slsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkei 83 
                                 L+ d + +ll+a ++    + rd++++el ++ glR sEl +l + d+dl  g +rv    +K++ker +p++ +++e+++++
  OBAL001.B.00005.C001_1793 113 LLDVDRATQLLDAPQGddFVQRRDHAMLELFYSSGLRLSELTGLNLPDIDLVAGQVRV---LGKGRKERFLPVGPQACEAIRQW 193
                                57888999999999996666789***********************************...*********************** PP

                                HHH.HHTTSTTSBSSBECTSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHHHHHH CS
            Phage_integrase  84 lsdrkkeaeerellfvskrgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHssisttkiYt 167
                                l  r++ a ++++lf+  rgk++ +stv+++++ a+++   + +l+pH+LRhsfa++lles  dl++vq+llGH++istt+iYt
  OBAL001.B.00005.C001_1793 194 LTVRGQVAAQDDALFLGVRGKRIAASTVRERVRLAGARELGQ-HLHPHMLRHSFASHLLESSQDLRAVQELLGHADISTTQIYT 276
                                **************************************5555.9**************************************** PP

                                CCSHH CS
            Phage_integrase 168 hvake 172
                                h++ +
  OBAL001.B.00005.C001_1793 277 HLDFQ 281
                                *9987 PP

>> OBAL001.B.00005.C001_388  # 418072 # 419283 # 1 # ;gc_cont=0.469
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   85.3   0.0   1.2e-27   5.4e-25       6     167 ..     212     375 ..     209     380 .. 0.92

  Alignments for each domain:
  == domain 1  score: 85.3 bits;  conditional E-value: 1.2e-27
                               HHHHHHHHCCCTHHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE.........CCSSSCC.EEEEE-HHHHHHHHH CS
           Phage_integrase   6 deverllaaleeslsirdrllvellleTglRisEllslrvkdldldngtirvp.......aretKtkke.rtvplseellevlke 82 
                               +eve++l+a++e    ++r +  +  +Tg+R++E+ +l+++ +d+++++i ++       ++ tKt  + r++p+   ++++lk+
  OBAL001.B.00005.C001_388 212 SEVETILQAVRE----DYRNYYLVRFYTGMRTGEIDGLKWEFVDFQRREILIRetliggeTEYTKTDGSqREIPMFGPVYDALKK 292
                               69**********....8899999999**********************************9999999888*************** PP

                               HHHH.HHTTSTTSBSSBECTSSB..HHHHHHHHH.HHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHH.HH CS
           Phage_integrase  83 ilsdrkkeaeerellfvskrgkplsdstvnrafk.ravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHssistt.ki 165
                               +++      + ++++f++++gkpl++++v +++   ++++ +++ + +p + Rh++at +l sG +++ v++ lGHs++++  k+
  OBAL001.B.00005.C001_388 293 QYEAT---GKLSKYVFCNREGKPLDHNNVTKRVWyPLLRALNLK-KRRPYQTRHTAATLFLASGENPEWVARTLGHSTTEMLfKV 373
                               **999...999*****************9877650555556666.9************************************9*9 PP

                               HH CS
           Phage_integrase 166 Yt 167
                               Y 
  OBAL001.B.00005.C001_388 374 YS 375
                               97 PP

>> OBAL001.B.00005.C001_399  # 434671 # 440118 # -1 # ;gc_cont=0.465
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !    9.8   0.0   0.00019     0.085      24     161 ..     599     756 ..     579     765 .. 0.55
   2 ?   -3.9   0.1       3.1   1.4e+03      87     118 ..     922     966 ..     908     969 .. 0.59
   3 !    8.5   0.0   0.00048      0.21      24     142 ..    1344    1480 ..    1316    1488 .. 0.69
   4 ?   -1.0   0.0      0.38   1.7e+02     149     170 ..    1520    1542 ..    1515    1544 .. 0.87
   5 ?   -0.1   0.0      0.21        92      41     156 ..    1614    1710 ..    1589    1717 .. 0.52

  Alignments for each domain:
  == domain 1  score: 9.8 bits;  conditional E-value: 0.00019
                               HHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCC.EEEEE-HHHHHHHHHHHHH.HHTTSTTSBSSBEC...... CS
           Phage_integrase  24 rllvellleTglRisEllslrvkdldldngtirvparetKtkke.rtvplseellevlkeilsdrkkeaeerellfvsk...... 101
                               ++++ l ++Tg+R       ++  ++++ ++++++ ++  + ++ r vpl++ +le l++++ +    a+  +l+ +++      
  OBAL001.B.00005.C001_399 599 YVVTALYAATGARPLRDPFESLTYFSFKYSCVFINDKSDDGLHSgRLVPLPKLALEILQRYVMHL---AKVADLIEPHRpelatk 680
                               55555555555555443333444456677777776666666666677777777777777766666...33333333333333332 PP

                               ....................TSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----SHH CS
           Phage_integrase 102 ....................rgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHssis 161
                                                    g +  + t  + ++  + +      l  + +Rh ++ +l ++Gv ++v++  +GH++  
  OBAL001.B.00005.C001_399 681 irlltqdnsdagmplfflldSGLQWHSMTSAEKLNCSLFD----WGLPANLFRHRYSQQLFKEGVHPEVIEGWMGHAERG 756
                               3333333444222222221122222222222222222222....23346779************************9865 PP

  == domain 2  score: -3.9 bits;  conditional E-value: 3.1
                               ..HHTTSTTSBSSBEC............TSSB..HHHHHHHHHHH CS
           Phage_integrase  87 r.kkeaeerellfvsk............rgkplsdstvnrafkra 118
                               r  ++ ee++l+ ++               ++ s++t+ + +k+a
  OBAL001.B.00005.C001_399 922 RlVQQSEEHSLITAEVaqamelyqelksWANRTSKQTFAGKIKKA 966
                               245556666666666655555555665556667777777777666 PP

  == domain 3  score: 8.5 bits;  conditional E-value: 0.00048
                                HHHHHHHHHHT--HHHHHC-BGGG.EECTTEEEEEE......CCSSSCC.EEEEE-HHHHHHHHHHHHH.........HHTTS CS
           Phage_integrase   24 rllvellleTglRisEllslrvkd.ldldngtirvp....aretKtkke.rtvplseellevlkeilsdr........kkeae 92  
                                 +++   +  glR +E+++l + d  ++ + t+ +     +r  K+  + r++pl  +l +  ++i+++         +k+ +
  OBAL001.B.00005.C001_399 1344 GFVLLATYRFGLRAQEAVGLLRRDwCQTTDYTWVLVqnsqYRTLKSASSrRAIPLLFRLSDIEQDIIERTlaryqsiaGKAIN 1426
                                3336667888********987777455555554443677677777776699*****999999999999887777777666666 PP

                                TTSBSSBECTSSB..HHHHHHHHHHHHHHTT--....CC-HHHHHHHHHHHHHH CS
           Phage_integrase   93 erellfvskrgkplsdstvnrafkravkeagie....keltpHtLRhsfatall 142 
                                ++ l+ +s  g++   + +  ++++a+ +  ++    +el  H  Rhsf  + +
  OBAL001.B.00005.C001_399 1427 RPILCEPSSTGNQPILTSIAPRISEALIQLLRNvignPELVLHHCRHSFYNRVA 1480
                                666777888888877888888888888885555556689999****99987765 PP

  == domain 4  score: -1.0 bits;  conditional E-value: 0.38
                                HHHHHH----SHHHH.HHHHCCS CS
           Phage_integrase  149 kvvqkllGHssistt.kiYthva 170 
                                +++++l+GH+  st  k+Y h+a
  OBAL001.B.00005.C001_399 1520 MALARLMGHKFPSTGlKNYFHLA 1542
                                6899**********999999986 PP

  == domain 5  score: -0.1 bits;  conditional E-value: 0.21
                                HC-BGGGEECTTEEEEEE..CCSSSCC....EEEEE-.....HHHHHHHHHHHHH.HHTTSTTSBSSBECTSSB..HHHHHHH CS
           Phage_integrase   41 lslrvkdldldngtirvparetKtkke....rtvpls.....eellevlkeilsdrkkeaeerellfvskrgkplsdstvnra 114 
                                ++l+ k + l ++++ +  ++t ++ +     +v ls     ++lle+                          + d++ +r 
  OBAL001.B.00005.C001_399 1614 MRLEPKYVALLQKVFET--TNTRMRFSassdKRVKLSgelcpNALLEA--------------------------IPDQAWQRL 1668
                                55555555555555555..33333333334444444433333333333..........................334444444 PP

                                HHHHHHHTT--.CC-HHHHHHHHHHHHHHH----HHHHHHH-- CS
           Phage_integrase  115 fkravkeagie.keltpHtLRhsfatallesGvdlkvvqkllG 156 
                                ++ a +++g +  +l   +LRh     +l s+ + +v++++lG
  OBAL001.B.00005.C001_399 1669 LHCAKENTGSTdLKLEDDELRHLQELPYLISQ-NRQVLMEQLG 1710
                                44555555544445555555555555555554.5555555555 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (173 nodes)
Target sequences:                         2218  (714377 residues searched)
Passed MSV filter:                        61  (0.0275023); expected 44.4 (0.02)
Passed bias filter:                       46  (0.0207394); expected 44.4 (0.02)
Passed Vit filter:                         7  (0.003156); expected 2.2 (0.001)
Passed Fwd filter:                         5  (0.00225428); expected 0.0 (1e-05)
Initial search space (Z):               2218  [actual number of targets]
Domain search space  (domZ):               5  [number of targets reported over threshold]
# CPU time: 0.01u 0.00s 00:00:00.01 Elapsed: 00:00:00.01
# Mc/sec: 12358.72
//
[ok]
