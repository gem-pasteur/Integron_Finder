# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/data/Models/phage-int.hmm
# target sequence database:        /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_ACBA.0917.00019/tmp_ACBA.0917.00019.0001/ACBA.0917.00019.0001.prt
# output directed to file:         /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_ACBA.0917.00019/tmp_ACBA.0917.00019.0001/ACBA.0917.00019.0001_phage_int.res
# per-seq hits tabular output:     /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_ACBA.0917.00019/tmp_ACBA.0917.00019.0001/ACBA.0917.00019.0001_phage_int_table.res
# number of worker threads:        1
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
    7.5e-41  137.7   0.8    2.6e-40  136.0   0.0    2.0  2  ACBA.0917.00019.i0001_03058  915 xerC_2 | Tyrosine recombinas
    1.4e-28   97.8   0.0    3.8e-28   96.4   0.0    1.8  2  ACBA.0917.00019.i0001_00980  1284 intS_1 | Prophage integrase
    2.1e-23   80.9   0.2    2.1e-23   80.9   0.2    2.0  2  ACBA.0917.00019.i0001_02762  987 xerD_3 | Tyrosine recombinas
    1.3e-17   62.1   0.0    2.3e-17   61.3   0.0    1.4  1  ACBA.0917.00019.i0001_03204  1245 intA | Prophage integrase I
    5.6e-17   60.0   0.0    9.4e-15   52.8   0.0    2.3  2  ACBA.0917.00019.i0001_01335  1206 intS_2 | Prophage integrase
    2.9e-15   54.4   0.0    7.9e-15   53.0   0.0    1.8  1  ACBA.0917.00019.i0001_01431  927 xerC_1 | Tyrosine recombinas
    2.1e-13   48.4   0.0    4.6e-13   47.3   0.0    1.6  1  ACBA.0917.00019.i0001_00654  1149 NA | hypothetical protein |
    2.8e-07   28.4   0.0    7.4e-07   27.1   0.0    1.7  1  ACBA.0917.00019.i0001_01756  1233 NA | hypothetical protein |


Domain annotation for each sequence (and alignments):
>> ACBA.0917.00019.i0001_00298  1035 xerD_1 | Tyrosine recombinase XerD | NA | similar to AA sequence:UniProtKB:P9WF33
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?   -2.3   0.0         2   7.7e+02       3      26 ..      56      83 ..      54      88 .. 0.66
   2 !  219.5   0.0   1.4e-68   5.5e-66       2     171 ..     116     317 ..     115     319 .. 0.99

  Alignments for each domain:
  == domain 1  score: -2.3 bits;  conditional E-value: 2
                                 HHHHHHHHHHHCCCT....HHHHHHHHH CS
              Phage_integrase  3 Ltedeverllaalee....slsirdrll 26
                                 L ++eve++l +l +    s s+++++l
  ACBA.0917.00019.i0001_00298 56 LGSSEVEAFLSWLANerkvSVSTHRQAL 83
                                 6678999999999985555555555555 PP

  == domain 2  score: 219.5 bits;  conditional E-value: 1.4e-68
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
   1 ?   -3.0   0.0       3.3   1.3e+03      59      83 ..      73      99 ..      59     107 .. 0.64
   2 !  171.4   0.1   8.7e-54   3.4e-51       3     172 ..     123     293 ..     121     294 .. 0.98

  Alignments for each domain:
  == domain 1  score: -3.0 bits;  conditional E-value: 3.3
                                 ..CCSSSCCEEEEE-..HHHHHHHHHH CS
              Phage_integrase 59 aretKtkkertvpls..eellevlkei 83
                                 ++ +K+ ++++  ls  ++++++l+e+
  ACBA.0917.00019.i0001_00338 73 TKVGKSPRSIARCLSalRQFYKFLREQ 99
                                 445666666777777777777777775 PP

  == domain 2  score: 171.4 bits;  conditional E-value: 8.7e-54
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

>> ACBA.0917.00019.i0001_03058  915 xerC_2 | Tyrosine recombinase XerC | NA | similar to AA sequence:UniProtKB:P44818
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?   -0.2   0.2      0.45   1.7e+02      77     102 ..       5      45 ..       2      83 .. 0.43
   2 !  136.0   0.0   6.6e-43   2.6e-40      18     172 ..     134     288 ..     113     289 .. 0.92

  Alignments for each domain:
  == domain 1  score: -0.2 bits;  conditional E-value: 0.45
                                  HHHHHHHHHH.HHTTSTTSBSSBEC...............T CS
              Phage_integrase  77 levlkeilsdrkkeaeerellfvsk...............r 102
                                  l++l  +l++rk ++++++++ +                 +
  ACBA.0917.00019.i0001_03058   5 LQLLSMWLKERKIQNQSEHTITAYErdvrsflefcelkkvD 45 
                                  55666666666444444444444432222222222222222 PP

  == domain 2  score: 136.0 bits;  conditional E-value: 6.6e-43
                                  HHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHHHHH....HHTTSTTSB CS
              Phage_integrase  18 slsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkeilsdr...kkeaeerel 96 
                                  +l +rd++++ell++ glR +El +l++kd+d++++ +r+    +K++k r+vp+++++ e+l ++l+     k   +++  
  ACBA.0917.00019.i0001_03058 134 QLWLRDKAMLELLYSSGLRLAELQGLTIKDIDFNRQLVRI---TGKGNKTRIVPFGKKAKESLLNWLKIYniwKGHFDQNAS 212
                                  56689***********************************...*************9999999999998889888899999* PP

                                  SSBECTSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHHHHHHCCSHH CS
              Phage_integrase  97 lfvskrgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHssisttkiYthvake 172
                                  +f+s+rg  l+ +++++++k  +++ag++ +l+pH LRh fa+++l s  dl+ vq++lGHs++stt+iYth + +
  ACBA.0917.00019.i0001_03058 213 VFISQRGGALTPRQIEKRVKLQAQRAGVNVDLHPHLLRHCFASHMLSSSGDLRSVQEMLGHSNLSTTQIYTHIDFD 288
                                  ************************************************************************9875 PP

>> ACBA.0917.00019.i0001_00980  1284 intS_1 | Prophage integrase IntS | NA | similar to AA sequence:UniProtKB:P37326
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?   -3.8   0.0       5.6   2.2e+03      62      87 ..      19      44 ..      15      72 .. 0.58
   2 !   96.4   0.0   9.9e-31   3.8e-28       4     171 ..     236     397 ..     233     399 .. 0.93

  Alignments for each domain:
  == domain 1  score: -3.8 bits;  conditional E-value: 5.6
                                 CSSSCCEEEEE-HHHHHHHHHHHHH. CS
              Phage_integrase 62 tKtkkertvplseellevlkeilsdr 87
                                 +K+ +  + pls+++l  l+   ++ 
  ACBA.0917.00019.i0001_00980 19 MKRTEIKRRPLSDTVLANLEPESKEY 44
                                 45555577788877777776655544 PP

  == domain 2  score: 96.4 bits;  conditional E-value: 9.9e-31
                                  HHHHHHHHHHCCCTHHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHHHH CS
              Phage_integrase   4 tedeverllaaleeslsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkeils 85 
                                   e+e+++ll a+++   ++ r+ ++ll+   +R  El++ +++++dl++g++++pa+++K+++e+ vpl+++++  l e   
  ACBA.0917.00019.i0001_00980 236 SEQELPALLRAINSYPTMDVRMGLQLLAMLFCRPTELREAKWQEFDLNQGIWNIPAERMKKRREHVVPLPRQAITILNELKT 317
                                  57899*********888899999********************************************************777 PP

                                  H.HHTTSTTSBSSBEC..TSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHHHH CS
              Phage_integrase  86 drkkeaeerellfvsk..rgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHssisttki 165
                                   +     ++e+lf+s+  + kp+sd+ +  a++r+++e     + tpH +Rh ++t l ++G+d + ++  l H +  ++ +
  ACBA.0917.00019.i0001_00980 318 YE----TNSEYLFPSRsdKSKPKSDTVFIMALRRMGYE----GRQTPHGFRHIASTLLNNRGFDERHIEAALAHVKDGVAGV 391
                                  77....99********9999999999999999999999....556************************************* PP

                                  HHCCSH CS
              Phage_integrase 166 Ythvak 171
                                  Y+++++
  ACBA.0917.00019.i0001_00980 392 YNKAQY 397
                                  **9876 PP

>> ACBA.0917.00019.i0001_02762  987 xerD_3 | Tyrosine recombinase XerD | NA | similar to AA sequence:UniProtKB:P0A8P8
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?    0.2   0.3      0.34   1.3e+02      59     110 ..      38      95 ..      17     125 .. 0.65
   2 !   80.9   0.2   5.4e-26   2.1e-23       2     171 ..     161     320 ..     160     322 .. 0.85

  Alignments for each domain:
  == domain 1  score: 0.2 bits;  conditional E-value: 0.34
                                  ..CCSSSCCEEEEE-.HHHHHHHHHHHHH......HHTTSTTSBSSBECTSSB..HHH CS
              Phage_integrase  59 aretKtkkertvpls.eellevlkeilsdr.....kkeaeerellfvskrgkplsdst 110
                                  +ret  +++    l+  ++le+++++ls+      + e+e + + f +++ k+l++++
  ACBA.0917.00019.i0001_02762  38 KRETQLREQSHGKLPdHSFLEAIERYLSEVsvkkkTHENEVKRMAFFKREYKKLCQKQ 95 
                                  5555666666666666788888888888886643333333333555555555555554 PP

  == domain 2  score: 80.9 bits;  conditional E-value: 5.4e-26
                                  HHHHHHHHHHHHCCCT....HHHHHHHHH..HHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHH CS
              Phage_integrase   2 vLtedeverllaalee....slsirdrll..vellleTglRisEllslrvkdldldngtirvparetKtkkertvplseell 77 
                                  ++ +de++rl  a++     + +  ++++  + ++ eT++R +E+++l+++ + l++++  +   etK++++r+vpls++++
  ACBA.0917.00019.i0001_02762 161 RIAQDEIDRLCLAANWdnnvPVNSTQQIIiaFLFAIETAMRAGEIVGLTWDRVYLKDRYLVL--NETKNGTKRNVPLSKRAV 240
                                  56777888877777663445223333333458999***************************..****************** PP

                                  HHHHHHHHH.HHTTSTTSBSSBECTSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----S CS
              Phage_integrase  78 evlkeilsdrkkeaeerellfvskrgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHss 159
                                  e+l   l+       +++ +f+       +++ + + ++++  + +i+ +l++H+ Rh+++t+l+++   +  ++++ GH++
  ACBA.0917.00019.i0001_02762 241 ELLTL-LKGL-----DKKQVFTC------NSQSFDTLWRKLRDRCQIT-DLHFHDTRHEACTRLARKL-EVLDLARMIGHKD 308
                                  **987.4444.....58899955......59*****************.****************986.9************ PP

                                  HHHHHHHHCCSH CS
              Phage_integrase 160 isttkiYthvak 171
                                  +++ ++Y ++++
  ACBA.0917.00019.i0001_02762 309 LRSLMVYYNATA 320
                                  *******99876 PP

>> ACBA.0917.00019.i0001_03204  1245 intA | Prophage integrase IntA | NA | similar to AA sequence:UniProtKB:P32053
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   61.3   0.0   5.9e-20   2.3e-17      20     169 ..     231     383 ..     211     387 .. 0.86

  Alignments for each domain:
  == domain 1  score: 61.3 bits;  conditional E-value: 5.9e-20
                                  HHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE....CCSSSCCEEEEE-HHHHHHHHHHHHH.HHTTSTTSBSSB CS
              Phage_integrase  20 sirdrllvellleTglRisEllslrvkdldldngtirvp..aretKtkkertvplseellevlkeilsdrkkeaeerellfv 99 
                                     ++++++++ +  +R +El+  ++ d+dl+   +r    +++ Kt+ ++ vp+++++  +l +i +  ++ + e+e++f+
  ACBA.0917.00019.i0001_03204 231 TFITQMALKIAPYVFVRPGELRYAKWPDIDLEIDLWRYTppKTKNKTGVQHLVPIPRQVKVLLLKIKELTYDPDGESEYVFP 312
                                  5557889999999************************9855666677778******88888888866666999********* PP

                                  EC..TSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH.----HHHHHHH----SHHHH.HHHHCC CS
              Phage_integrase 100 sk..rgkplsdstvnrafkravkeagiekeltpHtLRhsfatalle.sGvdlkvvqkllGHssistt.kiYthv 169
                                  s   + kp+s++t+n+a++r+++     ++++ H +R s+ t l e  +++++ ++++l H+   +  ++Y++ 
  ACBA.0917.00019.i0001_03204 313 SMtsKLKPMSENTINQALRRLGYT--S-EQVCGHGFRASARTILEEvLNYPIEIIEQQLAHKVKDMHgRAYNRT 383
                                  *9999******************9..3.478*************995699***************999999876 PP

>> ACBA.0917.00019.i0001_01335  1206 intS_2 | Prophage integrase IntS | NA | similar to AA sequence:UniProtKB:P37326
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?    5.0   0.0     0.011       4.3     102     143 ..      91     132 ..      62     139 .. 0.87
   2 !   52.8   0.0   2.4e-17   9.4e-15       1     157 [.     205     366 ..     205     377 .. 0.82

  Alignments for each domain:
  == domain 1  score: 5.0 bits;  conditional E-value: 0.011
                                  TSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH CS
              Phage_integrase 102 rgkplsdstvnrafkravkeagiekeltpHtLRhsfatalle 143
                                  ++k + +st+++ f  +++ + i+k+ +  ++R+sf  +  +
  ACBA.0917.00019.i0001_01335  91 QQKYIEASTLKEVFDDWYESYCIKKKTSAKDIRRSFEHHVFD 132
                                  678899******************************988765 PP

  == domain 2  score: 52.8 bits;  conditional E-value: 2.4e-17
                                  -HHHHHHHHHHHHCCCT.HHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHH CS
              Phage_integrase   1 kvLtedeverllaalee.slsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlk 81 
                                  +vLt +e+  + +a++e ++ ++++++++l+l  g+R +El++   +d+dl++++++vp ++ K +k++   + + +l+ ++
  ACBA.0917.00019.i0001_01335 205 RVLTDEEITMVWKAIDEsKALLKNKIFLKLCLMYGCRNGELRKALKSDFDLKRKVWIVPVENNKVGKKTGREIVRPILPEME 286
                                  68***************88899*************************************99999999988888888888888 PP

                                  HHHHH.HHTTSTTSBSSBEC.....TSSB..HHHHHHHHHHHHHHTT--.CC-HHHHHHHHHHHHHHH----HHHHHHH--- CS
              Phage_integrase  82 eilsdrkkeaeerellfvsk.....rgkplsdstvnrafkravkeagie.keltpHtLRhsfatallesGvdlkvvqkllGH 157
                                  e + ++  + +++e+ +++       g+  s s + + ++r+ ++++   ++++ H+LR+++ t+++    + +v + + GH
  ACBA.0917.00019.i0001_01335 287 ELIVEA-MSLNDSEYFLTNDddvtpMGHGSSNSLAANVMERLRRHYNFHmPHWSLHDLRRTARTNFSAFT-SRDVAELMIGH 366
                                  888888.55556666666667777777777777777788888888888799**************97543.44444444455 PP

>> ACBA.0917.00019.i0001_01431  927 xerC_1 | Tyrosine recombinase XerC | NA | protein motif:HAMAP:MF_01808
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   53.0   0.0     2e-17   7.9e-15       5     171 ..     138     298 ..     134     300 .. 0.82

  Alignments for each domain:
  == domain 1  score: 53.0 bits;  conditional E-value: 2e-17
                                  HHHHHHHHHCCCT.....HHHHHHHHH..HHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-...HHH CS
              Phage_integrase   5 edeverllaalee.....slsirdrll..vellleTglRisEllslrvkdldldngtirvparetKtkkertvpls...eel 76 
                                  ++ ++++la l++     + + r+r++  + ++leT++R +E+ls + + +      ir+   +tK+++ r vpl+   +el
  ACBA.0917.00019.i0001_01431 138 QNYIDKVLAGLDYewgkvPVQPRHRVAwsFLFALETAIRKGEILSVEKSLI--FPDFIRL--LDTKNGTTRDVPLTtkaKEL 215
                                  55666666666666666666778888878999*************999986..6667777..7*************666666 PP

                                  HHHHHHHHHH.HHTTSTTSBSSBECTSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH---.-HHHHHHH--- CS
              Phage_integrase  77 levlkeilsdrkkeaeerellfvskrgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGv.dlkvvqkllGH 157
                                  l++l +  +d+              r  pl++++++ +++r +++ g++  +t+H+ Rh+++t+++     +++ ++k+ GH
  ACBA.0917.00019.i0001_01431 216 LSWLPDDPDDN--------------RMVPLTSNAFRLIWQRNLRRVGLDGVITFHDTRHEAITRFVHDYRlPVEILAKITGH 283
                                  66666655555..............6667889*********************************987544*********** PP

                                  -SHHHH.HHHHCCSH CS
              Phage_integrase 158 ssistt.kiYthvak 171
                                  ++is+  + Y + ++
  ACBA.0917.00019.i0001_01431 284 KTISVLvNTYYNPTA 298
                                  *****9999998776 PP

>> ACBA.0917.00019.i0001_00654  1149 NA | hypothetical protein | NA | NA
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   47.3   0.0   1.2e-15   4.6e-13      20     170 ..     213     364 ..     190     367 .. 0.86

  Alignments for each domain:
  == domain 1  score: 47.3 bits;  conditional E-value: 1.2e-15
                                  HHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSS.SCCEEEEE-HHHHHHHHHHHHH.HHTTSTTSBSSBE CS
              Phage_integrase  20 sirdrllvellleTglRisEllslrvkdldldngtirvparetKt.kkertvplseellevlkeilsdrkkeaeerellfvs 100
                                  ++ ++l++ ++     R +El +l ++d d  ++ ++v   + K+ + +   + s e+le  k+i++   +++ ++ +l + 
  ACBA.0917.00019.i0001_00654 213 KYPMHLIIWFAIFSCRREAELTRLWLQDYDSYHSAWKV--YDLKNpNGSKGNHKSFEVLEPCKRIIELLLDNEVRSRMLQLG 292
                                  566778888888888999********************..666652677888889999999999999998889999999999 PP

                                  C...TSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHHHHHHCCS CS
              Phage_integrase 101 k...rgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHssisttkiYthva 170
                                         pl+ +t+ + f+ a++  gie +l++H+LRh+ +t+l+e+ + + ++qk+  H s  + ++Y  v 
  ACBA.0917.00019.i0001_00654 293 YdkqLLLPLNPKTIGKEFRDACRMLGIE-DLHFHDLRHEGCTRLAEQSFTIPEIQKVSLHDSWGSLQRYVSVK 364
                                  844444599*******************.****************************************9876 PP

>> ACBA.0917.00019.i0001_01756  1233 NA | hypothetical protein | NA | NA
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   27.1   0.0   1.9e-09   7.4e-07      17     139 ..     236     350 ..     216     358 .. 0.78

  Alignments for each domain:
  == domain 1  score: 27.1 bits;  conditional E-value: 1.9e-09
                                  T.HHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHHHHH.HHTTSTTSBS CS
              Phage_integrase  17 e.slsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkeilsdrkkeaeerell 97 
                                  + +++  ++ ++ +l +Tgl  +E+ sl++k++dl++gt      +  ++   ++ +++ l +++k+   ++     + e++
  ACBA.0917.00019.i0001_01756 236 RgQENETNKDFLLTLILTGLFRNECESLHWKNIDLEEGTLSF--LNPYNNVYYKIYMGNFLWHLMKKRRIQN-----KGEWV 310
                                  3677888999999*****************************..8888888888888888888888855555.....9**** PP

                                  SBECTSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHH CS
              Phage_integrase  98 fvskrgkplsdstvnrafkravkeagiekeltpHtLRhsfat 139
                                  f+s + +     ++++  k++ ++ +++  +t+ +LR++f  
  ACBA.0917.00019.i0001_01756 311 FPSVKSESGHIINISKFRKKINEQCKLS--FTFQDLRRTFYF 350
                                  ****666666666666666666666667..*********976 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (173 nodes)
Target sequences:                         3870  (1172415 residues searched)
Passed MSV filter:                       129  (0.0333333); expected 77.4 (0.02)
Passed bias filter:                       95  (0.0245478); expected 77.4 (0.02)
Passed Vit filter:                        18  (0.00465116); expected 3.9 (0.001)
Passed Fwd filter:                        10  (0.00258398); expected 0.0 (1e-05)
Initial search space (Z):               3870  [actual number of targets]
Domain search space  (domZ):              10  [number of targets reported over threshold]
# CPU time: 0.02u 0.00s 00:00:00.02 Elapsed: 00:00:00.02
# Mc/sec: 10141.39
//
[ok]
