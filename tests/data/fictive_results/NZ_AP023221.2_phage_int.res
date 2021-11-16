# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/integron_finder/data/Models/phage-int.hmm
# target sequence database:        /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_NZ_AP023221/tmp_NZ_AP023221.2/NZ_AP023221.2.prt
# output directed to file:         /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_NZ_AP023221/tmp_NZ_AP023221.2/NZ_AP023221.2_phage_int.res
# per-seq hits tabular output:     /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_NZ_AP023221/tmp_NZ_AP023221.2/NZ_AP023221.2_phage_int_table.res
# model-specific thresholding:     GA cutoffs
# number of worker threads:        1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       Phage_integrase  [M=173]
Accession:   PF00589.16
Description: Phage integrase family
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence          Description
    ------- ------ -----    ------- ------ -----   ---- --  --------          -----------
      5e-30   97.9   0.0    7.3e-30   97.4   0.0    1.3  1  NZ_AP023221.2_109  # 93412 # 94152 # -1 # ;gc_cont=0.588
    5.5e-29   94.5   0.1    7.8e-29   94.1   0.1    1.2  1  NZ_AP023221.2_106  # 91461 # 92201 # -1 # ;gc_cont=0.591
    1.3e-06   21.7   0.0    1.7e-06   21.3   0.0    1.1  0  NZ_AP023221.2_85   # 73154 # 73432 # -1 # ;gc_cont=0.573


Domain annotation for each sequence (and alignments):
>> NZ_AP023221.2_109  # 93412 # 94152 # -1 # ;gc_cont=0.588
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   97.4   0.0   1.4e-31   7.3e-30       2     171 ..      44     217 ..      43     219 .. 0.90

  Alignments for each domain:
  == domain 1  score: 97.4 bits;  conditional E-value: 1.4e-31
                        HHHHHHHHHHHHCCCTHHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE.........CCSSSCC..EEEEE-.HHHHHHHHHH CS
    Phage_integrase   2 vLtedeverllaaleeslsirdrllvellleTglRisEllslrvkdldldngtirvp.......aretKtkke..rtvpls.eellevlkei 83 
                        +L + ev++ll ++ +    r+++l+ +l++Tg+Ri+E++ l+ + +dld  + +v+       ar++   k+  r vpl+ ++++ +++ +
  NZ_AP023221.2_109  44 YLLAPEVSALLFYMPD---QRHHMLFATLWNTGMRIGEARMLTPESFDLDGARPFVRvlsekvrARRGRPPKDevRLVPLTdASFVRQMESW 132
                        78899***********...9*************************************88888887777777777899999988889999998 PP

                        HHH.HHTTSTTSBSSBEC..TSSB..HHHHHHHHHHHHHHTT--.....CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHHHHHHC CS
    Phage_integrase  84 lsdrkkeaeerellfvsk..rgkplsdstvnrafkravkeagie.....keltpHtLRhsfatallesGvdlkvvqkllGHssisttkiYth 168
                        +              +++  + +p++d+t+++++k+avk+a+ +      ++tpHt+Rhs+++++l +  + kv+q l GH++ +++++Yt+
  NZ_AP023221.2_109 133 MVTT----------RPRRrePLWPVTDETMRNWLKQAVKRAEADgvhfsIPVTPHTFRHSYIMHMLYHRQPRKVIQALAGHRDPRSMEVYTR 214
                        8888..........233344667899**************9999999999*****************************************9 PP

                        CSH CS
    Phage_integrase 169 vak 171
                        v++
  NZ_AP023221.2_109 215 VFA 217
                        986 PP

>> NZ_AP023221.2_106  # 91461 # 92201 # -1 # ;gc_cont=0.591
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   94.1   0.1   1.5e-30   7.8e-29       2     171 ..      44     217 ..      43     219 .. 0.88

  Alignments for each domain:
  == domain 1  score: 94.1 bits;  conditional E-value: 1.5e-30
                        HHHHHHHHHHHHCCCTHHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE.........CCSSSCC..EEEEE-.HHHHHHHHHH CS
    Phage_integrase   2 vLtedeverllaaleeslsirdrllvellleTglRisEllslrvkdldldngtirvp.......aretKtkke..rtvpls.eellevlkei 83 
                        +L + ev++ll ++ +    r+++l+ +l++Tg+Ri+E++ l+ + +dld  + +v+       ar++   k+  r vpl+  +++ +++ +
  NZ_AP023221.2_106  44 YLLAPEVSALLFYMPD---QRHHMLFATLWNTGMRIGEARMLTPESFDLDGVRPFVRilsekvrARRGRPPKDevRLVPLTdISYVRQMESW 132
                        78899***********...9*************************************77777776666666667788888877788888888 PP

                        HHH.HHTTSTTSBSSBEC.....TSSB..HHHHHHHHHHHHHHTT--.....CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHHHH CS
    Phage_integrase  84 lsdrkkeaeerellfvsk.....rgkplsdstvnrafkravkeagie.....keltpHtLRhsfatallesGvdlkvvqkllGHssisttki 165
                        +                +     + + ++d+t+++++k+avk+a+ +      ++tpHt+Rhs+++++l +  + kv+q l GH++ +++++
  NZ_AP023221.2_106 133 MITT-------------RprrraPLWAVTDETMRNWLKQAVKRAEADgvhfsIPVTPHTFRHSYIMHMLYHRQPRKVIQALAGHRDPRSMEV 211
                        7777.............222233456689**************9999999999*************************************** PP

                        HHCCSH CS
    Phage_integrase 166 Ythvak 171
                        Yt+v++
  NZ_AP023221.2_106 212 YTRVFA 217
                        **9986 PP

>> NZ_AP023221.2_85  # 73154 # 73432 # -1 # ;gc_cont=0.573
   [No individual domains that satisfy reporting thresholds (although complete target did)]



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (173 nodes)
Target sequences:                          155  (36040 residues searched)
Passed MSV filter:                         7  (0.0451613); expected 3.1 (0.02)
Passed bias filter:                        5  (0.0322581); expected 3.1 (0.02)
Passed Vit filter:                         3  (0.0193548); expected 0.2 (0.001)
Passed Fwd filter:                         3  (0.0193548); expected 0.0 (1e-05)
Initial search space (Z):                155  [actual number of targets]
Domain search space  (domZ):               3  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# Mc/sec: inf
//
[ok]
