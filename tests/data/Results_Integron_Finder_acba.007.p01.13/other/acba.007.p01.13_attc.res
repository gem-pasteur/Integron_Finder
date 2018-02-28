# cmsearch :: search CM(s) against a sequence database
# INFERNAL 1.1.2 (July 2016)
# Copyright (C) 2016 Howard Hughes Medical Institute.
# Freely distributed under a BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query CM file:                         /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/data/Models/attc_4.cm
# target sequence database:              /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/tests/data/Replicons/acba.007.p01.13.fst
# output directed to file:               /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_attc.res
# tabular output of hits:                /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_attc_table.res
# sequence reporting threshold:          E-value <= 10
# number of worker threads:              1 [--cpu]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       attC_4  [CLEN=47]
Hit scores:
 rank     E-value  score  bias  sequence         start    end   mdl trunc   gc  description
 ----   --------- ------ -----  --------------- ------ ------   --- ----- ----  -----------
  (1) !     1e-09   46.4   0.0  ACBA.007.P01_13  17884  17825 -  cm    no 0.55  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid 
  (2) !   1.1e-07   38.4   0.0  ACBA.007.P01_13  19726  19618 -  cm    no 0.59  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid 
  (3) !    0.0001   26.6   0.0  ACBA.007.P01_13  19149  19080 -  cm    no 0.69  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid 


Hit alignments:
>> ACBA.007.P01_13  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid pMDR-ZJ06, complete
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (1) !     1e-09   46.4   0.0  cm        1       47 []       17884       17825 - .. 1.00    no 0.55

                                                                                            NC
                        <<<<<<<--------<<<-<<<<....................>>>>>>>---------->>>>>>> CS
           attC_4     1 GcCUAACAAgUCAUUGuUCAAGc....................gCUUAaCUCgGcCAUUCGUUAGgC 47   
                         :CUAACAA+UC   GUUCAAGC                    GCUUAACUC+++C    GUUAG: 
  ACBA.007.P01_13 17884 GUCUAACAAUUC---GUUCAAGCcgacgccgcuucgcggcgcgGCUUAACUCAAGC----GUUAGAU 17825
                        ************...*****************************************....******* PP
                        [Rsec=]========[=Lsec=]....................[Lprim]==========[Rprim] RF

>> ACBA.007.P01_13  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid pMDR-ZJ06, complete
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (2) !   1.1e-07   38.4   0.0  cm        1       47 []       19726       19618 - .. 1.00    no 0.59

                        v                                                                                          NC
                        <<<<<<<--------<<<-<<<<................................................................... CS
           attC_4     1 GcCUAACAAgUCAUUGuUCAAGc................................................................... 23   
                         CC AACAA+UC   GUUCAAGC                                                                   
  ACBA.007.P01_13 19726 ACCUAACAAUUC---GUUCAAGCcgagaucgcuucgcggccgcggaguuguucggaaaaauugucacaacgccgcggccgcaaagcgcuc 19640
                        ************...*************************************************************************** PP
                        [Rsec=]========[=Lsec=]................................................................... RF

                                                 v NC
                        ..>>>>>>>---------->>>>>>> CS
           attC_4    24 ..gCUUAaCUCgGcCAUUCGUUAGgC 47   
                          GCUUAACUC+G+C    GUU GG 
  ACBA.007.P01_13 19639 cgGCUUAACUCAGGC----GUUGGGC 19618
                        ***************....******* PP
                        ..[Lprim]==========[Rprim] RF

>> ACBA.007.P01_13  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid pMDR-ZJ06, complete
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (3) !    0.0001   26.6   0.0  cm        1       47 []       19149       19080 - .. 1.00    no 0.69

                        vv v                                                                     v vv NC
                        <<<<<<<--------<<<-<<<<..............................>>>>>>>---------->>>>>>> CS
           attC_4     1 GcCUAACAAgUCAUUGuUCAAGc..............................gCUUAaCUCgGcCAUUCGUUAGgC 47   
                          C AACA G+C   G:UCAAGC                              GCUUA:CU GG+C    GUU G  
  ACBA.007.P01_13 19149 GCCCAACAUGGC---GCUCAAGCcgaccggccagcccugcgggcuguccgucgGCUUAGCUAGGGC----GUUAGAG 19080
                        ************...***************************************************....******* PP
                        [Rsec=]========[=Lsec=]..............................[Lprim]==========[Rprim] RF



Internal CM pipeline statistics summary:
----------------------------------------
Query model(s):                                                  1  (47 consensus positions)
Target sequences:                                                1  (40602 residues searched)
Target sequences re-searched for truncated hits:                 1  (3580 residues re-searched)
Windows   passing  local HMM SSV           filter:              19  (0.3933); expected (0.35)
Windows   passing  local HMM Viterbi       filter:                  (off)
Windows   passing  local HMM Viterbi  bias filter:                  (off)
Windows   passing  local HMM Forward       filter:               5  (0.1125); expected (0.02)
Windows   passing  local HMM Forward  bias filter:               5  (0.1125); expected (0.02)
Windows   passing glocal HMM Forward       filter:               3  (0.09221); expected (0.02)
Windows   passing glocal HMM Forward  bias filter:               3  (0.09221); expected (0.02)
Envelopes passing glocal HMM envelope defn filter:               4  (0.009574); expected (0.02)
Envelopes passing  local CM  CYK           filter:               3  (0.005409); expected (0.0001)
Total CM hits reported:                                          3  (0.005409); includes 0 truncated hit(s)

# CPU time: 0.06u 0.00s 00:00:00.06 Elapsed: 00:00:00.05
//
[ok]
