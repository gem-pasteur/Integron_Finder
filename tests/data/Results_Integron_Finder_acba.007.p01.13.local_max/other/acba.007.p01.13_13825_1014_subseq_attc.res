# cmsearch :: search CM(s) against a sequence database
# INFERNAL 1.1.2 (July 2016)
# Copyright (C) 2016 Howard Hughes Medical Institute.
# Freely distributed under a BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query CM file:                         /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/data/Models/attc_4.cm
# target sequence database:              /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_subseq.fst
# database size is set to:               0.0 Mb
# output directed to file:               /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_13825_1014_subseq_attc.res
# tabular output of hits:                /home/bneron/Projects/GEM/Integron_Finder/src/Integron_Finder/Results_Integron_Finder_acba.007.p01.13/other/acba.007.p01.13_13825_1014_subseq_attc_table.res
# sequence reporting threshold:          E-value <= 10
# Max sensitivity mode:                  on [all heuristic filters off]
# search bottom-strand only:             on
# truncated hit detection:               off [due to --max]
# number of worker threads:              1 [--cpu]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       attC_4  [CLEN=47]
Hit scores:
 rank     E-value  score  bias  sequence         start    end   mdl trunc   gc  description
 ----   --------- ------ -----  --------------- ------ ------   --- ----- ----  -----------
  (1) !   5.1e-10   46.4   0.0  ACBA.007.P01_13   4059   4000 -  cm    no 0.55  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid 
  (2) !   5.5e-08   38.4   0.0  ACBA.007.P01_13   5901   5793 -  cm    no 0.59  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid 
  (3) !   3.3e-05   27.3   0.0  ACBA.007.P01_13   5324   5254 -  cm    no 0.69  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid 
 ------ inclusion threshold ------
  (4) ?     0.016   16.7   0.0  ACBA.007.P01_13   4924   4856 -  cm    no 0.61  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid 
  (5) ?       2.5    8.0   0.0  ACBA.007.P01_13   6388   6370 -  cm    no 0.37  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid 
  (6) ?         3    7.7   0.0  ACBA.007.P01_13   2523   2502 -  cm    no 0.55  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid 
  (7) ?       6.6    6.4   0.0  ACBA.007.P01_13   2449   2398 -  cm    no 0.54  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid 


Hit alignments:
>> ACBA.007.P01_13  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid pMDR-ZJ06, complete
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       cyksc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ------ ----- ----
  (1) !   5.1e-10   46.4   0.0  cm        1       47 []        4059        4000 - ..   46.4    no 0.55

                                                                                           NC
                       <<<<<<<--------<<<-<<<<....................>>>>>>>---------->>>>>>> CS
           attC_4    1 GcCUAACAAgUCAUUGuUCAAGc....................gCUUAaCUCgGcCAUUCGUUAGgC 47  
                        :CUAACAA+UC   GUUCAAGC                    GCUUAACUC+++C    GUUAG: 
  ACBA.007.P01_13 4059 GUCUAACAAUUC---GUUCAAGCcgacgccgcuucgcggcgcgGCUUAACUCAAGC----GUUAGAU 4000
                       [Rsec=]========[=Lsec=]....................[Lprim]==========[Rprim] RF

>> ACBA.007.P01_13  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid pMDR-ZJ06, complete
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       cyksc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ------ ----- ----
  (2) !   5.5e-08   38.4   0.0  cm        1       47 []        5901        5793 - ..   38.4    no 0.59

                       v                                                                                            NC
                       <<<<<<<--------<<<-<<<<..................................................................... CS
           attC_4    1 GcCUAACAAgUCAUUGuUCAAGc..................................................................... 23  
                        CC AACAA+UC   GUUCAAGC                                                                     
  ACBA.007.P01_13 5901 ACCUAACAAUUC---GUUCAAGCcgagaucgcuucgcggccgcggaguuguucggaaaaauugucacaacgccgcggccgcaaagcgcuccg 5813
                       [Rsec=]========[=Lsec=]..................................................................... RF

                                              v NC
                       >>>>>>>---------->>>>>>> CS
           attC_4   24 gCUUAaCUCgGcCAUUCGUUAGgC 47  
                       GCUUAACUC+G+C    GUU GG 
  ACBA.007.P01_13 5812 GCUUAACUCAGGC----GUUGGGC 5793
                       [Lprim]==========[Rprim] RF

>> ACBA.007.P01_13  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid pMDR-ZJ06, complete
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       cyksc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ------ ----- ----
  (3) !   3.3e-05   27.3   0.0  cm        1       47 []        5324        5254 - ..   27.0    no 0.69

                          v                                                                     v     NC
                       <<<<<<<--------<<<-<<<<..............................>>>>>>>---------->>>>>.>> CS
           attC_4    1 GcCUAACAAgUCAUUGuUCAAGc..............................gCUUAaCUCgGcCAUUCGUUAG.gC 47  
                       GCC AACA G+C   G:UCAAGC                              GCUUA:CU GG+C    GUU G GC
  ACBA.007.P01_13 5324 GCCCAACAUGGC---GCUCAAGCcgaccggccagcccugcgggcuguccgucgGCUUAGCUAGGGC----GUUAGaGC 5254
                       [Rsec=]========[=Lsec=]..............................[Lprim]==========[Rpri.m] RF

>> ACBA.007.P01_13  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid pMDR-ZJ06, complete
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       cyksc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ------ ----- ----
  (4) ?     0.016   16.7   0.0  cm        1       47 []        4924        4856 - ..   16.0    no 0.61

                       v v                   v                             v                    v v NC
                       <<<<<<<--------<<<-<<<<.............................>>>>>>>---------->>>>>>> CS
           attC_4    1 GcCUAACAAgUCAUUGuUCAAGc.............................gCUUAaCUCgGcCAUUCGUUAGgC 47  
                        : UAAC A U    :UUC:AG                               CU:AA:  GG+C    GUUA : 
  ACBA.007.P01_13 4924 CUCUAACCAAUG---CUUCCAGCcgacgcuacggccuuacggccuccgcgcgUCUGAAGCUGGGC----GUUAAAC 4856
                       [Rsec=]========[=Lsec=].............................[Lprim]==========[Rprim] RF

>> ACBA.007.P01_13  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid pMDR-ZJ06, complete
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       cyksc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ------ ----- ----
  (5) ?       2.5    8.0   0.0  cm        1       47 []        6388        6370 - ..    6.1    no 0.37

                       v                v NC
                       <<<<<<~~~~~~>>>>>> CS
           attC_4    1 GcCUAA*[35]*UUAGgC 47  
                        CCUAA      UUAGG 
  ACBA.007.P01_13 6388 CCCUAA*[ 7]*UUAGGU 6370
                       [Rsec=~~~~~~Rprim] RF

>> ACBA.007.P01_13  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid pMDR-ZJ06, complete
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       cyksc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ------ ----- ----
  (6) ?         3    7.7   0.0  cm        1       47 []        2523        2502 - ..    7.0    no 0.55

                                            NC
                       <<<<<<<~~~~~~>>>>>>> CS
           attC_4    1 GcCUAAC*[33]*GUUAGgC 47  
                        CC:AAC      GUU:GG 
  ACBA.007.P01_13 2523 CCCCAAC*[ 8]*GUUGGGG 2502
                       [Rsec=]~~~~~~[Rprim] RF

>> ACBA.007.P01_13  08-JUN-2013 20301 bp Acinetobacter baumannii MDR-ZJ06 plasmid pMDR-ZJ06, complete
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       cyksc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ------ ----- ----
  (7) ?       6.6    6.4   0.0  cm        1       47 []        2449        2398 - ..    4.0    no 0.54

                          v                v    NC
                       <<<<<<<----~~~~~~>>>>>>> CS
           attC_4    1 GcCUAACAAgU*[29]*GUUAGgC 47  
                       GC: AAC  GU      GUU :GC
  ACBA.007.P01_13 2449 GCUCAACCGGU*[34]*GUUAAGC 2398
                       [Rsec=]====~~~~~~[Rprim] RF



Internal CM pipeline statistics summary:
----------------------------------------
Query model(s):                                                  1  (47 consensus positions)
Target sequences:                                                1  (7490 residues searched)
Target sequences re-searched for truncated hits:                 0  (0 residues re-searched)
Windows   passing  local HMM SSV           filter:                  (off)
Windows   passing  local HMM Viterbi       filter:                  (off)
Windows   passing  local HMM Viterbi  bias filter:                  (off)
Windows   passing  local HMM Forward       filter:                  (off)
Windows   passing  local HMM Forward  bias filter:                  (off)
Windows   passing glocal HMM Forward       filter:                  (off)
Windows   passing glocal HMM Forward  bias filter:                  (off)
Envelopes passing glocal HMM envelope defn filter:                  (off)
Envelopes passing  local CM  CYK           filter:                  (off)
Total CM hits reported:                                          7  (0.05367); includes 0 truncated hit(s)

# CPU time: 31.42u 0.01s 00:00:31.43 Elapsed: 00:00:31.43
//
[ok]
