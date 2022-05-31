#!/bin/bash
#bed2diff:
/Users/chris/Desktop/Workshop_2/work/Lflavomaculatus/local_files/eems-master/bed2diffs/src-wout-openmp/bed2diffs_v1 --bfile /Users/chris/Desktop/Workshop_2/work/Lflavomaculatus/Inputs/EEMS/Lflav
/Users/chris/Desktop/Workshop_2/work/Lflavomaculatus/local_files/eems-master/bed2diffs/src-wout-openmp/bed2diffs_v2 --bfile /Users/chris/Desktop/Workshop_2/work/Lflavomaculatus/Inputs/EEMS/Lflav

#runeems_SNPs:
/Users/chris/Desktop/Workshop_2/work/Lflavomaculatus/local_files/eems-master/runeems_snps/src/runeems_snps --params /Users/chris/Desktop/Workshop_2/work/Lflavomaculatus/Inputs/EEMS/700-Lflav-chain1r.ini
/Users/chris/Desktop/Workshop_2/work/Lflavomaculatus/local_files/eems-master/runeems_snps/src/runeems_snps --params /Users/chris/Desktop/Workshop_2/work/Lflavomaculatus/Inputs/EEMS/700-Lflav-chain2r.ini
