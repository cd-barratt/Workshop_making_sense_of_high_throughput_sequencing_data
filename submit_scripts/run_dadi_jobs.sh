#!/bin/bash                                                                                                                            

cd /public/barratt/submit_scripts/dadi/

qsub dadi_Run_Optimizations_ancmig_1.sh
qsub dadi_Run_Optimizations_ancmig_2.sh
qsub dadi_Run_Optimizations_ancmig_3.sh
qsub dadi_Run_Optimizations_refugia_1.sh
qsub dadi_Run_Optimizations_refugia_2.sh
qsub dadi_Run_Optimizations_refugia_3.sh
qsub dadi_Run_Optimizations_sim_split_no_mig.sh
qsub dadi_Run_Optimizations_split_no_mig_human.sh
qsub dadi_Run_Optimizations_sim_split_refugia_sym_mig_adjacent.sh
qsub dadi_Run_Optimizations_sim_split_refugia_sym_mig_all.sh
qsub dadi_Run_Optimizations_sim_split_sym_mig_adjacent.sh
qsub dadi_Run_Optimizations_sim_split_sym_mig_all.sh
qsub dadi_Run_Optimizations_split_no_mig.sh
qsub dadi_Run_Optimizations_split_sym_mig_adjacent.sh
qsub dadi_Run_Optimizations_split_sym_mig_all.sh
