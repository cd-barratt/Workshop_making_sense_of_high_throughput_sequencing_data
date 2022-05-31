####### UFZ cluster job submit #########

#!/bin/bash                                                                                                                            
#$ -N Admixture
#$ -M christopher_david.barratt@uni-leipzig.de
#$ -m beas
#$ -wd /public/barratt/
#$ -S /bin/bash
#$ -o /public/barratt/work/Lflavomaculatus/job_logs/$JOB_NAME-$JOB_ID.log
#$ -j y
#$ -pe smp -1
#$ -l h_vmem=3G
#$ -l highmem
#$ -l h_rt=1:00:00 
 
#get slot information etc. (determine no. of slots with -pe above (can also be a range so that it takes the maximum available slots but will settle for less))
echo "Job name: $JOB_NAME"
echo "$NSLOTS slots have been dedicated to this job"
echo "On following node(s)"
cat $PE_HOSTFILE
echo "pe host file is $PE_HOSTFILE"

#load environment
module load miniconda/3
source activate conda_env_RADseq
source /public/barratt/virtual_env_RAD_seq/bin/activate

#load modules, start job
module load admixture/1.3.0-1

# Lflavomaculatus
# redo conversion using plink just to be sure...
plink --file /public/barratt/work/Lflavomaculatus/Inputs/Admixture/Lflav --make-bed --recode --out /public/barratt/work/Lflavomaculatus/Outputs/Admixture/Lflav
#run admixture
cd /public/barratt/work/Lflavomaculatus/Outputs/Admixture/
bash
for K in 1 2 3 4 5 6 7 8 9 10; \
do admixture --cv=10 ./Lflav.bed $K| tee log_Lflav.${K}.out; done
grep -h CV log_Lflav*.out > ./RESULTS_Lflavomaculatus_admixture.txt
