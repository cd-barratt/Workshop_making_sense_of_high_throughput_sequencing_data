####### UFZ cluster job submit #########

#!/bin/bash                                                                                                                            
#$ -N dadi_model
#$ -M christopher_david.barratt@uni-leipzig.de
#$ -m beas
#$ -wd /public/barratt/work/Lflavomaculatus/Outputs/dadi/
#$ -S /bin/bash
#$ -o /public/barratt/work/Lflavomaculatus/job_logs/$JOB_NAME-$JOB_ID.log
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=25G
#$ -l h_rt=96:00:00 
 
#get slot information etc. (determine no. of slots with -pe above (can also be a range so that it takes the maximum available slots but will settle for less))
echo "Job name: $JOB_NAME"
echo "$NSLOTS slots have been dedicated to this job"
echo "On following node(s)"
cat $PE_HOSTFILE
echo "pe host file is $PE_HOSTFILE"

#load environment
module load python/2/7.13-2
source /public/barratt/virtual_env_RADseq/bin/activate

#load modules, start job
cd /public/barratt/work/Lflavomaculatus/Inputs/dadi/
python /public/barratt/work/Lflavomaculatus/Inputs/dadi/dadi_Run_Optimizations_split_symmig_all.py
