####### UFZ cluster job submit #########

#!/bin/bash                                                                                                                            
#$ -N SNAPP
#$ -M christopher_david.barratt@uni-leipzig.de
#$ -m beas
#$ -wd /work/$USER
#$ -S /bin/bash
#$ -o /public/barratt/work/Lflavomaculatus/job_logs/$JOB_NAME-$JOB_ID.log
#$ -j y
#$ -pe smp 4
#$ -l h_vmem=10G
#$ -l h_rt=120:00:00 
 
#get slot information etc. (determine no. of slots with -pe above (can also be a range so that it takes the maximum available slots but will settle for less))
echo "Job name: $JOB_NAME"
echo "$NSLOTS slots have been dedicated to this job"
echo "On following node(s)"
cat $PE_HOSTFILE
echo "pe host file is $PE_HOSTFILE"

#load modules, start job
source ~/barratt_env/bin/activate
module load beast/2.4.6-1
module load java

cd /public/barratt/work/Lflavomaculatus/Outputs/SNAPP/
 
export JVM_ARGS="-Xms100M -Xmx25G"
beast -overwrite -beagle_GPU -threads ${NSLOTS:-1} /public/barratt/work/Lflavomaculatus/Inputs/SNAPP/Lflavomaculatus_thinned.xml
