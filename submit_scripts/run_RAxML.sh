####### UFZ cluster job submit #########

#!/bin/bash                                                                                                                            
#$ -N RAxML
#$ -M christopher_david.barratt@uni-leipzig.de
#$ -m beas
#$ -wd /public/barratt/work
#$ -S /bin/bash
#$ -o /public/barratt/work/Lflavomaculatus/job_logs/$JOB_NAME-$JOB_ID.log
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=72:00:00 
 
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
module load python/2/7.13-2

module load raxml/8.1.17-1

cd /public/barratt/work/Lflavomaculatus/Outputs/RAxML/

raxmlHPC -f a -m GTRCAT -p 12345 -x 12345 -# 100 -s /public/barratt/work/Lflavomaculatus/Inputs/RAxML/Lflavomaculatus.phy ­­asc­corr=lewis -n bootstrapped 
raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.bootstrapped -z RAxML_bootstrap.bootstrapped -n FINAL_bootstrapped.tre

