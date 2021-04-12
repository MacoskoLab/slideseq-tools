#!/bin/bash
#$ -l h_vmem=10g
#$ -l h_rt=100:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

# This script is to call run_mergebarcodes.py

source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3
reuse Java-1.8
source activate slideseq_pipeline_env

submission=$0
manifest=$1
scriptpath=$2
outputpath=$3

echo ${submission}

# in case of exit, set all permissions 
function finish {
    log_lib=${outputpath}/logs
    if [ -d "$log_lib" ]; then
echo "$log_lib" >> /broad/macosko/jlanglie/tmp/SLIDE_SEQ_GROUP/$(date +"%d-%m-%Y__%H_%M_%S")__$RANDOM
    fi
	status_lib=${outputpath}/status
    if [ -d "$status_lib" ]; then
echo "$status_lib" >> /broad/macosko/jlanglie/tmp/SLIDE_SEQ_GROUP/$(date +"%d-%m-%Y__%H_%M_%S")__$RANDOM
    fi
	tmp_lib=${outputpath}/tmp
    if [ -d "$tmp_lib" ]; then
echo "$tmp_lib" >> /broad/macosko/jlanglie/tmp/SLIDE_SEQ_GROUP/$(date +"%d-%m-%Y__%H_%M_%S")__$RANDOM
    fi
}
trap finish SIGUSR2 EXIT


python ${scriptpath}/run_mergebarcodes.py ${manifest}

