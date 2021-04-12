#!/bin/bash
#$ -l h_vmem=150G
#$ -l h_rt=80:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

# This script is to call run_barcodes2sam.py

source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3
reuse Java-1.8
source activate slideseq_pipeline_env

submission=$0
manifest=$1
commandStr=$2
lane=$3
slice=$4
scriptpath=$5
outputpath=$6
specpath=$7

echo ${submission}

# in case of exit, set all permissions 
function finish {
    if [ -d ${specpath} ]; then
echo "${specpath}" >> /broad/macosko/jlanglie/tmp/SLIDE_SEQ_GROUP/$(date +"%d-%m-%Y__%H_%M_%S")__$RANDOM
    fi
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


python ${scriptpath}/run_barcodes2sam.py ${manifest} "${commandStr}" ${lane} ${slice}

