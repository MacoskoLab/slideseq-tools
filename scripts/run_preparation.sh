#!/bin/bash

# This script is to call run_preparation.py

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
    if [ -d ${outputpath} ]; then
echo "${outputpath}" >> /broad/macosko/jlanglie/tmp/SLIDE_SEQ_GROUP/$(date +"%d-%m-%Y__%H_%M_%S")__$RANDOM
    fi

}
trap finish SIGUSR2 EXIT

python ${scriptpath}/run_preparation.py ${manifest}

