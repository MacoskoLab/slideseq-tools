#!/bin/bash

# This script is to call run_alignment.py

source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3
reuse Java-1.8
reuse .samtools-1.7
source activate slideseq_pipeline_env

submission=$0
manifest=$1
library=$2
lane=$3
slice=$4
barcode=$5
scriptpath=$6
outputpath=$7
specpath=$8

echo ${submission}

# in case of exit, set all permissions 
function finish {
    if [ -d ${specpath} ]; then
echo "${specpath}" >> /broad/macosko/jlanglie/tmp/SLIDE_SEQ_GROUP/$(date +"%d-%m-%Y__%H_%M_%S")__$RANDOM
    fi
    if [ -d ${outputpath} ]; then
echo "${outputpath}" >> /broad/macosko/jlanglie/tmp/SLIDE_SEQ_GROUP/$(date +"%d-%m-%Y__%H_%M_%S")__$RANDOM
    fi

}
trap finish SIGUSR2 EXIT

python ${scriptpath}/run_alignment.py ${manifest} ${library} ${lane} ${slice} ${barcode}

