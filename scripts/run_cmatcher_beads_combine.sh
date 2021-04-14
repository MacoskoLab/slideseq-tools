#!/bin/bash
#$ -l h_vmem=10G
#$ -l h_rt=40:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

# This script is to call run_cmatcher_beads_combine.py

source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3
reuse Java-1.8
source activate slideseq_pipeline_env

submission=$0
manifest=$1
library=$2
scriptpath=$3
outputpath=$4
specpath=$5

echo ${submission}

python ${scriptpath}/run_cmatcher_beads_combine.py ${manifest} ${library}

