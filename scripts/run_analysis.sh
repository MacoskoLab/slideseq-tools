#!/bin/bash
#$ -l h_vmem=30G
#$ -l h_rt=100:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

# This script is to call run_analysis.py

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

python ${scriptpath}/run_analysis.py ${manifest} ${library}

