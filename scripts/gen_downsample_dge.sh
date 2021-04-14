#!/bin/bash
#$ -l h_vmem=37G
#$ -l h_rt=15:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

# This script is to call gen_downsample_dge.py

source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3
reuse Java-1.8
source activate slideseq_pipeline_env

submission=$0
manifest=$1
library=$2
scriptpath=$3
locusfunclist=$4
ratio=$5
outputpath=$6
specpath=$7

echo ${submission}

python ${scriptpath}/gen_downsample_dge.py ${manifest} ${library} ${locusfunclist} ${ratio}

