#!/bin/bash
#$ -l h_vmem=10G
#$ -l h_rt=30:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

# This script is to call generate_plot_downsampling.py

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
outputpath=$5
specpath=$6

echo ${submission}

python ${scriptpath}/generate_plot_downsampling.py ${manifest} ${library} ${locusfunclist}

