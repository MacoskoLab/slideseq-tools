#!/bin/bash
#$ -l h_vmem=62G
#$ -l h_rt=23:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

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

python ${scriptpath}/run_alignment.py ${manifest} ${library} ${lane} ${slice} ${barcode}

