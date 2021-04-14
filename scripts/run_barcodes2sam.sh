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

python ${scriptpath}/run_barcodes2sam.py ${manifest} "${commandStr}" ${lane} ${slice}

