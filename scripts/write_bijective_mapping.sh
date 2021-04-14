#!/bin/bash
#$ -l h_vmem=45G
#$ -l h_rt=2:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

# This script is to call write_bijective_mapping.py

source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3
reuse Java-1.8
source activate slideseq_pipeline_env

submission=$0
manifest=$1
library=$2
scriptpath=$3
locus_function=$4
outputpath=$5
specpath=$6

echo ${submission}

python ${scriptpath}/write_bijective_mapping.py ${manifest} ${library} ${locus_function}

