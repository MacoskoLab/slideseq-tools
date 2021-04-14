#!/bin/bash
#$ -l h_vmem=70G
#$ -l h_rt=6:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

# This script is to call gen_sparse_matrix.py

source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3
reuse Java-1.8
source activate slideseq_pipeline_env

submission=$0
manifest=$1
library=$2
locus_function_list=$3
folder=$4
filename=$5
scriptpath=$6
outputpath=$7
specpath=$8

echo ${submission}

python ${scriptpath}/gen_sparse_matrix.py ${manifest} ${library} ${locus_function_list} ${folder} ${filename}

