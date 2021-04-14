#!/bin/bash

# This script is to call tag_matched_bam.py

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
locus_function_list=$6
scriptpath=$7
outputpath=$8
specpath=$9

echo ${submission}

python ${scriptpath}/tag_matched_bam.py ${manifest} ${library} ${lane} ${slice} ${barcode} ${locus_function_list}

