#!/bin/bash

# This script is to call run_analysis_UPdistance.py

source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3
reuse Java-1.8
source activate slideseq_pipeline_env

submission=$0
manifest=$1
library=$2
lane=$3
scriptpath=$4
outputpath=$5
specpath=$6

echo ${submission}

python ${scriptpath}/run_analysis_UPdistance.py ${manifest} ${library} ${lane}

