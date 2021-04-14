#!/bin/bash

# This script is to call gen_readme.py

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

python ${scriptpath}/gen_readme.py ${manifest} ${library} ${locusfunclist}

