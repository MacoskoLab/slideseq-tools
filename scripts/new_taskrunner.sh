#!/bin/bash

# This script is to call write_bijective_mapping.py

source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3
reuse Java-1.8
source activate slideseq_pipeline_env

scriptpath=$1
output_dir=$2

# -u so unbuffered and write to log
python -u ${scriptpath}/new_taskrunner.py $scriptpath $output_dir 
