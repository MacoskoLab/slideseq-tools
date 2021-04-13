#!/bin/bash
#$ -l h_vmem=27G
#$ -l h_rt=20:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

# This script is to call CMatcher

source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3
source activate slideseq_pipeline_env

submission=$0
scriptpath=$1
bead_barcode_file=$2
select_cell_file=$3
output_distance_file=$4
output_detail_file=$5
bead_type=$6
outputpath=$7
specpath=$8

echo ${submission}

# in case of exit, set all permissions 
function finish {
    if [ -d ${specpath} ]; then
echo "${specpath}" >> /broad/macosko/jlanglie/tmp/SLIDE_SEQ_GROUP/$(date +"%d-%m-%Y__%H_%M_%S")__$RANDOM
    fi
    log_lib=${outputpath}/logs
    if [ -d "$log_lib" ]; then
echo "$log_lib" >> /broad/macosko/jlanglie/tmp/SLIDE_SEQ_GROUP/$(date +"%d-%m-%Y__%H_%M_%S")__$RANDOM
    fi

}
trap finish SIGUSR2 EXIT

${scriptpath}/cmatcher ${bead_barcode_file} ${select_cell_file} ${output_distance_file} ${output_detail_file} ${bead_type}

