#!/bin/bash

# This script is to call cmatcher_beads

source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3
source activate slideseq_pipeline_env

submission=$0
scriptpath=$1
bead_barcode_file1=$2
bead_barcode_file2=$3
bead_location_file=$4
output_file_01=$5
output_file_2=$6
bead_type=$7
outputpath=$8
specpath=$9

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

${scriptpath}/cmatcher_beads ${bead_barcode_file1} ${bead_barcode_file2} ${bead_location_file} ${output_file_01} ${output_file_2} ${bead_type}

