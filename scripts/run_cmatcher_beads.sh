#!/bin/bash
#$ -l h_vmem=10G
#$ -l h_rt=5:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

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

${scriptpath}/cmatcher_beads ${bead_barcode_file1} ${bead_barcode_file2} ${bead_location_file} ${output_file_01} ${output_file_2} ${bead_type}

