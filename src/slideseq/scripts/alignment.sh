#!/bin/bash
#$ -l os=RedHat7
#$ -l h_vmem=8g
#$ -l h_rt=8:0:0
#$ -pe smp 8
#$ -binding linear:8
#$ -terse
#$ -notify
#$ -P macosko_lab
#$ -j y
#$ -m beas

source /broad/software/scripts/useuse
reuse Anaconda3
reuse Java-1.8

if [ -z "$CONDA_ENV" ]
then
  conda activate ${CONDA_ENV}
else
  echo "Error: conda environment is not set"
  exit 1
fi

# task array id is the lane
