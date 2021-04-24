#!/bin/bash
#$ -l os=RedHat7
#$ -l h_vmem=8g
#$ -l h_rt=1:0:0
#$ -pe smp 4
#$ -binding linear:4
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
  echo "Error: conda environment is not set"
  exit 1
else
  conda activate ${CONDA_ENV}
fi

# do whatever post-processing is needed
