#!/bin/bash
#$ -N downsampling
#$ -l os=RedHat7
#$ -l h_vmem=16G
#$ -l h_rt=96:0:0
#$ -pe smp 2
#$ -binding linear:2
#$ -terse
#$ -notify
#$ -R y
#$ -j y
#$ -m eas

source /broad/software/scripts/useuse
reuse Anaconda3
reuse Java-1.8
reuse Google-Cloud-SDK

set -e

if [ -z "${CONDA_ENV}" ]; then
  echo "Error: conda environment is not set"
  exit 1
else
  source activate ${CONDA_ENV}
fi

# do whatever post-processing is needed
downsample_library ${DEBUG} \
  --library-index ${SGE_TASK_ID} \
  --manifest-file ${MANIFEST}
