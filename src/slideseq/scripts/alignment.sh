#!/bin/bash
#$ -N alignment
#$ -l os=RedHat7
#$ -l h_vmem=8G
#$ -l h_rt=48:0:0
#$ -pe smp 8
#$ -binding linear:8
#$ -terse
#$ -notify
#$ -P macosko_lab
#$ -R y
#$ -j y
#$ -m beas

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

# this will run slideseq.pipeline.alignment:main
align_library ${DEBUG} \
  --flowcell ${FLOWCELL} \
  --lane ${LANE} \
  --library-index ${SGE_TASK_ID} \
  --manifest-file ${MANIFEST}
