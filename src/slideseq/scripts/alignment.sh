#!/bin/bash
#$ -l os=RedHat7
#$ -l h_vmem=8g
#$ -l h_rt=8:0:0
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

if [ -z "${CONDA_ENV}" ]; then
  echo "Error: conda environment is not set"
  exit 1
else
  source activate ${CONDA_ENV}
fi

# this will run slideseq.pipeline.alignment:main
align_sample ${DEBUG} \
  --lane ${LANE} \
  --sample-index ${SGE_TASK_ID} \
  --manifest-file ${MANIFEST}
