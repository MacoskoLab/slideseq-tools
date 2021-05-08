#!/bin/bash
#$ -N build_reference
#$ -l os=RedHat7
#$ -l h_vmem=8G
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

set -e

if [ -z "${CONDA_ENV}" ]; then
  echo "Error: conda environment is not set"
  exit 1
else
  source activate ${CONDA_ENV}
fi

OUTPUT_FASTA=${OUTPUT_DIR}/${GENOME_NAME}.fasta
SEQUENCE_DICTIONARY=${OUTPUT_DIR}/${GENOME_NAME}.dict
OUTPUT_GTF=${OUTPUT_DIR}/${GENOME_NAME}.gtf
REDUCED_GTF=${OUTPUT_DIR}/${GENOME_NAME}.reduced.gtf

java -jar ${PICARD_JAR} NormalizeFasta \
  --INPUT ${REFERENCE_FASTA} \
  --OUTPUT ${OUTPUT_FASTA}

java -jar ${PICARD_JAR} CreateSequenceDictionary \
  --REFERENCE ${OUTPUT_FASTA} \
  --OUTPUT ${SEQUENCE_DICTIONARY} \
  --SPECIES ${GENOME_NAME}

${DROPSEQ_DIR}/FilterGtf \
  GTF=${REFERENCE_GTF} \
  SEQUENCE_DICTIONARY=${SEQUENCE_DICTIONARY} \
  OUTPUT=${OUTPUT_GTF} \
  ${FILTERED_BIOTYPES}

${DROPSEQ_DIR}/ConvertToRefFlat \
  ANNOTATIONS_FILE=${OUTPUT_GTF} \
  SEQUENCE_DICTIONARY=${SEQUENCE_DICTIONARY} \
  OUTPUT=${OUTPUT_DIR}/${GENOME_NAME}.refFlat

${DROPSEQ_DIR}/ReduceGtf \
  GTF=${OUTPUT_GTF} \
  SEQUENCE_DICTIONARY=${SEQUENCE_DICTIONARY} \
  OUTPUT=${REDUCED_GTF}

if [ -z "${MT_SEQUENCE}" ]; then
  ${DROPSEQ_DIR}/CreateIntervalsFiles \
    SEQUENCE_DICTIONARY=${SEQUENCE_DICTIONARY} \
    REDUCED_GTF=${REDUCED_GTF} \
    PREFIX=${GENOME_NAME} \
    OUTPUT=${OUTPUT_DIR}
else
  ${DROPSEQ_DIR}/CreateIntervalsFiles \
    SEQUENCE_DICTIONARY=${SEQUENCE_DICTIONARY} \
    REDUCED_GTF=${REDUCED_GTF} \
    PREFIX=${GENOME_NAME} \
    OUTPUT=${OUTPUT_DIR} \
    MT_SEQUENCE=${MT_SEQUENCE}
fi

STAR --runMode genomeGenerate \
  --genomeDir ${OUTPUT_DIR}/STAR \
  --outFilenamePrefix ${OUTPUT_DIR}/ \
  --genomeFastaFiles ${OUTPUT_FASTA} \
  --sjdbGTFfile ${OUTPUT_GTF} \
  --sjdbOverhang 97 \
  --runThreadN 8

bgzip ${OUTPUT_FASTA}

samtools faidx ${OUTPUT_FASTA}.gz

gunzip ${OUTPUT_FASTA}.gz
