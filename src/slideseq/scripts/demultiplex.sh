#!/bin/bash
#$ -l os=RedHat7
#$ -l h_vmem=16g
#$ -l h_rt=4:0:0
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

# check that the data directory has all the needed files
java -Djava.io.tmpdir=${TMP_DIR} -XX:+UseParallelGC \
  -XX:GCTimeLimit=20 -XX:GCHeapFreeLimit=10 -Xms8g -Xmx124g \
  -jar ${PICARD_JAR} CheckIlluminaDirectory \
  --TMP_DIR ${TMP_DIR} \
  --LANE ${SGE_TASK_ID} \
  --BASECALLS_DIR ${BASECALLS_DIR} \
  --READ_STRUCTURE ${READ_STRUCTURE}


# extract barcodes to $BARCODE_FILE
java -Djava.io.tmpdir=${TMP_DIR} -XX:+UseParallelGC \
  -XX:GCTimeLimit=20 -XX:GCHeapFreeLimit=10 -Xms64g -Xmx124g \
  -jar ${PICARD_JAR} ExtractIlluminaBarcodes \
  --TMP_DIR ${TMP_DIR} \
  --LANE ${SGE_TASK_ID} \
  --BASECALLS_DIR ${BASECALLS_DIR} \
  --READ_STRUCTURE ${READ_STRUCTURE} \
  --OUTPUT_DIR ${OUTPUT_DIR}/L00${SGE_TASK_ID}/barcodes \
  --BARCODE_FILE ${BARCODE_FILE} \
  --METRICS_FILE ${METRICS_FILE} \
  --COMPRESS_OUTPUTS true \
  --NUM_PROCESSORS 8


# create uBAM files
java -Djava.io.tmpdir=${TMP_DIR} -XX:+UseParallelGC \
  -XX:GCTimeLimit=20 -XX:GCHeapFreeLimit=10 -Xms64g -Xmx124g \
  -jar ${PICARD_JAR} IlluminaBasecallsToSam \
  --TMP_DIR ${TMP_DIR} \
  --LANE ${SGE_TASK_ID} \
  --BASECALLS_DIR ${BASECALLS_DIR} \
  --READ_STRUCTURE ${READ_STRUCTURE} \
  --RUN_BARCODE ${RUN_BARCODE} \
  --BARCODES_DIR ${BARCODE_DIR} \
  --LIBRARY_PARAMS ${LIBRARY_PARAMS} \
  --INCLUDE_NON_PF_READS false \
  --APPLY_EAMSS_FILTER false \
  --ADAPTERS_TO_CHECK null \
  --IGNORE_UNEXPECTED_BARCODES true \
  --SEQUENCING_CENTER BI
