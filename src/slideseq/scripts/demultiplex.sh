#!/bin/bash
#$ -N demux
#$ -l os=RedHat7
#$ -l h_vmem=16G
#$ -l h_rt=12:0:0
#$ -pe smp 8
#$ -binding linear:8
#$ -terse
#$ -notify
#$ -P macosko_lab
#$ -R y
#$ -j y
#$ -m eas

source /broad/software/scripts/useuse
reuse Anaconda3
reuse Java-1.8

set -e

# check that the data directory has all the needed files
java -Djava.io.tmpdir=${TMP_DIR} -XX:+UseParallelGC \
  -XX:GCTimeLimit=20 -XX:GCHeapFreeLimit=10 -Xms8g -Xmx124g \
  -jar ${PICARD_JAR} CheckIlluminaDirectory \
  --TMP_DIR ${TMP_DIR} \
  --LANES ${SGE_TASK_ID} \
  --BASECALLS_DIR ${BASECALLS_DIR} \
  --READ_STRUCTURE ${READ_STRUCTURE}

# extract barcodes
java -Djava.io.tmpdir=${TMP_DIR} -XX:+UseParallelGC \
  -XX:GCTimeLimit=20 -XX:GCHeapFreeLimit=10 -Xms64g -Xmx124g \
  -jar ${PICARD_JAR} ExtractIlluminaBarcodes \
  --TMP_DIR ${TMP_DIR} \
  --LANE ${SGE_TASK_ID} \
  --BASECALLS_DIR ${BASECALLS_DIR} \
  --READ_STRUCTURE ${READ_STRUCTURE} \
  --OUTPUT_DIR ${OUTPUT_DIR}/${FLOWCELL}/L00${SGE_TASK_ID}/barcodes \
  --BARCODE_FILE ${OUTPUT_DIR}/${FLOWCELL}/L00${SGE_TASK_ID}/barcode_params.txt \
  --METRICS_FILE ${OUTPUT_DIR}/${FLOWCELL}/L00${SGE_TASK_ID}/${FLOWCELL}.barcode_metrics.txt \
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
  --RUN_BARCODE ${FLOWCELL} \
  --BARCODES_DIR ${OUTPUT_DIR}/${FLOWCELL}/L00${SGE_TASK_ID}/barcodes \
  --LIBRARY_PARAMS ${OUTPUT_DIR}/${FLOWCELL}/L00${SGE_TASK_ID}/library_params.txt \
  --INCLUDE_NON_PF_READS false \
  --APPLY_EAMSS_FILTER false \
  --ADAPTERS_TO_CHECK null \
  --IGNORE_UNEXPECTED_BARCODES true \
  --SEQUENCING_CENTER BI

# give group access to all files and directories
chmod -R --silent g+rwX ${OUTPUT_DIR}
