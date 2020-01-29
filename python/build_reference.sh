#!/bin/bash

# This script is to build genome reference for the Slide-seq pipeline

submission=$0
picard_folder=$1
dropseq_folder=$2
STAR_folder=$3
bgzipLocation=$4
outdir=$5
reference_fasta=$6
gtf=$7
name=$8
filtered_gene_biotypes=$9

echo ${submission}

output_fasta=${outdir}/${name}.fasta
sequence_dictionary=${outdir}/${name}.dict
output_gtf=${outdir}/${name}.gtf
reduced_gtf=${outdir}/${name}.reduced.gtf

java -jar ${picard_folder}/picard.jar NormalizeFasta INPUT=${reference_fasta} OUTPUT=$output_fasta

java -jar ${picard_folder}/picard.jar CreateSequenceDictionary REFERENCE=$output_fasta OUTPUT=$sequence_dictionary SPECIES=${name}

${dropseq_folder}/FilterGtf GTF=${gtf} SEQUENCE_DICTIONARY=$sequence_dictionary OUTPUT=$output_gtf ${filtered_gene_biotypes}

${dropseq_folder}/ConvertToRefFlat ANNOTATIONS_FILE=$output_gtf SEQUENCE_DICTIONARY=$sequence_dictionary OUTPUT=${outdir}/${name}.refFlat

${dropseq_folder}/ReduceGtf GTF=$output_gtf SEQUENCE_DICTIONARY=$sequence_dictionary OUTPUT=$reduced_gtf

${dropseq_folder}/CreateIntervalsFiles SEQUENCE_DICTIONARY=$sequence_dictionary REDUCED_GTF=$reduced_gtf PREFIX=${name} OUTPUT=${outdir} MT_SEQUENCE=chrM

${STAR_folder}/STAR --runMode genomeGenerate --genomeDir ${outdir}/STAR --genomeFastaFiles $output_fasta --sjdbGTFfile $output_gtf --sjdbOverhang 97

${bgzipLocation} $output_fasta

samtools faidx $output_fasta.gz

gunzip $output_fasta.gz

