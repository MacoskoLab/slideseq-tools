#!/bin/bash
#$ -l h_vmem=50G
#$ -l h_rt=10:0:0
#$ -l os=RedHat7
#$ -notify
#$ -P macosko_lab
#$ -j y

# This script is to build genome reference

source /broad/software/scripts/useuse
reuse Picard-Tools
reuse STAR

submission=$0

dropseq_folder=$2

bgzipLocation=$4
outdir=$5
reference_fasta=$6
gtf=$7
name=$8
filtered_gene_biotypes=$9
mt_sequence="${10}"

echo ${submission}

output_fasta=${outdir}/${name}.fasta
sequence_dictionary=${outdir}/${name}.dict
output_gtf=${outdir}/${name}.gtf
reduced_gtf=${outdir}/${name}.reduced.gtf

java -jar ${PICARD} NormalizeFasta INPUT=${reference_fasta} OUTPUT=$output_fasta

java -jar ${PICARD} CreateSequenceDictionary REFERENCE=$output_fasta OUTPUT=$sequence_dictionary SPECIES=${name}

${dropseq_folder}/FilterGtf GTF=${gtf} SEQUENCE_DICTIONARY=$sequence_dictionary OUTPUT=$output_gtf ${filtered_gene_biotypes}

${dropseq_folder}/ConvertToRefFlat ANNOTATIONS_FILE=$output_gtf SEQUENCE_DICTIONARY=$sequence_dictionary OUTPUT=${outdir}/${name}.refFlat

${dropseq_folder}/ReduceGtf GTF=$output_gtf SEQUENCE_DICTIONARY=$sequence_dictionary OUTPUT=$reduced_gtf

if [ "${mt_sequence}" == "" ]; then
	${dropseq_folder}/CreateIntervalsFiles SEQUENCE_DICTIONARY=$sequence_dictionary REDUCED_GTF=$reduced_gtf PREFIX=${name} OUTPUT=${outdir}
else
	${dropseq_folder}/CreateIntervalsFiles SEQUENCE_DICTIONARY=$sequence_dictionary REDUCED_GTF=$reduced_gtf PREFIX=${name} OUTPUT=${outdir} MT_SEQUENCE=${mt_sequence}
fi

STAR --runMode genomeGenerate --genomeDir ${outdir}/STAR --genomeFastaFiles $output_fasta --sjdbGTFfile $output_gtf --sjdbOverhang 97

${bgzipLocation} $output_fasta

samtools faidx $output_fasta.gz

gunzip $output_fasta.gz

