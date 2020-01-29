# slideseq-pipeline
A pipeline for processing Slide-seq data

# Requirement

The pipeline needs several public tools pre-installed:
1) Drop-seq tools: https://github.com/broadinstitute/Drop-seq
2) Picard: https://broadinstitute.github.io/picard/
3) STAR: https://github.com/alexdobin/STAR
4) Java
5) Samtools
6) gcc/g++
7) Python (prefer 3.6 or above)

Several Python packages need to be installed for calculation and ploting:
1) numpy
2) pandas
3) plotnine
4) matplotlib
    
# Build genome reference

The Slide-seq pipeline needs a genome reference in specific format for alignment and analysis. You could build a genome reference based on input fasta and gtf files. 

Command: build_reference.py manifest_file

Notice:
1) build_reference.py calls build_reference.sh
2) check example.buildreference.txt for manifest_file format

# Run the pipeline

Add below commands into run.sh and build_reference.sh or your bashrc file (more or less based on your machine):
1) use Java-1.8
2) use .samtools-1.7
3) use Python-3.6

Compile CMatcher (command might be different on your system): g++ -std=c++11 -o cmatcher cmatcher.cpp

Submit a request to the pipeline: python run_pipeline.py manifest_file

Notice: 
1) an email from slideseq@broadinstitute.org will be sent to you if email_address is specified in the manifest file when the submission is received, the workflow finishes, and/or any job fails.
2) in order to speed up the process of NovaSeq data and NovaSeq S4 data, the pipeline splits each lane into a few slices, and runs the alignment steps on the slices parallelly and combines the alignment outputs together. 

# Manifest

flowcell_directory: specify the directory of the Illumina BCL files that will be used as input

output_folder: specify the directory where the output files will be saved. Its parent directory must exist

library_folder: specify the directory where the alignment and barcode matching outputs are saved. Default value: output_folder/libraries

dropseq_folder: specify the directory of the Drop-seq tools

picard_folder: specify the directory of the Picard tool

STAR_folder: specify the directory of the STAR aligner

scripts_folder: specify the directory of the pipeline scripts

temp_folder: specify the directory where the temporary files will be saved. Default value: output_folder/tmp

flowcell_barcode: specify the flowcell barcode that will be used as part of the output folder and file names

metadata_file: specify the meta data file containing information including library, experiment_date, lane, sample_barcode, bead_structure, estimated_num_cells, estimated_num_beads, reference, locus_function_list, start_sequence, base_quality, min_transcripts_per_cell, run_barcodematching, bead_barcode_file (see example.metadata.txt)

option_file: specify the options for calling the Drop-seq tools (see options.txt)

illumina_platform: specify the Illumina platform. It could be MiniSeq, NextSeq, NovaSeq or NovaSeqS4. Default value is NextSeq

email_address: specify email address(es) for receiving message from the pipeline

# Metadata

library: specify the library ID

experiment_date: specify the experiment date that will be used as part of the output folder and file names

lane: specify the lane id for the library. It could be a single lane number like "1", a combination lane numbers like "1,2" or "{LANE}" for all of the lanes in the experiment. Default value is {LANE}

sample_barcode: specify sample barcode such as TAAGGCGA

bead_structure: specify bead structure such as 8C18X7C8M1X|*T

estimated_num_cells: specify the estimated number of cells

estimated_num_beads: specify the estimated number of beads

reference: specify the reference build to use for alignment

start_sequence: specify the starting sequence that is used to trim reads. Default value is AAGCAGTGGTATCAACGCAGAGTGAATGGG

base_quality: specify the minimum quality of reads that will be kept during the alignment process. Default value is 10

min_transcripts_per_cell: specify the minimum number of transcripts per cell that is used to filter barcodes during generating digital expression. Default value is 10

locus_function_list: specify a list of functional annotations that reads need to be completely contained by to be considered for analysis. This option is used to generate digital expression files. Possible values include exonic, exonic+intronic and intronic. Default value is exonic+intronic

run_barcodematching: specify whether to run barcode matching. Default value is False

bead_barcode_file: specify bead barcode file that contains a list of bead barcodes and related x coordinate and y coordinate as three columns. (See example.beadbarcodes.txt)

