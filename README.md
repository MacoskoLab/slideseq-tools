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

Notice: an email from slideseq@broadinstitute.org will be sent to you if email_address is specified in the manifest file when the submission is received, the workflow finishes, and/or any job fails.

