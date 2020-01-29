from __future__ import print_function

# This script is to build genome reference for the Slide-seq pipeline

import sys
import os
import getopt

import argparse
import glob
import re
from subprocess import call


if len(sys.argv) != 2:
    print("Please provide one argument: manifest file!")
    sys.exit()

manifest_file = sys.argv[1]

# Check if the manifest file exists
if not os.path.isfile(manifest_file):
    print("File {} does not exist. Exiting...".format(manifest_file))
    sys.exit()

# Read manifest file
options = {}
with open(manifest_file,"r") as fp:
    for line in fp:
        dict = line.rstrip().split("=")
        options[dict[0]] = dict[1]
fp.close()

reference_fasta = options['reference_fasta']
gtf = options['gtf']
name = options['name']
output_folder = options['output_folder']
dropseq_folder = options['dropseq_folder']
picard_folder = options['picard_folder']
STAR_folder = options['STAR_folder']
bgzip_location = options['bgzip_location']
submission_script = options['shell_script']
filtered_gene_biotypes = options['filtered_gene_biotypes'] if 'filtered_gene_biotypes' in options else 'processed_pseudogene,unprocessed_pseudogene,transcribed_unprocessed_pseudogene,pseudogene,IG_V_pseudogene,transcribed_processed_pseudogene,TR_J_pseudogene,TR_V_pseudogene,unitary_pseudogene,polymorphic_pseudogene,IG_D_pseudogene,translated_processed_pseudogene,translated_unprocessed_pseudogene,IG_C_pseudogene'

print("Creating Slide-seq genome reference...")

filtered_gene_biotypes2 = ''
if len(filtered_gene_biotypes) > 0:
    types = filtered_gene_biotypes.split(',')
    for type in types:
        filtered_gene_biotypes2 += 'G={} '.format(type)

if not os.path.isdir(output_folder):
    call(['mkdir', output_folder])
if not os.path.isdir('{}/STAR'.format(output_folder)):
    call(['mkdir', '{}/STAR'.format(output_folder)])

sequence_dictionary = '{}/{}.dict'.format(output_folder, name)
if os.path.isdir(sequence_dictionary):
    call(['rm', '-r', sequence_dictionary])

output_file = '{}/run.log'.format(output_folder)
call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=50g', '-l', 'h_rt=10:0:0', '-j', 'y', submission_script, 
             picard_folder, dropseq_folder, STAR_folder, bgzip_location, output_folder, reference_fasta, gtf, name, filtered_gene_biotypes2]
call(call_args)

