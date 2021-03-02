#!/usr/bin/python

# This script is to generate BijectiveMapping.mat

from __future__ import print_function

import sys
import os
import getopt
import csv

import argparse
import glob
import re
import time
from subprocess import call
from datetime import datetime

import numpy as np

from new_submit_to_taskrunner import call_to_taskrunner
import traceback

# Convert string to boolean
def str2bool(s):
    return s.lower() == "true"
	
	
# Write to log file
def write_log(log_file, flowcell_barcode, log_string):
    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as logfile:
        logfile.write(dt_string+" [Slide-seq Flowcell Alignment Workflow - "+flowcell_barcode+"]: "+log_string+"\n")
    logfile.close()
    

def main():
    if len(sys.argv) != 4:
        print("Please provide three arguments: manifest file, library ID and locus function!")
        sys.exit()
    
    manifest_file = sys.argv[1]
    library = sys.argv[2]
    locus_function_list = sys.argv[3]

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
    
    flowcell_directory = options['flowcell_directory']
    output_folder = options['output_folder']
    metadata_file = options['metadata_file']
    flowcell_barcode = options['flowcell_barcode']
    
    library_folder = options['library_folder'] if 'library_folder' in options else '{}/libraries'.format(output_folder)
    tmpdir = options['temp_folder'] if 'temp_folder' in options else '{}/tmp'.format(output_folder)
    dropseq_folder = options['dropseq_folder'] if 'dropseq_folder' in options else '/broad/macosko/bin/dropseq-tools'
    picard_folder = options['picard_folder'] if 'picard_folder' in options else '/broad/macosko/bin/dropseq-tools/3rdParty/picard'
    STAR_folder = options['STAR_folder'] if 'STAR_folder' in options else '/broad/macosko/bin/dropseq-tools/3rdParty/STAR-2.5.2a'
    scripts_folder = options['scripts_folder'] if 'scripts_folder' in options else '/broad/macosko/jilong/slideseq_pipeline/scripts'
    is_NovaSeq = str2bool(options['is_NovaSeq']) if 'is_NovaSeq' in options else False
    is_NovaSeq_S4 = str2bool(options['is_NovaSeq_S4']) if 'is_NovaSeq_S4' in options else False
    num_slice_NovaSeq = int(options['num_slice_NovaSeq']) if 'num_slice_NovaSeq' in options else 10
    num_slice_NovaSeq_S4 = int(options['num_slice_NovaSeq_S4']) if 'num_slice_NovaSeq_S4' in options else 40
    
    runinfo_file = '{}/RunInfo.xml'.format(flowcell_directory)
    log_file = '{}/logs/workflow.log'.format(output_folder)
    
    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    reference = ''
    run_barcodematching = False
    puckcaller_path = ''
    bead_type = ''
    sequence = 'AAGCAGTGGTATCAACGCAGAGTGAATGGG'
    base_quality = '10'
    min_transcripts_per_cell = '10'
    email_address = ''
    experiment_date = ''
    with open('{}/parsed_metadata.txt'.format(output_folder), 'r') as fin:
        reader = csv.reader(fin, delimiter='\t')
        rows = list(reader)
        row0 = rows[0]
        for i in range(1, len(rows)):
            row = rows[i]
            lanes.append(row[row0.index('lane')])
            if row[row0.index('lane')] not in lanes_unique:
                lanes_unique.append(row[row0.index('lane')])
            libraries.append(row[row0.index('library')])
            if row[row0.index('library')] not in libraries_unique:
                libraries_unique.append(row[row0.index('library')])
            barcodes.append(row[row0.index('sample_barcode')])
            bead_structures.append(row[row0.index('bead_structure')])
            if row[row0.index('library')] == library:
                reference = row[row0.index('reference')]
                sequence = row[row0.index('start_sequence')]
                base_quality = row[row0.index('base_quality')]
                min_transcripts_per_cell = row[row0.index('min_transcripts_per_cell')]
                email_address = row[row0.index('email')]
                run_barcodematching = str2bool(row[row0.index('run_barcodematching')])
                puckcaller_path = row[row0.index('puckcaller_path')]
                bead_type = row[row0.index('bead_type')]
                experiment_date = row[row0.index('date')]
    fin.close()
    
    reference_folder = reference[:reference.rfind('/')]
    referencePure = reference[reference.rfind('/') + 1:]
    if (referencePure.endswith('.gz')):
        referencePure = referencePure[:referencePure.rfind('.')]
    referencePure = referencePure[:referencePure.rfind('.')]
    genome_dir = '{}/STAR'.format(reference_folder)
    intervals = '{}/{}.genes.intervals'.format(reference_folder, referencePure)
    annotations_file = '{}/{}.gtf'.format(reference_folder, referencePure)
    ref_flat = '{}/{}.refFlat'.format(reference_folder, referencePure)
    ribosomal_intervals = '{}/{}.rRNA.intervals'.format(reference_folder, referencePure)
    reference2 = referencePure + '.' + locus_function_list
    
    folder_running = '{}/status/running.write_bijective_mapping_{}_{}'.format(output_folder, library, locus_function_list)
    folder_finished = '{}/status/finished.write_bijective_mapping_{}_{}'.format(output_folder, library, locus_function_list)
    folder_failed = '{}/status/failed.write_bijective_mapping_{}_{}'.format(output_folder, library, locus_function_list)
    
    alignment_folder = '{}/{}_{}/{}/alignment'.format(library_folder, experiment_date, library, reference2)
    barcode_matching_folder = '{}/{}_{}/{}/barcode_matching'.format(library_folder, experiment_date, library, reference2)
    dge_gzfile = '{}/{}.digital_expression.txt.gz'.format(alignment_folder, library)
    dge_file = '{}/{}.digital_expression2.txt'.format(alignment_folder, library)
    uniqueMappedDge_file = '{}/{}.UniqueMappedDge.txt'.format(alignment_folder, library)
    MappedDGEForR_file = '{}/MappedDGEForR.csv'.format(alignment_folder)
    
    call(['mkdir', '-p', folder_running])

    try:
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
        # UniqueMappedIlluminaBarcodes
        bci_file = '{}/{}_barcode_matching.txt'.format(barcode_matching_folder, library)
        unique_bci_file = '{}/{}_unique_matched_illumina_barcodes.txt'.format(barcode_matching_folder, library)
        
        if not os.path.isfile(dge_file):
            os.system('gunzip -c {} > {}'.format(dge_gzfile, dge_file))
        
        location_file = '{}/{}_matched_bead_locations.txt'.format(barcode_matching_folder, library)       
        genename_file = '{}/{}_genenames.txt'.format(barcode_matching_folder, library)
        bcb_file = '{}/{}_unique_matched_beads.txt'.format(barcode_matching_folder, library)
        commandStr = 'perl '+scripts_folder+'/get_unique_mapped_dge.pl '+dge_file+' '+uniqueMappedDge_file+' '+genename_file+' '+bcb_file
        os.system(commandStr)
        
        # Call run_WriteBijectiveMapping
        # 'UniqueMappedBeads','UniqueMappedDGE','UniqueMappedIlluminaBarcodes','GeneNames'
        output_file = '{}/logs/run_WriteBijectiveMapping_{}_{}.log'.format(output_folder, library, locus_function_list)
        submission_script = '/broad/macosko/jilong/slideseq_pipeline/scripts/run_WriteBijectiveMapping.sh'
        call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=65G', '-notify', '-l', 'h_rt=10:0:0', '-j', 'y', '-P', 'macosko_lab', '-l', 'os=RedHat7', submission_script, '/broad/software/nonfree/Linux/redhat_7_x86_64/pkgs/matlab_2019a', scripts_folder, bcb_file, uniqueMappedDge_file, unique_bci_file, genename_file, location_file, puckcaller_path, output_folder]
        call_to_taskrunner(output_folder, call_args)
        
        commandStr = 'perl '+scripts_folder+'/txt2csv.pl '+dge_file+' '+MappedDGEForR_file
        os.system(commandStr)
        
        if os.path.isfile(dge_file):
            call(['rm', dge_file])
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
        call(['mv', folder_running, folder_finished])
    except Exception as exp:
        print("EXCEPTION:!")
        print(exp)
        traceback.print_tb(exp.__traceback__, file=sys.stdout)
        if os.path.isdir(folder_running):
            call(['mv', folder_running, folder_failed])
        elif os.path.isdir(folder_waiting):
            call(['mv', folder_waiting, folder_failed])
        else:
            call(['mkdir', '-p', folder_failed])
            
        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = "The Slide-seq workflow for "+library+" "+locus_function_list+" failed at the step of generating BijectiveMapping.mat. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)

        sys.exit()
    

if __name__ == "__main__":
    main()


