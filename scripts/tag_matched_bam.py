#!/usr/bin/python

# This script is to tag bam using matched bead barcodes

from __future__ import print_function

import sys
import os
import getopt
import csv

import argparse
import glob
import re
from subprocess import call
from datetime import datetime


# Convert string to boolean
def str2bool(s):
    return s.lower() == "true"


# Write to log file
def write_log(log_file, flowcell_barcode, log_string):
    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as logfile:
        logfile.write('{} [Slide-seq Flowcell Alignment Workflow - {}]: {}\n'.format(dt_string, flowcell_barcode, log_string))
    logfile.close()


def main():
    if len(sys.argv) != 6:
        print("Please provide five arguments: manifest file, library ID, lane ID, slice ID and locus function list!")
        sys.exit()
    
    manifest_file = sys.argv[1]
    library = sys.argv[2]
    lane = sys.argv[3]
    slice = sys.argv[4]
    locus_function_list = sys.argv[5]

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
    scripts_folder = options['scripts_folder']
    output_folder = options['output_folder']
    metadata_file = options['metadata_file']
    option_file = options['option_file']
    flowcell_barcode = options['flowcell_barcode']
    
    library_folder = options['library_folder'] if 'library_folder' in options else '{}/libraries'.format(output_folder)
    tmpdir = options['temp_folder'] if 'temp_folder' in options else '{}/tmp'.format(output_folder)
    illumina_platform = options['illumina_platform'] if 'illumina_platform' in options else 'NextSeq'
    email_address = options['email_address'] if 'email_address' in options else ''
    
    # Read info from metadata file
    barcode = ''
    reference = ''
    experiment_date = ''
    with open('{}/parsed_metadata.txt'.format(output_folder), 'r') as fin:
        reader = csv.reader(fin, delimiter='\t')
        rows = list(reader)
        row0 = rows[0]
        for i in range(1, len(rows)):
            row = rows[i]
            if row[row0.index('library')] == library:
                barcode = row[row0.index('sample_barcode')]
                reference = row[row0.index('reference')]
                experiment_date = row[row0.index('experiment_date')]
                break
    fin.close()
    
    referencePure = reference[reference.rfind('/') + 1:]
    referencePure = referencePure[:referencePure.rfind('.')]
    reference2 = referencePure + '.' + locus_function_list

    log_file = '{}/logs/workflow.log'.format(output_folder)
    
    barcode_matching_folder = '{}/{}_{}/{}/barcode_matching'.format(library_folder, experiment_date, library, reference2)
    prefix_libraries = '{}/{}.{}.{}.{}'.format(barcode_matching_folder, flowcell_barcode, lane, slice, library)
    if (barcode):
        prefix_libraries += '.'+barcode
    mapped_bam = '{}.star_gene_exon_tagged2.bam'.format(prefix_libraries)
    mapped_sam = '{}/{}_{}_{}_aligned.sam'.format(barcode_matching_folder, library, lane, slice)
    tagged_sam = '{}/{}_{}_{}_tagged.sam'.format(barcode_matching_folder, library, lane, slice)
    tagged_bam = '{}/{}_{}_{}_tagged.bam'.format(barcode_matching_folder, library, lane, slice)
    combined_cmatcher_file = '{}/{}_barcode_matching.txt'.format(barcode_matching_folder, library)
    
    if not os.path.isfile(mapped_bam):
        write_log(log_file, flowcell_barcode, 'TagMatchedBam error: '+mapped_bam+' does not exist!')
        raise Exception('TagMatchedBam error: '+mapped_bam+' does not exist!')
    
    if not os.path.isfile(combined_cmatcher_file):
        write_log(log_file, flowcell_barcode, 'TagMatchedBam error: '+combined_cmatcher_file+' does not exist!')
        raise Exception('TagMatchedBam error: '+combined_cmatcher_file+' does not exist!')
        
    folder_running = '{}/status/running.tag_matched_bam_{}_{}_{}_{}'.format(output_folder, library, lane, slice, reference2)
    folder_finished = '{}/status/finished.tag_matched_bam_{}_{}_{}_{}'.format(output_folder, library, lane, slice, reference2)
    folder_failed = '{}/status/failed.tag_matched_bam_{}_{}_{}_{}'.format(output_folder, library, lane, slice, reference2)

    try:
        call(['mkdir', folder_running])
        
        write_log(log_file, flowcell_barcode, "Tag matched bam for "+library+" "+reference2+" in lane "+lane+" slice "+slice)
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
    
        commandStr = 'samtools view -h -o {} {}'.format(mapped_sam, mapped_bam)
        os.system(commandStr)
        
        call(['rm', mapped_bam])

        dict1 = {}
        with open(combined_cmatcher_file, 'r') as fin:
            j = 0
            for line in fin:
                j += 1
                if j > 1:
                    dict1[line.split('\t')[0]] = line.split('\t')[2]
        fin.close()

        with open(tagged_sam, 'w') as fout:
            with open(mapped_sam, 'r') as fin:
                for line in fin:
                    if line[0] == '@':
                        fout.write(line)
                    else:
                        items1 = line.split('\t')
                        bc1 = items1[11]
                        items2 = bc1.split(':')
                        bc2 = items2[2]
                        if bc2 in dict1:
                            items2[2] = dict1[bc2]
                            items1[11] = ':'.join(items2)
                            fout.write('\t'.join(items1))
            fin.close()
        fout.close()
        
        commandStr = 'samtools view -S -b {} > {}'.format(tagged_sam, tagged_bam)
        os.system(commandStr)
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
        if os.path.isfile(mapped_sam):
            call(['rm', mapped_sam])
        if os.path.isfile(tagged_sam):
            call(['rm', tagged_sam])
            
        write_log(log_file, flowcell_barcode, "Tag matched bam for "+library+" "+reference2+" in lane "+lane+" slice "+slice+" is done. ")
        
        call(['mv', folder_running, folder_finished])
    except:
        if os.path.isdir(folder_running):
            call(['mv', folder_running, folder_failed])
        else:
            call(['mkdir', folder_failed])
            
        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = "The Slide-seq workflow for "+library+" "+reference2+" in lane "+lane+" slice "+slice+" failed at the step of tagging matched bam. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)
        
        sys.exit()


if __name__ == "__main__":
    main()


