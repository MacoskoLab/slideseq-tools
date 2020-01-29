#!/usr/bin/python

# This script is to filter unmapped bam using matched barcodes

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


# Get bead structure range
def get_bead_structure_range(bs, type):
    #12C8M|*T
    #7C18X7C8M2X|*T
    l = re.split('C|X|M', re.split('\|', bs)[0])
    res = ''
    i = 1
    p = -1
    for it in l:
        if it:
            p += len(it) + 1
            if bs[p] == type:
                res += str(i) + '-' + str(i + int(it) - 1) + ':'
            i += int(it)
    return res[:-1]


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
    bead_structure = ''
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
                bead_structure = row[row0.index('bead_structure')]
                experiment_date = row[row0.index('experiment_date')]
                break
    fin.close()
    
    referencePure = reference[reference.rfind('/') + 1:]
    referencePure = referencePure[:referencePure.rfind('.')]
    reference2 = referencePure + '.' + locus_function_list

    log_file = '{}/logs/workflow.log'.format(output_folder)

    prefix_libraries = '{}/{}_{}/{}.{}.{}.{}'.format(library_folder, experiment_date, library, flowcell_barcode, lane, slice, library)
    if (barcode):
        prefix_libraries += '.'+barcode
    unmapped_bam = '{}.unmapped.bam'.format(prefix_libraries)
    
    barcode_matching_folder = '{}/{}_{}/{}/barcode_matching'.format(library_folder, experiment_date, library, reference2)
    unmapped_sam = '{}/{}_{}_{}_unmapped.sam'.format(barcode_matching_folder, library, lane, slice)
    filtered_sam = '{}/{}_{}_{}_filtered.sam'.format(barcode_matching_folder, library, lane, slice)
    filtered_bam = '{}/{}_{}_{}_filtered.bam'.format(barcode_matching_folder, library, lane, slice)
    combined_cmatcher_file = '{}/{}_barcode_matching.txt'.format(barcode_matching_folder, library)
    
    if not os.path.isfile(unmapped_bam):
        write_log(log_file, flowcell_barcode, 'TagMatchedBam error: '+unmapped_bam+' does not exist!')
        raise Exception('TagMatchedBam error: '+unmapped_bam+' does not exist!')
    
    if not os.path.isfile(combined_cmatcher_file):
        write_log(log_file, flowcell_barcode, 'TagMatchedBam error: '+combined_cmatcher_file+' does not exist!')
        raise Exception('TagMatchedBam error: '+combined_cmatcher_file+' does not exist!')
        
    folder_running = '{}/status/running.filter_unmapped_bam_{}_{}_{}_{}'.format(output_folder, library, lane, slice, reference2)
    folder_finished = '{}/status/finished.filter_unmapped_bam_{}_{}_{}_{}'.format(output_folder, library, lane, slice, reference2)
    folder_failed = '{}/status/failed.filter_unmapped_bam_{}_{}_{}_{}'.format(output_folder, library, lane, slice, reference2)

    try:
        call(['mkdir', folder_running])
        
        write_log(log_file, flowcell_barcode, "Filter unmapped bam for "+library+" "+reference2+" in lane "+lane+" slice "+slice)
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)

        commandStr = 'samtools view -h -o {} {}'.format(unmapped_sam, unmapped_bam)
        os.system(commandStr)

        dict1 = {}
        with open(combined_cmatcher_file, 'r') as fin:
            j = 0
            for line in fin:
                j += 1
                if j > 1:
                    dict1[line.split('\t')[0]] = line.split('\t')[2]
        fin.close()
        
        bs_range1 = get_bead_structure_range(bead_structure, 'C') # 1-7:26-32
        lists = re.split(':', bs_range1)
        with open(filtered_sam, 'w') as fout:
            with open(unmapped_sam, 'r') as fin:
                flag = False
                for line in fin:
                    if line[0] == '@':
                        fout.write(line)
                    else:
                        items1 = line.split('\t')
                        if items1[1] == '77':
                            read1 = items1[9]
                            r = ""
                            for s in lists:
                                v = re.split('-', s)
                                r += read1[int(v[0])-1:int(v[1])]
                            if r in dict1:
                                fout.write(line)
                                flag = True
                            else:
                                flag = False
                        elif items1[1] == '141':
                            if flag:
                                fout.write(line)
            fin.close()
        fout.close()
        
        commandStr = 'samtools view -S -b {} > {}'.format(filtered_sam, filtered_bam)
        os.system(commandStr)
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
        if os.path.isfile(unmapped_sam):
            call(['rm', unmapped_sam])
        if os.path.isfile(filtered_sam):
            call(['rm', filtered_sam])
            
        write_log(log_file, flowcell_barcode, "Filter unmapped bam for "+library+" "+reference2+" in lane "+lane+" slice "+slice+" is done. ")
        
        call(['mv', folder_running, folder_finished])
    except:
        if os.path.isdir(folder_running):
            call(['mv', folder_running, folder_failed])
        else:
            call(['mkdir', folder_failed])
            
        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = "The Slide-seq workflow for "+library+" "+reference2+" in lane "+lane+" slice "+slice+" failed at the step of filtering unmapped bam. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)
        
        sys.exit()


if __name__ == "__main__":
    main()


