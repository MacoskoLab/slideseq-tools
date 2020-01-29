#!/usr/bin/python

# This script is to call the alignment step

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


# Get read structure from RunInfo.xml
def get_read_structure(x):
    posts = 'TBT'
    i = 0
    str = ''
    with open(x, 'r') as fin:
        for line in fin:
            line = line.strip(' \t\n')
            if (line.startswith('<Read ', 0)):
                l = line.split('=')[2]
                l = l.split('\"')[1]
                str += l + posts[i]
                i += 1
                if (i >= 3):
                    break;
    fin.close()
    return str


# Get tile information from RunInfo.xml
def get_tiles(x, lane):
    tiles = []
    with open(x, 'r') as fin:
        for line in fin:
            line = line.strip(' \t\n')
            if (line.startswith('<Tile>', 0)):
                l = line[6:].split('<')[0]
                if (l.split('_')[0] == lane):
                    tiles.append(l.split('_')[1])
    fin.close()
    tiles.sort()
    return tiles


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
    
    flowcell_directory = options['flowcell_directory']
    dropseq_folder = options['dropseq_folder']
    picard_folder = options['picard_folder']
    STAR_folder = options['STAR_folder']
    scripts_folder = options['scripts_folder']
    output_folder = options['output_folder']
    metadata_file = options['metadata_file']
    option_file = options['option_file']
    flowcell_barcode = options['flowcell_barcode']
    
    library_folder = options['library_folder'] if 'library_folder' in options else '{}/libraries'.format(output_folder)
    tmpdir = options['temp_folder'] if 'temp_folder' in options else '{}/tmp'.format(output_folder)
    illumina_platform = options['illumina_platform'] if 'illumina_platform' in options else 'NextSeq'
    email_address = options['email_address'] if 'email_address' in options else ''
    
    is_NovaSeq = True if illumina_platform == 'NovaSeq' else False
    is_NovaSeq_S4 = True if illumina_platform == 'NovaSeq_S4' else False
    num_slice_NovaSeq = 10
    num_slice_NovaSeq_S4 = 40

    runinfo_file = '{}/RunInfo.xml'.format(flowcell_directory)
    log_file = '{}/logs/workflow.log'.format(output_folder)
    
    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    references_unique = []
    locus_function_list_unique = []
    experiment_date = []
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
                references_unique.append(row[row0.index('reference')])
                locus_function_list_unique.append(row[row0.index('locus_function_list')])
                experiment_date.append(row[row0.index('experiment_date')])
            barcodes.append(row[row0.index('sample_barcode')])
            bead_structures.append(row[row0.index('bead_structure')])
    fin.close()

    # Get tile information from RunInfo.xml
    slice_id = {}
    slice_first_tile = {}
    slice_tile_limit = {}
    for lane in lanes_unique:
        tile_nums = get_tiles(runinfo_file, lane)
        tile_cou = len(tile_nums)
        if ((not is_NovaSeq) and (not is_NovaSeq_S4)):
            slice_id[lane] = ['0']
            slice_first_tile[lane] = [str(tile_nums[0])]
            slice_tile_limit[lane] = [str(tile_cou)]
        else:
            slice_cou = num_slice_NovaSeq if is_NovaSeq else num_slice_NovaSeq_S4
            tile_cou_per_slice = (tile_cou // slice_cou) + 1
            slice_id[lane] = []
            slice_first_tile[lane] = []
            slice_tile_limit[lane] = []
            for i in range(slice_cou):
                if tile_cou_per_slice * i >= tile_cou:
                    break
                slice_id[lane].append(str(i))
                slice_first_tile[lane].append(str(tile_nums[tile_cou_per_slice * i]))
                slice_tile_limit[lane].append(str(tile_cou_per_slice))
    
    folder_waiting = '{}/status/waiting.mergebarcodes'.format(output_folder)
    folder_running = '{}/status/running.mergebarcodes'.format(output_folder)
    folder_finished = '{}/status/finished.mergebarcodes'.format(output_folder)
    folder_failed = '{}/status/failed.mergebarcodes'.format(output_folder)
    
    call(['mkdir', folder_waiting])
    
    # Wait till all of run_processbarcodes and run_barcodes2sam finish
    while 1:
        all_failed_1 = True
        all_failed_2 = True
        for lane in lanes_unique:
            failed_1 = '{}/status/failed.processbarcodes_lane_{}'.format(output_folder, lane)
            if not os.path.isdir(failed_1):
                all_failed_1 = False
            for i in range(len(slice_id[lane])):
                failed_2 = '{}/status/failed.barcodes2sam_lane_{}_{}'.format(output_folder, lane, slice_id[lane][i])
                if not os.path.isdir(failed_2):
                    all_failed_2 = False
                    break
        if all_failed_1:
            write_log(log_file, flowcell_barcode, "All run_processbarcodes failed. Exiting...")
            sys.exit()
        if all_failed_2:
            write_log(log_file, flowcell_barcode, "All run_barcodes2sam failed. Exiting...")
            sys.exit()
        
        f = True
        for lane in lanes_unique:
            for i in range(len(slice_id[lane])):
                fol1 = '{}/status/finished.barcodes2sam_lane_{}_{}'.format(output_folder, lane, slice_id[lane][i])
                fol2 = '{}/status/failed.barcodes2sam_lane_{}_{}'.format(output_folder, lane, slice_id[lane][i])
                if (not os.path.isdir(fol1)) and (not os.path.isdir(fol2)):
                    f = False
                    break
            if not f:
                break
        if f:
            break
        time.sleep(30)
    
    if os.path.isdir(folder_waiting):
        call(['mv', folder_waiting, folder_running])
    else:
        call(['mkdir', folder_running])
        
    try:       
        for j in range(len(libraries_unique)):
            library = libraries_unique[j]
            os.system('mkdir ' + '{}/{}_{}'.format(library_folder, experiment_date[j], library))
            for i in range(len(lanes)):
                if libraries[i] != library:
                    continue
                for slice in slice_id[lanes[i]]:
                    unmapped_bam = '{}/{}/{}/{}/'.format(output_folder, lanes[i], slice, library)
                    if (barcodes[i]):
                        unmapped_bam += '{}/{}.{}.{}.{}.{}.unmapped.bam'.format(barcodes[i], flowcell_barcode, lanes[i], slice, library, barcodes[i])
                    else:
                        unmapped_bam += '{}.{}.{}.{}.unmapped.bam'.format(flowcell_barcode, lanes[i], slice, library)
                    unmapped_bam2 = '{}/{}_{}/{}.{}.{}.{}'.format(library_folder, experiment_date[j], library, flowcell_barcode, lanes[i], slice, library)
                    if (barcodes[i]):
                        unmapped_bam2 += '.'+barcodes[i]
                    unmapped_bam2 += '.unmapped.bam'
                    if os.path.isfile(unmapped_bam):
                        os.system('mv ' + unmapped_bam + ' ' + unmapped_bam2)
    
                    # Call run_alignment
                    output_file = '{}/logs/run_alignment_{}_{}_{}.log'.format(output_folder, library, lanes[i], slice)
                    submission_script = '{}/run.sh'.format(scripts_folder)
                    call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=35g', '-notify', '-l', 'h_rt=10:0:0', '-j', 'y', submission_script, 'run_alignment', manifest_file, library, lanes[i], slice, scripts_folder]
                    call(call_args)
    
        # Call run_analysis
        for j in range(len(libraries_unique)):
            library = libraries_unique[j]
            output_file = '{}/logs/run_analysis_{}.log'.format(output_folder, library)
            submission_script = '{}/run.sh'.format(scripts_folder)
            call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=30g', '-notify', '-l', 'h_rt=30:0:0', '-j', 'y', submission_script, 'run_analysis', manifest_file, library, scripts_folder]
            call(call_args)
        
        call(['mv', folder_running, folder_finished])
    except:
        if os.path.isdir(folder_running):
            call(['mv', folder_running, folder_failed])
        elif os.path.isdir(folder_waiting):
            call(['mv', folder_waiting, folder_failed])
        else:
            call(['mkdir', folder_failed])
            
        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = "The Slide-seq workflow failed at the step of merging barcode matrics. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)

        sys.exit()
        
    
if __name__ == "__main__":
    main()


