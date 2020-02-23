#!/usr/bin/python

# This script is to check the Illumina directory, parse input data, and call the steps 
# of extracting Illumina barcodes and converting barcodes to bam files

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
                str += l.split('\"')[1] + posts[i]
                i += 1
                if (i >= 3):
                    break
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
    
    basecalls_dir = '{}/Data/Intensities/BaseCalls'.format(flowcell_directory)
    
    # Get read structure from RunInfo.xml
    runinfo_file = '{}/RunInfo.xml'.format(flowcell_directory)
    read_structure = get_read_structure(runinfo_file)

    log_file = '{}/logs/workflow.log'.format(output_folder)

    # Parse metadata file
    write_log(log_file, flowcell_barcode, "Parse metadata file. ")
    commandStr = 'python {}/parse_metadata.py -i {} -r {} -o {}/parsed_metadata.txt'.format(scripts_folder, metadata_file, runinfo_file, output_folder)
    os.system(commandStr)
    
    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    references_unique = []
    locus_function_list_unique = []
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
            barcodes.append(row[row0.index('sample_barcode')])
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
                if (tile_cou_per_slice * i >= tile_cou):
                    break
                slice_id[lane].append(str(i))
                slice_first_tile[lane].append(str(tile_nums[tile_cou_per_slice * i]))
                slice_tile_limit[lane].append(str(tile_cou_per_slice))
        
    folder_running = '{}/status/running.run_preparation'.format(output_folder)
    folder_finished = '{}/status/finished.run_preparation'.format(output_folder)
    folder_failed = '{}/status/failed.run_preparation'.format(output_folder)

    try:
        call(['mkdir', folder_running])
        
        # Check the input Illumina folder
        commandStr = 'java -Djava.io.tmpdir={} -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8192m '.format(tmpdir)
        commandStr += '-jar {}/picard.jar CheckIlluminaDirectory TMP_DIR={} VALIDATION_STRINGENCY=SILENT '.format(picard_folder, tmpdir)
        commandStr += 'BASECALLS_DIR={} READ_STRUCTURE={}'.format(basecalls_dir, read_structure)
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += ' LINK_LOCS=false'
        for lane in lanes_unique:
            commandStr += ' L=' + lane
        write_log(log_file, flowcell_barcode, "CheckIlluminaDirectory Command=" + commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "CheckIlluminaDirectory is done. ")
        
        # Create directories
        write_log(log_file, flowcell_barcode, "Creating directories. ")
        for lane in lanes_unique:
            call(['mkdir', '{}/{}'.format(output_folder, lane)])
            call(['mkdir', '{}/{}/barcodes'.format(output_folder, lane)])
            for slice in slice_id[lane]:
                call(['mkdir', '{}/{}/{}'.format(output_folder, lane, slice)])
        for i in range(len(lanes)):
            for slice in slice_id[lane]:
                if not os.path.isdir('{}/{}/{}/{}'.format(output_folder, lanes[i], slice, libraries[i])):
                    call(['mkdir', '{}/{}/{}/{}'.format(output_folder, lanes[i], slice, libraries[i])])
                if (barcodes[i]):
                    call(['mkdir', '{}/{}/{}/{}/{}'.format(output_folder, lanes[i], slice, libraries[i], barcodes[i])])
        
        # Generate barcode_params.txt that is needed by ExtractIlluminaBarcodes
        for lane in lanes_unique:
            write_log(log_file, flowcell_barcode, "Generating barcode_params.txt for Lane " + lane)
            commandStr = 'python {}/gen_barcode_params.py -i {}/parsed_metadata.txt -o {}/{}/barcode_params.txt -l {}'.format(scripts_folder, output_folder, output_folder, lane, lane)
            os.system(commandStr)

        # Generate library_params that is needed by IlluminaBasecallsToSam
        for lane in lanes_unique:
            write_log(log_file, flowcell_barcode, "Generating library_params.txt for Lane " + lane)
            for slice in slice_id[lane]:
                commandStr = 'python {}/gen_library_params.py -i {}/parsed_metadata.txt -o {}/{}/{}/library_params.txt -b '.format(scripts_folder, output_folder, output_folder, lane, slice)
                commandStr += '{}/{}/{}/ -n {}.{}.{} -l {}'.format(output_folder, lane, slice, flowcell_barcode, lane, slice, lane)
                os.system(commandStr)
        
        # Call run_processbarcodes
        for lane in lanes_unique:
            output_file = '{}/logs/run_processbarcodes_lane_{}.log'.format(output_folder, lane)
            submission_script = '{}/run.sh'.format(scripts_folder)
            call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=30g', '-notify', '-l', 'h_rt=10:0:0', '-j', 'y', submission_script, 'run_processbarcodes', manifest_file, lane, scripts_folder]
            call(call_args)
        
        # Call run_mergebarcodes
        output_file = '{}/logs/run_mergebarcodes.log'.format(output_folder)
        submission_script = '{}/run.sh'.format(scripts_folder)
        call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=4g', '-notify', '-l', 'h_rt=40:0:0', '-j', 'y', submission_script, 'run_mergebarcodes', manifest_file, scripts_folder]
        call(call_args)
    
        call(['mv', folder_running, folder_finished])
    except:
        if os.path.isdir(folder_running):
            call(['mv', folder_running, folder_failed])
        else:
            call(['mkdir', folder_failed])
            
        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = "The Slide-seq workflow failed at the step of preparation. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)
        
        sys.exit()
    

if __name__ == "__main__":
    main()

