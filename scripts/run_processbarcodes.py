#!/usr/bin/python

# This script is to extract Illumina barcodes

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
    if len(sys.argv) != 3:
        print("Please provide two arguments: manifest file and lane ID!")
        sys.exit()
    
    manifest_file = sys.argv[1]
    lane = sys.argv[2]

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
    
    # Get tile information from RunInfo.xml
    slice_id = {}
    slice_first_tile = {}
    slice_tile_limit = {}
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
    
    folder_running = '{}/status/running.processbarcodes_lane_{}'.format(output_folder, lane)
    folder_finished = '{}/status/finished.processbarcodes_lane_{}'.format(output_folder, lane)
    folder_failed = '{}/status/failed.processbarcodes_lane_{}'.format(output_folder, lane)
    
    try:
        call(['mkdir', folder_running])
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
    
        # Extract Illumina barcodes
        commandStr = 'java -Djava.io.tmpdir={} -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m '.format(tmpdir)
        commandStr += '-jar {}/picard.jar ExtractIlluminaBarcodes TMP_DIR={} VALIDATION_STRINGENCY=SILENT '.format(picard_folder, tmpdir)
        commandStr += 'BASECALLS_DIR={} OUTPUT_DIR={}/{}/barcodes LANE={} '.format(basecalls_dir, output_folder, lane, lane)
        commandStr += 'READ_STRUCTURE={} BARCODE_FILE={}/{}/barcode_params.txt '.format(read_structure, output_folder, lane)
        commandStr += 'METRICS_FILE={}/{}/{}.{}.barcode_metrics COMPRESS_OUTPUTS=true NUM_PROCESSORS=4'.format(output_folder, lane, flowcell_barcode, lane)
        write_log(log_file, flowcell_barcode, "ExtractIlluminaBarcodes for Lane " + lane + " Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "ExtractIlluminaBarcodes for Lane " + lane + " is done. ")
        
        # Convert Illumina base calls to sam (unmapped.bam)
        for i in range(len(slice_id[lane])):
            commandStr = 'java -Djava.io.tmpdir={} -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx10192m '.format(tmpdir)
            commandStr += '-jar {}/picard.jar IlluminaBasecallsToSam TMP_DIR={} VALIDATION_STRINGENCY=SILENT '.format(picard_folder, tmpdir)
            commandStr += 'BASECALLS_DIR={} LANE={} RUN_BARCODE={} NUM_PROCESSORS=4 '.format(basecalls_dir, lane, flowcell_barcode)
            commandStr += 'READ_STRUCTURE={} LIBRARY_PARAMS={}/{}/{}/library_params.txt INCLUDE_NON_PF_READS=false '.format(read_structure, output_folder, lane, slice_id[lane][i])
            commandStr += 'APPLY_EAMSS_FILTER=false MAX_READS_IN_RAM_PER_TILE=600000 ADAPTERS_TO_CHECK=null IGNORE_UNEXPECTED_BARCODES=true'
            commandStr += ' SEQUENCING_CENTER=BI BARCODES_DIR={}/{}/barcodes FIRST_TILE={} TILE_LIMIT={}'.format(output_folder, lane, slice_first_tile[lane][i], slice_tile_limit[lane][i])

            output_file = '{}/logs/run_barcodes2sam_lane_{}_{}.log'.format(output_folder, lane, slice_id[lane][i])
            submission_script = '{}/run.sh'.format(scripts_folder)
            call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=25g', '-notify', '-l', 'h_rt=10:0:0', '-j', 'y', submission_script, 'run_barcodes2sam', manifest_file, commandStr, lane, slice_id[lane][i], scripts_folder]
            call(call_args)
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
        call(['mv', folder_running, folder_finished])
    except:
        if os.path.isdir(folder_running):
            call(['mv', folder_running, folder_failed])
        else:
            call(['mkdir', folder_failed])
            
        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = "The Slide-seq workflow for lane "+lane+" failed at the step of processing barcodes. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)
        
        sys.exit()
    

if __name__ == "__main__":
    main()

