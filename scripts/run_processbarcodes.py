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

from new_submit_to_taskrunner import call_to_taskrunner
import traceback

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
        logfile.write(dt_string+" [Slide-seq Flowcell Alignment Workflow - "+flowcell_barcode+"]: "+log_string+"\n")
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
    email_address = options['email_address'] if 'email_address' in options else ''
    
    basecalls_dir = '{}/Data/Intensities/BaseCalls'.format(flowcell_directory)
    log_file = '{}/logs/workflow.log'.format(output_folder)
    
    # Get read structure from RunInfo.xml
    runinfo_file = '{}/RunInfo.xml'.format(flowcell_directory)
    read_structure = get_read_structure(runinfo_file)
    
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
        call(['mkdir', '-p', folder_running])
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
    
        # Extract Illumina barcodes
        commandStr = 'java -Djava.io.tmpdir='+tmpdir+' -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m '
        commandStr += '-jar '+picard_folder+'/picard.jar ExtractIlluminaBarcodes TMP_DIR='+tmpdir+' VALIDATION_STRINGENCY=SILENT '
        commandStr += 'BASECALLS_DIR='+basecalls_dir+' OUTPUT_DIR='+output_folder+'/'+lane+'/barcodes LANE='+lane+' '
        commandStr += 'READ_STRUCTURE='+read_structure+' BARCODE_FILE='+output_folder+'/'+lane+'/barcode_params.txt '
        commandStr += 'METRICS_FILE='+output_folder+'/'+lane+'/'+flowcell_barcode+'.'+lane+'.barcode_metrics COMPRESS_OUTPUTS=true NUM_PROCESSORS=4'
        write_log(log_file, flowcell_barcode, "ExtractIlluminaBarcodes for Lane "+lane+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "ExtractIlluminaBarcodes for Lane "+lane+" is done. ")
        
        # Convert Illumina base calls to sam (unmapped.bam)
        for i in range(len(slice_id[lane])):
            commandStr = 'java -Djava.io.tmpdir='+tmpdir+' -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx10192m '
            commandStr += '-jar '+picard_folder+'/picard.jar IlluminaBasecallsToSam TMP_DIR='+tmpdir+' VALIDATION_STRINGENCY=SILENT '
            commandStr += 'BASECALLS_DIR='+basecalls_dir+' LANE='+lane+' RUN_BARCODE='+flowcell_barcode+' NUM_PROCESSORS=4 '
            commandStr += 'READ_STRUCTURE='+read_structure+' LIBRARY_PARAMS='+output_folder+'/'+lane+'/'+slice_id[lane][i]+'/library_params.txt INCLUDE_NON_PF_READS=false '
            commandStr += 'APPLY_EAMSS_FILTER=false MAX_READS_IN_RAM_PER_TILE=600000 ADAPTERS_TO_CHECK=null IGNORE_UNEXPECTED_BARCODES=true'
            commandStr += ' SEQUENCING_CENTER=BI BARCODES_DIR='+output_folder+'/'+lane+'/barcodes FIRST_TILE='+slice_first_tile[lane][i]+' TILE_LIMIT='+slice_tile_limit[lane][i]

            output_file = '{}/logs/run_barcodes2sam_lane_{}_{}.log'.format(output_folder, lane, slice_id[lane][i])
            submission_script = '{}/run_barcodes2sam.sh'.format(scripts_folder)
            call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=150G', '-notify', '-l', 'h_rt=80:0:0', '-j', 'y', '-P', 'macosko_lab', '-l', 'os=RedHat7', submission_script, manifest_file, commandStr, lane, slice_id[lane][i], scripts_folder, output_folder, '{}/{}'.format(output_folder, lane)]
            call_to_taskrunner(output_folder, call_args)
        
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
        else:
            call(['mkdir', '-p', folder_failed])
            
        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = "The Slide-seq workflow for lane "+lane+" failed at the step of processing barcodes. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)
        
        sys.exit()
    

if __name__ == "__main__":
    main()


