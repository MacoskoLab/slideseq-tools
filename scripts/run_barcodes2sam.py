#!/usr/bin/python

# This script is to convert Illumina barcodes to unmapped.bam

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


# Write to log file
def write_log(log_file, flowcell_barcode, log_string):
    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as logfile:
        logfile.write(dt_string+" [Slide-seq Flowcell Alignment Workflow - "+flowcell_barcode+"]: "+log_string+"\n")
    logfile.close()
    

def main():
    if len(sys.argv) != 5:
        print("Please provide four arguments: manifest file, commandStr, lane ID and slice ID!")
        sys.exit()
    
    manifest_file = sys.argv[1]
    commandStr = sys.argv[2]
    lane = sys.argv[3]
    slice = sys.argv[4]

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
    
    output_folder = options['output_folder']
    flowcell_barcode = options['flowcell_barcode']
    scripts_folder = options['scripts_folder'] if 'scripts_folder' in options else '/broad/macosko/jilong/slideseq_pipeline/scripts'
    email_address = options['email_address'] if 'email_address' in options else ''
    log_file = '{}/logs/workflow.log'.format(output_folder)
    
    folder_running = '{}/status/running.barcodes2sam_lane_{}_{}'.format(output_folder, lane, slice)
    folder_finished = '{}/status/finished.barcodes2sam_lane_{}_{}'.format(output_folder, lane, slice)
    folder_failed = '{}/status/failed.barcodes2sam_lane_{}_{}'.format(output_folder, lane, slice)
    
    try:
        call(['mkdir', '-p', folder_running])
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
        # Convert Illumina base calls to sam (unmapped.bam)
        write_log(log_file, flowcell_barcode, "IlluminaBasecallsToSam for Lane "+lane+'_'+slice+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "IlluminaBasecallsToSam for Lane "+lane+'_'+slice+" is done. ")
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)

        call(['mv', folder_running, folder_finished])
    except:
        if os.path.isdir(folder_running):
            call(['mv', folder_running, folder_failed])
        else:
            call(['mkdir', '-p', folder_failed])
            
        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = "The Slide-seq workflow for lane "+lane+" slice "+slice+" failed at the step of converting barcodes to sam. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)

        sys.exit()
            

if __name__ == "__main__":
    main()


