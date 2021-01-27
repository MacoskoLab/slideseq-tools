#!/usr/bin/python

# This script is to run the Slide-seq flowcell alignment pipeline

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
        logfile.write(dt_string+" [Slide-seq Flowcell Alignment Workflow - "+flowcell_barcode+"]: "+log_string+"\n")
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
    
    if 'flowcell_directory' not in options:
        print('flowcell_directory is not specified in the manifest file. Exiting...')
        sys.exit()
    
    if 'output_folder' not in options:
        print('output_folder is not specified in the manifest file. Exiting...')
        sys.exit()
    
    if 'metadata_file' not in options:
        print('metadata_file is not specified in the manifest file. Exiting...')
        sys.exit()
        
    if 'flowcell_barcode' not in options:
        print('flowcell_barcode is not specified in the manifest file. Exiting...')
        sys.exit()

    flowcell_directory = options['flowcell_directory']
    output_folder = options['output_folder']
    metadata_file = options['metadata_file']
    flowcell_barcode = options['flowcell_barcode']
    
    dropseq_folder = options['dropseq_folder'] if 'dropseq_folder' in options else '/broad/macosko/bin/dropseq-tools'
    picard_folder = options['picard_folder'] if 'picard_folder' in options else '/broad/macosko/bin/dropseq-tools/3rdParty/picard'
    STAR_folder = options['STAR_folder'] if 'STAR_folder' in options else '/broad/macosko/bin/dropseq-tools/3rdParty/STAR-2.5.2a'
    scripts_folder = options['scripts_folder'] if 'scripts_folder' in options else '/broad/macosko/jilong/slideseq_pipeline/scripts'
    email_address = options['email_address'] if 'email_address' in options else ''
    
    if not os.path.isdir(flowcell_directory):
        print("Folder {} does not exist. Exiting...".format(flowcell_directory))
        sys.exit()
    
    if not os.path.isfile(metadata_file):
        print("File {} does not exist. Exiting...".format(metadata_file))
        sys.exit()
    
    if not os.path.isdir(dropseq_folder):
        print("Folder {} does not exist. Exiting...".format(dropseq_folder))
        sys.exit()
    
    if not os.path.isdir(picard_folder):
        print("Folder {} does not exist. Exiting...".format(picard_folder))
        sys.exit()
    
    if not os.path.isdir(STAR_folder):
        print("Folder {} does not exist. Exiting...".format(STAR_folder))
        sys.exit()
    
    if not os.path.isdir(scripts_folder):
        print("Folder {} does not exist. Exiting...".format(scripts_folder))
        sys.exit()
    
    library_folder = options['library_folder'] if 'library_folder' in options else '{}/libraries'.format(output_folder)
    
    runinfo_file = '{}/RunInfo.xml'.format(flowcell_directory)
    if not os.path.isfile(runinfo_file):
        print("File {} does not exist. Exiting...".format(runinfo_file))
        sys.exit()
        
    try:
        # Create directories
        if not os.path.isdir(output_folder):
            call(['mkdir', '-p', output_folder])
        if not os.path.isdir('{}/logs'.format(output_folder)):
            call(['mkdir', '-p', '{}/logs'.format(output_folder)])
        call(['mkdir', '-p', '{}/status'.format(output_folder)])
        if not os.path.isdir(library_folder):
            call(['mkdir', '-p', library_folder])
        if 'temp_folder' not in options:
            call(['mkdir', '-p', '{}/tmp'.format(output_folder)])
    except:
        print("Folder {} cannot be created. Exiting...".format(output_folder))
        sys.exit()
    
    log_file = '{}/logs/workflow.log'.format(output_folder)
    write_log(log_file, flowcell_barcode, "The Slide-seq alignment pipeline is starting to run. ")
    
    # Call run_preparation
    output_file = '{}/logs/run_preparation.log'.format(output_folder)
    submission_script = '{}/run_preparation.sh'.format(scripts_folder)
    call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=54g', '-notify', '-l', 'h_rt=15:0:0', '-j', 'y', '-P', 'macosko_lab', '-l', 'os=RedHat7', submission_script, manifest_file, scripts_folder, output_folder]
    call(call_args)
    
    if len(email_address) > 1:
        subject = "Submission received for " + flowcell_barcode
        content = "Thank you for your interest on the Slide-seq tools! We received your request. An email will be sent to you once the workflow finishes. "
        call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
        call(call_args)


if __name__ == "__main__":
    main()


