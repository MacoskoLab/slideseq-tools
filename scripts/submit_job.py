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
    
    if 'flowcell_directory' not in options:
        print('flowcell_directory is not specified in the manifest file. Exiting...')
        sys.exit()
    
    if 'dropseq_folder' not in options:
        print('dropseq_folder is not specified in the manifest file. Exiting...')
        sys.exit()
        
    if 'picard_folder' not in options:
        print('picard_folder is not specified in the manifest file. Exiting...')
        sys.exit()

    if 'STAR_folder' not in options:
        print('STAR_folder is not specified in the manifest file. Exiting...')
        sys.exit()

    if 'scripts_folder' not in options:
        print('scripts_folder is not specified in the manifest file. Exiting...')
        sys.exit()

    if 'output_folder' not in options:
        print('output_folder is not specified in the manifest file. Exiting...')
        sys.exit()
    
    if 'metadata_file' not in options:
        print('metadata_file is not specified in the manifest file. Exiting...')
        sys.exit()
        
    if 'option_file' not in options:
        print('option_file is not specified in the manifest file. Exiting...')
        sys.exit()
    
    if 'flowcell_barcode' not in options:
        print('flowcell_barcode is not specified in the manifest file. Exiting...')
        sys.exit()
        
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
    email_address = options['email_address'] if 'email_address' in options else ''
    
    if not os.path.isdir(flowcell_directory):
        print("Folder {} does not exist. Exiting...".format(flowcell_directory))
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
    
    if not os.path.isfile(metadata_file):
        print("File {} does not exist. Exiting...".format(metadata_file))
        sys.exit()
    
    runinfo_file = '{}/RunInfo.xml'.format(flowcell_directory)
    if not os.path.isfile(runinfo_file):
        print("File {} does not exist. Exiting...".format(runinfo_file))
        sys.exit()

    try:
        # Create directories
        if os.path.isdir(output_folder):
            print("Folder {} already exists. Exiting...".format(output_folder))
            sys.exit()
        call(['mkdir', output_folder])
        call(['mkdir', '{}/logs'.format(output_folder)])
        call(['mkdir', '{}/status'.format(output_folder)])
        if not os.path.isdir(library_folder):
            call(['mkdir', library_folder])
        if not os.path.isdir(tmpdir):
            call(['mkdir', tmpdir])
    except:
        print("Folder {} cannot be created. Exiting...".format(output_folder))
        sys.exit()
    
    log_file = '{}/logs/workflow.log'.format(output_folder)
    write_log(log_file, flowcell_barcode, "The Slide-seq alignment pipeline is starting to run. ")
    
    # Call run_preparation
    output_file = '{}/logs/run_preparation.log'.format(output_folder)
    submission_script = '{}/run.sh'.format(scripts_folder)
    call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=4g', '-notify', '-l', 'h_rt=1:0:0', '-j', 'y', submission_script, 'run_preparation', manifest_file, scripts_folder]
    call(call_args)
    
    if len(email_address) > 1:
        subject = "Submission received for " + flowcell_barcode
        content = "Thank you for your interest on the Slide-seq tools! We received your request. An email will be sent to you once the workflow finishes. "
        call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
        call(call_args)


if __name__ == "__main__":
    main()


