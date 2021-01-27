#!/usr/bin/python

# This script is to combine outputs from cmatcher_beads

from __future__ import print_function

import sys
import os
import getopt
import csv
import gzip
import shutil

import argparse
import glob
import re
import time
from subprocess import call
from datetime import datetime

import numpy as np
# silence warnings for pandas below
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import pandas as pd

import random
from random import sample

import itertools

import math
import numpy.polynomial.polynomial as poly

import plotnine as pn
from plotnine import *

import matplotlib as mpl
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
    if len(sys.argv) != 3:
        print("Please provide two arguments: manifest file and library ID!")
        sys.exit()
    
    manifest_file = sys.argv[1]
    library = sys.argv[2]
    
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
    
    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    reference = ''
    base_quality = '10'
    min_transcripts_per_cell = '10'
    email_address = ''
    bead_type = '180402'
    bead_structure = ''
    run_puckmatcher = False
    experiment_date = ''
    gen_read1_plot = False
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
                base_quality = row[row0.index('base_quality')]
                min_transcripts_per_cell = row[row0.index('min_transcripts_per_cell')]
                email_address = row[row0.index('email')]
                bead_type = row[row0.index('bead_type')]
                bead_structure = row[row0.index('bead_structure')]
                run_puckmatcher = str2bool(row[row0.index('run_barcodematching')])
                experiment_date = row[row0.index('date')]
                if 'gen_read1_plot' in row0:
                    gen_read1_plot = str2bool(row[row0.index('gen_read1_plot')])
    fin.close()

    runinfo_file = '{}/RunInfo.xml'.format(flowcell_directory)
    log_file = '{}/logs/workflow.log'.format(output_folder)

    analysis_folder = '{}/{}_{}'.format(library_folder, experiment_date, library)
    bead_barcode_file = '{}/BeadBarcodes.txt'.format(analysis_folder)
    bead_location_file = '{}/BeadLocations.txt'.format(analysis_folder)
    
    if not os.path.isfile(bead_barcode_file):
        write_log(log_file, flowcell_barcode, 'run_cmatcher_beads_combine error: '+bead_barcode_file+' does not exist!')
        raise Exception('run_cmatcher_beads_combine error: '+bead_barcode_file+' does not exist!')
    
    folder_running = '{}/status/running.cmatcher_beads_combine_{}'.format(output_folder, library)
    folder_finished = '{}/status/finished.cmatcher_beads_combine_{}'.format(output_folder, library)
    folder_failed = '{}/status/failed.cmatcher_beads_combine_{}'.format(output_folder, library)
    
    try:        
        call(['mkdir', '-p', folder_running])
        
        l = 0
        with open(bead_barcode_file, 'r') as fin:
            for line in fin:
                l += 1
        fin.close()
        k = 10000
        ls = l // k
        
        while 1:
            f = True
            for i in range(ls + 1):
                if i * k >= l:
                    break;
                file2 = '{}/{}_barcode_matching_01_{}.finished'.format(analysis_folder, library, str(i + 1))
                if not os.path.isfile(file2):
                    f = False
                    break
            if f:
                break
            time.sleep(30)

        print('combine cmatcher_beads outputs...')
        write_log(log_file, flowcell_barcode, "Combine cmatcher_beads outputs for "+library)
        combined_cmatcher_file = '{}/{}_barcode_matching_01.txt'.format(analysis_folder, library)
        with open(combined_cmatcher_file, 'w') as fout:
            for i in range(ls + 1):
                if i * k >= l:
                    break;
                file2 = '{}/{}_barcode_matching_01_{}.txt'.format(analysis_folder, library, str(i + 1))
                with open(file2, 'r') as fin:
                    for line in fin:
                        fout.write(line)
                fin.close()
        fout.close()
        
        combined_cmatcher_file2 = '{}/{}_barcode_matching_2.txt'.format(analysis_folder, library)
        with open(combined_cmatcher_file2, 'w') as fout:
            for i in range(ls + 1):
                if i * k >= l:
                    break;
                file2 = '{}/{}_barcode_matching_2_{}.txt'.format(analysis_folder, library, str(i + 1))
                with open(file2, 'r') as fin:
                    for line in fin:
                        fout.write(line)
                fin.close()
        fout.close()
        
        for i in range(ls + 1):
            if i * k >= l:
                break;
            file1 = '{}/{}_barcode_matching_01_{}.txt'.format(analysis_folder, library, str(i + 1))
            file2 = '{}/{}_barcode_matching_01_{}.finished'.format(analysis_folder, library, str(i + 1))
            file3 = '{}/BeadBarcodes_{}.txt'.format(analysis_folder, str(i + 1))
            file4 = '{}/{}_barcode_matching_2_{}.txt'.format(analysis_folder, library, str(i + 1))
            if os.path.isfile(file1):
                call(['rm', file1])
            if os.path.isfile(file2):
                call(['rm', file2])
            if os.path.isfile(file3):
                call(['rm', file3])
            if os.path.isfile(file4):
                call(['rm', file4])
        
        write_log(log_file, flowcell_barcode, "Combine cmatcher_beads outputs for "+library+" is done. ")
        
        # Create degenerate bead barcodes
        print('Create degenerate bead barcodes...')
        combined_cmatcher_file3 = '{}/BeadBarcodes_degenerate.txt'.format(analysis_folder)
        commandStr = scripts_folder+'/degenerate_beads '+bead_barcode_file+' '+bead_location_file+' '+combined_cmatcher_file+' '+combined_cmatcher_file2+' '+combined_cmatcher_file3
        os.system(commandStr)
        
        file = '{}/BeadBarcodes_degenerate.finished'.format(analysis_folder)
        with open(file, 'w') as fout:
            fout.write('finished')
        fout.close()
        
        write_log(log_file, flowcell_barcode, "Create degenerate bead barcodes for "+library+" is done. ")
        
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
            content = "The Slide-seq workflow for "+library+" failed at the step of running cmatcher beads combine. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)

        sys.exit()


if __name__ == "__main__":
    main()
    
    
