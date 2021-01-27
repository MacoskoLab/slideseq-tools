#!/usr/bin/python

# This script is to calculate UP distance and generate plot

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

import numpy as np
# silence warnings for pandas below
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import pandas as pd

import random
from random import sample

from plotnine import *

def levenshtein(seq1, seq2):
    size_x = len(seq1) + 1
    size_y = len(seq2) + 1
    matrix = np.zeros ((size_x, size_y))
    for x in range(size_x):
        matrix [x, 0] = x
    for y in range(size_y):
        matrix [0, y] = y

    for x in range(1, size_x):
        for y in range(1, size_y):
            if seq1[x-1] == seq2[y-1]:
                matrix [x,y] = min(
                    matrix[x-1, y] + 1,
                    matrix[x-1, y-1],
                    matrix[x, y-1] + 1
                )
            else:
                matrix [x,y] = min(
                    matrix[x-1,y] + 1,
                    matrix[x-1,y-1] + 1,
                    matrix[x,y-1] + 1
                )
    return (matrix[size_x - 1, size_y - 1])
    
    
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
    if len(sys.argv) != 4:
        print("Please provide two arguments: manifest file, library ID and lane ID!")
        sys.exit()
    
    manifest_file = sys.argv[1]
    library = sys.argv[2]
    this_lane = sys.argv[3]

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
    
    runinfo_file = '{}/RunInfo.xml'.format(flowcell_directory)
    log_file = '{}/logs/workflow.log'.format(output_folder)
    
    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    reference = ''
    locus_function_list = 'exonic+intronic'
    email_address = ''
    experiment_date = ''
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
                locus_function_list = row[row0.index('locus_function_list')]
                email_address = row[row0.index('email')]
                experiment_date = row[row0.index('date')]
    fin.close()
    
    folder_running = '{}/status/running.analysis_UPdistance_{}_{}'.format(output_folder, library, this_lane)
    folder_finished = '{}/status/finished.analysis_UPdistance_{}_{}'.format(output_folder, library, this_lane)
    folder_failed = '{}/status/failed.analysis_UPdistance_{}_{}'.format(output_folder, library, this_lane)

    try:
        call(['mkdir', '-p', folder_running])
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
        write_log(log_file, flowcell_barcode, "Start to calculate UP distance for "+library+" in lane "+this_lane)
        
        analysis_folder = '{}/{}_{}'.format(library_folder, experiment_date, library)
        read1_file = '{}/{}.{}.read1.fastq'.format(analysis_folder, library, this_lane)
        read1_gzfile = '{}/{}.{}.read1.fastq.gz'.format(analysis_folder, library, this_lane)
        summary_file = '{}/{}.UPsummary.txt'.format(analysis_folder, library)
        UPbases_file = '{}/{}.UPbases.txt'.format(analysis_folder, library)
        UPdistances_filename = '{}.UPdistances.pdf'.format(library)
        
        os.system('gzip -c '+read1_file+' > '+read1_gzfile)
        
        if os.path.isfile(read1_file):
            call(['rm', read1_file])
        
        #Total number of reads in read 1
        with open(summary_file, "a") as summaryfile:
            summaryfile.write("Total number of reads in read 1\n")
        summaryfile.close()
        os.system('zcat '+read1_gzfile+'| grep -cP "[AGTC]{42}" >> '+summary_file)

        #Total number of reads with full length UP
        with open(summary_file, "a") as summaryfile:
            summaryfile.write("\nTotal number of reads with full length UP\n")
        summaryfile.close()
        os.system('zcat '+read1_gzfile+' | grep -cP "TCTTCAGCGTTCCCGAGA" >> '+summary_file)

        #Total number of reads with full length UP starting at the right position
        with open(summary_file, "a") as summaryfile:
            summaryfile.write("\nTotal number of reads with full length UP starting at the right position\n")
        summaryfile.close()
        os.system('zcat '+read1_gzfile+' | grep -cP "[AGTC]{8}TCTTCAGCGTTCCCGAGA[AGTC]{16}" >> '+summary_file)

        #Total number of reads with full length UP starting at the -1 position
        with open(summary_file, "a") as summaryfile:
            summaryfile.write("\nTotal number of reads with full length UP starting at the -1 position\n")
        summaryfile.close()
        os.system('zcat '+read1_gzfile+' | grep -cP "[AGTC]{7}TCTTCAGCGTTCCCGAGA[AGTC]{17}" >> '+summary_file)

        #Extract bases at the UP position and save to file
        os.system('zcat '+read1_gzfile+'| grep -P "[AGTC]{42}" | cut -c 9-26 > '+UPbases_file)
        
        if os.path.isfile(read1_gzfile):
            call(['rm', read1_gzfile])

        # Plot
        with open(UPbases_file) as f:
            x = f.read().splitlines()
        f.close()
        
        #Randomly sample 20,000 sequences to get a good representation
        random.seed(1)

        #k = 20000
        k = min(len(x), 1000000)
        idx = sample(x, k)

        a = []
        for i in range(0, k):
            a.append(levenshtein(idx[i], "TCTTCAGCGTTCCCGAGA"))
        
        shuffles = []
        for i in range(0, k):
            l = list(idx[i])
            random.shuffle(l)
            shuffles.append(''.join(l))

        shuffdist = []
        for i in range(0, k):
            shuffdist.append(levenshtein(shuffles[i], "TCTTCAGCGTTCCCGAGA"))
        
        df = pd.DataFrame({
            "Real": a,
            "Shuffled": shuffdist
        })
        df = pd.melt(df)

        p = ggplot(aes(x='value', fill='variable'), data=df) +\
        geom_histogram(alpha=0.3, binwidth=1, position="identity") +\
        xlab("Distance") +\
        ylab("Count") +\
        ggtitle("Distance vs Count")
        ggsave(plot = p, height=10, width=10, filename=UPdistances_filename, path=analysis_folder)
        
        if os.path.isfile(UPbases_file):
            call(['rm', UPbases_file])
            
        write_log(log_file, flowcell_barcode, "Calculating UP distance for "+library+" in lane "+this_lane+" is done. ")
        
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
            content = "The Slide-seq workflow for "+library+" in lane "+this_lane+" failed at the step of running analysis UPdistance. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)

        sys.exit()
    

if __name__ == "__main__":
    main()


