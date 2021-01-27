#!/usr/bin/python

# This script is to generate PDF for downsampling and projection

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

def my_smoother(data, xseq, **params):
    x, y = data['x'], data['y']
    coefs = np.polynomial.polynomial.polyfit(np.log(x), y, 1)
    ffit = np.polynomial.polynomial.polyval(np.log(xseq), coefs)
    data = pd.DataFrame({'x': xseq, 'y': ffit})
    return data

    
# Write to log file
def write_log(log_file, flowcell_barcode, log_string):
    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as logfile:
        logfile.write(dt_string+" [Slide-seq Flowcell Alignment Workflow - "+flowcell_barcode+"]: "+log_string+"\n")
        

def main():
    if len(sys.argv) != 4:
        print("Please provide three arguments: manifest file, library ID and locus function!")
        sys.exit()
    
    manifest_file = sys.argv[1]
    library = sys.argv[2]
    locus_function_list = sys.argv[3]
    
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
    
    flowcell_directory = options['flowcell_directory']
    output_folder = options['output_folder']
    metadata_file = options['metadata_file']
    flowcell_barcode = options['flowcell_barcode']
    
    library_folder = options['library_folder'] if 'library_folder' in options else '{}/libraries'.format(output_folder)
    scripts_folder = options['scripts_folder'] if 'scripts_folder' in options else '/broad/macosko/jilong/slideseq_pipeline/scripts'
    
    # Read info of lane, library and barcode
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    reference = ''
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
                email_address = row[row0.index('email')]
                experiment_date = row[row0.index('date')]
    fin.close()
    
    referencePure = reference[reference.rfind('/') + 1:]
    if (referencePure.endswith('.gz')):
        referencePure = referencePure[:referencePure.rfind('.')]
    referencePure = referencePure[:referencePure.rfind('.')]
    reference2 = referencePure + '.' + locus_function_list

    log_file = '{}/logs/workflow.log'.format(output_folder)
    
    downsample_folder = '{}/{}_{}/{}/downsample/'.format(library_folder, experiment_date, library, reference2)
    alignment_folder = '{}/{}_{}/{}/alignment'.format(library_folder, experiment_date, library, reference2)

    folder_running = '{}/status/running.gen_downsampling_plot_{}_{}'.format(output_folder, library, locus_function_list)
    folder_finished = '{}/status/finished.gen_downsampling_plot_{}_{}'.format(output_folder, library, locus_function_list)
    folder_failed = '{}/status/failed.gen_downsampling_plot_{}_{}'.format(output_folder, library, locus_function_list)
    
    # Wait till all of gen_downsample_dge finish
    ratio = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    while 1:
        for r in ratio:
            failed = '{}/status/failed.gen_downsample_dge_{}_{}_{}'.format(output_folder, library, locus_function_list, str(r))
            if os.path.isdir(failed):
                print('All gen_downsample_dge failed. Exiting...')
                sys.exit()
        
        f = True
        for r in ratio:
            finished = '{}/status/finished.gen_downsample_dge_{}_{}_{}'.format(output_folder, library, locus_function_list, str(r))
            if not os.path.isdir(finished):
                f = False
                break
        if f:
            break
        time.sleep(30)
        
    try:
        call(['mkdir', '-p', folder_running])

        ratio = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
        df_y = []
        for i in range(0, 10, 1):
            file = downsample_folder+library+'_'+str(ratio[i])+'.digital_expression_summary.txt'
            mat = np.loadtxt(file, delimiter='\t', dtype='int', skiprows=7, usecols=(2))
            v = int(sum(mat[:10000]) / min(len(mat),10000))
            df_y.append(v)
        df_x = np.tile(np.arange(0.1, 1.1, 0.1), 1)
        df = pd.DataFrame({"x":df_x, "y":df_y})
        p = ggplot(aes(x='x', y='y'), df) +\
        geom_point(size=2.5, color='red', fill='white', show_legend=False) +\
        geom_smooth(method=my_smoother, se=False, fullrange=True, size=0.4, color='red', fill='white', show_legend=False) +\
        scale_x_continuous(limits=(0.1,10), breaks=(0.5,1.0,2.0,5.0,10.0), labels=[0.5,1.0,2.0,5.0,10.0]) +\
        theme(axis_text_y=element_text(rotation=90, hjust=1)) +\
        xlab("Reads generated (relative to this run)") +\
        ylab("Transcripts per cell") +\
        ggtitle("Return to sequencing coverage (downsampling + projection)")
        ggsave(plot=p, height=8, width=8, filename=library+'_'+reference2+"_downsampling.pdf", path=alignment_folder)
        
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
            content = "The Slide-seq workflow for "+library+" "+reference2+" failed at the step of generating plot for downsampling and projection. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)
        
        sys.exit()
    

if __name__ == "__main__":
    main()
    
    
