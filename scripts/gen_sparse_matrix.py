#!/usr/bin/python

# This script is to generate mtx files for sparse digital expression matrix

from __future__ import print_function

import sys
import os
import getopt
import csv
import gzip
import shutil
import scipy.io
import scipy.sparse
import io

import argparse
import glob
import re
from subprocess import call
from datetime import datetime

import numpy as np
# silence warnings for pandas below
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
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
    if len(sys.argv) != 6:
        print("Please provide five arguments: manifest file, library ID, locus function list, input folder and file name!")
        sys.exit()
    
    manifest_file = sys.argv[1]
    library = sys.argv[2]
    locus_function_list = sys.argv[3]
    input_folder = sys.argv[4]
    file_name = sys.argv[5]

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
    scripts_folder = options['scripts_folder']

    # Read info from metadata file
    reference = ''
    email_address = ''
    with open('{}/parsed_metadata.txt'.format(output_folder), 'r') as fin:
        reader = csv.reader(fin, delimiter='\t')
        rows = list(reader)
        row0 = rows[0]
        for i in range(1, len(rows)):
            row = rows[i]
            if row[row0.index('library')] == library:
                reference = row[row0.index('reference')]
                email_address = row[row0.index('email')]
                break
    fin.close()

    reference_folder = reference[:reference.rfind('/')]
    referencePure = reference[reference.rfind('/') + 1:]
    if (referencePure.endswith('.gz')):
        referencePure = referencePure[:referencePure.rfind('.')]
    referencePure = referencePure[:referencePure.rfind('.')]
    reference2 = referencePure + '.' + locus_function_list
    annotations_file = '{}/{}.gtf'.format(reference_folder, referencePure)

    log_file = '{}/logs/workflow.log'.format(output_folder)

    dge_gzfile = '{}/{}.txt.gz'.format(input_folder, file_name)
    dge_file = '{}/{}.txt'.format(input_folder, file_name)
    mat_file = '{}/{}_matrix.mtx'.format(input_folder, file_name)
    barcodes_file = '{}/{}_barcodes.tsv'.format(input_folder, file_name)
    genes_file = '{}/{}_features.tsv'.format(input_folder, file_name)
    
    folder_running = '{}/status/running.gen_sparse_matrix_{}_{}_{}'.format(output_folder, library, reference2, file_name)
    folder_finished = '{}/status/finished.gen_sparse_matrix_{}_{}_{}'.format(output_folder, library, reference2, file_name)
    folder_failed = '{}/status/failed.gen_sparse_matrix_{}_{}_{}'.format(output_folder, library, reference2, file_name)

    try:
        call(['mkdir', '-p', folder_running])
        
        write_log(log_file, flowcell_barcode, "Generating sparse matrix for "+library+" "+reference2+" "+file_name)
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
        os.system('gunzip -c {} > {}'.format(dge_gzfile, dge_file))
    
        # read column names
        with open(dge_file) as f:
            cols = f.readline().split('\t')
        f.close()
        
        # num of columns
        ncols = len(cols)
        
        # write barcodes file
        cols = cols[1:]
        with open(barcodes_file, 'w') as fout:
            for bc in cols:
                fout.write('{}-1\n'.format(bc.strip(' \t\n')))
        fout.close()

        # read gene id and name mapping from gtf file
        dict = {}
        with open(annotations_file, 'r') as fin:
            for line in fin:
                line = line.strip(' \t\n')
                if len(line) < 1 or line[0] == '#':
                    continue
                items = line.split('\t')[8]
                items = items.split(';')
                id = ''
                name = ''
                for item in items:
                    item = item.strip(' \t\n')
                    if item.split(' ')[0] == 'gene_id':
                        id = item.split(' ')[1]
                        id = id.strip(' \"')
                    elif item.split(' ')[0] == 'gene_name':
                        name = item.split(' ')[1]
                        name = name.strip(' \"')
                if not name in dict:
                    dict[name] = id
        fin.close()

        # write features (genes) file
        genes = np.loadtxt(dge_file, delimiter='\t', dtype='str', skiprows=1, usecols=(0))
        with open(genes_file, 'w') as fout:
            for gene in genes:
                if gene in dict:
                    fout.write('{}\t{}\n'.format(dict[gene], gene))
                else:
                    fout.write('{}\t{}\n'.format(gene, gene))
        fout.close()

        # load matrix
        data = np.loadtxt(dge_file, delimiter='\t', dtype='i', skiprows=1, usecols=range(1, ncols))
        
        data = scipy.sparse.csr_matrix(data)

        # write mtx file
        scipy.io.mmwrite(mat_file, data)

        if os.path.isfile(dge_file):
            call(['rm', dge_file])
            
        call(['gzip', mat_file])
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)

        write_log(log_file, flowcell_barcode, "Generating sparse matrix "+library+" "+reference2+" "+file_name+" is done. ")
        
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
            content = "The Slide-seq workflow for "+library+" "+reference2+" "+file_name+" failed at the step of generating sparse matrix. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)
        
        sys.exit()


if __name__ == "__main__":
    main()


