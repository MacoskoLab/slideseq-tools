#!/usr/bin/python

# This script is to combine outputs from cmatcher and call tag_matched_bam and filter_unmapped_bam

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


# Get bead structure range
def get_bead_structure_range(bs, type):
    #12C8M|*T
    #7C18X7C8M2X|*T
    l = re.split('C|X|M', re.split('\|', bs)[0])
    res = ''
    i = 1
    p = -1
    for it in l:
        if it:
            p += len(it) + 1
            if bs[p] == type:
                res += str(i) + '-' + str(i + int(it) - 1) + ':'
            i += int(it)
    return res[:-1]


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
    if len(sys.argv) != 4:
        print("Please provide three arguments: manifest file, library ID and locus function list!")
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
    bead_structure = ''
    experiment_date = ''
    run_puckmatcher = False
    bead_type = '180402'
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
                bead_structure = row[row0.index('bead_structure')]
                run_puckmatcher = str2bool(row[row0.index('run_barcodematching')])
                experiment_date = row[row0.index('experiment_date')]
    fin.close()
    
    reference_folder = reference[:reference.rfind('/')]
    referencePure = reference[reference.rfind('/') + 1:]
    if (referencePure.endswith('.gz')):
        referencePure = referencePure[:referencePure.rfind('.')]
    referencePure = referencePure[:referencePure.rfind('.')]
    genome_dir = '{}/STAR'.format(reference_folder)
    intervals = '{}/{}.genes.intervals'.format(reference_folder, referencePure)
    annotations_file = '{}/{}.gtf'.format(reference_folder, referencePure)
    ref_flat = '{}/{}.refFlat'.format(reference_folder, referencePure)
    ribosomal_intervals = '{}/{}.rRNA.intervals'.format(reference_folder, referencePure)
    reference2 = referencePure + '.' + locus_function_list

    runinfo_file = '{}/RunInfo.xml'.format(flowcell_directory)
    log_file = '{}/logs/workflow.log'.format(output_folder)

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
                if tile_cou_per_slice * i >= tile_cou:
                    break
                slice_id[lane].append(str(i))
                slice_first_tile[lane].append(str(tile_nums[tile_cou_per_slice * i]))
                slice_tile_limit[lane].append(str(tile_cou_per_slice))
    
    alignment_folder = '{}/{}_{}/{}/alignment/'.format(library_folder, experiment_date, library, reference2)
    barcode_matching_folder = '{}/{}_{}/{}/barcode_matching/'.format(library_folder, experiment_date, library, reference2)
    select_cell_file = '{}{}.{}_transcripts_mq_{}_selected_cells.txt'.format(alignment_folder, library, min_transcripts_per_cell, base_quality)
    bead_barcode_file = '{}/BeadBarcodes.txt'.format(barcode_matching_folder)
    
    if not os.path.isfile(select_cell_file):
        write_log(log_file, flowcell_barcode, 'run_cmatcher_combine error: '+select_cell_file+' does not exist!')
        raise Exception('run_cmatcher_combine error: '+select_cell_file+' does not exist!')
        
    folder_running = '{}/status/running.cmatcher_combine_{}_{}'.format(output_folder, library, reference2)
    folder_finished = '{}/status/finished.cmatcher_combine_{}_{}'.format(output_folder, library, reference2)
    folder_failed = '{}/status/failed.cmatcher_combine_{}_{}'.format(output_folder, library, reference2)
    
    try:        
        call(['mkdir', folder_running])
        
        l = 0
        with open(select_cell_file, 'r') as fin:
            for line in fin:
                l += 1
        fin.close()
        k = 50000
        ls = l // k
        
        print('# selected cells: ' + str(l))
        
        while 1:
            f = True
            for i in range(ls + 1):
                if i * k >= l:
                    break;
                file2 = '{}/{}_barcode_matching_{}.finished'.format(barcode_matching_folder, library, str(i + 1))
                if not os.path.isfile(file2):
                    f = False
                    break
            if f:
                break
            time.sleep(30)

        print('combine cmatcher outputs...')
        write_log(log_file, flowcell_barcode, "Combine CMatcher outputs for "+library+" "+reference2)
        combined_cmatcher_file = '{}/{}_barcode_matching.txt'.format(barcode_matching_folder, library)
        with open(combined_cmatcher_file, 'w') as fout:
            fout.write('IlluminaBarcodes\tProcessedIlluminaBarcodes\tBeadBarcodes\tDistance\tX\tY\n')
            for i in range(ls + 1):
                if i * k >= l:
                    break;
                file2 = '{}/{}_barcode_matching_{}.txt'.format(barcode_matching_folder, library, str(i + 1))
                with open(file2, 'r') as fin:
                    j = 0
                    for line in fin:
                        j += 1
                        if j > 1:
                            fout.write(line)
                fin.close()
        fout.close()
        
        for i in range(ls + 1):
            if i * k >= l:
                break;
            file1 = '{}/{}_barcode_matching_{}.txt'.format(barcode_matching_folder, library, str(i + 1))
            file2 = '{}/{}_barcode_matching_{}.finished'.format(barcode_matching_folder, library, str(i + 1))
            name = '{}.{}_transcripts_mq_{}_selected_cells'.format(library, min_transcripts_per_cell, base_quality)
            file3 = '{}/{}_{}.txt'.format(alignment_folder, name, str(i + 1))
            if os.path.isfile(file1):
                call(['rm', file1])
            if os.path.isfile(file2):
                call(['rm', file2])
            if os.path.isfile(file3):
                call(['rm', file3])

        combined_cmatcher_file2 = '{}/{}_barcode_matching_distance.txt'.format(barcode_matching_folder, library)
        with open(combined_cmatcher_file2, 'w') as fout:
            fout.write('IlluminaBarcodes\tProcessedIlluminaBarcodes\tBeadBarcodes\tDistance\n')
            for i in range(ls + 1):
                if i * k >= l:
                    break;
                file2 = '{}/{}_barcode_matching_distance_{}.txt'.format(barcode_matching_folder, library, str(i + 1))
                with open(file2, 'r') as fin:
                    j = 0
                    for line in fin:
                        j += 1
                        if j > 1:
                            fout.write(line)
                fin.close()
        fout.close()
        
        write_log(log_file, flowcell_barcode, "Combine CMatcher outputs for "+library+" "+reference2+" is done. ")
        
        # Get unique matched bead barcodes and locations
        print('Get unique matched bead barcodes and locations...')
        write_log(log_file, flowcell_barcode, "Get unique matched bead barcodes and locations for "+library+" "+reference2)
        dict = {}
        matched_bead_barcode_file = '{}/{}_matched_bead_barcodes.txt'.format(barcode_matching_folder, library)
        matched_bead_location_file = '{}/{}_matched_bead_locations.txt'.format(barcode_matching_folder, library)
        with open(matched_bead_barcode_file, 'w') as fout1:
            with open(matched_bead_location_file, 'w') as fout2:
                with open(combined_cmatcher_file, 'r') as fin:
                    j = 0
                    for line in fin:
                        j += 1
                        if j > 1:
                            bc = line.split('\t')[2]
                            dist = line.split('\t')[3]
                            x = line.split('\t')[4]
                            y = line.split('\t')[5]
                            if not bc in dict:
                                fout1.write(bc + '\n')
                                fout2.write(dist + '\t' + x + '\t' + y)
                                dict[bc] = 1
                fin.close()
            fout2.close()
        fout1.close()
        
        print('Gzip unique matched bead barcode file...')
        matched_bead_barcode_gzfile = '{}/{}_matched_bead_barcodes.txt.gz'.format(barcode_matching_folder, library)
        os.system('gzip -c {} > {}'.format(matched_bead_barcode_file, matched_bead_barcode_gzfile))
        write_log(log_file, flowcell_barcode, "Get unique matched bead barcodes and locations for "+library+" "+reference2+" is done. ")

        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for slice in slice_id[lanes[i]]:
                # Call tag_matched_bam
                output_file = '{}/logs/tag_matched_bam_{}_{}_{}_{}_{}.log'.format(output_folder, library, lanes[i], slice, barcodes[i], reference2)
                submission_script = '{}/run.sh'.format(scripts_folder)
                call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=20g', '-notify', '-l', 'h_rt=5:0:0', '-j', 'y', submission_script, 'tag_matched_bam', manifest_file, library, lanes[i], slice, barcodes[i], locus_function_list, scripts_folder]
                call(call_args)
                
                # Call filter_unmapped_bam
                output_file = '{}/logs/filter_unmapped_bam_{}_{}_{}_{}_{}.log'.format(output_folder, library, lanes[i], slice, barcodes[i], reference2)
                submission_script = '{}/run.sh'.format(scripts_folder)
                call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=20g', '-notify', '-l', 'h_rt=5:0:0', '-j', 'y', submission_script, 'filter_unmapped_bam', manifest_file, library, lanes[i], slice, barcodes[i], locus_function_list, scripts_folder]
                call(call_args)
                
        # Call generate_plots_cmatcher
        output_file = '{}/logs/generate_plots_cmatcher_{}_{}.log'.format(output_folder, library, reference2)
        submission_script = '{}/run.sh'.format(scripts_folder)
        call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=30g', '-notify', '-l', 'h_rt=20:0:0', '-j', 'y', submission_script, 'generate_plots_cmatcher', manifest_file, library, scripts_folder, locus_function_list]
        call(call_args)
        
        call(['mv', folder_running, folder_finished])
    except:
        if os.path.isdir(folder_running):
            call(['mv', folder_running, folder_failed])
        else:
            call(['mkdir', folder_failed])
        
        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = "The Slide-seq workflow for "+library+" "+reference2+" failed at the step of running cmatcher combine. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)

        sys.exit()


if __name__ == "__main__":
    main()
    
    
