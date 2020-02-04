#!/usr/bin/python

# This script is to generate analysis outputs from alignment and generate plots if not running barcode matching

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
    sequence = 'AAGCAGTGGTATCAACGCAGAGTGAATGGG'
    base_quality = '10'
    min_transcripts_per_cell = '10'
    bead_structure = ''
    experiment_date = ''
    run_barcodematching = False
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
                sequence = row[row0.index('start_sequence')]
                base_quality = row[row0.index('base_quality')]
                min_transcripts_per_cell = row[row0.index('min_transcripts_per_cell')]
                bead_structure = row[row0.index('bead_structure')]
                experiment_date = row[row0.index('experiment_date')]
                run_barcodematching = str2bool(row[row0.index('run_barcodematching')])
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
    
    folder_running = '{}/status/running.generate_plots_{}'.format(output_folder, library)
    folder_finished = '{}/status/finished.generate_plots_{}'.format(output_folder, library)
    folder_failed = '{}/status/failed.generate_plots_{}'.format(output_folder, library)
    
    alignment_folder = '{}/{}_{}/'.format(library_folder, experiment_date, library)
    combined_bamfile = '{}/{}.bam'.format(alignment_folder, library)

    try:
        call(['mkdir', folder_running])
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
        # Bam tag histogram
        commandStr = '{}/BamTagHistogram '.format(dropseq_folder)
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += '-m 15884m '
        else:
            commandStr += '-m 7692m '
        commandStr += 'I={} OUTPUT={}{}.numReads_perCell_XC_mq_{}.txt.gz '.format(combined_bamfile, alignment_folder, library, base_quality)
        commandStr += 'TAG=XC FILTER_PCR_DUPLICATES=false TMP_DIR={} READ_MQ={} VALIDATION_STRINGENCY=SILENT'.format(tmpdir, base_quality)
        write_log(log_file, flowcell_barcode, "BamTagHistogram for "+library+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "BamTagHistogram for "+library+" is done. ")
        
        # Collect RnaSeq metrics
        commandStr = 'java -Djava.io.tmpdir={} -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 '.format(tmpdir)
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += '-Xmx16384m '
        else:
            commandStr += '-Xmx8192m '
        commandStr += '-jar {}/picard.jar CollectRnaSeqMetrics TMP_DIR={}'.format(picard_folder, tmpdir)
        commandStr += ' VALIDATION_STRINGENCY=SILENT I={} REF_FLAT={} STRAND_SPECIFICITY=NONE '.format(combined_bamfile, ref_flat)
        commandStr += 'OUTPUT={}{}.fracIntronicExonic.txt RIBOSOMAL_INTERVALS={}'.format(alignment_folder, library, ribosomal_intervals)
        write_log(log_file, flowcell_barcode, "CollectRnaSeqMetrics for "+library+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "CollectRnaSeqMetrics for "+library+" is done. ")
        
        # Base distribution at read position for cellular
        commandStr = '{}/BaseDistributionAtReadPosition '.format(dropseq_folder)
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += '-m 15884m '
        else:
            commandStr += '-m 7692m '
        commandStr += 'I={} OUTPUT={}{}.barcode_distribution_XC.txt TMP_DIR={} TAG=XC VALIDATION_STRINGENCY=SILENT'.format(combined_bamfile, alignment_folder, library, tmpdir)
        write_log(log_file, flowcell_barcode, "BaseDistributionAtReadPosition Cellular for "+library+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "BaseDistributionAtReadPosition Cellular for "+library+" is done. ")
        
        # Base distribution at read position for molecular
        commandStr = '{}/BaseDistributionAtReadPosition '.format(dropseq_folder)
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += '-m 15884m '
        else:
            commandStr += '-m 7692m '
        commandStr += 'I={} OUTPUT={}{}.barcode_distribution_XM.txt TMP_DIR={} TAG=XM VALIDATION_STRINGENCY=SILENT'.format(combined_bamfile, alignment_folder, library, tmpdir)
        write_log(log_file, flowcell_barcode, "BaseDistributionAtReadPosition Molecular for "+library+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "BaseDistributionAtReadPosition Molecular for "+library+" is done. ")
        
        # Gather read quality metrics
        commandStr = '{}/GatherReadQualityMetrics '.format(dropseq_folder)
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += '-m 15884m '
        else:
            commandStr += '-m 7692m '
        commandStr += 'I={} TMP_DIR={} OUTPUT={}{}.ReadQualityMetrics.txt VALIDATION_STRINGENCY=SILENT'.format(combined_bamfile, tmpdir, alignment_folder, library)
        write_log(log_file, flowcell_barcode, "GatherReadQualityMetrics for "+library+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "GatherReadQualityMetrics for "+library+" is done. ")
        
        # Single cell RnaSeq metrics collector
        commandStr = '{}/SingleCellRnaSeqMetricsCollector '.format(dropseq_folder)
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += '-m 32268m '
        else:
            commandStr += '-m 15884m '
        commandStr += 'I={} ANNOTATIONS_FILE={}'.format(combined_bamfile, annotations_file)
        commandStr += ' OUTPUT={}{}.fracIntronicExonicPerCell.txt.gz RIBOSOMAL_INTERVALS={} CELL_BARCODE_TAG=XC READ_MQ={}'.format(alignment_folder, library, ribosomal_intervals, base_quality)
        commandStr += ' TMP_DIR={} CELL_BC_FILE={}{}.{}_transcripts_mq_{}_selected_cells.txt.gz MT_SEQUENCE=MT VALIDATION_STRINGENCY=SILENT'.format(tmpdir, alignment_folder, library, min_transcripts_per_cell, base_quality)
        #write_log(log_file, flowcell_barcode, "SingleCellRnaSeqMetricsCollector for "+library+" Command="+commandStr)
        #os.system(commandStr)
        #write_log(log_file, flowcell_barcode, "SingleCellRnaSeqMetricsCollector for "+library+" is done. ")
        
        if not run_barcodematching:
            pp1 = PdfPages('{}/{}.pdf'.format(alignment_folder, library))
            
            file = '{}/{}.ReadQualityMetrics.txt'.format(alignment_folder, library)
            mat = np.loadtxt(file, delimiter='\t', dtype='int', skiprows=3, max_rows=1, usecols=(1,2,3,4))
            totalreads = mat[0]
            df_z = [mat[0],mat[1],mat[2],mat[3]]
            if mat[3] >= 1000000:
                df_u = [int(mat[0]/1000000),int(mat[1]/1000000),int(mat[2]/1000000),int(mat[3]/1000000)]
                yl = "# Reads [millions]"
            else:
                df_u = [mat[0],mat[1],mat[2],mat[3]]
                yl = "# Reads"
            df_y = [mat[0]/mat[0]*100,mat[1]/mat[0]*100,mat[2]/mat[0]*100,mat[3]/mat[0]*100]
            df_v = ['{:,}'.format(mat[0]),'{:,}'.format(mat[1]),'{:,}'.format(mat[2]),'{:,}'.format(mat[3])]
            labels = []
            for i in range(4):
                labels.append('{0:.3g}%'.format(df_y[i]))
            df = pd.DataFrame({"x":['Total','Mapped','HQ','HQ No Dupes'], "z":df_z, "u":df_u})
            fig, ax = plt.subplots(figsize=(8, 8))
            bp = plt.bar(df['x'], df['u'], width=0.7, color='lightskyblue', edgecolor='black')
            for idx,rect in enumerate(bp):
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width()/2., 0.5*height, labels[idx], ha='center', va='bottom')
                ax.text(rect.get_x() + rect.get_width()/2., height, df_v[idx], ha='center', va='bottom')
            plt.yticks(rotation=90)
            plt.ylabel(yl)
            plt.title("Alignment quality for all reads")
            plt.savefig(pp1, format='pdf')
            
            file = '{}/{}.fracIntronicExonic.txt'.format(alignment_folder, library)
            mat = np.loadtxt(file, delimiter='\t', dtype='float', skiprows=7, max_rows=1, usecols=(15,16,17,18,19))
            mat[1] += mat[2]
            mat[2] = mat[1] + mat[3]
            df_x = ['ribosomal','exonic','genic','intronic','intergenic']
            df_y = mat * 100
            max_y = max(df_y)
            labels = []
            for i in range(5):
                labels.append('{0:.3g}'.format(df_y[i]))
            df = pd.DataFrame({"x":df_x, "y":df_y})
            fig, ax = plt.subplots(figsize=(8, 8))
            bp = plt.bar(df['x'], df['y'], width=0.7, color='lightskyblue', edgecolor='black')
            for idx,rect in enumerate(bp):
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width()/2., 0.5*height, labels[idx], ha='center', va='bottom')
            plt.yticks(rotation=90)
            plt.ylabel("Percentage")
            plt.title("All reads")
            plt.savefig(pp1, format='pdf')
            
            f1 = '{}/{}.numReads_perCell_XC_mq_{}.txt'.format(alignment_folder, library, base_quality)
            f2 = '{}/{}.numReads_perCell_XC_mq_{}.txt.gz'.format(alignment_folder, library, base_quality)
            if not os.path.isfile(f1):
                with gzip.open(f2, 'rb') as f_in:
                    with open(f1, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    f_out.close()
                f_in.close()

            mat = np.loadtxt(f1, delimiter='\t', dtype='int', skiprows=1, usecols=(0))
            l = int(len(mat) / 10)
            df_x = np.arange(1, l+1, 1)
            y = np.cumsum(mat)
            df_y = y / max(y)
            df = pd.DataFrame({"x":df_x, "y":df_y[:l]})
            plt.figure(figsize=(8, 8))
            plt.plot(df['x'], df['y'], color='green')
            plt.ylim((0, 1))
            plt.yticks(rotation=90)
            plt.xlabel("cell barcodes sorted by number of reads [descending]")
            plt.ylabel("cumulative fraction of reads")
            plt.title("Cumulative fraction of reads per cell barcode")
            plt.savefig(pp1, format='pdf')
            
            if os.path.isfile(f1):
                call(['rm', f1])
            
            file1 = '{}/{}.ReadQualityMetrics.txt'.format(alignment_folder, library)
            cou1 = np.loadtxt(file1, delimiter='\t', dtype='int', skiprows=3, max_rows=1, usecols=(1))
            l = 42
            f = False
            for i in range(len(lanes)):
                if libraries[i] != library:
                    continue
                for slice in slice_id[lanes[i]]:
                    file = '{}/{}.{}.{}.{}.{}.polyA_trimming_report.txt'.format(alignment_folder, flowcell_barcode, lanes[i], slice, library, barcodes[i])
                    if os.path.isfile(file):
                        mat = np.loadtxt(file, delimiter='\t', dtype='int', skiprows=7)
                        l = len(mat)
                        f = True
                        break
                if f:
                    break
            df_x = np.arange(0, l, 1)
            df_y = [0] * l
            for i in range(len(lanes)):
                if libraries[i] != library:
                    continue
                for slice in slice_id[lanes[i]]:
                    file = '{}/{}.{}.{}.{}.{}.polyA_trimming_report.txt'.format(alignment_folder, flowcell_barcode, lanes[i], slice, library, barcodes[i])
                    if os.path.isfile(file):
                        mat = np.loadtxt(file, delimiter='\t', dtype='int', skiprows=7)
                        for j in range(0, len(mat)):
                            df_y[mat[j,0]] += mat[j,1]
            max_y = max(df_y[1:])
            min_y = min(df_y[1:])
            cou = sum(df_y[1:])
            val = '{0:.3g}'.format(cou / cou1 * 100)
            df = pd.DataFrame({"x":df_x, "y":df_y})
            plt.figure(figsize=(8, 8))
            plt.plot(df['x'], df['y'], color='black')
            plt.xlim((1, l))
            plt.ylim((max(0,min_y-10000), max_y+10000))
            plt.yticks(rotation=90)
            plt.xlabel("first base of PolyA tail trimmed")
            plt.ylabel("number of reads")
            plt.title("% Reads trimmed by 3' PolyA trimmer: " + val)
            plt.savefig(pp1, format='pdf')
            
            file = '{}/{}.barcode_distribution_XC.txt'.format(alignment_folder, library)
            mat = np.loadtxt(file, delimiter='\t', dtype='int', skiprows=1)
            cou = mat[0,1] + mat[0,2] + mat[0,3] + mat[0,4]
            df_x = []
            df_x.extend(mat[:,0])
            df_x.extend(mat[:,0])
            df_x.extend(mat[:,0])
            df_x.extend(mat[:,0])
            df_y = []
            df_y.extend(mat[:,1] / cou * 100)
            df_y.extend(mat[:,2] / cou * 100)
            df_y.extend(mat[:,3] / cou * 100)
            df_y.extend(mat[:,4] / cou * 100)
            max_y = int(max(df_y))
            l = len(mat)
            fig, ax = plt.subplots(figsize=(8, 8))
            colors = ['red', 'blue', 'green', 'purple']
            labels = ['A', 'C', 'G', 'T']
            for i in range(4):
                ax.scatter(df_x[i*l:(i+1)*l], df_y[i*l:(i+1)*l], c=colors[i], s=20, label=labels[i])
            ax.legend(loc="lower right")                 
            plt.xlim((0, l+2))
            plt.ylim((0, max_y+2))
            plt.yticks(rotation=90)
            plt.xlabel("base position")
            plt.ylabel("fraction of reads")
            plt.title("Cell barcodes for all reads")
            plt.savefig(pp1, format='pdf')
            
            file = '{}/{}.barcode_distribution_XM.txt'.format(alignment_folder, library)
            mat = np.loadtxt(file, delimiter='\t', dtype='int', skiprows=1)
            cou = mat[0,1] + mat[0,2] + mat[0,3] + mat[0,4]
            df_x = []
            df_x.extend(mat[:,0])
            df_x.extend(mat[:,0])
            df_x.extend(mat[:,0])
            df_x.extend(mat[:,0])
            df_y = []
            df_y.extend(mat[:,1] / cou * 100)
            df_y.extend(mat[:,2] / cou * 100)
            df_y.extend(mat[:,3] / cou * 100)
            df_y.extend(mat[:,4] / cou * 100)
            max_y = int(max(df_y))
            l = len(mat)
            fig, ax = plt.subplots(figsize=(8, 8))
            colors = ['red', 'blue', 'green', 'purple']
            labels = ['A', 'C', 'G', 'T']
            for i in range(4):
                ax.scatter(df_x[i*l:(i+1)*l], df_y[i*l:(i+1)*l], c=colors[i], s=20, label=labels[i])
            ax.legend(loc="lower right")                 
            plt.xlim((0, l+1))
            plt.ylim((0, max_y+2))
            plt.yticks(rotation=90)
            plt.xlabel("base position")
            plt.ylabel("fraction of reads")
            plt.title("Molecular barcodes for all reads")
            plt.savefig(pp1, format='pdf')
            
            pp1.close()
        
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
            content = "The Slide-seq workflow for "+library+" failed at the step of generating plots. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)
        
        sys.exit()


if __name__ == "__main__":
    main()
    
    
