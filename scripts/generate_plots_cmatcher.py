#!/usr/bin/python

# This script is to generate analysis outputs and PDFs on matched barcodes

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


# Get read 1 length
def get_read1_len(bs):
    #12C8M|*T
    #7C18X7C8M2X|*T
    l = re.split('C|X|M', re.split('\|', bs)[0])
    i = 0
    for it in l:
        if it:
            i += int(it)
    return i
    

def hamming(seq1, seq2):
    c = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            c += 1
    return c


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
    
    analysis_folder = '{}/{}_{}/'.format(library_folder, experiment_date, library)
    alignment_folder = '{}/{}_{}/{}/alignment/'.format(library_folder, experiment_date, library, reference2)
    barcode_matching_folder = '{}/{}_{}/{}/barcode_matching/'.format(library_folder, experiment_date, library, reference2)
    
    bead_barcode_file = '{}/BeadBarcodes.txt'.format(barcode_matching_folder)

    folder_running = '{}/status/running.generate_plots_cmatcher_{}_{}'.format(output_folder, library, reference2)
    folder_finished = '{}/status/finished.generate_plots_cmatcher_{}_{}'.format(output_folder, library, reference2)
    folder_failed = '{}/status/failed.generate_plots_cmatcher_{}_{}'.format(output_folder, library, reference2)

    # Wait for all of tag_matched_bam finish
    while 1:
        f = True
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for slice in slice_id[lanes[i]]:
                fol1 = '{}/status/finished.tag_matched_bam_{}_{}_{}_{}'.format(output_folder, library, lanes[i], slice, reference2)
                fol2 = '{}/status/failed.tag_matched_bam_{}_{}_{}_{}'.format(output_folder, library, lanes[i], slice, reference2)
                if (not os.path.isdir(fol1)) and (not os.path.isdir(fol2)):
                    f = False
                    break
            if not f:
                break
        if f:
            break
        time.sleep(60)
    
    try:
        call(['mkdir', folder_running])
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
        print('Merge tagged matched bam files...')
        write_log(log_file, flowcell_barcode, "Merge tagged matched bam files for "+library+" "+reference2)
        
        matched_bam_file = '{}/{}_matched.bam'.format(barcode_matching_folder, library)
        commandStr = 'java -Djava.io.tmpdir={} -Dsamjdk.buffer_size=131072 -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8192m '.format(tmpdir)
        commandStr += '-jar {}/picard.jar MergeSamFiles TMP_DIR={} CREATE_INDEX=true CREATE_MD5_FILE=false VALIDATION_STRINGENCY=SILENT '.format(picard_folder, tmpdir)
        commandStr += 'OUTPUT={} SORT_ORDER=coordinate ASSUME_SORTED=true'.format(matched_bam_file)
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for slice in slice_id[lanes[i]]:
                filtered_bam = '{}/{}_{}_{}_tagged.bam'.format(barcode_matching_folder, library, lanes[i], slice)
                if os.path.isfile(filtered_bam):
                    commandStr += ' INPUT={}'.format(filtered_bam)
                else:
                    print(filtered_bam + ' not found!')
        write_log(log_file, flowcell_barcode, "MergeSamFiles for "+library+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "MergeSamFiles for "+library+" is done. ")
        
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for slice in slice_id[lanes[i]]:
                filtered_bam = '{}/{}_{}_{}_tagged.bam'.format(barcode_matching_folder, library, lanes[i], slice)
                if os.path.isfile(filtered_bam):
                    call(['rm', filtered_bam])
        
        write_log(log_file, flowcell_barcode, "Merge tagged matched bam files for "+library+" "+reference2+" is done. ")
        
        matched_bead_barcode_file = '{}/{}_matched_bead_barcodes.txt'.format(barcode_matching_folder, library)
        matched_bead_location_file = '{}/{}_matched_bead_locations.txt'.format(barcode_matching_folder, library)        
        matched_bead_barcode_gzfile = '{}/{}_matched_bead_barcodes.txt.gz'.format(barcode_matching_folder, library)

        # Generate digital expression files
        print('Generate digital expression files...')
        commandStr = '{}/DigitalExpression -m 7692m I={} O={}{}.digital_expression.txt.gz '.format(dropseq_folder, matched_bam_file, alignment_folder, library)
        commandStr += 'SUMMARY={}{}.digital_expression_summary.txt EDIT_DISTANCE=1 READ_MQ={} MIN_BC_READ_THRESHOLD=0 '.format(alignment_folder, library, base_quality)
        commandStr += 'CELL_BC_FILE={} TMP_DIR={} '.format(matched_bead_barcode_gzfile, tmpdir)
        commandStr += 'OUTPUT_HEADER=true UEI={} VALIDATION_STRINGENCY=SILENT'.format(library)
        if locus_function_list == 'exonic+intronic':
            commandStr += ' LOCUS_FUNCTION_LIST=INTRONIC'
        elif locus_function_list == 'intronic':
            commandStr += ' LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC'       
        write_log(log_file, flowcell_barcode, "DigitalExpression for "+library+" on matched barcodes Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "DigitalExpression for "+library+" on matched barcodes is done. ")

        # Bam tag histogram
        print('Bam tag histogram...')
        commandStr = '{}/BamTagHistogram -m 7692m '.format(dropseq_folder)
        commandStr += 'I={} OUTPUT={}{}.numReads_perCell_XC_mq_{}.txt.gz '.format(matched_bam_file, alignment_folder, library, base_quality)
        commandStr += 'TAG=XC FILTER_PCR_DUPLICATES=false TMP_DIR={} READ_MQ={} VALIDATION_STRINGENCY=SILENT'.format(tmpdir, base_quality)
        write_log(log_file, flowcell_barcode, "BamTagHistogram for "+library+" on matched barcodes Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "BamTagHistogram for "+library+" on matched barcodes is done. ")
        
        # Collect RnaSeq metrics
        print('Collect RnaSeq metrics...')
        commandStr = 'java -Djava.io.tmpdir={} -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8192m '.format(tmpdir)
        commandStr += '-jar {}/picard.jar CollectRnaSeqMetrics TMP_DIR={}'.format(picard_folder, tmpdir)
        commandStr += ' VALIDATION_STRINGENCY=SILENT I={} REF_FLAT={} STRAND_SPECIFICITY=NONE '.format(matched_bam_file, ref_flat)
        commandStr += 'OUTPUT={}{}.fracIntronicExonic.txt RIBOSOMAL_INTERVALS={}'.format(alignment_folder, library, ribosomal_intervals)
        write_log(log_file, flowcell_barcode, "CollectRnaSeqMetrics for "+library+" on matched barcodes Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "CollectRnaSeqMetrics for "+library+" on matched barcodes is done. ")
        
        # Base distribution at read position for cellular
        print('Base distribution at read position for cellular...')
        commandStr = '{}/BaseDistributionAtReadPosition -m 7692m I={}'.format(dropseq_folder, matched_bam_file)
        commandStr += ' OUTPUT={}{}.barcode_distribution_XC.txt TMP_DIR={} TAG=XC VALIDATION_STRINGENCY=SILENT'.format(alignment_folder, library, tmpdir)
        write_log(log_file, flowcell_barcode, "BaseDistributionAtReadPosition Cellular for "+library+" on matched barcodes Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "BaseDistributionAtReadPosition Cellular for "+library+" on matched barcodes is done. ")
        
        # Base distribution at read position for molecular
        print('Base distribution at read position for molecular...')
        commandStr = '{}/BaseDistributionAtReadPosition -m 7692m I={}'.format(dropseq_folder, matched_bam_file)
        commandStr += ' OUTPUT={}{}.barcode_distribution_XM.txt TMP_DIR={} TAG=XM VALIDATION_STRINGENCY=SILENT'.format(alignment_folder, library, tmpdir)
        write_log(log_file, flowcell_barcode, "BaseDistributionAtReadPosition Molecular for "+library+" on matched barcodes Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "BaseDistributionAtReadPosition Molecular for "+library+" on matched barcodes is done. ")

        # Gather read quality metrics
        print('Gather read quality metrics...')
        commandStr = '{}/GatherReadQualityMetrics -m 7692m I={}'.format(dropseq_folder, matched_bam_file)
        commandStr += ' TMP_DIR={} OUTPUT={}{}.ReadQualityMetrics.txt VALIDATION_STRINGENCY=SILENT'.format(tmpdir, alignment_folder, library)
        write_log(log_file, flowcell_barcode, "GatherReadQualityMetrics for "+library+" on matched barcodes Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "GatherReadQualityMetrics for "+library+" on matched barcodes is done. ")

        # Single cell RnaSeq metrics collector
        print('Single cell RnaSeq metrics collector...')
        commandStr = '{}/SingleCellRnaSeqMetricsCollector -m 15884m I={} ANNOTATIONS_FILE={}'.format(dropseq_folder, matched_bam_file, annotations_file)
        commandStr += ' OUTPUT={}{}.fracIntronicExonicPerCell.txt.gz RIBOSOMAL_INTERVALS={} CELL_BARCODE_TAG=XC READ_MQ={}'.format(alignment_folder, library, ribosomal_intervals, base_quality)
        commandStr += ' TMP_DIR={} CELL_BC_FILE={} MT_SEQUENCE=MT VALIDATION_STRINGENCY=SILENT'.format(tmpdir, matched_bead_barcode_gzfile)
        write_log(log_file, flowcell_barcode, "SingleCellRnaSeqMetricsCollector for "+library+" on matched barcodes Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "SingleCellRnaSeqMetricsCollector for "+library+" on matched barcodes is done. ")

        if os.path.isfile('{}/{}.SelectCellsByNumTranscripts_metrics'.format(alignment_folder, library)):
            call(['rm', '{}/{}.SelectCellsByNumTranscripts_metrics'.format(alignment_folder, library)])
        if os.path.isfile(matched_bead_barcode_gzfile):
            call(['rm', matched_bead_barcode_gzfile])

        print('generating plots... \n')
        write_log(log_file, flowcell_barcode, "Generate plots for matched barcodes for "+library+" "+reference2)
        
        matched_barcodes = np.loadtxt(matched_bead_barcode_file, delimiter='\t', dtype='str', usecols=(0))
        distances = np.loadtxt(matched_bead_location_file, delimiter='\t', dtype='int', usecols=(0))
        coordinatesx = np.loadtxt(matched_bead_location_file, delimiter='\t', dtype='float', usecols=(1))
        coordinatesy = np.loadtxt(matched_bead_location_file, delimiter='\t', dtype='float', usecols=(2))

        pp = PdfPages('{}/{}_{}.pdf'.format(alignment_folder, library, reference2))

        file = '{}/{}.ReadQualityMetrics.txt'.format(analysis_folder, library)
        print(file + '... \n')
        while 1:
            if os.path.isfile(file):
                break
            time.sleep(60)
        time.sleep(10)
        mat = np.loadtxt(file, delimiter='\t', dtype='int', skiprows=3, max_rows=1, usecols=(1,2,3,4))
        file1 = '{}/{}.ReadQualityMetrics.txt'.format(alignment_folder, library)
        mat1 = np.loadtxt(file1, delimiter='\t', dtype='int', skiprows=3, max_rows=1, usecols=(1,2,3,4))

        df_z = [mat[0],mat[1],mat[2],mat1[0]]
        if mat[2] >= 1000000 and mat1[0] >= 1000000:
            df_u = [int(mat[0]/1000000),int(mat[1]/1000000),int(mat[2]/1000000),int(mat1[0]/1000000)]
            yl = "# Reads [millions]"
        else:
            df_u = [mat[0],mat[1],mat[2],mat1[0]]
            yl = "# Reads"
        df_y = [mat[0]/mat[0]*100,mat[1]/mat[0]*100,mat[2]/mat[0]*100,mat1[0]/mat[0]*100]
        df_v = ['{:,}'.format(mat[0]),'{:,}'.format(mat[1]),'{:,}'.format(mat[2]),'{:,}'.format(mat1[0])]
        labels = []
        for i in range(4):
            labels.append('{0:.3g}%'.format(df_y[i]))
        df = pd.DataFrame({"x":['Total','Mapped','HQ','HQ matched'], "z":df_z, "u":df_u})
        fig, ax = plt.subplots(figsize=(8, 8))
        bp = plt.bar(df['x'], df['u'], width=0.7, color='lightskyblue', edgecolor='black')
        for idx,rect in enumerate(bp):
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 0.5*height, labels[idx], ha='center', va='bottom')
            ax.text(rect.get_x() + rect.get_width()/2., height, df_v[idx], ha='center', va='bottom')
        plt.yticks(rotation=90)
        plt.ylabel(yl)
        plt.title("Alignment quality for all reads")
        plt.savefig(pp, format='pdf')
        
        file = '{}/{}.fracIntronicExonic.txt'.format(analysis_folder, library)
        print(file + '... \n')
        while 1:
            if os.path.isfile(file):
                break
            time.sleep(60)
        time.sleep(10)
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
        plt.savefig(pp, format='pdf')
        
        f1 = '{}/{}.numReads_perCell_XC_mq_{}.txt'.format(analysis_folder, library, base_quality)
        f2 = '{}/{}.numReads_perCell_XC_mq_{}.txt.gz'.format(analysis_folder, library, base_quality)
        print(f2 + '... \n')
        while 1:
            if os.path.isfile(f1) or os.path.isfile(f2):
                break
            time.sleep(60)
        time.sleep(10)
        if not os.path.isfile(f1):
            with gzip.open(f2, 'rb') as f_in:
                with open(f1, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                f_out.close()
            f_in.close()
        
        print(f1 + '... \n')
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
        plt.savefig(pp, format='pdf')
        
        if os.path.isfile(f1):
            call(['rm', f1])
        
        file1 = '{}/{}.ReadQualityMetrics.txt'.format(analysis_folder, library)
        print(file1 + '... \n')
        cou1 = np.loadtxt(file1, delimiter='\t', dtype='int', skiprows=3, max_rows=1, usecols=(1))
        l = 42
        f = False
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for slice in slice_id[lanes[i]]:
                file = '{}/{}.{}.{}.{}.{}.polyA_trimming_report.txt'.format(analysis_folder, flowcell_barcode, lanes[i], slice, library, barcodes[i])
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
                file = '{}/{}.{}.{}.{}.{}.polyA_trimming_report.txt'.format(analysis_folder, flowcell_barcode, lanes[i], slice, library, barcodes[i])
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
        plt.savefig(pp, format='pdf')
        
        file = '{}/{}.barcode_distribution_XC.txt'.format(analysis_folder, library)
        print(file + '... \n')
        while 1:
            if os.path.isfile(file):
                break
            time.sleep(60)
        time.sleep(10)
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
        plt.savefig(pp, format='pdf')
        
        file = '{}/{}.barcode_distribution_XM.txt'.format(analysis_folder, library)
        print(file + '... \n')
        while 1:
            if os.path.isfile(file):
                break
            time.sleep(60)
        time.sleep(10)
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
        plt.savefig(pp, format='pdf')
        
        num_beads = 0
        if os.path.isfile(bead_barcode_file):
            with open(bead_barcode_file, 'r') as fin:
                for line in fin:
                    num_beads += 1
            fin.close()

        num_cells = 0
        selected_cells = '{}/{}.{}_transcripts_mq_{}_selected_cells.txt'.format(alignment_folder, library, min_transcripts_per_cell, base_quality)
        if os.path.isfile(selected_cells):
            with open(selected_cells, 'r') as fin:
                for line in fin:
                    num_cells += 1
            fin.close()
                
        num_matched_beads = 0
        dge_summary = '{}/{}.digital_expression_summary.txt'.format(alignment_folder, library)
        if os.path.isfile(dge_summary):
            dge_summary_reads = np.loadtxt(dge_summary, delimiter='\t', dtype='int', skiprows=7, usecols=(1))
            num_matched_beads = len(dge_summary_reads)

        num_matched_cells = -1
        combined_cmatcher_file = '{}/{}_barcode_matching.txt'.format(barcode_matching_folder, library)
        if os.path.isfile(combined_cmatcher_file):
            with open(combined_cmatcher_file, 'r') as fin:
                for line in fin:
                    num_matched_cells += 1
            fin.close()

        xtitle = ''
        if num_cells > 0 and num_beads > 0:
            ratio_cell = '{0:.3g}%'.format(num_matched_cells / num_cells * 100)
            ratio_bead = '{0:.3g}%'.format(num_matched_beads / num_beads * 100)
            xtitle = '{} Illumina barcodes matching to bead barcodes\n{} bead barcodes matching to Illumina barcodes'.format(ratio_cell, ratio_bead)
        
        file1 = '{}/{}.ReadQualityMetrics.txt'.format(alignment_folder, library)
        print(file1 + '... \n')
        mat1 = np.loadtxt(file1, delimiter='\t', dtype='int', skiprows=3, max_rows=1, usecols=(1,2,3,4))
        df_z = [mat1[0],mat1[1],mat1[2],mat1[3]]
        if mat1[3] >= 1000000:
            df_u = [mat1[0]/1000000,mat1[1]/1000000,mat1[2]/1000000,mat1[3]/1000000]
            yl = "# Reads [millions]"
        else:
            df_u = [mat1[0],mat1[1],mat1[2],mat1[3]]
            yl = "# Reads"
        df_y = [mat1[0]/mat1[0]*100,mat1[1]/mat1[0]*100,mat1[2]/mat1[0]*100,mat1[3]/mat1[0]*100]
        df_v = ['{:,}'.format(mat1[0]),'{:,}'.format(mat1[1]),'{:,}'.format(mat1[2]),'{:,}'.format(mat1[3])]
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
        plt.xlabel(xtitle)
        plt.ylabel(yl)
        plt.title("Alignment quality for matched barcodes")
        plt.savefig(pp, format='pdf')
        
        file = '{}/{}.fracIntronicExonic.txt'.format(alignment_folder, library)
        print(file + '... \n')
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
        plt.title("All reads that mapped to matched barcodes")
        plt.savefig(pp, format='pdf')

        df_y = [0,0,0]
        for i in range(len(distances)):
            df_y[distances[i]] += 1
        df = pd.DataFrame({"x":[0,1,2],"y":df_y})
        plt.figure(figsize=(8, 8))
        plt.bar(df['x'], df['y'], width=1, color='lightskyblue', edgecolor='black')
        plt.xlim((-1, 15))
        plt.yticks(rotation=90)
        plt.xlabel("hamming distance\nThreshold for choosing matched barcodes: hamming distance <= 1")
        plt.ylabel("number of matched barcodes")
        plt.title("Histogram of barcode matches by hamming distance")
        plt.savefig(pp, format='pdf')

        matched_dge_summary = '{}/{}.digital_expression_summary.txt'.format(alignment_folder, library)
        print(matched_dge_summary + '... \n')
        matched_dge_summary_barcodes = np.loadtxt(matched_dge_summary, delimiter='\t', dtype='str', skiprows=7, usecols=(0))
        matched_dge_summary_reads = np.loadtxt(matched_dge_summary, delimiter='\t', dtype='int', skiprows=7, usecols=(1))
        matched_distances = [0,0,0]
        for i in range(len(matched_dge_summary_barcodes)):
            if matched_dge_summary_barcodes[i] in matched_barcodes:
                matched_distances[int(distances[np.where(matched_barcodes==matched_dge_summary_barcodes[i])])] += matched_dge_summary_reads[i]
        for i in range(3):
            matched_distances[i] = matched_distances[i] / 1000000
        df = pd.DataFrame({"x":[0,1,2],"y":matched_distances})
        plt.figure(figsize=(8, 8))
        plt.bar(df['x'], df['y'], width=1, color='lightskyblue', edgecolor='black')
        plt.xlim((-1, 15))
        plt.yticks(rotation=90)
        plt.xlabel("hamming distance\nThreshold for choosing matched barcodes: hamming distance <= 1")
        plt.ylabel("# Reads [millions]")
        plt.title("Histogram of reads matches by hamming distance")
        plt.savefig(pp, format='pdf')

        matched_dge_summary_transcripts = np.loadtxt(matched_dge_summary, delimiter='\t', dtype='float', skiprows=7, usecols=(2))
        for i in range(0, len(matched_dge_summary_transcripts)):
            matched_dge_summary_transcripts[i] = math.log10(matched_dge_summary_transcripts[i])
        df = pd.DataFrame({"x":matched_dge_summary_transcripts})
        plt.figure(figsize=(8, 8))
        plt.hist(df['x'], bins=70, facecolor='lightskyblue', edgecolor='black')
        plt.yticks(rotation=90)
        plt.xlabel("log10-based number of UMIs")
        plt.title("Histogram of UMIs per matched barcode")
        plt.savefig(pp, format='pdf')

        matched_dge_summary_genes = np.loadtxt(matched_dge_summary, delimiter='\t', dtype='float', skiprows=7, usecols=(3))
        for i in range(0, len(matched_dge_summary_genes)):
            matched_dge_summary_genes[i] = math.log10(matched_dge_summary_genes[i])
        df = pd.DataFrame({"x":matched_dge_summary_genes})
        plt.figure(figsize=(8, 8))
        plt.hist(df['x'], bins=70, facecolor='lightskyblue', edgecolor='black')
        plt.yticks(rotation=90)
        plt.xlabel("log10-based number of genes")
        plt.title("Histogram of genes per matched barcode")
        plt.savefig(pp, format='pdf')

        file = '{}/{}.barcode_distribution_XC.txt'.format(alignment_folder, library)
        print(file + '... \n')
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
        plt.title("Matched cell barcodes")
        plt.savefig(pp, format='pdf')

        file = '{}/{}.barcode_distribution_XM.txt'.format(alignment_folder, library)
        print(file + '... \n')
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
        plt.title("Matched molecular barcodes")
        plt.savefig(pp, format='pdf')
        
        # Wait for all of filter_unmapped_bam finish
        while 1:
            f = True
            for i in range(len(lanes)):
                if libraries[i] != library:
                    continue
                for slice in slice_id[lanes[i]]:
                    fol1 = '{}/status/finished.filter_unmapped_bam_{}_{}_{}_{}'.format(output_folder, library, lanes[i], slice, reference2)
                    fol2 = '{}/status/failed.filter_unmapped_bam_{}_{}_{}_{}'.format(output_folder, library, lanes[i], slice, reference2)
                    if (not os.path.isdir(fol1)) and (not os.path.isdir(fol2)):
                        f = False
                        break
                if not f:
                    break
            if f:
                break
            time.sleep(60)

        print('Merge filtered unmapped bam files...')
        write_log(log_file, flowcell_barcode, "Merge filtered unmapped bam files for "+library+" "+reference2)
        
        unmapped_bam_file = '{}/{}_unmapped.bam'.format(barcode_matching_folder, library)
        read1_file = '{}/{}.base_distribution_read1.txt'.format(alignment_folder, library)
        commandStr = 'java -Djava.io.tmpdir={} -Dsamjdk.buffer_size=131072 -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8192m '.format(tmpdir)
        commandStr += '-jar {}/picard.jar MergeSamFiles TMP_DIR={} CREATE_INDEX=true CREATE_MD5_FILE=false VALIDATION_STRINGENCY=SILENT '.format(picard_folder, tmpdir)
        commandStr += 'OUTPUT={} SORT_ORDER=coordinate ASSUME_SORTED=true'.format(unmapped_bam_file)
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for slice in slice_id[lanes[i]]:
                filtered_bam = '{}/{}_{}_{}_filtered.bam'.format(barcode_matching_folder, library, lanes[i], slice)
                if os.path.isfile(filtered_bam):
                    commandStr += ' INPUT={}'.format(filtered_bam)
                else:
                    print(filtered_bam + ' not found!')
        write_log(log_file, flowcell_barcode, "MergeSamFiles for "+library+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "MergeSamFiles for "+library+" is done. ")
        
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for slice in slice_id[lanes[i]]:
                filtered_bam = '{}/{}_{}_{}_filtered.bam'.format(barcode_matching_folder, library, lanes[i], slice)
                if os.path.isfile(filtered_bam):
                    call(['rm', filtered_bam])
        
        write_log(log_file, flowcell_barcode, "Merge filtered unmapped bam files for "+library+" "+reference2+" is done. ")
        
        # Base distribution at read position for read 1
        commandStr = '{}/BaseDistributionAtReadPosition -m 7692m I={}'.format(dropseq_folder, unmapped_bam_file)
        commandStr += ' OUTPUT={} TMP_DIR={} READ_NUMBER=1 VALIDATION_STRINGENCY=SILENT'.format(read1_file, tmpdir)
        os.system(commandStr)
        
        if os.path.isfile(unmapped_bam_file):
            call(['rm', unmapped_bam_file])
        if os.path.isfile('{}/{}_unmapped.bai'.format(barcode_matching_folder, library)):
            call(['rm', '{}/{}_unmapped.bai'.format(barcode_matching_folder, library)])

        read1_len = get_read1_len(bead_structure) # 42
        mat = np.loadtxt(read1_file, delimiter='\t', dtype='int', skiprows=1)
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
        fig, ax = plt.subplots(figsize=(8, 8))
        colors = ['red', 'blue', 'green', 'purple']
        labels = ['A', 'C', 'G', 'T']
        for i in range(4):
            ax.scatter(df_x[i*read1_len:(i+1)*read1_len], df_y[i*read1_len:(i+1)*read1_len], c=colors[i], s=20, label=labels[i])
        ax.legend(loc="lower right")                 
        plt.xlim((0, read1_len+2))
        plt.ylim((0, max_y+2))
        plt.yticks(rotation=90)
        plt.xlabel("base position")
        plt.ylabel("fraction of reads")
        plt.title("Read 1 for matched barcodes")
        plt.savefig(pp, format='pdf')

        xs = []
        ys = []
        zs = []
        matched_dge_summary_transcripts = np.loadtxt(matched_dge_summary, delimiter='\t', dtype='float', skiprows=7, usecols=(2))
        for i in range(len(matched_dge_summary_barcodes)):
            barcode = matched_dge_summary_barcodes[i]
            if barcode in matched_barcodes:
                zs.append(math.log10(matched_dge_summary_transcripts[i]))
                xs.append(coordinatesx[np.where(matched_barcodes==barcode)])
                ys.append(coordinatesy[np.where(matched_barcodes==barcode)])
        df = pd.DataFrame({"x":xs, "y":ys, "z":zs})
        plt.figure(figsize=(8, 8))
        plt.set_cmap('viridis_r')
        plt.scatter(df['x'], df['y'], c=df['z'], s=0.5)
        plt.colorbar()
        plt.axis('equal')
        plt.yticks(rotation=90)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title("log10 based total # UMIs per matched bead")
        plt.savefig(pp, format='pdf')

        f1 = '{}/{}.fracIntronicExonicPerCell.txt'.format(alignment_folder, library)
        f2 = '{}/{}.fracIntronicExonicPerCell.txt.gz'.format(alignment_folder, library)
        if not os.path.isfile(f1):
            with gzip.open(f2, 'rb') as f_in:
                with open(f1, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                f_out.close()
            f_in.close()
        
        #PCT_RIBOSOMAL_BASES	PCT_CODING_BASES	PCT_UTR_BASES	PCT_INTRONIC_BASES	PCT_INTERGENIC_BASES
        #17	18	19	20	21
        #mat[1] += mat[2]
        #mat[2] = mat[1] + mat[3]
        #df_x = ['ribosomal','exonic','genic','intronic','intergenic']
        pct_mt = np.loadtxt(f1, delimiter='\t', dtype='float', skiprows=7, usecols=(1))
        pct_ribosomal = np.loadtxt(f1, delimiter='\t', dtype='float', skiprows=7, usecols=(17))
        pct_coding = np.loadtxt(f1, delimiter='\t', dtype='float', skiprows=7, usecols=(18))
        pct_utr = np.loadtxt(f1, delimiter='\t', dtype='float', skiprows=7, usecols=(19))
        barcodes = np.loadtxt(f1, delimiter='\t', dtype='str', skiprows=7, usecols=(29))
        xs = []
        ys = []
        zs = []
        us = []
        vs = []
        for i in range(len(barcodes)):
            barcode = barcodes[i]
            if barcode in matched_barcodes:
                xs.append(coordinatesx[np.where(matched_barcodes==barcode)])
                ys.append(coordinatesy[np.where(matched_barcodes==barcode)])
                zs.append(pct_mt[i])
                us.append(pct_coding[i]+pct_utr[i])
                vs.append(pct_ribosomal[i])
        
        df = pd.DataFrame({"x":xs, "y":ys, "z":zs})
        plt.figure(figsize=(8, 8))
        plt.set_cmap('viridis_r')
        plt.scatter(df['x'], df['y'], c=df['z'], s=0.5)
        plt.colorbar()
        plt.axis('equal')
        plt.yticks(rotation=90)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title("% mitochondrial reads per matched bead")
        plt.savefig(pp, format='pdf')

        df = pd.DataFrame({"x":xs, "y":ys, "z":us})
        plt.figure(figsize=(8, 8))
        plt.set_cmap('viridis_r')
        plt.scatter(df['x'], df['y'], c=df['z'], s=0.5)
        plt.colorbar()
        plt.axis('equal')
        plt.yticks(rotation=90)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title("% exonic reads per matched bead")
        plt.savefig(pp, format='pdf')
        
        df = pd.DataFrame({"x":xs, "y":ys, "z":vs})
        plt.figure(figsize=(8, 8))
        plt.set_cmap('viridis_r')
        plt.scatter(df['x'], df['y'], c=df['z'], s=0.5)
        plt.colorbar()
        plt.axis('equal')
        plt.yticks(rotation=90)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title("% ribosomal reads per matched bead")
        plt.savefig(pp, format='pdf')
        
        pp.close()
        write_log(log_file, flowcell_barcode, "Generate plots for matched barcodes for "+library+" "+reference2+" is done. ")
        
        if os.path.isfile(f1):
            call(['rm', f1])
        
        if len(email_address) > 1:
            subject = "Slide-seq workflow finished for " + flowcell_barcode
            content = "The Slide-seq workflow for "+library+"_"+reference2+" is finished. Please check the output folder for the results. Thank you for using the Slide-seq tools! "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)
        
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
            content = "The Slide-seq workflow for "+library+" "+reference2+" failed at the step of generating plots for matched barcodes. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)

        sys.exit()
    

if __name__ == "__main__":
    main()
    
    
