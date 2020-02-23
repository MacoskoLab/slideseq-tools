#!/usr/bin/python

# This script is to generate digital expression and other analysis outputs

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
    
    basecalls_dir = '{}/Data/Intensities/BaseCalls'.format(flowcell_directory)
    
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
    sequence = 'AAGCAGTGGTATCAACGCAGAGTGAATGGG'
    base_quality = '10'
    min_transcripts_per_cell = '10'
    experiment_date = ''
    run_barcodematching = False
    bead_barcode_file = ''
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
                sequence = row[row0.index('start_sequence')]
                base_quality = row[row0.index('base_quality')]
                min_transcripts_per_cell = row[row0.index('min_transcripts_per_cell')]
                run_barcodematching = str2bool(row[row0.index('run_barcodematching')])
                bead_barcode_file = row[row0.index('bead_barcode_file')]
                experiment_date = row[row0.index('experiment_date')]
    fin.close()
    
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
                if (tile_cou_per_slice * i >= tile_cou):
                    break
                slice_id[lane].append(str(i))
                slice_first_tile[lane].append(str(tile_nums[tile_cou_per_slice * i]))
                slice_tile_limit[lane].append(str(tile_cou_per_slice))
	
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
    
    folder_running = '{}/status/running.analysis_spec_{}_{}'.format(output_folder, library, locus_function_list)
    folder_finished = '{}/status/finished.analysis_spec_{}_{}'.format(output_folder, library, locus_function_list)
    folder_failed = '{}/status/failed.analysis_spec_{}_{}'.format(output_folder, library, locus_function_list)
    
    alignment_folder = '{}/{}_{}/{}/alignment/'.format(library_folder, experiment_date, library, reference2)
    combined_bamfile = '{}/{}_{}/{}.bam'.format(library_folder, experiment_date, library, library)
    barcode_matching_folder = '{}/{}_{}/{}/barcode_matching/'.format(library_folder, experiment_date, library, reference2)
    
    call(['mkdir', folder_running])

    try:
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
    
        # Select cells by num transcripts
        commandStr = '{}/SelectCellsByNumTranscripts '.format(dropseq_folder)
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += '-m 24076m I={} MIN_TRANSCRIPTS_PER_CELL={} READ_MQ={}'.format(combined_bamfile, min_transcripts_per_cell, base_quality)
        else:
            commandStr += '-m 7692m I={} MIN_TRANSCRIPTS_PER_CELL={} READ_MQ={}'.format(combined_bamfile, min_transcripts_per_cell, base_quality)
        commandStr += ' OUTPUT={}{}.{}_transcripts_mq_{}_selected_cells.txt.gz '.format(alignment_folder, library, min_transcripts_per_cell, base_quality)
        commandStr += 'TMP_DIR={} VALIDATION_STRINGENCY=SILENT'.format(tmpdir)
        if locus_function_list == 'exonic+intronic':
            commandStr += ' LOCUS_FUNCTION_LIST=INTRONIC'
        elif locus_function_list == 'intronic':
            commandStr += ' LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC'
        write_log(log_file, flowcell_barcode, "SelectCellsByNumTranscripts for "+library+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "SelectCellsByNumTranscripts for "+library+" is done. ")

        # Call run_cmatcher
        if run_barcodematching:
            if os.path.isfile(bead_barcode_file):               
                name = '{}.{}_transcripts_mq_{}_selected_cells'.format(library, min_transcripts_per_cell, base_quality)
                select_cell_file = '{}{}.txt'.format(alignment_folder, name)
                select_cell_gzfile = '{}.gz'.format(select_cell_file)
                os.system('gunzip -c {} > {}'.format(select_cell_gzfile, select_cell_file))

                l = 0
                with open(select_cell_file, 'r') as fin:
                    for line in fin:
                        l += 1
                fin.close()
                k = 50000
                ls = l // k

                for i in range(ls + 1):
                    if i * k >= l:
                        break;
                    
                    infile2 = '{}/{}_{}.txt'.format(alignment_folder, name, str(i + 1))
                    commandStr = 'awk \'NR >= {} && NR <= {}\' {} > {}'.format(str(i * k + 1), str((i+1) * k), select_cell_file, infile2)
                    os.system(commandStr)
                    
                    file4 = '{}/{}_barcode_matching_distance_{}.txt'.format(barcode_matching_folder, library, str(i + 1))
                    file5 = '{}/{}_barcode_matching_{}.txt'.format(barcode_matching_folder, library, str(i + 1))
                    output_file = '{}/logs/run_cmatcher_{}_{}_{}.log'.format(output_folder, library, locus_function_list, str(i + 1))
                    submission_script = '{}/run.sh'.format(scripts_folder)
                    call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=10g', '-notify', '-l', 'h_rt=5:0:0', '-j', 'y', submission_script, 'cmatcher', scripts_folder, bead_barcode_file, infile2, file4, file5, bead_type, '1']
                    call(call_args)
                    write_log(log_file, flowcell_barcode, "Run CMatcher for "+library+" "+reference2+" "+str(i + 1))

                # Call run_cmatcher_combine
                output_file = '{}/logs/run_cmatcher_combine_{}_{}.log'.format(output_folder, library, locus_function_list)
                submission_script = '{}/run.sh'.format(scripts_folder)
                call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=20g', '-notify', '-l', 'h_rt=12:0:0', '-j', 'y', submission_script, 'run_cmatcher_combine', manifest_file, library, scripts_folder, locus_function_list]
                call(call_args)
            else:
                run_barcodematching = False
                if os.path.isdir(barcode_matching_folder):
                    call(['rm', '-r', barcode_matching_folder])
                print("File {} does not exist. Do not run barcode matching...".format(bead_barcode_file))
        
        # Generate digital expression files for all Illumina barcodes
        commandStr = '{}/DigitalExpression '.format(dropseq_folder)
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += '-m 32268m '
        else:
            commandStr += '-m 7692m '
        commandStr += 'I={} O={}{}.AllIllumina.digital_expression.txt.gz '.format(combined_bamfile, alignment_folder, library)
        commandStr += 'SUMMARY={}{}.AllIllumina.digital_expression_summary.txt EDIT_DISTANCE=1 READ_MQ={} MIN_BC_READ_THRESHOLD=0 '.format(alignment_folder, library, base_quality)
        commandStr += 'CELL_BC_FILE={}{}.{}_transcripts_mq_{}_selected_cells.txt.gz TMP_DIR={} '.format(alignment_folder, library, min_transcripts_per_cell, base_quality, tmpdir)
        commandStr += 'OUTPUT_HEADER=false UEI={} VALIDATION_STRINGENCY=SILENT'.format(library)
        if locus_function_list == 'exonic+intronic':
            commandStr += ' LOCUS_FUNCTION_LIST=INTRONIC'
        elif locus_function_list == 'intronic':
            commandStr += ' LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC'
        write_log(log_file, flowcell_barcode, "DigitalExpression for "+library+" for all Illumina barcodes Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "DigitalExpression for "+library+" for all Illumina barcodes is done. ")
        
        if not run_barcodematching:
            if os.path.isdir(barcode_matching_folder):
                call(['rm', '-r', barcode_matching_folder])
            if len(email_address) > 1:
                subject = "Slide-seq workflow finished for " + flowcell_barcode
                content = "The Slide-seq workflow for "+library+"_"+locus_function_list+" is finished. Please check the output folder for the results. Thank you for using the Slide-seq tools! "
                call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
                call(call_args)
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
        call(['mv', folder_running, folder_finished])
    except:
        if os.path.isdir(folder_running):
            call(['mv', folder_running, folder_failed])
        elif os.path.isdir(folder_waiting):
            call(['mv', folder_waiting, folder_failed])
        else:
            call(['mkdir', folder_failed])
            
        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = "The Slide-seq workflow for "+library+" "+locus_function_list+" failed at the step of running specific analysis. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)

        sys.exit()
    

if __name__ == "__main__":
    main()

