#!/usr/bin/python

# This script is to call the alignment step

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

from new_submit_to_taskrunner import call_to_taskrunner

import traceback

# Get read structure from RunInfo.xml
def get_read_structure(x):
    posts = 'TBT'
    i = 0
    str = ''
    with open(x, 'r') as fin:
        for line in fin:
            line = line.strip(' \t\n')
            if (line.startswith('<Read ', 0)):
                l = line.split('=')[2]
                l = l.split('\"')[1]
                str += l + posts[i]
                i += 1
                if (i >= 3):
                    break;
    fin.close()
    return str


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
    
    flowcell_directory = options['flowcell_directory']
    output_folder = options['output_folder']
    metadata_file = options['metadata_file']
    flowcell_barcode = options['flowcell_barcode']
    
    library_folder = options['library_folder'] if 'library_folder' in options else '{}/libraries'.format(output_folder)
    tmpdir = options['temp_folder'] if 'temp_folder' in options else '{}/tmp'.format(output_folder)
    dropseq_folder = options['dropseq_folder'] if 'dropseq_folder' in options else '/broad/macosko/bin/dropseq-tools'
    scripts_folder = options['scripts_folder'] if 'scripts_folder' in options else '/broad/macosko/jilong/slideseq_pipeline/scripts'
    is_NovaSeq = str2bool(options['is_NovaSeq']) if 'is_NovaSeq' in options else False
    is_NovaSeq_S4 = str2bool(options['is_NovaSeq_S4']) if 'is_NovaSeq_S4' in options else False
    num_slice_NovaSeq = int(options['num_slice_NovaSeq']) if 'num_slice_NovaSeq' in options else 10
    num_slice_NovaSeq_S4 = int(options['num_slice_NovaSeq_S4']) if 'num_slice_NovaSeq_S4' in options else 40
    email_address = options['email_address'] if 'email_address' in options else ''
    resubmit = str2bool(options['resubmit']) if 'resubmit' in options else False

    runinfo_file = '{}/RunInfo.xml'.format(flowcell_directory)
    log_file = '{}/logs/workflow.log'.format(output_folder)
    
    # Parse metadata file
    if resubmit:
        write_log(log_file, flowcell_barcode, "Parse metadata file. ")
        commandStr = 'python '+scripts_folder+'/parse_metadata.py -i '+metadata_file+' -r '+runinfo_file+' -o '+'{}/parsed_metadata.txt'.format(output_folder)
        os.system(commandStr)
    
    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    references_unique = []
    locus_function_list_unique = []
    resubmit_unique = []
    experiment_date = []
    run_barcodematching = []
    puckcaller_path = []
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
                references_unique.append(row[row0.index('reference')])
                locus_function_list_unique.append(row[row0.index('locus_function_list')])
                resubmit_unique.append(row[row0.index('resubmit')])
                experiment_date.append(row[row0.index('date')])
                run_barcodematching.append(str2bool(row[row0.index('run_barcodematching')]))
                puckcaller_path.append(row[row0.index('puckcaller_path')])
            barcodes.append(row[row0.index('sample_barcode')])
            bead_structures.append(row[row0.index('bead_structure')])
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
                if tile_cou_per_slice * i >= tile_cou:
                    break
                slice_id[lane].append(str(i))
                slice_first_tile[lane].append(str(tile_nums[tile_cou_per_slice * i]))
                slice_tile_limit[lane].append(str(tile_cou_per_slice))
    
    folder_waiting = '{}/status/waiting.mergebarcodes'.format(output_folder)
    folder_running = '{}/status/running.mergebarcodes'.format(output_folder)
    folder_finished = '{}/status/finished.mergebarcodes'.format(output_folder)
    folder_failed = '{}/status/failed.mergebarcodes'.format(output_folder)
    
    call(['mkdir', '-p', folder_waiting])
    
    # Wait till all of run_processbarcodes and run_barcodes2sam finish
    while 1:
        all_failed_1 = True
        all_failed_2 = True
        for lane in lanes_unique:
            failed_1 = '{}/status/failed.processbarcodes_lane_{}'.format(output_folder, lane)
            if not os.path.isdir(failed_1):
                all_failed_1 = False
            for i in range(len(slice_id[lane])):
                failed_2 = '{}/status/failed.barcodes2sam_lane_{}_{}'.format(output_folder, lane, slice_id[lane][i])
                if not os.path.isdir(failed_2):
                    all_failed_2 = False
                    break
        if all_failed_1:
            write_log(log_file, flowcell_barcode, "All run_processbarcodes failed. Exiting...")
            sys.exit()
        if all_failed_2:
            write_log(log_file, flowcell_barcode, "All run_barcodes2sam failed. Exiting...")
            sys.exit()
        
        f = True
        for lane in lanes_unique:
            for i in range(len(slice_id[lane])):
                fol1 = '{}/status/finished.barcodes2sam_lane_{}_{}'.format(output_folder, lane, slice_id[lane][i])
                fol2 = '{}/status/failed.barcodes2sam_lane_{}_{}'.format(output_folder, lane, slice_id[lane][i])
                if (not os.path.isdir(fol1)) and (not os.path.isdir(fol2)):
                    f = False
                    break
            if not f:
                break
        if f:
            break
        time.sleep(30)
    
    if os.path.isdir(folder_waiting):
        call(['mv', folder_waiting, folder_running])
    else:
        call(['mkdir', '-p', folder_running])
    
    try:
        for j in range(len(libraries_unique)):
            if (not resubmit) or resubmit_unique[j] == 'TRUE':
                library = libraries_unique[j]
                
                if resubmit:
                    for i in range(len(lanes)):
                        if libraries[i] != library:
                            continue
                        for slice in slice_id[lanes[i]]:
                            unmapped_bam = '{}/{}/{}/{}/'.format(output_folder, lanes[i], slice, library)
                            if (barcodes[i]):
                                unmapped_bam += '{}/{}.{}.{}.{}.{}.unmapped.bam'.format(barcodes[i], flowcell_barcode, lanes[i], slice, library, barcodes[i])
                            else:
                                unmapped_bam += '{}.{}.{}.{}.unmapped.bam'.format(flowcell_barcode, lanes[i], slice, library)
                            unmapped_bam2 = '{}/{}_{}/{}.{}.{}.{}'.format(library_folder, experiment_date[j], library, flowcell_barcode, lanes[i], slice, library)
                            if (barcodes[i]):
                                unmapped_bam2 += '.'+barcodes[i]
                            unmapped_bam2 += '.unmapped.bam'
                            if os.path.isfile(unmapped_bam2):
                                os.system('mv ' + unmapped_bam2 + ' ' + unmapped_bam)
                    if os.path.isdir('{}/{}_{}'.format(library_folder, experiment_date[j], library)):
                        os.system('rm -r ' + '{}/{}_{}'.format(library_folder, experiment_date[j], library))
                    os.system('rm ' + '{}/logs/*{}*'.format(output_folder, library))
                    os.system('rm -r ' + '{}/status/*{}*'.format(output_folder, library))
                
                os.system('mkdir -p ' + '{}/{}_{}'.format(library_folder, experiment_date[j], library))
                for i in range(len(lanes)):
                    if libraries[i] != library:
                        continue
                    for slice in slice_id[lanes[i]]:
                        unmapped_bam = '{}/{}/{}/{}/'.format(output_folder, lanes[i], slice, library)
                        if (barcodes[i]):
                            unmapped_bam += '{}/{}.{}.{}.{}.{}.unmapped.bam'.format(barcodes[i], flowcell_barcode, lanes[i], slice, library, barcodes[i])
                        else:
                            unmapped_bam += '{}.{}.{}.{}.unmapped.bam'.format(flowcell_barcode, lanes[i], slice, library)
                        unmapped_bam2 = '{}/{}_{}/{}.{}.{}.{}'.format(library_folder, experiment_date[j], library, flowcell_barcode, lanes[i], slice, library)
                        if (barcodes[i]):
                            unmapped_bam2 += '.'+barcodes[i]
                        unmapped_bam2 += '.unmapped.bam'
                        if os.path.isfile(unmapped_bam):
                            os.system('mv ' + unmapped_bam + ' ' + unmapped_bam2)
        
                        # Call run_alignment
                        output_file = '{}/logs/run_alignment_{}_{}_{}_{}.log'.format(output_folder, library, lanes[i], slice, barcodes[i])
                        submission_script = '{}/run_alignment.sh'.format(scripts_folder)
                        call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=60G', '-notify', '-l', 'h_rt=24:0:0', '-j', 'y', '-P', 'macosko_lab', '-l', 'os=RedHat7', submission_script, manifest_file, library, lanes[i], slice, barcodes[i], scripts_folder, output_folder, '{}/{}_{}'.format(library_folder, experiment_date[j], library)]
                        call_to_taskrunner(output_folder, call_args)
                
                if run_barcodematching[j]:
                    puckcaller_path1 = puckcaller_path[j]
                    file1 = '{}/AnalysisOutputs-selected.mat'.format(puckcaller_path1)
                    file2 = '{}/BeadBarcodes.txt'.format(puckcaller_path1)
                    file3 = '{}/BeadLocations.txt'.format(puckcaller_path1)
                    if puckcaller_path1[-1] != '/':
                        puckcaller_path1 += '/'
                    if (not os.path.isfile(file2)) or (not os.path.isfile(file3)):
                        print('{} and/or {} are not found!'.format(file2, file3))
                        if os.path.isfile(file1):
                            output_file = '{}/logs/ExtractBeadBarcode_{}.log'.format(output_folder, library)
                            submission_script = '{}/puckcaller/run_ExtractBeadBarcode.sh'.format(scripts_folder)
                            call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=10g', '-notify', '-l', 'h_rt=5:0:0', '-j', 'y', '-P', 'macosko_lab', '-l', 'os=RedHat7', submission_script, '/broad/software/nonfree/Linux/redhat_7_x86_64/pkgs/matlab_2019a', puckcaller_path1, scripts_folder, output_folder]
                            call_to_taskrunner(output_folder, call_args)
                        else:
                            print('{} is not found!'.format(file1))
                
                if is_NovaSeq or is_NovaSeq_S4:
                    time.sleep(1800)
                else:
                    time.sleep(600)
                
                # Call run_analysis
                output_file = '{}/logs/run_analysis_{}.log'.format(output_folder, library)
                submission_script = '{}/run_analysis.sh'.format(scripts_folder)
                call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=30g', '-notify', '-l', 'h_rt=100:0:0', '-j', 'y', '-P', 'macosko_lab', '-l', 'os=RedHat7', submission_script, manifest_file, library, scripts_folder, output_folder, '{}/{}_{}'.format(library_folder, experiment_date[j], library)]
                call_to_taskrunner(output_folder, call_args)

        call(['mv', folder_running, folder_finished])
    except Exception as exp:
        print("EXCEPTION:!")
        print(exp)
        traceback.print_tb(exp.__traceback__, file=sys.stdout)
        if os.path.isdir(folder_running):
            call(['mv', folder_running, folder_failed])
        elif os.path.isdir(folder_waiting):
            call(['mv', folder_waiting, folder_failed])
        else:
            call(['mkdir', '-p', folder_failed])
            
        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = "The Slide-seq workflow failed at the step of merging barcode matrics. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)

        sys.exit()
        
    
if __name__ == "__main__":
    main()


