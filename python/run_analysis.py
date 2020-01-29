#!/usr/bin/python

# This script is to combine bam files from slice alignments

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
    locus_function_list = 'exonic+intronic'
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
                locus_function_list = row[row0.index('locus_function_list')]
                run_barcodematching = str2bool(row[row0.index('run_barcodematching')])
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
    
    folder_waiting = '{}/status/waiting.analysis_{}'.format(output_folder, library)
    folder_running = '{}/status/running.analysis_{}'.format(output_folder, library)
    folder_finished = '{}/status/finished.analysis_{}'.format(output_folder, library)
    folder_failed = '{}/status/failed.analysis_{}'.format(output_folder, library)
    
    analysis_folder = '{}/{}_{}'.format(library_folder, experiment_date, library)

    call(['mkdir', folder_waiting])

    # Wait for all of run_alignment finish
    failed_list = []
    while 1:
        f = True
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for slice in slice_id[lanes[i]]:
                fol1 = '{}/status/finished.alignment_{}_{}_{}'.format(output_folder, library, lanes[i], slice)
                fol2 = '{}/status/failed.alignment_{}_{}_{}'.format(output_folder, library, lanes[i], slice)
                if (not os.path.isdir(fol1)) and (not os.path.isdir(fol2)):
                    f = False
                prefix_libraries = '{}/{}.{}.{}.{}'.format(analysis_folder, flowcell_barcode, lanes[i], slice, library)
                if (barcodes[i]):
                    prefix_libraries += '.'+barcodes[i]
                star_bamfile = prefix_libraries + '.star_gene_exon_tagged2.bam'
                if (os.path.isdir(fol1) or os.path.isdir(fol2)) and (not os.path.isfile(star_bamfile)):
                    if star_bamfile not in failed_list:
                        failed_list.append(star_bamfile)
                        if os.path.isdir(fol1):
                            call(['rm', '-r', fol1])
                        if os.path.isdir(fol2):
                            call(['rm', '-r', fol2])
                        if os.path.isfile(prefix_libraries+'.star.Log.final.out'):
                            call(['rm', prefix_libraries+'.star.Log.final.out'])
                        if os.path.isfile(prefix_libraries+'.star.Log.out'):
                            call(['rm', prefix_libraries+'.star.Log.out'])
                        if os.path.isfile(prefix_libraries+'.star.Log.progress.out'):
                            call(['rm', prefix_libraries+'.star.Log.progress.out'])
                        if os.path.isfile(prefix_libraries+'.star.SJ.out.tab'):
                            call(['rm', prefix_libraries+'.star.SJ.out.tab'])
                        if os.path.isdir(prefix_libraries+'.star._STARtmp'):
                            call(['rm', '-r', prefix_libraries+'.star._STARtmp'])
                        output_file = '{}/logs/run_alignment_{}_{}_{}.log'.format(output_folder, library, lanes[i], slice)
                        submission_script = '{}/run.sh'.format(scripts_folder)
                        call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=35g', '-notify', '-l', 'h_rt=10:0:0', '-j', 'y', submission_script, 'run_alignment', manifest_file, library, lanes[i], slice, scripts_folder]
                        call(call_args)
                        f = False
                    else:
                        write_log(log_file, flowcell_barcode, 'MergeSamFiles error: '+star_bamfile+' does not exist!')
                        raise Exception(star_bamfile + ' does not exist!')
        if f:
            break
        time.sleep(60)
    
    if os.path.isdir(folder_waiting):
        call(['mv', folder_waiting, folder_running])
    else:
        call(['mkdir', folder_running])
    
    try:
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
    
        # Merge bam files
        combined_bamfile = '{}/{}.bam'.format(analysis_folder, library)
        commandStr = 'java -Djava.io.tmpdir={} -Dsamjdk.buffer_size=131072 -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8192m '.format(tmpdir)
        commandStr += '-jar {}/picard.jar MergeSamFiles TMP_DIR={} CREATE_INDEX=true CREATE_MD5_FILE=false VALIDATION_STRINGENCY=SILENT '.format(picard_folder, tmpdir)
        commandStr += 'OUTPUT={} SORT_ORDER=coordinate ASSUME_SORTED=true'.format(combined_bamfile)
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for slice in slice_id[lanes[i]]:
                star_bamfile = '{}/{}.{}.{}.{}'.format(analysis_folder, flowcell_barcode, lanes[i], slice, library)
                if (barcodes[i]):
                    star_bamfile += '.'+barcodes[i]
                star_bamfile += '.star_gene_exon_tagged2.bam'
                if not os.path.isfile(star_bamfile):
                    write_log(log_file, flowcell_barcode, 'MergeSamFiles error: '+star_bamfile+' does not exist!')
                    raise Exception(star_bamfile + ' does not exist!')
                commandStr += ' INPUT='+star_bamfile
        write_log(log_file, flowcell_barcode, "MergeSamFiles for "+library+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "MergeSamFiles for "+library+" is done. ")
        
        # Validate bam file
        commandStr = 'java -Djava.io.tmpdir={} -Dsamjdk.buffer_size=131072 -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx16384m '.format(tmpdir)
        commandStr += '-jar {}/picard.jar ValidateSamFile TMP_DIR={} VALIDATION_STRINGENCY=SILENT '.format(picard_folder, tmpdir)
        commandStr += 'INPUT={} MODE=SUMMARY'.format(combined_bamfile)
        if (not is_NovaSeq) and (not is_NovaSeq_S4):
            commandStr += ' IGNORE=MISSING_PLATFORM_VALUE IGNORE=INVALID_VERSION_NUMBER'
        write_log(log_file, flowcell_barcode, "ValidateSamFile for "+library+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "ValidateSamFile for "+library+" is done. ")
        
        # Call generate_plots
        output_file = '{}/logs/generate_plots_{}.log'.format(output_folder, library)
        submission_script = '{}/run.sh'.format(scripts_folder)
        call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=30g', '-notify', '-l', 'h_rt=10:0:0', '-j', 'y', submission_script, 'generate_plots', manifest_file, library, scripts_folder]
        call(call_args)
        
        lists = locus_function_list.split(',')
        referencePure = reference[reference.rfind('/') + 1:]
        referencePure = referencePure[:referencePure.rfind('.')]   
        for l in lists:
            call(['mkdir', '{}/{}.{}'.format(analysis_folder, referencePure, l)])
            call(['mkdir', '{}/{}.{}/alignment'.format(analysis_folder, referencePure, l)])
            
            if run_barcodematching:
                barcode_matching_folder = '{}/{}.{}/barcode_matching/'.format(analysis_folder, referencePure, l)
                call(['mkdir', barcode_matching_folder])
                for i in range(len(lanes)):
                    if libraries[i] != library:
                        continue
                    for slice in slice_id[lanes[i]]:
                        toCopyFile = '{}/{}.{}.{}.{}'.format(analysis_folder, flowcell_barcode, lanes[i], slice, library)
                        if (barcodes[i]):
                            toCopyFile += '.'+barcodes[i]
                        toCopyFile += '.star_gene_exon_tagged2.bam'
                        if os.path.isfile(toCopyFile):
                            call(['cp', toCopyFile, barcode_matching_folder])
            
            # Call run_analysis_spec
            output_file = '{}/logs/run_analysis_spec_{}_{}.log'.format(output_folder, library, l)
            submission_script = '{}/run.sh'.format(scripts_folder)
            call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=30g', '-notify', '-l', 'h_rt=12:0:0', '-j', 'y', submission_script, 'run_analysis_spec', manifest_file, library, scripts_folder, l]
            call(call_args)
        
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for slice in slice_id[lanes[i]]:
                toDeleteFile = '{}/{}.{}.{}.{}'.format(analysis_folder, flowcell_barcode, lanes[i], slice, library)
                if (barcodes[i]):
                    toDeleteFile += '.'+barcodes[i]
                toDeleteFile += '.star_gene_exon_tagged2.bam'
                if os.path.isfile(toDeleteFile):
                    call(['rm', toDeleteFile])
        
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
            content = "The Slide-seq workflow for "+library+" failed at the step of running analysis. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)
        
        sys.exit()
    

if __name__ == "__main__":
    main()


