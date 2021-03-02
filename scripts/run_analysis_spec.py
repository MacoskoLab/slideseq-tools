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
import random
from random import sample

from new_submit_to_taskrunner import call_to_taskrunner
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
    run_barcodematching = False
    puckcaller_path = ''
    bead_type = '180402'
    sequence = 'AAGCAGTGGTATCAACGCAGAGTGAATGGG'
    base_quality = '10'
    min_transcripts_per_cell = '10'
    email_address = ''
    experiment_date = ''
    gen_downsampling = False
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
                email_address = row[row0.index('email')]
                run_barcodematching = str2bool(row[row0.index('run_barcodematching')])
                puckcaller_path = row[row0.index('puckcaller_path')]
                bead_type = row[row0.index('bead_type')]
                experiment_date = row[row0.index('date')]
                if 'gen_downsampling' in row0:
                    gen_downsampling = str2bool(row[row0.index('gen_downsampling')])
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
    
    folder_running = '{}/status/running.analysis_spec_{}_{}'.format(output_folder, library, locus_function_list)
    folder_finished = '{}/status/finished.analysis_spec_{}_{}'.format(output_folder, library, locus_function_list)
    folder_failed = '{}/status/failed.analysis_spec_{}_{}'.format(output_folder, library, locus_function_list)
    
    analysis_folder = '{}/{}_{}'.format(library_folder, experiment_date, library)
    alignment_folder = '{}/{}/alignment/'.format(analysis_folder, reference2)
    barcode_matching_folder = '{}/{}/barcode_matching/'.format(analysis_folder, reference2)
    combined_bamfile = '{}/{}.bam'.format(analysis_folder, library)
    
    call(['mkdir', '-p', folder_running])

    try:
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
    
        # Select cells by num transcripts
        commandStr = dropseq_folder+'/SelectCellsByNumTranscripts '
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += '-m 24076m I='+combined_bamfile+' MIN_TRANSCRIPTS_PER_CELL='+min_transcripts_per_cell+' READ_MQ='+base_quality
        else:
            commandStr += '-m 7692m I='+combined_bamfile+' MIN_TRANSCRIPTS_PER_CELL='+min_transcripts_per_cell+' READ_MQ='+base_quality
        commandStr += ' OUTPUT='+alignment_folder+library+'.'+min_transcripts_per_cell+'_transcripts_mq_'+base_quality+'_selected_cells.txt.gz '
        commandStr += 'TMP_DIR='+tmpdir+' VALIDATION_STRINGENCY=SILENT'
        if locus_function_list == 'exonic+intronic':
            commandStr += ' LOCUS_FUNCTION_LIST=INTRONIC'
        elif locus_function_list == 'intronic':
            commandStr += ' LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC'
        write_log(log_file, flowcell_barcode, "SelectCellsByNumTranscripts for "+library+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "SelectCellsByNumTranscripts for "+library+" is done. ")

        # Call run_cmatcher
        if run_barcodematching:
            finish_file = '{}/BeadBarcodes_degenerate.finished'.format(analysis_folder)
            while 1:
                if os.path.isfile(finish_file):
                    call(['rm', finish_file])
                    break
                time.sleep(30)
            
            bead_barcode_file = '{}/BeadBarcodes_degenerate.txt'.format(analysis_folder)
            select_cell_gzfile = alignment_folder+library+'.'+min_transcripts_per_cell+'_transcripts_mq_'+base_quality+'_selected_cells.txt.gz'
            select_cell_file = alignment_folder+library+'.'+min_transcripts_per_cell+'_transcripts_mq_'+base_quality+'_selected_cells.txt'
            name = library+'.'+min_transcripts_per_cell+'_transcripts_mq_'+base_quality+'_selected_cells'
            name_shuffled = library+'.'+min_transcripts_per_cell+'_transcripts_mq_'+base_quality+'_selected_cells.shuffled'
            os.system('gunzip -c '+select_cell_gzfile+' > '+select_cell_file)
            
            select_cell_shuffled_file = alignment_folder+library+'.'+min_transcripts_per_cell+'_transcripts_mq_'+base_quality+'_selected_cells.shuffled.txt'
            with open(select_cell_shuffled_file, 'w') as fout:
                with open(select_cell_file, 'r') as fin:
                    for line in fin:
                        line = line.strip(' \t\n')
                        items = list(line)
                        random.shuffle(items)
                        bc = ''.join(items)
                        fout.write(bc + '\n')
                fin.close()
            fout.close()

            l = 0
            with open(select_cell_file, 'r') as fin:
                for line in fin:
                    l += 1
            fin.close()
            k = 10000
            ls = l // k

            for i in range(ls + 1):
                if i * k >= l:
                    break;
                
                # real barcodes
                infile2 = '{}/{}_{}.txt'.format(alignment_folder, name, str(i + 1))
                commandStr = 'awk \'NR >= {} && NR <= {}\' {} > {}'.format(str(i * k + 1), str((i+1) * k), select_cell_file, infile2)
                os.system(commandStr)
                
                file4 = '{}/{}_barcode_matching_distance_{}.txt'.format(barcode_matching_folder, library, str(i + 1))
                file5 = '{}/{}_barcode_matching_{}.txt'.format(barcode_matching_folder, library, str(i + 1))
                output_file = '{}/logs/run_cmatcher_{}_{}_{}.log'.format(output_folder, library, locus_function_list, str(i + 1))
                submission_script = '{}/run_cmatcher.sh'.format(scripts_folder)
                call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=30G', '-notify', '-l', 'h_rt=26:0:0', '-j', 'y', '-P', 'macosko_lab', '-l', 'os=RedHat7', submission_script, scripts_folder, bead_barcode_file, infile2, file4, file5, bead_type, output_folder, barcode_matching_folder]
                call_to_taskrunner(output_folder, call_args)
                write_log(log_file, flowcell_barcode, "Run CMatcher for "+library+" "+reference2+" "+str(i + 1))
                
                # shuffled barcodes
                infile2 = '{}/{}_{}.txt'.format(alignment_folder, name_shuffled, str(i + 1))
                commandStr = 'awk \'NR >= {} && NR <= {}\' {} > {}'.format(str(i * k + 1), str((i+1) * k), select_cell_shuffled_file, infile2)
                os.system(commandStr)
                
                file4 = '{}/{}_barcode_matching_distance_shuffled_{}.txt'.format(barcode_matching_folder, library, str(i + 1))
                file5 = '{}/{}_barcode_matching_shuffled_{}.txt'.format(barcode_matching_folder, library, str(i + 1))
                output_file = '{}/logs/run_cmatcher_{}_{}_shuffled_{}.log'.format(output_folder, library, locus_function_list, str(i + 1))
                submission_script = '{}/run_cmatcher.sh'.format(scripts_folder)
                call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=30G', '-notify', '-l', 'h_rt=26:0:0', '-j', 'y', '-P', 'macosko_lab', '-l', 'os=RedHat7', submission_script, scripts_folder, bead_barcode_file, infile2, file4, file5, bead_type, output_folder, barcode_matching_folder]
                call_to_taskrunner(output_folder, call_args)
                write_log(log_file, flowcell_barcode, "Run CMatcher for "+library+" "+reference2+" "+str(i + 1))
                
            # Call run_cmatcher_combine
            output_file = '{}/logs/run_cmatcher_combine_{}_{}.log'.format(output_folder, library, locus_function_list)
            submission_script = '{}/run_cmatcher_combine.sh'.format(scripts_folder)
            call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=10G', '-notify', '-l', 'h_rt=48:0:0', '-j', 'y', '-P', 'macosko_lab', '-l', 'os=RedHat7', submission_script, manifest_file, library, scripts_folder, locus_function_list, output_folder, '{}/{}'.format(analysis_folder, reference2)]
            call_to_taskrunner(output_folder, call_args)
        
        # Generate digital expression files for all Illumina barcodes
        commandStr = dropseq_folder+'/DigitalExpression '
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += '-m 32268m '
        else:
            commandStr += '-m 7692m '
        commandStr += 'I='+combined_bamfile+' O='+alignment_folder+library+'.AllIllumina.digital_expression.txt.gz '
        commandStr += 'SUMMARY='+alignment_folder+library+'.AllIllumina.digital_expression_summary.txt EDIT_DISTANCE=1 READ_MQ='+base_quality+' MIN_BC_READ_THRESHOLD=0 '
        commandStr += 'CELL_BC_FILE='+alignment_folder+library+'.'+min_transcripts_per_cell+'_transcripts_mq_'+base_quality+'_selected_cells.txt.gz TMP_DIR='+tmpdir+' '
        commandStr += 'OUTPUT_HEADER=false UEI='+library+' VALIDATION_STRINGENCY=SILENT'
        if locus_function_list == 'exonic+intronic':
            commandStr += ' LOCUS_FUNCTION_LIST=INTRONIC'
        elif locus_function_list == 'intronic':
            commandStr += ' LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC'
        write_log(log_file, flowcell_barcode, "DigitalExpression for "+library+" for all Illumina barcodes Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "DigitalExpression for "+library+" for all Illumina barcodes is done. ")
        
        if gen_downsampling:
            # Downsample bam
            downsample_folder = '{}/{}_{}/{}/downsample/'.format(library_folder, experiment_date, library, reference2)
            call(['mkdir', '-p', downsample_folder])
            f1 = '{}/{}.AllIllumina.digital_expression_summary.txt'.format(alignment_folder, library)
            f2 = '{}/{}_1.digital_expression_summary.txt'.format(downsample_folder, library)
            call(['cp', f1, f2])
            ratio = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
            for i in range(0, 9, 1):
                output_file = '{}/logs/gen_downsample_dge_{}_{}_{}.log'.format(output_folder, library, reference2, str(ratio[i]))
                submission_script = '{}/gen_downsample_dge.sh'.format(scripts_folder)
                call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=37G', '-notify', '-l', 'h_rt=14:0:0', '-j', 'y', '-P', 'macosko_lab', '-l', 'os=RedHat7', submission_script, manifest_file, library, scripts_folder, locus_function_list, str(ratio[i]), output_folder, downsample_folder]
                call_to_taskrunner(output_folder, call_args)
        
            # Call generate_plot_downsampling
            output_file = '{}/logs/generate_plot_downsampling_{}_{}.log'.format(output_folder, library, reference2)
            submission_script = '{}/generate_plot_downsampling.sh'.format(scripts_folder)
            call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=10G', '-notify', '-l', 'h_rt=40:0:0', '-j', 'y', '-P', 'macosko_lab', '-l', 'os=RedHat7', submission_script, manifest_file, library, scripts_folder, locus_function_list, output_folder, barcode_matching_folder]
            call_to_taskrunner(output_folder, call_args)
        
        if not run_barcodematching:
            if os.path.isdir(barcode_matching_folder):
                call(['rm', '-r', barcode_matching_folder])
            if len(email_address) > 1:
                subject = "Slide-seq workflow finished for " + flowcell_barcode
                content = "The Slide-seq workflow for "+library+"_"+locus_function_list+" is finished. Please check the output folder for the results. Thank you for using the Slide-seq tools! "
                call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
                call(call_args)

                output_file = '{}/logs/give_group_{}_{}.log'.format(output_folder, library, reference2)
                submission_script = '{}/give_all_group_write.sh'.format(scripts_folder)
                call_args = ['qsub', '-o', output_file, '-l', 'h_vmem=5G', '-notify', '-l', 'h_rt=10:0:0', '-j', 'y', '-P', 'macosko_lab', '-l', 'os=RedHat7', submission_script]
                call_to_taskrunner(output_folder, call_args)
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
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
            content = "The Slide-seq workflow for "+library+" "+locus_function_list+" failed at the step of running specific analysis. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)

        sys.exit()
    

if __name__ == "__main__":
    main()


