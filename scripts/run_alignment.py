#!/usr/bin/python

# This script is to run the alignment steps

from __future__ import print_function

import sys
import os
import getopt
import csv

import argparse
import glob
import re
from subprocess import call
from datetime import datetime


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
    if len(sys.argv) != 5:
        print("Please provide four arguments: manifest file, library ID, lane ID, and slice ID!")
        sys.exit()
    
    manifest_file = sys.argv[1]
    library = sys.argv[2]
    lane = sys.argv[3]
    slice = sys.argv[4]

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
    
    # Read info from metadata file
    barcode = ''
    bead_structure = ''
    reference = ''
    locus_function_list = 'exonic+intronic'
    sequence = 'AAGCAGTGGTATCAACGCAGAGTGAATGGG'
    base_quality = '10'
    experiment_date = ''
    with open('{}/parsed_metadata.txt'.format(output_folder), 'r') as fin:
        reader = csv.reader(fin, delimiter='\t')
        rows = list(reader)
        row0 = rows[0]
        for i in range(1, len(rows)):
            row = rows[i]
            if row[row0.index('library')] == library:
                barcode = row[row0.index('sample_barcode')]
                bead_structure = row[row0.index('bead_structure')]
                reference = row[row0.index('reference')]
                locus_function_list = row[row0.index('locus_function_list')]
                sequence = row[row0.index('start_sequence')]
                base_quality = row[row0.index('base_quality')]
                experiment_date = row[row0.index('experiment_date')]
                break
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
        
    prefix_libraries = '{}/{}_{}/{}.{}.{}.{}'.format(library_folder, experiment_date, library, flowcell_barcode, lane, slice, library)
    if (barcode):
        prefix_libraries += '.'+barcode
    
    unmapped_bam = prefix_libraries + '.unmapped.bam'
    if not os.path.isfile(unmapped_bam):
        unmapped_bam1 = '{}/{}/{}/{}/'.format(output_folder, lane, slice, library)
        if (barcode):
            unmapped_bam1 += '{}/{}.{}.{}.{}.{}.unmapped.bam'.format(barcode, flowcell_barcode, lane, slice, library, barcode)
        else:
            unmapped_bam1 += '{}.{}.{}.{}.unmapped.bam'.format(flowcell_barcode, lane, slice, library)
        if os.path.isfile(unmapped_bam1):
            os.system('mv ' + unmapped_bam1 + ' ' + unmapped_bam)
    
    bs_range1 = get_bead_structure_range(bead_structure, 'C')
    bs_range2 = get_bead_structure_range(bead_structure, 'M')
    
    folder_running = '{}/status/running.alignment_{}_{}_{}'.format(output_folder, library, lane, slice)
    folder_finished = '{}/status/finished.alignment_{}_{}_{}'.format(output_folder, library, lane, slice)
    folder_failed = '{}/status/failed.alignment_{}_{}_{}'.format(output_folder, library, lane, slice)

    try:
        call(['mkdir', folder_running])
        
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)
        
        # Tag bam with read sequence extended cellular
        commandStr = '{}/TagBamWithReadSequenceExtended O={}.unaligned_tagged_Cellular.bam COMPRESSION_LEVEL=0 TMP_DIR={}'.format(dropseq_folder, prefix_libraries, tmpdir)
        commandStr += ' SUMMARY={}.unaligned_tagged_Cellular.bam_summary.txt BASE_RANGE={} BASE_QUALITY={}'.format(prefix_libraries, bs_range1, base_quality)
        commandStr += ' BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 I={} VALIDATION_STRINGENCY=SILENT'.format(unmapped_bam)
        write_log(log_file, flowcell_barcode, "TagBamWithReadSequenceExtended Cellular for "+library+" in Lane "+lane+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "TagBamWithReadSequenceExtended Cellular for "+library+" in Lane "+lane+" is done. ")
        
        unaligned_cellular_file = '{}.unaligned_tagged_Cellular.bam'.format(prefix_libraries)
        if not os.path.isfile(unaligned_cellular_file):
            write_log(log_file, flowcell_barcode, 'TagBamWithReadSequenceExtended error: '+unaligned_cellular_file+' does not exist!')
            raise Exception('TagBamWithReadSequenceExtended error: '+unaligned_cellular_file+' does not exist!')
        
        # Tag bam with read sequence extended molecular
        commandStr = '{}/TagBamWithReadSequenceExtended O={}.unaligned_tagged_Molecular.bam COMPRESSION_LEVEL=0 TMP_DIR={}'.format(dropseq_folder, prefix_libraries, tmpdir)
        commandStr += ' SUMMARY={}.unaligned_tagged_Molecular.bam_summary.txt BASE_RANGE={} BASE_QUALITY={}'.format(prefix_libraries, bs_range2, base_quality)
        commandStr += ' BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 I={}.unaligned_tagged_Cellular.bam VALIDATION_STRINGENCY=SILENT'.format(prefix_libraries)
        write_log(log_file, flowcell_barcode, "TagBamWithReadSequenceExtended Molecular for "+library+" in Lane "+lane+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "TagBamWithReadSequenceExtended Molecular for "+library+" in Lane "+lane+" is done. ")
        
        unaligned_molecular_file = '{}.unaligned_tagged_Molecular.bam'.format(prefix_libraries)
        if not os.path.isfile(unaligned_molecular_file):
            write_log(log_file, flowcell_barcode, 'TagBamWithReadSequenceExtended error: '+unaligned_molecular_file+' does not exist!')
            raise Exception('TagBamWithReadSequenceExtended error: '+unaligned_molecular_file+' does not exist!')
        
        if os.path.isfile(unaligned_cellular_file):
            call(['rm', unaligned_cellular_file])
        
        # Filter low-quality reads
        commandStr = '{}/FilterBam TAG_REJECT=XQ I={}.unaligned_tagged_Molecular.bam '.format(dropseq_folder, prefix_libraries)
        commandStr += 'O={}.unaligned.filtered.bam COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=SILENT TMP_DIR={} OPTIONS_FILE={}'.format(prefix_libraries, tmpdir, option_file)
        write_log(log_file, flowcell_barcode, "FilterBam for "+library+" in Lane "+lane+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "FilterBam for "+library+" in Lane "+lane+" is done. ")
        
        unaligned_filtered_file = '{}.unaligned.filtered.bam'.format(prefix_libraries)
        if not os.path.isfile(unaligned_filtered_file):
            write_log(log_file, flowcell_barcode, 'FilterBam error: '+unaligned_filtered_file+' does not exist!')
            raise Exception('FilterBam error: '+unaligned_filtered_file+' does not exist!')
        
        if os.path.isfile(unaligned_molecular_file):
            call(['rm', unaligned_molecular_file])

        # Trim reads with starting sequence
        commandStr = '{}/TrimStartingSequence INPUT={}.unaligned.filtered.bam OUTPUT={}.unaligned_trimstartingsequence.filtered.bam COMPRESSION_LEVEL=0 TMP_DIR={} '.format(dropseq_folder, prefix_libraries, prefix_libraries, tmpdir)
        commandStr += 'OUTPUT_SUMMARY={}.adapter_trimming_report.txt SEQUENCE={} MISMATCHES=0 NUM_BASES=5 VALIDATION_STRINGENCY=SILENT'.format(prefix_libraries, sequence)
        write_log(log_file, flowcell_barcode, "TrimStartingSequence for "+library+" in Lane "+lane+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "TrimStartingSequence for "+library+" in Lane "+lane+" is done. ")
        
        adapter_trim_file = '{}.unaligned_trimstartingsequence.filtered.bam'.format(prefix_libraries)
        if not os.path.isfile(adapter_trim_file):
            write_log(log_file, flowcell_barcode, 'TrimStartingSequence error: '+adapter_trim_file+' does not exist!')
            raise Exception('TrimStartingSequence error: '+adapter_trim_file+' does not exist!')

        if os.path.isfile(unaligned_filtered_file):
            call(['rm', unaligned_filtered_file])

        # Adapter-aware poly A trimming
        commandStr = '{}/PolyATrimmer I={}.unaligned_trimstartingsequence.filtered.bam O={}.unaligned_mc_tagged_polyA_filtered.bam TMP_DIR={} '.format(dropseq_folder, prefix_libraries, prefix_libraries, tmpdir)
        commandStr += 'OUTPUT_SUMMARY={}.polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6 VALIDATION_STRINGENCY=SILENT USE_NEW_TRIMMER=true'.format(prefix_libraries)
        write_log(log_file, flowcell_barcode, "PolyATrimmer for "+library+" in Lane "+lane+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "PolyATrimmer for "+library+" in Lane "+lane+" is done. ")
        
        polyA_trim_file = '{}.unaligned_mc_tagged_polyA_filtered.bam'.format(prefix_libraries)
        if not os.path.isfile(polyA_trim_file):
            write_log(log_file, flowcell_barcode, 'PolyATrimmer error: '+polyA_trim_file+' does not exist!')
            raise Exception('PolyATrimmer error: '+polyA_trim_file+' does not exist!')
        
        if os.path.isfile(adapter_trim_file):
            call(['rm', adapter_trim_file])
        
        # Convert bam to fastq
        commandStr = 'java -Djava.io.tmpdir={} -Xmx500m -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 '.format(tmpdir)
        commandStr += '-jar {}/picard.jar SamToFastq I={}.unaligned_mc_tagged_polyA_filtered.bam F={}.fastq VALIDATION_STRINGENCY=SILENT'.format(picard_folder, prefix_libraries, prefix_libraries)
        write_log(log_file, flowcell_barcode, "SamToFastq for "+library+" in Lane "+lane+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "SamToFastq for "+library+" in Lane "+lane+" is done. ")
        
        fastq_file = '{}.fastq'.format(prefix_libraries)
        if not os.path.isfile(fastq_file):
            write_log(log_file, flowcell_barcode, 'SamToFastq error: '+fastq_file+' does not exist!')
            raise Exception('SamToFastq error: '+fastq_file+' does not exist!')
        
        # Map reads to genome sequence using STAR
        commandStr = '{}/STAR --genomeDir {} --readFilesIn {}.fastq '.format(STAR_folder, genome_dir, prefix_libraries)
        commandStr += '--outFileNamePrefix {}.star. --outStd Log --outSAMtype BAM Unsorted --outBAMcompression 0'.format(prefix_libraries)
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += ' --limitOutSJcollapsed 5000000'
        write_log(log_file, flowcell_barcode, "Mapping using STAR for "+library+" in Lane "+lane+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "Mapping using STAR for "+library+" in Lane "+lane+" is done. ")
        
        star_file = '{}.star.Aligned.out.bam'.format(prefix_libraries)
        if not os.path.isfile(star_file):
            write_log(log_file, flowcell_barcode, 'STAR error: '+star_file+' does not exist!')
            raise Exception('STAR error: '+star_file+' does not exist!')
        
        if os.path.isfile(fastq_file):
            call(['rm', fastq_file])
        
        # Sort aligned bam
        commandStr = 'java -Djava.io.tmpdir={} -Xmx4000m -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 '.format(tmpdir)
        commandStr += '-jar {}/picard.jar SortSam I={}.star.Aligned.out.bam '.format(picard_folder, prefix_libraries)
        commandStr += 'O={}.aligned.sorted.bam SORT_ORDER=queryname VALIDATION_STRINGENCY=SILENT TMP_DIR={}'.format(prefix_libraries, tmpdir)
        write_log(log_file, flowcell_barcode, "SortSam for "+library+" in Lane "+lane+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "SortSam for "+library+" in Lane "+lane+" is done. ")
        
        sortsam_file = '{}.aligned.sorted.bam'.format(prefix_libraries)
        if not os.path.isfile(sortsam_file):
            write_log(log_file, flowcell_barcode, 'SortSam error: '+sortsam_file+' does not exist!')
            raise Exception('SortSam error: '+sortsam_file+' does not exist!')
        
        if os.path.isfile(star_file):
            call(['rm', star_file])
        
        # Merge unmapped bam and aligned bam
        commandStr = 'java -Djava.io.tmpdir={} -Xmx8192m -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 '.format(tmpdir)
        commandStr += '-jar {}/picard.jar MergeBamAlignment R={} UNMAPPED={}.unaligned_mc_tagged_polyA_filtered.bam '.format(picard_folder, reference, prefix_libraries)
        commandStr += 'ALIGNED={}.aligned.sorted.bam O={}.merged.bam COMPRESSION_LEVEL=0 INCLUDE_SECONDARY_ALIGNMENTS=false CLIP_ADAPTERS=false '.format(prefix_libraries, prefix_libraries)
        commandStr += 'VALIDATION_STRINGENCY=SILENT TMP_DIR={}'.format(tmpdir)
        write_log(log_file, flowcell_barcode, "MergeBamAlignment for "+library+" in Lane "+lane+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "MergeBamAlignment for "+library+" in Lane "+lane+" is done. ")
        
        mergedbam_file = '{}.merged.bam'.format(prefix_libraries)
        if not os.path.isfile(mergedbam_file):
            write_log(log_file, flowcell_barcode, 'MergeBamAlignment error: '+mergedbam_file+' does not exist!')
            raise Exception('MergeBamAlignment error: '+mergedbam_file+' does not exist!')
        
        if os.path.isfile(polyA_trim_file):
            call(['rm', polyA_trim_file])
        if os.path.isfile(sortsam_file):
            call(['rm', sortsam_file])
        
        # Tag read with interval
        commandStr = '{}/TagReadWithInterval I={}.merged.bam O={}.merged.TagReadWithInterval.bam '.format(dropseq_folder, prefix_libraries, prefix_libraries)
        commandStr += 'COMPRESSION_LEVEL=0 TMP_DIR={} INTERVALS={} TAG=XG VALIDATION_STRINGENCY=SILENT'.format(tmpdir, intervals)
        write_log(log_file, flowcell_barcode, "TagReadWithInterval for "+library+" in Lane "+lane+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "TagReadWithInterval for "+library+" in Lane "+lane+" is done. ")
        
        merged_taginterval_file = '{}.merged.TagReadWithInterval.bam'.format(prefix_libraries)
        if not os.path.isfile(merged_taginterval_file):
            write_log(log_file, flowcell_barcode, 'TagReadWithInterval error: '+merged_taginterval_file+' does not exist!')
            raise Exception('TagReadWithInterval error: '+merged_taginterval_file+' does not exist!')
        
        if os.path.isfile(mergedbam_file):
            call(['rm', mergedbam_file])
        
        # Tag read with gene function
        commandStr = '{}/TagReadWithGeneFunction I={}.merged.TagReadWithInterval.bam O={}.star_gene_exon_tagged2.bam '.format(dropseq_folder, prefix_libraries, prefix_libraries)
        commandStr += 'ANNOTATIONS_FILE={} TMP_DIR={} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=false'.format(annotations_file, tmpdir)
        write_log(log_file, flowcell_barcode, "TagReadWithGeneFunction for "+library+" in Lane "+lane+" Command="+commandStr)
        os.system(commandStr)
        write_log(log_file, flowcell_barcode, "TagReadWithGeneFunction for "+library+" in Lane "+lane+" is done. ")
        
        merged_taggenefunc_file = '{}.star_gene_exon_tagged2.bam'.format(prefix_libraries)
        if not os.path.isfile(merged_taggenefunc_file):
            write_log(log_file, flowcell_barcode, 'TagReadWithGeneFunction error: '+merged_taggenefunc_file+' does not exist!')
            raise Exception('TagReadWithGeneFunction error: '+merged_taggenefunc_file+' does not exist!')

        if os.path.isfile(merged_taginterval_file):
            call(['rm', merged_taginterval_file])
        
        ToFolder = '{}/{}/{}/{}/'.format(output_folder, lane, slice, library)
        if (barcode):
            ToFolder += barcode + '/'
        if os.path.isfile(prefix_libraries+'.star.Log.final.out'):
            call(['mv', prefix_libraries+'.star.Log.final.out', ToFolder])
        if os.path.isfile(prefix_libraries+'.star.Log.out'):
            call(['mv', prefix_libraries+'.star.Log.out', ToFolder])
        if os.path.isfile(prefix_libraries+'.star.Log.progress.out'):
            call(['mv', prefix_libraries+'.star.Log.progress.out', ToFolder])
        if os.path.isfile(prefix_libraries+'.star.SJ.out.tab'):
            call(['mv', prefix_libraries+'.star.SJ.out.tab', ToFolder])
        if os.path.isdir(prefix_libraries+'.star._STARtmp'):
            call(['mv', prefix_libraries+'.star._STARtmp', ToFolder])
        
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
            content = "The Slide-seq workflow for "+library+" in lane "+lane+" slice "+slice+" failed at the step of running alignment. Please check the log file for the issues. "
            call_args = ['python', '{}/send_email.py'.format(scripts_folder), email_address, subject, content]
            call(call_args)
        
        sys.exit()


if __name__ == "__main__":
    main()


