#!/usr/bin/python

# This script is to generate readme.txt

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
    
    # Read info from metadata file
    reference = ''
    base_quality = '10'
    min_transcripts_per_cell = '10'
    email_address = ''
    experiment_date = ''
    with open('{}/parsed_metadata.txt'.format(output_folder), 'r') as fin:
        reader = csv.reader(fin, delimiter='\t')
        rows = list(reader)
        row0 = rows[0]
        for i in range(1, len(rows)):
            row = rows[i]
            if row[row0.index('library')] == library:                
                reference = row[row0.index('reference')]
                base_quality = row[row0.index('base_quality')]
                min_transcripts_per_cell = row[row0.index('min_transcripts_per_cell')]
                email_address = row[row0.index('email')]
                experiment_date = row[row0.index('date')]
    fin.close()
    
    reference_folder = reference[:reference.rfind('/')]
    referencePure = reference[reference.rfind('/') + 1:]
    if (referencePure.endswith('.gz')):
        referencePure = referencePure[:referencePure.rfind('.')]
    referencePure = referencePure[:referencePure.rfind('.')]
    reference2 = referencePure + '.' + locus_function_list
    
    analysis_folder = '{}/{}_{}'.format(library_folder, experiment_date, library)
    alignment_folder = '{}/{}/alignment'.format(analysis_folder, reference2)
    barcode_matching_folder = '{}/{}/barcode_matching'.format(analysis_folder, reference2)
    
    readme_file = '{}/readme.txt'.format(alignment_folder)
    with open(readme_file, 'w') as fout:
        fout.write('\n')
        
        f1 = '{}/BeadBarcodes.txt'.format(analysis_folder)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-bead barcodes from PuckCaller folder\n')
            fout.write('\n')
        
        f1 = '{}/BeadLocations.txt'.format(analysis_folder)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-bead locations from PuckCaller folder\n')
            fout.write('\n')

        f1 = '{}/{}_barcode_matching_01.txt'.format(analysis_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-barcode matching results (hamming distance <= 1) within beads for degenerate barcodes\n')
            fout.write('\n')

        f1 = '{}/{}_barcode_matching_2.txt'.format(analysis_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-barcode matching results (hamming distance >= 2) within beads\n')
            fout.write('\n')

        f1 = '{}/BeadBarcodes_degenerate.txt'.format(analysis_folder)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-the list of degenerate bead barcodes (HD <= 1) and other bead barcodes (HD >= 2)\n')
            fout.write('\n')

        f1 = '{}/{}.bam'.format(analysis_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-tagged aligned bam from STAR\n')
            fout.write('\n')

        f1 = '{}/{}_alignment_quality.pdf'.format(analysis_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-mapping quality plots\n')
            fout.write('\n')

        f1 = '{}/{}_mapping_rate.txt'.format(analysis_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-mapping rate summary\n')
            fout.write('\n')

        f1 = '{}/{}_barcode_matching.txt'.format(barcode_matching_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-barcode matching results between raw Illumina barcodes and bead barcodes in BeadBarcodes_degenerate.txt\n')
            fout.write('\n')
        
        f1 = '{}/{}_barcode_matching_shuffled.txt'.format(barcode_matching_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-barcode matching results between shuffled Illumina barcodes and bead barcodes in BeadBarcodes_degenerate.txt\n')
            fout.write('\n')
            
        f1 = '{}/{}_matched_bead_barcodes.txt'.format(barcode_matching_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-the list of matched bead barcodes\n')
            fout.write('\n')
            
        f1 = '{}/{}_matched_bead_locations.txt'.format(barcode_matching_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-hamming distance and XY coordinates of matched bead barcodes\n')
            fout.write('\n')
        
        f1 = '{}/{}_matched.bam'.format(barcode_matching_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-bam tagged (XC) with matched bead barcodes\n')
            fout.write('\n')

        f1 = '{}/BeadLocationsForR.csv'.format(barcode_matching_folder)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-matched bead barcodes and XY coordinates\n')
            fout.write('\n')
        
        f1 = '{}/{}_XYUMIs.txt'.format(barcode_matching_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-XY coordinates and the number of UMIs\n')
            fout.write('\n')

        f1 = '{}/{}.{}_transcripts_mq_{}_selected_cells.txt.gz'.format(alignment_folder, library, min_transcripts_per_cell, base_quality)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-selected top cells based on min_transcripts_per_cell and base_quality\n')
            fout.write('\n')

        f1 = '{}/{}.AllIllumina.digital_expression.txt.gz'.format(alignment_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-digital expression matrix on all of selected Illumina cell barcodes\n')
            fout.write('\n')

        f1 = '{}/{}.digital_expression_raw.txt.gz'.format(alignment_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-digital expression matrix on all of matched Illumina cell barcodes\n')
            fout.write('\n')

        f1 = '{}/{}.digital_expression_shuffled.txt.gz'.format(alignment_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-digital expression matrix on all of matched shuffled Illumina cell barcodes\n')
            fout.write('\n')
            
        f1 = '{}/{}.digital_expression.txt.gz'.format(alignment_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-digital expression matrix on all of matched bead barcodes\n')
            fout.write('\n')

        f1 = '{}/{}_{}.pdf'.format(alignment_folder, library, reference2)
        if os.path.isfile(f1):
            fout.write(f1 + '\n')
            fout.write('-plots of alignments and barcode matching\n')
            fout.write('\n')
    fout.close()


if __name__ == "__main__":
    main()
    
    
