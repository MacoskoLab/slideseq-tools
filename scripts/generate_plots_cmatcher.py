#!/usr/bin/python

# This script is to generate PDFs for matched barcodes from CMatcher

import csv
import gzip
import logging
import math
import os
import re
import shutil
import sys
import time
from subprocess import call

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from plotnine import (
    aes,
    geom_histogram,
    ggplot,
    ggsave,
    ggtitle,
    scale_x_continuous,
    xlab,
    ylab,
)

from slideseq.util import get_tiles, str2bool

log = logging.getLogger(__name__)


# Get read 1 length
def get_read1_len(bs):
    # 12C8M|*T
    # 7C18X7C8M2X|*T
    b = re.split("[CXM]", bs.split("|")[0])
    i = 0
    for it in b:
        if it:
            i += int(it)
    return i


def main():
    if len(sys.argv) != 4:
        print(
            "Please provide three arguments: manifest file, library ID and locus function list!"
        )
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
    with open(manifest_file, "r") as fp:
        for line in fp:
            key, value = line.rstrip().split("=")
            options[key] = value

    flowcell_directory = options["flowcell_directory"]
    output_folder = options["output_folder"]
    flowcell_barcode = options["flowcell_barcode"]

    library_folder = (
        options["library_folder"]
        if "library_folder" in options
        else "{}/libraries".format(output_folder)
    )
    tmpdir = (
        options["temp_folder"]
        if "temp_folder" in options
        else "{}/tmp".format(output_folder)
    )
    dropseq_folder = (
        options["dropseq_folder"]
        if "dropseq_folder" in options
        else "/broad/macosko/bin/dropseq-tools"
    )
    picard_folder = (
        options["picard_folder"]
        if "picard_folder" in options
        else "/broad/macosko/bin/dropseq-tools/3rdParty/picard"
    )

    scripts_folder = (
        options["scripts_folder"]
        if "scripts_folder" in options
        else "/broad/macosko/jilong/slideseq_pipeline/scripts"
    )
    is_NovaSeq = str2bool(options["is_NovaSeq"]) if "is_NovaSeq" in options else False
    is_NovaSeq_S4 = (
        str2bool(options["is_NovaSeq_S4"]) if "is_NovaSeq_S4" in options else False
    )
    num_slice_NovaSeq = (
        int(options["num_slice_NovaSeq"]) if "num_slice_NovaSeq" in options else 10
    )
    num_slice_NovaSeq_S4 = (
        int(options["num_slice_NovaSeq_S4"])
        if "num_slice_NovaSeq_S4" in options
        else 40
    )

    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    reference = ""
    base_quality = "10"
    min_transcripts_per_cell = "10"

    bead_structure = ""
    experiment_date = ""
    gen_read1_plot = False
    thres_fdr = 1
    with open("{}/parsed_metadata.txt".format(output_folder), "r") as fin:
        reader = csv.reader(fin, delimiter="\t")
        rows = list(reader)
        row0 = rows[0]
        for i in range(1, len(rows)):
            row = rows[i]
            lanes.append(row[row0.index("lane")])
            if row[row0.index("lane")] not in lanes_unique:
                lanes_unique.append(row[row0.index("lane")])
            libraries.append(row[row0.index("library")])
            if row[row0.index("library")] not in libraries_unique:
                libraries_unique.append(row[row0.index("library")])
            barcodes.append(row[row0.index("sample_barcode")])
            bead_structures.append(row[row0.index("bead_structure")])
            if row[row0.index("library")] == library:
                reference = row[row0.index("reference")]
                base_quality = row[row0.index("base_quality")]
                min_transcripts_per_cell = row[row0.index("min_transcripts_per_cell")]

                bead_structure = row[row0.index("bead_structure")]
                experiment_date = row[row0.index("date")]
                if "gen_read1_plot" in row0:
                    gen_read1_plot = str2bool(row[row0.index("gen_read1_plot")])
                if "fdr_threshold" in row0:
                    thres_fdr = int(row[row0.index("fdr_threshold")])

    reference_folder = reference[: reference.rfind("/")]
    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]

    annotations_file = f"{reference_folder}/{referencePure}.gtf"
    ref_flat = f"{reference_folder}/{referencePure}.refFlat"
    ribosomal_intervals = f"{reference_folder}/{referencePure}.rRNA.intervals"
    reference2 = referencePure + "." + locus_function_list

    runinfo_file = f"{flowcell_directory}/RunInfo.xml"

    # Get tile information from RunInfo.xml
    slice_id = {}
    slice_first_tile = {}
    slice_tile_limit = {}
    for lane in lanes_unique:
        tile_nums = get_tiles(runinfo_file, lane)
        tile_cou = len(tile_nums)
        if (not is_NovaSeq) and (not is_NovaSeq_S4):
            slice_id[lane] = ["0"]
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

    analysis_folder = "{}/{}_{}/".format(library_folder, experiment_date, library)
    alignment_folder = "{}/{}/alignment/".format(analysis_folder, reference2)
    barcode_matching_folder = "{}/{}/barcode_matching/".format(
        analysis_folder, reference2
    )
    bead_barcode_file = "{}/BeadBarcodes_degenerate.txt".format(analysis_folder)

    log.info("Merge tagged matched bam files...")
    log.info(
        f"{flowcell_barcode} - Merge tagged matched bam files for {library} {reference2}"
    )

    # matched bam
    matched_bam_file = "{}/{}_matched.bam".format(barcode_matching_folder, library)
    commandStr = (
        f"java -Djava.io.tmpdir={tmpdir}"
        " -Dsamjdk.buffer_size=131072 -XX:+UseParallelGC -XX:GCTimeLimit=20"
        " -XX:GCHeapFreeLimit=10 -Xmx8192m "
    )
    commandStr += (
        "-jar "
        + picard_folder
        + "/picard.jar MergeSamFiles TMP_DIR="
        + tmpdir
        + " CREATE_INDEX=true CREATE_MD5_FILE=false VALIDATION_STRINGENCY=SILENT "
    )
    commandStr += (
        "OUTPUT=" + matched_bam_file + " SORT_ORDER=coordinate ASSUME_SORTED=true"
    )
    for i in range(len(lanes)):
        if libraries[i] != library:
            continue
        for lane_slice in slice_id[lanes[i]]:
            filtered_bam = f"{barcode_matching_folder}/{library}_{lanes[i]}_{lane_slice}_{barcodes[i]}_tagged.bam"
            if os.path.isfile(filtered_bam):
                commandStr += " INPUT=" + filtered_bam
            else:
                log.error(f"{filtered_bam} not found!")

    log.info(f"{flowcell_barcode} MergeSamFiles for {library}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - MergeSamFiles for {library} is done.")

    for i in range(len(lanes)):
        if libraries[i] != library:
            continue
        for lane_slice in slice_id[lanes[i]]:
            filtered_bam = f"{barcode_matching_folder}/{library}_{lanes[i]}_{lane_slice}_{barcodes[i]}_tagged.bam"
            if os.path.isfile(filtered_bam):
                call(["rm", filtered_bam])

    # matched_raw_bam_file
    matched_raw_bam_file = f"{barcode_matching_folder}/{library}_matched_raw.bam"
    commandStr = (
        f"java -Djava.io.tmpdir={tmpdir}"
        " -Dsamjdk.buffer_size=131072 -XX:+UseParallelGC -XX:GCTimeLimit=20"
        " -XX:GCHeapFreeLimit=10 -Xmx8192m "
    )
    commandStr += (
        f"-jar {picard_folder}/picard.jar MergeSamFiles TMP_DIR={tmpdir}"
        " CREATE_INDEX=true CREATE_MD5_FILE=false VALIDATION_STRINGENCY=SILENT "
    )
    commandStr += (
        f"OUTPUT={matched_raw_bam_file} SORT_ORDER=coordinate ASSUME_SORTED=true"
    )
    for i in range(len(lanes)):
        if libraries[i] != library:
            continue
        for lane_slice in slice_id[lanes[i]]:
            filtered_bam = "{}/{}_{}_{}_{}_raw.bam".format(
                barcode_matching_folder, library, lanes[i], lane_slice, barcodes[i]
            )
            if os.path.isfile(filtered_bam):
                commandStr += " INPUT=" + filtered_bam
            else:
                log.error(f"{filtered_bam} not found!")
    os.system(commandStr)

    for i in range(len(lanes)):
        if libraries[i] != library:
            continue
        for lane_slice in slice_id[lanes[i]]:
            filtered_bam = "{}/{}_{}_{}_{}_raw.bam".format(
                barcode_matching_folder, library, lanes[i], lane_slice, barcodes[i]
            )
            if os.path.isfile(filtered_bam):
                call(["rm", filtered_bam])

    # matched_shuffled_bam_file
    matched_shuffled_bam_file = "{}/{}_matched_shuffled.bam".format(
        barcode_matching_folder, library
    )
    commandStr = (
        f"java -Djava.io.tmpdir={tmpdir}"
        " -Dsamjdk.buffer_size=131072 -XX:+UseParallelGC -XX:GCTimeLimit=20"
        " -XX:GCHeapFreeLimit=10 -Xmx8192m "
    )
    commandStr += (
        f"-jar {picard_folder}/picard.jar MergeSamFiles TMP_DIR={tmpdir}"
        " CREATE_INDEX=true CREATE_MD5_FILE=false VALIDATION_STRINGENCY=SILENT "
    )
    commandStr += (
        f"OUTPUT={matched_shuffled_bam_file} SORT_ORDER=coordinate ASSUME_SORTED=true"
    )
    for i in range(len(lanes)):
        if libraries[i] != library:
            continue
        for lane_slice in slice_id[lanes[i]]:
            filtered_bam = f"{barcode_matching_folder}/{library}_{lanes[i]}_{lane_slice}_{barcodes[i]}_shuffled.bam"
            if os.path.isfile(filtered_bam):
                commandStr += " INPUT=" + filtered_bam
            else:
                log.error(f"{filtered_bam} not found!")
    os.system(commandStr)

    for i in range(len(lanes)):
        if libraries[i] != library:
            continue
        for lane_slice in slice_id[lanes[i]]:
            filtered_bam = f"{barcode_matching_folder}/{library}_{lanes[i]}_{lane_slice}_{barcodes[i]}_shuffled.bam"
            if os.path.isfile(filtered_bam):
                call(["rm", filtered_bam])

    log.info(
        f"{flowcell_barcode} - Merge tagged matched bam files for "
        + library
        + " "
        + reference2
        + " is done. ",
    )

    # Generate digital expression files on shuffled Illumina barcodes
    log.info("Generate digital expression files on shuffled Illumina barcodes...")
    shuffled_bci_file = "{}/{}_unique_shuffled_illumina_barcodes.txt.gz".format(
        barcode_matching_folder, library
    )
    matched_dge_summary_shuffled = (
        "{}/{}.digital_expression_summary_shuffled.txt".format(
            alignment_folder, library
        )
    )
    commandStr = (
        dropseq_folder
        + "/DigitalExpression -m 7692m I="
        + matched_shuffled_bam_file
        + " O="
        + alignment_folder
        + library
        + ".digital_expression_shuffled.txt.gz "
    )
    commandStr += (
        "SUMMARY="
        + matched_dge_summary_shuffled
        + " EDIT_DISTANCE=1 READ_MQ="
        + base_quality
        + " MIN_BC_READ_THRESHOLD=0 "
    )
    commandStr += "CELL_BC_FILE=" + shuffled_bci_file + " TMP_DIR=" + tmpdir + " "
    commandStr += "OUTPUT_HEADER=false UEI=" + library + " VALIDATION_STRINGENCY=SILENT"
    if locus_function_list == "exonic+intronic":
        commandStr += " LOCUS_FUNCTION_LIST=INTRONIC"
    elif locus_function_list == "intronic":
        commandStr += " LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
    os.system(commandStr)
    if os.path.isfile(matched_shuffled_bam_file):
        call(["rm", matched_shuffled_bam_file])
    if os.path.isfile(shuffled_bci_file):
        call(["rm", shuffled_bci_file])
    matched_shuffled_bai_file = "{}/{}_matched_shuffled.bai".format(
        barcode_matching_folder, library
    )
    if os.path.isfile(matched_shuffled_bai_file):
        call(["rm", matched_shuffled_bai_file])

    # Generate digital expression files on raw Illumina barcodes
    log.info("Generate digital expression files on raw Illumina barcodes...")
    unique_bci_file = "{}/{}_unique_matched_illumina_barcodes.txt.gz".format(
        barcode_matching_folder, library
    )
    matched_dge_summary_raw = "{}/{}.digital_expression_summary_raw.txt".format(
        alignment_folder, library
    )
    commandStr = (
        dropseq_folder
        + "/DigitalExpression -m 7692m I="
        + matched_raw_bam_file
        + " O="
        + alignment_folder
        + library
        + ".digital_expression_raw.txt.gz "
    )
    commandStr += (
        "SUMMARY="
        + matched_dge_summary_raw
        + " EDIT_DISTANCE=1 READ_MQ="
        + base_quality
        + " MIN_BC_READ_THRESHOLD=0 "
    )
    commandStr += "CELL_BC_FILE=" + unique_bci_file + " TMP_DIR=" + tmpdir + " "
    commandStr += "OUTPUT_HEADER=false UEI=" + library + " VALIDATION_STRINGENCY=SILENT"
    if locus_function_list == "exonic+intronic":
        commandStr += " LOCUS_FUNCTION_LIST=INTRONIC"
    elif locus_function_list == "intronic":
        commandStr += " LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
    os.system(commandStr)
    if os.path.isfile(unique_bci_file):
        call(["rm", unique_bci_file])

    # remove matched beads with high FDR
    dict1 = {}
    dict2 = {}
    dge_barcodes_raw = np.loadtxt(
        matched_dge_summary_raw,
        delimiter="\t",
        dtype="str",
        skiprows=7,
        usecols=0,
    )
    dge_transcripts_raw = np.loadtxt(
        matched_dge_summary_raw,
        delimiter="\t",
        dtype="int",
        skiprows=7,
        usecols=2,
    )
    for i in range(len(dge_barcodes_raw)):
        dict1[dge_barcodes_raw[i]] = dge_transcripts_raw[i]
    bm_raw = f"{barcode_matching_folder}/{library}_barcode_matching.txt"
    bm_raw1 = np.loadtxt(bm_raw, delimiter="\t", dtype="str", skiprows=1, usecols=0)
    bm_raw2 = np.loadtxt(bm_raw, delimiter="\t", dtype="str", skiprows=1, usecols=2)
    for i in range(len(bm_raw1)):
        if bm_raw1[i] in dict1:
            if bm_raw2[i] in dict2:
                dict2[bm_raw2[i]] += dict1[bm_raw1[i]]
            else:
                dict2[bm_raw2[i]] = dict1[bm_raw1[i]]

    dict3 = {}
    dict4 = {}
    dge_barcodes_shuffled = np.loadtxt(
        matched_dge_summary_shuffled,
        delimiter="\t",
        dtype="str",
        skiprows=7,
        usecols=0,
    )
    dge_transcripts_shuffled = np.loadtxt(
        matched_dge_summary_shuffled,
        delimiter="\t",
        dtype="int",
        skiprows=7,
        usecols=2,
    )
    for i in range(len(dge_barcodes_shuffled)):
        dict3[dge_barcodes_shuffled[i]] = dge_transcripts_shuffled[i]
    bm_shuffled = f"{barcode_matching_folder}/{library}_barcode_matching_shuffled.txt"
    bm_shuffled1 = np.loadtxt(
        bm_shuffled, delimiter="\t", dtype="str", skiprows=1, usecols=0
    )
    bm_shuffled2 = np.loadtxt(
        bm_shuffled, delimiter="\t", dtype="str", skiprows=1, usecols=2
    )
    for i in range(len(bm_shuffled1)):
        if bm_shuffled1[i] in dict3:
            if bm_shuffled2[i] in dict4:
                dict4[bm_shuffled2[i]] += dict3[bm_shuffled1[i]]
            else:
                dict4[bm_shuffled2[i]] = dict3[bm_shuffled1[i]]

    matched_bead_barcode_file = (
        f"{barcode_matching_folder}/{library}_matched_bead_barcodes.txt"
    )
    matched_barcodes = np.loadtxt(
        matched_bead_barcode_file, delimiter="\t", dtype="str", usecols=0
    )
    non_selected_bc = []
    for i in range(len(matched_barcodes)):
        bc = matched_barcodes[i]
        v_raw = 0
        v_shuffled = 0
        if bc in dict2:
            v_raw = dict2[bc]
        if bc in dict4:
            v_shuffled = dict4[bc]
        ratio = 100 * v_shuffled / (v_shuffled + v_raw)
        if ratio > thres_fdr:
            non_selected_bc.append(bc)

    # modified files based on new matched beads
    matched_bead_location_file = (
        f"{barcode_matching_folder}/{library}_matched_bead_locations.txt"
    )
    matched_bead_barcode_file2 = (
        f"{barcode_matching_folder}/{library}_matched_bead_barcodes_revised.txt"
    )
    matched_bead_location_file2 = (
        f"{barcode_matching_folder}/{library}_matched_bead_locations_revised.txt"
    )
    with open(matched_bead_barcode_file2, "w") as fout1:
        with open(matched_bead_location_file2, "w") as fout2:
            with open(matched_bead_location_file, "r") as fin:
                j = -1
                for line in fin:
                    j += 1
                    if not matched_barcodes[j] in non_selected_bc:
                        fout1.write(matched_barcodes[j] + "\n")
                        fout2.write(line)

    if os.path.isfile(matched_bead_barcode_file):
        call(["rm", matched_bead_barcode_file])
    if os.path.isfile(matched_bead_location_file):
        call(["rm", matched_bead_location_file])

    os.system(f"mv {matched_bead_barcode_file2} {matched_bead_barcode_file}")
    os.system(f"mv {matched_bead_location_file2} {matched_bead_location_file}")

    matched_bead_barcode_gzfile = (
        f"{barcode_matching_folder}/{library}_matched_bead_barcodes.txt.gz"
    )
    os.system(
        "gzip -c " + matched_bead_barcode_file + " > " + matched_bead_barcode_gzfile
    )

    # Generate digital expression files on matched bead barcodes
    log.info("Generate digital expression files on matched bead barcodes...")
    commandStr = (
        dropseq_folder
        + "/DigitalExpression -m 7692m I="
        + matched_bam_file
        + " O="
        + alignment_folder
        + library
        + ".digital_expression.txt.gz "
    )
    commandStr += (
        "SUMMARY="
        + alignment_folder
        + library
        + ".digital_expression_summary.txt EDIT_DISTANCE=1 READ_MQ="
        + base_quality
        + " MIN_BC_READ_THRESHOLD=0 "
    )
    commandStr += (
        "CELL_BC_FILE=" + matched_bead_barcode_gzfile + " TMP_DIR=" + tmpdir + " "
    )
    commandStr += "OUTPUT_HEADER=false UEI=" + library + " VALIDATION_STRINGENCY=SILENT"
    if locus_function_list == "exonic+intronic":
        commandStr += " LOCUS_FUNCTION_LIST=INTRONIC"
    elif locus_function_list == "intronic":
        commandStr += " LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
    log.info(
        f"{flowcell_barcode} - DigitalExpression for {library} on matched barcodes"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - DigitalExpression for {library} on matched barcodes is done."
    )

    # Call gen_sparse_matrix
    file_name = f"{library}.digital_expression"
    output_file = (
        f"{output_folder}/logs/gen_sparse_matrix_{library}_{locus_function_list}.log"
    )
    submission_script = f"{scripts_folder}/gen_sparse_matrix.sh"
    call_args = [
        "qsub",
        "-o",
        output_file,
        submission_script,
        manifest_file,
        library,
        locus_function_list,
        alignment_folder,
        file_name,
        scripts_folder,
        output_folder,
        f"{library_folder}/{experiment_date}_{library}/{reference2}",
    ]
    call(call_args)

    # # Call write_bijective_mapping.sh...but what we just didn't
    # output_file = f"{output_folder}/logs/write_bijective_mapping_{library}_{locus_function_list}.log"
    # submission_script = f"{scripts_folder}/write_bijective_mapping.sh"
    # call_args = [
    #     "qsub",
    #     "-o",
    #     output_file,
    #     submission_script,
    #     manifest_file,
    #     library,
    #     scripts_folder,
    #     locus_function_list,
    #     output_folder,
    #     f"{analysis_folder}/{reference2}",
    # ]
    # call(call_args)

    # Generate libraryID_XYUMIs.txt
    matched_barcodes = np.loadtxt(
        matched_bead_barcode_file, delimiter="\t", dtype="str", usecols=0
    )
    distances = np.loadtxt(
        matched_bead_location_file, delimiter="\t", dtype="int", usecols=0
    )
    coordinatesx = np.loadtxt(
        matched_bead_location_file, delimiter="\t", dtype="float", usecols=1
    )
    coordinatesy = np.loadtxt(
        matched_bead_location_file, delimiter="\t", dtype="float", usecols=2
    )
    matched_dge_summary = f"{alignment_folder}/{library}.digital_expression_summary.txt"
    matched_dge_summary_barcodes = np.loadtxt(
        matched_dge_summary, delimiter="\t", dtype="str", skiprows=7, usecols=0
    )
    matched_dge_summary_transcripts = np.loadtxt(
        matched_dge_summary, delimiter="\t", dtype="int", skiprows=7, usecols=2
    )
    XYUMIs_file = f"{barcode_matching_folder}/{library}_XYUMIs.txt"
    dictX = {}
    dictY = {}
    for i in range(len(matched_barcodes)):
        dictX[matched_barcodes[i]] = coordinatesx[i]
        dictY[matched_barcodes[i]] = coordinatesy[i]
    with open(XYUMIs_file, "w") as fout:
        for i in range(len(matched_dge_summary_barcodes)):
            barcode = matched_dge_summary_barcodes[i]
            umi = matched_dge_summary_transcripts[i]
            fout.write(f"{dictX[barcode]}\t{dictY[barcode]}\t{umi}\n")

    # Bam tag histogram
    commandStr = (
        f"{dropseq_folder}/BamTagHistogram -m 7692m"
        f" I={matched_bam_file} OUTPUT={alignment_folder}{library}.numReads_perCell_XC_mq_{base_quality}.txt.gz"
        f" TAG=XC FILTER_PCR_DUPLICATES=false TMP_DIR={tmpdir} READ_MQ={base_quality} VALIDATION_STRINGENCY=SILENT"
    )
    log.info(f"{flowcell_barcode} - BamTagHistogram for {library}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - BamTagHistogram for {library} on matched barcodes is done."
    )

    # Bam tag histogram on UMI
    commandStr = (
        f"{dropseq_folder}/BamTagHistogram -m 7692m"
        f" I={matched_bam_file} OUTPUT={alignment_folder}{library}.numReads_perUMI_XM_mq_{base_quality}.txt.gz"
        f" TAG=XM FILTER_PCR_DUPLICATES=false TMP_DIR={tmpdir} READ_MQ={base_quality} VALIDATION_STRINGENCY=SILENT"
    )
    log.info(
        f"{flowcell_barcode} - BamTagHistogram on UMI for {library} on matched barcodes"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - BamTagHistogram on UMI for {library} on matched barcodes is done."
    )

    # Collect RnaSeq metrics
    commandStr = (
        f"java -Djava.io.tmpdir={tmpdir} -XX:+UseParallelGC"
        f" -XX:GCTimeLimit=20 -XX:GCHeapFreeLimit=10 -Xmx8192m"
        f" -jar {picard_folder}/picard.jar CollectRnaSeqMetrics TMP_DIR={tmpdir}"
        f" VALIDATION_STRINGENCY=SILENT I={matched_bam_file} REF_FLAT={ref_flat} STRAND_SPECIFICITY=NONE"
        f" OUTPUT={alignment_folder}{library}.fracIntronicExonic.txt RIBOSOMAL_INTERVALS={ribosomal_intervals}"
    )
    log.info(
        f"{flowcell_barcode} - CollectRnaSeqMetrics for {library} on matched barcodes"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - CollectRnaSeqMetrics for {library} on matched barcodes is done."
    )

    # Base distribution at read position for cellular
    commandStr = (
        f"{dropseq_folder}/BaseDistributionAtReadPosition -m 7692m I={matched_raw_bam_file}"
        f" OUTPUT={alignment_folder}{library}.barcode_distribution_XC.txt TMP_DIR={tmpdir}"
        f" TAG=XC VALIDATION_STRINGENCY=SILENT"
    )
    log.info(
        f"{flowcell_barcode} - BaseDistributionAtReadPosition Cellular for {library}"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - BaseDistributionAtReadPosition Cellular for {library} on matched barcodes is done."
    )

    # Base distribution at read position for molecular
    commandStr = (
        f"{dropseq_folder}/BaseDistributionAtReadPosition -m 7692m I={matched_raw_bam_file}"
        f" OUTPUT={alignment_folder}{library}.barcode_distribution_XM.txt TMP_DIR={tmpdir}"
        f" TAG=XM VALIDATION_STRINGENCY=SILENT"
    )
    log.info(
        f"{flowcell_barcode} - BaseDistributionAtReadPosition Molecular for {library} on matched barcodes"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - BaseDistributionAtReadPosition Molecular for {library} on matched barcodes is done."
    )

    if os.path.isfile(matched_raw_bam_file):
        call(["rm", matched_raw_bam_file])
    matched_raw_bai_file = f"{barcode_matching_folder}/{library}_matched_raw.bai"
    if os.path.isfile(matched_raw_bai_file):
        call(["rm", matched_raw_bai_file])

    # Gather read quality metrics
    commandStr = (
        f"{dropseq_folder}/GatherReadQualityMetrics -m 7692m I={matched_bam_file}"
        f" TMP_DIR={tmpdir} OUTPUT={alignment_folder}{library}.ReadQualityMetrics.txt VALIDATION_STRINGENCY=SILENT"
    )
    log.info(
        f"{flowcell_barcode} - GatherReadQualityMetrics for {library} on matched barcodes"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - GatherReadQualityMetrics for {library} on matched barcodes is done.",
    )

    # Single cell RnaSeq metrics collector
    commandStr = (
        f"{dropseq_folder}/SingleCellRnaSeqMetricsCollector -m 15884m I={matched_bam_file}"
        f" ANNOTATIONS_FILE={annotations_file} OUTPUT={alignment_folder}{library}.fracIntronicExonicPerCell.txt.gz"
        f" RIBOSOMAL_INTERVALS={ribosomal_intervals} CELL_BARCODE_TAG=XC READ_MQ={base_quality}"
        f" TMP_DIR={tmpdir} CELL_BC_FILE={matched_bead_barcode_gzfile} VALIDATION_STRINGENCY=SILENT"
    )
    log.info(f"{flowcell_barcode} - SingleCellRnaSeqMetricsCollector for {library}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - SingleCellRnaSeqMetricsCollector for {library} on matched barcodes is done.",
    )

    if os.path.isfile(
        f"{alignment_folder}/{library}.SelectCellsByNumTranscripts_metrics"
    ):
        call(
            [
                "rm",
                f"{alignment_folder}/{library}.SelectCellsByNumTranscripts_metrics",
            ]
        )
    if os.path.isfile(matched_bead_barcode_gzfile):
        call(["rm", matched_bead_barcode_gzfile])

    log.info("generating plots...")
    log.info(
        f"{flowcell_barcode} - Generate plots for matched barcodes for {library} {reference2}"
    )
    pp = PdfPages(f"{alignment_folder}/{library}_{reference2}.pdf")

    file = f"{analysis_folder}/{library}.ReadQualityMetrics.txt"
    while 1:
        if os.path.isfile(file):
            break
        time.sleep(60)
    time.sleep(10)
    mat = np.loadtxt(
        file,
        delimiter="\t",
        dtype="int",
        skiprows=3,
        max_rows=1,
        usecols=(1, 2, 3, 4),
    )
    file1 = f"{alignment_folder}/{library}.ReadQualityMetrics.txt"
    mat1 = np.loadtxt(
        file1,
        delimiter="\t",
        dtype="int",
        skiprows=3,
        max_rows=1,
        usecols=(1, 2, 3, 4),
    )

    df_z = [mat[0], mat[1], mat[2], mat1[2]]
    if mat[2] >= 1000000 and mat1[2] >= 1000000:
        df_u = [
            int(mat[0] / 1000000),
            int(mat[1] / 1000000),
            int(mat[2] / 1000000),
            int(mat1[2] / 1000000),
        ]
        yl = "# Reads [millions]"
    else:
        df_u = [mat[0], mat[1], mat[2], mat1[2]]
        yl = "# Reads"
    df_y = [
        mat[0] / mat[0] * 100,
        mat[1] / mat[0] * 100,
        mat[2] / mat[0] * 100,
        mat1[2] / mat[0] * 100,
    ]
    df_v = [
        "{:,}".format(mat[0]),
        "{:,}".format(mat[1]),
        "{:,}".format(mat[2]),
        "{:,}".format(mat1[2]),
    ]
    labels = []
    for i in range(4):
        labels.append(f"{df_y[i]:.3g}%")
    df = pd.DataFrame(
        {"x": ["Total", "Mapped", "HQ", "HQ matched"], "z": df_z, "u": df_u}
    )
    fig, ax = plt.subplots(figsize=(8, 8))
    bp = plt.bar(df["x"], df["u"], width=0.7, color="lightskyblue", edgecolor="black")
    for idx, rect in enumerate(bp):
        height = rect.get_height()
        ax.text(
            rect.get_x() + rect.get_width() / 2.0,
            0.5 * height,
            labels[idx],
            ha="center",
            va="bottom",
        )
        ax.text(
            rect.get_x() + rect.get_width() / 2.0,
            height,
            df_v[idx],
            ha="center",
            va="bottom",
        )
    plt.yticks(rotation=90)
    plt.ylabel(yl)
    plt.title("Alignment quality for all reads")
    plt.savefig(pp, format="pdf")

    file = f"{analysis_folder}/{library}.fracIntronicExonic.txt"
    while 1:
        if os.path.isfile(file):
            break
        time.sleep(60)
    time.sleep(10)
    mat = np.loadtxt(
        file,
        delimiter="\t",
        dtype="float",
        skiprows=7,
        max_rows=1,
        usecols=(15, 16, 17, 18, 19),
    )
    mat[1] += mat[2]
    mat[2] = mat[1] + mat[3]
    df_x = ["ribosomal", "exonic", "genic", "intronic", "intergenic"]
    df_y = mat * 100

    labels = []
    for i in range(5):
        labels.append(f"{df_y[i]:.3g}")
    df = pd.DataFrame({"x": df_x, "y": df_y})
    fig, ax = plt.subplots(figsize=(8, 8))
    bp = plt.bar(df["x"], df["y"], width=0.7, color="lightskyblue", edgecolor="black")
    for idx, rect in enumerate(bp):
        height = rect.get_height()
        ax.text(
            rect.get_x() + rect.get_width() / 2.0,
            0.5 * height,
            labels[idx],
            ha="center",
            va="bottom",
        )
    plt.yticks(rotation=90)
    plt.ylabel("Percentage")
    plt.title("All reads")
    plt.savefig(pp, format="pdf")

    f1 = f"{analysis_folder}/{library}.numReads_perCell_XC_mq_{base_quality}.txt"
    f2 = f"{analysis_folder}/{library}.numReads_perCell_XC_mq_{base_quality}.txt.gz"
    while 1:
        if os.path.isfile(f1) or os.path.isfile(f2):
            break
        time.sleep(60)
    time.sleep(10)
    if not os.path.isfile(f1):
        with gzip.open(f2, "rb") as f_in:
            with open(f1, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    mat = np.loadtxt(f1, delimiter="\t", dtype="int", skiprows=1, usecols=0)
    j = int(len(mat) / 10)
    df_x = np.arange(1, j + 1, 1)
    y = np.cumsum(mat)
    df_y = y / max(y)
    df = pd.DataFrame({"x": df_x, "y": df_y[:j]})
    plt.figure(figsize=(8, 8))
    plt.plot(df["x"], df["y"], color="green")
    plt.ylim((0, 1))
    plt.yticks(rotation=90)
    plt.xlabel("cell barcodes sorted by number of reads [descending]")
    plt.ylabel("cumulative fraction of reads")
    plt.title("Cumulative fraction of reads per cell barcode")
    plt.savefig(pp, format="pdf")

    if os.path.isfile(f1):
        call(["rm", f1])

    file1 = f"{analysis_folder}/{library}.ReadQualityMetrics.txt"
    cou1 = np.loadtxt(
        file1, delimiter="\t", dtype="int", skiprows=3, max_rows=1, usecols=1
    )
    j = 42
    f = False
    for i in range(len(lanes)):
        if libraries[i] != library:
            continue
        for lane_slice in slice_id[lanes[i]]:
            file = f"{analysis_folder}/{flowcell_barcode}.{lanes[i]}.{lane_slice}.{library}.{barcodes[i]}.polyA_trimming_report.txt"
            if os.path.isfile(file):
                mat = np.loadtxt(file, delimiter="\t", dtype="int", skiprows=7)
                j = len(mat)
                f = True
                break
        if f:
            break
    df_x = np.arange(0, j, 1)
    df_y = [0] * j
    for i in range(len(lanes)):
        if libraries[i] != library:
            continue
        for lane_slice in slice_id[lanes[i]]:
            file = f"{analysis_folder}/{flowcell_barcode}.{lanes[i]}.{lane_slice}.{library}.{barcodes[i]}.polyA_trimming_report.txt"
            if os.path.isfile(file):
                mat = np.loadtxt(file, delimiter="\t", dtype="int", skiprows=7)
                for j in range(0, len(mat)):
                    df_y[mat[j, 0]] += mat[j, 1]
    max_y = max(df_y[1:])
    min_y = min(df_y[1:])
    cou = sum(df_y[1:])
    val = f"{cou / cou1 * 100:.3g}"
    df = pd.DataFrame({"x": df_x, "y": df_y})
    plt.figure(figsize=(8, 8))
    plt.plot(df["x"], df["y"], color="black")
    plt.xlim((1, j))
    plt.ylim((max(0, min_y - 10000), max_y + 10000))
    plt.yticks(rotation=90)
    plt.xlabel("first base of PolyA tail trimmed")
    plt.ylabel("number of reads")
    plt.title("% Reads trimmed by 3' PolyA trimmer: " + val)
    plt.savefig(pp, format="pdf")

    file = f"{analysis_folder}/{library}.barcode_distribution_XC.txt"
    while 1:
        if os.path.isfile(file):
            break
        time.sleep(60)
    time.sleep(10)
    mat = np.loadtxt(file, delimiter="\t", dtype="int", skiprows=1)
    cou = mat[0, 1] + mat[0, 2] + mat[0, 3] + mat[0, 4]
    df_x = []
    df_x.extend(mat[:, 0])
    df_x.extend(mat[:, 0])
    df_x.extend(mat[:, 0])
    df_x.extend(mat[:, 0])
    df_y = []
    df_y.extend(mat[:, 1] / cou * 100)
    df_y.extend(mat[:, 2] / cou * 100)
    df_y.extend(mat[:, 3] / cou * 100)
    df_y.extend(mat[:, 4] / cou * 100)
    max_y = int(max(df_y))
    j = len(mat)
    fig, ax = plt.subplots(figsize=(8, 8))
    colors = ["red", "blue", "green", "purple"]
    labels = ["A", "C", "G", "T"]
    for i in range(4):
        ax.scatter(
            df_x[i * j : (i + 1) * j],
            df_y[i * j : (i + 1) * j],
            c=colors[i],
            s=20,
            label=labels[i],
        )
    ax.legend(loc="lower right")
    plt.xlim((0, j + 2))
    plt.ylim((0, max_y + 2))
    plt.yticks(rotation=90)
    plt.xlabel("base position")
    plt.ylabel("fraction of reads")
    plt.title("Cell barcodes for all reads")
    plt.savefig(pp, format="pdf")

    file = f"{analysis_folder}/{library}.barcode_distribution_XM.txt"
    while 1:
        if os.path.isfile(file):
            break
        time.sleep(60)
    time.sleep(10)
    mat = np.loadtxt(file, delimiter="\t", dtype="int", skiprows=1)
    cou = mat[0, 1] + mat[0, 2] + mat[0, 3] + mat[0, 4]
    df_x = []
    df_x.extend(mat[:, 0])
    df_x.extend(mat[:, 0])
    df_x.extend(mat[:, 0])
    df_x.extend(mat[:, 0])
    df_y = []
    df_y.extend(mat[:, 1] / cou * 100)
    df_y.extend(mat[:, 2] / cou * 100)
    df_y.extend(mat[:, 3] / cou * 100)
    df_y.extend(mat[:, 4] / cou * 100)
    max_y = int(max(df_y))
    j = len(mat)
    fig, ax = plt.subplots(figsize=(8, 8))
    colors = ["red", "blue", "green", "purple"]
    labels = ["A", "C", "G", "T"]
    for i in range(4):
        ax.scatter(
            df_x[i * j : (i + 1) * j],
            df_y[i * j : (i + 1) * j],
            c=colors[i],
            s=20,
            label=labels[i],
        )
    ax.legend(loc="lower right")
    plt.xlim((0, j + 1))
    plt.ylim((0, max_y + 2))
    plt.yticks(rotation=90)
    plt.xlabel("base position")
    plt.ylabel("fraction of reads")
    plt.title("Molecular barcodes for all reads")
    plt.savefig(pp, format="pdf")

    num_beads = 1
    if os.path.isfile(bead_barcode_file):
        num_beads = len(open(bead_barcode_file).readlines())

    num_cells = 1
    selected_cells = f"{alignment_folder}/{library}.{min_transcripts_per_cell}_transcripts_mq_{base_quality}_selected_cells.txt"
    if os.path.isfile(selected_cells):
        num_cells = len(open(selected_cells).readlines())
        call(["rm", selected_cells])

    num_matched_beads = len(matched_dge_summary_barcodes)

    num_matched_cells = -1
    combined_cmatcher_file = f"{barcode_matching_folder}/{library}_barcode_matching.txt"
    if os.path.isfile(combined_cmatcher_file):
        num_matched_cells = len(open(combined_cmatcher_file).readlines())

    ratio_cell = "{0:.3g}%".format(num_matched_cells / num_cells * 100)
    ratio_bead = "{0:.3g}%".format(num_matched_beads / num_beads * 100)
    xtitle = (
        f"{ratio_cell} Illumina barcodes matching to bead barcodes\n"
        f"{ratio_bead} bead barcodes matching to Illumina barcodes"
    )

    dict5 = {}
    dict6 = {}
    for i in range(len(bm_raw1)):
        dict5[bm_raw1[i]] = bm_raw2[i]
    for i in range(len(bm_shuffled1)):
        dict6[bm_shuffled1[i]] = bm_shuffled2[i]
    sum_raw = 0
    for i in range(len(dge_barcodes_raw)):
        if not dict5[dge_barcodes_raw[i]] in non_selected_bc:
            sum_raw += dge_transcripts_raw[i]
    sum_shuffled = 0
    for i in range(len(dge_barcodes_shuffled)):
        if (not dict6[dge_barcodes_shuffled[i]] in non_selected_bc) and (
            dict6[dge_barcodes_shuffled[i]] in dge_barcodes_raw
        ):
            sum_shuffled += dge_transcripts_shuffled[i]
    ratio_UMIs = f"{100 * sum_shuffled / sum_raw:.3g}%"

    file1 = f"{alignment_folder}/{library}.ReadQualityMetrics.txt"
    mat1 = np.loadtxt(
        file1,
        delimiter="\t",
        dtype="int",
        skiprows=3,
        max_rows=1,
        usecols=(1, 2, 3, 4),
    )
    df_z = [mat1[0], mat1[1], mat1[2], mat1[3]]
    if mat1[3] >= 1000000:
        df_u = [
            mat1[0] / 1000000,
            mat1[1] / 1000000,
            mat1[2] / 1000000,
            mat1[3] / 1000000,
        ]
        yl = "# Reads [millions]"
    else:
        df_u = [mat1[0], mat1[1], mat1[2], mat1[3]]
        yl = "# Reads"
    df_y = [
        mat1[0] / mat1[0] * 100,
        mat1[1] / mat1[0] * 100,
        mat1[2] / mat1[0] * 100,
        mat1[3] / mat1[0] * 100,
    ]
    df_v = [
        "{:,}".format(mat1[0]),
        "{:,}".format(mat1[1]),
        "{:,}".format(mat1[2]),
        "{:,}".format(mat1[3]),
    ]
    labels = []
    for i in range(4):
        labels.append("{0:.3g}%".format(df_y[i]))
    df = pd.DataFrame(
        {"x": ["Total", "Mapped", "HQ", "HQ No Dupes"], "z": df_z, "u": df_u}
    )
    fig, ax = plt.subplots(figsize=(8, 8))
    bp = plt.bar(df["x"], df["u"], width=0.7, color="lightskyblue", edgecolor="black")
    for idx, rect in enumerate(bp):
        height = rect.get_height()
        ax.text(
            rect.get_x() + rect.get_width() / 2.0,
            0.5 * height,
            labels[idx],
            ha="center",
            va="bottom",
        )
        ax.text(
            rect.get_x() + rect.get_width() / 2.0,
            height,
            df_v[idx],
            ha="center",
            va="bottom",
        )
    plt.yticks(rotation=90)
    plt.xlabel(xtitle)
    plt.ylabel(yl)
    plt.title("Alignment quality for matched barcodes (FDR = " + str(ratio_UMIs) + ")")
    plt.savefig(pp, format="pdf")

    file = "{}/{}.fracIntronicExonic.txt".format(alignment_folder, library)
    mat = np.loadtxt(
        file,
        delimiter="\t",
        dtype="float",
        skiprows=7,
        max_rows=1,
        usecols=(15, 16, 17, 18, 19),
    )
    mat[1] += mat[2]
    mat[2] = mat[1] + mat[3]
    df_x = ["ribosomal", "exonic", "genic", "intronic", "intergenic"]
    df_y = mat * 100

    labels = []
    for i in range(5):
        labels.append("{0:.3g}".format(df_y[i]))
    df = pd.DataFrame({"x": df_x, "y": df_y})
    fig, ax = plt.subplots(figsize=(8, 8))
    bp = plt.bar(df["x"], df["y"], width=0.7, color="lightskyblue", edgecolor="black")
    for idx, rect in enumerate(bp):
        height = rect.get_height()
        ax.text(
            rect.get_x() + rect.get_width() / 2.0,
            0.5 * height,
            labels[idx],
            ha="center",
            va="bottom",
        )
    plt.yticks(rotation=90)
    plt.ylabel("Percentage")
    plt.title("All reads that mapped to matched barcodes")
    plt.savefig(pp, format="pdf")

    df_y = [0, 0, 0]
    for i in range(len(distances)):
        df_y[distances[i]] += 1
    df = pd.DataFrame({"x": [0, 1, 2], "y": df_y})
    plt.figure(figsize=(8, 8))
    plt.bar(df["x"], df["y"], width=1, color="lightskyblue", edgecolor="black")
    plt.xlim((-1, 15))
    plt.yticks(rotation=90)
    plt.ylabel("number of matched barcodes")
    plt.title("Histogram of barcode matches by hamming distance")
    plt.savefig(pp, format="pdf")

    matched_dge_summary_reads = np.loadtxt(
        matched_dge_summary, delimiter="\t", dtype="int", skiprows=7, usecols=1
    )
    matched_distances = [0, 0, 0]
    for i in range(len(matched_dge_summary_barcodes)):
        if matched_dge_summary_barcodes[i] in matched_barcodes:
            matched_distances[
                int(
                    distances[
                        np.where(matched_barcodes == matched_dge_summary_barcodes[i])
                    ]
                )
            ] += matched_dge_summary_reads[i]
    for i in range(3):
        matched_distances[i] = matched_distances[i] / 1000000
    df = pd.DataFrame({"x": [0, 1, 2], "y": matched_distances})
    plt.figure(figsize=(8, 8))
    plt.bar(df["x"], df["y"], width=1, color="lightskyblue", edgecolor="black")
    plt.xlim((-1, 15))
    plt.yticks(rotation=90)
    plt.ylabel("# Reads [millions]")
    plt.title("Histogram of reads matches by hamming distance")
    plt.savefig(pp, format="pdf")

    matched_dge_summary_transcripts = np.loadtxt(
        matched_dge_summary, delimiter="\t", dtype="float", skiprows=7, usecols=2
    )
    for i in range(0, len(matched_dge_summary_transcripts)):
        matched_dge_summary_transcripts[i] = math.log10(
            matched_dge_summary_transcripts[i]
        )
    df = pd.DataFrame({"x": matched_dge_summary_transcripts})
    plt.figure(figsize=(8, 8))
    plt.hist(df["x"], bins=70, facecolor="lightskyblue", edgecolor="black")
    plt.yticks(rotation=90)
    plt.xlabel("log10-based number of UMIs")
    plt.title("Histogram of UMIs per matched barcode")
    plt.savefig(pp, format="pdf")

    matched_dge_summary_genes = np.loadtxt(
        matched_dge_summary, delimiter="\t", dtype="float", skiprows=7, usecols=3
    )
    for i in range(0, len(matched_dge_summary_genes)):
        matched_dge_summary_genes[i] = math.log10(matched_dge_summary_genes[i])
    df = pd.DataFrame({"x": matched_dge_summary_genes})
    plt.figure(figsize=(8, 8))
    plt.hist(df["x"], bins=70, facecolor="lightskyblue", edgecolor="black")
    plt.yticks(rotation=90)
    plt.xlabel("log10-based number of genes")
    plt.title("Histogram of genes per matched barcode")
    plt.savefig(pp, format="pdf")

    file = "{}/{}.barcode_distribution_XC.txt".format(alignment_folder, library)
    mat = np.loadtxt(file, delimiter="\t", dtype="int", skiprows=1)
    cou = mat[0, 1] + mat[0, 2] + mat[0, 3] + mat[0, 4]
    df_x = []
    df_x.extend(mat[:, 0])
    df_x.extend(mat[:, 0])
    df_x.extend(mat[:, 0])
    df_x.extend(mat[:, 0])
    df_y = []
    df_y.extend(mat[:, 1] / cou * 100)
    df_y.extend(mat[:, 2] / cou * 100)
    df_y.extend(mat[:, 3] / cou * 100)
    df_y.extend(mat[:, 4] / cou * 100)
    max_y = int(max(df_y))
    j = len(mat)
    fig, ax = plt.subplots(figsize=(8, 8))
    colors = ["red", "blue", "green", "purple"]
    labels = ["A", "C", "G", "T"]
    for i in range(4):
        ax.scatter(
            df_x[i * j : (i + 1) * j],
            df_y[i * j : (i + 1) * j],
            c=colors[i],
            s=20,
            label=labels[i],
        )
    ax.legend(loc="lower right")
    plt.xlim((0, j + 2))
    plt.ylim((0, max_y + 2))
    plt.yticks(rotation=90)
    plt.xlabel("base position")
    plt.ylabel("fraction of reads")
    plt.title("Matched cell barcodes")
    plt.savefig(pp, format="pdf")

    file = "{}/{}.barcode_distribution_XM.txt".format(alignment_folder, library)
    mat = np.loadtxt(file, delimiter="\t", dtype="int", skiprows=1)
    cou = mat[0, 1] + mat[0, 2] + mat[0, 3] + mat[0, 4]
    df_x = []
    df_x.extend(mat[:, 0])
    df_x.extend(mat[:, 0])
    df_x.extend(mat[:, 0])
    df_x.extend(mat[:, 0])
    df_y = []
    df_y.extend(mat[:, 1] / cou * 100)
    df_y.extend(mat[:, 2] / cou * 100)
    df_y.extend(mat[:, 3] / cou * 100)
    df_y.extend(mat[:, 4] / cou * 100)
    max_y = int(max(df_y))
    j = len(mat)
    fig, ax = plt.subplots(figsize=(8, 8))
    colors = ["red", "blue", "green", "purple"]
    labels = ["A", "C", "G", "T"]
    for i in range(4):
        ax.scatter(
            df_x[i * j : (i + 1) * j],
            df_y[i * j : (i + 1) * j],
            c=colors[i],
            s=20,
            label=labels[i],
        )
    ax.legend(loc="lower right")
    plt.xlim((0, j + 1))
    plt.ylim((0, max_y + 2))
    plt.yticks(rotation=90)
    plt.xlabel("base position")
    plt.ylabel("fraction of reads")
    plt.title("Matched molecular barcodes")
    plt.savefig(pp, format="pdf")

    if gen_read1_plot:
        # Wait for all of filter_unmapped_bam finish

        log.info(
            f"{flowcell_barcode} - Merge filtered unmapped bam files for {library} {reference2}"
        )

        unmapped_bam_file = "{}/{}_unmapped.bam".format(
            barcode_matching_folder, library
        )
        read1_file = "{}/{}.base_distribution_read1.txt".format(
            alignment_folder, library
        )
        commandStr = (
            f"java -Djava.io.tmpdir={tmpdir}"
            " -Dsamjdk.buffer_size=131072 -XX:+UseParallelGC -XX:GCTimeLimit=20"
            " -XX:GCHeapFreeLimit=10 -Xmx8192m "
        )
        commandStr += (
            "-jar "
            + picard_folder
            + "/picard.jar MergeSamFiles TMP_DIR="
            + tmpdir
            + " CREATE_INDEX=true CREATE_MD5_FILE=false VALIDATION_STRINGENCY=SILENT "
        )
        commandStr += (
            "OUTPUT=" + unmapped_bam_file + " SORT_ORDER=coordinate ASSUME_SORTED=true"
        )
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for lane_slice in slice_id[lanes[i]]:
                filtered_bam = f"{barcode_matching_folder}/{library}_{lanes[i]}_{lane_slice}_{barcodes[i]}_filtered.bam"
                if os.path.isfile(filtered_bam):
                    commandStr += " INPUT=" + filtered_bam
                else:
                    log.error(f"{filtered_bam} not found!")
        log.info(f"{flowcell_barcode} MergeSamFiles for {library}")
        log.debug(f"Command = {commandStr}")
        os.system(commandStr)
        log.info(f"{flowcell_barcode} MergeSamFiles for {library} is done.")

        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for lane_slice in slice_id[lanes[i]]:
                filtered_bam = "{}/{}_{}_{}_{}_filtered.bam".format(
                    barcode_matching_folder,
                    library,
                    lanes[i],
                    lane_slice,
                    barcodes[i],
                )
                if os.path.isfile(filtered_bam):
                    call(["rm", filtered_bam])

        log.info(
            f"{flowcell_barcode} - Merge filtered unmapped bam files for "
            + library
            + " "
            + reference2
            + " is done. ",
        )

        # Base distribution at read position for read 1
        commandStr = (
            dropseq_folder
            + "/BaseDistributionAtReadPosition -m 7692m I="
            + unmapped_bam_file
        )
        commandStr += (
            " OUTPUT="
            + read1_file
            + " TMP_DIR="
            + tmpdir
            + " READ_NUMBER=1 VALIDATION_STRINGENCY=SILENT"
        )
        os.system(commandStr)

        if os.path.isfile(unmapped_bam_file):
            call(["rm", unmapped_bam_file])
        if os.path.isfile(
            "{}/{}_unmapped.bai".format(barcode_matching_folder, library)
        ):
            call(
                [
                    "rm",
                    "{}/{}_unmapped.bai".format(barcode_matching_folder, library),
                ]
            )

        read1_len = get_read1_len(bead_structure)  # 42
        mat = np.loadtxt(read1_file, delimiter="\t", dtype="int", skiprows=1)
        cou = mat[0, 1] + mat[0, 2] + mat[0, 3] + mat[0, 4]
        df_x = []
        df_x.extend(mat[:, 0])
        df_x.extend(mat[:, 0])
        df_x.extend(mat[:, 0])
        df_x.extend(mat[:, 0])
        df_y = []
        df_y.extend(mat[:, 1] / cou * 100)
        df_y.extend(mat[:, 2] / cou * 100)
        df_y.extend(mat[:, 3] / cou * 100)
        df_y.extend(mat[:, 4] / cou * 100)
        max_y = int(max(df_y))
        fig, ax = plt.subplots(figsize=(8, 8))
        colors = ["red", "blue", "green", "purple"]
        labels = ["A", "C", "G", "T"]
        for i in range(4):
            ax.scatter(
                df_x[i * read1_len : (i + 1) * read1_len],
                df_y[i * read1_len : (i + 1) * read1_len],
                c=colors[i],
                s=20,
                label=labels[i],
            )
        ax.legend(loc="lower right")
        plt.xlim((0, read1_len + 2))
        plt.ylim((0, max_y + 2))
        plt.yticks(rotation=90)
        plt.xlabel("base position")
        plt.ylabel("fraction of reads")
        plt.title("Read 1 for matched barcodes")
        plt.savefig(pp, format="pdf")

    xs = []
    ys = []
    zs = []
    matched_dge_summary_transcripts = np.loadtxt(
        matched_dge_summary, delimiter="\t", dtype="float", skiprows=7, usecols=2
    )
    k = int(len(matched_dge_summary_transcripts) / 20)
    for i in range(k):
        matched_dge_summary_transcripts[i] = matched_dge_summary_transcripts[k]
    for i in range(len(matched_dge_summary_barcodes)):
        barcode = matched_dge_summary_barcodes[i]
        if barcode in matched_barcodes:
            zs.append(matched_dge_summary_transcripts[i])
            xs.append(coordinatesx[np.where(matched_barcodes == barcode)])
            ys.append(coordinatesy[np.where(matched_barcodes == barcode)])

    j = len(zs)
    for i in range(j - 1):
        for j in range(i + 1, j):
            if zs[i] > zs[j]:
                xs[i], xs[j] = xs[j], xs[i]
                ys[i], ys[j] = ys[j], ys[i]
                zs[i], zs[j] = zs[j], zs[i]

    df = pd.DataFrame({"x": xs, "y": ys, "z": zs})
    plt.figure(figsize=(8, 8))
    plt.set_cmap("viridis_r")
    plt.scatter(df["x"], df["y"], c=df["z"], s=0.5)
    plt.colorbar()
    plt.axis("equal")
    plt.yticks(rotation=90)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("total # UMIs per matched bead (95% perecntile)")
    plt.savefig(pp, format="pdf")

    f1 = "{}/{}.fracIntronicExonicPerCell.txt".format(alignment_folder, library)
    f2 = "{}/{}.fracIntronicExonicPerCell.txt.gz".format(alignment_folder, library)
    if not os.path.isfile(f1):
        with gzip.open(f2, "rb") as f_in:
            with open(f1, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    # PCT_RIBOSOMAL_BASES	PCT_CODING_BASES	PCT_UTR_BASES	PCT_INTRONIC_BASES	PCT_INTERGENIC_BASES
    # 17	18	19	20	21
    # mat[1] += mat[2]
    # mat[2] = mat[1] + mat[3]
    # df_x = ['ribosomal','exonic','genic','intronic','intergenic']
    pct_mt = np.loadtxt(f1, delimiter="\t", dtype="float", skiprows=7, usecols=1)
    pct_ribosomal = np.loadtxt(
        f1, delimiter="\t", dtype="float", skiprows=7, usecols=17
    )
    pct_coding = np.loadtxt(f1, delimiter="\t", dtype="float", skiprows=7, usecols=18)
    pct_utr = np.loadtxt(f1, delimiter="\t", dtype="float", skiprows=7, usecols=19)
    barcodes = np.loadtxt(f1, delimiter="\t", dtype="str", skiprows=7, usecols=29)
    xs = []
    ys = []
    zs = []
    us = []
    vs = []
    for i in range(len(barcodes)):
        barcode = barcodes[i]
        if barcode in matched_barcodes:
            xs.append(coordinatesx[np.where(matched_barcodes == barcode)])
            ys.append(coordinatesy[np.where(matched_barcodes == barcode)])
            zs.append(pct_mt[i])
            us.append(pct_coding[i] + pct_utr[i])
            vs.append(pct_ribosomal[i])

    j = len(zs)
    for i in range(j - 1):
        for j in range(i + 1, j):
            if zs[i] > zs[j]:
                xs[i], xs[j] = xs[j], xs[i]
                ys[i], ys[j] = ys[j], ys[i]
                zs[i], zs[j] = zs[j], zs[i]
                us[i], us[j] = us[j], us[i]
                vs[i], vs[j] = vs[j], vs[i]

    df = pd.DataFrame({"x": xs, "y": ys, "z": zs})
    plt.figure(figsize=(8, 8))
    plt.set_cmap("viridis_r")
    plt.scatter(df["x"], df["y"], c=df["z"], s=0.5)
    plt.colorbar()
    plt.axis("equal")
    plt.yticks(rotation=90)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("% mitochondrial reads per matched bead")
    plt.savefig(pp, format="pdf")

    for i in range(j - 1):
        for j in range(i + 1, j):
            if us[i] > us[j]:
                xs[i], xs[j] = xs[j], xs[i]
                ys[i], ys[j] = ys[j], ys[i]
                us[i], us[j] = us[j], us[i]
                vs[i], vs[j] = vs[j], vs[i]

    df = pd.DataFrame({"x": xs, "y": ys, "z": us})
    plt.figure(figsize=(8, 8))
    plt.set_cmap("viridis_r")
    plt.scatter(df["x"], df["y"], c=df["z"], s=0.5)
    plt.colorbar()
    plt.axis("equal")
    plt.yticks(rotation=90)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("% exonic reads per matched bead")
    plt.savefig(pp, format="pdf")

    for i in range(j - 1):
        for j in range(i + 1, j):
            if vs[i] > vs[j]:
                xs[i], xs[j] = xs[j], xs[i]
                ys[i], ys[j] = ys[j], ys[i]
                vs[i], vs[j] = vs[j], vs[i]

    df = pd.DataFrame({"x": xs, "y": ys, "z": vs})
    plt.figure(figsize=(8, 8))
    plt.set_cmap("viridis_r")
    plt.scatter(df["x"], df["y"], c=df["z"], s=0.5)
    plt.colorbar()
    plt.axis("equal")
    plt.yticks(rotation=90)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("% ribosomal reads per matched bead")
    plt.savefig(pp, format="pdf")

    pp.close()

    if os.path.isfile(f1):
        call(["rm", f1])

    combined_cmatcher_file = "{}/{}_barcode_matching_distance.txt".format(
        barcode_matching_folder, library
    )
    combined_cmatcher_shuffled_file = (
        "{}/{}_barcode_matching_distance_shuffled.txt".format(
            barcode_matching_folder, library
        )
    )
    distances = np.loadtxt(
        combined_cmatcher_file, delimiter="\t", dtype="int", skiprows=1, usecols=3
    )
    shuffdists = np.loadtxt(
        combined_cmatcher_shuffled_file,
        delimiter="\t",
        dtype="int",
        skiprows=1,
        usecols=3,
    )
    if len(distances) > 100000:
        distances = distances[:100000]
        shuffdists = shuffdists[:100000]
    df = pd.DataFrame({"Real": distances, "Shuffled": shuffdists})
    df = pd.melt(df)
    p = (
        ggplot(aes(x="value", fill="variable"), data=df)
        + geom_histogram(alpha=0.3, binwidth=1, position="identity")
        + scale_x_continuous(
            limits=(-1, 15), breaks=np.arange(0, 16, 2), labels=np.arange(0, 16, 2)
        )
        + xlab("hamming distance")
        + ylab("# beads")
        + ggtitle("hamming distance (real vs shuffled)")
    )
    ggsave(
        plot=p,
        height=10,
        width=10,
        filename=library + "_" + reference2 + "_shuffled.pdf",
        path=alignment_folder,
        verbose=False,
    )

    log.info(
        f"{flowcell_barcode} - Generate plots for matched barcodes for "
        + library
        + " "
        + reference2
        + " is done. ",
    )

    # Call gen_readme
    output_file = "{}/logs/gen_readme_{}_{}.log".format(
        output_folder, library, reference2
    )
    submission_script = "{}/gen_readme.sh".format(scripts_folder)
    call_args = [
        "qsub",
        "-o",
        output_file,
        "-l",
        "h_vmem=2g",
        "-notify",
        "-l",
        "h_rt=1:0:0",
        "-j",
        "y",
        "-P",
        "macosko_lab",
        "-l",
        "os=RedHat7",
        submission_script,
        manifest_file,
        library,
        scripts_folder,
        locus_function_list,
        output_folder,
        alignment_folder,
    ]
    call(call_args)


if __name__ == "__main__":
    main()
