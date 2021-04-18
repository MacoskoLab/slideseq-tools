#!/usr/bin/python

# This script is to tag bam using matched bead barcodes

import csv
import logging
import os
import sys
from subprocess import call

import numpy as np

log = logging.getLogger(__name__)


def main():
    if len(sys.argv) != 7:
        print(
            "Please provide six arguments: "
            "manifest file, library ID, lane ID, slice ID, sample barcode and locus function list!"
        )
        sys.exit()

    manifest_file = sys.argv[1]
    library = sys.argv[2]
    lane = sys.argv[3]
    lane_slice = sys.argv[4]
    barcode = sys.argv[5]
    locus_function_list = sys.argv[6]

    # Check if the manifest file exists
    if not os.path.isfile(manifest_file):
        print(f"File {manifest_file} does not exist. Exiting...")
        sys.exit()

    # Read manifest file
    options = {}
    with open(manifest_file, "r") as fp:
        for line in fp:
            key, value = line.rstrip().split("=")
            options[key] = value

    output_folder = options["output_folder"]
    flowcell_barcode = options["flowcell_barcode"]
    library_folder = (
        options["library_folder"]
        if "library_folder" in options
        else f"{output_folder}/libraries"
    )

    # Read info from metadata file
    reference = ""
    experiment_date = ""
    base_quality = "10"
    min_transcripts_per_cell = "10"
    with open("{}/parsed_metadata.txt".format(output_folder), "r") as fin:
        reader = csv.reader(fin, delimiter="\t")
        rows = list(reader)
        row0 = rows[0]
        for i in range(1, len(rows)):
            row = rows[i]
            if row[row0.index("library")] == library:
                reference = row[row0.index("reference")]
                experiment_date = row[row0.index("date")]
                base_quality = row[row0.index("base_quality")]
                min_transcripts_per_cell = row[row0.index("min_transcripts_per_cell")]
                break

    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]
    reference2 = f"{referencePure}.{locus_function_list}"

    alignment_folder = (
        f"{library_folder}/{experiment_date}_{library}/{reference2}/alignment/"
    )
    barcode_matching_folder = (
        f"{library_folder}/{experiment_date}_{library}/{reference2}/barcode_matching"
    )
    prefix_libraries = (
        f"{barcode_matching_folder}/{flowcell_barcode}.{lane}.{lane_slice}.{library}"
    )
    if barcode:
        prefix_libraries += "." + barcode
    mapped_bam = prefix_libraries + ".star_gene_exon_tagged2.bam"
    mapped_sam = (
        f"{barcode_matching_folder}/{library}_{lane}_{lane_slice}_{barcode}_aligned.sam"
    )
    tagged_sam = (
        f"{barcode_matching_folder}/{library}_{lane}_{lane_slice}_{barcode}_tagged.sam"
    )
    tagged_bam = (
        f"{barcode_matching_folder}/{library}_{lane}_{lane_slice}_{barcode}_tagged.bam"
    )
    raw_sam = (
        f"{barcode_matching_folder}/{library}_{lane}_{lane_slice}_{barcode}_raw.sam"
    )
    raw_bam = (
        f"{barcode_matching_folder}/{library}_{lane}_{lane_slice}_{barcode}_raw.bam"
    )
    shuffled_sam = f"{barcode_matching_folder}/{library}_{lane}_{lane_slice}_{barcode}_shuffled.sam"
    shuffled_bam = f"{barcode_matching_folder}/{library}_{lane}_{lane_slice}_{barcode}_shuffled.bam"
    combined_cmatcher_file = f"{barcode_matching_folder}/{library}_barcode_matching.txt"
    combined_cmatcher_shuffled_file = (
        f"{barcode_matching_folder}/{library}_barcode_matching_shuffled.txt"
    )

    if not os.path.isfile(mapped_bam):
        log.error(
            f"{flowcell_barcode} - TagMatchedBam error: {mapped_bam} does not exist!"
        )
        raise Exception(
            f"{flowcell_barcode} - TagMatchedBam error: {mapped_bam} does not exist!"
        )

    if not os.path.isfile(combined_cmatcher_file):
        log.error(
            f"{flowcell_barcode} - TagMatchedBam error: {combined_cmatcher_file} does not exist!"
        )
        raise Exception(
            f"{flowcell_barcode} - TagMatchedBam error: {combined_cmatcher_file} does not exist!"
        )

    if not os.path.isfile(combined_cmatcher_shuffled_file):
        log.error(
            f"{flowcell_barcode} - TagMatchedBam error: {combined_cmatcher_shuffled_file} does not exist!"
        )
        raise Exception(
            f"{flowcell_barcode} - TagMatchedBam error: {combined_cmatcher_shuffled_file} does not exist!"
        )

    log.info(
        f"{flowcell_barcode} - Tag matched bam for {library} {reference2} in lane {lane} slice {lane_slice}"
    )

    log.info("mapped bam to sam")
    commandStr = "samtools view -h -o " + mapped_sam + " " + mapped_bam
    os.system(commandStr)

    call(["rm", mapped_bam])

    log.info("read combined_cmatcher_file into dict1")
    dict1 = {}
    with open(combined_cmatcher_file, "r") as fin:
        j = 0
        for line in fin:
            j += 1
            if j > 1:
                dict1[line.split("\t")[0]] = line.split("\t")[2]

    log.info("read raw2shuffle into dict2")
    dict2 = {}
    select_cell_file = f"{alignment_folder}{library}.{min_transcripts_per_cell}_transcripts_mq_{base_quality}_selected_cells.txt"
    select_cell_shuffled_file = f"{alignment_folder}{library}.{min_transcripts_per_cell}_transcripts_mq_{base_quality}_selected_cells.shuffled.txt"
    bc1 = np.loadtxt(select_cell_file, delimiter="\t", dtype="str", usecols=0)
    bc2 = np.loadtxt(select_cell_shuffled_file, delimiter="\t", dtype="str", usecols=0)
    for i in range(len(bc1)):
        dict2[bc1[i].strip("\n")] = bc2[i].strip("\n")

    log.info("read combined_cmatcher_shuffled_file into dict3")
    dict3 = {}
    with open(combined_cmatcher_shuffled_file, "r") as fin:
        j = 0
        for line in fin:
            j += 1
            if j > 1:
                dict3[line.split("\t")[0]] = line.split("\t")[2]

    log.info("gen tagged_sam")
    with open(tagged_sam, "w") as fout1:
        with open(raw_sam, "w") as fout2:
            with open(shuffled_sam, "w") as fout3:
                with open(mapped_sam, "r") as fin:
                    for line in fin:
                        if line[0] == "@":
                            fout1.write(line)
                            fout2.write(line)
                            fout3.write(line)
                        else:
                            items1 = line.split("\t")
                            bc1 = items1[11]
                            items2 = bc1.split(":")
                            bc2 = items2[2]
                            if bc2 in dict1:
                                fout2.write(line)
                                items2[2] = dict1[bc2]
                                items1[11] = ":".join(items2)
                                fout1.write("\t".join(items1))
                            if (bc2 in dict2) and (dict2[bc2] in dict3):
                                items2[2] = dict2[bc2]
                                items1[11] = ":".join(items2)
                                fout3.write("\t".join(items1))

    if os.path.isfile(mapped_sam):
        call(["rm", mapped_sam])

    commandStr = "samtools view -S -b " + tagged_sam + " > " + tagged_bam
    os.system(commandStr)

    if os.path.isfile(tagged_sam):
        call(["rm", tagged_sam])

    commandStr = "samtools view -S -b " + raw_sam + " > " + raw_bam
    os.system(commandStr)

    if os.path.isfile(raw_sam):
        call(["rm", raw_sam])

    commandStr = "samtools view -S -b " + shuffled_sam + " > " + shuffled_bam
    os.system(commandStr)

    if os.path.isfile(shuffled_sam):
        call(["rm", shuffled_sam])

    log.info(
        f"{flowcell_barcode} - Tag matched bam for {library} {reference2} in lane {lane} slice {lane_slice} is done"
    )


if __name__ == "__main__":
    main()
