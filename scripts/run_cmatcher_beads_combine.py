#!/usr/bin/python

# This script is to combine outputs from cmatcher_beads

import csv
import logging
import os
import sys
import time
from subprocess import call

log = logging.getLogger(__name__)


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
    with open(manifest_file, "r") as fp:
        for line in fp:
            key, value = line.rstrip().split("=")
            options[key] = value

    output_folder = options["output_folder"]
    flowcell_barcode = options["flowcell_barcode"]

    library_folder = (
        options["library_folder"]
        if "library_folder" in options
        else "{}/libraries".format(output_folder)
    )
    scripts_folder = (
        options["scripts_folder"]
        if "scripts_folder" in options
        else "/broad/macosko/jilong/slideseq_pipeline/scripts"
    )

    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    experiment_date = ""
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
                experiment_date = row[row0.index("date")]

    analysis_folder = "{}/{}_{}".format(library_folder, experiment_date, library)
    bead_barcode_file = "{}/BeadBarcodes.txt".format(analysis_folder)
    bead_location_file = "{}/BeadLocations.txt".format(analysis_folder)

    if not os.path.isfile(bead_barcode_file):
        log.error(
            f"{flowcell_barcode} - TagMatchedBam error: {bead_barcode_file} does not exist!",
        )
        raise FileNotFoundError(
            f"TagMatchedBam error: {bead_barcode_file} does not exist!"
        )

    with open(bead_barcode_file, "r") as fin:
        j = sum(1 for _ in fin)

    k = 10000
    ls = j // k

    while 1:
        f = True
        for i in range(ls + 1):
            if i * k >= j:
                break
            file2 = (
                f"{analysis_folder}/{library}_barcode_matching_01_{str(i + 1)}.finished"
            )
            if not os.path.isfile(file2):
                f = False
                break
        if f:
            break
        time.sleep(30)

    log.info("combine cmatcher_beads outputs...")
    combined_cmatcher_file = f"{analysis_folder}/{library}_barcode_matching_01.txt"
    with open(combined_cmatcher_file, "w") as fout:
        for i in range(ls + 1):
            if i * k >= j:
                break
            file2 = f"{analysis_folder}/{library}_barcode_matching_01_{str(i + 1)}.txt"
            with open(file2, "r") as fin:
                for line in fin:
                    fout.write(line)

    combined_cmatcher_file2 = f"{analysis_folder}/{library}_barcode_matching_2.txt"
    with open(combined_cmatcher_file2, "w") as fout:
        for i in range(ls + 1):
            if i * k >= j:
                break
            file2 = f"{analysis_folder}/{library}_barcode_matching_2_{str(i + 1)}.txt"
            with open(file2, "r") as fin:
                for line in fin:
                    fout.write(line)

    for i in range(ls + 1):
        if i * k >= j:
            break
        file1 = f"{analysis_folder}/{library}_barcode_matching_01_{str(i + 1)}.txt"
        file2 = f"{analysis_folder}/{library}_barcode_matching_01_{str(i + 1)}.finished"
        file3 = f"{analysis_folder}/BeadBarcodes_{str(i + 1)}.txt"
        file4 = f"{analysis_folder}/{library}_barcode_matching_2_{str(i + 1)}.txt"
        if os.path.isfile(file1):
            call(["rm", file1])
        if os.path.isfile(file2):
            call(["rm", file2])
        if os.path.isfile(file3):
            call(["rm", file3])
        if os.path.isfile(file4):
            call(["rm", file4])

    log.info(
        f"{flowcell_barcode} - Combine cmatcher_beads outputs for {library} is done."
    )

    # Create degenerate bead barcodes
    log.info("Create degenerate bead barcodes...")
    combined_cmatcher_file3 = "{}/BeadBarcodes_degenerate.txt".format(analysis_folder)
    commandStr = (
        scripts_folder
        + "/degenerate_beads "
        + bead_barcode_file
        + " "
        + bead_location_file
        + " "
        + combined_cmatcher_file
        + " "
        + combined_cmatcher_file2
        + " "
        + combined_cmatcher_file3
    )
    os.system(commandStr)

    file = f"{analysis_folder}/BeadBarcodes_degenerate.finished"
    with open(file, "w") as fout:
        fout.write("finished")

    log.info(
        f"{flowcell_barcode} - Create degenerate bead barcodes for {library} is done."
    )


if __name__ == "__main__":
    main()
