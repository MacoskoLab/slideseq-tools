#!/usr/bin/python

# This script is to combine outputs from cmatcher_beads

import csv
import os
import sys
import time
import traceback

# silence warnings for pandas below
from datetime import datetime
from subprocess import call


# Write to log file
def write_log(log_file, flowcell_barcode, log_string):
    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as logfile:
        logfile.write(
            dt_string
            + " [Slide-seq Flowcell Alignment Workflow - "
            + flowcell_barcode
            + "]: "
            + log_string
            + "\n"
        )


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
    email_address = ""
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
                email_address = row[row0.index("email")]
                experiment_date = row[row0.index("date")]

    log_file = "{}/logs/workflow.log".format(output_folder)

    analysis_folder = "{}/{}_{}".format(library_folder, experiment_date, library)
    bead_barcode_file = "{}/BeadBarcodes.txt".format(analysis_folder)
    bead_location_file = "{}/BeadLocations.txt".format(analysis_folder)

    if not os.path.isfile(bead_barcode_file):
        write_log(
            log_file,
            flowcell_barcode,
            "run_cmatcher_beads_combine error: "
            + bead_barcode_file
            + " does not exist!",
        )
        raise Exception(
            "run_cmatcher_beads_combine error: "
            + bead_barcode_file
            + " does not exist!"
        )

    folder_running = "{}/status/running.cmatcher_beads_combine_{}".format(
        output_folder, library
    )
    folder_finished = "{}/status/finished.cmatcher_beads_combine_{}".format(
        output_folder, library
    )
    folder_failed = "{}/status/failed.cmatcher_beads_combine_{}".format(
        output_folder, library
    )

    try:
        call(["mkdir", "-p", folder_running])

        with open(bead_barcode_file, "r") as fin:
            j = sum(1 for _ in fin)

        k = 10000
        ls = j // k

        while 1:
            f = True
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = "{}/{}_barcode_matching_01_{}.finished".format(
                    analysis_folder, library, str(i + 1)
                )
                if not os.path.isfile(file2):
                    f = False
                    break
            if f:
                break
            time.sleep(30)

        print("combine cmatcher_beads outputs...")
        write_log(
            log_file, flowcell_barcode, "Combine cmatcher_beads outputs for " + library
        )
        combined_cmatcher_file = "{}/{}_barcode_matching_01.txt".format(
            analysis_folder, library
        )
        with open(combined_cmatcher_file, "w") as fout:
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = "{}/{}_barcode_matching_01_{}.txt".format(
                    analysis_folder, library, str(i + 1)
                )
                with open(file2, "r") as fin:
                    for line in fin:
                        fout.write(line)

        combined_cmatcher_file2 = "{}/{}_barcode_matching_2.txt".format(
            analysis_folder, library
        )
        with open(combined_cmatcher_file2, "w") as fout:
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = "{}/{}_barcode_matching_2_{}.txt".format(
                    analysis_folder, library, str(i + 1)
                )
                with open(file2, "r") as fin:
                    for line in fin:
                        fout.write(line)

        for i in range(ls + 1):
            if i * k >= j:
                break
            file1 = "{}/{}_barcode_matching_01_{}.txt".format(
                analysis_folder, library, str(i + 1)
            )
            file2 = "{}/{}_barcode_matching_01_{}.finished".format(
                analysis_folder, library, str(i + 1)
            )
            file3 = "{}/BeadBarcodes_{}.txt".format(analysis_folder, str(i + 1))
            file4 = "{}/{}_barcode_matching_2_{}.txt".format(
                analysis_folder, library, str(i + 1)
            )
            if os.path.isfile(file1):
                call(["rm", file1])
            if os.path.isfile(file2):
                call(["rm", file2])
            if os.path.isfile(file3):
                call(["rm", file3])
            if os.path.isfile(file4):
                call(["rm", file4])

        write_log(
            log_file,
            flowcell_barcode,
            "Combine cmatcher_beads outputs for " + library + " is done. ",
        )

        # Create degenerate bead barcodes
        print("Create degenerate bead barcodes...")
        combined_cmatcher_file3 = "{}/BeadBarcodes_degenerate.txt".format(
            analysis_folder
        )
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

        file = "{}/BeadBarcodes_degenerate.finished".format(analysis_folder)
        with open(file, "w") as fout:
            fout.write("finished")

        write_log(
            log_file,
            flowcell_barcode,
            "Create degenerate bead barcodes for " + library + " is done. ",
        )

        call(["mv", folder_running, folder_finished])
    except Exception as exp:
        print("EXCEPTION:!")
        print(exp)
        traceback.print_tb(exp.__traceback__, file=sys.stdout)
        if os.path.isdir(folder_running):
            call(["mv", folder_running, folder_failed])
        else:
            call(["mkdir", "-p", folder_failed])

        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = (
                "The Slide-seq workflow for "
                + library
                + " failed at the step of running cmatcher beads combine. Please check the log file for the issues. "
            )
            call_args = [
                "python",
                "{}/send_email.py".format(scripts_folder),
                email_address,
                subject,
                content,
            ]
            call(call_args)

        sys.exit()


if __name__ == "__main__":
    main()
