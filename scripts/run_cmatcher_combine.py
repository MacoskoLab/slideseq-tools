#!/usr/bin/python

# This script is to combine outputs from cmatcher
# and call tag_matched_bam

import csv
import os
import sys
import time
import traceback
from datetime import datetime
from subprocess import call

import numpy as np
from new_submit_to_taskrunner import call_to_taskrunner


# Get tile information from RunInfo.xml
def get_tiles(x, lane):
    tiles = []
    with open(x, "r") as fin:
        for line in fin:
            line = line.strip(" \t\n")
            if line.startswith("<Tile>", 0):
                lane_split = line[6:].split("<")[0]
                if lane_split.split("_")[0] == lane:
                    tiles.append(lane_split.split("_")[1])

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
        logfile.write(
            dt_string
            + " [Slide-seq Flowcell Alignment Workflow - "
            + flowcell_barcode
            + "]: "
            + log_string
            + "\n"
        )


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
    email_address = ""
    experiment_date = ""
    gen_read1_plot = False
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
                email_address = row[row0.index("email")]
                experiment_date = row[row0.index("date")]
                if "gen_read1_plot" in row0:
                    gen_read1_plot = str2bool(row[row0.index("gen_read1_plot")])

    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]
    reference2 = referencePure + "." + locus_function_list

    runinfo_file = "{}/RunInfo.xml".format(flowcell_directory)
    log_file = "{}/logs/workflow.log".format(output_folder)

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

    analysis_folder = "{}/{}_{}".format(library_folder, experiment_date, library)
    alignment_folder = "{}/{}/alignment/".format(analysis_folder, reference2)
    barcode_matching_folder = "{}/{}/barcode_matching/".format(
        analysis_folder, reference2
    )
    select_cell_file = (
        alignment_folder
        + library
        + "."
        + min_transcripts_per_cell
        + "_transcripts_mq_"
        + base_quality
        + "_selected_cells.txt"
    )

    if not os.path.isfile(select_cell_file):
        write_log(
            log_file,
            flowcell_barcode,
            "run_cmatcher_combine error: " + select_cell_file + " does not exist!",
        )
        raise Exception(
            "run_cmatcher_combine error: " + select_cell_file + " does not exist!"
        )

    folder_running = "{}/status/running.cmatcher_combine_{}_{}".format(
        output_folder, library, reference2
    )
    folder_finished = "{}/status/finished.cmatcher_combine_{}_{}".format(
        output_folder, library, reference2
    )
    folder_failed = "{}/status/failed.cmatcher_combine_{}_{}".format(
        output_folder, library, reference2
    )

    try:
        call(["mkdir", "-p", folder_running])

        with open(select_cell_file, "r") as fin:
            j = sum(1 for _ in fin)

        k = 10000
        ls = j // k

        print("# selected cells: " + str(j))

        while 1:
            f = True
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = "{}/{}_barcode_matching_{}.finished".format(
                    barcode_matching_folder, library, str(i + 1)
                )
                if not os.path.isfile(file2):
                    f = False
                    break
            if f:
                break
            time.sleep(30)

        while 1:
            f = True
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = "{}/{}_barcode_matching_shuffled_{}.finished".format(
                    barcode_matching_folder, library, str(i + 1)
                )
                if not os.path.isfile(file2):
                    f = False
                    break
            if f:
                break
            time.sleep(30)

        print("combine cmatcher outputs...")
        write_log(
            log_file,
            flowcell_barcode,
            "Combine CMatcher outputs for " + library + " " + reference2,
        )
        combined_cmatcher_file = "{}/{}_barcode_matching.txt".format(
            barcode_matching_folder, library
        )
        with open(combined_cmatcher_file, "w") as fout:
            fout.write(
                "IlluminaBarcodes\tProcessedIlluminaBarcodes\tBeadBarcodes\tDistance\tX\tY\n"
            )
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = "{}/{}_barcode_matching_{}.txt".format(
                    barcode_matching_folder, library, str(i + 1)
                )
                with open(file2, "r") as fin:
                    j = 0
                    for line in fin:
                        j += 1
                        if j > 1:
                            fout.write(line)

        # Combine CMatcher logs
        combined_cmatcher_summary = "{}/{}_barcode_matching_summary.txt".format(
            barcode_matching_folder, library
        )
        total = 0
        unique = 0
        multi = 0
        for i in range(ls + 1):
            if i * k >= j:
                break
            file2 = "{}/{}_barcode_matching_{}.txt.log".format(
                barcode_matching_folder, library, str(i + 1)
            )
            if not os.path.isfile(file2):
                continue
            j = 0
            with open(file2, "r") as fin:
                for line in fin:
                    j += 1
                    s = line.split(":")[1]
                    s = s.strip(" \t\n")
                    if j == 1:
                        total += int(s)
                    elif j == 2:
                        unique += int(s)
                    elif j == 3:
                        multi += int(s)

        with open(combined_cmatcher_summary, "w") as fout:
            fout.write("Total # barcodes: {}\n".format(str(total)))
            fout.write(
                "# unique matched barcodes: {}, {}%\n".format(
                    str(unique), str(unique * 100 / total)
                )
            )
            fout.write(
                "# multiple matched barcodes: {}, {}%\n".format(
                    str(multi), str(multi * 100 / total)
                )
            )

        for i in range(ls + 1):
            if i * k >= j:
                break
            file1 = "{}/{}_barcode_matching_{}.txt".format(
                barcode_matching_folder, library, str(i + 1)
            )
            file2 = "{}/{}_barcode_matching_{}.finished".format(
                barcode_matching_folder, library, str(i + 1)
            )
            name = (
                library
                + "."
                + min_transcripts_per_cell
                + "_transcripts_mq_"
                + base_quality
                + "_selected_cells"
            )
            file3 = "{}/{}_{}.txt".format(alignment_folder, name, str(i + 1))
            file4 = "{}/{}_barcode_matching_{}.txt.log".format(
                barcode_matching_folder, library, str(i + 1)
            )
            if os.path.isfile(file1):
                call(["rm", file1])
            if os.path.isfile(file2):
                call(["rm", file2])
            if os.path.isfile(file3):
                call(["rm", file3])
            if os.path.isfile(file4):
                call(["rm", file4])

        combined_cmatcher_file2 = "{}/{}_barcode_matching_distance.txt".format(
            barcode_matching_folder, library
        )
        with open(combined_cmatcher_file2, "w") as fout:
            fout.write(
                "IlluminaBarcodes\tProcessedIlluminaBarcodes\tBeadBarcodes\tDistance\tX\tY\n"
            )
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = "{}/{}_barcode_matching_distance_{}.txt".format(
                    barcode_matching_folder, library, str(i + 1)
                )
                with open(file2, "r") as fin:
                    j = 0
                    for line in fin:
                        j += 1
                        if j > 1:
                            fout.write(line)
                call(["rm", file2])

        # UniqueMappedIlluminaBarcodes
        bci = np.loadtxt(
            combined_cmatcher_file, delimiter="\t", dtype="str", skiprows=1, usecols=1
        )
        bci = np.unique(bci)
        unique_bci_file = "{}/{}_unique_matched_illumina_barcodes.txt".format(
            barcode_matching_folder, library
        )
        with open(unique_bci_file, "w") as f1:
            for bc in bci:
                f1.write("%s\n" % bc)

        os.system("gzip -c " + unique_bci_file + " > " + unique_bci_file + ".gz")

        write_log(
            log_file,
            flowcell_barcode,
            "Combine CMatcher outputs for " + library + " " + reference2 + " is done. ",
        )

        # Get unique matched bead barcodes and locations
        print("Get unique matched bead barcodes and locations...")
        write_log(
            log_file,
            flowcell_barcode,
            "Get unique matched bead barcodes and locations for "
            + library
            + " "
            + reference2,
        )
        bead_dict = {}
        matched_bead_barcode_file = "{}/{}_matched_bead_barcodes.txt".format(
            barcode_matching_folder, library
        )
        matched_bead_location_file = "{}/{}_matched_bead_locations.txt".format(
            barcode_matching_folder, library
        )
        bead_location_forR = "{}/BeadLocationsForR.csv".format(barcode_matching_folder)
        with open(matched_bead_barcode_file, "w") as fout1:
            with open(matched_bead_location_file, "w") as fout2:
                with open(bead_location_forR, "w") as fout3:
                    fout3.write("barcodes,xcoord,ycoord\n")
                    with open(combined_cmatcher_file, "r") as fin:
                        j = 0
                        for line in fin:
                            j += 1
                            if j > 1:
                                bc = line.split("\t")[2]
                                dist = line.split("\t")[3]
                                x = line.split("\t")[4]
                                y = line.split("\t")[5]
                                if bc not in bead_dict:
                                    fout1.write(bc + "\n")
                                    fout2.write(dist + "\t" + x + "\t" + y)
                                    fout3.write(bc + "," + x + "," + y)
                                    bead_dict[bc] = 1

        write_log(
            log_file,
            flowcell_barcode,
            "Get unique matched bead barcodes and locations for "
            + library
            + " "
            + reference2
            + " is done. ",
        )

        combined_cmatcher_file = "{}/{}_barcode_matching_shuffled.txt".format(
            barcode_matching_folder, library
        )
        with open(combined_cmatcher_file, "w") as fout:
            fout.write(
                "IlluminaBarcodes\tProcessedIlluminaBarcodes\tBeadBarcodes\tDistance\tX\tY\n"
            )
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = "{}/{}_barcode_matching_shuffled_{}.txt".format(
                    barcode_matching_folder, library, str(i + 1)
                )
                with open(file2, "r") as fin:
                    j = 0
                    for line in fin:
                        j += 1
                        if j > 1:
                            fout.write(line)

        for i in range(ls + 1):
            if i * k >= j:
                break
            file1 = "{}/{}_barcode_matching_shuffled_{}.txt".format(
                barcode_matching_folder, library, str(i + 1)
            )
            file2 = "{}/{}_barcode_matching_shuffled_{}.finished".format(
                barcode_matching_folder, library, str(i + 1)
            )
            name = (
                library
                + "."
                + min_transcripts_per_cell
                + "_transcripts_mq_"
                + base_quality
                + "_selected_cells.shuffled"
            )
            file3 = "{}/{}_{}.txt".format(alignment_folder, name, str(i + 1))
            file4 = "{}/{}_barcode_matching_shuffled_{}.txt.log".format(
                barcode_matching_folder, library, str(i + 1)
            )
            if os.path.isfile(file1):
                call(["rm", file1])
            if os.path.isfile(file2):
                call(["rm", file2])
            if os.path.isfile(file3):
                call(["rm", file3])
            if os.path.isfile(file4):
                call(["rm", file4])

        combined_cmatcher_file2 = "{}/{}_barcode_matching_distance_shuffled.txt".format(
            barcode_matching_folder, library
        )
        with open(combined_cmatcher_file2, "w") as fout:
            fout.write(
                "IlluminaBarcodes\tProcessedIlluminaBarcodes\tBeadBarcodes\tDistance\tX\tY\n"
            )
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = "{}/{}_barcode_matching_distance_shuffled_{}.txt".format(
                    barcode_matching_folder, library, str(i + 1)
                )
                with open(file2, "r") as fin:
                    j = 0
                    for line in fin:
                        j += 1
                        if j > 1:
                            fout.write(line)

                call(["rm", file2])

        # UniqueMappedIlluminaBarcodes
        bci = np.loadtxt(
            combined_cmatcher_file, delimiter="\t", dtype="str", skiprows=1, usecols=1
        )
        bci = np.unique(bci)
        shuffled_bci_file = "{}/{}_unique_shuffled_illumina_barcodes.txt".format(
            barcode_matching_folder, library
        )
        with open(shuffled_bci_file, "w") as f1:
            for bc in bci:
                f1.write("%s\n" % bc)

        os.system("gzip -c " + shuffled_bci_file + " > " + shuffled_bci_file + ".gz")

        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for lane_slice in slice_id[lanes[i]]:
                # Call tag_matched_bam
                output_file = "{}/logs/tag_matched_bam_{}_{}_{}_{}_{}.log".format(
                    output_folder,
                    library,
                    lanes[i],
                    lane_slice,
                    barcodes[i],
                    reference2,
                )
                submission_script = "{}/tag_matched_bam.sh".format(scripts_folder)
                call_args = [
                    "qsub",
                    "-o",
                    output_file,
                    "-l",
                    "h_vmem=10g",
                    "-notify",
                    "-l",
                    "h_rt=10:0:0",
                    "-j",
                    "y",
                    "-P",
                    "macosko_lab",
                    "-l",
                    "os=RedHat7",
                    submission_script,
                    manifest_file,
                    library,
                    lanes[i],
                    lane_slice,
                    barcodes[i],
                    locus_function_list,
                    scripts_folder,
                    output_folder,
                    analysis_folder,
                ]
                call_to_taskrunner(output_folder, call_args)

                # Call filter_unmapped_bam
                if gen_read1_plot:
                    output_file = (
                        "{}/logs/filter_unmapped_bam_{}_{}_{}_{}_{}.log".format(
                            output_folder,
                            library,
                            lanes[i],
                            lane_slice,
                            barcodes[i],
                            reference2,
                        )
                    )
                    submission_script = "{}/filter_unmapped_bam.sh".format(
                        scripts_folder
                    )
                    call_args = [
                        "qsub",
                        "-o",
                        output_file,
                        "-l",
                        "h_vmem=10g",
                        "-notify",
                        "-l",
                        "h_rt=10:0:0",
                        "-j",
                        "y",
                        "-P",
                        "macosko_lab",
                        "-l",
                        "os=RedHat7",
                        submission_script,
                        manifest_file,
                        library,
                        lanes[i],
                        lane_slice,
                        barcodes[i],
                        locus_function_list,
                        scripts_folder,
                        output_folder,
                        analysis_folder,
                    ]
                    call_to_taskrunner(output_folder, call_args)

        # Call generate_plots_cmatcher
        output_file = "{}/logs/generate_plots_cmatcher_{}_{}.log".format(
            output_folder, library, reference2
        )
        submission_script = "{}/generate_plots_cmatcher.sh".format(scripts_folder)
        call_args = [
            "qsub",
            "-o",
            output_file,
            "-l",
            "h_vmem=33g",
            "-notify",
            "-l",
            "h_rt=40:0:0",
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
            "{}/{}".format(analysis_folder, reference2),
        ]
        call_to_taskrunner(output_folder, call_args)

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
                + " "
                + reference2
                + " failed at the step of running cmatcher combine. Please check the log file for the issues. "
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
