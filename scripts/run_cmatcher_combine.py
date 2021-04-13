#!/usr/bin/python

# This script is to combine outputs from cmatcher
# and call tag_matched_bam

import csv
import logging
import os
import sys
import time
from subprocess import call

import numpy as np

from slideseq.logging import create_logger
from slideseq.util import get_tiles, str2bool

log = logging.getLogger(__name__)


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
        log.critical(f"File {manifest_file} does not exist. Exiting.")
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
        else f"{output_folder}/libraries"
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
    with open(f"{output_folder}/parsed_metadata.txt", "r") as fin:
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

    runinfo_file = f"{flowcell_directory}/RunInfo.xml"
    log_file = f"{output_folder}/logs/workflow.log"
    create_logger(log_file)

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

    analysis_folder = f"{library_folder}/{experiment_date}_{library}"
    alignment_folder = f"{analysis_folder}/{reference2}/alignment/"
    barcode_matching_folder = f"{analysis_folder}/{reference2}/barcode_matching/"
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
        log.error(f"{flowcell_barcode} - {select_cell_file} does not exist!")
        raise FileNotFoundError(f"{select_cell_file} does not exist!")

    try:
        with open(select_cell_file, "r") as fin:
            j = sum(1 for _ in fin)

        k = 10000
        ls = j // k

        log.info("# selected cells: " + str(j))

        # wait for barcode matching and shuffled matching

        log.info(
            f"{flowcell_barcode} - Combine CMatcher outputs for {library} + {reference2}"
        )

        combined_cmatcher_file = (
            f"{barcode_matching_folder}/{library}_barcode_matching.txt"
        )
        combined_cmatcher_header = (
            "IlluminaBarcodes",
            "ProcessedIlluminaBarcodes",
            "BeadBarcodes",
            "Distance",
            "X",
            "Y",
        )
        with open(combined_cmatcher_file, "w") as fout:
            print("\t".join(combined_cmatcher_header), file=fout)
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = f"{barcode_matching_folder}/{library}_barcode_matching_{str(i + 1)}.txt"
                with open(file2, "r") as fin:
                    j = 0
                    for line in fin:
                        j += 1
                        if j > 1:
                            fout.write(line)

        # Combine CMatcher logs
        combined_cmatcher_summary = (
            f"{barcode_matching_folder}/{library}_barcode_matching_summary.txt"
        )
        total = 0
        unique = 0
        multi = 0
        for i in range(ls + 1):
            if i * k >= j:
                break
            file2 = f"{barcode_matching_folder}/{library}_barcode_matching_{str(i + 1)}.txt.log"
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
            print(f"Total # barcodes: {str(total)}", file=fout)
            print(
                f"# unique matched barcodes: {str(unique)}, {str(unique * 100 / total)}%",
                file=fout,
            )
            print(
                f"# multiple matched barcodes: {str(multi)}, {str(multi * 100 / total)}%",
                file=fout,
            )

        for i in range(ls + 1):
            if i * k >= j:
                break
            file1 = (
                f"{barcode_matching_folder}/{library}_barcode_matching_{str(i + 1)}.txt"
            )
            file2 = f"{barcode_matching_folder}/{library}_barcode_matching_{str(i + 1)}.finished"

            name = f"{library}.{min_transcripts_per_cell}_transcripts_mq_{base_quality}_selected_cells"
            file3 = f"{alignment_folder}/{name}_{str(i + 1)}.txt"
            file4 = f"{barcode_matching_folder}/{library}_barcode_matching_{str(i + 1)}.txt.log"
            if os.path.isfile(file1):
                call(["rm", file1])
            if os.path.isfile(file2):
                call(["rm", file2])
            if os.path.isfile(file3):
                call(["rm", file3])
            if os.path.isfile(file4):
                call(["rm", file4])

        combined_cmatcher_file2 = (
            f"{barcode_matching_folder}/{library}_barcode_matching_distance.txt"
        )

        with open(combined_cmatcher_file2, "w") as fout:
            print("\t".join(combined_cmatcher_header), file=fout)
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = f"{barcode_matching_folder}/{library}_barcode_matching_distance_{str(i + 1)}.txt"
                with open(file2) as fin:
                    next(fin)
                    for line in fin:
                        print(line, file=fout)
                call(["rm", file2])

        # UniqueMappedIlluminaBarcodes
        bci = np.loadtxt(
            combined_cmatcher_file, delimiter="\t", dtype="str", skiprows=1, usecols=1
        )
        bci = np.unique(bci)
        unique_bci_file = (
            f"{barcode_matching_folder}/{library}_unique_matched_illumina_barcodes.txt"
        )
        with open(unique_bci_file, "w") as f1:
            for bc in bci:
                f1.write("%s\n" % bc)

        os.system("gzip -c " + unique_bci_file + " > " + unique_bci_file + ".gz")

        log.info(
            f"{flowcell_barcode} - Combine CMatcher outputs for {library} {reference2} completed"
        )

        # Get unique matched bead barcodes and locations
        log.info(
            f"{flowcell_barcode} - Get unique matched bead barcodes and locations for"
            f" {library} {reference2}"
        )
        bead_dict = set()
        matched_bead_barcode_file = (
            f"{barcode_matching_folder}/{library}_matched_bead_barcodes.txt"
        )
        matched_bead_location_file = (
            f"{barcode_matching_folder}/{library}_matched_bead_locations.txt"
        )
        bead_location_forR = f"{barcode_matching_folder}/BeadLocationsForR.csv"
        with open(matched_bead_barcode_file, "w") as fout1:
            with open(matched_bead_location_file, "w") as fout2:
                with open(bead_location_forR, "w") as fout3:
                    print("barcodes,xcoord,ycoord", file=fout3)
                    with open(combined_cmatcher_file, "r") as fin:
                        rdr = csv.reader(fin, delimiter="\t")
                        next(rdr)
                        for row in rdr:
                            bc = row[2]
                            dist = row[3]
                            x = row[4]
                            y = row[5]
                            if bc not in bead_dict:
                                print(bc, file=fout1)
                                print(f"{dist}\t{x}\t{y}", file=fout2)
                                print(f"{bc},{x},{y}", file=fout3)
                                bead_dict.add(bc)

        log.info(
            f"{flowcell_barcode} - Finished unique matched bead barcodes and locations"
            f" for {library} {reference2}"
        )

        combined_cmatcher_file = (
            f"{barcode_matching_folder}/{library}_barcode_matching_shuffled.txt"
        )
        with open(combined_cmatcher_file, "w") as fout:
            print("\t".join(combined_cmatcher_header), file=fout)
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = f"{barcode_matching_folder}/{library}_barcode_matching_shuffled_{str(i + 1)}.txt"
                with open(file2, "r") as fin:
                    next(fin)
                    for line in fin:
                        print(line[:-1], file=fout)

        for i in range(ls + 1):
            if i * k >= j:
                break
            file1 = f"{barcode_matching_folder}/{library}_barcode_matching_shuffled_{str(i + 1)}.txt"
            file2 = f"{barcode_matching_folder}/{library}_barcode_matching_shuffled_{str(i + 1)}.finished"
            name = (
                library
                + "."
                + min_transcripts_per_cell
                + "_transcripts_mq_"
                + base_quality
                + "_selected_cells.shuffled"
            )
            file3 = f"{alignment_folder}/{name}_{str(i + 1)}.txt"
            file4 = f"{barcode_matching_folder}/{library}_barcode_matching_shuffled_{str(i + 1)}.txt.log"
            if os.path.isfile(file1):
                call(["rm", file1])
            if os.path.isfile(file2):
                call(["rm", file2])
            if os.path.isfile(file3):
                call(["rm", file3])
            if os.path.isfile(file4):
                call(["rm", file4])

        combined_cmatcher_file2 = f"{barcode_matching_folder}/{library}_barcode_matching_distance_shuffled.txt"
        with open(combined_cmatcher_file2, "w") as fout:
            print("\t".join(combined_cmatcher_header), file=fout)
            for i in range(ls + 1):
                if i * k >= j:
                    break
                file2 = f"{barcode_matching_folder}/{library}_barcode_matching_distance_shuffled_{str(i + 1)}.txt"
                with open(file2, "r") as fin:
                    next(fin)
                    for line in fin:
                        print(line[:-1], file=fout)

                call(["rm", file2])

        # UniqueMappedIlluminaBarcodes
        bci = np.loadtxt(
            combined_cmatcher_file, delimiter="\t", dtype="str", skiprows=1, usecols=1
        )
        bci = np.unique(bci)
        shuffled_bci_file = (
            f"{barcode_matching_folder}/{library}_unique_shuffled_illumina_barcodes.txt"
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
                output_file = f"{output_folder}/logs/tag_matched_bam_{library}_{lanes[i]}_{lane_slice}_{barcodes[i]}_{reference2}.log"
                submission_script = f"{scripts_folder}/tag_matched_bam.sh"
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
                call(call_args)

                # Call filter_unmapped_bam
                if gen_read1_plot:
                    output_file = f"{output_folder}/logs/filter_unmapped_bam_{library}_{lanes[i]}_{lane_slice}_{barcodes[i]}_{reference2}.log"
                    submission_script = f"{scripts_folder}/filter_unmapped_bam.sh"
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
                    call(call_args)

        # Call generate_plots_cmatcher
        output_file = (
            f"{output_folder}/logs/generate_plots_cmatcher_{library}_{reference2}.log"
        )
        submission_script = f"{scripts_folder}/generate_plots_cmatcher.sh"
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
            f"{analysis_folder}/{reference2}",
        ]
        call(call_args)
    except:
        log.exception("Exception!")

        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = (
                f"The Slide-seq workflow for {library} {reference2} failed at the step "
                "of running cmatcher combine. Please check the log file for the issues."
            )
            call_args = [
                "python",
                f"{scripts_folder}/send_email.py",
                email_address,
                subject,
                content,
            ]
            call(call_args)

        sys.exit(1)


if __name__ == "__main__":
    main()
