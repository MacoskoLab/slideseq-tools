#!/usr/bin/python

# This script is to call the alignment step

import csv
import logging
import os
import sys
import time

from subprocess import call

from new_submit_to_taskrunner import call_to_taskrunner

from slideseq.logging import create_logger
from slideseq.util import get_tiles, str2bool


log = logging.getLogger(__name__)


def main():
    if len(sys.argv) != 2:
        print("Please provide one argument: manifest file!")
        sys.exit(1)

    manifest_file = sys.argv[1]

    # Check if the manifest file exists
    if not os.path.isfile(manifest_file):
        print(f"File {manifest_file} does not exist. Exiting...")
        sys.exit(1)

    # Read manifest file
    options = {}
    with open(manifest_file, "r") as fp:
        for line in fp:
            key, value = line.rstrip().split("=")
            options[key] = value

    flowcell_directory = options["flowcell_directory"]
    output_folder = options["output_folder"]
    metadata_file = options["metadata_file"]
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
    email_address = options["email_address"] if "email_address" in options else ""
    resubmit = str2bool(options["resubmit"]) if "resubmit" in options else False

    runinfo_file = f"{flowcell_directory}/RunInfo.xml"
    log_file = f"{output_folder}/logs/workflow.log"
    create_logger(log_file, logging.INFO)

    # Parse metadata file
    # if resubmit:
    #     write_log(log_file, flowcell_barcode, "Parse metadata file. ")
    #     commandStr = (
    #         "python "
    #         + scripts_folder
    #         + "/parse_metadata.py -i "
    #         + metadata_file
    #         + " -r "
    #         + runinfo_file
    #         + " -o "
    #         + "{}/parsed_metadata.txt".format(output_folder)
    #     )
    #     os.system(commandStr)

    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    references_unique = []
    locus_function_list_unique = []
    resubmit_unique = []
    experiment_date = []
    run_barcodematching = []
    puckcaller_path = []
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
                references_unique.append(row[row0.index("reference")])
                locus_function_list_unique.append(
                    row[row0.index("locus_function_list")]
                )
                resubmit_unique.append(row[row0.index("resubmit")])
                experiment_date.append(row[row0.index("date")])
                run_barcodematching.append(
                    str2bool(row[row0.index("run_barcodematching")])
                )
                puckcaller_path.append(row[row0.index("puckcaller_path")])
            barcodes.append(row[row0.index("sample_barcode")])
            bead_structures.append(row[row0.index("bead_structure")])

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

    folder_waiting = f"{output_folder}/status/waiting.mergebarcodes"
    folder_running = f"{output_folder}/status/running.mergebarcodes"
    folder_finished = f"{output_folder}/status/finished.mergebarcodes"
    folder_failed = f"{output_folder}/status/failed.mergebarcodes"

    call(["mkdir", "-p", folder_waiting])

    # Wait till all of run_processbarcodes and run_barcodes2sam finish
    while 1:
        all_failed_1 = True
        all_failed_2 = True
        for lane in lanes_unique:
            failed_1 = f"{output_folder}/status/failed.processbarcodes_lane_{lane}"
            if not os.path.isdir(failed_1):
                all_failed_1 = False
            for i in range(len(slice_id[lane])):
                failed_2 = f"{output_folder}/status/failed.barcodes2sam_lane_{lane}_{slice_id[lane][i]}"
                if not os.path.isdir(failed_2):
                    all_failed_2 = False
                    break
        if all_failed_1:
            log.error(
                f"{flowcell_barcode} - All run_processbarcodes failed. Exiting..."
            )
            sys.exit(1)
        if all_failed_2:
            log.error(f"{flowcell_barcode} - All run_barcodes2sam failed. Exiting...")
            sys.exit(1)

        f = True
        for lane in lanes_unique:
            for i in range(len(slice_id[lane])):
                fol1 = f"{output_folder}/status/finished.barcodes2sam_lane_{lane}_{slice_id[lane][i]}"
                fol2 = f"{output_folder}/status/failed.barcodes2sam_lane_{lane}_{slice_id[lane][i]}"
                if (not os.path.isdir(fol1)) and (not os.path.isdir(fol2)):
                    f = False
                    break
            if not f:
                break
        if f:
            break
        time.sleep(30)

    if os.path.isdir(folder_waiting):
        call(["mv", folder_waiting, folder_running])
    else:
        call(["mkdir", "-p", folder_running])

    try:
        for j in range(len(libraries_unique)):
            if (not resubmit) or resubmit_unique[j] == "TRUE":
                library = libraries_unique[j]

                # if resubmit:
                #     for i in range(len(lanes)):
                #         if libraries[i] != library:
                #             continue
                #         for lane_slice in slice_id[lanes[i]]:
                #             unmapped_bam = f"{output_folder}/{lanes[i]}/{lane_slice}/{library}/"
                #             if barcodes[i]:
                #                 unmapped_bam += f"{barcodes[i]}/{flowcell_barcode}.{lanes[i]}.{lane_slice}.{library}.{barcodes[i]}.unmapped.bam"
                #             else:
                #                 unmapped_bam += f"{flowcell_barcode}.{lanes[i]}.{lane_slice}.{library}.unmapped.bam"
                #             unmapped_bam2 = f"{library_folder}/{experiment_date[j]}_{library}/{flowcell_barcode}.{lanes[i]}.{lane_slice}.{library}"
                #             if barcodes[i]:
                #                 unmapped_bam2 += "." + barcodes[i]
                #             unmapped_bam2 += ".unmapped.bam"
                #             if os.path.isfile(unmapped_bam2):
                #                 os.system("mv " + unmapped_bam2 + " " + unmapped_bam)
                #     if os.path.isdir(f"{library_folder}/{experiment_date[j]}_{library}"):
                #         os.system(f"rm -r {library_folder}/{experiment_date[j]}_{library}")
                #     os.system(f"rm {output_folder}/logs/*{library}*")
                #     os.system(f"rm -r {output_folder}/status/*{library}*")

                os.system(f"mkdir -p {library_folder}/{experiment_date[j]}_{library}")
                for i in range(len(lanes)):
                    if libraries[i] != library:
                        continue
                    for lane_slice in slice_id[lanes[i]]:
                        unmapped_bam = (
                            f"{output_folder}/{lanes[i]}/{lane_slice}/{library}/"
                        )
                        if barcodes[i]:
                            unmapped_bam += f"{barcodes[i]}/{flowcell_barcode}.{lanes[i]}.{lane_slice}.{library}.{barcodes[i]}.unmapped.bam"
                        else:
                            unmapped_bam += f"{flowcell_barcode}.{lanes[i]}.{lane_slice}.{library}.unmapped.bam"
                        unmapped_bam2 = f"{library_folder}/{experiment_date[j]}_{library}/{flowcell_barcode}.{lanes[i]}.{lane_slice}.{library}"
                        if barcodes[i]:
                            unmapped_bam2 += "." + barcodes[i]
                        unmapped_bam2 += ".unmapped.bam"
                        if os.path.isfile(unmapped_bam):
                            os.system("mv " + unmapped_bam + " " + unmapped_bam2)

                        # Call run_alignment
                        output_file = f"{output_folder}/logs/run_alignment_{library}_{lanes[i]}_{lane_slice}_{barcodes[i]}.log"
                        submission_script = f"{scripts_folder}/run_alignment.sh"
                        call_args = [
                            "qsub",
                            "-o",
                            output_file,
                            submission_script,
                            manifest_file,
                            library,
                            lanes[i],
                            lane_slice,
                            barcodes[i],
                            scripts_folder,
                            output_folder,
                            f"{library_folder}/{experiment_date[j]}_{library}",
                        ]
                        call_to_taskrunner(output_folder, call_args)

                if run_barcodematching[j]:
                    puckcaller_path1 = puckcaller_path[j]
                    file1 = f"{puckcaller_path1}/AnalysisOutputs-selected.mat"
                    file2 = f"{puckcaller_path1}/BeadBarcodes.txt"
                    file3 = f"{puckcaller_path1}/BeadLocations.txt"
                    if puckcaller_path1[-1] != "/":
                        puckcaller_path1 += "/"
                    if (not os.path.isfile(file2)) or (not os.path.isfile(file3)):
                        log.error(f"{file2} and/or {file3} are not found!")
                        if os.path.isfile(file1):
                            output_file = (
                                f"{output_folder}/logs/ExtractBeadBarcode_{library}.log"
                            )
                            submission_script = (
                                f"{scripts_folder}/puckcaller/run_ExtractBeadBarcode.sh"
                            )
                            call_args = [
                                "qsub",
                                "-o",
                                output_file,
                                submission_script,
                                "/broad/software/nonfree/Linux/redhat_7_x86_64/pkgs/matlab_2019a",
                                puckcaller_path1,
                                scripts_folder,
                                output_folder,
                            ]
                            call_to_taskrunner(output_folder, call_args)
                        else:
                            log.error(f"{file1} is not found!")

                if is_NovaSeq or is_NovaSeq_S4:
                    time.sleep(1800)
                else:
                    time.sleep(600)

                # Call run_analysis
                output_file = f"{output_folder}/logs/run_analysis_{library}.log"
                submission_script = f"{scripts_folder}/run_analysis.sh"
                call_args = [
                    "qsub",
                    "-o",
                    output_file,
                    submission_script,
                    manifest_file,
                    library,
                    scripts_folder,
                    output_folder,
                    f"{library_folder}/{experiment_date[j]}_{library}",
                ]
                call_to_taskrunner(output_folder, call_args)

        call(["mv", folder_running, folder_finished])
    except:
        log.exception("EXCEPTION!")

        if os.path.isdir(folder_running):
            call(["mv", folder_running, folder_failed])
        elif os.path.isdir(folder_waiting):
            call(["mv", folder_waiting, folder_failed])
        else:
            call(["mkdir", "-p", folder_failed])

        if len(email_address) > 1:
            subject = f"Slide-seq workflow failed for {flowcell_barcode}"
            content = (
                "The Slide-seq workflow failed at the step of merging barcode matrics."
                " Please check the log file for the issues. "
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
