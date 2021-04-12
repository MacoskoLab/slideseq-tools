#!/usr/bin/python

# This script is to check the Illumina directory, parse input data,
# and call the steps of extracting Illumina barcodes and
# converting barcodes to bam files

import csv
import logging
import os
import sys
from subprocess import call

from new_submit_to_taskrunner import call_to_taskrunner

from slideseq.logging import create_logger
from slideseq.util import get_tiles, get_read_structure, str2bool


log = logging.getLogger(__name__)


def main():
    if len(sys.argv) != 2:
        print("Please provide one argument: manifest file!")
        sys.exit()

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

    tmpdir = (
        options["temp_folder"] if "temp_folder" in options else f"{output_folder}/tmp"
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
    email_address = options["email_address"] if "email_address" in options else ""

    basecalls_dir = f"{flowcell_directory}/Data/Intensities/BaseCalls"
    log_file = f"{output_folder}/logs/workflow.log"
    create_logger(log_file, logging.INFO)

    # Get read structure from RunInfo.xml
    runinfo_file = f"{flowcell_directory}/RunInfo.xml"
    read_structure = get_read_structure(runinfo_file)

    # Parse metadata file
    log.info(f"{flowcell_barcode} - Parsing metadata file.")
    commandStr = (
        f"python {scripts_folder}/parse_metadata.py -i {metadata_file}"
        f" -r {runinfo_file} -o {output_folder}/parsed_metadata.txt"
    )
    os.system(commandStr)

    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    references_unique = []
    locus_function_list_unique = []
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

    folder_running = f"{output_folder}/status/running.run_preparation"
    folder_finished = f"{output_folder}/status/finished.run_preparation"
    folder_failed = f"{output_folder}/status/failed.run_preparation"

    try:
        call(["mkdir", "-p", folder_running])

        # Check if the input Illumina folder is in correct format
        commandStr = (
            f"java -Djava.io.tmpdir={tmpdir} -XX:+UseParallelOldGC"
            f" -XX:ParallelGCThreads=1 -XX:GCTimeLimit=5"
            f" -XX:GCHeapFreeLimit=10 -Xmx8192m "
            f"-jar {picard_folder}/picard.jar CheckIlluminaDirectory TMP_DIR={tmpdir}"
            f" VALIDATION_STRINGENCY=SILENT"
            f" BASECALLS_DIR={basecalls_dir} READ_STRUCTURE={read_structure}"
        )
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += " LINK_LOCS=false"
        for lane in lanes_unique:
            commandStr += " L=" + lane
        log.info(f"{flowcell_barcode} - CheckIlluminaDirectory Command - {commandStr}")
        os.system(commandStr)
        log.info(f"{flowcell_barcode} - CheckIlluminaDirectory done")

        # Create directories
        log.info(f"{flowcell_barcode} - Creating directories.")
        for lane in lanes_unique:
            call(["mkdir", "-p", f"{output_folder}/{lane}"])
            call(["mkdir", "-p", f"{output_folder}/{lane}/barcodes"])
            for lane_slice in slice_id[lane]:
                call(["mkdir", "-p", f"{output_folder}/{lane}/{lane_slice}"])
        for i, lane in enumerate(lanes):
            for lane_slice in slice_id[lane]:
                if not os.path.isdir(
                    f"{output_folder}/{lane}/{lane_slice}/{libraries[i]}"
                ):
                    call(
                        [
                            "mkdir",
                            "-p",
                            f"{output_folder}/{lane}/{lane_slice}/{libraries[i]}",
                        ]
                    )
                if barcodes[i]:
                    call(
                        [
                            "mkdir",
                            "-p",
                            f"{output_folder}/{lane}/{lane_slice}/{libraries[i]}/{barcodes[i]}",
                        ]
                    )

        # Generate barcode_params.txt that is needed by ExtractIlluminaBarcodes
        for lane in lanes_unique:
            log.info(
                f"{flowcell_barcode} - Generating barcode_params.txt for Lane {lane}"
            )
            commandStr = (
                f"python {scripts_folder}/gen_barcode_params.py"
                f" -i {output_folder}/parsed_metadata.txt"
                f" -o {output_folder}/{lane}/barcode_params.txt -l {lane}"
            )
            os.system(commandStr)

        # Generate library_params that is needed by IlluminaBasecallsToSam
        for lane in lanes_unique:
            log.info(
                f"{flowcell_barcode} - Generating library_params.txt for Lane {lane}"
            )
            for lane_slice in slice_id[lane]:
                commandStr = (
                    f"python {scripts_folder}/gen_library_params.py"
                    f" -i {output_folder}/parsed_metadata.txt"
                    f" -o {output_folder}/{lane}/{lane_slice}/library_params.txt"
                    f" -b {output_folder}/{lane}/{lane_slice}/"
                    f" -n {flowcell_barcode}.{lane}.{lane_slice}"
                    f" -l {lane}"
                )
                os.system(commandStr)

        # Call run_processbarcodes
        for lane in lanes_unique:
            output_file = f"{output_folder}/logs/run_processbarcodes_lane_{lane}.log"
            submission_script = f"{scripts_folder}/run_processbarcodes.sh"
            call_args = [
                "qsub",
                "-o",
                output_file,
                submission_script,
                manifest_file,
                lane,
                scripts_folder,
                output_folder,
                f"{output_folder}/{lane}",
            ]
            call_to_taskrunner(output_folder, call_args)

        # Call run_mergebarcodes
        output_file = f"{output_folder}/logs/run_mergebarcodes.log"
        submission_script = f"{scripts_folder}/run_mergebarcodes.sh"
        call_args = [
            "qsub",
            "-o",
            output_file,
            submission_script,
            manifest_file,
            scripts_folder,
            output_folder,
        ]
        call_to_taskrunner(output_folder, call_args)

        call(["mv", folder_running, folder_finished])
    except:
        log.exception("EXCEPTION")

        if os.path.isdir(folder_running):
            call(["mv", folder_running, folder_failed])
        else:
            call(["mkdir", "-p", folder_failed])

        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = (
                "The Slide-seq workflow failed at the step of preparation."
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
