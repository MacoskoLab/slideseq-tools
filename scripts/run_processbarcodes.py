#!/usr/bin/python

# This script is to extract Illumina barcodes

import logging
import os
import sys
from subprocess import call

from slideseq.logging import create_logger
from slideseq.util import get_read_structure, get_tiles, str2bool

log = logging.getLogger(__name__)


def main():
    if len(sys.argv) != 3:
        print("Please provide two arguments: manifest file and lane ID!")
        sys.exit(1)

    manifest_file = sys.argv[1]
    lane = sys.argv[2]

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

    # Get read structure from RunInfo.xml
    runinfo_file = f"{flowcell_directory}/RunInfo.xml"
    read_structure = get_read_structure(runinfo_file)

    # Get tile information from RunInfo.xml
    slice_id = {}
    slice_first_tile = {}
    slice_tile_limit = {}
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

    try:
        log.info(f"{flowcell_barcode} - Running ExtractIlluminaBarcodes")
        # Extract Illumina barcodes
        commandStr = (
            f"java -Djava.io.tmpdir={tmpdir} -XX:+UseParallelOldGC"
            " -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
            " -Xmx4000m -jar {picard_folder}/picard.jar ExtractIlluminaBarcodes"
            f" TMP_DIR={tmp_dir} VALIDATION_STRINGENCY=SILENT"
            f" BASECALLS_DIR={basecalls_dir} OUTPUT_DIR={output_folder}/{lane}/barcodes"
            f" LANE={lane} READ_STRUCTURE={read_structure}"
            f" BARCODE_FILE={output_folder}/{lane}/barcode_params.txt"
            f" METRICS_FILE={output_folder}/{lane}/{flowcell_barcode}.{lane}.barcode_metrics"
            " COMPRESS_OUTPUTS=true NUM_PROCESSORS=4"
        )
        log.info(f"{flowcell_barcode} - ExtractIlluminaBarcodes for Lane {lane}")
        log.debug(f"Command = {commandStr}")

        os.system(commandStr)

        log.info(
            f"{flowcell_barcode} - ExtractIlluminaBarcodes for Lane {lane} is done."
        )

        # Convert Illumina base calls to sam (unmapped.bam)
        for i in range(len(slice_id[lane])):
            commandStr = (
                f"java -Djava.io.tmpdir={tmpdir}"
                " -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50"
                " -XX:GCHeapFreeLimit=10 -Xmx10192m "
                f" -jar {picard_folder} /picard.jar IlluminaBasecallsToSam"
                f" TMP_DIR={tmpdir} VALIDATION_STRINGENCY=SILENT"
                f" BASECALLS_DIR={basecalls_dir} LANE={lane}"
                f" RUN_BARCODE={flowcell_barcode} NUM_PROCESSORS=4"
                f" READ_STRUCTURE={read_structure}"
                f" LIBRARY_PARAMS={output_folder}/{lane}/{slice_id[lane][i]}/library_params.txt"
                " INCLUDE_NON_PF_READS=false APPLY_EAMSS_FILTER=false"
                " MAX_READS_IN_RAM_PER_TILE=600000 ADAPTERS_TO_CHECK=null"
                " IGNORE_UNEXPECTED_BARCODES=true SEQUENCING_CENTER=BI"
                f" BARCODES_DIR={output_folder}/{lane}/barcodes"
                f" FIRST_TILE={slice_first_tile[lane][i]}"
                f" TILE_LIMIT={slice_tile_limit[lane][i]}"
            )

            output_file = f"{output_folder}/logs/run_barcodes2sam_lane_{lane}_{slice_id[lane][i]}.log"
            submission_script = f"{scripts_folder}/run_barcodes2sam.sh"
            call_args = [
                "qsub",
                "-o",
                output_file,
                submission_script,
                manifest_file,
                commandStr,
                lane,
                slice_id[lane][i],
                scripts_folder,
                output_folder,
                f"{output_folder}/{lane}",
            ]
            call(call_args)
    except:
        log.exception("EXCEPTION!")

        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = (
                f"The Slide-seq workflow for lane {lane} failed at the step of"
                f" processing barcodes. Please check the log file for the issues."
            )
            call_args = [
                "python",
                f"{scripts_folder}/send_email.py",
                email_address,
                subject,
                content,
            ]
            call(call_args)

        sys.exit()


if __name__ == "__main__":
    main()
