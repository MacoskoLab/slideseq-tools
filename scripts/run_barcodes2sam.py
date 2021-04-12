#!/usr/bin/python

# This script is to convert Illumina barcodes to unmapped.bam

import logging
import os
import sys
from subprocess import call

from slideseq.logging import create_logger


log = logging.getLogger(__name__)


def main():
    if len(sys.argv) != 5:
        print(
            "Please provide four arguments: manifest file, commandStr, lane ID and slice ID!"
        )
        sys.exit(1)

    manifest_file = sys.argv[1]
    commandStr = sys.argv[2]
    lane = sys.argv[3]
    lane_slice = sys.argv[4]

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

    output_folder = options["output_folder"]
    flowcell_barcode = options["flowcell_barcode"]
    scripts_folder = (
        options["scripts_folder"]
        if "scripts_folder" in options
        else "/broad/macosko/jilong/slideseq_pipeline/scripts"
    )
    email_address = options["email_address"] if "email_address" in options else ""
    log_file = f"{output_folder}/logs/workflow.log"
    create_logger(log_file, logging.INFO)

    folder_running = (
        f"{output_folder}/status/running.barcodes2sam_lane_{lane}_{lane_slice}"
    )
    folder_finished = (
        f"{output_folder}/status/finished.barcodes2sam_lane_{lane}_{lane_slice}"
    )
    folder_failed = (
        f"{output_folder}/status/failed.barcodes2sam_lane_{lane}_{lane_slice}"
    )

    try:
        call(["mkdir", "-p", folder_running])

        # Convert Illumina base calls to sam (unmapped.bam)
        log.info(
            f"{flowcell_barcode} - IlluminaBasecallsToSam for Lane {lane}_{lane_slice}"
        )
        log.info(f"{flowcell_barcode} - Command = {commandStr}")
        os.system(commandStr)
        log.info(f"{flowcell_barcode} - IlluminaBasecallsToSam is done.")

        call(["mv", folder_running, folder_finished])
    except:
        log.exception("EXCEPTION")

        if os.path.isdir(folder_running):
            call(["mv", folder_running, folder_failed])
        else:
            call(["mkdir", "-p", folder_failed])

        if len(email_address) > 1:
            subject = f"Slide-seq workflow failed for {flowcell_barcode}"
            content = (
                f"The Slide-seq workflow for lane {lane} slice {lane_slice} failed at"
                f" the step of converting barcodes to sam."
                f" Please check the log file for the issues."
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
