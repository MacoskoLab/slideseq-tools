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

    try:
        # Convert Illumina base calls to sam (unmapped.bam)
        log.info(
            f"{flowcell_barcode} - IlluminaBasecallsToSam for Lane {lane}_{lane_slice}"
        )
        log.debug(f"Command = {commandStr}")
        os.system(commandStr)
        log.info(f"{flowcell_barcode} - IlluminaBasecallsToSam is done.")
    except:
        log.exception("EXCEPTION")

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
