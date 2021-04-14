#!/usr/bin/python

# This script is to convert Illumina barcodes to unmapped.bam

import logging
import os
import sys

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

    flowcell_barcode = options["flowcell_barcode"]

    # Convert Illumina base calls to sam (unmapped.bam)
    log.info(
        f"{flowcell_barcode} - IlluminaBasecallsToSam for Lane {lane}_{lane_slice}"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - IlluminaBasecallsToSam is done.")


if __name__ == "__main__":
    main()
