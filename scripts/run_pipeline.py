#!/usr/bin/python

# This script is to run the Slide-seq flowcell alignment pipeline

import logging
import os
import sys
from subprocess import call

from slideseq.logging import create_logger

log = logging.getLogger(__name__)


def main():
    if len(sys.argv) != 2:
        print("Please provide one argument: manifest file!")
        sys.exit()

    manifest_file = sys.argv[1]

    # Check if the manifest file exists
    if not os.path.isfile(manifest_file):
        print(f"File {manifest_file} does not exist. Exiting...")
        sys.exit()

    # Read manifest file
    options = {}
    with open(manifest_file, "r") as fp:
        for line in fp:
            key, value = line.rstrip().split("=")
            options[key] = value

    if "flowcell_directory" not in options:
        print("flowcell_directory is not specified in the manifest file. Exiting...")
        sys.exit()

    if "output_folder" not in options:
        print("output_folder is not specified in the manifest file. Exiting...")
        sys.exit()

    if "metadata_file" not in options:
        print("metadata_file is not specified in the manifest file. Exiting...")
        sys.exit()

    if "flowcell_barcode" not in options:
        print("flowcell_barcode is not specified in the manifest file. Exiting...")
        sys.exit()

    flowcell_directory = options["flowcell_directory"]
    output_folder = options["output_folder"]
    metadata_file = options["metadata_file"]
    flowcell_barcode = options["flowcell_barcode"]

    dropseq_folder = (
        options["dropseq_folder"]
        if "dropseq_folder" in options
        else "/broad/macosko/bin/dropseq-tools"
    )
    picard_folder = (
        options["picard_folder"]
        if "picard_folder" in options
        else "/broad/macosko/bin/dropseq-tools/3rdParty/picard"
    )
    STAR_folder = (
        options["STAR_folder"]
        if "STAR_folder" in options
        else "/broad/macosko/bin/dropseq-tools/3rdParty/STAR-2.5.2a"
    )
    scripts_folder = (
        options["scripts_folder"]
        if "scripts_folder" in options
        else "/broad/macosko/jilong/slideseq_pipeline/scripts"
    )
    email_address = options["email_address"] if "email_address" in options else ""

    if not os.path.isdir(flowcell_directory):
        print(f"Folder {flowcell_directory} does not exist. Exiting...")
        sys.exit()

    if not os.path.isfile(metadata_file):
        print(f"File {metadata_file} does not exist. Exiting...")
        sys.exit()

    if not os.path.isdir(dropseq_folder):
        print(f"Folder {dropseq_folder} does not exist. Exiting...")
        sys.exit()

    if not os.path.isdir(picard_folder):
        print(f"Folder {picard_folder} does not exist. Exiting...")
        sys.exit()

    if not os.path.isdir(STAR_folder):
        print(f"Folder {STAR_folder} does not exist. Exiting...")
        sys.exit()

    if not os.path.isdir(scripts_folder):
        print(f"Folder {scripts_folder} does not exist. Exiting...")
        sys.exit()

    library_folder = (
        options["library_folder"]
        if "library_folder" in options
        else f"{output_folder}/libraries"
    )

    runinfo_file = f"{flowcell_directory}/RunInfo.xml"
    if not os.path.isfile(runinfo_file):
        print(f"File {runinfo_file} does not exist. Exiting...")
        sys.exit()

    try:
        # Create directories
        if not os.path.isdir(output_folder):
            call(["mkdir", "-p", output_folder])
        if not os.path.isdir(f"{output_folder}/logs"):
            call(["mkdir", "-p", f"{output_folder}/logs"])
        call(["mkdir", "-p", f"{output_folder}/status"])
        if not os.path.isdir(library_folder):
            call(["mkdir", "-p", library_folder])
        if "temp_folder" not in options:
            call(["mkdir", "-p", f"{output_folder}/tmp"])
    except:
        log.exception("EXCEPTION")
        log.error(f"Folder {output_folder} cannot be created. Exiting...")
        sys.exit(1)

    log.info(f"{flowcell_barcode} - starting SlideSeq alignment pipeline")

    # Call run_preparation
    output_file = f"{output_folder}/logs/run_preparation.log"
    submission_script = f"{scripts_folder}/run_preparation.sh"
    call_args = [
        "qsub",
        "-o",
        output_file,
        submission_script,
        manifest_file,
        scripts_folder,
        output_folder,
    ]
    call(call_args)

    if len(email_address) > 1:
        subject = f"Submission received for {flowcell_barcode}"
        content = (
            "Thank you for your interest on the Slide-seq tools! "
            "We received your request. An email will be sent to you "
            "once the workflow finishes. "
        )
        call_args = [
            "python",
            f"{scripts_folder}/send_email.py",
            email_address,
            subject,
            content,
        ]
        call(call_args)


if __name__ == "__main__":
    main()
