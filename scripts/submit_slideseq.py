#!/usr/bin/python

# This script is to submit a request to the Slide-seq flowcell alignment pipeline

import argparse
import os
import re
import warnings
from subprocess import call

import gspread
import numpy as np
import pandas as pd
from new_submit_to_taskrunner import call_to_taskrunner
from oauth2client.service_account import ServiceAccountCredentials

scope = [
    "https://spreadsheets.google.com/feeds",
    "https://www.googleapis.com/auth/drive",
]

json_file = "/broad/macosko/jilong/slideseq_pipeline/slideseq-6e435328493d.json"
credentials = ServiceAccountCredentials.from_json_keyfile_name(json_file, scope)

gc = gspread.authorize(credentials)

parser = argparse.ArgumentParser(
    description="Submit flowcells to the Slide-seq flowcell alignment pipeline."
)

parser.add_argument(
    "flowcells",
    metavar="flowcells",
    type=str,
    nargs="+",
    help="Names of flowcells to submit for processing (separated by spaces)",
)
parser.add_argument(
    "--spreadsheet",
    dest="spreadsheet",
    metavar="spreadsheet",
    default="1kwnKrkbl80LyE9lND0UZZJXipL4yfBbGjkTe6hcwJic",
    action="store",
    help="Optional spreadsheet key to open and populate (default is for Macosko Slide-seq Flowcell Alignment)",
)
parser.add_argument(
    "--worksheet_ind",
    dest="wks_ind",
    metavar="worksheet_ind",
    type=int,
    default=0,
    help="Which worksheet to open in spreadsheet (0-indexed) (default 0)",
)
parser.add_argument(
    "--resubmit",
    dest="resubmit",
    default=False,
    action="store_true",
    help="Whether these flowcells are being resubmitted",
)
parser.add_argument(
    "--dryrun",
    dest="dryrun",
    default=False,
    action="store_true",
    help="Show commands without running them",
)

args = parser.parse_args()

if os.path.isfile(args.flowcells[0]):
    with open(args.flowcells[0], "r") as name_file:
        lines = name_file.read().splitlines()

    flowcells = lines
else:
    flowcells = args.flowcells

if args.dryrun:
    print("DRY RUN ONLY -- No files will be written and no submissions made.")

print(f"Beginning submission for {len(flowcells)} flowcells")

# convert non-decimals
non_decimal = re.compile(r"[^\d.]+")

wks = gc.open_by_key(args.spreadsheet).get_worksheet(args.wks_ind)

# IMPORTANT NOTE: row and columns are 1-indexed for range and other utilities
# these column indices are hardcoded for now based on spreadsheet organization
lib_i = 1
gen_downsampling_i = 26

workflow_dir = "/broad/macosko/data/workflows/flowcell"
library_dir = "/broad/macosko/data/libraries"

scripts_folder = "/broad/macosko/jilong/slideseq_pipeline/scripts"
platforms = ["MiniSeq", "NextSeq", "NovaSeq", "NovaSeqS4"]

submitted = []
for flowcell in flowcells:
    if "-" in flowcell and not args.diff_align:
        warnings.warn(
            (
                "Flowcell name follows syntax for different secondary alignment;"
                + " please set --diff_align flag before running."
            ).format(flowcell)
        )
        continue
    flow_cells = wks.findall(flowcell)
    if len(np.unique([cell.col for cell in flow_cells])) > 1:
        warnings.warn(
            (
                "Flowcell {} found in multiple columns;"
                + " not running until resolved."
            ).format(flowcell)
        )
        continue
    elif len(flow_cells) < 1:
        warnings.warn(
            "Flowcell {} not found in spreadsheet; please add to sheet.".format(
                flowcell
            )
        )
        continue
    flow_rows = [cell.row for cell in flow_cells]

    # collect submission metadata
    submit_meta = [
        [cell.value for cell in wks.range(row, lib_i, row, gen_downsampling_i)]
        for row in flow_rows
    ]
    submit_df = pd.DataFrame(
        submit_meta,
        columns=[
            cell.value.lower() for cell in wks.range(1, lib_i, 1, gen_downsampling_i)
        ],
    )
    submit_df.drop(submit_df.columns[[1, 2]], axis=1, inplace=True)

    if len(np.unique(submit_df.library)) != len(submit_df.library):
        warnings.warn(
            (
                "Flowcell {} does not have unique library names; please fill out"
                + " naming metadata correctly or add suffixes."
            ).format(flowcell),
            stacklevel=2,
        )
        continue

    # check for spaces
    if any([" " in name for name in submit_df.library.values]):
        warnings.warn(
            (
                "Flowcell {} has naming metadata containing spaces;"
                + " please remove all spaces before running."
            ).format(flowcell),
            stacklevel=2,
        )
        continue

    if "" in submit_df.values:
        warnings.warn(
            (
                "Flowcell {} does not have complete submission metadata (orange and blue cols);"
                + " please fill out before running."
            ).format(flowcell),
            stacklevel=2,
        )
        continue

    # if bcl paths don't match
    if len(np.unique(submit_df.BCLPath)) > 1:
        warnings.warn(
            (
                "Flowcell {} has multiple BCLPaths associated with it;"
                + " please correct to single path before running."
            ).format(flowcell),
            stacklevel=2,
        )
        continue

    if not os.path.isdir(submit_df.BCLPath.values[0]):
        warnings.warn(
            (
                "Flowcell {} has incorrect BCLPath associated with it;"
                + " please correct path before running."
            ).format(flowcell),
            stacklevel=2,
        )
        continue

    if len(np.unique(submit_df.IlluminaPlatform)) > 1:
        warnings.warn(
            (
                "Flowcell {} has multiple Illumina platforms associated with it;"
                + " please correct to single platform before running."
            ).format(flowcell),
            stacklevel=2,
        )
        continue

    # check if platform supported
    flowcell_platform = submit_df.IlluminaPlatform.values[0]
    if flowcell_platform not in platforms:
        warnings.warn(
            (
                "Flowcell {} has unsupported platform; please correct to one of "
                + " {} before running."
            ).format(flowcell, " ".join(platforms)),
            stacklevel=2,
        )
        continue

    # check if references exist
    build_bool = [os.path.isfile(build) for build in submit_df.reference.values]
    if not all(build_bool):
        warnings.warn(
            (
                "Flowcell {} has reference which does not exist;"
                + " please correct reference values before running."
            ).format(flowcell),
            stacklevel=2,
        )
        continue

    flowcell_name = flowcell
    # without A or B prefix
    flowcell_short = flowcell_name[1:]

    # data locations
    output_dir = "{}/{}".format(workflow_dir, flowcell_short)
    manifest_file = "{}/me.manifest".format(output_dir)
    metadata_file = "{}/me.metadata".format(output_dir)

    submission_script = "{}/run_pipeline.sh".format(scripts_folder)
    if args.resubmit:
        submission_script = "{}/run_mergebarcodes.sh".format(scripts_folder)
        if not os.path.isdir(output_dir):
            warnings.warn(
                (
                    "Folder {} does not exist; please do not use"
                    + " --resubmit flag for this flowcell."
                ).format(output_dir),
                stacklevel=2,
            )
            continue
    else:
        if os.path.isdir(output_dir):
            warnings.warn(
                (
                    "Folder {} exists; please remove the folder or use --resubmit flag."
                ).format(output_dir),
                stacklevel=2,
            )
            continue

    if not args.dryrun:
        if not args.resubmit:
            call(["mkdir", "-p", output_dir])
            call(["mkdir", "-p", "{}/logs".format(output_dir)])
            call(["mkdir", "-p", "{}/tmp_taskrunner".format(output_dir)])
            call(["mkdir", "-p", "{}/tmp_taskrunner/done".format(output_dir)])

        # write to me.manifest
        with open(manifest_file, "w") as f:
            f.write("flowcell_directory=" + submit_df.BCLPath[0] + "\n")
            f.write("scripts_folder=" + scripts_folder + "\n")
            f.write("dropseq_folder=/broad/macosko/bin/dropseq-tools" + "\n")
            f.write(
                "picard_folder=/broad/macosko/bin/dropseq-tools/3rdParty/picard" + "\n"
            )
            f.write(
                "STAR_folder=/broad/macosko/bin/dropseq-tools/3rdParty/STAR-2.5.2a"
                + "\n"
            )
            f.write("output_folder=" + output_dir + "\n")
            f.write("library_folder=" + library_dir + "\n")
            f.write("metadata_file=" + metadata_file + "\n")
            f.write("flowcell_barcode=" + flowcell_short + "\n")
            f.write("experiment_date=" + submit_df.date.values[0] + "\n")
            if submit_df.IlluminaPlatform[0] == "NovaSeq":
                f.write("is_NovaSeq=True" + "\n")
            elif submit_df.IlluminaPlatform[0] == "NovaSeqS4":
                f.write("is_NovaSeq_S4=True" + "\n")
            f.write("email_address=" + ",".join(np.unique(submit_df.email)) + "\n")
            if args.resubmit:
                f.write("resubmit=True" + "\n")

    # round num_expected_cells to integers
    submit_df = submit_df.round({"num_expected_cells": 0})

    if args.dryrun:
        print("Would write following in {}".format(metadata_file))
        print(submit_df)
    else:
        submit_df.to_csv(metadata_file, sep="\t", header=True, index=False)

    output_file = "{}/logs/taskrunner.log".format(output_dir)
    # yes 20 days is huge but needs to stay alive for everything else, also 2g only

    taskrunner_script = "{}/new_taskrunner.sh".format(scripts_folder)
    call_args = [
        "qsub",
        "-o",
        output_file,
        "-l",
        "h_vmem=2g",
        "-notify",
        "-l",
        "h_rt=480:00:0",
        "-j",
        "y",
        "-P",
        "macosko_lab",
        "-l",
        "os=RedHat7",
        taskrunner_script,
        scripts_folder,
        output_dir,
    ]
    call(call_args)

    # Call run pipeline
    if args.resubmit:
        # TODO: this is probably broken
        output_file = f"{output_dir}/logs/run_mergebarcodes.log"
        call_args = [
            "qsub",
            "-o",
            output_file,
            submission_script,
            manifest_file,
            scripts_folder,
            output_dir,
        ]
    else:
        output_file = f"{output_dir}/logs/run_pipeline.log"
        call_args = [
            "qsub",
            "-o",
            output_file,
            submission_script,
            manifest_file,
            scripts_folder,
            output_dir,
        ]

    print("Command issued:")
    print(" ".join(call_args))

    if not args.dryrun:
        call_to_taskrunner(output_dir, call_args)

    submitted.append(flowcell)

    print("Flowcells {} submitted for processing".format(" ".join(submitted)))
    skipped = [flowcell for flowcell in flowcells if flowcell not in submitted]
    if skipped:
        print(
            "\nFlowcells {} were skipped -- please see warnings above.".format(
                "_".join(skipped)
            )
        )
