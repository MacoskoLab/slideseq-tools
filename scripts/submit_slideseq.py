#!/usr/bin/python

# This script is to submit a request to the Slide-seq flowcell alignment pipeline

import argparse
import logging
import os
import warnings
from subprocess import call

import gspread
import numpy as np
import pandas as pd
from oauth2client.service_account import ServiceAccountCredentials


log = logging.getLogger(__name__)

scope = [
    "https://spreadsheets.google.com/feeds",
    "https://www.googleapis.com/auth/drive",
]

json_file = "/broad/macosko/jilong/slideseq_pipeline/slideseq-6e435328493d.json"
credentials = ServiceAccountCredentials.from_json_keyfile_name(json_file, scope)

gc = gspread.authorize(credentials)


def main():
    parser = argparse.ArgumentParser(
        description="Submit flowcells to the Slide-seq flowcell alignment pipeline."
    )

    parser.add_argument(
        "flowcells",
        nargs="+",
        help="Names of flowcells to submit for processing (separated by spaces)",
    )
    parser.add_argument(
        "--spreadsheet",
        default="1kwnKrkbl80LyE9lND0UZZJXipL4yfBbGjkTe6hcwJic",
        action="store",
        help="Optional spreadsheet key to open and populate (default is for Macosko Slide-seq Flowcell Alignment)",
    )
    parser.add_argument(
        "--worksheet_ind",
        type=int,
        default=0,
        help="Which worksheet to open in spreadsheet (0-indexed) (default 0)",
    )
    parser.add_argument(
        "--dryrun",
        action="store_true",
        help="Show commands without running them",
    )

    args = parser.parse_args()

    flowcells = args.flowcells

    if args.dryrun:
        log.info("DRY RUN ONLY -- No files will be written and no submissions made.")

    log.info(f"Beginning submission for {len(flowcells)} flowcells")

    wks = gc.open_by_key(args.spreadsheet).get_worksheet(args.worksheet_ind)

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
            log.warning(
                (
                    "Flowcell name follows syntax for different secondary alignment;"
                    " please set --diff_align flag before running."
                )
            )
            continue
        flow_cells = wks.findall(flowcell)
        if len(np.unique([cell.col for cell in flow_cells])) > 1:
            log.warning(
                f"Flowcell {flowcell} found in multiple columns; not running until resolved."
            )
            continue
        elif len(flow_cells) < 1:
            log.warning(
                f"Flowcell {flowcell} not found in spreadsheet; please add to sheet."
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
            log.warning(
                (
                    f"Flowcell {flowcell} does not have unique library names;"
                    f" please fill out naming metadata correctly or add suffixes."
                ),
                stacklevel=2,
            )
            continue

        # check for spaces
        if any([" " in name for name in submit_df.library.values]):
            log.warning(
                (
                    f"Flowcell {flowcell} has naming metadata containing spaces;"
                    " please remove all spaces before running."
                ),
                stacklevel=2,
            )
            continue

        if "" in submit_df.values:
            log.warning(
                (
                    f"Flowcell {flowcell} does not have complete submission metadata (orange and blue cols);"
                    " please fill out before running."
                ),
                stacklevel=2,
            )
            continue

        # if bcl paths don't match
        if len(np.unique(submit_df.BCLPath)) > 1:
            log.warning(
                (
                    f"Flowcell {flowcell} has multiple BCLPaths associated with it;"
                    " please correct to single path before running."
                ),
                stacklevel=2,
            )
            continue

        if not os.path.isdir(submit_df.BCLPath.values[0]):
            log.warning(
                (
                    f"Flowcell {flowcell} has incorrect BCLPath associated with it;"
                    " please correct path before running."
                ),
                stacklevel=2,
            )
            continue

        if len(np.unique(submit_df.IlluminaPlatform)) > 1:
            log.warning(
                (
                    f"Flowcell {flowcell} has multiple Illumina platforms associated with it;"
                    " please correct to single platform before running."
                ),
                stacklevel=2,
            )
            continue

        # check if platform supported
        flowcell_platform = submit_df.IlluminaPlatform.values[0]
        if flowcell_platform not in platforms:
            log.warning(
                (
                    f"Flowcell {flowcell} has unsupported platform; please correct to one of"
                    f" {' '.join(platforms)} before running."
                ),
                stacklevel=2,
            )
            continue

        # check if references exist
        build_bool = [os.path.isfile(build) for build in submit_df.reference.values]
        if not all(build_bool):
            warnings.warn(
                f"Reference for {flowcell} does not exist; please correct reference values before running.",
                stacklevel=2,
            )
            continue

        # data locations
        output_dir = f"{workflow_dir}/{flowcell}"
        manifest_file = f"{output_dir}/me.manifest"
        metadata_file = f"{output_dir}/me.metadata"

        submission_script = f"{scripts_folder}/run_pipeline.sh"

        if not args.dryrun:
            call(["mkdir", "-p", output_dir])
            call(["mkdir", "-p", f"{output_dir}/logs"])

            # write to me.manifest
            with open(manifest_file, "w") as out:
                print(f"flowcell_directory={submit_df.BCLPath[0]}", file=out)
                print(f"scripts_folder={scripts_folder}", file=out)
                # print("dropseq_folder=/broad/macosko/bin/dropseq-tools", file=out)
                # print(
                #     "picard_folder=/broad/macosko/bin/dropseq-tools/3rdParty/picard",
                #     file=out
                # )
                # print(
                #     "STAR_folder=/broad/macosko/bin/dropseq-tools/3rdParty/STAR-2.5.2a",
                #     file=out,
                # )
                print(f"output_folder={output_dir}", file=out)
                print(f"library_folder={library_dir}", file=out)
                print(f"metadata_file={metadata_file}", file=out)
                print(f"flowcell_barcode={flowcell}" + flowcell, file=out)
                print(f"experiment_date={submit_df.data.values[0]}", file=out)
                print(f"is_NovaSeq={submit_df.IlluminaPlatform[0].startswith('NovaSeq')}", file=out)
                print(f"is_NovaSeq_S4={submit_df.IlluminaPlatform[0] == 'NovaSeqS4'}", file=out)
                print(f"email_address={','.join(np.unique(submit_df.email))}", file=out)

        # round num_expected_cells to integers
        submit_df = submit_df.round({"num_expected_cells": 0})

        if args.dryrun:
            log.info(f"Would write following in {metadata_file}")
            log.info(submit_df)
        else:
            submit_df.to_csv(metadata_file, sep="\t", header=True, index=False)

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

        log.debug("Command issued:")
        log.debug(" ".join(call_args))

        if not args.dryrun:
            call(call_args)

        submitted.append(flowcell)

        log.info(f"Flowcells {' '.join(submitted)} submitted for processing")
        skipped = [flowcell for flowcell in flowcells if flowcell not in submitted]
        if skipped:
            log.info(
                f"Flowcells {'_'.join(skipped)} were skipped -- please see warnings above."
            )
