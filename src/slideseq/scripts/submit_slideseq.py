#!/usr/bin/python

# This script is to submit a request to the Slide-seq flowcell alignment pipeline

import importlib.resources
import logging
import os
import pathlib
from subprocess import run

import click
import pandas as pd

import slideseq.pipeline.metadata
import slideseq.util.google as gutil
import slideseq.util.constants as constants
from slideseq.logging import create_logger
from slideseq.util import get_env_name


log = logging.getLogger(__name__)


def validate_flowcell_df(flowcell: str, flowcell_df: pd.DataFrame) -> bool:
    warning_logs = []

    if len(set(flowcell_df.library)) != len(flowcell_df.library):
        warning_logs.append(
            f"Flowcell {flowcell} does not have unique library names;"
            f" please fill out naming metadata correctly or add suffixes."
        )

    if any(" " in name for name in flowcell_df.library):
        warning_logs.append(
            f"The 'library' column for {flowcell} contains spaces;"
            " please remove all spaces before running."
        )

    if any(flowcell_df.isna()):
        warning_logs.append(
            f"Flowcell {flowcell} does not have complete submission metadata (orange and blue cols);"
            " please fill out before running."
        )

    if len(set(flowcell_df.BCLPath)) > 1:
        warning_logs.append(
            f"Flowcell {flowcell} has multiple BCLPaths associated with it;"
            " please correct to single path before running."
        )

    if not os.path.isdir(flowcell_df.BCLPath.values[0]):
        warning_logs.append(
            f"Flowcell {flowcell} has incorrect BCLPath associated with it;"
            " please correct path before running."
        )

    if len(set(flowcell_df.IlluminaPlatform)) > 1:
        warning_logs.append(
            f"Flowcell {flowcell} has multiple Illumina platforms associated with it;"
            " please correct to single platform before running."
        )

    # check if platform supported
    if flowcell_df.IlluminaPlatform[0] not in constants.PLATFORMS:
        warning_logs.append(
            f"Flowcell {flowcell} has unsupported platform; please correct to one of"
            f" {' '.join(constants.PLATFORMS)} before running."
        )

    # check if references exist
    if not all(os.path.isfile(build) for build in flowcell_df.reference):
        warning_logs.append(
            f"Reference for {flowcell} does not exist; please correct reference values before running.",
        )

    if warning_logs:
        for msg in warning_logs:
            log.warning(msg)

        return False
    else:
        return True


@click.command(name="submit_slideseq", no_args_is_help=True)
@click.argument("flowcell", nargs=-1, name="flowcells")
@click.option(
    "-s",
    "--spreadsheet",
    default="1kwnKrkbl80LyE9lND0UZZJXipL4yfBbGjkTe6hcwJic",
    help="ID of the Google Sheet to use (default is Macosko Slide-seq Flowcell Alignment)",
)
@click.option("-w", "--worksheet", default="Experiment Log", help="Worksheet to open")
@click.option("-d", "--dryrun", is_flag=True, help="Show the plan but don't execute")
@click.option("--debug", is_flag=True, help="Turn on debug logging")
@click.option("--log_file", type=click.Path(exists=False, writable=True))
def main(flowcells, spreadsheet, worksheet, dryrun=False, debug=False, log_file=None):
    """
    Submit flowcells to the Slide-Seq alignment pipeline.

    See README.md for instructions and requirements: github.com/MacoskoLab/slideseq-tools
    """

    create_logger(debug=debug, log_file=log_file)

    if dryrun:
        log.info("DRY RUN ONLY -- No files will be written and no submissions made.")

    log.debug("Fetching Google credentials")
    google_creds = gutil.get_secrets_manager_credentials()

    log.debug("Setting up Google Drive service")
    drive_service = gutil.get_service(google_creds)

    # get a pandas DataFrame of the worksheet
    log.debug(f"Downloading Google workshet id={spreadsheet} worksheet={worksheet}")
    worksheet_df = gutil.GoogleSheet(drive_service, spreadsheet)[worksheet]
    worksheet_df = worksheet_df.dropna(axis=0, how="all")

    log.debug(f"Retreived worksheet {worksheet} with {len(worksheet_df)} rows")

    log.info(f"Beginning submission for {len(flowcells)} flowcells")
    submitted = []
    skipped = []

    for flowcell in flowcells:
        flowcell_df = worksheet_df.loc[worksheet_df["flowcell"] == flowcell]

        if not len(flowcell_df):
            log.warning(f"Flowcell {flowcell} not found in spreadsheet; please add to sheet.")
            skipped.append(flowcell)
            continue

        log.debug(f"Found {len(flowcell_df)} in worksheet")

        # subset to metadata columns
        flowcell_df = flowcell_df[constants.METADATA_COLS]
        flowcell_df.columns = [c.lower() for c in constants.METADATA_COLS]

        if not validate_flowcell_df(flowcell, flowcell_df):
            skipped.append(flowcell)
            continue

        # data locations
        output_dir = constants.WORKFLOW_DIR / flowcell
        manifest_file = output_dir / "me.manifest"
        metadata_file = output_dir / "me.metadata"

        # this just makes a bunch of directories. skipping / doing it all in here
        # submission_script = run_pipeline.sh
        # it would also launch run_preparation.sh

        # TODO: implement this script as a function or something
        # run_pipeline()

        manifest = slideseq.pipeline.metadata.Manifest(
            flowcell_directory=pathlib.Path(flowcell_df.BCLPath[0]),
            output_directory=output_dir,
            metadata_file=metadata_file,
            flowcell=flowcell,
            experiment_date=flowcell_df.date[0].date(),
            is_novaseq=flowcell_df.IlluminaPlatform[0].startswith('NovaSeq'),
            email_addresses=sorted(set(flowcell_df.email)),
        )

        if not dryrun:
            os.makedirs(output_dir / "logs", exist_ok=True)

            # # write to me.manifest
            # with open(manifest_file, "w") as out:
            #     print(f"flowcell_directory={flowcell_df.BCLPath[0]}", file=out)
            #     print(f"output_folder={output_dir}", file=out)
            #     print(f"metadata_file={metadata_file}", file=out)
            #     print(f"flowcell_barcode={flowcell}", file=out)
            #     print(f"experiment_date={flowcell_df.date[0].date()}", file=out)
            #     print(f"is_NovaSeq={flowcell_df.IlluminaPlatform[0].startswith('NovaSeq')}", file=out)
            #     print(f"email_address={','.join(sorted(set(flowcell_df.email)))}", file=out)

        # round num_expected_cells to integers
        flowcell_df = flowcell_df.round({"num_expected_cells": 0})

        if dryrun:
            log.info(f"Would write following in {metadata_file}")
            log.info(flowcell_df)
        else:
            flowcell_df.to_csv(metadata_file, sep="\t", header=True, index=False)

        output_file = output_dir / "logs" / "demultiplex.log"

        # this script will check the sequencing directory, extract barcodes, and demultiplex to BAM files
        with importlib.resources.path(slideseq.scripts, "demultiplex.sh") as qsub_script:
            call_args = [
                "qsub",
                "-terse",
                "-o",
                f"{output_file}",
                f"{qsub_script.absolute()}",
                manifest_file,
                output_dir,
            ]

        log.debug("Command issued:")
        log.debug(" ".join(call_args))

        if not dryrun:
            proc = run(call_args, capture_output=True, text=True)

            # using -terse, proc.stdout is the job id
            jid = proc.stdout.strip()

        submitted.append(flowcell)

    log.info(f"Flowcells {' '.join(submitted)} submitted for processing")
    if skipped:
        log.info(
            f"Flowcells {'_'.join(skipped)} were skipped -- please see warnings above."
        )
