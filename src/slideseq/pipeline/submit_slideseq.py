#!/usr/bin/python

# This script is to submit a request to the Slide-seq flowcell alignment pipeline

import importlib.resources
import logging
import pathlib
from subprocess import run

import click

import slideseq.scripts
import slideseq.util.google as gutil
import slideseq.util.constants as constants

from slideseq.logging import create_logger
from slideseq.pipeline.metadata import Manifest, validate_flowcell_df
from slideseq.util import get_env_name, qsub_args, get_lanes, get_read_structure
from slideseq.util.constants import MAX_QSUB

log = logging.getLogger(__name__)


def qsub_attempt(qsub_arg_list: list[str], flowcell: str, job_name: str):
    # qsub may sporadically fail due to network issues
    for _ in range(MAX_QSUB):
        proc = run(qsub_arg_list, capture_output=True, text=True)
        if int(proc.returncode) == 0:
            return proc
        else:
            log.warning("qsub failed, retrying")
    else:
        log.error(f"Unable to launch {job_name} job for {flowcell}")
        return None


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

    env_name = get_env_name()
    log.debug(f"Running in Conda env {env_name}")

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
    flowcell_errors = []

    for flowcell in flowcells:
        flowcell_df = worksheet_df.loc[worksheet_df["flowcell"] == flowcell]

        if not len(flowcell_df):
            log.warning(f"Flowcell {flowcell} not found in spreadsheet; please add to sheet.")
            flowcell_errors.append(flowcell)
            continue

        log.debug(f"Found {len(flowcell_df)} in worksheet")

        # subset to metadata columns
        flowcell_df = flowcell_df[constants.METADATA_COLS]
        flowcell_df.columns = [c.lower() for c in constants.METADATA_COLS]

        if not validate_flowcell_df(flowcell, flowcell_df):
            flowcell_errors.append(flowcell)
            continue

        # data locations
        output_dir = (constants.WORKFLOW_DIR / flowcell).absolute().resolve()
        log_dir = output_dir / "logs"
        tmp_dir = output_dir / "tmp"
        flowcell_dir = pathlib.Path(flowcell_df.BCLPath[0])

        runinfo_file = flowcell_dir / "RunInfo.xml"
        if not runinfo_file.exists():
            log.warning(f"{runinfo_file} does not exist, skipping")
            flowcell_errors.append(flowcell)
            continue
        else:
            lanes = get_lanes(run_info_file=runinfo_file)
            read_structure = get_read_structure(run_info_file=runinfo_file)

        manifest_file = output_dir / "manifest.yaml"
        metadata_file = output_dir / "metadata.csv"

        manifest = Manifest(
            flowcell_directory=flowcell_dir,
            output_directory=output_dir,
            metadata_file=metadata_file,
            flowcell=flowcell,
            experiment_date=flowcell_df.date[0].date(),
            is_novaseq=flowcell_df.IlluminaPlatform[0].startswith('NovaSeq'),
            email_addresses=sorted(set(flowcell_df.email)),
        )

        if not dryrun:
            log.debug("Creating output directories")
            output_dir.mkdir(exist_ok=True)
            log_dir.mkdir(exist_ok=True)
            tmp_dir.mkdir(exist_ok=True)

            manifest.to_file(manifest_file)

        # round num_expected_cells to integers
        flowcell_df = flowcell_df.round({"num_expected_cells": 0})

        if dryrun:
            log.info(f"Would write following in {metadata_file}")
            log.info(flowcell_df)
        else:
            # TODO parse metadata first
            flowcell_df.to_csv(metadata_file, sep="\t", header=True, index=False)

        demux_log = log_dir / "demultiplex.log"
        alignment_log = log_dir / "alignment.log"
        processing_log = log_dir / "processing.log"

        # barcode_file, metrics_file, library_params = prepare_demux_files()

        # this script will check the sequencing directory, extract barcodes, and demultiplex to BAM files
        with importlib.resources.path(slideseq.scripts, "demultiplex.sh") as qsub_script:
            # request a high-cpu, high-mem machine for this step
            # TODO: tweak this if it doesn't queue fast enough, maybe put in a config file
            demux_args = qsub_args(
                log_file=demux_log,
                include_env=True,
                PICARD_JAR=constants.PICARD,
                TMP_DIR=tmp_dir,
                BASECALLS_DIR=flowcell_dir / "Data" / "Intensities" / "BaseCalls",
                READ_STRUCTURE=read_structure,
                OUTPUT_DIR=output_dir,
                RUN_BARCODE=flowcell,
                # BARCODE_FILE=barcode_file,
                # METRICS_FILE=metrics_file,
                # LIBRARY_PARAMS=library_params,
            )
            demux_args.extend(
                [
                    "-t",
                    f"{min(lanes)-max(lanes)}",
                    f"{qsub_script.absolute()}",
                    manifest_file,
                ]
            )

        log.debug(f"Demux command:\n\t{' '.join(demux_args)}")

        # this script processes/filters the extracted uBAMs and aligns them
        # depends on previous jobs per lane, so we use -hold_jid_ad for fun
        with importlib.resources.path(slideseq.scripts, "alignment.sh") as qsub_script:
            # request a high-cpu, high-mem machine for this step
            alignment_args = qsub_args(cpu=8, mem=8, hours=8, log_file=alignment_log)
            alignment_args.extend(
                [
                    "-hold_jid_ad",
                    "HOLD_JID",
                    "-t",
                    f"{min(lanes) - max(lanes)}",
                    f"{qsub_script.absolute()}",
                    f"{manifest_file}",
                    f"{metadata_file}",
                ]
            )

        log.debug(f"Alignment command:\n\t{' '.join(alignment_args)}")

        # this script analyzes the alignment output, generates plots, matches to puck, etc
        # need to wait on the whole alignment job to finish, so we use -hold_jid
        with importlib.resources.path(slideseq.scripts, "processing.sh") as qsub_script:
            # this step requires considerably fewer resources
            processing_args = qsub_args(cpu=4, mem=8, hours=1, log_file=processing_log)
            processing_args.extend(
                [
                    "-hold_jid",
                    "HOLD_JID",
                    f"{qsub_script.absolute()}",
                    manifest_file,
                    output_dir,
                ]
            )

        log.debug(f"Processing command:\n\t{' '.join(processing_args)}")

        if not dryrun:
            jid = None
            for arg_list, job_name in (
                (demux_args, "demultiplex"),
                (alignment_args, "alignment"),
                (processing_args, "processing")
            ):
                if jid is not None:
                    arg_list[arg_list.index("HOLD_JID")] = jid

                proc = qsub_attempt(arg_list, flowcell, job_name)
                if proc is None:
                    flowcell_errors.append(flowcell)
                    continue

                # using -terse, proc.stdout is the job id
                jid = proc.stdout.strip()

        submitted.append(flowcell)

    log.info(f"Flowcells {', '.join(submitted)} submitted for processing")
    if flowcell_errors:
        log.info(
            f"Flowcells {', '.join(flowcell_errors)} had errors -- please see warnings above."
        )
