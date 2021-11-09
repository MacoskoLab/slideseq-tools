#!/usr/bin/env python

# This script is to submit a request to the Slide-seq flowcell alignment pipeline

import importlib.resources
import logging
from pathlib import Path
from subprocess import run

import click

import slideseq.scripts
import slideseq.util.constants as constants
import slideseq.util.gutil as gutil
from slideseq.config import get_config
from slideseq.metadata import Manifest, split_sample_lanes, validate_run_df
from slideseq.pipeline.preparation import (
    prepare_demux,
    validate_alignment,
    validate_demux,
)
from slideseq.util import get_env_name, qsub_args
from slideseq.util.logger import create_logger
from slideseq.util.run_info import get_run_info

log = logging.getLogger(__name__)


def attempt_qsub(qsub_arg_list: list[str], run_name: str, job_name: str, dryrun: bool):
    log.debug(f"{job_name} command:\n\t{' '.join(qsub_arg_list)}")
    if dryrun:
        return "DRYRUN"

    # qsub may sporadically fail due to network issues
    for _ in range(constants.MAX_QSUB):
        proc = run(qsub_arg_list, capture_output=True, text=True)
        if int(proc.returncode) == 0:
            # using -terse, proc.stdout is [job id].1-N:1' for array jobs
            return proc.stdout.strip().split(".")[0]
        else:
            log.warning("qsub failed, retrying")
            log.debug(f"Error: {proc.stderr}")
    else:
        log.error(f"Unable to launch {job_name} job for {run_name}")
        return None


@click.command(name="submit_slideseq", no_args_is_help=True)
@click.argument("runs", nargs=-1, metavar="RUN [RUN ...]")
@click.option("--demux/--no-demux", default=True, help="Whether to run demultiplexing")
@click.option("--align/--no-align", default=True, help="Whether to run alignment")
@click.option("--dryrun", is_flag=True, help="Show the plan but don't execute")
@click.option("--debug", is_flag=True, help="Turn on debug logging")
@click.option("--log-file", type=click.Path(exists=False), help="File to write logs")
def main(
    runs: list[str],
    demux: bool = True,
    align: bool = True,
    dryrun: bool = False,
    debug: bool = False,
    log_file: str = None,
):
    """
    Submit each RUN to the Slide-seq alignment pipeline.

    See README.md for instructions and requirements: github.com/MacoskoLab/slideseq-tools
    """

    create_logger(debug=debug, dryrun=dryrun, log_file=log_file)
    env_name = get_env_name()
    log.debug(f"Running in conda env {env_name}")
    config = get_config()

    # you shouldn't demux without aligning, that's weird
    if demux and not align:
        log.debug("Assuming --no-align should also set --no-demux")
        demux = False

    log.debug("Fetching Google credentials")
    google_creds = gutil.get_secrets_manager_credentials(config.gsecret_name)

    log.debug("Setting up Google Drive service")
    drive_service = gutil.get_service(google_creds)

    # get a pandas DataFrame of the worksheet
    log.debug(
        f"Downloading Google Sheet, id={config.gsheet_id} worksheet={config.worksheet}"
    )
    worksheet_df = gutil.GoogleSheet(drive_service, config.gsheet_id)[config.worksheet]
    worksheet_df = worksheet_df.dropna(axis=0, how="all")

    log.debug(f"Retreived worksheet {config.worksheet} with {len(worksheet_df)} rows")

    log.info(f"Beginning submission for {len(runs)} runs")
    submitted = set()
    run_errors = set()

    for run_name in runs:
        run_df = worksheet_df.loc[worksheet_df["run_name"] == run_name]

        if not len(run_df):
            log.warning(f"{run_name} not found in worksheet; please add to sheet.")
            run_errors.add(run_name)
            continue

        log.debug(f"Found {len(run_df)} libraries in worksheet")

        # subset to metadata columns
        run_df = run_df[constants.METADATA_COLS]
        run_df.columns = [c.lower() for c in constants.METADATA_COLS]

        if not validate_run_df(run_name, run_df):
            run_errors.add(run_name)
            continue

        # convert columns to desired types
        run_df = run_df.astype(constants.METADATA_TYPES)

        # data locations
        output_dir = config.workflow_dir / run_name
        flowcell_dirs = sorted(Path(fd) for fd in set(run_df.bclpath))

        manifest_file = output_dir / "manifest.yaml"
        metadata_file = output_dir / "metadata.csv"

        run_info_list = []

        for flowcell_dir in flowcell_dirs:
            run_info_list.append(get_run_info(flowcell_dir))

        # break rows out per flowcell/lane for processing
        run_df = split_sample_lanes(run_df, run_info_list)

        manifest = Manifest(
            run_name=run_name,
            flowcell_dirs=flowcell_dirs,
            workflow_dir=output_dir,
            library_dir=config.library_dir,
            metadata_file=metadata_file,
            metadata=split_sample_lanes(run_df, run_info_list),
            email_addresses=sorted(
                set(e.strip() for v in run_df.email for e in v.split(","))
            ),
        )

        n_libraries = len(list(manifest.libraries))

        if not dryrun:
            log.debug("Creating output directories")
            output_dir.mkdir(exist_ok=True)

            manifest.log_dir.mkdir(exist_ok=True)
            if list(manifest.log_dir.glob("*.log")):
                log.warning(
                    "Log files already exist for this job, new output will be appended"
                )

            manifest.library_dir.mkdir(exist_ok=True)
            manifest.tmp_dir.mkdir(exist_ok=True)

            log.debug(f"Writing manifest to {manifest_file}")
            if metadata_file.exists():
                log.info(f"Overwriting metadata file {metadata_file}")
            manifest.to_file(manifest_file)

            if demux:
                # make various directories
                prepare_demux(run_info_list, manifest)
            elif not validate_demux(manifest):
                # appears that demux was not run previously
                run_errors.add(run_name)
                continue
            elif not (align or validate_alignment(manifest, n_libraries)):
                # appears that alignment was not run and wasn't requested
                run_errors.add(run_name)
                continue

        demux_jids = dict()

        # this script will check the sequencing directory, extract barcodes,
        # and demultiplex to BAM files
        with importlib.resources.path(
            slideseq.scripts, "demultiplex.sh"
        ) as qsub_script:
            for run_info in run_info_list:
                demux_args = qsub_args(
                    log_file=manifest.log_dir / run_info.demux_log,
                    email=",".join(manifest.email_addresses),
                    PICARD_JAR=config.picard,
                    TMP_DIR=manifest.tmp_dir,
                    FLOWCELL=run_info.flowcell,
                    BASECALLS_DIR=run_info.basecall_dir,
                    READ_STRUCTURE=run_info.read_structure,
                    OUTPUT_DIR=output_dir,
                )
                demux_args.extend(
                    [
                        "-t",
                        f"{min(run_info.lanes)}-{max(run_info.lanes)}",
                        f"{qsub_script.absolute()}",
                    ]
                )

                if demux:
                    demux_jids[run_info.flowcell] = attempt_qsub(
                        demux_args, run_name, "demultiplex", dryrun
                    )
                    if demux_jids[run_info.flowcell] is None:
                        run_errors.add(run_name)
                        continue
                else:
                    demux_jids[run_info.flowcell] = None
                    log.debug("Skipping demux step")

        if run_name in run_errors:
            continue

        alignment_jids = dict()

        # this script processes/filters the extracted uBAMs and aligns them to the
        # specified reference. this depends on previous jobs per lane, so we use
        # -hold_jid on the lane-specific demux job
        with importlib.resources.path(slideseq.scripts, "alignment.sh") as qsub_script:
            for run_info in run_info_list:
                for lane in run_info.lanes:
                    if demux and demux_jids[run_info.flowcell] is None:
                        log.debug(
                            f"Not aligning {lane} because demux was not submitted"
                        )
                        continue

                    alignment_args = qsub_args(
                        log_file=manifest.log_dir / run_info.alignment_log(lane),
                        email=",".join(manifest.email_addresses),
                        debug=debug,
                        CONDA_ENV=env_name,
                        FLOWCELL=run_info.flowcell,
                        LANE=lane,
                        MANIFEST=manifest_file,
                    )

                    if demux:
                        alignment_args.extend(
                            ["-hold_jid", f"{demux_jids[run_info.flowcell]}[{lane}]"]
                        )

                    alignment_args.extend(
                        [
                            "-t",
                            f"1-{n_libraries}",
                            f"{qsub_script.absolute()}",
                        ]
                    )

                    if align:
                        alignment_jids[run_info.flowcell, lane] = attempt_qsub(
                            alignment_args, run_name, "alignment", dryrun
                        )
                        if alignment_jids[run_info.flowcell, lane] is None:
                            run_errors.add(run_name)
                    else:
                        alignment_jids[run_info.flowcell, lane] = None
                        log.info("Skipping alignment step")

        if run_name in run_errors:
            continue

        # this script analyzes the alignment output, generates plots, matches to puck, etc
        # this is per-library, which means each library needs to wait on the alignment jobs
        # for the relevant lane(s) that contain that library. We're going to wait on all the
        # lanes, but use hold_jid_ad to wait on only that library's alignments
        with importlib.resources.path(slideseq.scripts, "processing.sh") as qsub_script:
            if align and any(
                alignment_jids[run_info.flowcell, lane] is None
                for run_info in run_info_list
                for lane in run_info.lanes
            ):
                log.debug("Not processing because some alignments were not submitted")
                continue

            processing_args = qsub_args(
                log_file=manifest.log_dir / "processing.$TASK_ID.log",
                email=",".join(manifest.email_addresses),
                debug=debug,
                CONDA_ENV=env_name,
                MANIFEST=manifest_file,
            )

            if align:
                for run_info in run_info_list:
                    for lane in run_info.lanes:
                        processing_args.extend(
                            [
                                "-hold_jid_ad",
                                f"{alignment_jids[run_info.flowcell, lane]}",
                            ]
                        )

            processing_args.extend(
                [
                    "-t",
                    f"1-{n_libraries}",
                    f"{qsub_script.absolute()}",
                ]
            )

            processing_jid = attempt_qsub(
                processing_args, run_name, "processing", dryrun
            )
            if processing_jid is None:
                run_errors.add(run_name)
            else:
                submitted.add(run_name)

    if submitted and not dryrun:
        log.info(f"Flowcells {', '.join(submitted)} submitted for processing")
    else:
        log.info("No flowcells submitted for processing")

    if run_errors:
        log.info(f"Flowcells {', '.join(run_errors)} had errors -- see warnings above.")
