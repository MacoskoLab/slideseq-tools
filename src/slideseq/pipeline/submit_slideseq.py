#!/usr/bin/env python

# This script is to submit a request to the Slide-seq flowcell alignment pipeline

import importlib.resources
import logging
import pathlib
from subprocess import run

import click

import slideseq.scripts
import slideseq.util.constants as constants
import slideseq.util.gutil as gutil
from slideseq.metadata import Manifest, split_sample_lanes, validate_flowcell_df
from slideseq.pipeline.preparation import prepare_demux, validate_demux
from slideseq.util import get_env_name, get_lanes, get_read_structure, qsub_args
from slideseq.util.constants import MAX_QSUB
from slideseq.util.logger import create_logger

log = logging.getLogger(__name__)


def attempt_qsub(qsub_arg_list: list[str], flowcell: str, job_name: str, dryrun: bool):
    log.debug(f"{job_name} command:\n\t{' '.join(qsub_arg_list)}")
    if dryrun:
        return "DRYRUN"

    # qsub may sporadically fail due to network issues
    for _ in range(MAX_QSUB):
        proc = run(qsub_arg_list, capture_output=True, text=True)
        if int(proc.returncode) == 0:
            # using -terse, proc.stdout is [job id].1-N:1' for array jobs
            return proc.stdout.strip().split(".")[0]
        else:
            log.warning("qsub failed, retrying")
            log.debug(f"Error: {proc.stderr}")
    else:
        log.error(f"Unable to launch {job_name} job for {flowcell}")
        return None


@click.command(name="submit_slideseq", no_args_is_help=True)
@click.argument("flowcells", nargs=-1)
@click.option(
    "--spreadsheet",
    default="1kwnKrkbl80LyE9lND0UZZJXipL4yfBbGjkTe6hcwJic",
    help="ID of the Google Sheet (default is Macosko Slide-seq Flowcell Alignment)",
)
@click.option(
    "--worksheet", default="Experiment Log", help="Worksheet to open", show_default=True
)
@click.option("--demux/--no-demux", default=True, help="Whether to run demultiplexing")
@click.option("--align/--no-align", default=True, help="Whether to run alignment")
@click.option("--dryrun", is_flag=True, help="Show the plan but don't execute")
@click.option("--debug", is_flag=True, help="Turn on debug logging")
@click.option("--log-file", type=click.Path(exists=False), help="File to write logs")
def main(
    flowcells: list[str],
    spreadsheet: str,
    worksheet: str,
    demux: bool = True,
    align: bool = True,  # TODO validate this
    dryrun: bool = False,
    debug: bool = False,
    log_file: str = None,
):
    """
    Submit FLOWCELLS to the Slide-Seq alignment pipeline.

    See README.md for instructions and requirements: github.com/MacoskoLab/slideseq-tools
    """

    create_logger(debug=debug, dryrun=dryrun, log_file=log_file)
    env_name = get_env_name()
    log.debug(f"Running in conda env {env_name}")

    log.debug("Fetching Google credentials")
    google_creds = gutil.get_secrets_manager_credentials()

    log.debug("Setting up Google Drive service")
    drive_service = gutil.get_service(google_creds)

    # get a pandas DataFrame of the worksheet
    log.debug(f"Downloading Google worksheet id={spreadsheet} worksheet={worksheet}")
    worksheet_df = gutil.GoogleSheet(drive_service, spreadsheet)[worksheet]
    worksheet_df = worksheet_df.dropna(axis=0, how="all")

    log.debug(f"Retreived worksheet {worksheet} with {len(worksheet_df)} rows")

    log.info(f"Beginning submission for {len(flowcells)} flowcells")
    submitted = set()
    flowcell_errors = set()

    for flowcell in flowcells:
        flowcell_df = worksheet_df.loc[worksheet_df["flowcell"] == flowcell]

        if not len(flowcell_df):
            log.warning(
                f"Flowcell {flowcell} not found in worksheet; please add to sheet."
            )
            flowcell_errors.add(flowcell)
            continue

        log.debug(f"Found {len(flowcell_df)} libraries in worksheet")

        # subset to metadata columns
        flowcell_df = flowcell_df[constants.METADATA_COLS]
        flowcell_df.columns = [c.lower() for c in constants.METADATA_COLS]

        if not validate_flowcell_df(flowcell, flowcell_df):
            flowcell_errors.add(flowcell)
            continue

        # convert columns to desired types
        flowcell_df = flowcell_df.astype(constants.METADATA_TYPES)

        # data locations
        output_dir = constants.WORKFLOW_DIR / flowcell
        flowcell_dir = pathlib.Path(flowcell_df.bclpath.values[0])

        run_info_file = flowcell_dir / "RunInfo.xml"
        lanes = get_lanes(run_info_file=run_info_file)
        read_structure = get_read_structure(run_info_file=run_info_file)

        manifest_file = output_dir / "manifest.yaml"
        metadata_file = output_dir / "metadata.csv"

        manifest = Manifest(
            flowcell_directory=flowcell_dir,
            output_directory=output_dir,
            metadata_file=metadata_file,
            flowcell=flowcell,
            email_addresses=sorted(
                set(e.strip() for v in flowcell_df.email for e in v.split(","))
            ),
        )

        if not dryrun:
            log.debug("Creating output directories")
            output_dir.mkdir(exist_ok=True)

            manifest.log_dir.mkdir(exist_ok=True)
            if list(manifest.log_dir.glob("*.log")):
                log.warning(
                    "Log files already exist for this job, new output will be appended"
                )

            manifest.tmp_dir.mkdir(exist_ok=True)

            log.debug(f"Writing manifest to {manifest_file}")
            manifest.to_file(manifest_file)

        # break rows out per lane for processing
        flowcell_df = split_sample_lanes(flowcell_df, lanes)

        libraries = sorted(set(flowcell_df.library))
        n_libraries = len(libraries)

        if dryrun:
            log.info(f"Would write following in {metadata_file}")
            log.info(flowcell_df)
        elif demux:
            # make various directories
            prepare_demux(flowcell_df, manifest)
            flowcell_df.to_csv(metadata_file, header=True, index=False)
        elif not validate_demux(manifest):
            # appears that demux was not run previously
            flowcell_errors.add(flowcell)
            continue

        # this script will check the sequencing directory, extract barcodes,
        # and demultiplex to BAM files
        with importlib.resources.path(
            slideseq.scripts, "demultiplex.sh"
        ) as qsub_script:
            # request a high-cpu, high-mem machine for this step
            demux_args = qsub_args(
                log_file=manifest.log_dir / "demultiplex.L00$TASK_ID.log",
                PICARD_JAR=constants.PICARD,
                TMP_DIR=manifest.tmp_dir,
                BASECALLS_DIR=flowcell_dir / "Data" / "Intensities" / "BaseCalls",
                READ_STRUCTURE=read_structure,
                OUTPUT_DIR=output_dir,
                RUN_BARCODE=flowcell,
            )
            demux_args.extend(
                [
                    "-t",
                    f"{min(lanes)}-{max(lanes)}",
                    f"{qsub_script.absolute()}",
                ]
            )

            if demux:
                demux_jid = attempt_qsub(demux_args, flowcell, "demultiplex", dryrun)
                if demux_jid is None:
                    flowcell_errors.add(flowcell)
                    continue
            else:
                demux_jid = None
                log.debug("Skipping demux step")

        if flowcell in flowcell_errors:
            continue

        alignment_jids = dict()

        # this script processes/filters the extracted uBAMs and aligns them to the
        # specified reference. this depends on previous jobs per lane, so we use
        # -hold_jid on the lane-specific demux job
        with importlib.resources.path(slideseq.scripts, "alignment.sh") as qsub_script:
            for lane in lanes:
                if demux and demux_jid is None:
                    log.debug(f"Not aligning {lane} because demux was not submitted")
                    continue

                # request a high-cpu, high-mem machine for this step
                alignment_args = qsub_args(
                    log_file=manifest.log_dir / f"alignment.L00{lane}.$TASK_ID.log",
                    debug=debug,
                    CONDA_ENV=env_name,
                    LANE=lane,
                    MANIFEST=manifest_file,
                )

                if demux:
                    alignment_args.extend(["-hold_jid", f"{demux_jid}[{lane}]"])

                alignment_args.extend(
                    [
                        "-t",
                        f"1-{n_libraries}",
                        f"{qsub_script.absolute()}",
                    ]
                )

                if align:
                    alignment_jids[lane] = attempt_qsub(
                        alignment_args, flowcell, "alignment", dryrun
                    )
                    if alignment_jids[lane] is None:
                        flowcell_errors.add(flowcell)
                else:
                    alignment_jids[lane] = None
                    log.info("Skipping alignment step")

        if flowcell in flowcell_errors:
            continue

        # this script analyzes the alignment output, generates plots, matches to puck, etc
        # this is per-library, which means each library needs to wait on the alignment jobs
        # for the relevant lane(s) that contain that library. We're going to wait on all the
        # lanes, but use hold_jid_ad to wait on only that library's alignments
        with importlib.resources.path(slideseq.scripts, "processing.sh") as qsub_script:
            if align and any(alignment_jids[lane] is None for lane in lanes):
                log.debug("Not processing because some alignments were not submitted")
                continue

            # this step requires considerably fewer resources
            processing_args = qsub_args(
                log_file=manifest.log_dir / "processing.$TASK_ID.log",
                debug=debug,
                CONDA_ENV=env_name,
                MANIFEST=manifest_file,
            )

            if align:
                for lane in lanes:
                    processing_args.extend(["-hold_jid_ad", f"{alignment_jids[lane]}"])

            processing_args.extend(
                [
                    "-t",
                    f"1-{n_libraries}",
                    f"{qsub_script.absolute()}",
                ]
            )

            processing_jid = attempt_qsub(
                processing_args, flowcell, "processing", dryrun
            )
            if processing_jid is None:
                flowcell_errors.add(flowcell)
            else:
                submitted.add(flowcell)

    if submitted and not dryrun:
        log.info(f"Flowcells {', '.join(submitted)} submitted for processing")
    else:
        log.info("No flowcells submitted for processing")

    if flowcell_errors:
        log.info(
            f"Flowcells {', '.join(flowcell_errors)} had errors -- see warnings above."
        )
