# This script is to build genome reference

import importlib.resources
import logging
import os
import pathlib
import shutil
import sys
from subprocess import run

import click

import slideseq.scripts
import slideseq.util.constants as constants
from slideseq.util import qsub_args
from slideseq.util.logger import create_logger

log = logging.getLogger(__name__)


@click.command(name="build_ref", no_args_is_help=True)
@click.option("-n", "--genome-name", help="Name for the reference")
@click.option("-f", "--reference-fasta", type=click.Path(dir_okay=False, exists=True))
@click.option("-g", "--reference-gtf", type=click.Path(dir_okay=False, exists=True))
@click.option(
    "--reference-dir",
    type=click.Path(dir_okay=True, file_okay=False),
    default=constants.REFERENCE_DIR,
    show_default=True,
    help="Location where genome references are stored",
)
@click.option("--mt-sequence", help="Name prefix used in GTF for mitochrondrial genes")
@click.option(
    "-F",
    "--filter-biotypes",
    multiple=True,
    default=constants.FILTERED_BIOTYPES,
    show_default=True,
    help="Gene biotypes to filter out",
)
@click.option("--overwrite", is_flag=True, help="Overwrite any existing reference")
@click.option("--debug", is_flag=True, help="Turn on debug logging")
@click.option("--dryrun", is_flag=True, help="Show the plan but don't execute")
@click.option("--log-file", type=click.Path(exists=False))
def main(
    genome_name: str,
    reference_fasta: str,
    reference_gtf: str,
    reference_dir: str,
    mt_sequence: str,
    filter_biotypes: list[str],
    overwrite: bool = False,
    debug: bool = False,
    dryrun: bool = False,
    log_file: str = None,
):
    create_logger(debug=debug, dryrun=dryrun, log_file=log_file)
    env_name = slideseq.util.get_env_name()
    log.debug(f"Running in Conda env {env_name}")

    reference_fasta = pathlib.Path(reference_fasta)
    reference_gtf = pathlib.Path(reference_gtf)
    output_dir = pathlib.Path(reference_dir) / genome_name
    star_dir = output_dir / "STAR"

    log.info(f"Building reference for genome {genome_name}")

    if output_dir.exists():
        if overwrite:
            log.warning(f"{output_dir} exists, overwriting existing reference")
            if not dryrun:
                if star_dir.exists():
                    log.debug(f"Removing and remaking {star_dir}")
                    shutil.rmtree(star_dir)
                    star_dir.mkdir()
                if (output_dir / f"{genome_name}.dict").exists():
                    log.debug(f"Removing {output_dir}/{genome_name}.dict")
                    os.remove(output_dir / f"{genome_name}.dict")
        else:
            log.error(f"{star_dir} already exists and overwrite=False, aborting")
            sys.exit(1)
    else:
        log.info(f"Creating output directory {star_dir}")
        if not dryrun:
            star_dir.mkdir(parents=True)

    log.info(f"Creating genome reference for {reference_fasta}")

    # this script will create a reference for slideseq
    with importlib.resources.path(
        slideseq.scripts, "build_reference.sh"
    ) as qsub_script:
        mkref_args = qsub_args(
            log_file=output_dir / "build_reference.log",
            CONDA_ENV=env_name,
            PICARD_JAR=constants.PICARD,
            DROPSEQ_DIR=constants.DROPSEQ_DIR,
            GENOME_NAME=genome_name,
            REFERENCE_FASTA=reference_fasta,
            REFERENCE_GTF=reference_gtf,
            OUTPUT_DIR=output_dir,
            MT_SEQUENCE=mt_sequence,
            FILTERED_BIOTYPES=" ".join(f"G={biotype}" for biotype in filter_biotypes),
        )
        mkref_args.append(f"{qsub_script.absolute()}")

        log.debug(f"Build-reference command:\n\t{' '.join(mkref_args)}")
        if dryrun:
            return

        # qsub may sporadically fail due to network issues
        for _ in range(constants.MAX_QSUB):
            proc = run(mkref_args, capture_output=True, text=True)
            if int(proc.returncode) != 0:
                log.warning("qsub failed, retrying")
                log.debug(f"Error: {proc.stderr}")
            else:
                break
        else:
            log.error(f"Unable to launch build_ref job for {reference_fasta}")
