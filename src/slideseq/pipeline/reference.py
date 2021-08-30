# This script is to build genome reference

import csv
import importlib.resources
import logging
import os
import shutil
import sys
from pathlib import Path
from subprocess import run

import click

import slideseq.scripts
import slideseq.util.constants as constants
from slideseq.config import get_config
from slideseq.util import qsub_args
from slideseq.util.logger import create_logger

log = logging.getLogger(__name__)


def check_gtf(reference_gtf: Path, mt_sequence: str):
    warnings = {"header": False, "mt_missing": True}
    fields = ["gene_name", "gene_id"]
    for field in fields:
        warnings[field] = False

    with reference_gtf.open() as fh:
        rdr = csv.reader(fh, delimiter="\t")
        row = next(rdr)
        if row[0][0] != "#":
            warnings["header"] = True
        for row in rdr:
            if row[0][0] == "#":
                continue
            elif row[0].startswith(mt_sequence):
                warnings["mt_missing"] = False

            for field in fields:
                if row[8].find(field) == -1:
                    warnings[field] = True

    if not any(warnings.values()):
        return True
    else:
        if warnings["header"]:
            log.warning("No header detected in GTF file")
        if warnings["mt_missing"]:
            log.warning(f"No entries matched the MT string '{mt_sequence}'")
        for field in fields:
            if warnings[field]:
                log.warning(f"Some rows are missing the field '{field}'")

        return False


@click.command(name="build_ref", no_args_is_help=True)
@click.option("-n", "--genome-name", required=True, help="Name for the reference")
@click.option(
    "-f",
    "--reference-fasta",
    required=True,
    type=click.Path(dir_okay=False, exists=True),
)
@click.option(
    "-g", "--reference-gtf", required=True, type=click.Path(dir_okay=False, exists=True)
)
@click.option(
    "-m", "--mt-sequence", required=True, help="Name in GTF for mitochrondrial sequence"
)
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
    config = get_config()

    reference_fasta = Path(reference_fasta)
    reference_gtf = Path(reference_gtf)
    output_dir = config.reference_dir / genome_name
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

    if check_gtf(reference_gtf, mt_sequence):
        log.info("GTF has all required fields")
    else:
        log.error("Need to fix GTF errors")
        return

    log.info(f"Creating genome reference for {reference_fasta}")

    # this script will create a reference for slideseq
    with importlib.resources.path(
        slideseq.scripts, "build_reference.sh"
    ) as qsub_script:
        mkref_args = qsub_args(
            log_file=output_dir / "build_reference.log",
            CONDA_ENV=env_name,
            PICARD_JAR=config.picard,
            DROPSEQ_DIR=config.dropseq_dir,
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
