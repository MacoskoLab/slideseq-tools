#!/usr/bin/env python

import csv
import logging
from pathlib import Path

import pandas as pd

from slideseq.metadata import Manifest
from slideseq.util import get_lanes

log = logging.getLogger(__name__)


def gen_barcode_file(flowcell_df: pd.DataFrame, lane: int, output_file: Path):
    with output_file.open("w") as out:
        wtr = csv.writer(out, delimiter="\t")
        wtr.writerow(("barcode_sequence_1", "library_name", "barcode_name"))

        for _, row in flowcell_df.loc[flowcell_df["lane"] == lane].iterrows():
            # we don't write out barcode_name but the column is required
            wtr.writerow((row.sample_barcode, row.library, ""))


def gen_library_params(
    flowcell_df: pd.DataFrame, lane: int, output_file: Path, library_dir: Path
):
    with output_file.open("w") as out:
        wtr = csv.writer(out, delimiter="\t")
        wtr.writerow(("OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME", "BARCODE_1"))

        for _, row in flowcell_df.loc[flowcell_df["lane"] == lane].iterrows():
            # output the uBAM directly to library directory
            lane_dir = library_dir / f"{row.date}_{row.library}" / f"L{lane:03d}"
            lane_dir.mkdir(exist_ok=True)
            output_bam = f"{row.library}.unmapped.bam"

            wtr.writerow(
                (lane_dir / output_bam, row.library, row.library, row.sample_barcode)
            )


def prepare_demux(
    flowcell_df: pd.DataFrame, lanes: range, output_dir: Path, library_dir: Path
):
    """create a bunch of directories, and write some input files for picard"""
    # Create directories
    log.info(f"Creating directories in {output_dir} and {library_dir}")
    for lane in lanes:
        output_lane_dir = output_dir / f"L{lane:03d}"
        output_lane_dir.mkdir(exist_ok=True)
        (output_lane_dir / "barcodes").mkdir(exist_ok=True)

        # Generate barcode_params.txt that is needed by ExtractIlluminaBarcodes
        gen_barcode_file(flowcell_df, lane, output_lane_dir / "barcode_params.txt")

        # Generate library_params that is needed by IlluminaBasecallsToSam
        gen_library_params(
            flowcell_df, lane, output_lane_dir / "library_params.txt", library_dir
        )


def validate_demux(manifest: Manifest):
    """verify that `prepare_demux` was run previously"""
    if not manifest.workflow_dir.exists():
        log.error(f"{manifest.workflow_dir} does not exist")
        return False

    if not manifest.metadata_file.exists():
        log.error(f"{manifest.metadata_file} does not exist")
        return False

    runinfo_file = manifest.flowcell_dir / "RunInfo.xml"

    lanes = get_lanes(runinfo_file)

    # Create directories
    log.info(f"Checking directories in {manifest.workflow_dir}")
    for lane in lanes:
        output_lane_dir = manifest.workflow_dir / f"L{lane:03d}"

        for p in (
            output_lane_dir,
            output_lane_dir / "barcodes",
            output_lane_dir / "barcode_params.txt",
            output_lane_dir / "library_params.txt",
        ):
            if not p.exists():
                log.error(f"{p} does not exist, demux looks incomplete")
                return False
    else:
        return True


def validate_alignment(manifest: Manifest, n_libraries: int):
    """verify that alignment was run and output is present"""

    for i in range(n_libraries):
        library = manifest.get_library(i)

        for p_list in (
            library.polya_filtering_summaries,
            library.star_logs,
            library.alignment_pickles,
            library.processed_bams,
        ):
            for p in p_list:
                if not p.exists():
                    log.error(f"{p} does not exist, alignment looks incomplete")
                    return False
    else:
        return True
