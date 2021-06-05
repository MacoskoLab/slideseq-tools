#!/usr/bin/env python

# This script is to check the Illumina directory, parse input data,
# and call the steps of extracting Illumina barcodes and
# converting barcodes to bam files

import logging

import pandas as pd

from slideseq.metadata import Manifest
from slideseq.util import get_lanes

log = logging.getLogger(__name__)


def gen_barcode_file(flowcell_df: pd.DataFrame, manifest: Manifest, lane: int):
    output_file = manifest.workflow_dir / f"L{lane:03d}" / "barcode_params.txt"

    with output_file.open("w") as out:
        print("barcode_sequence_1\tlibrary_name\tbarcode_name", file=out)

        for _, row in flowcell_df.loc[flowcell_df["lane"] == lane].iterrows():
            # we don't write out barcode_name but the column is required
            print(
                f"{row.sample_barcode}\t{row.library}\t",
                file=out,
            )


def gen_library_params(flowcell_df: pd.DataFrame, manifest: Manifest, lane: int):
    output_file = manifest.workflow_dir / f"L{lane:03d}" / "library_params.txt"

    with output_file.open("w") as out:
        print("OUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\tBARCODE_1", file=out)

        for _, row in flowcell_df.loc[flowcell_df["lane"] == lane].iterrows():
            # output the uBAM directly to library directory
            lane_dir = (
                manifest.library_dir / f"{row.date}_{row.library}" / f"L{lane:03d}"
            )
            lane_dir.mkdir(exist_ok=True, parents=True)
            output_bam = f"{row.library}.unmapped.bam"

            print(
                f"{lane_dir / output_bam}\t{row.library}\t{row.library}\t{row.sample_barcode}",
                file=out,
            )


def prepare_demux(flowcell_df: pd.DataFrame, manifest: Manifest):
    """create a bunch of directories, and write some input files for picard"""
    runinfo_file = manifest.flowcell_dir / "RunInfo.xml"

    lanes = get_lanes(runinfo_file)

    # Create directories
    log.info(f"Creating directories in {manifest.workflow_dir}")
    for lane in lanes:
        output_lane_dir = manifest.workflow_dir / f"L{lane:03d}"
        output_lane_dir.mkdir(exist_ok=True)
        (output_lane_dir / "barcodes").mkdir(exist_ok=True)

        # Generate barcode_params.txt that is needed by ExtractIlluminaBarcodes
        gen_barcode_file(flowcell_df, manifest, lane)

        # Generate library_params that is needed by IlluminaBasecallsToSam
        gen_library_params(flowcell_df, manifest, lane)


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
        for p in (
            manifest.workflow_dir / f"L{lane:03d}",
            manifest.workflow_dir / f"L{lane:03d}" / "barcodes",
            manifest.workflow_dir / f"L{lane:03d}" / "barcode_params.txt",
            manifest.workflow_dir / f"L{lane:03d}" / "library_params.txt",
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
