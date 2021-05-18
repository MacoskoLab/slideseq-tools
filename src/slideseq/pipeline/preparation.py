#!/usr/bin/env python

# This script is to check the Illumina directory, parse input data,
# and call the steps of extracting Illumina barcodes and
# converting barcodes to bam files

import logging

import pandas as pd

import slideseq.util.constants as constants
from slideseq.metadata import Manifest
from slideseq.util import get_lanes

log = logging.getLogger(__name__)


def gen_barcode_file(flowcell_df: pd.DataFrame, manifest: Manifest, lane: int):
    output_file = manifest.output_directory / f"L{lane:03d}" / "barcode_params.txt"

    with output_file.open("w") as out:
        print("barcode_sequence_1\tlibrary_name\tbarcode_name", file=out)

        for _, row in flowcell_df.loc[flowcell_df["lane"] == lane].iterrows():
            # we don't write out barcode_name but the column is required
            print(
                f"{row.sample_barcode}\t{row.library}\t",
                file=out,
            )


def gen_library_params(flowcell_df: pd.DataFrame, manifest: Manifest, lane: int):
    output_file = manifest.output_directory / f"L{lane:03d}" / "library_params.txt"

    with output_file.open("w") as out:
        print("OUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\tBARCODE_1", file=out)

        for _, row in flowcell_df.loc[flowcell_df["lane"] == lane].iterrows():
            # output the uBAM directly to library directory
            library_dir = (
                constants.LIBRARY_DIR / f"{row.date}_{row.library}" / f"L{lane:03d}"
            )
            library_dir.mkdir(exist_ok=True, parents=True)
            output_bam = f"{manifest.flowcell}.L{lane:03d}.{row.library}.{row.sample_barcode}.unmapped.bam"

            print(
                f"{library_dir / output_bam}\t{row.library}\t{row.library}\t{row.sample_barcode}",
                file=out,
            )


def prepare_demux(flowcell_df: pd.DataFrame, manifest: Manifest):
    """create a bunch of directories, and write some input files for picard"""
    output_dir = manifest.output_directory

    runinfo_file = manifest.flowcell_directory / "RunInfo.xml"

    lanes = get_lanes(runinfo_file)

    # Create directories
    log.info(f"Creating directories in {output_dir}")
    for lane in lanes:
        output_lane_dir = output_dir / f"L{lane:03d}"
        output_lane_dir.mkdir(exist_ok=True)
        (output_lane_dir / "barcodes").mkdir(exist_ok=True)

        # Generate barcode_params.txt that is needed by ExtractIlluminaBarcodes
        gen_barcode_file(flowcell_df, manifest, lane)

        # Generate library_params that is needed by IlluminaBasecallsToSam
        gen_library_params(flowcell_df, manifest, lane)


def validate_demux(manifest: Manifest):
    """verify that `prepare_demux` was run previously"""
    output_dir = manifest.output_directory

    if not output_dir.exists():
        log.error(f"{output_dir} does not exist")
        return False

    if not manifest.metadata_file.exists():
        log.error(f"{manifest.metadata_file} does not exist")
        return False

    runinfo_file = manifest.flowcell_directory / "RunInfo.xml"

    lanes = get_lanes(runinfo_file)

    # Create directories
    log.info(f"Checking directories in {output_dir}")
    for lane in lanes:
        for p in (
            output_dir / f"L{lane:03d}",
            output_dir / f"L{lane:03d}" / "barcodes",
            output_dir / f"L{lane:03d}" / "barcode_params.txt",
            output_dir / f"L{lane:03d}" / "library_params.txt",
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
