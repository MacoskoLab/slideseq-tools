#!/usr/bin/env python

import csv
import logging
from pathlib import Path

from slideseq.metadata import Manifest
from slideseq.util.run_info import RunInfo, get_run_info

log = logging.getLogger(__name__)


def gen_barcode_file(manifest: Manifest, flowcell: str, lane: int, output_file: Path):
    with output_file.open("w") as out:
        wtr = csv.writer(out, delimiter="\t")
        wtr.writerow(("barcode_sequence_1", "library_name", "barcode_name"))

        for library in manifest.libraries:
            if (flowcell, lane) in library.samples:
                for barcode in library.samples[flowcell, lane]:
                    # we don't write out barcode_name but the column is required
                    wtr.writerow((barcode, library.name, ""))


def gen_library_params(manifest: Manifest, flowcell: str, lane: int, output_file: Path):
    with output_file.open("w") as out:
        wtr = csv.writer(out, delimiter="\t")
        wtr.writerow(("OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME", "BARCODE_1"))

        for sample in manifest.samples:
            if sample.flowcell == flowcell and sample.lane == lane:
                # output the uBAM directly to library directory
                sample.lane_dir.mkdir(exist_ok=True, parents=True)

                for barcode, output_ubam in zip(sample.barcodes, sample.barcode_ubams):
                    wtr.writerow((output_ubam, sample.name, sample.name, barcode))


def prepare_demux(run_info_list: list[RunInfo], manifest: Manifest):
    """create a bunch of directories, and write some input files for picard"""
    # Create directories
    log.info(
        f"Creating directories in {manifest.workflow_dir} and {manifest.library_dir}"
    )

    for run_info in run_info_list:
        for lane in run_info.lanes:
            output_lane_dir = manifest.workflow_dir / run_info.flowcell / f"L{lane:03d}"

            output_lane_dir.mkdir(exist_ok=True, parents=True)
            (output_lane_dir / "barcodes").mkdir(exist_ok=True)

            # Generate barcode_params.txt that is needed by ExtractIlluminaBarcodes
            gen_barcode_file(
                manifest,
                run_info.flowcell,
                lane,
                output_lane_dir / "barcode_params.txt",
            )

            # Generate library_params that is needed by IlluminaBasecallsToSam
            gen_library_params(
                manifest,
                run_info.flowcell,
                lane,
                output_lane_dir / "library_params.txt",
            )


def validate_demux(manifest: Manifest):
    """verify that `prepare_demux` was run previously"""
    if not manifest.workflow_dir.exists():
        log.error(f"{manifest.workflow_dir} does not exist")
        return False

    for flowcell_dir in manifest.flowcell_dirs:
        run_info = get_run_info(flowcell_dir)

        # Create directories
        log.info(f"Checking directories in {manifest.workflow_dir / run_info.flowcell}")
        for lane in run_info.lanes:
            output_lane_dir = manifest.workflow_dir / run_info.flowcell / f"L{lane:03d}"

            for p in (
                output_lane_dir,
                output_lane_dir / "barcodes",
                output_lane_dir / "barcode_params.txt",
                output_lane_dir / "library_params.txt",
            ):
                if not p.exists():
                    log.error(f"{p} does not exist, demux looks incomplete")
                    return False

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
