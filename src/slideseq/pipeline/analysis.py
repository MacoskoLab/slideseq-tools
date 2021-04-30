#!/usr/bin/python

# This script is to generate digital expression and other analysis outputs

import logging
import pathlib

import pandas as pd

from slideseq.util import dropseq_cmd, picard_cmd, start_popen
from slideseq.util.bead_matching import match_barcodes

log = logging.getLogger(__name__)


def run_barcodematching(
    selected_cells: pathlib.Path, row: pd.Series, library_dir: pathlib.Path
):
    reference = pathlib.Path(row.reference)
    puckcaller_dir = pathlib.Path(row.puckcaller_path)

    reference2 = f"{reference.stem}.{row.locus_function_list}"

    barcode_matching_folder = library_dir / reference2 / "barcode_matching"
    barcode_matching_folder.mkdir(exist_ok=True, parents=True)

    barcode_mapping = match_barcodes(
        selected_cells,
        puckcaller_dir / "BeadBarcodes.txt",
        puckcaller_dir / "BeadLocations.txt",
    )

    with (barcode_matching_folder / f"{row.library}_barcode_matching.txt").open(
        "w"
    ) as out:
        print(
            "\n".join(
                f"{bead_bc}\t{seq_bc}" for bead_bc, seq_bc in barcode_mapping.items()
            ),
            file=out,
        )

    # TODO: standalone script to analyze FDR by shuffling barcodes.
    # Doesn't seem necessary for the main pipeline?


def downsample_dge(
    bam_file: pathlib.Path,
    downsample_dir: pathlib.Path,
    row: pd.Series,
    ratio: float,
    tmp_dir: pathlib.Path,
):
    digital_expression_summary = (
        downsample_dir / f"{row.library}_{ratio}.digital_expression_summary.txt"
    )

    procs = []

    # Downsample reads
    cmd = picard_cmd("DownsampleSam", tmp_dir)
    cmd.extend(["--INPUT", bam_file, "--OUTPUT", "/dev/stdout", "-P", ratio])
    procs.append(start_popen(cmd, "DownsampleSam", row.library))

    # output to /dev/null because we don't want to keep the DGE matrix
    cmd = dropseq_cmd("DigitalExpression", "/dev/stdin", "/dev/null")
    cmd.extend(
        [
            f"SUMMARY={digital_expression_summary}",
            "EDIT_DISTANCE=1",
            f"MIN_NUM_TRANSCRIPTS_PER_CELL={row.min_transcripts_per_cell}",
            f"READ_MQ={row.base_quality}",
            "MIN_BC_READ_THRESHOLD=0",
            "OUTPUT_HEADER=true",
            f"UEI={row.library}",
        ]
    )
    if row.locus_function_list == "intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=null", "LOCUS_FUNCTION_LIST=INTRONIC"])
    elif row.locus_function_list == "exonic+intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=INTRONIC"])

    procs.append(start_popen(cmd, "DigitalExpression", row.library, procs[-1]))

    # close intermediate streams
    for p in procs[:-1]:
        p.stdout.close()

    # wait for final process to finish
    procs[-1].communicate()
    log.debug(f"Finished with downsampling at ratio {ratio}")
