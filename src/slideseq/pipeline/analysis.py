#!/usr/bin/python

# This script is to generate digital expression and other analysis outputs

import logging
import pathlib

import pandas as pd

import slideseq.util.constants as constants
from slideseq.util.bead_matching import match_barcodes

log = logging.getLogger(__name__)


def run_barcodematching(selected_cells: pathlib.Path, row: pd.Series):
    library_dir = constants.LIBRARY_DIR / f"{row.date}_{row.library}"

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


def gen_downsampling(downsample_dir: pathlib.Path, row: pd.Series):
    # Downsample bam
    downsample_dir.mkdir(exist_ok=True, parents=True)

    # call(["mkdir", "-p", downsample_folder])
    # f1 = f"{alignment_folder}/{library}.AllIllumina.digital_expression_summary.txt"
    # f2 = f"{downsample_folder}/{library}_1.digital_expression_summary.txt"
    # call(["cp", f1, f2])
    # ratio = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    # for i in range(0, 9, 1):
    #     output_file = f"{output_folder}/logs/gen_downsample_dge_{library}_{reference2}_{str(ratio[i])}.log"
    #     submission_script = f"{scripts_folder}/gen_downsample_dge.sh"
    #     call_args = [
    #         "qsub",
    #         "-o",
    #         output_file,
    #         submission_script,
    #         manifest_file,
    #         library,
    #         scripts_folder,
    #         locus_function_list,
    #         str(ratio[i]),
    #         output_folder,
    #         downsample_folder,
    #     ]
    #     call(call_args)

    # Call generate_plot_downsampling
    # output_file = f"{output_folder}/logs/generate_plot_downsampling_{library}_{reference2}.log"
    # submission_script = f"{scripts_folder}/generate_plot_downsampling.sh"
    # call_args = [
    #     "qsub",
    #     "-o",
    #     output_file,
    #     submission_script,
    #     manifest_file,
    #     library,
    #     scripts_folder,
    #     locus_function_list,
    #     output_folder,
    #     barcode_matching_folder,
    # ]
    # call(call_args)
