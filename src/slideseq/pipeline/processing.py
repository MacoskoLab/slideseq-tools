#!/usr/bin/python

import gzip
import logging
import os
from pathlib import Path

import click
import numpy as np
import pandas as pd

import slideseq.alignment_quality as alignment_quality
import slideseq.util.constants as constants
from slideseq.bead_matching import match_barcodes
from slideseq.metadata import Manifest
from slideseq.pipeline.downsampling import downsample_dge
from slideseq.plot.plot_library_metrics import make_library_plots
from slideseq.retag_bam import write_retagged_bam
from slideseq.util import dropseq_cmd, picard_cmd, run_command
from slideseq.util.logger import create_logger

log = logging.getLogger(__name__)


def validate_library_df(library: str, library_df: pd.DataFrame):
    """Verify that all of the columns in the dataframe are constant, except
    for the lane which was expanded out earlier"""

    for col in constants.METADATA_COLS:
        if col == "lane":
            continue

        if len(set(library_df[col])) != 1:
            raise ValueError(f"Library {library} has multiple values in column {col}")


def calc_alignment_metrics(
    input_bam: Path,
    reference: Path,
    row: pd.Series,
    manifest: Manifest,
    matched_barcodes: Path = None,
    cell_tag: str = "XC",
):
    reference_dir = reference.parent

    # assumes that reference always ends in .fasta, not .fasta.gz
    ref_flat = reference_dir / f"{reference.stem}.refFlat"
    ribosomal_intervals = reference_dir / f"{reference.stem}.rRNA.intervals"

    reads_per_cell = input_bam.with_suffix(
        f".numReads_perCell_{cell_tag}_mq_{row.base_quality}.txt.gz"
    )
    reads_per_umi = input_bam.with_suffix(
        f".numReads_perUMI_XM_mq_{row.base_quality}.txt.gz"
    )
    frac_intronic_exonic = input_bam.with_suffix(".fracIntronicExonic.txt")
    xc_barcode_distribution = input_bam.with_suffix(".barcode_distribution_XC.txt")
    xm_barcode_distribution = input_bam.with_suffix(".barcode_distribution_XM.txt")
    read_quality_metrics = input_bam.with_suffix(".ReadQualityMetrics.txt")
    scnra_metrics = input_bam.with_suffix(".fracIntronicExonicPerCell.txt.gz")

    selected_cells = input_bam.with_suffix(
        f".{row.min_transcripts_per_cell}_transcripts_mq_{row.base_quality}_selected_cells.txt.gz"
    )
    digital_expression = input_bam.with_suffix(".digital_expression.txt.gz")
    digital_expression_summary = input_bam.with_suffix(
        ".digital_expression_summary.txt"
    )

    # Bam tag histogram (cells)
    cmd = dropseq_cmd("BamTagHistogram", input_bam, reads_per_cell)
    cmd.extend(
        [
            f"TAG={cell_tag}",
            "FILTER_PCR_DUPLICATES=false",
            f"READ_MQ={row.base_quality}",
            f"TMP_DIR={manifest.tmp_dir}",
        ]
    )

    run_command(cmd, "BamTagHistogram (cells)", row.library)

    # Bam tag histogram (UMIs)
    cmd = dropseq_cmd("BamTagHistogram", input_bam, reads_per_umi)
    cmd.extend(
        [
            "TAG=XM",
            "FILTER_PCR_DUPLICATES=false",
            f"READ_MQ={row.base_quality}",
            f"TMP_DIR={manifest.tmp_dir}",
        ]
    )

    run_command(cmd, "BamTagHistogram (UMIs)", row.library)

    # Collect RnaSeq metrics
    cmd = picard_cmd("CollectRnaSeqMetrics", manifest.tmp_dir)
    cmd.extend(
        [
            "--INPUT",
            input_bam,
            "--REF_FLAT",
            ref_flat,
            "--OUTPUT",
            frac_intronic_exonic,
            "--STRAND_SPECIFICITY",
            "NONE",
            "--RIBOSOMAL_INTERVALS",
            ribosomal_intervals,
        ]
    )

    run_command(cmd, "CollectRnaSeqMetrics", row.library)

    # Base distribution at read position for raw cellular barcode
    cmd = dropseq_cmd(
        "BaseDistributionAtReadPosition", input_bam, xc_barcode_distribution
    )
    cmd.extend(["TAG=XC"])
    run_command(cmd, "BaseDistributionAtReadPosition (Cellular)", row.library)

    # Base distribution at read position for molecular barcode
    cmd = dropseq_cmd(
        "BaseDistributionAtReadPosition", input_bam, xm_barcode_distribution
    )
    cmd.extend(["TAG=XM"])
    run_command(cmd, "BaseDistributionAtReadPosition (Molecular)", row.library)

    # Gather read quality metrics
    cmd = dropseq_cmd("GatherReadQualityMetrics", input_bam, read_quality_metrics)
    run_command(cmd, "GatherReadQualityMetrics", row.library)

    if matched_barcodes is not None:
        cmd = dropseq_cmd("SingleCellRnaSeqMetricsCollector", input_bam, scnra_metrics)
        cmd.extend(
            [
                f"CELL_BARCODE_TAG={cell_tag}",
                f"READ_MQ={row.base_quality}",
                f"CELL_BC_FILE={matched_barcodes}",
                f"RIBOSOMAL_INTERVALS={ribosomal_intervals}",
            ]
        )
        run_command(cmd, "SingleCellRnaSeqMetricsCollector", row.library)

    # Select cells by num transcripts
    cmd = dropseq_cmd("SelectCellsByNumTranscripts", input_bam, selected_cells)
    cmd.extend(
        [
            f"MIN_TRANSCRIPTS_PER_CELL={row.min_transcripts_per_cell}",
            f"READ_MQ={row.base_quality}",
        ]
    )
    if row.locus_function_list == "intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=null", "LOCUS_FUNCTION_LIST=INTRONIC"])
    elif row.locus_function_list == "exonic+intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=INTRONIC"])
    run_command(cmd, "SelectCellsByNumTranscripts", row.library)

    # Generate digital expression files for all Illumina barcodes
    cmd = dropseq_cmd("DigitalExpression", input_bam, digital_expression, compression=6)
    cmd.extend(
        [
            f"SUMMARY={digital_expression_summary}",
            "EDIT_DISTANCE=1",
            f"READ_MQ={row.base_quality}",
            "MIN_BC_READ_THRESHOLD=0",
            f"CELL_BC_FILE={selected_cells}",
            "OUTPUT_HEADER=false",
            f"UEI={row.library}",
        ]
    )
    if row.locus_function_list == "intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=null", "LOCUS_FUNCTION_LIST=INTRONIC"])
    elif row.locus_function_list == "exonic+intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=INTRONIC"])
    run_command(cmd, "DigitalExpression", row.library)

    return selected_cells


@click.command("process_library")
@click.option(
    "-i",
    "--library-index",
    type=int,
    required=True,
    help="Which sample from the metadata to align",
)
@click.option(
    "--manifest-file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="YAML file containing the flowcell manifest",
)
@click.option("-d", "--dryrun", is_flag=True, help="Show the plan but don't execute")
@click.option("--debug", is_flag=True, help="Turn on debug logging")
@click.option("--log-file", type=click.Path(exists=False))
def main(
    library_index: int,
    manifest_file: str,
    dryrun: bool = False,
    debug: bool = False,
    log_file: str = None,
):
    create_logger(debug=debug, dryrun=dryrun, log_file=log_file)

    log.debug(f"Reading manifest from {manifest_file}")
    manifest = Manifest.from_file(Path(manifest_file))
    metadata_df = pd.read_csv(manifest.metadata_file)

    # task array is 1-indexed
    library = sorted(set(metadata_df.library))[library_index - 1]
    library_df = metadata_df.loc[metadata_df.library == library]
    log.debug(f"Processing alignments for library {library}")

    validate_library_df(library, library_df)
    row = library_df.iloc[0]

    lanes = sorted(set(library_df.lanes))

    library_dir = constants.LIBRARY_DIR / f"{row.date}_{library}"
    reference = Path(row.reference)
    reference_dir = library_dir / f"{reference.stem}.{row.locus_function_list}"

    library_bases = [
        f"{manifest.flowcell}.L{lane:03d}.{library}.{row.sample_barcode}.$"
        for lane in lanes
    ]
    alignment_bases = [
        library_dir / f"L{lane:03d}" / "alignment" / base
        for lane, base in zip(lanes, library_bases)
    ]

    bam_files = [
        (library_dir / f"L{lane:03d}" / base).with_suffix(".final.bam")
        for lane, base in zip(lanes, library_bases)
    ]
    alignment_stats = [
        base.with_suffix(".alignment_statistics.pickle") for base in alignment_bases
    ]
    star_logs = [base.with_suffix(".star.Log.final.out") for base in alignment_bases]

    # Combine check_alignments_quality files, and plot histograms
    alignment_quality.combine_alignment_stats(
        alignment_stats, star_logs, library_dir / f"{library}.$"
    )

    # paths for merged full BAM and BAM filtererd to matched barcodes
    combined_bam = library_dir / f"{library}.bam"
    matched_bam = combined_bam.with_suffix(".matched.bam")

    # Merge bam files
    cmd = picard_cmd("MergeSamFiles", manifest.tmp_dir)
    cmd.extend(
        [
            "--CREATE_INDEX",
            "true",
            "--CREATE_MD5_FILE",
            "false",
            "--OUTPUT",
            combined_bam,
            "--SORT_ORDER",
            "coordinate",
            "--ASSUMED_SORTED",
            "true",
        ]
    )

    for bam_file in bam_files:
        cmd.extend(["--INPUT", bam_file])

    run_command(cmd, "MergeBamFiles", library)

    # remove the individual bam files after merging
    for bam_file in bam_files:
        log.debug(f"Removing {bam_file}")
        os.remove(bam_file)

    # Validate bam file
    cmd = picard_cmd("ValidateSamFile", manifest.tmp_dir)
    cmd.extend(["--INPUT", combined_bam, "--MODE", "SUMMARY"])
    # is this necessary?
    if not row.illuminaplatform.startswith("NovaSeq"):
        cmd.extend(
            ["--IGNORE", "MISSING_PLATFORM_VALUE", "--IGNORE", "INVALID_VERSION_NUMBER"]
        )

    run_command(cmd, "ValidateSamFile", library)

    # generate various metrics files, including digital expression matrix
    selected_cells = calc_alignment_metrics(combined_bam, reference, row, manifest)

    if row.run_barcodematching:
        puckcaller_dir = Path(row.puckcaller_path)

        barcode_matching_folder = reference_dir / "barcode_matching"
        barcode_matching_folder.mkdir(exist_ok=True, parents=True)
        barcode_matching_file = (
            barcode_matching_folder / f"{row.library}_barcode_matching.txt.gz"
        )
        matched_barcodes_file = (
            barcode_matching_file / f"{row.library}_matched_barcodes.txt.gz"
        )

        barcode_mapping, bead_xy, bead_graph = match_barcodes(
            selected_cells,
            puckcaller_dir / "BeadBarcodes.txt",
            puckcaller_dir / "BeadLocations.txt",
        )

        with gzip.open(barcode_matching_file, "wt") as out:
            with gzip.open(matched_barcodes_file, "wt") as out2:
                for seq_bc, bead_bc in barcode_mapping.items():
                    x, y = bead_xy[bead_bc]
                    print(f"{seq_bc}\t{bead_bc}\t{x:.1f}\t{y:.1f}", file=out)
                    print(bead_bc, file=out2)

        # subset to the matched beads and add combined barcode as XB tag
        write_retagged_bam(combined_bam, matched_bam, barcode_mapping)

        # do it all again, but should be faster on the smaller file
        calc_alignment_metrics(
            matched_bam, reference, row, manifest, matched_barcodes_file, "XB"
        )
        os.remove(matched_barcodes_file)
    else:
        matched_bam = None
        bead_xy = None

    make_library_plots(row, lanes, manifest, matched_bam, bead_xy)

    if row.gen_downsampling:
        downsample_dir = reference_dir / "downsample"
        downsample_dir.mkdir(exist_ok=True, parents=True)

        # this might take a long time, we'll see...
        for ratio in np.linspace(0.1, 0.9, 9):
            downsample_dge(
                bam_file=combined_bam,
                downsample_dir=downsample_dir,
                row=row,
                ratio=ratio,
                tmp_dir=manifest.tmp_dir,
            )
