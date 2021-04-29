#!/usr/bin/python

import logging
import os
import pathlib

import click
import pandas as pd

import slideseq.pipeline.analysis as analysis
import slideseq.util.alignment_quality as alignment_quality
import slideseq.util.constants as constants
from slideseq.pipeline.metadata import Manifest
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
    manifest = Manifest.from_file(pathlib.Path(manifest_file))
    metadata_df = pd.read_csv(manifest.metadata_file)
    tmp_dir = manifest.output_directory / "tmp"

    # task array is 1-indexed
    library = sorted(set(metadata_df.library))[library_index - 1]
    library_df = metadata_df.loc[metadata_df.library == library]
    log.debug(f"Processing alignments for library {library}")

    validate_library_df(library, library_df)
    row = library_df.iloc[0]

    lanes = sorted(set(library_df.lanes))

    library_dir = constants.LIBRARY_DIR / f"{row.date}_{library}"

    reference = pathlib.Path(row.reference)
    reference_dir = reference.parent

    # assumes that reference always ends in .fasta, not .fasta.gz
    ref_flat = reference_dir / f"{reference.stem}.refFlat"
    ribosomal_intervals = reference_dir / f"{reference.stem}.rRNA.intervals"

    bam_base = [
        f"{manifest.flowcell}.L{lane:03d}.{row.library}.{row.sample_barcode}"
        for lane in lanes
    ]
    bam_files = [library_dir / f"{base}.final.bam" for base in bam_base]
    alignment_stats = [
        library_dir / f"{base}.alignment_statistics.pickle" for base in bam_base
    ]
    star_logs = [library_dir / f"{base}.star.Log.final.out" for base in bam_base]

    # Merge bam files
    combined_bam = library_dir / f"{library}.bam"

    reads_per_cell = combined_bam.with_suffix(
        f".numReads_perCell_XC_mq_{row.base_quality}.txt.gz"
    )
    frac_intronic_exonic = combined_bam.with_suffix(".fracIntronicExonic.txt")
    xc_barcode_distribution = combined_bam.with_suffix(".barcode_distribution_XC.txt")
    xm_barcode_distribution = combined_bam.with_suffix(".barcode_distribution_XM.txt")
    read_quality_metrics = combined_bam.with_suffix(".ReadQualityMetrics.txt")

    selected_cells = combined_bam.with_suffix(
        f".{row.min_transcripts_per_cell}_transcripts_mq_{row.base_quality}_selected_cells.txt.gz"
    )
    digital_expression = combined_bam.with_suffix(
        ".AllIllumina.digital_expression.txt.gz"
    )
    digital_expression_summary = combined_bam.with_suffix(
        ".AllIllumina.digital_expression_summary.txt"
    )

    cmd = picard_cmd("MergeSamFiles", tmp_dir)
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

    run_command(cmd, "MergeBamFiles", manifest.flowcell, library)

    # remove the individual bam files after merging
    for bam_file in bam_files:
        log.debug(f"Removing {bam_file}")
        os.remove(bam_file)

    # Validate bam file
    cmd = picard_cmd("ValidateSamFile", tmp_dir)
    cmd.extend(["--INPUT", combined_bam, "--MODE", "SUMMARY"])
    # is this necessary?
    if not row.illuminaplatform.startswith("NovaSeq"):
        cmd.extend(
            ["--IGNORE", "MISSING_PLATFORM_VALUE", "--IGNORE", "INVALID_VERSION_NUMBER"]
        )

    run_command(cmd, "ValidateSamFile", manifest.flowcell, library)

    # Bam tag histogram
    cmd = dropseq_cmd("BamTagHistogram", combined_bam, reads_per_cell)
    cmd.extend(
        [
            "TAG=XC",
            "FILTER_PCR_DUPLICATES=false",
            f"READ_MQ={row.base_quality}",
            f"TMP_DIR={tmp_dir}",
        ]
    )

    run_command(cmd, "BamTagHistogram", manifest.flowcell, library)

    # Collect RnaSeq metrics
    cmd = picard_cmd("CollectRnaSeqMetrics", tmp_dir)
    cmd.extend(
        [
            "--INPUT",
            combined_bam,
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

    run_command(cmd, "CollectRnaSeqMetrics", manifest.flowcell, library)

    # Base distribution at read position for cellular barcode
    cmd = dropseq_cmd(
        "BaseDistributionAtReadPosition", combined_bam, xc_barcode_distribution
    )
    cmd.extend(["TAG=XC"])

    run_command(
        cmd, "BaseDistributionAtReadPosition (Cellular)", manifest.flowcell, library
    )

    # Base distribution at read position for molecular barcode
    cmd = dropseq_cmd(
        "BaseDistributionAtReadPosition", combined_bam, xm_barcode_distribution
    )
    cmd.extend(["TAG=XM"])

    run_command(
        cmd, "BaseDistributionAtReadPosition (Molecular)", manifest.flowcell, library
    )

    # Gather read quality metrics
    cmd = dropseq_cmd("GatherReadQualityMetrics", combined_bam, read_quality_metrics)

    run_command(cmd, "GatherReadQualityMetrics", manifest.flowcell, library)

    # Select cells by num transcripts
    cmd = dropseq_cmd("SelectCellsByNumTranscripts", combined_bam, selected_cells)
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

    run_command(cmd, "SelectCellsByNumTranscripts", manifest.flowcell, library)

    # Generate digital expression files for all Illumina barcodes
    cmd = dropseq_cmd("DigitalExpression", combined_bam, digital_expression)
    cmd.extend(
        [
            f"SUMMARY={digital_expression_summary}",
            "EDIT_DISTANCE=1",
            f"READ_MQ={row.base_quality}",
            "MIN_BC_READ_THRESHOLD=0",
            f"CELL_BC_FILE={selected_cells}",
            "OUTPUT_HEADER=false",
            f"UEI={library}",
        ]
    )
    if row.locus_function_list == "intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=null", "LOCUS_FUNCTION_LIST=INTRONIC"])
    elif row.locus_function_list == "exonic+intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=INTRONIC"])

    run_command(cmd, "DigitalExpression", library)

    # TODO: run_analysis_spec
    # barcode matching and downsampling

    if row.run_barcodematching:
        analysis.run_barcodematching(selected_cells, row, library_dir)

    # gen_downsampling
    if row.gen_downsampling:
        analysis.gen_downsampling()

    # Combine check_alignments_quality files
    alignment_quality.combine_alignment_stats(
        alignment_stats, star_logs, library_dir / f"{library}.$"
    )

    # TODO: plot alignment histogram
    # commandStr = "python {}/plot_alignment_histogram.py {} {} {}".format(
    #     scripts_folder, analysis_folder, library, library
    # )
    # os.system(commandStr)
