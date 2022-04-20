#!/usr/bin/python

import gzip
import logging
import os
import shutil
from pathlib import Path

import click

import slideseq.alignment_quality as alignment_quality
import slideseq.bead_matching as bead_matching
from slideseq.config import Config, get_config
from slideseq.library import Base, Library
from slideseq.metadata import Manifest
from slideseq.pipeline.write_matrix import write_sparse_matrix
from slideseq.plot.plot_library_metrics import make_library_plots
from slideseq.retag_bam import write_retagged_bam
from slideseq.util import give_group_access, rsync_to_google, run_command
from slideseq.util.logger import create_logger

log = logging.getLogger(__name__)


def calc_alignment_metrics(
    config: Config,
    input_base: Base,
    library: Library,
    tmp_dir: Path,
    matched_barcodes: Path = None,
    cell_tag: str = "XC",
):
    # Bam tag histogram (cells)
    cmd = config.dropseq_cmd(
        "BamTagHistogram", input_base.bam, input_base.reads_per_cell(cell_tag), tmp_dir
    )
    cmd.extend(
        [
            f"TAG={cell_tag}",
            "FILTER_PCR_DUPLICATES=false",
            f"READ_MQ={library.base_quality}",
        ]
    )

    run_command(cmd, "BamTagHistogram (cells)", library)

    # Bam tag histogram (UMIs)
    cmd = config.dropseq_cmd(
        "BamTagHistogram", input_base.bam, input_base.reads_per_umi, tmp_dir
    )
    cmd.extend(
        [
            "TAG=XM",
            "FILTER_PCR_DUPLICATES=false",
            f"READ_MQ={library.base_quality}",
        ]
    )

    run_command(cmd, "BamTagHistogram (UMIs)", library)

    # Collect RnaSeq metrics
    cmd = config.picard_cmd("CollectRnaSeqMetrics", tmp_dir)
    cmd.extend(
        [
            "--INPUT",
            input_base.bam,
            "--REF_FLAT",
            library.reference.ref_flat,
            "--OUTPUT",
            input_base.frac_intronic_exonic,
            "--STRAND_SPECIFICITY",
            "NONE",
            "--RIBOSOMAL_INTERVALS",
            library.reference.ribosomal_intervals,
        ]
    )

    run_command(cmd, "CollectRnaSeqMetrics", library)

    # Base distribution at read position for raw cellular barcode
    cmd = config.dropseq_cmd(
        "BaseDistributionAtReadPosition",
        input_base.bam,
        input_base.xc_distribution,
        tmp_dir,
    )
    cmd.extend(["TAG=XC"])
    run_command(cmd, "BaseDistributionAtReadPosition (Cellular)", library)

    # Base distribution at read position for molecular barcode
    cmd = config.dropseq_cmd(
        "BaseDistributionAtReadPosition",
        input_base.bam,
        input_base.xm_distribution,
        tmp_dir,
    )
    cmd.extend(["TAG=XM"])
    run_command(cmd, "BaseDistributionAtReadPosition (Molecular)", library)

    # Gather read quality metrics
    cmd = config.dropseq_cmd(
        "GatherReadQualityMetrics",
        input_base.bam,
        input_base.read_quality_metrics,
        tmp_dir,
    )
    run_command(cmd, "GatherReadQualityMetrics", library)

    if matched_barcodes is not None:
        cmd = config.dropseq_cmd(
            "SingleCellRnaSeqMetricsCollector",
            input_base.bam,
            input_base.frac_intronic_exonic_per_cell,
            tmp_dir,
            compression=6,
        )
        cmd.extend(
            [
                f"CELL_BARCODE_TAG={cell_tag}",
                f"READ_MQ={library.base_quality}",
                f"CELL_BC_FILE={matched_barcodes}",
                f"RIBOSOMAL_INTERVALS={library.reference.ribosomal_intervals}",
                f"ANNOTATIONS_FILE={library.reference.annotations}",
            ]
        )
        run_command(cmd, "SingleCellRnaSeqMetricsCollector", library)

    # Select cells by num transcripts
    cmd = config.dropseq_cmd(
        "SelectCellsByNumTranscripts",
        input_base.bam,
        input_base.selected_cells,
        tmp_dir,
        compression=6,
    )
    cmd.extend(
        [
            f"CELL_BARCODE_TAG={cell_tag}",
            f"MIN_TRANSCRIPTS_PER_CELL={library.min_transcripts_per_cell}",
            f"READ_MQ={library.base_quality}",
        ]
    )
    if library.locus_function_list == "intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=null", "LOCUS_FUNCTION_LIST=INTRONIC"])
    elif library.locus_function_list == "exonic+intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=INTRONIC"])
    run_command(cmd, "SelectCellsByNumTranscripts", library)

    # Generate digital expression files for all Illumina barcodes
    cmd = config.dropseq_cmd(
        "DigitalExpression",
        input_base.bam,
        input_base.digital_expression,
        tmp_dir,
        compression=6,
    )
    cmd.extend(
        [
            f"CELL_BARCODE_TAG={cell_tag}",
            f"CELL_BC_FILE={input_base.selected_cells}",
            f"SUMMARY={input_base.digital_expression_summary}",
            f"READ_MQ={library.base_quality}",
            "OUTPUT_HEADER=false",
            f"UEI={library}",
        ]
    )
    if library.locus_function_list == "intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=null", "LOCUS_FUNCTION_LIST=INTRONIC"])
    elif library.locus_function_list == "exonic+intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=INTRONIC"])
    run_command(cmd, "DigitalExpression", library)


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
    type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=Path),
    help="YAML file containing the manifest",
)
@click.option("--debug", is_flag=True, help="Turn on debug logging")
@click.option("--log-file", type=click.Path(path_type=Path))
def main(
    library_index: int, manifest_file: Path, debug: bool = False, log_file: Path = None
):
    create_logger(debug=debug, log_file=log_file)
    config = get_config()

    log.debug(f"Reading manifest from {manifest_file}")
    manifest = Manifest.from_file(manifest_file)

    # task array is 1-indexed
    library = manifest.get_library(library_index - 1)

    log.debug(f"Processing alignments for library {library.name}")

    # define barcode matching files, if needed
    barcode_matching_file = (
        library.barcode_matching_dir / f"{library}_barcode_matching.txt.gz"
    )
    barcode_coordinate_file = (
        library.barcode_matching_dir / f"{library}_barcode_xy.txt.gz"
    )
    matched_barcodes_file = (
        library.barcode_matching_dir / f"{library}_matched_barcodes.txt.gz"
    )

    # Combine check_alignments_quality files, and plot histograms
    alignment_quality.combine_alignment_stats(library)

    # Merge bam files
    cmd = config.picard_cmd("MergeSamFiles", manifest.tmp_dir)
    cmd.extend(
        [
            "--CREATE_INDEX",
            "true",
            "--CREATE_MD5_FILE",
            "false",
            "--OUTPUT",
            library.merged.bam,
            "--SORT_ORDER",
            "coordinate",
            "--ASSUME_SORTED",
            "true",
        ]
    )

    for bam_file in library.processed_bams:
        cmd.extend(["--INPUT", bam_file])

    run_command(cmd, "MergeBamFiles", library)

    # generate various metrics files, including digital expression matrix
    calc_alignment_metrics(config, library.merged, library, manifest.tmp_dir)

    if library.run_barcodematching:
        library.barcode_matching_dir.mkdir(exist_ok=True, parents=True)
        shutil.copy(library.bead_barcodes, library.barcode_matching_dir)
        shutil.copy(library.bead_locations, library.barcode_matching_dir)

        barcode_list, barcode_mapping, bead_xy, _ = bead_matching.match_barcodes(
            library.merged.selected_cells, library.bead_barcodes, library.bead_locations
        )

        bead_matching.write_barcode_mapping(
            barcode_mapping, bead_xy, barcode_matching_file
        )
        bead_matching.write_barcode_xy(barcode_list, bead_xy, barcode_coordinate_file)

        with gzip.open(matched_barcodes_file, "wt") as out:
            for bead_bc in sorted(set(barcode_mapping.values())):
                print(bead_bc, file=out)

        # subset to the matched beads and add combined barcode as XB tag
        write_retagged_bam(library.merged.bam, library.matched.bam, barcode_mapping)

        # do it all again, but should be faster on the smaller file
        calc_alignment_metrics(
            config,
            library.matched,
            library,
            manifest.tmp_dir,
            matched_barcodes_file,
            "XB",
        )

        write_sparse_matrix(library)

        make_library_plots(library, bead_xy)
    else:
        make_library_plots(library)

    # remove unneeded files now that we're done
    for bam_file in library.processed_bams:
        log.debug(f"Removing {bam_file}")
        os.remove(bam_file)

    if matched_barcodes_file.exists():
        log.debug(f"Removing {matched_barcodes_file}")
        os.remove(matched_barcodes_file)

    log.debug("Setting group permissions")
    give_group_access(library.dir)
    if config.gs_path is not None:
        log.debug("Copying data to google storage")
        rsync_to_google(library.dir, config.gs_path / library.date_name)

    log.info(f"Processing for {library} complete")
