#!/usr/bin/python

import logging
import os
import pathlib

import click
import pandas as pd

import slideseq.util.alignment_quality as alignment_quality
import slideseq.util.constants as constants
from slideseq.pipeline.metadata import Manifest
from slideseq.util import picard_cmd, run_command
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
    lane_str = f"{min(lanes)}-{max(lanes)}"

    output_dir = constants.LIBRARY_DIR / f"{row.date}_{library}"

    bam_base = [
        f"{manifest.flowcell}.L{lane:03d}.{row.library}.{row.sample_barcode}"
        for lane in lanes
    ]
    bam_files = [output_dir / f"{base}.final.bam" for base in bam_base]
    alignment_stats = [
        output_dir / f"{base}.alignment_statistics.pickle" for base in bam_base
    ]
    star_logs = [output_dir / f"{base}.star.Log.final.out" for base in bam_base]

    # Merge bam files
    combined_bam = output_dir / f"{library}.bam"

    cmd = picard_cmd("MergeSamFiles", tmp_dir)
    cmd.extend(
        [
            "--CREATE_INDEX",
            "true",
            "--CREATE_MD5_FILE",
            "false",
            "--OUTPUT",
            f"{combined_bam}",
            "--SORT_ORDER",
            "coordinate",
            "--ASSUMED_SORTED",
            "true",
        ]
    )

    for bam_file in bam_files:
        cmd.extend(["--INPUT", f"{bam_file}"])

    run_command(cmd, "MergeBamFiles", manifest.flowcell, library, lane_str)

    # Validate bam file
    cmd = picard_cmd("ValidateSamFile", tmp_dir)
    cmd.extend(
        [
            "--INPUT",
            f"{combined_bam}",
            "--MODE",
            "SUMMARY",
        ]
    )
    # is this necessary?
    if not row.illuminaplatform.startswith("NovaSeq"):
        cmd.extend(
            ["--IGNORE", "MISSING_PLATFORM_VALUE", "--IGNORE", "INVALID_VERSION_NUMBER"]
        )

    run_command(cmd, "ValidateSamFile", manifest.flowcell, library, lane_str)

    # TODO: generate_plots
    # output_file = f"{output_dir}/logs/generate_plots_{library}.log"
    # submission_script = f"{scripts_folder}/generate_plots.sh"
    # call_args = [
    #     "qsub",
    #     "-o",
    #     output_file,
    #     submission_script,
    #     manifest_file,
    #     library,
    #     scripts_folder,
    #     output_dir,
    #     analysis_folder,
    # ]
    # call(call_args)

    # reference = pathlib.Path(row.reference)
    # reference_dir = reference.parent

    # assumes that reference always ends in .fasta, not .fasta.gz
    # genome_dir = reference_dir / "STAR"
    # intervals = reference_dir / f"{reference.stem}.genes.intervals"
    # annotations_file = reference_dir / f"{reference.stem}.gtf"

    # not gonna support multiple values in this column, for now
    # lists = locus_function_list.split(",")

    # TODO: make output dirs?
    # for j in lists:
    #     call(["mkdir", "-p", f"{analysis_folder}/{reference.stem}.{j}"])
    #     call(
    #         [
    #             "mkdir",
    #             "-p",
    #             f"{analysis_folder}/{reference.stem}.{j}/alignment",
    #         ]
    #     )
    #
    #     if run_barcodematching:
    #         barcode_matching_folder = (
    #             f"{analysis_folder}/{reference.stem}.{j}/barcode_matching/"
    #         )
    #         call(["mkdir", "-p", barcode_matching_folder])
    #         for i in range(len(lanes)):
    #             if libraries[i] != library:
    #                 continue
    #             for lane_slice in slice_id[lanes[i]]:
    #                 toCopyFile = f"{analysis_folder}/{flowcell_barcode}.{lanes[i]}.{lane_slice}.{library}"
    #                 if barcodes[i]:
    #                     toCopyFile += "." + barcodes[i]
    #                 toCopyFile += ".star_gene_exon_tagged2.bam"
    #                 if os.path.isfile(toCopyFile):
    #                     call(["cp", toCopyFile, barcode_matching_folder])

    # TODO: run_analysis_spec
    # output_file = f"{output_dir}/logs/run_analysis_spec_{library}_{j}.log"
    # submission_script = f"{scripts_folder}/run_analysis_spec.sh"
    # call_args = [
    #     "qsub",
    #     "-o",
    #     output_file,
    #     submission_script,
    #     manifest_file,
    #     library,
    #     scripts_folder,
    #     j,
    #     output_dir,
    #     f"{analysis_folder}/{reference.stem}.{j}",
    # ]
    # call(call_args)

    for bam_file in bam_files:
        log.debug(f"Removing {bam_file}")
        os.remove(bam_file)

    # Combine check_alignments_quality files
    alignment_quality.combine_alignment_stats(
        alignment_stats,
        star_logs,
        output_dir / f"{library}.$",
    )

    # TODO: plot alignment histogram
    # commandStr = "python {}/plot_alignment_histogram.py {} {} {}".format(
    #     scripts_folder, analysis_folder, library, library
    # )
    # os.system(commandStr)
