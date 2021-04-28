#!/usr/bin/env python

import logging
import os
import pathlib
import re

import click
import pandas as pd

import slideseq.util.constants as constants
from slideseq.pipeline.metadata import Manifest
from slideseq.util import picard_cmd, run_command
from slideseq.util.alignment_quality import write_alignment_stats
from slideseq.util.logger import create_logger

log = logging.getLogger(__name__)


# Get bead structure range
def get_bead_structure_range(bs, structure_type):
    # 12C8M|*T
    # 7C18X7C8M2X|*T
    ell = re.split("[CXM]", bs.split("|")[0])
    res = ""
    i = 1
    p = -1
    for it in ell:
        if it:
            p += len(it) + 1
            if bs[p] == structure_type:
                res += str(i) + "-" + str(i + int(it) - 1) + ":"
            i += int(it)
    return res[:-1]


@click.command(name="align_library", no_args_is_help=True)
@click.option(
    "--lane", type=int, required=True, help="Lane of the flowcell being aligned"
)
@click.option(
    "-i",
    "--library-index",
    type=int,
    required=True,
    help="Which library from the metadata to align",
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
    lane: int,
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

    # task array is 1-indexed
    library = sorted(set(metadata_df.library))[library_index - 1]
    row = metadata_df.loc[(metadata_df.lane == lane) & (metadata_df.library == library)]
    if len(row) == 0:
        log.debug("Library not present in this lane, nothing to do here")

    assert len(row) == 1, f"More than one library specified by {lane} and {library}"
    row = row.iloc[0]  # convert single-row DataFrame to a Series

    tmp_dir = manifest.output_directory / "tmp"

    output_dir = constants.LIBRARY_DIR / f"{row.date}_{row.library}" / f"L{lane:03d}"

    reference = pathlib.Path(row.reference)
    reference_dir = reference.parent

    # assumes that reference always ends in .fasta, not .fasta.gz
    genome_dir = reference_dir / "STAR"
    intervals = reference_dir / f"{reference.stem}.genes.intervals"
    annotations_file = reference_dir / f"{reference.stem}.gtf"

    # define all the intermediate files we need
    # TODO: figure out if we can combine some of these steps, maybe
    bam_base = f"{manifest.flowcell}.L{lane:03d}.{row.library}.{row.sample_barcode}.$"
    bam_base = output_dir / bam_base

    unmapped_bam = bam_base.with_suffix(".unmapped.bam")

    cellular_tagged_bam = bam_base.with_suffix(".unmapped_tagged_cellular.bam")
    cellular_tagged_summary = bam_base.with_suffix(".cellular_tagging.summary.txt")

    molecular_tagged_bam = bam_base.with_suffix(".unmapped_tagged_molecular.bam")
    molecular_tagged_summary = bam_base.with_suffix(".molecular_tagging.summary.txt")

    filtered_ubam = bam_base.with_suffix(".unmapped.filtered.bam")
    trimmed_ubam = bam_base.with_suffix(".unmapped_trimstartingsequence.filtered.bam")

    trimming_summary = bam_base.with_suffix(".adapter_trimming.summary.txt")

    polya_filtered_ubam = bam_base.with_suffix(
        ".unaligned_mc_tagged_polyA_filtered.bam"
    )
    polya_filtered_summary = bam_base.with_suffix(".polyA_filtering.summary.txt")
    polya_filtered_fastq = polya_filtered_ubam.with_suffix(".fastq.gz")

    # intermediate files
    aligned_bam = bam_base.with_suffix(".star.Aligned.out.bam")
    alignment_statistics = bam_base.with_suffix(".alignment_statistics.pickle")
    aligned_sorted_bam = bam_base.with_suffix(".aligned.sorted.bam")
    aligned_merged_bam = bam_base.with_suffix(".merged.bam")
    aligned_tagged_bam = bam_base.with_suffix(".merged.TagReadWithInterval.bam")

    # the final file
    final_aligned_bam = bam_base.with_suffix(".final.bam")

    bs_range1 = get_bead_structure_range(row.bead_structure, "C")
    bs_range2 = get_bead_structure_range(row.bead_structure, "M")

    cmd = [
        f"{constants.DROPSEQ_DIR / 'TagBamWithReadSequenceExtended'}",
        f"I={unmapped_bam}",
        f"O={cellular_tagged_bam}",
        f"SUMMARY={cellular_tagged_summary}",
        f"BASE_RANGE={bs_range1}",
        f"BASE_QUALITY={row.base_quality}",
        "BARCODED_READ=1",
        "DISCARD_READ=false",
        "TAG_NAME=XC",
        "NUM_BASES_BELOW_QUALITY=1",
    ]

    run_command(
        cmd,
        "TagBamWithReadSequenceExtended (Cellular)",
        manifest.flowcell,
        row.library,
        lane,
    )

    # Tag bam with read sequence extended molecular
    cmd = [
        f"{constants.DROPSEQ_DIR / 'TagBamWithReadSequenceExtended'}",
        f"I={cellular_tagged_bam}",
        f"O={molecular_tagged_bam}",
        f"SUMMARY={molecular_tagged_summary}",
        f"BASE_RANGE={bs_range2}",
        f"BASE_QUALITY={row.base_quality}",
        "BARCODED_READ=1",
        "DISCARD_READ=true",
        "TAG_NAME=XM",
        "NUM_BASES_BELOW_QUALITY=1",
    ]

    run_command(
        cmd,
        "TagBamWithReadSequenceExtended (Molecular)",
        manifest.flowcell,
        row.library,
        lane,
    )
    os.remove(cellular_tagged_bam)

    # Filter low-quality reads
    cmd = [
        f"{constants.DROPSEQ_DIR / 'FilterBam'}",
        f"I={molecular_tagged_bam}",
        f"O={filtered_ubam}",
        "PASSING_READ_THRESHOLD=0.1",
        "TAG_REJECT=XQ",
    ]

    run_command(cmd, "FilterBam", manifest.flowcell, row.library, lane)
    os.remove(molecular_tagged_bam)

    # Trim reads with starting sequence
    cmd = [
        f"{constants.DROPSEQ_DIR / 'TrimStartingSequence'}",
        f"I={filtered_ubam}",
        f"O={trimmed_ubam}",
        f"OUTPUT_SUMMARY={trimming_summary}",
        f"SEQUENCE={row.start_sequence}",
        "MISMATCHES=0",
        "NUM_BASES=5",
    ]

    run_command(cmd, "TrimStartingSequence", manifest.flowcell, row.library, lane)
    os.remove(filtered_ubam)

    # Adapter-aware poly A trimming
    cmd = [
        f"{constants.DROPSEQ_DIR / 'PolyATrimmer'}",
        f"I={trimmed_ubam}",
        f"O={polya_filtered_ubam}",
        f"OUTPUT_SUMMARY={polya_filtered_summary}",
        "MISMATCHES=0",
        "NUM_BASES=6",
        "USE_NEW_TRIMMER=true",
    ]

    run_command(cmd, "PolyATrimmer", manifest.flowcell, row.library, lane)
    os.remove(trimmed_ubam)

    # convert to fastq like a loser
    cmd = picard_cmd("SamToFastq", tmp_dir)
    cmd.extend(
        [
            "-I",
            f"{polya_filtered_ubam}",
            "-F",
            f"{polya_filtered_fastq}",
        ]
    )
    run_command(cmd, "SamToFastq", manifest.flowcell, row.library, lane)

    # Map reads to genome sequence using STARsolo
    cmd = [
        "STAR",
        "--genomeDir",
        f"{genome_dir}",
        "--readFilesIn",
        f"{polya_filtered_fastq}",
        "--readFilesCommand",
        "zcat",
        "--outFileNamePrefix",
        f"{bam_base.with_suffix('.star.')}",
        "--outStd",
        "Log",
        "--outSAMtype",
        "BAM",
        "Unsorted",
        "--outBAMcompression",
        "0",
        "--limitOutSJcollapsed",
        "5000000",
        "--runThreadN",
        "8",
    ]

    run_command(cmd, "STAR", manifest.flowcell, row.library, lane)

    # Check alignments quality
    log.debug(
        f"Writing alignment statistics for {aligned_bam} to {alignment_statistics}"
    )
    write_alignment_stats(aligned_bam, alignment_statistics)

    # Sort aligned bam
    cmd = picard_cmd("SortSam", tmp_dir)
    cmd.extend(
        [
            "-I",
            f"{aligned_bam}",
            "-O",
            f"{aligned_sorted_bam}",
            "--SORT_ORDER",
            "queryname",
        ]
    )
    run_command(cmd, "SortSam", manifest.flowcell, row.library, lane)
    os.remove(aligned_bam)

    # Merge unmapped bam and aligned bam
    cmd = picard_cmd("MergeBamAlignment", tmp_dir)
    cmd.extend(
        [
            "-R",
            f"{row.reference}",
            "--UNMAPPED",
            f"{polya_filtered_ubam}",
            "--ALIGNED",
            f"{aligned_sorted_bam}",
            "-O",
            f"{aligned_merged_bam}",
            "--COMPRESSION_LEVEL",
            "0",
            "--INCLUDE_SECONDARY_ALIGNMENTS",
            "false",
            "--CLIP_ADAPTERS",
            "false",
        ]
    )

    run_command(cmd, "MergeBamAlignment", manifest.flowcell, row.library, lane)
    os.remove(polya_filtered_ubam)
    os.remove(aligned_sorted_bam)

    # Tag read with interval
    cmd = [
        f"{constants.DROPSEQ_DIR / 'TagReadWithInterval'}",
        f"I={aligned_merged_bam}",
        f"O={aligned_tagged_bam}",
        f"INTERVALS={intervals}",
        "TAG=XG",
    ]

    run_command(cmd, "TagReadWithInterval", manifest.flowcell, row.library, lane)
    os.remove(aligned_merged_bam)

    # Tag read with gene function
    cmd = [
        f"{constants.DROPSEQ_DIR / 'TagReadWithGeneFunction'}",
        f"I={aligned_tagged_bam}",
        f"O={final_aligned_bam}",
        f"ANNOTATIONS_FILE={annotations_file}",
        "CREATE_INDEX=false",
    ]

    run_command(cmd, "TagReadWithGeneFunction", manifest.flowcell, row.library, lane)
    os.remove(aligned_tagged_bam)

    log.info(f"Alignment for {row.library} completed")
