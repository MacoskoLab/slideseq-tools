#!/usr/bin/env python

import logging
import os
import pathlib
import re
import sys
from subprocess import run

import click
import pandas as pd

import slideseq.util.constants as constants
from slideseq.logger import create_logger
from slideseq.pipeline.metadata import Manifest
from slideseq.util import picard_cmd

log = logging.getLogger(__name__)


def run_command(cmd: list[str], name: str, flowcell: str, library: str, lane: int):
    log.info(f"{flowcell} - {name} for {library} in lane {lane}")
    log.debug(f"Command = {' '.join(cmd)}")
    proc = run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        log.error(f"Error running {name}:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info(f"{flowcell} - {name} completed")


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


@click.command(name="align_sample", no_args_is_help=True)
@click.option(
    "--lane", type=int, required=True, help="Lane of the flowcell being aligned"
)
@click.option(
    "-i",
    "--sample_index",
    type=int,
    required=True,
    help="Which sample from the metadata to align",
)
@click.option(
    "--manifest_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="YAML file containing the flowcell manifest",
)
@click.option("-d", "--dryrun", is_flag=True, help="Show the plan but don't execute")
@click.option("--debug", is_flag=True, help="Turn on debug logging")
@click.option("--log_file", type=click.Path(exists=False))
def main(
    lane: int,
    sample_index: int,
    manifest_file: str,
    dryrun: bool = False,
    debug: bool = False,
    log_file: str = None,
):
    create_logger(debug=debug, log_file=log_file)

    if dryrun:
        log.info("DRY RUN ONLY -- No files will be written and no jobs run")

    log.debug(f"Reading manifest from {manifest_file}")
    manifest = Manifest.from_file(pathlib.Path(manifest_file))
    metadata_df = pd.read_csv(manifest.metadata_file)

    tmp_dir = manifest.output_directory / "tmp"

    row = metadata_df.iloc[sample_index - 1]

    output_dir = constants.LIBRARY_DIR / f"{row.date}_{row.library}" / f"L{lane:03d}"

    reference = pathlib.Path(row.reference)
    reference_folder = reference.parent

    # assumes that reference always ends in .fasta, not .fasta.gz
    genome_dir = reference_folder / "STAR"
    intervals = reference_folder / f"{reference.stem}.genes.intervals"
    annotations_file = reference_folder / f"{reference.stem}.gtf"

    # define all the intermediate files we need
    # TODO: figure out if we can combine some of these steps, maybe
    bam_base = f"{manifest.flowcell}.L{lane:03d}.{row.library}.{row.sample_barcode}"
    unmapped_bam = output_dir / f"{bam_base}.unmapped.bam"

    cellular_tagged_bam = output_dir / f"{bam_base}.unmapped_tagged_cellular.bam"
    cellular_tagged_summary = output_dir / f"{bam_base}.cellular_tagging.summary.txt"

    molecular_tagged_bam = output_dir / f"{bam_base}.unmapped_tagged_molecular.bam"
    molecular_tagged_summary = output_dir / f"{bam_base}.molecular_tagging.summary.txt"

    filtered_ubam = output_dir / f"{bam_base}.unmapped.filtered.bam"
    trimmed_ubam = output_dir / f"{bam_base}.unmapped_trimstartingsequence.filtered.bam"

    trimming_summary = output_dir / f"{bam_base}.adapter_trimming.summary.txt"

    polya_filtered_ubam = (
        output_dir / f"{bam_base}.unaligned_mc_tagged_polyA_filtered.bam"
    )
    polya_filtered_summary = output_dir / f"{bam_base}.polyA_filtering.summary.txt"
    polya_filtered_fastq = polya_filtered_ubam.with_suffix("fastq")

    # prefix for aligned bam file
    aligned_bam = output_dir / f"{bam_base}"

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
        "--outFileNamePrefix",
        f"{aligned_bam}.star.",
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
    os.remove(polya_filtered_fastq)

    # TODO: this thing

    # Check alignments quality
    # star_file2 = f"{prefix_libraries}.star.Aligned.out.sam"
    # commandStr = f"samtools view -h -o {star_file2} {star_file}"
    # os.system(commandStr)
    # commandStr = f"check_alignments_quality {star_file2}"
    # log.info(
    # f"{manifest.flowcell} - Check alignments quality for {library} in Lane {lane} is done."
    # )
    # log.debug(f"Command = {commandStr}")
    # os.system(commandStr)
    # log.info(
    #     f"{manifest.flowcell} - Check alignments quality for"
    #     f" {sample_row.library} in Lane {lane} is done."
    # )
    # call(["rm", star_file2])

    # Sort aligned bam
    cmd = picard_cmd("SortSam", tmp_dir)
    cmd.extend(
        [
            "-I",
            f"{aligned_bam}.star.Aligned.out.bam",
            "-O",
            f"{aligned_bam}.aligned.sorted.bam",
            "--SORT_ORDER",
            "queryname",
        ]
    )
    run_command(cmd, "SortSam", manifest.flowcell, row.library, lane)
    os.remove(f"{aligned_bam}.star.Aligned.out.bam")

    # Merge unmapped bam and aligned bam
    cmd = picard_cmd("MergeBamAlignment", tmp_dir)
    cmd.extend(
        [
            "-R",
            f"{row.reference}",
            "--UNMAPPED",
            f"{polya_filtered_ubam}",
            "--ALIGNED",
            f"{aligned_bam}.aligned.sorted.bam",
            "-O",
            f"{aligned_bam}.merged.bam",
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
    os.remove(f"{aligned_bam}.aligned.sorted.bam")

    # Tag read with interval
    cmd = [
        f"{constants.DROPSEQ_DIR / 'TagReadWithInterval'}",
        f"I={aligned_bam}.merged.bam",
        f"O={aligned_bam}.merged.TagReadWithInterval.bam",
        f"INTERVALS={intervals}",
        "TAG=XG",
    ]

    run_command(cmd, "TagReadWithInterval", manifest.flowcell, row.library, lane)
    os.remove(f"{aligned_bam}.merged.bam")

    # Tag read with gene function
    cmd = [
        f"{constants.DROPSEQ_DIR / 'TagReadWithGeneFunction'}",
        f"I={aligned_bam}.merged.TagReadWithInterval.bam",
        f"O={aligned_bam}.star_gene_exon_tagged2.bam",
        f"ANNOTATIONS_FILE={annotations_file}",
        "CREATE_INDEX=false",
    ]

    run_command(cmd, "TagReadWithGeneFunction", manifest.flowcell, row.library, lane)
    os.remove(f"{aligned_bam}.merged.TagReadWithInterval.bam")

    log.info(f"Alignment for {row.library} completed")
