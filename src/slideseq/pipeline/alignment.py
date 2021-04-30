#!/usr/bin/env python

import logging
import os
import pathlib
import re

import click
import pandas as pd

import slideseq.util.constants as constants
from slideseq.pipeline.metadata import Manifest
from slideseq.util import dropseq_cmd, picard_cmd, run_command, start_popen
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
    flowcell = manifest.flowcell
    metadata_df = pd.read_csv(manifest.metadata_file)

    # task array is 1-indexed
    library = sorted(set(metadata_df.library))[library_index - 1]
    row = metadata_df.loc[(metadata_df.lane == lane) & (metadata_df.library == library)]
    if len(row) == 0:
        log.debug("Library not present in this lane, nothing to do here")

    assert len(row) == 1, f"More than one library specified by {lane} and {library}"
    row = row.iloc[0]  # convert single-row DataFrame to a Series

    tmp_dir = manifest.output_directory / "tmp"

    output_dir = constants.LIBRARY_DIR / f"{row.date}_{library}" / f"L{lane:03d}"

    reference = pathlib.Path(row.reference)
    reference_dir = reference.parent

    # assumes that reference always ends in .fasta, not .fasta.gz
    genome_dir = reference_dir / "STAR"
    intervals = reference_dir / f"{reference.stem}.genes.intervals"
    annotations_file = reference_dir / f"{reference.stem}.gtf"

    # define all the files we will produce
    bam_base = f"{flowcell}.L{lane:03d}.{library}.{row.sample_barcode}.$"
    bam_base = output_dir / bam_base

    # unmapped input. Keep this for posterity
    unmapped_bam = bam_base.with_suffix(".unmapped.bam")

    cellular_tagged_summary = bam_base.with_suffix(".cellular_tagging.summary.txt")
    molecular_tagged_summary = bam_base.with_suffix(".molecular_tagging.summary.txt")
    trimming_summary = bam_base.with_suffix(".adapter_trimming.summary.txt")

    polya_filtered_ubam = bam_base.with_suffix(
        ".unaligned_mc_tagged_polyA_filtered.bam"
    )
    polya_filtered_summary = bam_base.with_suffix(".polyA_filtering.summary.txt")

    # needed for STAR alignment, might as well keep it for later
    polya_filtered_fastq = polya_filtered_ubam.with_suffix(".fastq.gz")

    aligned_bam = bam_base.with_suffix(".star.Aligned.out.bam")
    alignment_statistics = bam_base.with_suffix(".alignment_statistics.pickle")

    # the final bam file
    final_aligned_bam = bam_base.with_suffix(".final.bam")

    bs_range1 = get_bead_structure_range(row.bead_structure, "C")
    bs_range2 = get_bead_structure_range(row.bead_structure, "M")

    procs = []

    cmd = dropseq_cmd("TagBamWithReadSequenceExtended", unmapped_bam, "/dev/stdout")
    cmd.extend(
        [
            f"SUMMARY={cellular_tagged_summary}",
            f"BASE_RANGE={bs_range1}",
            f"BASE_QUALITY={row.base_quality}",
            "BARCODED_READ=1",
            "DISCARD_READ=false",
            "TAG_NAME=XC",
            "NUM_BASES_BELOW_QUALITY=1",
        ]
    )

    procs.append(
        start_popen(
            cmd,
            "TagBamWithReadSequenceExtended (Cellular)",
            flowcell,
            library,
            lane,
        )
    )

    # Tag bam with read sequence extended molecular
    cmd = dropseq_cmd("TagBamWithReadSequenceExtended", "/dev/stdin", "/dev/stdout")
    cmd.extend(
        [
            f"SUMMARY={molecular_tagged_summary}",
            f"BASE_RANGE={bs_range2}",
            f"BASE_QUALITY={row.base_quality}",
            "BARCODED_READ=1",
            "DISCARD_READ=true",
            "TAG_NAME=XM",
            "NUM_BASES_BELOW_QUALITY=1",
        ]
    )
    procs.append(
        start_popen(
            cmd,
            "TagBamWithReadSequenceExtended (Molecular)",
            flowcell,
            library,
            lane,
            procs[-1],
        )
    )

    # Filter low-quality reads
    cmd = dropseq_cmd("FilterBam", "/dev/stdin", "/dev/stdout")
    cmd.extend(["PASSING_READ_THRESHOLD=0.1", "TAG_REJECT=XQ"])

    procs.append(start_popen(cmd, "FilterBam", flowcell, library, lane, procs[-1]))

    # Trim reads with starting sequence
    cmd = dropseq_cmd("TrimStartingSequence", "/dev/stdin", "/dev/stdout")
    cmd.extend(
        [
            f"OUTPUT_SUMMARY={trimming_summary}",
            f"SEQUENCE={row.start_sequence}",
            "MISMATCHES=0",
            "NUM_BASES=5",
        ]
    )

    procs.append(
        start_popen(cmd, "TrimStartingSequence", flowcell, library, lane, procs[-1])
    )

    # Adapter-aware poly A trimming
    cmd = dropseq_cmd("PolyATrimmer", "/dev/stdin", polya_filtered_ubam)
    cmd.extend(
        [
            f"OUTPUT_SUMMARY={polya_filtered_summary}",
            "MISMATCHES=0",
            "NUM_BASES=6",
            "USE_NEW_TRIMMER=true",
        ]
    )

    procs.append(start_popen(cmd, "PolyATrimmer", flowcell, library, lane, procs[-1]))

    # close intermediate streams
    for p in procs[:-1]:
        p.stdout.close()
    # wait for final process to finish
    log.debug(f"Finished with pre-alignment: {procs[-1].communicate()[0]}")

    # convert to fastq like a loser
    cmd = picard_cmd("SamToFastq", tmp_dir)
    cmd.extend(["-I", polya_filtered_ubam, "-F", polya_filtered_fastq])
    run_command(cmd, "SamToFastq", flowcell, library, lane)

    # Map reads to genome sequence using STARsolo
    cmd = [
        "STAR",
        "--genomeDir",
        genome_dir,
        "--readFilesIn",
        polya_filtered_fastq,
        "--readFilesCommand",
        "zcat",
        "--outFileNamePrefix",
        bam_base.with_suffix(".star."),
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

    run_command(cmd, "STAR", flowcell, library, lane)

    # Check alignments quality
    log.debug(
        f"Writing alignment statistics for {aligned_bam} to {alignment_statistics}"
    )
    write_alignment_stats(aligned_bam, alignment_statistics)

    procs = []

    # Sort aligned bam
    cmd = picard_cmd("SortSam", tmp_dir, mem="24g")
    cmd.extend(["-I", aligned_bam, "-O", "/dev/stdout", "--SORT_ORDER", "queryname"])
    procs.append(start_popen(cmd, "SortSam", flowcell, library, lane))

    # Merge unmapped bam and aligned bam
    cmd = picard_cmd("MergeBamAlignment", tmp_dir, mem="24g")
    cmd.extend(
        [
            "-R",
            row.reference,
            "--UNMAPPED",
            polya_filtered_ubam,
            "--ALIGNED",
            "/dev/stdin",
            "-O",
            "/dev/stdout",
            "--COMPRESSION_LEVEL",
            "0",
            "--INCLUDE_SECONDARY_ALIGNMENTS",
            "false",
            "--CLIP_ADAPTERS",
            "false",
        ]
    )

    procs.append(
        start_popen(cmd, "MergeBamAlignment", flowcell, library, lane, procs[-1])
    )

    # Tag read with interval
    cmd = dropseq_cmd("TagReadWithInterval", "/dev/stdin", "/dev/stdout")
    cmd.extend([f"INTERVALS={intervals}", "TAG=XG"])

    procs.append(
        start_popen(cmd, "TagReadWithInterval", flowcell, library, lane, procs[-1])
    )

    # Tag read with gene function
    cmd = dropseq_cmd("TagReadWithGeneFunction", "/dev/stdin", final_aligned_bam)
    cmd.extend([f"ANNOTATIONS_FILE={annotations_file}", "CREATE_INDEX=false"])

    procs.append(
        start_popen(cmd, "TagReadWithGeneFunction", flowcell, library, lane, procs[-1])
    )

    # close intermediate streams
    for p in procs[:-1]:
        p.stdout.close()

    # wait for final process to finish
    log.debug(f"Finished with post-alignment: {procs[-1].communicate()[0]}")

    os.remove(polya_filtered_ubam)
    os.remove(aligned_bam)

    log.info(f"Alignment for {library} completed")
