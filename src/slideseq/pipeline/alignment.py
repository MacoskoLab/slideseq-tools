#!/usr/bin/env python

import logging
import os
import pathlib
import re

import click
import pandas as pd

import slideseq.util.constants as constants
from slideseq.alignment_quality import write_alignment_stats
from slideseq.metadata import Manifest
from slideseq.util import dropseq_cmd, picard_cmd, run_command, start_popen
from slideseq.util.logger import create_logger

log = logging.getLogger(__name__)


def get_bead_structure(bead_structure: str):
    """Extract bead structure ranges from a format string. Take input like:
    8C18X6C9M1X and return:
      [('C', 1, 8), ('X', 9, 26), ('C', 27, 32), ('M', 33, 41), ('X', 42, 42)]

    :param bead_structure: the structure to parse, using C for cell barcode, M
                           for UMI, and X for spacers/fixed sequences
    """
    i = 1
    intervals = []
    # this is obfuscated to hell but I can't help myself
    #   - regex splits on C, X or M but retains the character found
    #   - 2 * (iter(re),) makes a second reference (not copy) to the iterator
    #   - zip(*(2iter)) zips together pairs of elements
    # then, add up the lengths to get ranges
    for j, c in zip(*(2 * (iter(re.split("([CXM])", bead_structure)),))):
        j = int(j)
        intervals.append((c, i, i + j - 1))
        i += j

    return intervals


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
        log.warning("Library not present in this lane, nothing to do here")
        return

    assert len(row) == 1, f"More than one library specified by {lane} and {library}"
    row = row.iloc[0]  # convert single-row DataFrame to a Series

    library_dir = constants.LIBRARY_DIR / f"{row.date}_{library}" / f"L{lane:03d}"

    reference = pathlib.Path(row.reference)

    # assumes that reference always ends in .fasta, not .fasta.gz
    genome_dir = reference.parent / "STAR"
    intervals = reference.with_suffix(".genes.intervals")
    annotations_file = reference.with_suffix(".gtf")

    # define the base name for the library we're aligning here
    library_base = f"{manifest.flowcell}.L{lane:03d}.{library}.{row.sample_barcode}.$"

    # unmapped input. we'll keep this file for posterity
    unmapped_bam = (library_dir / library_base).with_suffix(".unmapped.bam")

    # create a subdirectory for alignment
    (library_dir / "alignment").mkdir(exist_ok=True)
    alignment_base = library_dir / "alignment" / library_base

    cellular_tagged_summary = alignment_base.with_suffix(
        ".cellular_tagging.summary.txt"
    )
    molecular_tagged_summary = alignment_base.with_suffix(
        ".molecular_tagging.summary.txt"
    )
    trimming_summary = alignment_base.with_suffix(".adapter_trimming.summary.txt")

    polya_filtered_ubam = alignment_base.with_suffix(
        ".unaligned_mc_tagged_polyA_filtered.bam"
    )
    polya_filtered_summary = alignment_base.with_suffix(".polyA_filtering.summary.txt")

    # needed for STAR alignment
    polya_filtered_fastq = polya_filtered_ubam.with_suffix(".fastq.gz")

    aligned_bam = alignment_base.with_suffix(".star.Aligned.out.bam")
    alignment_statistics = alignment_base.with_suffix(".alignment_statistics.pickle")

    # the final bam file
    final_aligned_bam = (library_dir / library_base).with_suffix(".final.bam")

    bead_structure = get_bead_structure(row.bead_structure)

    xc_range = ":".join(f"{i}-{j}" for c, i, j in bead_structure if c == "C")
    xm_range = ":".join(f"{i}-{j}" for c, i, j in bead_structure if c == "M")

    procs = []

    cmd = dropseq_cmd("TagBamWithReadSequenceExtended", unmapped_bam, "/dev/stdout")
    cmd.extend(
        [
            f"SUMMARY={cellular_tagged_summary}",
            f"BASE_RANGE={xc_range}",
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
            library,
            lane,
        )
    )

    # Tag bam with read sequence extended molecular
    cmd = dropseq_cmd("TagBamWithReadSequenceExtended", "/dev/stdin", "/dev/stdout")
    cmd.extend(
        [
            f"SUMMARY={molecular_tagged_summary}",
            f"BASE_RANGE={xm_range}",
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
            library,
            lane,
            procs[-1],
        )
    )

    # Filter low-quality reads
    cmd = dropseq_cmd("FilterBam", "/dev/stdin", "/dev/stdout")
    cmd.extend(["PASSING_READ_THRESHOLD=0.1", "TAG_REJECT=XQ"])

    procs.append(start_popen(cmd, "FilterBam", library, lane, procs[-1]))

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

    procs.append(start_popen(cmd, "TrimStartingSequence", library, lane, procs[-1]))

    # Adapter-aware poly A trimming
    cmd = dropseq_cmd("PolyATrimmer", "/dev/stdin", polya_filtered_ubam, compression=5)
    cmd.extend(
        [
            f"OUTPUT_SUMMARY={polya_filtered_summary}",
            "MISMATCHES=0",
            "NUM_BASES=6",
            "USE_NEW_TRIMMER=true",
        ]
    )

    procs.append(start_popen(cmd, "PolyATrimmer", library, lane, procs[-1]))

    # close intermediate streams
    for p in procs[:-1]:
        p.stdout.close()

    # wait for final process to finish
    procs[-1].communicate()
    log.debug("Finished with pre-alignment processing")

    # convert to fastq like a loser
    cmd = picard_cmd("SamToFastq", manifest.tmp_dir)
    cmd.extend(["-I", polya_filtered_ubam, "-F", polya_filtered_fastq])
    run_command(cmd, "SamToFastq", library, lane)

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
        alignment_base.with_suffix(".star."),
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

    run_command(cmd, "STAR", library, lane)

    # Check alignments quality
    log.debug(
        f"Writing alignment statistics for {aligned_bam} to {alignment_statistics}"
    )
    write_alignment_stats(aligned_bam, alignment_statistics)

    procs = []

    # Sort aligned bam
    cmd = picard_cmd("SortSam", manifest.tmp_dir, mem="24g")
    cmd.extend(["-I", aligned_bam, "-O", "/dev/stdout", "--SORT_ORDER", "queryname"])
    procs.append(start_popen(cmd, "SortSam", library, lane))

    # Merge unmapped bam and aligned bam
    cmd = picard_cmd("MergeBamAlignment", manifest.tmp_dir, mem="24g")
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

    procs.append(start_popen(cmd, "MergeBamAlignment", library, lane, procs[-1]))

    # Tag read with interval
    cmd = dropseq_cmd("TagReadWithInterval", "/dev/stdin", "/dev/stdout")
    cmd.extend([f"INTERVALS={intervals}", "TAG=XG"])

    procs.append(start_popen(cmd, "TagReadWithInterval", library, lane, procs[-1]))

    # Tag read with gene function
    cmd = dropseq_cmd(
        "TagReadWithGeneFunction", "/dev/stdin", final_aligned_bam, compression=5
    )
    cmd.extend([f"ANNOTATIONS_FILE={annotations_file}", "CREATE_INDEX=false"])

    procs.append(start_popen(cmd, "TagReadWithGeneFunction", library, lane, procs[-1]))

    # close intermediate streams
    for p in procs[:-1]:
        p.stdout.close()

    # wait for final process to finish
    procs[-1].communicate()
    log.debug("Finished with post-alignment processing")

    # removed unneeded files
    os.remove(polya_filtered_ubam)
    os.remove(polya_filtered_fastq)
    os.remove(aligned_bam)

    log.info(f"Alignment for {library} completed")
