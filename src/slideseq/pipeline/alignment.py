#!/usr/bin/env python

import logging
import os
from pathlib import Path

import click

import slideseq.util.constants as constants
from slideseq.alignment_quality import write_alignment_stats
from slideseq.config import get_config
from slideseq.metadata import Manifest
from slideseq.util import give_group_access, rsync_to_google, run_command, start_popen
from slideseq.util.logger import create_logger

log = logging.getLogger(__name__)


@click.command(name="align_library", no_args_is_help=True)
@click.option("--flowcell", required=True, help="Flowcell being aligned")
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
    help="YAML file containing the manifest",
)
@click.option("--debug", is_flag=True, help="Turn on debug logging")
@click.option("--log-file", type=click.Path(exists=False))
def main(
    flowcell: str,
    lane: int,
    library_index: int,
    manifest_file: str,
    debug: bool = False,
    log_file: str = None,
):
    create_logger(debug=debug, log_file=log_file)
    config = get_config()

    log.debug(f"Reading manifest from {manifest_file}")
    manifest = Manifest.from_file(Path(manifest_file))

    # task array is 1-indexed
    sample = manifest.get_sample(library_index - 1, flowcell, lane)
    if sample is None:
        return
    else:
        for barcode_ubam in sample.barcode_ubams:
            if not barcode_ubam.exists():
                raise ValueError(f"{barcode_ubam} does not exist")

    # create a subdirectory for alignment
    (sample.lane_dir / "alignment").mkdir(exist_ok=True)

    bead_structure = sample.get_bead_structure()
    xc_range = ":".join(f"{i}-{j}" for c, i, j in bead_structure if c == "C")
    xm_range = ":".join(f"{i}-{j}" for c, i, j in bead_structure if c == "M")

    # if multiple barcodes correspond to one sample, merge them
    if len(sample.barcode_ubams) > 1:
        cmd = config.picard_cmd("MergeSamFiles", manifest.tmp_dir)
        cmd.extend(
            [
                "--OUTPUT",
                sample.raw_ubam,
                "--SORT_ORDER",
                "unsorted",
                "--ASSUME_SORTED",
                "true",
            ]
        )

        for ubam_file in sample.barcode_ubams:
            cmd.extend(["--INPUT", ubam_file])

        run_command(cmd, "MergeBamFiles", sample)
    else:
        # otherwise, just rename the file
        os.rename(sample.barcode_ubams[0], sample.raw_ubam)

    procs = []

    cmd = config.dropseq_cmd(
        "TagBamWithReadSequenceExtended",
        sample.raw_ubam,
        "/dev/stdout",
        manifest.tmp_dir,
    )
    cmd.extend(
        [
            f"SUMMARY={sample.cellular_tagged_summary}",
            f"BASE_RANGE={xc_range}",
            f"BASE_QUALITY={sample.base_quality}",
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
            sample,
            lane,
        )
    )

    # Tag bam with read sequence extended molecular
    cmd = config.dropseq_cmd(
        "TagBamWithReadSequenceExtended", "/dev/stdin", "/dev/stdout", manifest.tmp_dir
    )
    cmd.extend(
        [
            f"SUMMARY={sample.molecular_tagged_summary}",
            f"BASE_RANGE={xm_range}",
            f"BASE_QUALITY={sample.base_quality}",
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
            sample,
            lane,
            procs[-1],
        )
    )

    # Filter low-quality reads
    cmd = config.dropseq_cmd("FilterBam", "/dev/stdin", "/dev/stdout", manifest.tmp_dir)
    cmd.extend(["PASSING_READ_THRESHOLD=0.1", "TAG_REJECT=XQ"])

    procs.append(start_popen(cmd, "FilterBam", sample, lane, procs[-1]))

    if sample.start_sequence != constants.NO_START_SEQUENCE:
        # Trim reads with starting sequence
        cmd = config.dropseq_cmd(
            "TrimStartingSequence", "/dev/stdin", "/dev/stdout", manifest.tmp_dir
        )
        cmd.extend(
            [
                f"OUTPUT_SUMMARY={sample.trimming_summary}",
                f"SEQUENCE={sample.start_sequence}",
                "MISMATCHES=0",
                "NUM_BASES=5",
            ]
        )

        procs.append(start_popen(cmd, "TrimStartingSequence", sample, lane, procs[-1]))

    # Adapter-aware poly A trimming
    cmd = config.dropseq_cmd(
        "PolyATrimmer", "/dev/stdin", sample.polya_filtered_ubam, manifest.tmp_dir
    )
    cmd.extend(
        [
            f"OUTPUT_SUMMARY={sample.polya_filtering_summary}",
            "MISMATCHES=0",
            "NUM_BASES=6",
            "USE_NEW_TRIMMER=true",
        ]
    )

    procs.append(start_popen(cmd, "PolyATrimmer", sample, lane, procs[-1]))

    # close intermediate streams
    for p in procs[:-1]:
        p.stdout.close()

    # wait for final process to finish
    procs[-1].communicate()
    log.debug("Finished with pre-alignment processing")

    # convert to fastq like a loser
    cmd = config.picard_cmd("SamToFastq", manifest.tmp_dir)
    cmd.extend(
        [
            "-I",
            sample.polya_filtered_ubam,
            "-F",
            sample.polya_filtered_fastq,
        ]
    )
    run_command(cmd, "SamToFastq", sample, lane)

    # Map reads to genome sequence using STAR
    cmd = [
        "STAR",
        "--genomeDir",
        sample.reference.genome_dir,
        "--readFilesIn",
        sample.polya_filtered_fastq,
        "--readFilesCommand",
        "zcat",
        "--outFileNamePrefix",
        sample.star_prefix,
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

    run_command(cmd, "STAR", sample, lane)

    # Check alignments quality
    log.debug(
        f"Writing alignment statistics for {sample.aligned_bam} to"
        f" {sample.alignment_pickle}"
    )
    write_alignment_stats(sample.aligned_bam, sample.alignment_pickle)

    procs = []

    # Sort aligned bam
    cmd = config.picard_cmd("SortSam", manifest.tmp_dir, mem="24g")
    cmd.extend(
        [
            "-I",
            sample.aligned_bam,
            "-O",
            "/dev/stdout",
            "--SORT_ORDER",
            "queryname",
        ]
    )
    procs.append(start_popen(cmd, "SortSam", sample, lane))

    # Merge unmapped bam and aligned bam
    cmd = config.picard_cmd("MergeBamAlignment", manifest.tmp_dir, mem="24g")
    cmd.extend(
        [
            "-R",
            sample.reference.fasta,
            "--UNMAPPED",
            sample.polya_filtered_ubam,
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

    procs.append(start_popen(cmd, "MergeBamAlignment", sample, lane, procs[-1]))

    # Tag read with interval
    cmd = config.dropseq_cmd(
        "TagReadWithInterval", "/dev/stdin", "/dev/stdout", manifest.tmp_dir
    )
    cmd.extend([f"INTERVALS={sample.reference.intervals}", "TAG=XG"])

    procs.append(start_popen(cmd, "TagReadWithInterval", sample, lane, procs[-1]))

    # Tag read with gene function
    cmd = config.dropseq_cmd(
        "TagReadWithGeneFunction",
        "/dev/stdin",
        sample.processed_bam,
        manifest.tmp_dir,
        compression=5,
    )
    cmd.extend(
        [f"ANNOTATIONS_FILE={sample.reference.annotations}", "CREATE_INDEX=false"]
    )

    procs.append(start_popen(cmd, "TagReadWithGeneFunction", sample, lane, procs[-1]))

    # close intermediate streams
    for p in procs[:-1]:
        p.stdout.close()

    # wait for final process to finish
    procs[-1].communicate()
    log.debug("Finished with post-alignment processing")

    # removed unneeded files
    os.remove(sample.polya_filtered_ubam)
    os.remove(sample.polya_filtered_fastq)
    os.remove(sample.aligned_bam)

    log.debug("Setting group permissions")
    give_group_access(sample.dir)
    log.debug("Copying data to google storage")
    rsync_to_google(
        sample.lane_dir,
        config.gs_path / sample.lane_dir.relative_to(config.library_dir),
    )

    log.info(f"Alignment for {sample} completed")
