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
from slideseq.util import dropseq_cmd, picard_cmd

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

    sample_row = metadata_df.iloc[sample_index]

    library_dir = constants.LIBRARY_DIR / f"{sample_row.date}_{sample_row.library}"

    reference = sample_row.reference

    reference_folder = reference[: reference.rfind("/")]
    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]
    genome_dir = f"{reference_folder}/STAR"
    intervals = f"{reference_folder}/{referencePure}.genes.intervals"
    annotations_file = f"{reference_folder}/{referencePure}.gtf"

    output_dir = (
        library_dir / f"{lane:03d}" / sample_row.library / sample_row.sample_barcode
    )

    # define all the intermediate files we need
    # TODO: figure out if we can combine some of these steps, maybe
    bam_base = (
        f"{manifest.flowcell}.{lane}.{sample_row.library}.{sample_row.sample_barcode}"
    )
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

    # prefix for aligned bam file
    aligned_bam = output_dir / f"{bam_base}"

    bs_range1 = get_bead_structure_range(sample_row.bead_structure, "C")
    bs_range2 = get_bead_structure_range(sample_row.bead_structure, "M")

    cmd = dropseq_cmd("TagBamWithReadSequenceExtended", tmp_dir)
    cmd.extend(
        [
            "-I",
            f"{unmapped_bam}",
            "-O",
            f"{cellular_tagged_bam}",
            "--SUMMARY",
            f"{cellular_tagged_summary}",
            "--BASE_RANGE",
            f"{bs_range1}",
            "--BASE_QUALITY",
            f"{sample_row.base_quality}",
            "--BARCODED_READ",
            "1",
            "--DISCARD_READ",
            "false",
            "--TAG_NAME",
            "XC",
            "--NUM_BASES_BELOW_QUALITY",
            "1",
        ]
    )

    log.info(
        f"{manifest.flowcell} - TagBamWithReadSequenceExtended (Cellular)"
        f" for {sample_row.library} in lane {lane}"
    )
    log.debug(f"Command = {' '.join(cmd)}")
    proc = run(cmd, capture_output=True)
    if proc.returncode != 0:
        log.error(f"Error running TagBamWithReadSequenceExtended:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info(
            f"{manifest.flowcell} - TagBamWithReadSequenceExtended (Cellular) completed"
        )

    # Tag bam with read sequence extended molecular
    cmd = dropseq_cmd("TagBamWithReadSequenceExtended", tmp_dir)
    cmd.extend(
        [
            "-I",
            f"{cellular_tagged_bam}",
            "-O",
            f"{molecular_tagged_bam}",
            "--SUMMARY",
            f"{molecular_tagged_summary}",
            "--BASE_RANGE",
            f"{bs_range2}",
            "--BASE_QUALITY",
            f"{sample_row.base_quality}",
            "--BARCODED_READ",
            "1",
            "--DISCARD_READ",
            "true",
            "--TAG_NAME",
            "XM",
            "--NUM_BASES_BELOW_QUALITY",
            "1",
        ]
    )
    log.info(
        f"{manifest.flowcell} - TagBamWithReadSequenceExtended (Molecular)"
        f" for {sample_row.library} in lane {lane}"
    )
    log.debug(f"Command = {' '.join(cmd)}")
    proc = run(cmd, capture_output=True)
    if proc.returncode != 0:
        log.error(f"Error running TagBamWithReadSequenceExtended:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info(
            f"{manifest.flowcell} - TagBamWithReadSequenceExtended (Molecular) completed"
        )
        os.remove(cellular_tagged_bam)

    # Filter low-quality reads
    cmd = dropseq_cmd("FilterBam", tmp_dir)
    cmd.extend(
        [
            "-I",
            f"{molecular_tagged_bam}",
            "-O",
            f"{filtered_ubam}",
            "--PASSING_READ_THRESHOLD",
            "0.1",
            "--REPAIR_BARCODES",
            "false",
            "--TAG_REJECT",
            "XQ",
        ]
    )
    log.info(f"{manifest.flowcell} - FilterBam for {sample_row.library} in Lane {lane}")
    log.debug(f"Command = {' '.join(cmd)}")
    proc = run(cmd, capture_output=True)
    if proc.returncode != 0:
        log.error(f"Error running FilterBam:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info(f"{manifest.flowcell} - FilterBam complete")
        os.remove(molecular_tagged_bam)

    # Trim reads with starting sequence
    cmd = dropseq_cmd("TrimStartingSequence", tmp_dir)
    cmd.extend(
        [
            "-I",
            f"{filtered_ubam}",
            "-O",
            f"{trimmed_ubam}",
            "--OUTPUT_SUMMARY",
            f"{trimming_summary}",
            "--SEQUENCE",
            f"{sample_row.start_sequence}",
            "--MISMATCHES",
            "0",
            "--NUM_BASES",
            "5",
        ]
    )
    log.info(
        f"{manifest.flowcell} - TrimStartingSequence for {sample_row.library} in Lane {lane}"
    )
    log.debug(f"Command = {' '.join(cmd)}")
    proc = run(cmd, capture_output=True)
    if proc.returncode != 0:
        log.error(f"Error running TrimStartingSequence:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info(f"{manifest.flowcell} - TrimStartingSequence complete")
        os.remove(filtered_ubam)

    # Adapter-aware poly A trimming
    cmd = dropseq_cmd("PolyATrimmer", tmp_dir)
    cmd.extend(
        [
            "-I",
            f"{trimmed_ubam}",
            "-O",
            f"{polya_filtered_ubam}",
            "--OUTPUT_SUMMARY",
            f"{polya_filtered_summary}",
            "--MISMATCHES",
            "0",
            "--NUM_BASES",
            "6",
            "--USE_NEW_TRIMMER",
            "true",
        ]
    )
    log.info(
        f"{manifest.flowcell} - PolyATrimmer for {sample_row.library} in Lane {lane}"
    )

    log.debug(f"Command = {' '.join(cmd)}")
    proc = run(cmd, capture_output=True)
    if proc.returncode != 0:
        log.error(f"Error running PolyATrimmer:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info(f"{manifest.flowcell} - PolyATrimmer complete")
        os.remove(trimmed_ubam)

    # Map reads to genome sequence using STARsolo
    cmd = [
        "STAR",
        "--soloType",
        "CB_samTagOut",
        "--genomeDir",
        f"{genome_dir}",
        "--readFilesIn",
        f"{polya_filtered_ubam}",
        "--readFilesCommand",
        "samtools",
        "view",
        "-F",
        "0x100",
        "--soloInputSAMattrBarcodeSeq",
        "XC",
        "XM",
        "--soloFeatures",
        "Gene",
        # "GeneFull",
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
    log.info(
        f"{manifest.flowcell} - Mapping using STAR for {sample_row.library} in Lane {lane}"
    )
    log.debug(f"Command = {' '.join(cmd)}")
    proc = run(cmd, capture_output=True)
    if proc.returncode != 0:
        log.error(f"Error running STAR:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info(f"{manifest.flowcell} - STAR alignment complete")
        os.remove(polya_filtered_ubam)

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
    log.info(f"{manifest.flowcell} - SortSam for {sample_row.library} in lane {lane}")
    log.debug(f"Command = {' '.join(cmd)}")
    proc = run(cmd, capture_output=True)
    if proc.returncode != 0:
        log.error(f"Error running SortSam:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info(f"{manifest.flowcell} - SortSam complete")
        os.remove(f"{aligned_bam}.star.Aligned.out.bam")

    # Merge unmapped bam and aligned bam
    cmd = picard_cmd("MergeBamAlignment", tmp_dir)
    cmd.extend(
        [
            "-R",
            f"{sample_row.reference}",
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

    log.info(
        f"{manifest.flowcell} - MergeBamAlignment for {sample_row.library} in lane {lane}"
    )
    log.debug(f"Command = {' '.join(cmd)}")
    proc = run(cmd, capture_output=True)
    if proc.returncode != 0:
        log.error(f"Error running MergeBamAlignment:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info(f"{manifest.flowcell} - MergeBamAlignment complete")
        os.remove(polya_filtered_ubam)
        os.remove(f"{aligned_bam}.aligned.sorted.bam")

    # Tag read with interval
    cmd = dropseq_cmd("TagReadWithInterval", tmp_dir)
    cmd.extend(
        [
            "-I",
            f"{aligned_bam}.merged.bam",
            "-O",
            f"{aligned_bam}.merged.TagReadWithInterval.bam",
            "--INTERVALS",
            f"{intervals}",
            "--TAG",
            "XG",
        ]
    )

    log.info(
        f"{manifest.flowcell} - TagReadWithInterval for {sample_row.library} in lane {lane}"
    )
    log.debug(f"Command = {' '.join(cmd)}")
    proc = run(cmd, capture_output=True)
    if proc.returncode != 0:
        log.error(f"Error running TagReadWithInterval:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info(f"{manifest.flowcell} - TagReadWithInterval complete")
        os.remove(f"{aligned_bam}.merged.bam")

    # Tag read with gene function
    cmd = dropseq_cmd("TagReadWithGeneFunction", tmp_dir, compression=5)
    cmd.extend(
        [
            "-I",
            f"{aligned_bam}.merged.TagReadWithInterval.bam",
            "-O",
            f"{aligned_bam}.star_gene_exon_tagged2.bam",
            "--ANNOTATIONS_FILE",
            f"{annotations_file}",
            "--CREATE_INDEX",
            "false",
        ]
    )

    log.info(
        f"{manifest.flowcell} - TagReadWithGeneFunction for {sample_row.library} in lane {lane}"
    )
    log.debug(f"Command = {' '.join(cmd)}")
    proc = run(cmd, capture_output=True)
    if proc.returncode != 0:
        log.error(f"Error running TagReadWithGeneFunction:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info(f"{manifest.flowcell} - TagReadWithGeneFunction complete")
        os.remove(f"{aligned_bam}.merged.TagReadWithInterval.bam")

    log.info(f"Alignment for {sample_row.library} completed")
