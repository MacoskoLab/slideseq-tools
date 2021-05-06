import logging
from pathlib import Path

import pysam

log = logging.getLogger(__name__)


def write_retagged_bam(
    bam_file: Path, retagged_bam_file: Path, barcode_matching: dict[str, str]
):
    """
    Takes a mapped BAM file and translates the original sequence cell barcodes
    into mapped barcodes computed by the bead_matching procedure, filtering out
    all others.

    :param bam_file: Input bam file, aligned and tagged with barcodes as XC
    :param retagged_bam_file: Output bam file with translated XC tags
    :param barcode_matching: a dictionary of old to new barcodes
    """

    log.debug(f"Reading {bam_file}")
    log.debug(f"Writing to {retagged_bam_file}")
    with pysam.AlignmentFile(bam_file, mode="rb") as mapped_bam:
        with pysam.AlignmentFile(
            retagged_bam_file, mode="wb", template=mapped_bam
        ) as tagged_bam:
            for a in mapped_bam:
                if a.get_tag("XC") in barcode_matching:
                    a.set_tag("XC", barcode_matching[a.get_tag("XC")])
                    tagged_bam.write(a)

    log.debug(f"Finished writing to {retagged_bam_file}")
