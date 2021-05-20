import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from slideseq.util import constants as constants


def validate_library_df(library: str, library_df: pd.DataFrame):
    """Verify that all of the columns in the dataframe are constant, except
    for the lane which was expanded out earlier"""

    for col in constants.METADATA_COLS:
        if col.lower() == "lane":
            continue

        if len(set(library_df[col.lower()])) != 1:
            raise ValueError(f"Library {library} has multiple values in column {col}")


class Reference:
    """Represents a genome fasta and associated reference files"""

    def __init__(self, path):
        self.fasta = Path(path)

    @property
    def base(self) -> Path:
        return self.fasta.parent

    @property
    def genome_dir(self) -> Path:
        return self.base / "STAR"

    @property
    def annotations(self) -> Path:
        return self.fasta.with_suffix(".gtf")

    @property
    def intervals(self) -> Path:
        return self.fasta.with_suffix(".genes.intervals")

    @property
    def ribosomal_intervals(self) -> Path:
        return self.fasta.with_suffix(".rRNA.intervals")

    @property
    def ref_flat(self) -> Path:
        return self.fasta.with_suffix(".refFlat")


class Base:
    """Represents a BAM file and the associated metrics files"""

    def __init__(self, path: Path, base_quality: int, min_transcripts_per_cell: int):
        self.path = path
        self.base_quality = base_quality
        self.min_transcripts_per_cell = min_transcripts_per_cell

    def reads_per_cell(self, tag) -> Path:
        return self.path.with_suffix(
            f".numReads_perCell_{tag}_mq_{self.base_quality}.txt.gz"
        )

    def downsampled_bam(self, ratio) -> Path:
        return self.path.with_suffix(f".downsample_{ratio:.1f}.bam")

    @property
    def bam(self) -> Path:
        return self.path.with_suffix(".bam")

    @property
    def pdf(self) -> Path:
        return self.path.with_suffix(".pdf")

    @property
    def downsampling(self) -> Path:
        return self.path.with_suffix(".downsampling.pdf")

    @property
    def alignment_pdf(self) -> Path:
        return self.path.with_suffix(".alignment_quality.pdf")

    @property
    def mapping_rate(self) -> Path:
        return self.path.with_suffix(".mapping_rate.txt")

    @property
    def read_quality_metrics(self) -> Path:
        return self.path.with_suffix(".ReadQualityMetrics.txt")

    @property
    def digital_expression(self) -> Path:
        return self.path.with_suffix(".digital_expression.txt.gz")

    @property
    def digital_expression_summary(self) -> Path:
        return self.path.with_suffix(".digital_expression_summary.txt")

    @property
    def mtx(self) -> Path:
        return self.path.with_suffix(".digital_expression_matrix.mtx.gz")

    @property
    def barcodes(self) -> Path:
        return self.path.with_suffix(".digital_expression_barcodes.mtx.gz")

    @property
    def genes(self) -> Path:
        return self.path.with_suffix(".digital_expression_features.mtx.gz")

    @property
    def frac_intronic_exonic(self) -> Path:
        return self.path.with_suffix(".fracIntronicExonic.txt")

    @property
    def frac_intronic_exonic_per_cell(self) -> Path:
        return self.path.with_suffix(".fracIntronicExonicPerCell.txt.gz")

    @property
    def xc_distribution(self) -> Path:
        return self.path.with_suffix(".barcode_distribution_XC.txt")

    @property
    def xm_distribution(self) -> Path:
        return self.path.with_suffix(".barcode_distribution_XM.txt")

    @property
    def reads_per_umi(self) -> Path:
        return self.path.with_suffix(
            f".numReads_perUMI_XM_mq_{self.base_quality}.txt.gz"
        )

    @property
    def selected_cells(self) -> Path:
        return self.path.with_suffix(
            f".{self.min_transcripts_per_cell}_transcripts_mq_{self.base_quality}_selected_cells.txt.gz"
        )


@dataclass
class Library:
    # a dataclass to represent a single library that defines all the files
    # generated during processing
    name: str
    row: pd.Series
    lanes: list[int]
    reference: Reference

    def get_bead_structure(self) -> list[tuple[str, int, int]]:
        """
        Extract bead structure ranges from a format string. Takes an entry like:
            8C18X6C9M1X
        and returns:
            [('C', 1, 8), ('X', 9, 26), ('C', 27, 32), ('M', 33, 41), ('X', 42, 42)]
        """
        i = 1
        intervals = []

        # this is obfuscated to hell but I can't help myself
        #   - regex splits on C, X or M but retains the character found
        #   - 2 * (iter(re),) makes a second reference (not copy) to the iterator
        #   - zip(*(2iter)) zips together pairs of elements
        # then, add up the lengths to get ranges
        for j, c in zip(*(2 * (iter(re.split("([CXM])", self.row.bead_structure)),))):
            j = int(j)
            intervals.append((c, i, i + j - 1))
            i += j

        return intervals

    @property
    def base_quality(self) -> int:
        return self.row.base_quality

    @property
    def start_sequence(self) -> str:
        return self.row.start_sequence

    @property
    def min_transcripts_per_cell(self) -> int:
        return self.row.min_transcripts_per_cell

    @property
    def locus_function_list(self) -> str:
        return self.row.locus_function_list

    @property
    def puckcaller_path(self) -> Path:
        return Path(self.row.puckcaller_path)

    @property
    def bead_barcodes(self) -> Path:
        return self.puckcaller_path / "BeadBarcodes.txt"

    @property
    def bead_locations(self) -> Path:
        return self.puckcaller_path / "BeadLocations.txt"

    @property
    def run_barcodematching(self) -> bool:
        return self.row.run_barcodematching

    @property
    def gen_downsampling(self) -> bool:
        return self.row.gen_downsampling

    @property
    def dir(self) -> Path:
        """Base directory for the library data"""
        return constants.LIBRARY_DIR / f"{self.row.date}_{self.name}"

    @property
    def barcode_matching_dir(self) -> Path:
        return self.dir / "barcode_matching"

    @property
    def downsample_dir(self) -> Path:
        return self.dir / "downsample"

    @property
    def base(self) -> str:
        """base filename for the library. The '.$' is for with_suffix()"""
        return f"{self.name}.$"

    @property
    def pdf(self) -> Path:
        """PDF output"""
        return (self.dir / self.base).with_suffix(".pdf")

    @property
    def merged(self) -> Base:
        """Base for processed BAM merged across all lanes"""
        return Base(
            (self.dir / self.base).with_suffix(".all_illumina.$"),
            base_quality=self.base_quality,
            min_transcripts_per_cell=self.min_transcripts_per_cell,
        )

    @property
    def matched(self) -> Base:
        """Base for processed BAM filtered and retagged to matched barcodes"""
        return Base(
            (self.dir / self.base).with_suffix(".matched.$"),
            base_quality=self.base_quality,
            min_transcripts_per_cell=self.min_transcripts_per_cell,
        )

    def per_lane(self, *args, suffix) -> list[Path]:
        return [
            ((self.dir / f"L{lane:03d}").joinpath(*args) / self.base).with_suffix(
                suffix
            )
            for lane in self.lanes
        ]

    @property
    def polya_filtering_summaries(self) -> list[Path]:
        """polyA filtering summary"""
        return self.per_lane("alignment", suffix=".polyA_filtering.summary.txt")

    @property
    def star_logs(self) -> list[Path]:
        """Log files from STAR alignment"""
        return self.per_lane("alignment", suffix=".star.Log.final.out")

    @property
    def alignment_pickles(self) -> list[Path]:
        """Pickle files containing statistics about the alignment"""
        return self.per_lane("alignment", suffix=".alignment_statistics.pickle")

    @property
    def processed_bams(self) -> list[Path]:
        """Final aligned+unmapped output BAMs"""
        return self.per_lane(suffix=".final.bam")

    def __str__(self):
        return self.name


@dataclass
class LibraryLane(Library):
    lane: int

    @property
    def lane_dir(self) -> Path:
        return self.dir / f"L{self.lane:03d}"

    @property
    def raw_ubam(self) -> Path:
        """Raw unmapped BAM from demux. Keep this one for posterity"""
        return (self.lane_dir / self.base).with_suffix(".unmapped.bam")

    @property
    def processed_bam(self) -> Path:
        """Final merged, aligned output BAM"""
        return (self.lane_dir / self.base).with_suffix(".final.bam")

    @property
    def alignment_base(self) -> Path:
        return self.lane_dir / "alignment" / self.base

    @property
    def cellular_tagged_summary(self) -> Path:
        """summary file containing XC tag histograms"""
        return self.alignment_base.with_suffix(".cellular_tagging.summary.txt")

    @property
    def molecular_tagged_summary(self) -> Path:
        """summary file containing XM tag histograms"""
        return self.alignment_base.with_suffix(".molecular_tagging.summary.txt")

    @property
    def trimming_summary(self) -> Path:
        """summary file from adapter trimming"""
        return self.alignment_base.with_suffix(".adapter_trimming.summary.txt")

    @property
    def polya_filtering_summary(self) -> Path:
        """polyA filtering summary"""
        return self.alignment_base.with_suffix(".polyA_filtering.summary.txt")

    @property
    def polya_filtered_ubam(self) -> Path:
        """polyA filtered uBAM, used transiently"""
        return self.alignment_base.with_suffix(".polyA_filtered.bam")

    @property
    def polya_filtered_fastq(self) -> Path:
        """polyA filtered fastq.gz, input to STAR"""
        return self.alignment_base.with_suffix(".polyA_filtered.fastq.gz")

    @property
    def star_prefix(self) -> Path:
        """Prefix for STAR output"""
        return self.alignment_base.with_suffix(".star.")

    @property
    def aligned_bam(self) -> Path:
        """BAM file output from STAR alignment, used transiently"""
        return self.alignment_base.with_suffix(".star.Aligned.out.bam")

    @property
    def alignment_pickle(self) -> Path:
        """Pickle files containing statistics about the alignment"""
        return self.alignment_base.with_suffix(".alignment_statistics.pickle")
