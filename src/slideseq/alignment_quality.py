import logging
import pickle
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.figure
import pysam
from matplotlib.backends.backend_pdf import PdfPages

from slideseq.library import Library

log = logging.getLogger(__name__)


def write_alignment_stats(bam_file: Path, out_file: Path):
    """
    Read through a BAM file and collect some distributions:
      - alignment scores for uniquely-mapped reads
      - number of mismatches per uniquely-mapped read
      - ratio of matches / total for uniquely-mapped reads
      - the above three stats but for multi-mapped reads

    And write them into a pickle

    :param bam_file: aligned BAM file to check (output from STAR)
    :param out_file: file to output distributions (as a pickle)
    """
    aligned_bam = pysam.AlignmentFile(bam_file, mode="rb", threads=8)

    mp = Counter()
    for a in aligned_bam:
        mp[a.query_name] += 1

    mp = {qn for qn, count in mp.items() if count > 1}

    aligned_bam.reset()

    unique_score = Counter()
    unique_mismatch = Counter()
    unique_ratio = Counter()
    multi_score = Counter()
    multi_mismatch = Counter()
    multi_ratio = Counter()

    for a in aligned_bam:
        if a.query_name in mp:
            if a.has_tag("AS"):
                multi_score[a.get_tag("AS")] += 1
            if a.has_tag("nM"):
                multi_mismatch[a.get_tag("nM")] += 1
            r = round(100 * a.get_cigar_stats()[0][0] / a.query_length)
            multi_ratio[r] += 1
        else:
            if a.has_tag("AS"):
                unique_score[a.get_tag("AS")] += 1
            if a.has_tag("nM"):
                unique_mismatch[a.get_tag("nM")] += 1
            r = round(100 * a.get_cigar_stats()[0][0] / a.query_length)
            unique_ratio[r] += 1

    log.debug(f"Writing alignment stats for {bam_file} to {out_file}")

    with out_file.open("wb") as out:
        pickle.dump(
            (
                unique_score,
                unique_mismatch,
                unique_ratio,
                multi_score,
                multi_mismatch,
                multi_ratio,
            ),
            out,
        )


def plot_hist(kind: str, metric: str, dist: Counter, pdf_pages: PdfPages):
    """
    Bar plot of a given alignment stat distribution and adds the plot to
    a PDF

    :param kind: Either "unique" or "multi"
    :param metric: Either "ratio", "score", or "mismatch"
    :param dist: The distribution of counts
    :param pdf_pages: PDF to add the figure to
    """

    x, y = zip(*sorted(dist.items()))
    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.bar(x, y, width=0.7, color="lightskyblue", edgecolor="black")
    if metric == "ratio":
        ax.set_xlabel(f"alignment {metric} %")
    else:
        ax.set_xlabel(f"alignment {metric}")
    ax.set_title(f"Histogram of {kind} alignment {metric}")

    pdf_pages.savefig(fig)


def parse_star_log(log_file: Path) -> dict[str, str]:
    """
    Read the Log.final.out file from STAR and return the parsed data

    :param log_file: Path to STAR's final output log
    :return: Dictionary containing the run stats
    """
    with log_file.open() as fh:
        rows = [line.strip().split("|\t") for line in fh if "|" in line]
        log_data = {r[0].strip(): r[1].strip() for r in rows}

    return log_data


def combine_alignment_stats(library: Library):
    """
    Combines the alignment statistics files from multiple lanes into one aggregate,
    and outputs aggregated versions into text files. Also aggregates the STAR log
    reports into one combined file.

    :param library: object that contains metadata about this library
    """
    dists = defaultdict(Counter)

    log.debug(
        "Combining"
        f" {', '.join(stat_file.name for stat_file in library.alignment_pickles)}"
    )

    for stat_file in library.alignment_pickles:
        with stat_file.open("rb") as fh:
            us, um, ur, ms, mm, mr = pickle.load(fh)
            dists["unique", "ratio"] += ur
            dists["multi", "ratio"] += mr
            dists["unique", "score"] += us
            dists["multi", "score"] += ms
            dists["unique", "mismatch"] += um
            dists["multi", "mismatch"] += mm

    log.debug(f"Plotting alignment stats for {library.merged.bam}")
    alignment_quality_pdf = PdfPages(library.merged.alignment_pdf)

    # make PDF with combined stats
    for metric in ("ratio", "score", "mismatch"):
        for kind in ("unique", "multi"):
            plot_hist(kind, metric, dists[kind, metric], alignment_quality_pdf)

    alignment_quality_pdf.close()

    total_reads = 0
    unique_reads = 0
    multi_reads = 0
    too_many_reads = 0

    # collect some stats from star log files
    for star_log in library.star_logs:
        log_data = parse_star_log(star_log)
        total_reads += int(log_data["Number of input reads"])
        unique_reads += int(log_data["Uniquely mapped reads number"])
        multi_reads += int(log_data["Number of reads mapped to multiple loci"])
        too_many_reads += int(log_data["Number of reads mapped to too many loci"])

    mismatch = dists["unique", "mismatch"] + dists["multi", "mismatch"]

    with library.merged.mapping_rate.open("w") as out:
        print(f"total_reads\t{total_reads}", file=out)
        print(f"unique_aligned_reads\t{unique_reads}", file=out)
        print(f"unique_aligned_ratio\t{unique_reads / total_reads:.3%}", file=out)
        print(f"multi_aligned_reads\t{multi_reads}", file=out)
        print(f"multi_aligned_ratio\t{multi_reads / total_reads:.3%}", file=out)
        print(f"too_many_aligned_reads\t{too_many_reads}", file=out)
        print(
            f"too_many_aligned_ratio\t{too_many_reads / total_reads:.3%}",
            file=out,
        )
        print(f"mismatch1_rate\t{mismatch[1] / total_reads:.3%}", file=out)
        print(f"mismatch2_rate\t{mismatch[2] / total_reads:.3%}", file=out)
        print(f"mismatch3_rate\t{mismatch[3] / total_reads:.3%}", file=out)
