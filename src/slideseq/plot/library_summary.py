# This script is to generate PDFs for the alignment outputs

import csv
import gzip
import logging
from collections import Counter
from pathlib import Path

import matplotlib.figure
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

import slideseq.util.constants as constants
from slideseq.pipeline.metadata import Manifest

log = logging.getLogger(__name__)


def read_dropseq_metrics(metrics_file: Path, key_translate: dict[str, str] = None):
    """Reads the output from a dropseq analysis command, which is in a standard
    output. This works for ReadQualityMetrics.txt, fracExonicIntronic.txt and
    polyA_filtering.summary.txt

    :param metrics_file: File to read. Should be text with # as comment char,
                         metric keys on first non-comment row and values on second,
                         then a histogra,
    :param key_translate: Dictionary of strings to translate the keys for metrics
    """

    if key_translate is None:
        key_translate = dict()  # don't translate any names

    with metrics_file.open() as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and r[0][0] != "#"]

        # translate keys to nicer names
        metrics = {
            key_translate.get(k, k): int(v) for k, v in zip(rows[0][1:], rows[1][1:])
        }
        # histogram is of form `{bin: num reads}`
        histogram = Counter({int(r[0]): int(r[1]) for r in rows[3:]})

    return metrics, histogram


def read_quality_metrics(metrics_file: Path):
    """Reads the ReadQualityMetrics.txt file"""

    return read_dropseq_metrics(metrics_file, constants.READ_QUALITY_METRICS)


def read_frac_intronic_exonic(metrics_file: Path):
    """Reads the fracIntronicExonic.txt file"""

    return read_dropseq_metrics(metrics_file, constants.FRAC_INTRONIC_EXONIC)


def plot_mapping_quality(quality_metrics: Path, pdf_pages: PdfPages):
    # plot the overall mapping quality
    qm, _ = read_quality_metrics(quality_metrics)

    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    keys = ["Total", "Mapped", "HQ", "HQ No Dupes"]
    ax.bar(
        [f"{k}\n{qm[k]:,}\n({qm[k] / qm['Total']:.1%})" for k in keys],
        [qm[k] for k in keys],
        width=0.7,
        color="lightskyblue",
        edgecolor="black",
    )

    ax.set_ylabel("# Reads")
    ax.set_title("Alignment quality for all reads")
    pdf_pages.savefig(fig)

    return qm


def plot_frac_intronic_exonic(frac_intron_exon: Path, pdf_pages: PdfPages):
    # plot the fraction of reads mapping to different regions
    frac_ie, _ = read_frac_intronic_exonic(frac_intron_exon)

    # calculate some derivative measures
    frac_ie["exonic"] = frac_ie["coding_bases"] + frac_ie["utr_bases"]
    frac_ie["genic"] = frac_ie["exonic"] + frac_ie["intronic_bases"]

    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    keys = ["ribosomal", "exonic", "genic", "intronic", "intergenic"]

    ax.bar(
        [f"{k}\n{frac_ie[k]:.1%}" for k in keys],
        [100 * frac_ie[k] for k in keys],
        width=0.7,
        color="lightskyblue",
        edgecolor="black",
    )
    ax.set_ylim(0, 100)

    ax.set_ylabel("Percentage")
    ax.set_title("All reads")
    pdf_pages.savefig(fig)


def plot_reads_per_barcode(barcode_count, pdf_pages: PdfPages):
    with gzip.open(barcode_count, "rt") as fh:
        # this file has form `num_reads    barcodes` in descending order
        read_counts = [
            int(r[0]) for r in csv.reader(fh, delimiter="\t") if r[0][0] != "#"
        ]

    total_reads = sum(read_counts)
    # truncate to first 10% of barcodes
    read_counts = read_counts[: len(read_counts) // 10]

    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.plot(np.cumsum(read_counts) / total_reads, color="g")

    ax.set_ylim((0, 1))
    ax.set_xlabel("Top 10% of cell barcodes, by number of reads")
    ax.set_ylabel("Cumulative fraction of reads")
    ax.set_title("Cumulative fraction of reads per cell barcode")

    pdf_pages.savefig(fig)


def plot_poly_a_trimming(
    qm: dict[str, int], poly_a_summary_files: list[Path], pdf_pages: PdfPages
):
    total_hist = Counter()

    for summary_file in poly_a_summary_files:
        _, poly_a_hist = read_dropseq_metrics(summary_file)
        total_hist += poly_a_hist

    trimmed_count = sum(total_hist[k] for k in total_hist if k > 0)

    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.plot(
        np.arange(1, max(total_hist) + 1),
        [total_hist[k] for k in total_hist if k > 0],
        marker="o",
        color="k",
    )

    ax.set_xlabel("First base of PolyA tail trimmed")
    ax.set_ylabel("Number of reads")
    ax.set_title(
        f"% Reads trimmed by 3' PolyA trimmer: {trimmed_count / qm['Total']:.3%}"
    )
    ax.set_xlim(0, max(total_hist) + 1)
    pdf_pages.savefig(fig)


def plot_base_distribution(base_dist_file: Path, title: str, pdf_pages: PdfPages):
    with base_dist_file.open() as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
        barcode_distribution = np.array([[float(v) for v in r[1:-1]] for r in rows[1:]])
        barcode_distribution /= barcode_distribution.sum(axis=1, keepdims=True)

    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_prop_cycle("color", ["red", "blue", "green", "purple"])
    ax.plot(
        np.arange(1, barcode_distribution.shape[0] + 1),
        barcode_distribution,
        linewidth=0,
        marker="o",
        markersize=20,
        label=["A", "C", "G", "T"],
    )
    ax.legend(loc="lower right")
    ax.set_xlim(0, barcode_distribution.shape[0] + 2)
    ax.set_ylim(0, np.max(barcode_distribution) + 0.02)

    ax.set_xlabel("base position")
    ax.set_ylabel("fraction of reads")
    ax.set_title(title)
    pdf_pages.savefig(fig)


def make_library_summary(row: pd.Series, lanes: list[int], manifest: Manifest):
    """
    Creates a PDF with various library QA plots. This version is for slideseq runs
    that don't perform barcode matching

    :param row: contains metadata about the library
    :param lanes: the lanes that the library was sequenced across
    :param manifest: flowcell metadata
    """
    library_dir = constants.LIBRARY_DIR / f"{row.date}_{row.library}"
    library_base = library_dir / f"{row.library}.$"

    pdf_pages = PdfPages(library_base.with_suffix(".pdf"))

    qm = plot_mapping_quality(
        library_base.with_suffix(".ReadQualityMetrics.txt"), pdf_pages
    )

    plot_frac_intronic_exonic(
        library_base.with_suffix(".fracIntronicExonic.txt"), pdf_pages
    )

    plot_reads_per_barcode(
        library_base.with_suffix(f".numReads_perCell_XC_mq_{row.base_quality}.txt.gz"),
        pdf_pages,
    )

    poly_a_summaries = [
        library_dir
        / f"L{lane:03d}"
        / f"{manifest.flowcell}.L{lane:03d}.{row.library}.{row.barcode}.polyA_filtering.summary.txt"
        for lane in lanes
    ]

    # old pipeline put this in two places (I think) because he recalculated trimming after matching?
    plot_poly_a_trimming(qm, poly_a_summaries, pdf_pages)

    plot_base_distribution(
        library_base.with_suffix(".barcode_distribution_XC.txt"),
        "Cell barcodes for all reads",
        pdf_pages,
    )

    plot_base_distribution(
        library_base.with_suffix(".barcode_distribution_XM.txt"),
        "Cell barcodes for all reads",
        pdf_pages,
    )

    pdf_pages.close()
