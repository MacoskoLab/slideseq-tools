# This script is to generate PDFs for the alignment outputs

import csv
import gzip
import logging
from collections import Counter
from pathlib import Path

import matplotlib.colors
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from slideseq.library import Library
from slideseq.plot import (
    new_ax,
    read_dge_summary,
    read_dropseq_metrics,
    read_frac_intronic_exonic,
    read_quality_metrics,
)

log = logging.getLogger(__name__)

BeadXY = dict[str, tuple[float, float]]


def plot_combined_mapping_quality(
    pdf_pages: PdfPages, quality_metrics: Path, matched_quality_metrics: Path
):
    keys = ["Total", "Mapped", "HQ", "HQ No Dupes"]
    qm, _ = read_quality_metrics(quality_metrics)
    mm, _ = read_quality_metrics(matched_quality_metrics)
    ms = [qm, mm]

    x = np.arange(len(keys))

    with new_ax(pdf_pages) as ax:
        ax.set_prop_cycle("color", ["lightskyblue", "goldenrod"])
        for i, m in enumerate(ms):
            ax.bar(
                x + 0.1 + 0.4 * i,
                [m[k] for k in keys],
                width=0.4,
                align="edge",
                label=["All", "Matched"][i],
                edgecolor="black",
            )

        ax.set_xticks(x + 0.5)
        ax.set_xticklabels(keys)

        ax.set_xticks(
            [v + 0.3 + i for v in x for i in (0.0, 0.4)],
            minor=True,
        )
        ax.set_xticklabels(
            [f"{int(m[k]):,}\n({m[k] / m['Total']:.1%})" for k in keys for m in ms],
            minor=True,
            rotation=90,
        )

        ax.set_ylabel("# Reads")
        ax.set_title("Alignment Quality")

        ax.tick_params(axis="x", which="major", bottom=False, length=80)
        ax.tick_params(axis="x", which="minor")

        ax.legend()

    return qm


def plot_mapping_quality(pdf_pages: PdfPages, quality_metrics: Path):
    keys = ["Total", "Mapped", "HQ", "HQ No Dupes"]
    qm, _ = read_quality_metrics(quality_metrics)

    with new_ax(pdf_pages) as ax:
        ax.bar(
            [f"{k}\n{int(qm[k]):,}\n({qm[k] / qm['Total']:.1%})" for k in keys],
            [qm[k] for k in keys],
            width=0.8,
            color="lightskyblue",
            edgecolor="black",
        )

        ax.set_ylabel("# Reads")
        ax.set_title("Alignment Quality for All Reads")

    return qm


def plot_frac_intronic_exonic(pdf_pages: PdfPages, frac_intron_exon: Path):
    # plot the fraction of reads mapping to different regions
    frac_ie, _ = read_frac_intronic_exonic(frac_intron_exon)

    # calculate some derivative measures
    frac_ie["exonic"] = frac_ie["coding"] + frac_ie["utr"]
    frac_ie["genic"] = frac_ie["exonic"] + frac_ie["intronic"]
    keys = ["ribosomal", "exonic", "genic", "intronic", "intergenic"]

    with new_ax(pdf_pages) as ax:
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


def plot_reads_per_barcode(pdf_pages: PdfPages, barcode_count):
    with gzip.open(barcode_count, "rt") as fh:
        # this file has form `num_reads    barcodes` in descending order
        read_counts = [
            int(r[0]) for r in csv.reader(fh, delimiter="\t") if r[0][0] != "#"
        ]

    total_reads = sum(read_counts)
    # truncate to first 10% of barcodes
    read_counts = read_counts[: len(read_counts) // 10]

    with new_ax(pdf_pages) as ax:
        ax.plot(np.cumsum(read_counts) / total_reads, color="g")

        ax.set_ylim((0, 1))
        ax.set_xlabel("Top 10% of cell barcodes, by number of reads")
        ax.set_ylabel("Cumulative fraction of reads")
        ax.set_title("Cumulative fraction of reads per cell barcode")


def plot_poly_a_trimming(
    pdf_pages: PdfPages, qm: dict[str, int], poly_a_summary_files: list[Path]
):
    total_hist = Counter()

    for summary_file in poly_a_summary_files:
        _, poly_a_hist = read_dropseq_metrics(summary_file)
        total_hist += poly_a_hist

    trimmed_count = sum(total_hist[k] for k in total_hist if k > 0)

    with new_ax(pdf_pages) as ax:
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


def plot_base_distribution(pdf_pages: PdfPages, base_dist_file: Path, title: str):
    with base_dist_file.open() as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
        base_distribution = np.array([[float(v) for v in r[1:-1]] for r in rows[1:]])
        base_distribution /= base_distribution.sum(axis=1, keepdims=True)

    with new_ax(pdf_pages) as ax:
        ax.set_prop_cycle("color", ["red", "blue", "green", "purple"])
        ax.plot(
            np.arange(1, base_distribution.shape[0] + 1),
            base_distribution,
            linewidth=0,
            marker="o",
            markersize=10,
            alpha=0.8,
            label=["A", "C", "G", "T"],
        )
        ax.legend(loc="lower right")
        ax.set_xlim(0, base_distribution.shape[0] + 2)
        ax.set_ylim(0, np.max(base_distribution) + 0.02)

        ax.set_xlabel("base position")
        ax.set_ylabel("fraction of reads")
        ax.set_title(title)


def plot_spatial_distribution(
    pdf_pages: PdfPages, bead_xy: np.ndarray, dist: list[float], title: str
):
    with new_ax(pdf_pages, include_fig=True) as (fig, ax):
        c = ax.scatter(
            bead_xy[:, 0],
            bead_xy[:, 1],
            c=dist,
            s=0.5,
            cmap="viridis_r",
            norm=matplotlib.colors.Normalize(0, np.percentile(dist, 95), clip=True),
        )
        c.set_rasterized(True)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.axis("equal")
        ax.set_title(f"{title} per matched bead")
        fig.colorbar(c, ax=ax)


def plot_dge_summary(pdf_pages: PdfPages, summary_file: Path, bead_xy: BeadXY):
    """
    Plot histograms of UMIs and genes per barcode
    """

    barcodes, umis_per_bc, genes_per_bc = read_dge_summary(summary_file)

    bead_xy = np.vstack([bead_xy[bc] for bc in barcodes])

    for dist, title in ((umis_per_bc, "UMIs"), (genes_per_bc, "genes")):
        with new_ax(pdf_pages) as ax:
            # can safely assume no zeros in this distribution
            # 10000 should be good enough for anyone
            ax.hist(
                dist,
                bins=np.logspace(0, 4, 41),
                facecolor="lightskyblue",
                edgecolor="black",
            )
            ax.set_xscale("log")
            ax.set_xlabel(f"Number of {title} (log10)")
            ax.set_title(f"Histogram of {title} per matched barcode")

    plot_spatial_distribution(pdf_pages, bead_xy, umis_per_bc, "UMIs")


def plot_scrna_metrics(pdf_pages: PdfPages, metrics_file: Path, bead_xy: BeadXY):
    with gzip.open(metrics_file, "rt") as fh:
        rows = list(
            csv.DictReader(
                (line for line in fh if len(line) > 1 and line[0] != "#"),
                delimiter="\t",
            )
        )

    barcodes = [r["LIBRARY"] for r in rows]
    bead_xy = np.vstack([bead_xy[bc] for bc in barcodes])

    for k, title in (
        ("PCT_MT_BASES", "% mitochondrial"),
        ("PCT_CODING_BASES", "% exonic"),
        ("PCT_RIBOSOMAL_BASES", "% ribosomal"),
    ):
        dist = [float(r[k]) for r in rows]
        plot_spatial_distribution(pdf_pages, bead_xy, dist, title)


def make_library_plots(library: Library, bead_xy: BeadXY = None):
    """
    Creates a PDF with various library QA plots. This version is for slideseq runs
    that don't perform barcode matching

    :param library: contains metadata about the library
    :param bead_xy: xy coordinates for matched beads, if matching was performed
    """

    log.info(f"Creating report file {library.pdf}")
    pdf_pages = PdfPages(library.pdf)

    if library.matched.read_quality_metrics.exists():
        qm = plot_combined_mapping_quality(
            pdf_pages,
            library.merged.read_quality_metrics,
            library.matched.read_quality_metrics,
        )
    else:
        qm = plot_mapping_quality(pdf_pages, library.merged.read_quality_metrics)

    plot_frac_intronic_exonic(pdf_pages, library.merged.frac_intronic_exonic)

    plot_reads_per_barcode(pdf_pages, library.merged.reads_per_cell("XC"))

    plot_poly_a_trimming(pdf_pages, qm, library.polya_filtering_summaries)

    for file_path, title in (
        (library.merged.xc_distribution, "Cell"),
        (library.merged.xm_distribution, "Molecular"),
    ):
        plot_base_distribution(pdf_pages, file_path, f"{title} barcodes for all reads")

    if bead_xy is not None:
        assert library.matched.bam.exists(), f"{library.matched.bam} does not exist"

        plot_frac_intronic_exonic(pdf_pages, library.matched.frac_intronic_exonic)

        plot_dge_summary(pdf_pages, library.matched.digital_expression_summary, bead_xy)

        plot_scrna_metrics(pdf_pages, library.matched.frac_intronic_exonic, bead_xy)

        for file_path, title in (
            (library.matched.xc_distribution, "Cell"),
            (library.matched.xm_distribution, "Molecular"),
        ):
            plot_base_distribution(
                pdf_pages, file_path, f"{title} barcodes for matched reads"
            )

    pdf_pages.close()
