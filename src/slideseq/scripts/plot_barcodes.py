import gzip
import itertools
import logging
from collections import Counter, defaultdict
from pathlib import Path

import click
import matplotlib
import networkx as nx
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.neighbors import radius_neighbors_graph

from slideseq.bead_matching import bipartite_matching, degen_barcode, hamming1_adjacency
from slideseq.plot import new_ax
from slideseq.util.logger import create_logger

log = logging.getLogger("plot_barcodes")


def plot_log_hist(dist, title, pdf_pages):
    max_d = np.ceil(np.log10(max(dist)))
    with new_ax(pdf_pages) as ax:
        ax.hist(dist, bins=np.logspace(0, max_d, int(max_d * 10 + 1)), log=True)
        ax.set_xscale("log")
        ax.set_title(title)


def spatial_plot(bead_xy_a, dist, title, pdf_pages, pct: float = 95.0):
    with new_ax(pdf_pages, include_fig=True) as (fig, ax):
        # version of 'Blues' colormap that is pure white at the bottom
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "BluesW", [(1.0, 1.0, 1.0), (0.0314, 0.188, 0.450)]
        )

        c = ax.scatter(
            bead_xy_a[:, 0],
            bead_xy_a[:, 1],
            c=dist,
            s=0.5,
            cmap=cmap,
            norm=matplotlib.colors.Normalize(0, np.percentile(dist, pct), clip=True),
        )
        c.set_rasterized(True)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.axis("equal")
        ax.set_title(title)
        fig.colorbar(c, ax=ax)


@click.command()
@click.argument("fastq_r1", type=click.Path(exists=True, path_type=Path))
@click.argument("fastq_r2", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--barcodes",
    type=click.Path(exists=True, path_type=Path),
    help="Path to BeadBarcodes.txt",
    required=True,
)
@click.option(
    "--locations",
    type=click.Path(exists=True, path_type=Path),
    help="Path to BeadLocations.txt",
    required=True,
)
@click.option(
    "--tag_sequence", default="TCGAGAGATCTACGGGTGGC", help="sequence to match against"
)
@click.option(
    "--output_pdf",
    type=click.Path(path_type=Path),
    help="Path to output figure",
    required=True,
)
@click.option(
    "--extra_pdf",
    type=click.Path(path_type=Path),
    help="Optional path for a bunch of other plots",
)
@click.option("--debug", is_flag=True, help="Turn on debug logging")
@click.option(
    "--percentile", type=float, default=95.0, help="Percentile for scaling plots"
)
def main(
    fastq_r1,
    fastq_r2,
    barcodes,
    locations,
    tag_sequence,
    output_pdf,
    extra_pdf=None,
    debug=False,
    percentile=95.0,
):
    """
    This script generates some plots for mapping barcoded reads.

    Reads sequences from FASTQ_R1 and FASTQ_R2. Assumes that the first read
    contains a 15bp barcode split across two locations, along with an 8bp UMI.
    The second read is assumed to have TAG_SEQUENCE in bases 20-40.
    """
    create_logger(debug, dryrun=False)

    log.debug(f"Reading from {fastq_r1}")
    with gzip.open(fastq_r1, "rt") as fh:
        r1_reads = [line.strip() for line in itertools.islice(fh, 1, None, 4)]

    log.debug(f"Reading from {fastq_r2}")
    with gzip.open(fastq_r2, "rt") as fh:
        r2_reads = [line.strip() for line in itertools.islice(fh, 1, None, 4)]

    log.debug(f"Reading {barcodes}")
    with open(barcodes) as fh:
        raw_bcs = ["".join(line.strip().split(",")) for line in fh]

    log.debug(f"Reading {locations}")
    with open(locations) as fh:
        x = np.array([float(v) for v in fh.readline().strip().split(",")])
        y = np.array([float(v) for v in fh.readline().strip().split(",")])
        xy = np.vstack((x, y)).T

    if extra_pdf is not None:
        extra_pdf_pages = PdfPages(extra_pdf)

        umi_counts = Counter(r[32:41] for r in r1_reads)
        log.debug(
            f"Found {len(umi_counts)} UMIs with {sum(umi_counts.values())} total counts"
        )

        plot_log_hist(umi_counts.values(), "Reads per UMI", extra_pdf_pages)
    else:
        extra_pdf_pages = None

    # pre-emptively remove poly-T/N sequences
    ok_barcodes = [not set(bc).issubset({"T", "N"}) for bc in raw_bcs]
    xy = xy[ok_barcodes, :]
    bead_barcodes = [bc for ok, bc in zip(ok_barcodes, raw_bcs) if ok]

    log.info(f"Read {len(raw_bcs)} barcodes and filtered to {len(bead_barcodes)}")

    seq_barcodes = sorted(r1[:8] + r1[26:32] for r1 in r1_reads)
    # remove poly-T sequence if present
    seq_barcodes = [seq for seq in seq_barcodes if set(seq) != {"T"}]

    log.info(f"Found {len(set(seq_barcodes))} unique barcodes in sequencing data")

    log.info("Computing barcode matching")

    log.debug("Computing radius neighbor graph")
    # adjacency matrix for all beads within radius of each other
    radius_matrix = radius_neighbors_graph(xy, radius=10.0)

    log.debug("Computing hamming neighbor graph")
    # adjacency matrix for all barcodes within hamming distance 1
    hamming_matrix = hamming1_adjacency(bead_barcodes)

    # just multiply together to get the combined adjacency matrix!
    combined_graph = nx.from_scipy_sparse_matrix(radius_matrix.multiply(hamming_matrix))

    # add xy coordinates to graph so we can analyze later
    for n, (x, y) in zip(combined_graph.nodes, xy):
        combined_graph.nodes[n]["x"] = x
        combined_graph.nodes[n]["y"] = y

    # get connected components to find groups of similar/close barcodes
    bead_groups = list(nx.connected_components(combined_graph))

    # calculate degenerate (ambiguous bases -> N) barcodes
    degen_bead_barcodes = [
        degen_barcode({bead_barcodes[j] for j in bg}) for bg in bead_groups
    ]

    log.debug(
        f"Collapsed {len(bead_groups)} bead groups into"
        f" {len(set(degen_bead_barcodes))} barcodes"
    )

    # average xy for grouped beads to get centroids
    bead_xy = dict()
    for bg, degen_bc in zip(bead_groups, degen_bead_barcodes):
        bg_graph = combined_graph.subgraph(bg)
        mean_x, mean_y = np.array(
            [[nd["x"], nd["y"]] for _, nd in bg_graph.nodes(data=True)]
        ).mean(0)
        bead_xy[degen_bc] = (mean_x, mean_y)

    barcode_matching = bipartite_matching(
        bead_barcodes, degen_bead_barcodes, bead_groups, seq_barcodes
    )

    if extra_pdf is not None:
        tag_barcodes = [r2[20:40] for r2 in r2_reads]
        tag_counts = Counter(tag_barcodes)

        sum(1 for r1 in r1_reads if (r1[:8] + r1[26:32]) in barcode_matching)

        umis_per_tag = defaultdict(set)
        for r1, r2 in zip(r1_reads, r2_reads):
            umis_per_tag[r2[20:40]].add(r1[32:41])

        plot_log_hist(tag_counts.values(), "Reads per tag", extra_pdf_pages)
        plot_log_hist(
            list(map(len, umis_per_tag.values())), "UMIs per tag", extra_pdf_pages
        )

    log.debug(f"Counting UMIs and reads per bead for sequence {tag_sequence}")
    reads_per_umi_per_bead = defaultdict(Counter)
    umis_per_bead = defaultdict(set)
    reads_per_bead = Counter()

    for r1, r2 in zip(r1_reads, r2_reads):
        seq_bc = r1[:8] + r1[26:32]

        if seq_bc not in barcode_matching:
            continue
        if r2[20:40] != tag_sequence:
            continue

        bead_bc = barcode_matching[seq_bc]
        umi = r1[32:41]

        reads_per_umi_per_bead[bead_bc][umi] += 1
        umis_per_bead[bead_bc].add(umi)
        reads_per_bead[bead_bc] += 1

    filtered_barcodes = [bc for bc in degen_bead_barcodes if umis_per_bead[bc]]
    bead_xy_a = np.vstack([bead_xy[dbc] for dbc in filtered_barcodes])

    with gzip.open(output_pdf.with_suffix(".reads_per_umi.txt.gz"), "wt") as out:
        print("bead_barcodes\tumi\treads", file=out)
        for bc in filtered_barcodes:
            for umi in reads_per_umi_per_bead[bc][umi]:
                print(f"{bc}\t{umi}\t{reads_per_umi_per_bead[bc][umi]}", file=out)

    with output_pdf.with_suffix(".txt").open("w") as out:
        print("bead_barcode\tumis\treads", file=out)
        for bc in filtered_barcodes:
            print(f"{bc}\t{len(umis_per_bead[bc])}\t{reads_per_bead[bc]}", file=out)

    if extra_pdf is not None:
        plot_log_hist(
            [len(umis_per_bead[bc]) for bc in filtered_barcodes],
            "UMIs per bead",
            extra_pdf_pages,
        )

        extra_pdf_pages.close()

    pdf_pages = PdfPages(output_pdf)

    log.info("Making plots")
    spatial_plot(
        bead_xy_a,
        [len(umis_per_bead[bc]) for bc in filtered_barcodes],
        "UMIs per bead",
        pdf_pages,
        pct=percentile,
    )

    spatial_plot(
        bead_xy_a,
        [np.log10(len(umis_per_bead[bc])) for bc in filtered_barcodes],
        "log10 UMIs per bead",
        pdf_pages,
        pct=percentile,
    )

    spatial_plot(
        bead_xy_a,
        [reads_per_bead[bc] for bc in filtered_barcodes],
        "Reads per bead",
        pdf_pages,
        pct=percentile,
    )

    spatial_plot(
        bead_xy_a,
        [np.log10(1 + reads_per_bead[bc]) for bc in filtered_barcodes],
        "log10 reads per bead",
        pdf_pages,
        pct=percentile,
    )

    pdf_pages.close()
    log.info("Done!")


if __name__ == "__main__":
    main()
