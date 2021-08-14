import gzip
import itertools
import logging
from collections import Counter, defaultdict
from pathlib import Path

import click
import matplotlib
import matplotlib.colors
import networkx as nx
import numpy as np
import scipy.io
import scipy.sparse
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.neighbors import radius_neighbors_graph

import slideseq.bead_matching
from slideseq.plot import new_ax
from slideseq.util.logger import create_logger

log = logging.getLogger("barcode_matrix")


DEGENERATE_BASE_DICT = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "M": {"A", "C"},
    "K": {"G", "T"},
    "S": {"C", "G"},
    "W": {"A", "T"},
    "H": {"A", "C", "T"},
    "B": {"C", "G", "T"},
    "V": {"A", "C", "G"},
    "D": {"A", "G", "T"},
    "N": {"A", "C", "G", "T"},
}


def plot_log_hist(dist, title: str, pdf_pages: PdfPages):
    max_d = np.ceil(np.log10(max(dist)))
    with new_ax(pdf_pages) as ax:
        ax.hist(dist, bins=np.logspace(0, max_d, int(max_d * 10 + 1)), log=True)
        ax.set_xscale("log")
        ax.set_title(title)


def plot_hist(dist, title: str, pdf_pages: PdfPages):
    with new_ax(pdf_pages) as ax:
        ax.hist(dist, bins=100, log=True)
        ax.set_title(title)


def spatial_plot(bead_xy_a, dist, title: str, pdf_pages: PdfPages, pct: float = 95.0):
    with new_ax(pdf_pages, include_fig=True) as (fig, ax):
        # version of 'Blues' colormap that is pure white at the bottom
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "BluesW", [(1.0, 1.0, 1.0), (0.0314, 0.188, 0.450)]
        )
        norm = matplotlib.colors.Normalize(0, np.percentile(dist, pct), clip=True)

        c = ax.scatter(
            bead_xy_a[:, 0],
            bead_xy_a[:, 1],
            c=dist,
            s=0.5,
            cmap=cmap,
            norm=norm,
        )
        c.set_rasterized(True)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.axis("equal")
        ax.set_title(title)
        fig.colorbar(c, ax=ax)


def spatial_plot_plus_links(
    bead_xy_a, bead_pairs, umi_dist, title: str, pdf_pages: PdfPages, pct: float = 95
):
    with new_ax(pdf_pages, include_fig=True) as (fig, ax):
        # version of 'Blues' colormap that is pure white at the bottom
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "BluesW", [(1.0, 1.0, 1.0), (0.0314, 0.188, 0.450)]
        )
        norm = matplotlib.colors.Normalize(0, np.percentile(umi_dist, pct), clip=True)

        c = ax.scatter(
            bead_xy_a[:, 0],
            bead_xy_a[:, 1],
            c=umi_dist,
            s=0.5,
            cmap=cmap,
            norm=norm,
        )
        c.set_rasterized(True)

        xs = bead_xy_a[bead_pairs, 0].T
        ys = bead_xy_a[bead_pairs, 1].T

        ax.plot(xs, ys, alpha=0.1, color="g", linewidth=0.1)

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.axis("equal")
        ax.set_title(title)
        fig.colorbar(c, ax=ax)


def get_reads(seq_dir: Path, sample_n: str):
    fastq_r1 = next(seq_dir.glob(f"Sample{sample_n}_S?_L001_R1_001.fastq.gz"))
    fastq_r2 = next(seq_dir.glob(f"Sample{sample_n}_S?_L001_R2_001.fastq.gz"))

    with gzip.open(fastq_r1, "rt") as fh:
        r1_reads = [line.strip() for line in itertools.islice(fh, 1, None, 4)]

    with gzip.open(fastq_r2, "rt") as fh:
        r2_reads = [line.strip() for line in itertools.islice(fh, 1, None, 4)]

    return r1_reads, r2_reads


def get_barcodes(puck_dir: Path):
    log.info(f"Reading bead data from {puck_dir}")
    log.debug(f"Reading {puck_dir / 'BeadBarcodes.txt'}")
    with open(puck_dir / "BeadBarcodes.txt") as fh:
        raw_bcs = ["".join(line.strip().split(",")) for line in fh]

    log.debug(f"Reading {puck_dir / 'BeadLocations.txt'}")
    with open(puck_dir / "BeadLocations.txt") as fh:
        x = np.array([float(v) for v in fh.readline().strip().split(",")])
        y = np.array([float(v) for v in fh.readline().strip().split(",")])
        xy = np.vstack((x, y)).T

    ok_barcodes = [not set(bc).issubset({"T", "N"}) for bc in raw_bcs]
    xy = xy[ok_barcodes, :]
    bead_barcodes = [bc for ok, bc in zip(ok_barcodes, raw_bcs) if ok]

    log.info(f"Read {len(raw_bcs)} barcodes and filtered to {len(bead_barcodes)}")

    return bead_barcodes, xy


def bead_network(bead_barcodes, seq_barcodes, xy):
    # adjacency matrix for all beads within radius of each other
    radius_matrix = radius_neighbors_graph(xy, radius=10.0)

    # adjacency matrix for all barcodes within hamming distance 1
    hamming_matrix = slideseq.bead_matching.hamming1_adjacency(bead_barcodes)

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
        slideseq.bead_matching.degen_barcode({bead_barcodes[j] for j in bg})
        for bg in bead_groups
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

    barcode_matching = slideseq.bead_matching.bipartite_matching(
        bead_barcodes, degen_bead_barcodes, bead_groups, seq_barcodes
    )

    return degen_bead_barcodes, bead_xy, barcode_matching


def tag_network(r2_reads, barcode_codes, cutoff=30):
    tag_sequences = [
        r2[25:57]
        for r2 in r2_reads
        if sum(b in s for b, s in zip(r2[25:57], barcode_codes)) > cutoff
    ]
    tag_sequences = sorted(set(tag_sequences))

    tag_graph = nx.from_scipy_sparse_matrix(
        slideseq.bead_matching.hamming1_adjacency(tag_sequences)
    )

    tag_groups = list(nx.connected_components(tag_graph))
    log.debug(f"{len(tag_groups)} tag components")

    degen_tag_barcodes = [
        slideseq.bead_matching.degen_barcode({tag_sequences[j] for j in tg})
        for tg in tag_groups
    ]

    degen_tag_counter = Counter(degen_tag_barcodes)

    tag_to_degen_tag = {
        tag_sequences[j]: dtag
        for tg, dtag in zip(tag_groups, degen_tag_barcodes)
        if degen_tag_counter[dtag] == 1
        for j in tg
    }

    degen_tag_barcodes = [
        dtag for dtag in degen_tag_barcodes if degen_tag_counter[dtag] == 1
    ]
    log.debug(f"{len(degen_tag_barcodes)} unambiguous tags")

    return tag_sequences, tag_to_degen_tag


def match_tags(
    r1_reads, r2_reads, barcode_matching, constant_seq_hd1, tag_to_degen_tag
):
    matched_but_no_const = 0
    match_and_const_but_no_tag = 0

    umis_per_bead = defaultdict(lambda: defaultdict(set))
    reads_per_bead = defaultdict(Counter)

    for r1, r2 in zip(r1_reads, r2_reads):
        seq_bc = r1[:8] + r1[26:32]

        if seq_bc not in barcode_matching:
            continue

        if r2[57:85] not in constant_seq_hd1:
            matched_but_no_const += 1
            continue

        tag = r2[25:57]

        if tag not in tag_to_degen_tag:
            match_and_const_but_no_tag += 1
            continue

        bead_bc = barcode_matching[seq_bc]
        umi = r1[32:41]
        dtag = tag_to_degen_tag[tag]

        umis_per_bead[bead_bc][dtag].add(umi)
        reads_per_bead[bead_bc][dtag] += 1

    log.debug(f"bead matched, no constant seq: {matched_but_no_const}")
    log.debug(f"matched w/ const seq, no tag: {match_and_const_but_no_tag}")

    umis_per_bead = {
        bead_bc: {dtag: len(v) for dtag, v in dtags_per_bead.items()}
        for bead_bc, dtags_per_bead in umis_per_bead.items()
    }

    return umis_per_bead, reads_per_bead


def write_matrix(upb, beads, tags, output_file):
    m = scipy.sparse.dok_matrix((len(beads), len(tags)), dtype=np.int32)
    b2i = {b: i for i, b in enumerate(beads)}
    t2j = {t: j for j, t in enumerate(tags)}

    for b, bd in upb.items():
        for t, v in bd.items():
            m[b2i[b], t2j[t]] = v

    m = m.tocsr()

    with gzip.open(output_file, "wb") as out:
        scipy.io.mmwrite(out, m)

    return m


@click.command()
@click.argument("fastq-r1", type=click.Path(exists=True))
@click.argument("fastq-r2", type=click.Path(exists=True))
@click.option(
    "--puck-dir",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    help="Path containing BeadBarcodes.txt and BeadLocations.txt",
    required=True,
)
@click.option(
    "--output-dir",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    help="Path for output",
    required=True,
)
@click.option(
    "--tag-sequence",
    default="WSBWSDWSHWSVWSBWSWSHWSVWSBWSDWSH",
    help="Format for the tag sequence",
)
@click.option(
    "--constant-sequence",
    default="TCGAGAGATCTACGGGTGGCATCCCTGT",
    help="Constant sequence that should follow the tag",
)
@click.option("--extra-pdf", is_flag=True, help="Generate extra plots")
@click.option("--debug", is_flag=True, help="Turn on debug logging")
@click.option(
    "--percentile",
    type=float,
    default=95.0,
    help="Percentile for scaling plots",
    show_default=True,
)
@click.option(
    "--tag-mismatch",
    type=int,
    default=1,
    help="Number of mismatches to allow in tag",
    show_default=True,
)
@click.option(
    "--min-dist",
    type=int,
    default=1000,
    help="Minimum distance (in pixels) for plotting connections",
    show_default=True,
)
def main(
    fastq_r1,
    fastq_r2,
    puck_dir,
    output_dir,
    tag_sequence,
    constant_sequence,
    extra_pdf=False,
    debug=False,
    percentile=95.0,
    tag_mismatch=1,
    min_dist=1000,
):
    """
    This script generates some plots for mapping barcoded reads.

    Reads sequences from FASTQ_R1 and FASTQ_R2. Assumes that the first read
    contains a 15bp barcode split across two locations, along with an 8bp UMI.
    The second read is assumed to have TAG_SEQUENCE in bases 20-40.
    """
    create_logger(debug, dryrun=False)

    output_dir = Path(output_dir)
    output_pdf = output_dir / "plots.pdf"

    log.info(f"Saving output to {output_dir}")

    log.debug(f"Reading from {fastq_r1}")
    with gzip.open(fastq_r1, "rt") as fh:
        r1_reads = [line.strip() for line in itertools.islice(fh, 1, None, 4)]

    log.debug(f"Reading from {fastq_r2}")
    with gzip.open(fastq_r2, "rt") as fh:
        r2_reads = [line.strip() for line in itertools.islice(fh, 1, None, 4)]

    assert len(r1_reads) == len(r2_reads), "read different number of reads"
    log.info(f"Total of {len(r1_reads)} reads")

    bead_barcodes, xy = get_barcodes(Path(puck_dir))

    constant_sequence_hset = slideseq.bead_matching.hamming_set(
        slideseq.bead_matching.initial_h_set(constant_sequence), include_N=False
    )

    barcode_codes = [DEGENERATE_BASE_DICT[b] for b in tag_sequence]

    # get unique barcodes
    seq_barcodes = sorted({r1[:8] + r1[26:32] for r1 in r1_reads})
    # remove poly-T sequence if present
    seq_barcodes = [seq for seq in seq_barcodes if set(seq) != {"T"}]

    log.info(f"{len(seq_barcodes)} unique seq barcodes")

    log.debug("calculating bead network")
    degen_bead_barcodes, bead_xy, barcode_mapping = bead_network(
        bead_barcodes, seq_barcodes, xy
    )

    log.debug("calculating tag network")
    tag_sequences, tag_to_degen_tag = tag_network(
        r2_reads, barcode_codes, len(tag_sequence) - tag_mismatch
    )

    with gzip.open(output_dir / "tag_mapping.txt.gz", "wt") as out:
        for tag, dtag in tag_to_degen_tag.items():
            print(f"{tag}\t{dtag}", file=out)

    with gzip.open(output_dir / "raw_sequences.txt.gz", "wt") as out:
        print("\n".join(tag_sequences), file=out)

    umis_per_bead, reads_per_bead = match_tags(
        r1_reads, r2_reads, barcode_mapping, constant_sequence_hset, tag_to_degen_tag
    )

    slideseq.bead_matching.write_barcode_mapping(
        barcode_mapping, bead_xy, output_dir / "barcode_matching.txt.gz"
    )

    slideseq.bead_matching.write_barcode_xy(
        degen_bead_barcodes, bead_xy, output_dir / "barcode_coordinates.txt.gz"
    )

    filtered_barcodes = [bc for bc in degen_bead_barcodes if bc in umis_per_bead]
    bead_xy_a = np.vstack([bead_xy[dbc] for dbc in filtered_barcodes])

    log.info("Writing output files")
    beads = sorted(umis_per_bead)
    tags = sorted({t for b in beads for t in umis_per_bead[b]})

    with gzip.open(output_dir / "beads.txt.gz", "wt") as out:
        print("\n".join(beads), file=out)

    with gzip.open(output_dir / "tags.txt.gz", "wt") as out:
        print("\n".join(tags), file=out)

    log.debug("Writing umi matrix")
    m = write_matrix(umis_per_bead, beads, tags, output_dir / "umi_matrix.mtx.gz")

    log.debug("Writing read matrix")
    write_matrix(reads_per_bead, beads, tags, output_dir / "read_matrix.mtx.gz")

    bead_dist = []
    bead_pairs = []

    # tags that appear in >1 bead
    two_beads = np.asarray((m > 0).sum(0) > 1).squeeze().nonzero()[0]
    for i in two_beads:
        # beads with this tag
        nz_beads = np.asarray(m[:, i].todense()).squeeze().nonzero()[0]
        for bi, bj in itertools.combinations(nz_beads, 2):
            bead_pairs.append((bi, bj))
            bead_dist.append(
                np.sqrt(((bead_xy_a[bi, :] - bead_xy_a[bj, :]) ** 2).sum())
            )

    long_pairs = [b_p for b_p, b_d in zip(bead_pairs, bead_dist) if b_d > min_dist]

    log.debug(f"Identified {len(long_pairs)} pairs with distance > {min_dist}")

    if extra_pdf:
        extra_pdf = output_dir / "extra_plots.pdf"
        log.info(f"Saving extra plots to {extra_pdf}")
        extra_pdf = PdfPages(extra_pdf)

        umi_counts = Counter(r[32:41] for r in r1_reads)
        log.debug(
            f"Found {len(umi_counts)} UMIs with {sum(umi_counts.values())} total counts"
        )

        plot_log_hist(umi_counts.values(), "Reads per UMI", extra_pdf)

        plot_log_hist(
            [sum(umis_per_bead[bc].values()) for bc in umis_per_bead],
            "UMIs per bead",
            extra_pdf,
        )

        plot_hist(bead_dist, "Distance distribution for paired tags", extra_pdf)

        extra_pdf.close()

    pdf_pages = PdfPages(output_pdf)

    log.info("Making plots")

    umi_dist = [sum(umis_per_bead[bc].values()) for bc in filtered_barcodes]
    read_dist = [sum(reads_per_bead[bc].values()) for bc in filtered_barcodes]

    spatial_plot_plus_links(
        bead_xy_a,
        long_pairs,
        umi_dist,
        f"UMIs per bead and pairs > {min_dist} pixels",
        pdf_pages,
        pct=percentile,
    )

    spatial_plot(
        bead_xy_a,
        umi_dist,
        "UMIs per bead",
        pdf_pages,
        pct=percentile,
    )

    spatial_plot(
        bead_xy_a,
        np.log10(umi_dist),
        "log10 UMIs per bead",
        pdf_pages,
        pct=percentile,
    )

    spatial_plot(
        bead_xy_a,
        read_dist,
        "Reads per bead",
        pdf_pages,
        pct=percentile,
    )

    spatial_plot(
        bead_xy_a,
        np.log10(read_dist),
        "log10 reads per bead",
        pdf_pages,
        pct=percentile,
    )

    pdf_pages.close()
    log.info("Done!")


if __name__ == "__main__":
    main()
