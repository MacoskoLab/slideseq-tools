import gzip
import logging
from collections import Counter
from pathlib import Path

import click
import networkx as nx
import numpy as np
import scipy.sparse
from sklearn.neighbors import radius_neighbors_graph

log = logging.getLogger(__name__)


def degen_barcode(barcode_set: set[str]):
    """Given a set of barcodes return the degenerate sequence that matches all of them,
    where N is used for any base that isn't uniform across the set

    :param barcode_set: A degenerate barcode set, e.g. a group of barcodes that
                        are considered equivalent and should be merged
    :return: a single barcode with the consensus nucleotide when possible or N if not
    """
    return "".join(s.pop() if len(s) == 1 else "N" for s in map(set, zip(*barcode_set)))


def initial_h_set(*barcodes: str):
    base_d = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}

    return {tuple(base_d[c] for c in barcode) for barcode in barcodes}


def hamming_set(h_set: set[tuple[int]], d: int = 1, include_N: bool = True):
    """Given an index of bases in {ACGTN}, generate all indexes within hamming
    distance d of the input

    :param h_set: initial set of barcodes to file the hamming set for, encoded as tuples
                  of integers that index into the ACGTN ordering.
    :param d: maximum distance to allow
    :param include_N: include N when generating hamming ball
    :return: set of indexes within hamming distance d
    """

    h_set = h_set.copy()  # so as not to modify in-place, which might startle people
    bc_len = len(next(iter(h_set)))

    # method for generating a hamming ball: turn your barcode into an array. Make
    # new_base array for each position that consists of base i on the diagonal and
    # the other_bases = (1 - eye) to grab the rest

    # use array broadcasting to add new_base + a * other_bases, i.e. the set of arrays
    # with the new base in each position plus the set of arrays with the other bases
    # in their positions
    # e.g. if starting from [2 3 0 1 1]
    # [0 0 3 0 0] + [1 1 0 1 1] * [2 3 0 1 1] = [2 3 3 1 1] which is distance 1 away

    # because we collapse everything in a set, some redundancy is not a problem
    # and this is a lot simpler than checking everything up front (and pretty fast)

    new_base = [i * np.eye(bc_len, dtype=np.uint8) for i in range(4 + include_N)]
    other_bases = 1 - np.eye(bc_len, dtype=np.uint8)

    for _ in range(d):
        for a in list(map(np.array, h_set)):
            h_set.update(t for i in new_base for t in map(tuple, a * other_bases + i))

    h_set = {"".join("ACGTN"[i] for i in h) for h in h_set}

    return h_set


def hamming1_adjacency(barcodes: list[str]):
    """Constructs a sparse adjacency matrix for all pairs of barcodes within Hamming
    distance 1. Does not include self-edges.

    :param barcodes: List of barcodes to compute. The ordering will determine the index
                     of the resulting matrix
    :return: Adjacency matrix in sparse CSR format
    """
    assert len(barcodes) == len(set(barcodes)), "barcodes should be unique"

    bc_to_i = {bc: i for i, bc in enumerate(barcodes)}
    adjacency = scipy.sparse.dok_matrix((len(barcodes), len(barcodes)), dtype=int)

    for i, bead_bc in enumerate(barcodes):
        for bc1 in hamming_set(initial_h_set(bead_bc), d=1, include_N=True):
            if bc1 in bc_to_i and bead_bc != bc1:
                adjacency[i, bc_to_i[bc1]] = 1

    return adjacency.tocsr()


def bipartite_matching(
    bead_barcodes: list[str],
    degen_barcodes: list[str],
    bead_groups: list[set[int]],
    seq_barcodes: list[str],
):
    bead_lens = set(map(len, bead_barcodes))
    seq_lens = set(map(len, seq_barcodes))
    assert (
        bead_lens == seq_lens
    ), f"Beads have length {bead_lens} but sequenced lengths {seq_lens}"
    log.debug(msg=f"Barcodes all have length {bead_lens}")

    assert {c for bc in bead_barcodes for c in bc}.issubset("ACGTN")

    seq_nset = {c for bc in seq_barcodes for c in bc}
    assert seq_nset.issubset("ACGTN")

    # if N does not appear in the sequence data we can save a little time
    # when we make our graph, because none of the queries will have it
    include_N = "N" in seq_nset

    # construct bipartite graph: barcode to hamming set
    matching_graph = nx.Graph()
    for i, bead_bg in enumerate(bead_groups):
        matching_graph.add_node(i, bipartite=0)
        h_set = initial_h_set(*(bead_barcodes[j] for j in bead_bg))

        for bc1 in hamming_set(h_set, d=1, include_N=include_N):
            matching_graph.add_node(bc1, bipartite=1)
            matching_graph.add_edge(i, bc1)

    log.debug(msg=f"Created bipartite graph with {matching_graph.size()} edges")

    barcode_mapping = dict()

    # find unambiguous matches from sequencing: <= hamming distance 1 of exactly one group
    for seq_bc in seq_barcodes:
        if seq_bc in matching_graph:
            sample_set = set(matching_graph[seq_bc])

            if len(sample_set) == 1:
                barcode_mapping[seq_bc] = degen_barcodes[sample_set.pop()]

    return barcode_mapping


def match_barcodes(
    seq_barcode_file: Path,
    bead_barcode_file: Path,
    bead_location_file: Path,
    radius: float = 10.0,
):
    with bead_barcode_file.open() as fh:
        bead_barcodes = ["".join(line.strip().split(",")) for line in fh]

    with bead_location_file.open() as fh:
        # this file has x and y on two super long lines
        x_line = np.array([float(v) for v in fh.readline().strip().split(",")])
        y_line = np.array([float(v) for v in fh.readline().strip().split(",")])
        xy = np.vstack((x_line, y_line)).T

    log.debug(
        msg=f"Read {len(bead_barcodes)} ({len(set(bead_barcodes))}) bead barcodes"
    )

    assert xy.shape[0] == len(
        bead_barcodes
    ), f"Got {xy.shape[0]} bead locations for {len(bead_barcodes)} beads"

    # pre-emptively remove poly-T/N sequences
    ok_barcodes = [not set(bc).issubset({"T", "N"}) for bc in bead_barcodes]
    xy = xy[ok_barcodes, :]
    bead_barcodes = [bc for ok, bc in zip(ok_barcodes, bead_barcodes) if ok]

    assert xy.shape[0] == len(
        bead_barcodes
    ), "Removed different numbers of T/N barcodes somehow"

    with gzip.open(seq_barcode_file, mode="rt") as fh:
        seq_barcodes = sorted(line.strip().split("-")[0] for line in fh)
        # remove poly-T sequence if present
        seq_barcodes = [seq for seq in seq_barcodes if set(seq) != {"T"}]

    log.debug(
        msg=f"Read {len(seq_barcodes)} ({len(set(seq_barcodes))}) sequencing barcodes"
    )

    # adjacency matrix for all beads within radius of each other
    radius_matrix = radius_neighbors_graph(xy, radius=radius)
    log.debug(f"radius matrix (<{radius}) with {radius_matrix.nnz // 2} edges")
    # adjacency matrix for all barcodes within hamming distance 1
    hamming_matrix = hamming1_adjacency(bead_barcodes)
    log.debug(f"hamming matrix (<=1) with {hamming_matrix.nnz // 2} edges")

    # just multiply together to get the combined adjacency matrix!
    combined_graph = nx.from_scipy_sparse_matrix(radius_matrix.multiply(hamming_matrix))
    log.debug(f"combined graph with {combined_graph.number_of_edges()} edges")

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

    # check if we have new collisions
    log.debug(
        msg=(
            f"Collapsed {len(bead_groups)} bead groups into"
            f" {len(set(degen_bead_barcodes))} barcodes"
        )
    )

    # just in case, we'll add integer tags to each one so they are unique
    barcode_counter = Counter()
    for i, barcode in enumerate(degen_bead_barcodes):
        barcode_counter[barcode] += 1
        degen_bead_barcodes[i] = f"{barcode}-{barcode_counter[barcode]}"

    # average xy for grouped beads to get centroids
    bead_xy = dict()
    for bg, degen_bc in zip(bead_groups, degen_bead_barcodes):
        bg_graph = combined_graph.subgraph(bg)
        mean_x, mean_y = np.array(
            [[nd["x"], nd["y"]] for _, nd in bg_graph.nodes(data=True)]
        ).mean(0)
        bead_xy[degen_bc] = (mean_x, mean_y)

    log.info(f"Found {len(bead_groups)} groups of connected beads")
    log.debug(f"Size distribution: {Counter(map(len, bead_groups))}")

    barcode_matching = bipartite_matching(
        bead_barcodes, degen_bead_barcodes, bead_groups, seq_barcodes
    )

    # filter degen_bead_barcodes to only barcodes that were observed
    matched_set = set(barcode_matching.values())
    degen_bead_barcodes = [bc for bc in degen_bead_barcodes if bc in matched_set]

    return degen_bead_barcodes, barcode_matching, bead_xy, combined_graph


def write_barcode_mapping(barcode_mapping, bead_xy, output_file: Path):
    """Write the mapping of sequenced barcodes to collapsed bead barcodes,
    along with xy coordinates"""
    with gzip.open(output_file, "wt") as out:
        for seq_bc, bead_bc in barcode_mapping.items():
            x, y = bead_xy[bead_bc]
            print(f"{seq_bc}\t{bead_bc}\t{x:.1f}\t{y:.1f}", file=out)


def write_barcode_xy(bead_list, bead_xy, output_file: Path):
    """Write the collapsed bead barcodes with their coordinates, in the specified order"""
    with gzip.open(output_file, "wt") as out:
        for bead_bc in bead_list:
            x, y = bead_xy[bead_bc]
            print(f"{bead_bc}\t{x:.1f}\t{y:.1f}", file=out)


@click.command(name="bead_matching")
@click.option(
    "-s",
    "--sequence-barcodes",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "-b",
    "--bead-barcodes",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "-l",
    "--bead-locations",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "-o",
    "--output-file",
    required=True,
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
)
@click.option("--radius", default=10.0, type=float, show_default=True)
def main(
    sequence_barcodes: str,
    bead_barcodes: str,
    bead_locations: str,
    output_file: str,
    radius: float = 10.0,
):
    log.setLevel(logging.DEBUG)
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    stream_handler.setFormatter(formatter)

    log.addHandler(stream_handler)
    log.info(msg="Added stream handler")

    sequence_barcodes = Path(sequence_barcodes)
    bead_barcodes = Path(bead_barcodes)
    bead_locations = Path(bead_locations)
    output_file = Path(output_file)

    barcode_list, barcode_mapping, bead_xy, bead_graph = match_barcodes(
        sequence_barcodes, bead_barcodes, bead_locations, radius
    )

    log.info(msg=f"Found {len(set(barcode_mapping.values()))} matches")

    write_barcode_mapping(barcode_mapping, bead_xy, output_file)

    nx.write_gpickle(bead_graph, output_file.with_suffix(".bead_graph.pickle.gz"))


if __name__ == "__main__":
    main()
