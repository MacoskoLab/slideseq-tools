#!/usr/bin/python

# This script is to generate mtx files for sparse digital expression matrix

import csv
import gzip
import logging

import scipy.io
import scipy.sparse

from slideseq.library import Library

log = logging.getLogger(__name__)


def write_sparse_matrix(library: Library):
    # read gene id and name mapping from gtf file
    gene_dict = dict()
    with library.reference.annotations.open("r") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        for row in rdr:
            if not row or row[0][0] == "#":
                continue

            md = [v.strip().split(maxsplit=1) for v in row[8].strip().split(";") if v]
            md = dict((k, v[1:-1]) for k, v in md)

            # this is going to pick the first gene id for a gene name and ignore others
            # that's kinda weird but not going to mess with it for now
            if md["gene_name"] not in gene_dict:
                gene_dict[md["gene_name"]] = md["gene_id"]

    # get cols and rows from dge file
    with gzip.open(library.matched.digital_expression, "rt") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        cols = next(rdr)[1:]
        rows = [r[0] for r in rdr]

    with library.matched.barcodes.open("w") as out:
        for bc in cols:
            print(bc, file=out)

    # write features (genes) file
    with library.matched.genes.open("w") as out:
        for gene in rows:
            print(f"{gene_dict.get(gene, gene)}\t{gene}", file=out)

    # read in DGE as sparse matrix
    data = scipy.sparse.dok_matrix((len(rows), len(cols)), dtype=int)
    with gzip.open(library.matched.digital_expression, "rt") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        _ = next(rdr)
        for i, row in enumerate(rdr):
            for j, val in enumerate(row[1:]):
                data[i, j] = int(val)

    # write mtx file
    with gzip.open(library.matched.mtx, "wb") as out:
        scipy.io.mmwrite(out, data.tocsr())
