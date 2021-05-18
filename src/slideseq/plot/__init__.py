import csv
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

import matplotlib.figure
from matplotlib.axes import Axes
from matplotlib.backends.backend_pdf import PdfPages

from slideseq.util import constants as constants


@contextmanager
def new_ax(pdf_pages: PdfPages, include_fig=False):
    """
    Simple helper to make a new figure and provide the axes for plotting,
    then save to an open PDF when finished.
    """
    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax: Axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    if include_fig:
        yield fig, ax
    else:
        yield ax

    pdf_pages.savefig(fig)


def read_dge_summary(summary_file):
    with summary_file.open() as fh:
        rdr = csv.DictReader(
            (line for line in fh if len(line) > 1 and line[0] != "#"), delimiter="\t"
        )

        barcodes, umis_per_bc, genes_per_bc = zip(
            *(
                (r["CELL_BARCODE"], int(r["NUM_TRANSCRIPTS"]), int(r["NUM_GENES"]))
                for r in rdr
            )
        )
    return barcodes, umis_per_bc, genes_per_bc


def read_dropseq_metrics(metrics_file: Path, key_dict: dict[str, str] = None):
    """Reads the output from a dropseq analysis command, which is in a standard
    output. This works for ReadQualityMetrics.txt, fracExonicIntronic.txt and
    polyA_filtering.summary.txt

    :param metrics_file: File to read. Should be text with # as comment char,
                         metric keys on first non-comment row and values on second,
                         then a histogram
    :param key_dict: Dictionary of strings to translate the keys for metrics
    """

    if key_dict is None:
        key_dict = dict()  # don't translate any names

    with metrics_file.open() as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and r[0][0] != "#"]

        # translate keys to nicer names
        metrics = {
            key_dict[k]: float(v)
            for k, v in zip(rows[0][1:], rows[1][1:])
            if k in key_dict
        }
        # histogram is of form `{bin: num reads}`
        histogram = Counter({int(r[0]): float(r[1]) for r in rows[3:]})

    return metrics, histogram


def read_quality_metrics(metrics_file: Path):
    """Reads the ReadQualityMetrics.txt file"""

    return read_dropseq_metrics(metrics_file, constants.READ_QUALITY_METRICS)


def read_frac_intronic_exonic(metrics_file: Path):
    """Reads the fracIntronicExonic.txt file"""

    return read_dropseq_metrics(metrics_file, constants.FRAC_INTRONIC_EXONIC)
