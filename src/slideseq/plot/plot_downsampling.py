#!/usr/bin/python

# edit1 Ali Qutab
# This script is to generate PDF for downsampling

import logging
from pathlib import Path

import matplotlib.figure
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg

from slideseq.plot import read_dge_summary

log = logging.getLogger(__name__)


def plot_downsampling(downsampling_output: list[tuple[float, Path]], figure_path: Path):
    xy = []

    for ratio, downsample_summary in downsampling_output:
        _, umis_per_bc, _ = read_dge_summary(downsample_summary)
        # take the first 10000 barcodes as representative of real cells
        v = np.mean(umis_per_bc[:10000])

        xy.append((ratio, v))

    xy.sort()
    x, y = zip(*xy)

    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    ax.plot(x, y, marker="o", markersize=10, alpha=0.8)
    ax.set_xlabel("Subsampling Ratio")
    ax.set_ylabel("Transcripts per barcode")
    ax.set_title("Average transcripts for top 10,000 barcodes")

    ax.set_xlim(0.0, 1.1)

    FigureCanvasAgg(fig).print_figure(figure_path)
if __name__ == "__main__":
    # call the function here with inpu
    plot_downsampling(downsampling_output = [( 0.1, Path("/Users/aqutab/slideseq-tools/aqutab/AQ_downsampling/aqutab_files/Puck_210203_04_0.1.digital_expression_su    mmary.txt"))], figure_path = "aq_edit1_plot_downsampling.png"),
