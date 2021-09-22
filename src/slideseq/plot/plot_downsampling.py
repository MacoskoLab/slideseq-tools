#!/usr/bin/python

# edit5 Ali Qutab
# This script is to generate PDF for downsampling with subsampling code for better predictions.
# optimizing a function called model, which goes from ratio to UMI count.

import logging
from pathlib import Path

import matplotlib.figure
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg

from slideseq.plot import read_dge_summary

log = logging.getLogger(__name__)
# pycharmedit

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

    def model(r, params):
        """
        y = alpha * exp(-r) + beta
        y is average transcript count for the 0.1,0.2,... files
        r = 0.1, 0.2,...,1.0 - for each file of data given
        params is an array of two values [alpha, beta] which we will optimize
        """
        return params[0] * np.exp(-r) + params[1]

    """
    model function takes 0.1, 0.2 ... 1.0 for r and outputs the mean counts values which are y,
    based on some parameters alpha and beta, which we need to find to predict the result for any value of r.
    """

    """
    The least_squares function is to optimize the model function, we want to minimize the error, which is the difference between the model prediction and the data.
    """

    def model_least_squares(params, *, r, data):
        # computes model and compares it to the data
        return model(r, params) - data

    import scipy.optimize

    scipy.optimize.least_squares(
        model_least_squares,
        [-10.0, 10.],  # initial values for params
        bounds=([-np.inf, 0], [0, np.inf]),  # some bounds on params. alpha should be negative, beta is positive
        kwargs={"data": y, "r": x},
        method="dogbox",  # I found this method to work well for this problem
    )

if __name__ == "__main__":
    # call the plot_downsampling() function here with input
    plot_downsampling(downsampling_output=[
        (0.1, Path(
            "/Users/aqutab/slideseq-tools/aqutab/aq_downsampling/aqutab_files/Puck_210203_04_0.1.digital_expression_summary.txt")),
        (0.2, Path(
            "/Users/aqutab/slideseq-tools/aqutab/aq_downsampling/aqutab_files/Puck_210203_04_0.2.digital_expression_summary.txt")),
        (0.3, Path(
            "/Users/aqutab/slideseq-tools/aqutab/aq_downsampling/aqutab_files/Puck_210203_04_0.3.digital_expression_summary.txt")),
        (0.4, Path(
            "/Users/aqutab/slideseq-tools/aqutab/aq_downsampling/aqutab_files/Puck_210203_04_0.4.digital_expression_summary.txt")),
        (0.5, Path(
            "/Users/aqutab/slideseq-tools/aqutab/aq_downsampling/aqutab_files/Puck_210203_04_0.5.digital_expression_summary.txt")),
        (0.6, Path(
            "/Users/aqutab/slideseq-tools/aqutab/aq_downsampling/aqutab_files/Puck_210203_04_0.6.digital_expression_summary.txt")),
        (0.7, Path(
            "/Users/aqutab/slideseq-tools/aqutab/aq_downsampling/aqutab_files/Puck_210203_04_0.7.digital_expression_summary.txt")),
        (0.8, Path(
            "/Users/aqutab/slideseq-tools/aqutab/aq_downsampling/aqutab_files/Puck_210203_04_0.8.digital_expression_summary.txt")),
        (0.9, Path(
            "/Users/aqutab/slideseq-tools/aqutab/aq_downsampling/aqutab_files/Puck_210203_04_0.9.digital_expression_summary.txt"))
    ], figure_path="aq_edit5_plot_downsampling.png")
