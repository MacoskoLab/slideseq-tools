#!/usr/bin/python

# edit4 Ali Qutab
# This script is to generate PDF for downsampling

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

    """
    optimizing a function we care about is model, which goes from ratio to UMI count. 
    """

    def model(r, params):
        # y = alpha * exp(-r) + beta
        # y is average transcript count for the 0.1,0.2,... files
        # y is average # of UMIs for top 10%,20%,40%,60%,80%,100% of barcodes
        # r is a ratio - the input containing different subsampling fractions
        # r = 0.1, 0.2,...,1.0 - for each file of data given
        # params is an array of two values [alpha, beta] which we will optimize
        # alpha params[0] and beta params[1] - want to find the best ones to fit the equation to those data points
        # alpha params[0] and beta params[1] - use the best ones to predict results for other values of r (>1)
        # for example fit our 10 data points (0.1 through 1.0) and found the best fit was
        # alpha = -1000 and beta = 2000,
        # predict the result for 3x read depth as -1000 * exp( -3) + 2000 ~= 1950
        # predict the result for any value of r, so we can draw those lines on the plot
        return params[0] * np.exp(-r) + params[1]

    """
    that takes 0.1, 0.2 ... 1.0 for r and outputs the mean counts values which are y,
    based on some parameters alpha and beta, which we need to find.
    """

    """
    The least_squares function is to optimize that model function, we want to minimize the error, which is the difference between the model prediction and the data.
    """

    def model_least_squares(params, *, r, data):
        # computes the squared error for a given set of parameters
        # and the observed data
        return model(r, params) - data

    """
    computes model and compares it to the data.
    will find optimal values by repeatedly tweaking them until it gets the best fit, and then it returns them to you.
    The input to least_squares should be a simple function that just does some math and returns the result
    """

    import scipy.optimize

    scipy.optimize.least_squares(
        model_least_squares,
        [-10.0, 10.],  # initial values for params
        bounds=([-np.inf, 0], [0, np.inf]),  # some bounds on params. alpha should be negative, beta is positive
        kwargs={"data": y, "r": x},
        method="dogbox",  # I found this method to work well for this problem
    )

    """
    will optimize the parameters to try to minimize the error.
    when you are passing the input to the optimizer, it should be as an np.array of input values and an array of output values.
    like np.array([0.1, 0.2, 0.3, ... 1.0]) and similar for the corresponding outputs.
    """

    """
    Once you have the values, you can get the model prediction for any input by using the same values and calculating the output for r = 2 or whatever you like.
    """

    """
    the plot_downsampling function should be the wrapper for everything
    1 reading the data files
    2 finding the coefficients by optimization
    3 plotting the results.
    """

if __name__ == "__main__":
    # call the function here with input
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
    ], figure_path="aq_edit4_plot_downsampling.png")


