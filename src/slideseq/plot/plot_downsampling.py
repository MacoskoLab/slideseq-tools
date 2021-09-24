#!/usr/bin/python

# edit29 Ali Qutab
# plot real data and model data with output.x for optimal alpha and beta values
# model plot is only a line now, no marker

import logging
from pathlib import Path

import matplotlib.figure
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg

from slideseq.plot import read_dge_summary

import scipy.optimize

log = logging.getLogger(__name__)
# pycharmedit

def plot_downsampling(downsampling_output: list[tuple[float, Path]], figure_path: Path):
    xy = []

    # right now barcodes is a list
    bc_list, full_umis_per_bc, _ = read_dge_summary(Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04.matched.digital_expression_summary.txt"))
    # this is a set comprehension, so we can remove the -1 from the matched barcodes
    # bc.split("-") will split it into two parts, and we take the first one
    bc_set = {bc.split("-")[0] for bc in bc_list}

    for r, downsample_summary in downsampling_output:
        # read the barcodes and counts from this downsampled file
        barcodes, umis_per_bc, _ = read_dge_summary(downsample_summary)
        filtered_barcodes = []
        filtered_umis_per_bc = []
        # we zip together those lists so that we have each value as a pair
        for bc, umis in zip(barcodes, umis_per_bc):
            # checking for this is fast because it's a set
            if bc in bc_set:
                filtered_barcodes.append(bc)
                filtered_umis_per_bc.append(umis)
                # take all barcodes as representative of real cells
        data = np.mean(filtered_umis_per_bc)
        xy.append((r, data))

    xy.sort()
    x, y = zip(*xy)
    x = np.array(x)
    y = np.array(y)

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

    output = scipy.optimize.least_squares(
        model_least_squares,
        [-10.0, 10.],  # initial values for params
        bounds=([-np.inf, 0], [0, np.inf]),  # some bounds on params. alpha should be negative, beta is positive
        kwargs={"data": y, "r": x},
        method="dogbox",  # I found this method to work well for this problem
    )

    params = output.x

    x_values = np.linspace(0.1, 3.0, 30)  # this function creates a linear space of points: 30 points from 0.1 to 3.0 (0.1, 0.2, ... 2.9, 3.0)
    predicted_y = model(x_values, params)

    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    ax.scatter(x, y, marker="o", alpha=0.8, color="r")  # red scatter plot for actual data r = 0.1...1.0
    ax.plot(x_values, predicted_y, alpha=0.8, color="b")  # blue line plot for model data r = 0.1...3.0
    ax.set_xlabel("Subsampling Ratio")
    ax.set_ylabel("Transcripts per barcode")
    ax.set_title("Average transcripts for matched barcodes")

    ax.set_xlim(0.0, 3.1)  # was 1.1, but x_values for model go up to 3.0

    FigureCanvasAgg(fig).print_figure(figure_path)

if __name__ == "__main__":
    # call the function here with input
    plot_downsampling(downsampling_output =
                      [(0.1, Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04_0.1.digital_expression_summary.txt")),
                       (0.2, Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04_0.2.digital_expression_summary.txt")),
                       (0.3, Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04_0.3.digital_expression_summary.txt")),
                       (0.4, Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04_0.4.digital_expression_summary.txt")),
                       (0.5, Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04_0.5.digital_expression_summary.txt")),
                       (0.6, Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04_0.6.digital_expression_summary.txt")),
                       (0.7, Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04_0.7.digital_expression_summary.txt")),
                       (0.8, Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04_0.8.digital_expression_summary.txt")),
                       (0.9, Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04_0.9.digital_expression_summary.txt"))],
                      figure_path = Path("/Users/aqutab/aq/aq_downsampling/aq_plots/aq_edit29_plot_downsampling.png"))