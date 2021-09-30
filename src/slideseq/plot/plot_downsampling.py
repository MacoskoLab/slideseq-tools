#!/usr/bin/python

# edit32 Ali Qutab
# plot real data and model data
# now have five quantiles:
# this script edit plots five quantiles by fitting model for the top 20%, the top 40%, 60%, 80%, 100%
# coneptually have five different values of data and fit the model five times, plot five lines and five sets of points

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
    x_100y_100 = []
    x_80y_80 = []

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
        data_100 = np.mean(filtered_umis_per_bc)
        data_80 = np.mean(filtered_umis_per_bc[:(int(0.8 * (len(filtered_umis_per_bc))))]) # for example, top 80% is from top item to item at 0.8 * number of items in list filtered_umis_per_bc
        """
        data_60 = 
        data_40 = 
        data_20 = 
        """
        x_100y_100.append((r, data_100))
        x_80y_80.append((r, data_80))

    x_100y_100.sort()
    x_100, y_100 = zip(*x_100y_100)
    x_100 = np.array(x_100)
    y_100 = np.array(y_100)

    x_80y_80.sort()
    x_80, y_80 = zip(*x_80y_80)
    x_80 = np.array(x_80)
    y_80 = np.array(y_80)

    def model_100(r, params_100):
        """
        y = alpha * exp(-r) + beta
        y is average transcript count for the 0.1,0.2,... files
        r = 0.1, 0.2,...,1.0 - for each file of data given
        params is an array of two values [alpha, beta] which we will optimize
        """
        return params_100[0] * np.exp(-r) + params_100[1]

    def model_80(r, params_80):
        """
        y = alpha * exp(-r) + beta
        y is average transcript count for the 0.1,0.2,... files
        r = 0.1, 0.2,...,1.0 - for each file of data given
        params is an array of two values [alpha, beta] which we will optimize
        """
        return params_80[0] * np.exp(-r) + params_80[1]

    """
    model function takes 0.1, 0.2 ... 1.0 for r and outputs the mean counts values which are y,
    based on some parameters alpha and beta, which we need to find to predict the result for any value of r.
    """

    """
    The least_squares function is to optimize the model function, we want to minimize the error, which is the difference between the model prediction and the data.
    """

    def model_least_squares_100(params_100, *, r, data_100):
        # computes model and compares it to the data
        return model_100(r, params_100) - data_100

    def model_least_squares_80(params_80, *, r, data_80):
        # computes model and compares it to the data
        return model_80(r, params_80) - data_80

    output_100 = scipy.optimize.least_squares(
        model_least_squares_100,
        [-10.0, 10.],  # initial values for params
        bounds=([-np.inf, 0], [0, np.inf]),  # some bounds on params. alpha should be negative, beta is positive
        kwargs={"data_100": y_100, "r": x_100},
        method="dogbox",  # I found this method to work well for this problem
    )

    output_80 = scipy.optimize.least_squares(
        model_least_squares_80,
        [-10.0, 10.],  # initial values for params
        bounds=([-np.inf, 0], [0, np.inf]),  # some bounds on params. alpha should be negative, beta is positive
        kwargs={"data_80": y_80, "r": x_80},
        method="dogbox",  # I found this method to work well for this problem
    )

    params_100 = output_100.x
    params_80 = output_80.x

    x_values = np.linspace(0.1, 3.0, 30)  # this function creates a linear space of points: 30 points from 0.1 to 3.0 (0.1, 0.2, ... 2.9, 3.0)
    predicted_y_100 = model_100(x_values, params_100)
    predicted_y_80 = model_80(x_values, params_80)

    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    ax.scatter(x_100, y_100, marker="o", alpha=0.8, color="r")  # red scatter plot for actual data r = 0.1...1.0
    ax.plot(x_values, predicted_y_100, alpha=0.8, color="b")  # blue line plot for model data r = 0.1...3.0
    ax.plot(x_values, predicted_y_80, alpha=0.8, color="g")

    ax.set_xlabel("Subsampling Ratio")
    ax.set_ylabel("Transcripts per matched barcode")
    ax.set_title("Average transcripts for matched barcodes")

    ax.set_xlim(0.0, 3.1)  # r went upto 1.0 for actual data, but x_values for model go up to 3.0

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
                      figure_path = Path("/Users/aqutab/aq/aq_downsampling/aq_plots/aq_edit32_plot_downsampling.png"))