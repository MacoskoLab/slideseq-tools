#!/usr/bin/python

# edit51 Ali Qutab
# fit the model once for each of the y variables, in the same for-loop as the plotting
# use the parameters immediately to get pred_y

import logging

import matplotlib.figure
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg

from slideseq.plot import read_dge_summary

import scipy.optimize

import argparse
import os
from pathlib import Path

log = logging.getLogger(__name__)

def plot_downsampling(downsampling_output: list[tuple[float, Path]], figure_path: Path):
    xy = []

    # right now barcodes is a list
    bc_list, full_umis_per_bc, _ = read_dge_summary(matched_path)
    # this is a set comprehension, so we can remove the -1 from the matched barcodes
    # bc.split("-") will split it into two parts, and we take the first one
    bc_set = {bc.split("-")[0] for bc in bc_list}

    for r, downsample_summary in downsampling_output:
        # read the barcodes and counts from this downsampled file
        barcodes, umis_per_bc, _ = read_dge_summary(downsample_summary)
        filtered_barcodes = [] # will append barcodes that match barcodes in matched expression sumarry file
        filtered_umis_per_bc = [] # will append UMIs that match barcodes in matched expression sumarry file
        # we zip together those lists so that we have each value as a pair
        for bc, umis in zip(barcodes, umis_per_bc):
            # checking for this is fast because it's a set
            if bc in bc_set:
                filtered_barcodes.append(bc)
                filtered_umis_per_bc.append(umis)
                # take all barcodes as representative of real cells
        data_100 = np.mean(filtered_umis_per_bc)
        data_80 = np.mean(filtered_umis_per_bc[:(int(0.8 * (len(filtered_umis_per_bc))))])  # for example, top 80% is from top item to item at 0.8 * number of items in list filtered_umis_per_bc
        data_60 = np.mean(filtered_umis_per_bc[:(int(0.6 * (len(filtered_umis_per_bc))))])
        data_40 = np.mean(filtered_umis_per_bc[:(int(0.4 * (len(filtered_umis_per_bc))))])
        data_20 = np.mean(filtered_umis_per_bc[:(int(0.2 * (len(filtered_umis_per_bc))))])

        xy.append((r, data_100, data_80, data_60, data_40, data_20))

    xy.sort()
    x, y_100, y_80, y_60, y_40, y_20 = zip(*xy)
    # convert to arrays
    x = np.array(x)
    y_100 = np.array(y_100)
    y_80 = np.array(y_80)
    y_60 = np.array(y_60)
    y_40 = np.array(y_40)
    y_20 = np.array(y_20)

    def model(r, params):

        # y = alpha * exp(-r) + beta
        # y is average transcript count for the 0.1,0.2,... files
        # r = 0.1, 0.2,...,1.0 - for each file of data given
        # params is an array of two values [alpha, beta] which we will optimize

        return params[0] * np.exp(-r) + params[1]

    # model function takes 0.1, 0.2 ... 1.0 for r and outputs the mean counts values which are y,
    # based on some parameters alpha and beta, which we need to find to predict the result for any value of r.

    # The least_squares function is to optimize the model function, we want to minimize the error, which is the difference between the model prediction and the data.


    def model_least_squares(params, *, r, data):
        # computes model and compares it to the data
        return model(r, params) - data

    x_values = np.linspace(0.1, 3.0,30)  # this function creates a linear space of points: 30 points from 0.1 to 3.0 (0.1, 0.2, ... 2.9, 3.0)

    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    for y, label in zip(
            (y_100, y_80, y_60, y_40, y_20),
            ("100", "80", "60", "40", "20")
    ):
        output = scipy.optimize.least_squares(
            model_least_squares,
            [-10.0, 10.],  # initial values for params
            bounds=([-np.inf, 0], [0, np.inf]),  # some bounds on params. alpha should be negative, beta is positive
            kwargs={"data": y, "r": x},
            method="dogbox",  # I found this method to work well for this problem
        )

        params = output.x

        pred_y = model(x_values, params)
        # predicted_y_80 = model_80(x_values, params_80)
        # predicted_y_60 = model_80(x_values, params_60)
        # predicted_y_40 = model_80(x_values, params_40)
        # predicted_y_20 = model_80(x_values, params_20)

        ax.scatter(x, y, marker="o", alpha=0.8, color="r")
        ax.plot(x_values, pred_y, label=f"top {label}%", alpha=0.8, color="r")

    ax.set_xlabel("Subsampling Ratio")
    ax.set_ylabel("Transcripts per matched barcode")
    ax.set_title("Average transcripts for matched barcodes")

    ax.set_xlim(0.0, 3.1)  # r went upto 1.0 for actual data, but x_values for model go up to 3.0

    ax.legend()

    FigureCanvasAgg(fig).print_figure(figure_path)

if __name__ == "__main__":
    # use argparse to get a list of files from the command line, as strings
    parser = argparse.ArgumentParser(description='Read in downsample_summary text files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('path', nargs='+', help='Path of a downsample_summary text file')
    parser.add_argument("--output", help="output filename")
    args = parser.parse_args()

    downsampling_list = []  # empty list outside the for loops
    # Parse paths
    # downsample_summary is a string containing the full path of a filename
    # for matched expression summary file
    for matched_expression_summary in args.path:
        if matched_expression_summary.find("matched") > -1:
            matched_path = Path(matched_expression_summary)
            print(matched_path)
            basename = os.path.basename(matched_path)
            print(basename)
            # downsampling_list.append((1.0, matched_path)) # to get 1.0 point in the plot, using the matched data

    # for 0.1,0.2,... files
    for downsample_summary in args.path:
        split_path = downsample_summary.split("_")  # split on _
        characters = split_path[5] # grab the 5th item in the list split by '_'
        if downsample_summary.find("matched") == -1:
            r = float(characters[0:3]) # first three characters are 0.1,0.2,...
            # make the list of (float, Path) to pass into the main plotting function
            downsample_summary_path = Path(downsample_summary) # convert downsample_summary to path object
            downsampling_list.append((r, downsample_summary_path))

    # trying to use arguments for ratio and path instead of hardcoding them into the script
    plot_downsampling(downsampling_list, figure_path=Path("/Users/aqutab/aq/aq_downsampling/aq_plots/aq_edit51_plot_downsampling.png"))
