#!/usr/bin/python

import argparse
import logging
import os
from pathlib import Path

import matplotlib.figure
import numpy as np
import scipy.optimize
from matplotlib.backends.backend_agg import FigureCanvasAgg

from slideseq.plot import read_dge_summary

import matplotlib.pyplot as plt

from matplotlib.offsetbox import AnchoredText

log = logging.getLogger(__name__)


def plot_downsampling(downsampling_output: list[tuple[float, Path]], matched_path: Path, figure_path: Path):
    xy = []

    # right now barcodes is a list
    bc_list, full_umis_per_bc, _ = read_dge_summary(matched_path)
    # this is a set comprehension, so we can remove the -1 from the matched barcodes
    # bc.split("-") will split it into two parts, and we take the first one
    bc_set = {bc.split("-")[0] for bc in bc_list}
    data_100 = np.mean(full_umis_per_bc)
    data_80 = np.mean(
        full_umis_per_bc[: (int(0.8 * (len(full_umis_per_bc))))]
    )  # for example, top 80% is from top item to item at 0.8 * number of items in list filtered_umis_per_bc
    data_60 = np.mean(full_umis_per_bc[: (int(0.6 * (len(full_umis_per_bc))))])
    data_40 = np.mean(full_umis_per_bc[: (int(0.4 * (len(full_umis_per_bc))))])
    data_20 = np.mean(full_umis_per_bc[: (int(0.2 * (len(full_umis_per_bc))))])
    xy.append((1.0, data_100, data_80, data_60, data_40, data_20))

    for r, downsample_summary in downsampling_output:
        if r == 1.0:
            continue
        else:
            # read the barcodes and counts from this downsampled file
            barcodes, umis_per_bc, _ = read_dge_summary(downsample_summary)
            filtered_barcodes = (
                []
            )  # will append barcodes that match barcodes in matched expression sumarry file
            filtered_umis_per_bc = (
                []
            )  # will append UMIs that match barcodes in matched expression sumarry file
            # we zip together those lists so that we have each value as a pair
            for bc, umis in zip(barcodes, umis_per_bc):
                # checking for this is fast because it's a set
                if bc in bc_set:
                    filtered_barcodes.append(bc)
                    filtered_umis_per_bc.append(umis)
                    # take all barcodes as representative of real cells
            data_100 = np.mean(filtered_umis_per_bc)
            data_80 = np.mean(
                filtered_umis_per_bc[: (int(0.8 * (len(filtered_umis_per_bc))))]
            )  # for example, top 80% is from top item to item at 0.8 * number of items in list filtered_umis_per_bc
            data_60 = np.mean(
                filtered_umis_per_bc[: (int(0.6 * (len(filtered_umis_per_bc))))]
            )
            data_40 = np.mean(
                filtered_umis_per_bc[: (int(0.4 * (len(filtered_umis_per_bc))))]
            )
            data_20 = np.mean(
                filtered_umis_per_bc[: (int(0.2 * (len(filtered_umis_per_bc))))]
            )

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

    x_values = np.linspace(
        0.1, 3.0, 30
    )  # this function creates a linear space of points: 30 points from 0.1 to 3.0 (0.1, 0.2, ... 2.9, 3.0)

    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    for y, label, color in zip(
            (y_100, y_80, y_60, y_40, y_20),
            ("100", "80", "60", "40", "20"),
            ("r", "b", "g", "k", "p"),
    ):
        output = scipy.optimize.least_squares(
            model_least_squares,
            [-10.0, 10.0],  # initial values for params
            bounds=(
                [-np.inf, 0],
                [0, np.inf],
            ),  # some bounds on params. alpha should be negative, beta is positive
            kwargs={"data": y, "r": x},
            method="dogbox",  # I found this method to work well for this problem
        )

        params = output.x

        pred_y = model(x_values, params)

        ax.scatter(x, y, marker="o", alpha=0.8)
        ax.plot(x_values, pred_y, label=f"top {label}%", alpha=0.8)

    ax.set_xlabel("Subsampling Ratio")
    ax.set_ylabel("Transcripts per matched barcode")
    ax.set_title("Average transcripts for matched barcodes")

    ax.set_xlim(
        0.0, 3.1
    )  # r went upto 1.0 for actual data, but x_values for model go up to 3.0

    ax.legend()

    # calculate 2x and 10x depth for the 100% model with the ratios model(2) / model(1) and model(10) / model(1)
    r_2 = (model(2.0, params)/model(1.0, params)) # model(2) / model(1)
    r_10 = (model(10.0, params)/model(1.0, params)) # model(10) / model(1)

    # text box, bottom right, for a summary of the return for 2x and 10x depth for the 100% model
    textstr = '\n'.join(((f"{r_2:.1%}"), (f"{r_10:.1%}"))) # used string formatting to use the number from the code.
    # ax.text(0.99, 0.05, textstr, transform=ax.transAxes, fontsize=9,
    # verticalalignment='top', horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    at = AnchoredText(textstr,
                      loc='lower right', prop=dict(size=8), frameon=True,
                      )
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)

    FigureCanvasAgg(fig).print_figure(figure_path)


if __name__ == "__main__":
    # use argparse to get a list of files from the command line, as strings
    parser = argparse.ArgumentParser(
        description="Read in downsample_summary text files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--m_path", type=Path, help="Path of match text file")
    parser.add_argument(
        "ds_path", nargs="+", help="Path of a downsample_summary text file"
    )
    parser.add_argument("--output", help="output filename", required=True)
    args = parser.parse_args()

    downsampling_list = []  # empty list outside the for loops
    # Parse paths
    # downsample_summary is a string containing the full path of a filename

    # for 0.1,0.2,... files
    for downsample_summary in args.ds_path:
        split_path = downsample_summary.split("_")  # split on _
        characters = split_path[5]  # grab the 5th item in the list split by '_'
        r = float(characters[0:3])  # first three characters are 0.1,0.2,...
        # make the list of (float, Path) to pass into the main plotting function
        downsample_summary_path = Path(
            downsample_summary
        )  # convert downsample_summary to path object
        downsampling_list.append((r, downsample_summary_path))

    # trying to use arguments for ratio and path instead of hardcoding them into the script
    plot_downsampling(downsampling_list, matched_path=args.m_path, figure_path=args.output)
