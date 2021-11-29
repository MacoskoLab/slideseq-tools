#!/usr/bin/python

import argparse
import logging
from pathlib import Path

import matplotlib.figure
import numpy as np
import scipy.optimize
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.offsetbox import AnchoredText

from slideseq.plot import read_dge_summary

log = logging.getLogger(__name__)


def plot_downsampling(
    downsampling_output: list[tuple[float, Path]], matched_path: Path, figure_path: Path
):
    # percentiles of top barcodes to compute
    percentiles = np.array([80, 60, 40, 20])

    bc_list, full_umis_per_bc, _ = read_dge_summary(matched_path)
    # set comprehension to remove the -1 from the matched barcodes
    bc_set = {bc.split("-")[0] for bc in bc_list}

    total_umis = sum(full_umis_per_bc)
    n_bcs = len(full_umis_per_bc)

    quintiles = total_umis - np.percentile(
        np.cumsum(full_umis_per_bc[::-1]), percentiles[::-1]
    )
    quintiles /= n_bcs * percentiles / 100

    xy = [(1.0, total_umis / n_bcs, *quintiles)]

    for ratio, downsample_file in downsampling_output:
        # read the barcodes and counts from this downsampled file
        barcodes, umis_per_bc, _ = read_dge_summary(downsample_file)

        # filter to barcodes from the matched expression summary file
        r_umis_per_bc = [
            umis for bc, umis in zip(barcodes, umis_per_bc) if bc in bc_set
        ]

        total_umis = sum(r_umis_per_bc)
        n_bcs = len(r_umis_per_bc)

        # how this works: r_umis_per_bc is in descending order. We reverse it,
        # then calculate the cumulative sum. The percentiles of that sum are
        # equal to the total reads from the bottom 20%, 40%, etc. We subtract
        # from the overall total to get the top 80%, 60%, etc. Then divide to
        # get the mean
        quintiles = total_umis - np.percentile(
            np.cumsum(r_umis_per_bc[::-1]), percentiles[::-1]
        )
        quintiles /= n_bcs * percentiles / 100

        xy.append((ratio, total_umis / n_bcs, *quintiles))

    xy.sort()
    x, *ys = zip(*xy)

    # convert to arrays
    x = np.array(x)
    ys = list(map(np.array, ys))

    def model(r, params):
        # y = alpha * exp(-r) + beta
        # y is average transcript count for the 0.1,0.2,... files
        # r = 0.1, 0.2,...,1.0 - for each file of data given
        # params is an array of two values [alpha, beta] which we will optimize

        return params[0] * np.exp(-r) + params[1]

    def model_least_squares(params, *, r, data):
        # computes model and compares it to the data
        return model(r, params) - data

    x_values = np.linspace(0.1, 3.0, 30)
    parameters = dict()

    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    for y, label in zip(ys, (100, *percentiles)):
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

        parameters[label] = output.x
        pred_y = model(x_values, parameters[label])

        ax.scatter(x, y, marker="o", alpha=0.8)
        ax.plot(x_values, pred_y, label=f"top {label}%", alpha=0.8)

    ax.set_xlabel("Subsampling Ratio")
    ax.set_ylabel("Transcripts per matched barcode")
    ax.set_title("Average transcripts for matched barcodes")

    # r went up to 1.0 for actual data, but x_values for model go up to 3.0
    ax.set_xlim(0.0, 3.1)
    ax.legend()

    # calculate 2x and 10x depth for the 100% model with the
    # ratios model(2) / model(1) and model(10) / model(1)
    r_2 = model(2.0, parameters[20]) / model(1.0, parameters[20])
    r_10 = model(10.0, parameters[20]) / model(1.0, parameters[20])

    # text box, bottom right, for a summary of the return for 2x and 10x depth for the 100% model
    textstr = f"top 20%, 2x depth: {r_2:8.1%}\ntop 20%, 10x depth: {r_10:6.1%}"
    # use AnchoredText to position text box to bottom right
    at = AnchoredText(textstr, loc="lower right", prop=dict(size=10), frameon=True)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    at.patch.set_edgecolor((0.8, 0.8, 0.8))
    ax.add_artist(at)

    FigureCanvasAgg(fig).print_figure(figure_path)


if __name__ == "__main__":
    # use argparse to get a list of files from the command line, as strings
    parser = argparse.ArgumentParser(
        description="Read in downsample_summary text files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "ds_path", nargs="+", type=Path, help="Downsampled summary files"
    )
    parser.add_argument("--m_path", type=Path, help="Matched summary file")
    parser.add_argument("--output", help="output filename", required=True)
    args = parser.parse_args()

    # make the list of (float, Path) to pass into the main plotting function
    downsampling_list = []  # empty list outside the for loops

    # for 0.1, 0.2,... files
    for downsample_summary in args.ds_path:
        # extract ratio from file name
        downsample_ratio = float(downsample_summary.name.split("_")[3][:3])

        downsampling_list.append((downsample_ratio, downsample_summary))

    plot_downsampling(
        downsampling_list, matched_path=args.m_path, figure_path=args.output
    )
