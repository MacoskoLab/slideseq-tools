#!/usr/bin/python

# edit40 Ali Qutab
# plot real data and model data
# this script edit plots five quantiles by fitting model for the top 20%, 40%, 60%, 80%, 100%
# conceptually have five different values of data and fit the model five times, plot five lines and five sets of points
# all edits upto edit 34, code is reading specifically from the file on my computer, instead;
# trying to use arguments for ratio and path instead of hardcoding them into the script
# the downsampled files are all going to be [Puck name]_[ratio].digital_expression_summary.txt
# the matched file will be [Puck name].matched.digital_expression_summary.txt
# use argparse to get a list of files from the command line, as strings
# split each file and the first 3 characters to get the ratio value, and convert that with float
# make the list of (float, Path) to pass into the main plotting function

import logging
from pathlib import Path

import matplotlib.figure
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg

from slideseq.plot import read_dge_summary

import scipy.optimize

import argparse
import os
import glob
import pathlib

log = logging.getLogger(__name__)


# pycharmedit

def plot_downsampling(downsampling_output: list[tuple[float, Path]], figure_path: Path):
   x_100y_100 = []
   x_80y_80 = []
   x_60y_60 = []
   x_40y_40 = []
   x_20y_20 = []

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
       data_80 = np.mean(filtered_umis_per_bc[:(int(0.8 * (len(filtered_umis_per_bc))))])  # for example, top 80% is from top item to item at 0.8 * number of items in list filtered_umis_per_bc
       data_60 = np.mean(filtered_umis_per_bc[:(int(0.6 * (len(filtered_umis_per_bc))))])
       data_40 = np.mean(filtered_umis_per_bc[:(int(0.4 * (len(filtered_umis_per_bc))))])
       data_20 = np.mean(filtered_umis_per_bc[:(int(0.2 * (len(filtered_umis_per_bc))))])

       x_100y_100.append((r, data_100))
       x_80y_80.append((r, data_80))
       x_60y_60.append((r, data_60))
       x_40y_40.append((r, data_40))
       x_20y_20.append((r, data_20))

   x_100y_100.sort()
   x_100, y_100 = zip(*x_100y_100)
   x_100 = np.array(x_100)
   y_100 = np.array(y_100)

   x_80y_80.sort()
   x_80, y_80 = zip(*x_80y_80)
   x_80 = np.array(x_80)
   y_80 = np.array(y_80)

   x_60y_60.sort()
   x_60, y_60 = zip(*x_60y_60)
   x_60 = np.array(x_60)
   y_60 = np.array(y_60)

   x_40y_40.sort()
   x_40, y_40 = zip(*x_40y_40)
   x_40 = np.array(x_40)
   y_40 = np.array(y_40)

   x_20y_20.sort()
   x_20, y_20 = zip(*x_20y_20)
   x_20 = np.array(x_20)
   y_20 = np.array(y_20)

   def model_100(r, params_100):
       """
       y = alpha * exp(-r) + beta
       y is average transcript count for the 0.1,0.2,... files
       r = 0.1, 0.2,...,1.0 - for each file of data given
       params is an array of two values [alpha, beta] which we will optimize
       """
       return params_100[0] * np.exp(-r) + params_100[1]

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

   output_100 = scipy.optimize.least_squares(
       model_least_squares_100,
       [-10.0, 10.],  # initial values for params
       bounds=([-np.inf, 0], [0, np.inf]),  # some bounds on params. alpha should be negative, beta is positive
       kwargs={"data_100": y_100, "r": x_100},
       method="dogbox",  # I found this method to work well for this problem
   )

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

   def model_least_squares_80(params_80, *, r, data_80):
       # computes model and compares it to the data
       return model_80(r, params_80) - data_80

   output_80 = scipy.optimize.least_squares(
       model_least_squares_80,
       [-10.0, 10.],  # initial values for params
       bounds=([-np.inf, 0], [0, np.inf]),  # some bounds on params. alpha should be negative, beta is positive
       kwargs={"data_80": y_80, "r": x_80},
       method="dogbox",  # I found this method to work well for this problem
   )

   def model_60(r, params_60):
       """
       y = alpha * exp(-r) + beta
       y is average transcript count for the 0.1,0.2,... files
       r = 0.1, 0.2,...,1.0 - for each file of data given
       params is an array of two values [alpha, beta] which we will optimize
       """
       return params_60[0] * np.exp(-r) + params_60[1]

   """
   model function takes 0.1, 0.2 ... 1.0 for r and outputs the mean counts values which are y,
   based on some parameters alpha and beta, which we need to find to predict the result for any value of r.
   """

   """
   The least_squares function is to optimize the model function, we want to minimize the error, which is the difference between the model prediction and the data.
   """

   def model_least_squares_60(params_60, *, r, data_60):
       # computes model and compares it to the data
       return model_60(r, params_60) - data_60

   output_60 = scipy.optimize.least_squares(
       model_least_squares_60,
       [-10.0, 10.],  # initial values for params
       bounds=([-np.inf, 0], [0, np.inf]),  # some bounds on params. alpha should be negative, beta is positive
       kwargs={"data_60": y_60, "r": x_60},
       method="dogbox",  # I found this method to work well for this problem
   )

   def model_40(r, params_40):
       """
       y = alpha * exp(-r) + beta
       y is average transcript count for the 0.1,0.2,... files
       r = 0.1, 0.2,...,1.0 - for each file of data given
       params is an array of two values [alpha, beta] which we will optimize
       """
       return params_40[0] * np.exp(-r) + params_40[1]

   """
   model function takes 0.1, 0.2 ... 1.0 for r and outputs the mean counts values which are y,
   based on some parameters alpha and beta, which we need to find to predict the result for any value of r.
   """

   """
   The least_squares function is to optimize the model function, we want to minimize the error, which is the difference between the model prediction and the data.
   """

   def model_least_squares_40(params_40, *, r, data_40):
       # computes model and compares it to the data
       return model_40(r, params_40) - data_40

   output_40 = scipy.optimize.least_squares(
       model_least_squares_40,
       [-10.0, 10.],  # initial values for params
       bounds=([-np.inf, 0], [0, np.inf]),  # some bounds on params. alpha should be negative, beta is positive
       kwargs={"data_40": y_40, "r": x_40},
       method="dogbox",  # I found this method to work well for this problem
   )

   def model_20(r, params_20):
       """
       y = alpha * exp(-r) + beta
       y is average transcript count for the 0.1,0.2,... files
       r = 0.1, 0.2,...,1.0 - for each file of data given
       params is an array of two values [alpha, beta] which we will optimize
       """
       return params_20[0] * np.exp(-r) + params_20[1]

   """
   model function takes 0.1, 0.2 ... 1.0 for r and outputs the mean counts values which are y,
   based on some parameters alpha and beta, which we need to find to predict the result for any value of r.
   """

   """
   The least_squares function is to optimize the model function, we want to minimize the error, which is the difference between the model prediction and the data.
   """

   def model_least_squares_20(params_20, *, r, data_20):
       # computes model and compares it to the data
       return model_20(r, params_20) - data_20

   output_20 = scipy.optimize.least_squares(
       model_least_squares_20,
       [-10.0, 10.],  # initial values for params
       bounds=([-np.inf, 0], [0, np.inf]),  # some bounds on params. alpha should be negative, beta is positive
       kwargs={"data_20": y_20, "r": x_20},
       method="dogbox",  # I found this method to work well for this problem
   )

   params_100 = output_100.x
   params_80 = output_80.x
   params_60 = output_60.x
   params_40 = output_40.x
   params_20 = output_20.x

   x_values = np.linspace(0.1, 3.0,30)  # this function creates a linear space of points: 30 points from 0.1 to 3.0 (0.1, 0.2, ... 2.9, 3.0)
   predicted_y_100 = model_100(x_values, params_100)
   predicted_y_80 = model_80(x_values, params_80)
   predicted_y_60 = model_80(x_values, params_60)
   predicted_y_40 = model_80(x_values, params_40)
   predicted_y_20 = model_80(x_values, params_20)

   fig = matplotlib.figure.Figure(figsize=(8, 8))
   ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

   ax.scatter(x_100, y_100, marker="o", alpha=0.8, color="r")  # red scatter plot for actual data r = 0.1...1.0
   ax.plot(x_values, predicted_y_100, alpha=0.8, color="r",
           label='top 100%')  # blue line plot for model data r = 0.1...3.0

   ax.scatter(x_80, y_80, marker="o", alpha=0.8, color="g")
   ax.plot(x_values, predicted_y_80, alpha=0.8, color="g", label='top 80%')

   ax.scatter(x_60, y_60, marker="o", alpha=0.8, color="m")
   ax.plot(x_values, predicted_y_60, alpha=0.8, color="m", label='top 60%')

   ax.scatter(x_40, y_40, marker="o", alpha=0.8, color="k")
   ax.plot(x_values, predicted_y_40, alpha=0.8, color="k", label='top 40%')

   ax.scatter(x_20, y_20, marker="o", alpha=0.8, color="c")
   ax.plot(x_values, predicted_y_20, alpha=0.8, color="c", label='top 20%')

   ax.set_xlabel("Subsampling Ratio")
   ax.set_ylabel("Transcripts per matched barcode")
   ax.set_title("Average transcripts for matched barcodes")

   ax.set_xlim(0.0, 3.1)  # r went upto 1.0 for actual data, but x_values for model go up to 3.0

   ax.legend()

   FigureCanvasAgg(fig).print_figure(figure_path)

if __name__ == "__main__":
   # use argparse to get a list of files from the command line, as strings
   parser = argparse.ArgumentParser(description='Read in downsample_summary text files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('path', type=pathlib.Path, nargs='+', help='Path of a downsample_summary text file')
   args = parser.parse_args()
   print(args)

   # Parse paths
   for downsample_summary in args.path:
       # trying to read each file and convert to string
       read_downsample_summary = open(downsample_summary, 'r')
       string_downsample_summary = str(read_downsample_summary)
       split_path = string_downsample_summary.split("_")  # split each file and the first 3 characters to get the ratio value
       characters = (str(split_path[6:7]) + "\n")
       print(characters[2:5])
       # r = float(characters[2:5])

   # make the list of (float, Path) to pass into the main plotting function
   downsampling_list = []
   downsampling_list.append((r, downsample_summary))

   # trying to use arguments for ratio and path instead of hardcoding them into the script
   plot_downsampling(downsampling_list, figure_path=Path("/Users/aqutab/aq/aq_downsampling/aq_plots/aq_edit40_plot_downsampling.png"))

