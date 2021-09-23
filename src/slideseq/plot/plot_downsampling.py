#!/usr/bin/python

# edit13 Ali Qutab
# This script is to generate PDF for downsampling
# for each file r = 0.1...1.0, read umi_per_barcode for barcodes that match barcodes in match file
# use Path()
# edit 12 only plotted r = 1.0, indent "data" and "xy.append inside" "if bc in bc_set:"

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

   # right now barcodes is a list
   bc_list, full_umis_per_bc, _ = read_dge_summary(Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04.matched.digital_expression_summary.txt"))
   # make a set out of the list
   bc_set = set(bc_list)

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

   fig = matplotlib.figure.Figure(figsize=(8, 8))
   ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

   ax.scatter(x, y, marker="o", alpha=0.8, color='r')
   ax.set_xlabel("Subsampling Ratio")
   ax.set_ylabel("Transcripts per filtered barcode")
   ax.set_title("Average transcripts for all filtered barcodes")

   ax.set_xlim(0.0, 1.1)

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
   (0.9, Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04_0.9.digital_expression_summary.txt")),
   (1.0, Path("/Users/aqutab/aq/aq_downsampling/aq_files/Puck_210203_04.matched.digital_expression_summary.txt"))],
   figure_path = Path("aq_edit13_plot_downsampling.png"))