#!/usr/bin/python

# This script is to generate PDF for downsampling and projection

import csv
import logging
import os
import sys

import numpy as np
import pandas as pd
import plotnine

log = logging.getLogger(__name__)


def my_smoother(data, xseq, **params):
    x, y = data["x"], data["y"]
    coefs = np.polynomial.polynomial.polyfit(np.log(x), y, 1)
    ffit = np.polynomial.polynomial.polyval(np.log(xseq), coefs)
    data = pd.DataFrame({"x": xseq, "y": ffit})
    return data


def main():
    if len(sys.argv) != 4:
        print(
            "Please provide three arguments: manifest file, library ID and locus function!"
        )
        sys.exit()

    manifest_file = sys.argv[1]
    library = sys.argv[2]
    locus_function_list = sys.argv[3]

    # Check if the manifest file exists
    if not os.path.isfile(manifest_file):
        print("File {} does not exist. Exiting...".format(manifest_file))
        sys.exit()

    # Read manifest file
    options = {}
    with open(manifest_file, "r") as fp:
        for line in fp:
            key, value = line.rstrip().split("=")
            options[key] = value

    output_folder = options["output_folder"]

    library_folder = (
        options["library_folder"]
        if "library_folder" in options
        else "{}/libraries".format(output_folder)
    )

    # Read info of lane, library and barcode
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    reference = ""
    experiment_date = ""
    with open("{}/parsed_metadata.txt".format(output_folder), "r") as fin:
        reader = csv.reader(fin, delimiter="\t")
        rows = list(reader)
        row0 = rows[0]
        for i in range(1, len(rows)):
            row = rows[i]
            lanes.append(row[row0.index("lane")])
            if row[row0.index("lane")] not in lanes_unique:
                lanes_unique.append(row[row0.index("lane")])
            libraries.append(row[row0.index("library")])
            if row[row0.index("library")] not in libraries_unique:
                libraries_unique.append(row[row0.index("library")])
            barcodes.append(row[row0.index("sample_barcode")])
            bead_structures.append(row[row0.index("bead_structure")])
            if row[row0.index("library")] == library:
                reference = row[row0.index("reference")]
                experiment_date = row[row0.index("date")]

    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]
    reference2 = referencePure + "." + locus_function_list

    downsample_folder = "{}/{}_{}/{}/downsample/".format(
        library_folder, experiment_date, library, reference2
    )
    alignment_folder = "{}/{}_{}/{}/alignment".format(
        library_folder, experiment_date, library, reference2
    )

    ratio = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    df_y = []
    for i in range(0, 10, 1):
        file = (
            downsample_folder
            + library
            + "_"
            + str(ratio[i])
            + ".digital_expression_summary.txt"
        )
        mat = np.loadtxt(file, delimiter="\t", dtype="int", skiprows=7, usecols=2)
        v = int(sum(mat[:10000]) / min(len(mat), 10000))
        df_y.append(v)
    df_x = np.tile(np.arange(0.1, 1.1, 0.1), 1)
    df = pd.DataFrame({"x": df_x, "y": df_y})
    p = (
        plotnine.ggplot(plotnine.aes(x="x", y="y"), df)
        + plotnine.geom_point(size=2.5, color="red", fill="white", show_legend=False)
        + plotnine.geom_smooth(
            method=my_smoother,
            se=False,
            fullrange=True,
            size=0.4,
            color="red",
            fill="white",
            show_legend=False,
        )
        + plotnine.scale_x_continuous(
            limits=(0.1, 10),
            breaks=(0.5, 1.0, 2.0, 5.0, 10.0),
            labels=[0.5, 1.0, 2.0, 5.0, 10.0],
        )
        + plotnine.theme(axis_text_y=plotnine.element_text(rotation=90, hjust=1))
        + plotnine.xlab("Reads generated (relative to this run)")
        + plotnine.ylab("Transcripts per cell")
        + plotnine.ggtitle("Return to sequencing coverage (downsampling + projection)")
    )
    plotnine.ggsave(
        plot=p,
        height=9,
        width=8,
        filename=library + "_" + reference2 + "_downsampling.pdf",
        path=alignment_folder,
        verbose=False,
    )


if __name__ == "__main__":
    main()
