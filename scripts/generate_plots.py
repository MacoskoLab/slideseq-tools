#!/usr/bin/python

# This script is to generate PDFs for the alignment outputs

import csv
import gzip
import logging
import os
import shutil
import sys
from subprocess import call

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from slideseq.util import get_tiles, str2bool

log = logging.getLogger(__name__)


def main():
    if len(sys.argv) != 3:
        print("Please provide two arguments: manifest file and library ID!")
        sys.exit()

    manifest_file = sys.argv[1]
    library = sys.argv[2]

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

    flowcell_directory = options["flowcell_directory"]
    output_folder = options["output_folder"]
    flowcell_barcode = options["flowcell_barcode"]

    library_folder = (
        options["library_folder"]
        if "library_folder" in options
        else "{}/libraries".format(output_folder)
    )
    tmpdir = (
        options["temp_folder"]
        if "temp_folder" in options
        else "{}/tmp".format(output_folder)
    )
    dropseq_folder = (
        options["dropseq_folder"]
        if "dropseq_folder" in options
        else "/broad/macosko/bin/dropseq-tools"
    )
    picard_folder = (
        options["picard_folder"]
        if "picard_folder" in options
        else "/broad/macosko/bin/dropseq-tools/3rdParty/picard"
    )
    is_NovaSeq = str2bool(options["is_NovaSeq"]) if "is_NovaSeq" in options else False
    is_NovaSeq_S4 = (
        str2bool(options["is_NovaSeq_S4"]) if "is_NovaSeq_S4" in options else False
    )
    num_slice_NovaSeq = (
        int(options["num_slice_NovaSeq"]) if "num_slice_NovaSeq" in options else 10
    )
    num_slice_NovaSeq_S4 = (
        int(options["num_slice_NovaSeq_S4"])
        if "num_slice_NovaSeq_S4" in options
        else 40
    )

    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    reference = ""
    base_quality = "10"
    experiment_date = ""
    run_barcodematching = False
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
                base_quality = row[row0.index("base_quality")]
                experiment_date = row[row0.index("date")]
                run_barcodematching = str2bool(row[row0.index("run_barcodematching")])

    reference_folder = reference[: reference.rfind("/")]
    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]
    ref_flat = "{}/{}.refFlat".format(reference_folder, referencePure)
    ribosomal_intervals = "{}/{}.rRNA.intervals".format(reference_folder, referencePure)

    runinfo_file = "{}/RunInfo.xml".format(flowcell_directory)

    # Get tile information from RunInfo.xml
    slice_id = {}
    slice_first_tile = {}
    slice_tile_limit = {}
    for lane in lanes_unique:
        tile_nums = get_tiles(runinfo_file, lane)
        tile_cou = len(tile_nums)
        if (not is_NovaSeq) and (not is_NovaSeq_S4):
            slice_id[lane] = ["0"]
            slice_first_tile[lane] = [str(tile_nums[0])]
            slice_tile_limit[lane] = [str(tile_cou)]
        else:
            slice_cou = num_slice_NovaSeq if is_NovaSeq else num_slice_NovaSeq_S4
            tile_cou_per_slice = (tile_cou // slice_cou) + 1
            slice_id[lane] = []
            slice_first_tile[lane] = []
            slice_tile_limit[lane] = []
            for i in range(slice_cou):
                if tile_cou_per_slice * i >= tile_cou:
                    break
                slice_id[lane].append(str(i))
                slice_first_tile[lane].append(str(tile_nums[tile_cou_per_slice * i]))
                slice_tile_limit[lane].append(str(tile_cou_per_slice))

    alignment_folder = "{}/{}_{}/".format(library_folder, experiment_date, library)
    combined_bamfile = "{}/{}.bam".format(alignment_folder, library)

    # Bam tag histogram
    commandStr = dropseq_folder + "/BamTagHistogram "
    if is_NovaSeq or is_NovaSeq_S4:
        commandStr += "-m 15884m "
    else:
        commandStr += "-m 7692m "
    commandStr += (
        "I="
        + combined_bamfile
        + " OUTPUT="
        + alignment_folder
        + library
        + ".numReads_perCell_XC_mq_"
        + base_quality
        + ".txt.gz "
    )
    commandStr += (
        "TAG=XC FILTER_PCR_DUPLICATES=false TMP_DIR="
        + tmpdir
        + " READ_MQ="
        + base_quality
        + " VALIDATION_STRINGENCY=SILENT"
    )
    log.info(f"{flowcell_barcode} - BamTagHistogram for {library}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - BamTagHistogram for {library} is done.")

    # Collect RnaSeq metrics
    commandStr = (
        "java -Djava.io.tmpdir="
        + tmpdir
        + " -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 "
    )
    if is_NovaSeq or is_NovaSeq_S4:
        commandStr += "-Xmx16384m "
    else:
        commandStr += "-Xmx8192m "
    commandStr += (
        "-jar " + picard_folder + "/picard.jar CollectRnaSeqMetrics TMP_DIR=" + tmpdir
    )
    commandStr += (
        " VALIDATION_STRINGENCY=SILENT I="
        + combined_bamfile
        + " REF_FLAT="
        + ref_flat
        + " STRAND_SPECIFICITY=NONE "
    )
    commandStr += (
        "OUTPUT="
        + alignment_folder
        + library
        + ".fracIntronicExonic.txt RIBOSOMAL_INTERVALS="
        + ribosomal_intervals
    )
    log.info(f"{flowcell_barcode} - CollectRnaSeqMetrics for {library}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - CollectRnaSeqMetrics for {library} is done.")

    # Base distribution at read position for cellular
    commandStr = dropseq_folder + "/BaseDistributionAtReadPosition "
    if is_NovaSeq or is_NovaSeq_S4:
        commandStr += "-m 15884m "
    else:
        commandStr += "-m 7692m "
    commandStr += (
        "I="
        + combined_bamfile
        + " OUTPUT="
        + alignment_folder
        + library
        + ".barcode_distribution_XC.txt TMP_DIR="
        + tmpdir
        + " TAG=XC VALIDATION_STRINGENCY=SILENT"
    )
    log.info(
        f"{flowcell_barcode} - BaseDistributionAtReadPosition Cellular for {library}"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - BaseDistributionAtReadPosition Cellular for {library} is done."
    )

    # Base distribution at read position for molecular
    commandStr = dropseq_folder + "/BaseDistributionAtReadPosition "
    if is_NovaSeq or is_NovaSeq_S4:
        commandStr += "-m 15884m "
    else:
        commandStr += "-m 7692m "
    commandStr += (
        "I="
        + combined_bamfile
        + " OUTPUT="
        + alignment_folder
        + library
        + ".barcode_distribution_XM.txt TMP_DIR="
        + tmpdir
        + " TAG=XM VALIDATION_STRINGENCY=SILENT"
    )
    log.info(
        f"{flowcell_barcode} - BaseDistributionAtReadPosition Molecular for {library}"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - BaseDistributionAtReadPosition Molecular for {library} is done."
    )

    # Gather read quality metrics
    commandStr = dropseq_folder + "/GatherReadQualityMetrics "
    if is_NovaSeq or is_NovaSeq_S4:
        commandStr += "-m 15884m "
    else:
        commandStr += "-m 7692m "
    commandStr += (
        "I="
        + combined_bamfile
        + " TMP_DIR="
        + tmpdir
        + " OUTPUT="
        + alignment_folder
        + library
        + ".ReadQualityMetrics.txt VALIDATION_STRINGENCY=SILENT"
    )
    log.info(f"{flowcell_barcode} - GatherReadQualityMetrics for {library}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - GatherReadQualityMetrics for {library} is done.")

    if not run_barcodematching:
        pp1 = PdfPages(f"{alignment_folder}/{library}.pdf")

        file = f"{alignment_folder}/{library}.ReadQualityMetrics.txt"
        mat = np.loadtxt(
            file,
            delimiter="\t",
            dtype="int",
            skiprows=3,
            max_rows=1,
            usecols=(1, 2, 3, 4),
        )
        df_z = [mat[0], mat[1], mat[2], mat[3]]
        if mat[3] >= 1000000:
            df_u = [
                int(mat[0] / 1000000),
                int(mat[1] / 1000000),
                int(mat[2] / 1000000),
                int(mat[3] / 1000000),
            ]
            yl = "# Reads [millions]"
        else:
            df_u = [mat[0], mat[1], mat[2], mat[3]]
            yl = "# Reads"
        df_y = [
            mat[0] / mat[0] * 100,
            mat[1] / mat[0] * 100,
            mat[2] / mat[0] * 100,
            mat[3] / mat[0] * 100,
        ]
        df_v = [
            "{:,}".format(mat[0]),
            "{:,}".format(mat[1]),
            "{:,}".format(mat[2]),
            "{:,}".format(mat[3]),
        ]
        labels = []
        for i in range(4):
            labels.append("{0:.3g}%".format(df_y[i]))
        df = pd.DataFrame(
            {"x": ["Total", "Mapped", "HQ", "HQ No Dupes"], "z": df_z, "u": df_u}
        )
        fig, ax = plt.subplots(figsize=(8, 8))
        bp = plt.bar(
            df["x"], df["u"], width=0.7, color="lightskyblue", edgecolor="black"
        )
        for idx, rect in enumerate(bp):
            height = rect.get_height()
            ax.text(
                rect.get_x() + rect.get_width() / 2.0,
                0.5 * height,
                labels[idx],
                ha="center",
                va="bottom",
            )
            ax.text(
                rect.get_x() + rect.get_width() / 2.0,
                height,
                df_v[idx],
                ha="center",
                va="bottom",
            )
        plt.yticks(rotation=90)
        plt.ylabel(yl)
        plt.title("Alignment quality for all reads")
        plt.savefig(pp1, format="pdf")

        file = "{}/{}.fracIntronicExonic.txt".format(alignment_folder, library)
        mat = np.loadtxt(
            file,
            delimiter="\t",
            dtype="float",
            skiprows=7,
            max_rows=1,
            usecols=(15, 16, 17, 18, 19),
        )
        mat[1] += mat[2]
        mat[2] = mat[1] + mat[3]
        df_x = ["ribosomal", "exonic", "genic", "intronic", "intergenic"]
        df_y = mat * 100
        labels = []
        for i in range(5):
            labels.append("{0:.3g}".format(df_y[i]))
        df = pd.DataFrame({"x": df_x, "y": df_y})
        fig, ax = plt.subplots(figsize=(8, 8))
        bp = plt.bar(
            df["x"], df["y"], width=0.7, color="lightskyblue", edgecolor="black"
        )
        for idx, rect in enumerate(bp):
            height = rect.get_height()
            ax.text(
                rect.get_x() + rect.get_width() / 2.0,
                0.5 * height,
                labels[idx],
                ha="center",
                va="bottom",
            )
        plt.yticks(rotation=90)
        plt.ylabel("Percentage")
        plt.title("All reads")
        plt.savefig(pp1, format="pdf")

        f1 = "{}/{}.numReads_perCell_XC_mq_{}.txt".format(
            alignment_folder, library, base_quality
        )
        f2 = "{}/{}.numReads_perCell_XC_mq_{}.txt.gz".format(
            alignment_folder, library, base_quality
        )
        if not os.path.isfile(f1):
            with gzip.open(f2, "rb") as f_in:
                with open(f1, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

        mat = np.loadtxt(f1, delimiter="\t", dtype="int", skiprows=1, usecols=0)
        j = int(len(mat) / 10)
        df_x = np.arange(1, j + 1, 1)
        y = np.cumsum(mat)
        df_y = y / max(y)
        df = pd.DataFrame({"x": df_x, "y": df_y[:j]})
        plt.figure(figsize=(8, 8))
        plt.plot(df["x"], df["y"], color="green")
        plt.ylim((0, 1))
        plt.yticks(rotation=90)
        plt.xlabel("cell barcodes sorted by number of reads [descending]")
        plt.ylabel("cumulative fraction of reads")
        plt.title("Cumulative fraction of reads per cell barcode")
        plt.savefig(pp1, format="pdf")

        if os.path.isfile(f1):
            call(["rm", f1])

        file1 = "{}/{}.ReadQualityMetrics.txt".format(alignment_folder, library)
        cou1 = np.loadtxt(
            file1, delimiter="\t", dtype="int", skiprows=3, max_rows=1, usecols=1
        )
        j = 42
        f = False
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for lane_slice in slice_id[lanes[i]]:
                file = "{}/{}.{}.{}.{}.{}.polyA_trimming_report.txt".format(
                    alignment_folder,
                    flowcell_barcode,
                    lanes[i],
                    lane_slice,
                    library,
                    barcodes[i],
                )
                if os.path.isfile(file):
                    mat = np.loadtxt(file, delimiter="\t", dtype="int", skiprows=7)
                    j = len(mat)
                    f = True
                    break
            if f:
                break
        df_x = np.arange(0, j, 1)
        df_y = [0] * j
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue
            for lane_slice in slice_id[lanes[i]]:
                file = "{}/{}.{}.{}.{}.{}.polyA_trimming_report.txt".format(
                    alignment_folder,
                    flowcell_barcode,
                    lanes[i],
                    lane_slice,
                    library,
                    barcodes[i],
                )
                if os.path.isfile(file):
                    mat = np.loadtxt(file, delimiter="\t", dtype="int", skiprows=7)
                    for j in range(0, len(mat)):
                        df_y[mat[j, 0]] += mat[j, 1]
        max_y = max(df_y[1:])
        min_y = min(df_y[1:])
        cou = sum(df_y[1:])
        val = "{0:.3g}".format(cou / cou1 * 100)
        df = pd.DataFrame({"x": df_x, "y": df_y})
        plt.figure(figsize=(8, 8))
        plt.plot(df["x"], df["y"], color="black")
        plt.xlim((1, j))
        plt.ylim((max(0, min_y - 10000), max_y + 10000))
        plt.yticks(rotation=90)
        plt.xlabel("first base of PolyA tail trimmed")
        plt.ylabel("number of reads")
        plt.title("% Reads trimmed by 3' PolyA trimmer: " + val)
        plt.savefig(pp1, format="pdf")

        file = "{}/{}.barcode_distribution_XC.txt".format(alignment_folder, library)
        mat = np.loadtxt(file, delimiter="\t", dtype="int", skiprows=1)
        cou = mat[0, 1] + mat[0, 2] + mat[0, 3] + mat[0, 4]
        df_x = []
        df_x.extend(mat[:, 0])
        df_x.extend(mat[:, 0])
        df_x.extend(mat[:, 0])
        df_x.extend(mat[:, 0])
        df_y = []
        df_y.extend(mat[:, 1] / cou * 100)
        df_y.extend(mat[:, 2] / cou * 100)
        df_y.extend(mat[:, 3] / cou * 100)
        df_y.extend(mat[:, 4] / cou * 100)
        max_y = int(max(df_y))
        j = len(mat)
        fig, ax = plt.subplots(figsize=(8, 8))
        colors = ["red", "blue", "green", "purple"]
        labels = ["A", "C", "G", "T"]
        for i in range(4):
            ax.scatter(
                df_x[i * j : (i + 1) * j],
                df_y[i * j : (i + 1) * j],
                c=colors[i],
                s=20,
                label=labels[i],
            )
        ax.legend(loc="lower right")
        plt.xlim((0, j + 2))
        plt.ylim((0, max_y + 2))
        plt.yticks(rotation=90)
        plt.xlabel("base position")
        plt.ylabel("fraction of reads")
        plt.title("Cell barcodes for all reads")
        plt.savefig(pp1, format="pdf")

        file = "{}/{}.barcode_distribution_XM.txt".format(alignment_folder, library)
        mat = np.loadtxt(file, delimiter="\t", dtype="int", skiprows=1)
        cou = mat[0, 1] + mat[0, 2] + mat[0, 3] + mat[0, 4]
        df_x = []
        df_x.extend(mat[:, 0])
        df_x.extend(mat[:, 0])
        df_x.extend(mat[:, 0])
        df_x.extend(mat[:, 0])
        df_y = []
        df_y.extend(mat[:, 1] / cou * 100)
        df_y.extend(mat[:, 2] / cou * 100)
        df_y.extend(mat[:, 3] / cou * 100)
        df_y.extend(mat[:, 4] / cou * 100)
        max_y = int(max(df_y))
        j = len(mat)
        fig, ax = plt.subplots(figsize=(8, 8))
        colors = ["red", "blue", "green", "purple"]
        labels = ["A", "C", "G", "T"]
        for i in range(4):
            ax.scatter(
                df_x[i * j : (i + 1) * j],
                df_y[i * j : (i + 1) * j],
                c=colors[i],
                s=20,
                label=labels[i],
            )
        ax.legend(loc="lower right")
        plt.xlim((0, j + 1))
        plt.ylim((0, max_y + 2))
        plt.yticks(rotation=90)
        plt.xlabel("base position")
        plt.ylabel("fraction of reads")
        plt.title("Molecular barcodes for all reads")
        plt.savefig(pp1, format="pdf")

        pp1.close()


if __name__ == "__main__":
    main()
