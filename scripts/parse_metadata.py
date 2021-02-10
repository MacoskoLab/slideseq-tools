#!/usr/bin/python

# This script is to parse the original metadata file
# and generate the new one with extended lane IDs and missed columns

import csv
import getopt
import sys


# Get lane information from RunInfo.xml
def get_lanes(x):
    lanes = []
    with open(x, "r") as fin:
        for line in fin:
            line = line.strip(" \t\n")
            if line.startswith("<Tile>", 0):
                lane_split = line[6:].split("<")[0]
                lane = lane_split.split("_")[0]
                if lane not in lanes:
                    lanes.append(lane)

    return lanes


def main(argv):
    inputfile = ""
    runinfo_file = ""
    outputfile = ""
    try:
        opts, args = getopt.getopt(argv, "hi:r:o:", ["ifile=", "rfile=", "ofile="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-r", "--rfile"):
            runinfo_file = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    with open(inputfile, "r") as fin, open(outputfile, "w") as fout:
        reader = csv.reader(fin, delimiter="\t")
        idx_LANE = -1
        idx_SAMPLE = -1
        idx_LIBRARY = -1
        idx_SAMPLE_BARCODE = -1
        i = 1
        for row in reader:
            if i == 1:
                row1 = [x.lower() for x in row]
                if "lane" in row1:
                    idx_LANE = row1.index("lane")
                if "library" in row1:
                    idx_LIBRARY = row1.index("library")
                if "sample" in row1:
                    idx_SAMPLE = row1.index("sample")
                if "sample_barcode" in row1:
                    idx_SAMPLE_BARCODE = row1.index("sample_barcode")
                if idx_SAMPLE < 0:
                    row1.append("sample")
                fout.write("\t".join(row1) + "\n")
            else:
                lane = ""
                sample = ""
                sample_barcode = ""

                if idx_LANE >= 0:
                    lane = row[idx_LANE]
                if idx_SAMPLE >= 0:
                    sample = row[idx_SAMPLE]
                elif idx_LIBRARY >= 0:
                    sample = row[idx_LIBRARY]
                if idx_SAMPLE_BARCODE >= 0:
                    sample_barcode = row[idx_SAMPLE_BARCODE]

                if idx_SAMPLE >= 0:
                    row[idx_SAMPLE] = sample
                else:
                    row.append(sample)

                if lane == "{LANE}":
                    lanes = get_lanes(runinfo_file)
                else:
                    lanes = lane.split(",")

                barcodes = sample_barcode.strip(" \t\n").split(",")

                for lane in lanes:
                    if sample_barcode:
                        for b in barcodes:
                            row[idx_LANE] = lane
                            row[idx_SAMPLE_BARCODE] = b
                            fout.write("\t".join(row) + "\n")
                    else:
                        row[idx_LANE] = lane
                        fout.write("\t".join(row) + "\n")

            i = i + 1


if __name__ == "__main__":
    main(sys.argv[1:])
