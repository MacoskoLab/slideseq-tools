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
        reader = csv.DictReader(fin, delimiter="\t")
        h = list(reader.fieldnames)
        if "sample" not in h:
            h.append("sample")

        print("\t".join(h), file=fout)

        for row in reader:
            lane = row.get("lane", "")
            if lane == "{LANE}":
                lanes = get_lanes(runinfo_file)
            else:
                lanes = lane.split(",")

            sample = row.get("sample", row.get("library", ""))
            row["sample"] = sample

            sample_barcode = row.get("sample_barcode", "")
            barcodes = sample_barcode.strip().split(",")

            for lane in lanes:
                if sample_barcode:
                    for b in barcodes:
                        row["lane"] = lane
                        row["sample_barcode"] = b
                        print("\t".join(row[c] for c in h), file=fout)
                else:
                    row["lane"] = lane
                    print("\t".join(row[c] for c in h), file=fout)


if __name__ == "__main__":
    main(sys.argv[1:])
