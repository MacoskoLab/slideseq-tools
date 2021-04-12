#!/usr/bin/python

# This script is to generate barcode_params.txt
# that is needed by extracting Illumina barcodes

import csv
import getopt
import sys


def main(argv):
    inputfile = ""
    outputfile = ""
    lane = ""
    try:
        opts, args = getopt.getopt(argv, "hi:o:l:", ["ifile=", "ofile=", "lane="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-l", "--lane"):
            lane = arg

    with open(inputfile, "r") as fin, open(outputfile, "w") as fout:
        print("barcode_sequence_1\tlibrary_name\tbarcode_name", file=fout)

        reader = csv.DictReader(fin, delimiter="\t")
        for row in reader:
            if row["lane"] != lane:
                continue

            print(
                "\t".join(
                    (
                        row.get("sample_barcode", ""),
                        row.get("library", ""),
                        row.get("barcode_name", ""),
                    )
                ),
                file=fout,
            )


if __name__ == "__main__":
    main(sys.argv[1:])
