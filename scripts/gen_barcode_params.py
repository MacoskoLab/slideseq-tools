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

    title = "barcode_sequence_1\tlibrary_name\tbarcode_name\n"

    # barcode_sequence_1 => SAMPLE_BARCODE
    # library_name => LIBRARY
    # barcode_name => BARCODE_NAME

    with open(inputfile, "r") as fin, open(outputfile, "w") as fout:
        fout.write(title)

        reader = csv.reader(fin, delimiter="\t")
        idx_LANE = -1
        idx_SAMPLE_BARCODE = -1
        idx_LIBRARY = -1
        idx_BARCODE_NAME = -1
        i = 1
        for row in reader:
            if i == 1:
                if "lane" in row:
                    idx_LANE = row.index("lane")
                if "sample_barcode" in row:
                    idx_SAMPLE_BARCODE = row.index("sample_barcode")
                if "library" in row:
                    idx_LIBRARY = row.index("library")
                if "barcode_name" in row:
                    idx_BARCODE_NAME = row.index("barcode_name")
            else:
                if row[idx_LANE] != lane:
                    continue
                s = ""
                if idx_SAMPLE_BARCODE >= 0:
                    s += row[idx_SAMPLE_BARCODE] + "\t"
                else:
                    s += "\t"
                if idx_LIBRARY >= 0:
                    s += row[idx_LIBRARY] + "\t"
                else:
                    s += "\t"
                if idx_BARCODE_NAME >= 0:
                    s += row[idx_BARCODE_NAME] + "\n"
                else:
                    s += "\n"
                fout.write(s)
            i = i + 1


if __name__ == "__main__":
    main(sys.argv[1:])
