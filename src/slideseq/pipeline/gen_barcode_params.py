#!/usr/bin/python

# This script is to generate barcode_params.txt
# that is needed by extracting Illumina barcodes

import csv


def gene_barcode_params(input_file, output_file, lane: str):
    with open(input_file, "r") as fin, open(output_file, "w") as fout:
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
