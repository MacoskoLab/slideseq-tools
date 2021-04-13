#!/usr/bin/python

# This script is to generate readme.txt

import csv
import os
import sys


def main():
    if len(sys.argv) != 4:
        print(
            "Please provide three arguments: manifest file, library ID and locus function list!"
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

    # Read info from metadata file
    reference = ""
    base_quality = "10"
    min_transcripts_per_cell = "10"
    experiment_date = ""
    with open("{}/parsed_metadata.txt".format(output_folder), "r") as fin:
        reader = csv.reader(fin, delimiter="\t")
        rows = list(reader)
        row0 = rows[0]
        for i in range(1, len(rows)):
            row = rows[i]
            if row[row0.index("library")] == library:
                reference = row[row0.index("reference")]
                base_quality = row[row0.index("base_quality")]
                min_transcripts_per_cell = row[row0.index("min_transcripts_per_cell")]
                experiment_date = row[row0.index("date")]

    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]
    reference2 = referencePure + "." + locus_function_list

    analysis_folder = f"{library_folder}/{experiment_date}_{library}"
    alignment_folder = f"{analysis_folder}/{reference2}/alignment"
    barcode_matching_folder = f"{analysis_folder}/{reference2}/barcode_matching"

    readme_file = "{}/readme.txt".format(alignment_folder)
    with open(readme_file, "w") as fout:
        print(file=fout)

        f1 = f"{analysis_folder}/BeadBarcodes.txt"
        if os.path.isfile(f1):
            print(f"{f1}\n-bead barcodes from PuckCaller folder\n", file=fout)

        f1 = f"{analysis_folder}/BeadLocations.txt"
        if os.path.isfile(f1):
            print(f"{f1}\n-bead locations from PuckCaller folder\n", file=fout)

        f1 = f"{analysis_folder}/{library}_barcode_matching_01.txt"
        if os.path.isfile(f1):
            print(
                f"{f1}\n-barcode matching results (hamming distance <= 1) within beads for degenerate barcodes\n",
                file=fout,
            )

        f1 = f"{analysis_folder}/{library}_barcode_matching_2.txt"
        if os.path.isfile(f1):
            print(
                f"{f1}\n-barcode matching results (hamming distance >= 2) within beads\n",
                file=fout,
            )

        f1 = f"{analysis_folder}/BeadBarcodes_degenerate.txt"
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write(
                "-the list of degenerate bead barcodes (HD <= 1) and other bead barcodes (HD >= 2)\n"
            )
            fout.write("\n")

        f1 = "{}/{}.bam".format(analysis_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write("-tagged aligned bam from STAR\n")
            fout.write("\n")

        f1 = "{}/{}_alignment_quality.pdf".format(analysis_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write("-mapping quality plots\n")
            fout.write("\n")

        f1 = "{}/{}_mapping_rate.txt".format(analysis_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write("-mapping rate summary\n")
            fout.write("\n")

        f1 = "{}/{}_barcode_matching.txt".format(barcode_matching_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write(
                "-barcode matching results between raw Illumina barcodes and bead barcodes"
                " in BeadBarcodes_degenerate.txt\n"
            )
            fout.write("\n")

        f1 = "{}/{}_barcode_matching_shuffled.txt".format(
            barcode_matching_folder, library
        )
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write(
                "-barcode matching results between shuffled Illumina barcodes and bead barcodes"
                " in BeadBarcodes_degenerate.txt\n"
            )
            fout.write("\n")

        f1 = "{}/{}_matched_bead_barcodes.txt".format(barcode_matching_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write("-the list of matched bead barcodes\n")
            fout.write("\n")

        f1 = "{}/{}_matched_bead_locations.txt".format(barcode_matching_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write(
                "-hamming distance and XY coordinates of matched bead barcodes\n"
            )
            fout.write("\n")

        f1 = "{}/{}_matched.bam".format(barcode_matching_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write("-bam tagged (XC) with matched bead barcodes\n")
            fout.write("\n")

        f1 = "{}/BeadLocationsForR.csv".format(barcode_matching_folder)
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write("-matched bead barcodes and XY coordinates\n")
            fout.write("\n")

        f1 = "{}/{}_XYUMIs.txt".format(barcode_matching_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write("-XY coordinates and the number of UMIs\n")
            fout.write("\n")

        f1 = "{}/{}.{}_transcripts_mq_{}_selected_cells.txt.gz".format(
            alignment_folder, library, min_transcripts_per_cell, base_quality
        )
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write(
                "-selected top cells based on min_transcripts_per_cell and base_quality\n"
            )
            fout.write("\n")

        f1 = "{}/{}.AllIllumina.digital_expression.txt.gz".format(
            alignment_folder, library
        )
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write(
                "-digital expression matrix on all of selected Illumina cell barcodes\n"
            )
            fout.write("\n")

        f1 = "{}/{}.digital_expression_raw.txt.gz".format(alignment_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write(
                "-digital expression matrix on all of matched Illumina cell barcodes\n"
            )
            fout.write("\n")

        f1 = "{}/{}.digital_expression_shuffled.txt.gz".format(
            alignment_folder, library
        )
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write(
                "-digital expression matrix on all of matched shuffled Illumina cell barcodes\n"
            )
            fout.write("\n")

        f1 = "{}/{}.digital_expression.txt.gz".format(alignment_folder, library)
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write("-digital expression matrix on all of matched bead barcodes\n")
            fout.write("\n")

        f1 = "{}/{}_{}.pdf".format(alignment_folder, library, reference2)
        if os.path.isfile(f1):
            fout.write(f1 + "\n")
            fout.write("-plots of alignments and barcode matching\n")
            fout.write("\n")


if __name__ == "__main__":
    main()
