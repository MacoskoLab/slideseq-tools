#!/usr/bin/python

# This script is to generate BijectiveMapping.mat

import csv
import logging
import os
import sys
from subprocess import call

log = logging.getLogger(__name__)


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
    flowcell_barcode = options["flowcell_barcode"]

    library_folder = (
        options["library_folder"]
        if "library_folder" in options
        else "{}/libraries".format(output_folder)
    )
    scripts_folder = (
        options["scripts_folder"]
        if "scripts_folder" in options
        else "/broad/macosko/jilong/slideseq_pipeline/scripts"
    )

    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    reference = ""
    puckcaller_path = ""
    email_address = ""
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
                email_address = row[row0.index("email")]
                puckcaller_path = row[row0.index("puckcaller_path")]
                experiment_date = row[row0.index("date")]

    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]
    reference2 = referencePure + "." + locus_function_list

    alignment_folder = "{}/{}_{}/{}/alignment".format(
        library_folder, experiment_date, library, reference2
    )
    barcode_matching_folder = "{}/{}_{}/{}/barcode_matching".format(
        library_folder, experiment_date, library, reference2
    )
    dge_gzfile = "{}/{}.digital_expression.txt.gz".format(alignment_folder, library)
    dge_file = "{}/{}.digital_expression2.txt".format(alignment_folder, library)
    uniqueMappedDge_file = "{}/{}.UniqueMappedDge.txt".format(alignment_folder, library)
    MappedDGEForR_file = "{}/MappedDGEForR.csv".format(alignment_folder)

    try:
        # UniqueMappedIlluminaBarcodes
        unique_bci_file = f"{barcode_matching_folder}/{library}_unique_matched_illumina_barcodes.txt"

        if not os.path.isfile(dge_file):
            os.system(f"gunzip -c {dge_gzfile} > {dge_file}")

        location_file = f"{barcode_matching_folder}/{library}_matched_bead_locations.txt"
        genename_file = f"{barcode_matching_folder}/{library}_genenames.txt"
        bcb_file = f"{barcode_matching_folder}/{library}_unique_matched_beads.txt"

        # split dge_file into three files, plus a csv for R

        with open(dge_file) as fh:
            rdr = csv.reader(fh, delimiter="\t")
            rows = [row for row in rdr if row and row[0][0] != "#"]

            # write bead barcodes (columns of the tsv)
            with open(bcb_file, "w") as out:
                print("\n".join(rows[0][1:]), file=out)

            # write gene names (row labels)
            with open(genename_file, "w") as out:
                print("\n".join(row[0] for row in rows[1:]), file=out)

            # print out the values separately (???)
            with open(uniqueMappedDge_file, "w") as out:
                for row in rows[1:]:
                    print("\t".join(row[1:]), file=out)

            # convert the tsv to csv for R (but R can read tsvs...?)
            with open(MappedDGEForR_file, "w") as out:
                for row in rows:
                    print(",".join(row), file=out)

        # this matlab script

        # Call run_WriteBijectiveMapping
        # 'UniqueMappedBeads','UniqueMappedDGE','UniqueMappedIlluminaBarcodes','GeneNames'
        output_file = f"{output_folder}/logs/run_WriteBijectiveMapping_{library}_{locus_function_list}.log"
        submission_script = "/broad/macosko/jilong/slideseq_pipeline/scripts/run_WriteBijectiveMapping.sh"
        call_args = [
            "qsub",
            "-o",
            output_file,
            submission_script,
            "/broad/software/nonfree/Linux/redhat_7_x86_64/pkgs/matlab_2019a",
            scripts_folder,
            bcb_file,
            uniqueMappedDge_file,
            unique_bci_file,
            genename_file,
            location_file,
            puckcaller_path,
            output_folder,
        ]
        call(call_args)

        os.remove(dge_file)
    except:
        log.exception("EXCEPTION!")

        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = (
                "The Slide-seq workflow for "
                + library
                + " "
                + locus_function_list
                + " failed at the step of generating BijectiveMapping.mat. Please check the log file for the issues. "
            )
            call_args = [
                "python",
                "{}/send_email.py".format(scripts_folder),
                email_address,
                subject,
                content,
            ]
            call(call_args)

        raise


if __name__ == "__main__":
    main()
