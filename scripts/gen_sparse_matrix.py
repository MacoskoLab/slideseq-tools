#!/usr/bin/python

# This script is to generate mtx files for sparse digital expression matrix

import csv
import logging
import os
import sys
from subprocess import call

import numpy as np
import scipy.io
import scipy.sparse

log = logging.getLogger(__name__)


def main():
    if len(sys.argv) != 6:
        print(
            "Please provide five arguments: manifest file, library ID, locus function list, input folder and file name!"
        )
        sys.exit()

    manifest_file = sys.argv[1]
    library = sys.argv[2]
    locus_function_list = sys.argv[3]
    input_folder = sys.argv[4]
    file_name = sys.argv[5]

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
    scripts_folder = options["scripts_folder"]

    # Read info from metadata file
    reference = ""
    email_address = ""
    with open("{}/parsed_metadata.txt".format(output_folder), "r") as fin:
        reader = csv.reader(fin, delimiter="\t")
        rows = list(reader)
        row0 = rows[0]
        for i in range(1, len(rows)):
            row = rows[i]
            if row[row0.index("library")] == library:
                reference = row[row0.index("reference")]
                email_address = row[row0.index("email")]
                break

    reference_folder = reference[: reference.rfind("/")]
    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]
    reference2 = referencePure + "." + locus_function_list
    annotations_file = "{}/{}.gtf".format(reference_folder, referencePure)

    dge_gzfile = "{}/{}.txt.gz".format(input_folder, file_name)
    dge_file = "{}/{}.txt".format(input_folder, file_name)
    mat_file = "{}/{}_matrix.mtx".format(input_folder, file_name)
    barcodes_file = "{}/{}_barcodes.tsv".format(input_folder, file_name)
    genes_file = "{}/{}_features.tsv".format(input_folder, file_name)

    try:
        log.info(
            f"{flowcell_barcode} - Generating sparse matrix for {library} {reference2} {file_name}"
        )

        os.system("gunzip -c {} > {}".format(dge_gzfile, dge_file))

        # read column names
        with open(dge_file) as f:
            cols = f.readline().split("\t")

        # num of columns
        ncols = len(cols)

        # write barcodes file
        cols = cols[1:]
        with open(barcodes_file, "w") as fout:
            for bc in cols:
                fout.write("{}-1\n".format(bc.strip(" \t\n")))

        # read gene id and name mapping from gtf file
        gene_dict = {}
        with open(annotations_file, "r") as fin:
            for line in fin:
                line = line.strip(" \t\n")
                if len(line) < 1 or line[0] == "#":
                    continue
                items = line.split("\t")[8]
                items = items.split(";")
                gene_id = ""
                name = ""
                for item in items:
                    item = item.strip(" \t\n")
                    if item.split(" ")[0] == "gene_id":
                        gene_id = item.split(" ")[1]
                        gene_id = gene_id.strip(' "')
                    elif item.split(" ")[0] == "gene_name":
                        name = item.split(" ")[1]
                        name = name.strip(' "')
                if name not in gene_dict:
                    gene_dict[name] = gene_id

        # write features (genes) file
        genes = np.loadtxt(dge_file, delimiter="\t", dtype="str", skiprows=1, usecols=0)
        with open(genes_file, "w") as fout:
            for gene in genes:
                if gene in gene_dict:
                    fout.write("{}\t{}\n".format(gene_dict[gene], gene))
                else:
                    fout.write("{}\t{}\n".format(gene, gene))

        # load matrix
        data = np.loadtxt(
            dge_file, delimiter="\t", dtype="i", skiprows=1, usecols=range(1, ncols)
        )

        data = scipy.sparse.csr_matrix(data)

        # write mtx file
        scipy.io.mmwrite(mat_file, data)

        if os.path.isfile(dge_file):
            call(["rm", dge_file])

        call(["gzip", mat_file])

        log.info(
            f"{flowcell_barcode} - Generating sparse matrix {library} {reference2} {file_name} is done"
        )
    except:
        log.exception("EXCEPTION!")

        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = (
                "The Slide-seq workflow for "
                + library
                + " "
                + reference2
                + " "
                + file_name
                + " failed at the step of generating sparse matrix. Please check the log file for the issues. "
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
