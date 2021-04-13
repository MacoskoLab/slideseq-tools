#!/usr/bin/python

# This script is to filter unmapped bam using matched barcodes

import csv
import logging
import os
import re
import sys
from subprocess import call

log = logging.getLogger(__name__)


# Get bead structure range
def get_bead_structure_range(bs, structure_type):
    # 12C8M|*T
    # 7C18X7C8M2X|*T
    ell = re.split("[CXM]", bs.split("|")[0])
    res = ""
    i = 1
    p = -1
    for it in ell:
        if it:
            p += len(it) + 1
            if bs[p] == structure_type:
                res += str(i) + "-" + str(i + int(it) - 1) + ":"
            i += int(it)
    return res[:-1]


def main():
    if len(sys.argv) != 7:
        print(
            "Please provide six arguments: "
            "manifest file, library ID, lane ID, slice ID, sample barcode and locus function list!"
        )
        sys.exit()

    manifest_file = sys.argv[1]
    library = sys.argv[2]
    lane = sys.argv[3]
    slice_id = sys.argv[4]
    barcode = sys.argv[5]
    locus_function_list = sys.argv[6]

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
    reference = ""
    email_address = ""
    bead_structure = ""
    experiment_date = ""
    with open("{}/parsed_metadata.txt".format(output_folder), "r") as fin:
        reader = csv.reader(fin, delimiter="\t")
        rows = list(reader)
        row0 = rows[0]
        for i in range(1, len(rows)):
            row = rows[i]
            if row[row0.index("library")] == library:
                reference = row[row0.index("reference")]
                email_address = row[row0.index("email")]
                bead_structure = row[row0.index("bead_structure")]
                experiment_date = row[row0.index("date")]
                break

    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]
    reference2 = referencePure + "." + locus_function_list

    prefix_libraries = "{}/{}_{}/{}.{}.{}.{}".format(
        library_folder,
        experiment_date,
        library,
        flowcell_barcode,
        lane,
        slice_id,
        library,
    )
    if barcode:
        prefix_libraries += "." + barcode
    unmapped_bam = prefix_libraries + ".unmapped.bam"

    barcode_matching_folder = "{}/{}_{}/{}/barcode_matching".format(
        library_folder, experiment_date, library, reference2
    )
    unmapped_sam = "{}/{}_{}_{}_{}_unmapped.sam".format(
        barcode_matching_folder, library, lane, slice_id, barcode
    )
    filtered_sam = "{}/{}_{}_{}_{}_filtered.sam".format(
        barcode_matching_folder, library, lane, slice_id, barcode
    )
    filtered_bam = "{}/{}_{}_{}_{}_filtered.bam".format(
        barcode_matching_folder, library, lane, slice_id, barcode
    )
    combined_cmatcher_file = "{}/{}_barcode_matching.txt".format(
        barcode_matching_folder, library
    )

    if not os.path.isfile(unmapped_bam):
        log.error(
            f"{flowcell_barcode} - TagMatchedBam error: {unmapped_sam} does not exist!",
        )
        raise FileNotFoundError(f"TagMatchedBam error: {unmapped_bam} does not exist!")

    if not os.path.isfile(combined_cmatcher_file):
        log.error(
            f"{flowcell_barcode} - TagMatchedBam error: {combined_cmatcher_file} does not exist!",
        )
        raise FileNotFoundError(
            f"TagMatchedBam error: {combined_cmatcher_file} does not exist!"
        )

    try:
        log.info(
            f"{flowcell_barcode} - Filter unmapped bam for {library} {reference2} in lane {lane} slice {slice_id}"
        )

        commandStr = f"samtools view -h -o {unmapped_sam} {unmapped_bam}"
        os.system(commandStr)

        dict1 = {}
        with open(combined_cmatcher_file, "r") as fin:
            j = 0
            for line in fin:
                j += 1
                if j > 1:
                    dict1[line.split("\t")[0]] = line.split("\t")[2]

        bs_range1 = get_bead_structure_range(bead_structure, "C")  # 1-7:26-32
        lists = re.split(":", bs_range1)
        with open(filtered_sam, "w") as fout:
            with open(unmapped_sam, "r") as fin:
                flag = False
                for line in fin:
                    if line[0] == "@":
                        fout.write(line)
                    else:
                        items1 = line.split("\t")
                        if items1[1] == "77":
                            read1 = items1[9]
                            r = ""
                            for s in lists:
                                v = re.split("-", s)
                                r += read1[int(v[0]) - 1 : int(v[1])]
                            if r in dict1:
                                fout.write(line)
                                flag = True
                            else:
                                flag = False
                        elif items1[1] == "141":
                            if flag:
                                fout.write(line)

        if os.path.isfile(unmapped_sam):
            call(["rm", unmapped_sam])

        commandStr = "samtools view -S -b " + filtered_sam + " > " + filtered_bam
        os.system(commandStr)

        if os.path.isfile(filtered_sam):
            call(["rm", filtered_sam])

        log.info(
            f"{flowcell_barcode} - Filter unmapped bam for {library} {reference2} in lane {lane} slice {slice_id} done"
        )
    except:
        log.exception("EXCEPTION!")

        if len(email_address) > 1:
            subject = f"Slide-seq workflow failed for {flowcell_barcode}"
            content = (
                f"The Slide-seq workflow for {library} {reference2} in lane {lane} slice {slice_id}"
                " failed at the step of filtering unmapped bam. Please check the log file for the issues."
            )
            call_args = [
                "python",
                f"{scripts_folder}/send_email.py",
                email_address,
                subject,
                content,
            ]
            call(call_args)

        raise


if __name__ == "__main__":
    main()
