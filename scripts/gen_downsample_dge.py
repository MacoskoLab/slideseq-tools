#!/usr/bin/python

# This script is to downsample bam file and generate dge

import csv
import os
import sys
import traceback
from datetime import datetime
from subprocess import call


# Convert string to boolean
def str2bool(s):
    return s.lower() == "true"


# Write to log file
def write_log(log_file, flowcell_barcode, log_string):
    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as logfile:
        logfile.write(
            dt_string
            + " [Slide-seq Flowcell Alignment Workflow - "
            + flowcell_barcode
            + "]: "
            + log_string
            + "\n"
        )


def main():
    if len(sys.argv) != 5:
        print(
            "Please provide four arguments: manifest file, library ID, locus function and ratio!"
        )
        sys.exit()

    manifest_file = sys.argv[1]
    library = sys.argv[2]
    locus_function_list = sys.argv[3]
    ratio = sys.argv[4]

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
    scripts_folder = (
        options["scripts_folder"]
        if "scripts_folder" in options
        else "/broad/macosko/jilong/slideseq_pipeline/scripts"
    )
    is_NovaSeq = str2bool(options["is_NovaSeq"]) if "is_NovaSeq" in options else False
    is_NovaSeq_S4 = (
        str2bool(options["is_NovaSeq_S4"]) if "is_NovaSeq_S4" in options else False
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
    min_transcripts_per_cell = "10"
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
                base_quality = row[row0.index("base_quality")]
                min_transcripts_per_cell = row[row0.index("min_transcripts_per_cell")]
                email_address = row[row0.index("email")]
                experiment_date = row[row0.index("date")]

    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]
    reference2 = referencePure + "." + locus_function_list

    analysis_folder = "{}/{}_{}".format(library_folder, experiment_date, library)
    downsample_folder = "{}/{}/downsample/".format(analysis_folder, reference2)
    bam_file = "{}/{}.bam".format(analysis_folder, library)

    folder_running = "{}/status/running.gen_downsample_dge_{}_{}_{}".format(
        output_folder, library, locus_function_list, ratio
    )
    folder_finished = "{}/status/finished.gen_downsample_dge_{}_{}_{}".format(
        output_folder, library, locus_function_list, ratio
    )
    folder_failed = "{}/status/failed.gen_downsample_dge_{}_{}_{}".format(
        output_folder, library, locus_function_list, ratio
    )

    try:
        call(["mkdir", "-p", folder_running])

        # Down sample reads
        sampled_bam_file = downsample_folder + library + "_" + ratio + ".bam"
        commandStr = (
            "java -Djava.io.tmpdir="
            + tmpdir
            + " -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8192m "
        )
        commandStr += (
            "-jar " + picard_folder + "/picard.jar DownsampleSam TMP_DIR=" + tmpdir
        )
        commandStr += " I=" + bam_file + " O=" + sampled_bam_file + " P=" + ratio
        os.system(commandStr)

        # Select cells by num transcripts
        commandStr = dropseq_folder + "/SelectCellsByNumTranscripts "
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += (
                "-m 24076m I="
                + sampled_bam_file
                + " MIN_TRANSCRIPTS_PER_CELL="
                + min_transcripts_per_cell
                + " READ_MQ="
                + base_quality
            )
        else:
            commandStr += (
                "-m 7692m I="
                + sampled_bam_file
                + " MIN_TRANSCRIPTS_PER_CELL="
                + min_transcripts_per_cell
                + " READ_MQ="
                + base_quality
            )
        commandStr += (
            " OUTPUT="
            + downsample_folder
            + library
            + "_"
            + ratio
            + "."
            + min_transcripts_per_cell
            + "_transcripts_mq_"
            + base_quality
            + "_selected_cells.txt.gz TMP_DIR="
            + tmpdir
            + " VALIDATION_STRINGENCY=SILENT"
        )
        if locus_function_list == "exonic+intronic":
            commandStr += " LOCUS_FUNCTION_LIST=INTRONIC"
        elif locus_function_list == "intronic":
            commandStr += " LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
        os.system(commandStr)

        # Generate digital expression files
        commandStr = dropseq_folder + "/DigitalExpression -m 7692m "
        commandStr += (
            "I="
            + sampled_bam_file
            + " O="
            + downsample_folder
            + library
            + "_"
            + ratio
            + ".digital_expression.txt.gz "
        )
        commandStr += (
            "SUMMARY="
            + downsample_folder
            + library
            + "_"
            + ratio
            + ".digital_expression_summary.txt EDIT_DISTANCE=1 READ_MQ="
            + base_quality
            + " MIN_BC_READ_THRESHOLD=0 "
        )
        commandStr += (
            "CELL_BC_FILE="
            + downsample_folder
            + library
            + "_"
            + ratio
            + "."
            + min_transcripts_per_cell
            + "_transcripts_mq_"
            + base_quality
            + "_selected_cells.txt.gz TMP_DIR="
            + tmpdir
            + " "
        )
        commandStr += (
            "OUTPUT_HEADER=true UEI=" + library + " VALIDATION_STRINGENCY=SILENT"
        )
        if locus_function_list == "exonic+intronic":
            commandStr += " LOCUS_FUNCTION_LIST=INTRONIC"
        elif locus_function_list == "intronic":
            commandStr += " LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
        os.system(commandStr)

        if os.path.isfile(
            downsample_folder
            + library
            + "_"
            + ratio
            + "."
            + min_transcripts_per_cell
            + "_transcripts_mq_"
            + base_quality
            + "_selected_cells.txt.gz"
        ):
            os.system(
                "rm "
                + downsample_folder
                + library
                + "_"
                + ratio
                + "."
                + min_transcripts_per_cell
                + "_transcripts_mq_"
                + base_quality
                + "_selected_cells.txt.gz"
            )
        if os.path.isfile(
            downsample_folder
            + library
            + "_"
            + ratio
            + ".SelectCellsByNumTranscripts_metrics"
        ):
            os.system(
                "rm "
                + downsample_folder
                + library
                + "_"
                + ratio
                + ".SelectCellsByNumTranscripts_metrics"
            )
        if os.path.isfile(
            downsample_folder + library + "_" + ratio + ".digital_expression.txt.gz"
        ):
            os.system(
                "rm "
                + downsample_folder
                + library
                + "_"
                + ratio
                + ".digital_expression.txt.gz"
            )
        if os.path.isfile(sampled_bam_file):
            os.system("rm " + sampled_bam_file)

        call(["mv", folder_running, folder_finished])
    except Exception as exp:
        print("EXCEPTION:!")
        print(exp)
        traceback.print_tb(exp.__traceback__, file=sys.stdout)
        if os.path.isdir(folder_running):
            call(["mv", folder_running, folder_failed])
        else:
            call(["mkdir", "-p", folder_failed])

        if len(email_address) > 1:
            subject = "Slide-seq workflow failed for " + flowcell_barcode
            content = (
                "The Slide-seq workflow for "
                + library
                + " "
                + reference2
                + " failed at the step of generating downsampled dge. Please check the log file for the issues. "
            )
            call_args = [
                "python",
                "{}/send_email.py".format(scripts_folder),
                email_address,
                subject,
                content,
            ]
            call(call_args)

        sys.exit()


if __name__ == "__main__":
    main()
