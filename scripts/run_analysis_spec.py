#!/usr/bin/python

# This script is to generate digital expression and other analysis outputs

import csv
import logging
import os
import random
import sys
from subprocess import call

from slideseq.util import str2bool

log = logging.getLogger(__name__)


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
    run_barcodematching = False

    bead_type = "180402"

    base_quality = "10"
    min_transcripts_per_cell = "10"
    experiment_date = ""
    gen_downsampling = False
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
                run_barcodematching = str2bool(row[row0.index("run_barcodematching")])

                bead_type = row[row0.index("bead_type")]
                experiment_date = row[row0.index("date")]
                if "gen_downsampling" in row0:
                    gen_downsampling = str2bool(row[row0.index("gen_downsampling")])

    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]

    reference2 = referencePure + "." + locus_function_list

    analysis_folder = f"{library_folder}/{experiment_date}_{library}"
    alignment_folder = f"{analysis_folder}/{reference2}/alignment/"
    barcode_matching_folder = f"{analysis_folder}/{reference2}/barcode_matching/"
    combined_bamfile = f"{analysis_folder}/{library}.bam"

    # Select cells by num transcripts
    commandStr = dropseq_folder + "/SelectCellsByNumTranscripts "
    if is_NovaSeq or is_NovaSeq_S4:
        commandStr += (
            "-m 24076m I="
            + combined_bamfile
            + " MIN_TRANSCRIPTS_PER_CELL="
            + min_transcripts_per_cell
            + " READ_MQ="
            + base_quality
        )
    else:
        commandStr += (
            "-m 7692m I="
            + combined_bamfile
            + " MIN_TRANSCRIPTS_PER_CELL="
            + min_transcripts_per_cell
            + " READ_MQ="
            + base_quality
        )
    commandStr += (
        " OUTPUT="
        + alignment_folder
        + library
        + "."
        + min_transcripts_per_cell
        + "_transcripts_mq_"
        + base_quality
        + "_selected_cells.txt.gz "
    )
    commandStr += "TMP_DIR=" + tmpdir + " VALIDATION_STRINGENCY=SILENT"
    if locus_function_list == "exonic+intronic":
        commandStr += " LOCUS_FUNCTION_LIST=INTRONIC"
    elif locus_function_list == "intronic":
        commandStr += " LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
    log.info(f"{flowcell_barcode} - SelectCellsByNumTranscripts for {library}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - SelectCellsByNumTranscripts for {library} is done")

    # Call run_cmatcher
    if run_barcodematching:
        # wait for BeadBarcodes_degenerate for finish...

        bead_barcode_file = "{}/BeadBarcodes_degenerate.txt".format(analysis_folder)
        select_cell_gzfile = (
            alignment_folder
            + library
            + "."
            + min_transcripts_per_cell
            + "_transcripts_mq_"
            + base_quality
            + "_selected_cells.txt.gz"
        )
        select_cell_file = (
            alignment_folder
            + library
            + "."
            + min_transcripts_per_cell
            + "_transcripts_mq_"
            + base_quality
            + "_selected_cells.txt"
        )
        name = (
            library
            + "."
            + min_transcripts_per_cell
            + "_transcripts_mq_"
            + base_quality
            + "_selected_cells"
        )
        name_shuffled = (
            library
            + "."
            + min_transcripts_per_cell
            + "_transcripts_mq_"
            + base_quality
            + "_selected_cells.shuffled"
        )
        os.system("gunzip -c " + select_cell_gzfile + " > " + select_cell_file)

        select_cell_shuffled_file = (
            alignment_folder
            + library
            + "."
            + min_transcripts_per_cell
            + "_transcripts_mq_"
            + base_quality
            + "_selected_cells.shuffled.txt"
        )
        with open(select_cell_shuffled_file, "w") as fout:
            with open(select_cell_file, "r") as fin:
                for line in fin:
                    line = line.strip(" \t\n")
                    items = list(line)
                    random.shuffle(items)
                    bc = "".join(items)
                    fout.write(bc + "\n")

        with open(select_cell_file, "r") as fin:
            j = sum(1 for _ in fin)

        k = 10000
        ls = j // k

        for i in range(ls + 1):
            if i * k >= j:
                break

            # real barcodes
            infile2 = f"{alignment_folder}/{name}_{str(i + 1)}.txt"
            commandStr = f"awk 'NR >= {str(i * k + 1)} && NR <= {str((i + 1) * k)}' {select_cell_file} > {infile2}"
            os.system(commandStr)

            file4 = f"{barcode_matching_folder}/{library}_barcode_matching_distance_{str(i + 1)}.txt"
            file5 = f"{barcode_matching_folder}/{library}_barcode_matching_{str(i + 1)}.txt"
            output_file = f"{output_folder}/logs/run_cmatcher_{library}_{locus_function_list}_{str(i + 1)}.log"
            submission_script = f"{scripts_folder}/run_cmatcher.sh"
            call_args = [
                "qsub",
                "-o",
                output_file,
                submission_script,
                scripts_folder,
                bead_barcode_file,
                infile2,  # select_cell_file
                file4,  # output_distance_file
                file5,  # output_detail_file
                bead_type,
                output_folder,
                barcode_matching_folder,
            ]
            call(call_args)
            log.info(
                f"{flowcell_barcode} - Run CMatcher for {library} {reference2} {i+1}"
            )

            # shuffled barcodes
            infile2 = f"{alignment_folder}/{name_shuffled}_{str(i + 1)}.txt"
            commandStr = "awk 'NR >= {} && NR <= {}' {} > {}".format(
                str(i * k + 1), str((i + 1) * k), select_cell_shuffled_file, infile2
            )
            os.system(commandStr)

            file4 = f"{barcode_matching_folder}/{library}_barcode_matching_distance_shuffled_{str(i + 1)}.txt"
            file5 = f"{barcode_matching_folder}/{library}_barcode_matching_shuffled_{str(i + 1)}.txt"
            output_file = f"{output_folder}/logs/run_cmatcher_{library}_{locus_function_list}_shuffled_{str(i + 1)}.log"
            submission_script = f"{scripts_folder}/run_cmatcher.sh"
            call_args = [
                "qsub",
                "-o",
                output_file,
                submission_script,
                scripts_folder,
                bead_barcode_file,
                infile2,
                file4,
                file5,
                bead_type,
                output_folder,
                barcode_matching_folder,
            ]
            call(call_args)
            log.info(
                f"{flowcell_barcode} - Run CMatcher for {library} {reference2} {i + 1}"
            )

        # Call run_cmatcher_combine
        output_file = f"{output_folder}/logs/run_cmatcher_combine_{library}_{locus_function_list}.log"
        submission_script = f"{scripts_folder}/run_cmatcher_combine.sh"
        call_args = [
            "qsub",
            "-o",
            output_file,
            submission_script,
            manifest_file,
            library,
            scripts_folder,
            locus_function_list,
            output_folder,
            f"{analysis_folder}/{reference2}",
        ]
        call(call_args)

    # Generate digital expression files for all Illumina barcodes
    commandStr = dropseq_folder + "/DigitalExpression "
    if is_NovaSeq or is_NovaSeq_S4:
        commandStr += "-m 32268m "
    else:
        commandStr += "-m 7692m "
    commandStr += (
        "I="
        + combined_bamfile
        + " O="
        + alignment_folder
        + library
        + ".AllIllumina.digital_expression.txt.gz "
    )
    commandStr += (
        "SUMMARY="
        + alignment_folder
        + library
        + ".AllIllumina.digital_expression_summary.txt EDIT_DISTANCE=1 READ_MQ="
        + base_quality
        + " MIN_BC_READ_THRESHOLD=0 "
    )
    commandStr += (
        "CELL_BC_FILE="
        + alignment_folder
        + library
        + "."
        + min_transcripts_per_cell
        + "_transcripts_mq_"
        + base_quality
        + "_selected_cells.txt.gz TMP_DIR="
        + tmpdir
        + " "
    )
    commandStr += "OUTPUT_HEADER=false UEI=" + library + " VALIDATION_STRINGENCY=SILENT"
    if locus_function_list == "exonic+intronic":
        commandStr += " LOCUS_FUNCTION_LIST=INTRONIC"
    elif locus_function_list == "intronic":
        commandStr += " LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
    log.info(
        f"{flowcell_barcode} - DigitalExpression for {library} for all Illumina barcodes"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - DigitalExpression for {library} for all Illumina barcodes is done."
    )

    if gen_downsampling:
        # Downsample bam
        downsample_folder = "{}/{}_{}/{}/downsample/".format(
            library_folder, experiment_date, library, reference2
        )
        call(["mkdir", "-p", downsample_folder])
        f1 = "{}/{}.AllIllumina.digital_expression_summary.txt".format(
            alignment_folder, library
        )
        f2 = "{}/{}_1.digital_expression_summary.txt".format(downsample_folder, library)
        call(["cp", f1, f2])
        ratio = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        for i in range(0, 9, 1):
            output_file = "{}/logs/gen_downsample_dge_{}_{}_{}.log".format(
                output_folder, library, reference2, str(ratio[i])
            )
            submission_script = "{}/gen_downsample_dge.sh".format(scripts_folder)
            call_args = [
                "qsub",
                "-o",
                output_file,
                submission_script,
                manifest_file,
                library,
                scripts_folder,
                locus_function_list,
                str(ratio[i]),
                output_folder,
                downsample_folder,
            ]
            call(call_args)

        # Call generate_plot_downsampling
        output_file = "{}/logs/generate_plot_downsampling_{}_{}.log".format(
            output_folder, library, reference2
        )
        submission_script = "{}/generate_plot_downsampling.sh".format(scripts_folder)
        call_args = [
            "qsub",
            "-o",
            output_file,
            submission_script,
            manifest_file,
            library,
            scripts_folder,
            locus_function_list,
            output_folder,
            barcode_matching_folder,
        ]
        call(call_args)


if __name__ == "__main__":
    main()
