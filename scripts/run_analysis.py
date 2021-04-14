#!/usr/bin/python

# This script is to combine bam files from slice alignments

import csv
import logging
import os
import sys
import time
from subprocess import call

from slideseq.util import get_tiles, str2bool

log = logging.getLogger(__name__)


# Number of input reads |	3979435
def get_key(line):
    line = line.strip()
    res = line.split("|")[0]
    res = res.strip()
    return res


def get_val(line):
    line = line.strip()
    res = line.split("|")[1]
    res = res.strip()
    return res


def main():
    if len(sys.argv) != 3:
        print("Please provide two arguments: manifest file and library ID!")
        sys.exit()

    manifest_file = sys.argv[1]
    library = sys.argv[2]

    # Check if the manifest file exists
    if not os.path.isfile(manifest_file):
        print(f"File {manifest_file} does not exist. Exiting...")
        sys.exit(1)

    # Read manifest file
    options = {}
    with open(manifest_file, "r") as fp:
        for line in fp:
            key, value = line.rstrip().split("=")
            options[key] = value

    flowcell_directory = options["flowcell_directory"]
    output_folder = options["output_folder"]

    flowcell_barcode = options["flowcell_barcode"]

    library_folder = (
        options["library_folder"]
        if "library_folder" in options
        else f"{output_folder}/libraries"
    )
    tmpdir = (
        options["temp_folder"] if "temp_folder" in options else f"{output_folder}/tmp"
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
    num_slice_NovaSeq = (
        int(options["num_slice_NovaSeq"]) if "num_slice_NovaSeq" in options else 10
    )
    num_slice_NovaSeq_S4 = (
        int(options["num_slice_NovaSeq_S4"])
        if "num_slice_NovaSeq_S4" in options
        else 40
    )

    runinfo_file = f"{flowcell_directory}/RunInfo.xml"

    # Read info from metadata file
    lanes = []
    lanes_unique = []
    libraries = []
    libraries_unique = []
    barcodes = []
    bead_structures = []
    reference = ""
    locus_function_list = "exonic+intronic"
    run_barcodematching = False
    puckcaller_path = ""
    bead_type = "180402"
    experiment_date = ""
    gen_updistance_plot = False
    with open(f"{output_folder}/parsed_metadata.txt", "r") as fin:
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
                locus_function_list = row[row0.index("locus_function_list")]
                run_barcodematching = str2bool(row[row0.index("run_barcodematching")])
                puckcaller_path = row[row0.index("puckcaller_path")]
                bead_type = row[row0.index("bead_type")]
                experiment_date = row[row0.index("date")]
                if "gen_updistance_plot" in row0:
                    gen_updistance_plot = str2bool(
                        row[row0.index("gen_updistance_plot")]
                    )

    # Get tile information from RunInfo.xml
    slice_id = {}
    slice_first_tile = {}
    slice_tile_limit = {}
    for lane in lanes_unique:
        tile_nums = get_tiles(runinfo_file, lane)
        tile_cou = len(tile_nums)
        if (not is_NovaSeq) and (not is_NovaSeq_S4):
            slice_id[lane] = ["0"]
            slice_first_tile[lane] = [str(tile_nums[0])]
            slice_tile_limit[lane] = [str(tile_cou)]
        else:
            slice_cou = num_slice_NovaSeq if is_NovaSeq else num_slice_NovaSeq_S4
            tile_cou_per_slice = (tile_cou // slice_cou) + 1
            slice_id[lane] = []
            slice_first_tile[lane] = []
            slice_tile_limit[lane] = []
            for i in range(slice_cou):
                if tile_cou_per_slice * i >= tile_cou:
                    break
                slice_id[lane].append(str(i))
                slice_first_tile[lane].append(str(tile_nums[tile_cou_per_slice * i]))
                slice_tile_limit[lane].append(str(tile_cou_per_slice))

    analysis_folder = f"{library_folder}/{experiment_date}_{library}"

    if run_barcodematching:
        file2 = f"{puckcaller_path}/BeadBarcodes.txt"
        file3 = f"{puckcaller_path}/BeadLocations.txt"
        while 1:
            if os.path.isfile(file2) and os.path.isfile(file3):
                break
            time.sleep(30)

        call(["cp", file2, analysis_folder + "/"])
        call(["cp", file3, analysis_folder + "/"])
        bead_barcode_file = f"{analysis_folder}/BeadBarcodes.txt"
        bead_location_file = f"{analysis_folder}/BeadLocations.txt"

        with open(bead_barcode_file, "r") as fin:
            j = sum(1 for _ in fin)

        k = 10000
        ls = j // k

        for i in range(ls + 1):
            if i * k >= j:
                break

            infile2 = f"{analysis_folder}/BeadBarcodes_{str(i + 1)}.txt"
            commandStr = "awk 'NR >= {} && NR <= {}' {} > {}".format(
                str(i * k + 1), str((i + 1) * k), bead_barcode_file, infile2
            )
            os.system(commandStr)

            file4 = f"{analysis_folder}/{library}_barcode_matching_01_{str(i + 1)}.txt"
            file5 = f"{analysis_folder}/{library}_barcode_matching_2_{str(i + 1)}.txt"
            output_file = f"{output_folder}/logs/run_cmatcher_beads_{str(i + 1)}.log"
            submission_script = f"{scripts_folder}/run_cmatcher_beads.sh"
            call_args = [
                "qsub",
                "-o",
                output_file,
                submission_script,
                scripts_folder,
                infile2,
                bead_barcode_file,
                bead_location_file,
                file4,
                file5,
                bead_type,
                output_folder,
                analysis_folder,
            ]
            call(call_args)

        # Call run_cmatcher_beads_combine
        output_file = f"{output_folder}/logs/run_cmatcher_beads_combine_{library}.log"
        submission_script = f"{scripts_folder}/run_cmatcher_beads_combine.sh"
        call_args = [
            "qsub",
            "-o",
            output_file,
            submission_script,
            manifest_file,
            library,
            scripts_folder,
            output_folder,
            analysis_folder,
        ]
        call(call_args)

    # Wait for all of run_alignment finish

    # Merge bam files
    combined_bamfile = "{}/{}.bam".format(analysis_folder, library)
    commandStr = (
        f"java -Djava.io.tmpdir={tmpdir}"
        " -Dsamjdk.buffer_size=131072 -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50"
        " -XX:GCHeapFreeLimit=10 -Xmx8192m "
    )
    commandStr += (
        "-jar "
        + picard_folder
        + "/picard.jar MergeSamFiles TMP_DIR="
        + tmpdir
        + " CREATE_INDEX=true CREATE_MD5_FILE=false VALIDATION_STRINGENCY=SILENT "
    )
    commandStr += (
        "OUTPUT=" + combined_bamfile + " SORT_ORDER=coordinate ASSUME_SORTED=true"
    )
    for i in range(len(lanes)):
        if libraries[i] != library:
            continue
        for lane_slice in slice_id[lanes[i]]:
            star_bamfile = "{}/{}.{}.{}.{}".format(
                analysis_folder, flowcell_barcode, lanes[i], lane_slice, library
            )
            if barcodes[i]:
                star_bamfile += "." + barcodes[i]
            star_bamfile += ".star_gene_exon_tagged2.bam"
            if not os.path.isfile(star_bamfile):
                log.error(
                    f"{flowcell_barcode} - MergeSamFiles error: {star_bamfile} does not exist!",
                )
                raise Exception(f"MergeSamFiles error: {star_bamfile} does not exist!")
            commandStr += " INPUT=" + star_bamfile
    log.info(f"{flowcell_barcode} - MergeSamFiles for {library}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - MergeSamFiles for {library} is done. ")

    # Validate bam file
    commandStr = (
        f"java -Djava.io.tmpdir={tmpdir}"
        " -Dsamjdk.buffer_size=131072 -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50"
        " -XX:GCHeapFreeLimit=10 -Xmx16384m "
    )
    commandStr += f"-jar {picard_folder}/picard.jar ValidateSamFile TMP_DIR={tmpdir} VALIDATION_STRINGENCY=SILENT "
    commandStr += "INPUT=" + combined_bamfile + " MODE=SUMMARY"
    if (not is_NovaSeq) and (not is_NovaSeq_S4):
        commandStr += " IGNORE=MISSING_PLATFORM_VALUE IGNORE=INVALID_VERSION_NUMBER"
    log.info(f"{flowcell_barcode} - ValidateSamFile for " + library)
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - ValidateSamFile for {library} is done.")

    # Call generate_plots
    output_file = "{}/logs/generate_plots_{}.log".format(output_folder, library)
    submission_script = "{}/generate_plots.sh".format(scripts_folder)
    call_args = [
        "qsub",
        "-o",
        output_file,
        submission_script,
        manifest_file,
        library,
        scripts_folder,
        output_folder,
        analysis_folder,
    ]
    call(call_args)

    lists = locus_function_list.split(",")
    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]
    for j in lists:
        call(["mkdir", "-p", "{}/{}.{}".format(analysis_folder, referencePure, j)])
        call(
            [
                "mkdir",
                "-p",
                "{}/{}.{}/alignment".format(analysis_folder, referencePure, j),
            ]
        )

        if run_barcodematching:
            barcode_matching_folder = "{}/{}.{}/barcode_matching/".format(
                analysis_folder, referencePure, j
            )
            call(["mkdir", "-p", barcode_matching_folder])
            for i in range(len(lanes)):
                if libraries[i] != library:
                    continue
                for lane_slice in slice_id[lanes[i]]:
                    toCopyFile = "{}/{}.{}.{}.{}".format(
                        analysis_folder,
                        flowcell_barcode,
                        lanes[i],
                        lane_slice,
                        library,
                    )
                    if barcodes[i]:
                        toCopyFile += "." + barcodes[i]
                    toCopyFile += ".star_gene_exon_tagged2.bam"
                    if os.path.isfile(toCopyFile):
                        call(["cp", toCopyFile, barcode_matching_folder])

        # Call run_analysis_spec
        output_file = "{}/logs/run_analysis_spec_{}_{}.log".format(
            output_folder, library, j
        )
        submission_script = "{}/run_analysis_spec.sh".format(scripts_folder)
        call_args = [
            "qsub",
            "-o",
            output_file,
            submission_script,
            manifest_file,
            library,
            scripts_folder,
            j,
            output_folder,
            "{}/{}.{}".format(analysis_folder, referencePure, j),
        ]
        call(call_args)

    for i in range(len(lanes)):
        if libraries[i] != library:
            continue
        for lane_slice in slice_id[lanes[i]]:
            toDeleteFile = "{}/{}.{}.{}.{}".format(
                analysis_folder, flowcell_barcode, lanes[i], lane_slice, library
            )
            if barcodes[i]:
                toDeleteFile += "." + barcodes[i]
            toDeleteFile += ".star_gene_exon_tagged2.bam"
            if os.path.isfile(toDeleteFile):
                call(["rm", toDeleteFile])

    # Combine check_alignments_quality files
    dict_unique_score = {}
    dict_multi_score = {}
    dict_unique_mismatch = {}
    dict_multi_mismatch = {}
    dict_unique_ratio = {}
    dict_multi_ratio = {}
    for i in range(len(lanes)):
        if libraries[i] != library:
            continue
        for lane_slice in slice_id[lanes[i]]:
            star_samfile = "{}/{}.{}.{}.{}".format(
                analysis_folder, flowcell_barcode, lanes[i], lane_slice, library
            )
            if barcodes[i]:
                star_samfile += "." + barcodes[i]
            star_samfile += ".star.Aligned.out.sam"
            file1 = star_samfile + ".unique.score"
            file2 = star_samfile + ".multi.score"
            file3 = star_samfile + ".unique.mismatch"
            file4 = star_samfile + ".multi.mismatch"
            file5 = star_samfile + ".unique.ratio"
            file6 = star_samfile + ".multi.ratio"
            if os.path.isfile(file1):
                with open(file1, "r") as fin:
                    for line in fin:
                        c1 = line.split("\t")[0].strip(" \t\n")
                        c2 = int(line.split("\t")[1].strip(" \t\n"))
                        if c1 not in dict_unique_score:
                            dict_unique_score[c1] = c2
                        else:
                            dict_unique_score[c1] += c2

            if os.path.isfile(file2):
                with open(file2, "r") as fin:
                    for line in fin:
                        c1 = line.split("\t")[0].strip(" \t\n")
                        c2 = int(line.split("\t")[1].strip(" \t\n"))
                        if c1 not in dict_multi_score:
                            dict_multi_score[c1] = c2
                        else:
                            dict_multi_score[c1] += c2

            if os.path.isfile(file3):
                with open(file3, "r") as fin:
                    for line in fin:
                        c1 = line.split("\t")[0].strip(" \t\n")
                        c2 = int(line.split("\t")[1].strip(" \t\n"))
                        if c1 not in dict_unique_mismatch:
                            dict_unique_mismatch[c1] = c2
                        else:
                            dict_unique_mismatch[c1] += c2

            if os.path.isfile(file4):
                with open(file4, "r") as fin:
                    for line in fin:
                        c1 = line.split("\t")[0].strip(" \t\n")
                        c2 = int(line.split("\t")[1].strip(" \t\n"))
                        if c1 not in dict_multi_mismatch:
                            dict_multi_mismatch[c1] = c2
                        else:
                            dict_multi_mismatch[c1] += c2

            if os.path.isfile(file5):
                with open(file5, "r") as fin:
                    for line in fin:
                        c1 = line.split("\t")[0].strip(" \t\n")
                        c2 = int(line.split("\t")[1].strip(" \t\n"))
                        if c1 not in dict_unique_ratio:
                            dict_unique_ratio[c1] = c2
                        else:
                            dict_unique_ratio[c1] += c2

            if os.path.isfile(file6):
                with open(file6, "r") as fin:
                    for line in fin:
                        c1 = line.split("\t")[0].strip(" \t\n")
                        c2 = int(line.split("\t")[1].strip(" \t\n"))
                        if c1 not in dict_multi_ratio:
                            dict_multi_ratio[c1] = c2
                        else:
                            dict_multi_ratio[c1] += c2

            call(["rm", file1])
            call(["rm", file2])
            call(["rm", file3])
            call(["rm", file4])
            call(["rm", file5])
            call(["rm", file6])

    outfile1 = "{}/{}.unique.score".format(analysis_folder, library)
    outfile2 = "{}/{}.multi.score".format(analysis_folder, library)
    outfile3 = "{}/{}.unique.mismatch".format(analysis_folder, library)
    outfile4 = "{}/{}.multi.mismatch".format(analysis_folder, library)
    outfile5 = "{}/{}.unique.ratio".format(analysis_folder, library)
    outfile6 = "{}/{}.multi.ratio".format(analysis_folder, library)
    with open(outfile1, "w") as fout:
        for k in dict_unique_score:
            fout.write(k + "\t" + str(dict_unique_score[k]) + "\n")

    with open(outfile2, "w") as fout:
        for k in dict_multi_score:
            fout.write(k + "\t" + str(dict_multi_score[k]) + "\n")

    with open(outfile3, "w") as fout:
        for k in dict_unique_mismatch:
            fout.write(k + "\t" + str(dict_unique_mismatch[k]) + "\n")

    with open(outfile4, "w") as fout:
        for k in dict_multi_mismatch:
            fout.write(k + "\t" + str(dict_multi_mismatch[k]) + "\n")

    with open(outfile5, "w") as fout:
        for k in dict_unique_ratio:
            fout.write(k + "\t" + str(dict_unique_ratio[k]) + "\n")

    with open(outfile6, "w") as fout:
        for k in dict_multi_ratio:
            fout.write(k + "\t" + str(dict_multi_ratio[k]) + "\n")

    # plot
    commandStr = "python {}/plot_alignment_histogram.py {} {} {}".format(
        scripts_folder, analysis_folder, library, library
    )
    os.system(commandStr)

    # Summary mapping rate
    totalreads = 0
    uniquereads = 0
    multireads = 0
    toomanyreads = 0
    for i in range(len(lanes)):
        if libraries[i] != library:
            continue
        for lane_slice in slice_id[lanes[i]]:
            log_file = "{}/{}/{}/{}/{}/{}.{}.{}.{}.{}.star.Log.final.out".format(
                output_folder,
                lanes[i],
                lane_slice,
                library,
                barcodes[i],
                flowcell_barcode,
                lanes[i],
                lane_slice,
                library,
                barcodes[i],
            )
            if not os.path.isfile(log_file):
                continue
            with open(log_file, "r") as f3:
                for line3 in f3:
                    if get_key(line3) == "Number of input reads":
                        totalreads += int(get_val(line3))
                    if get_key(line3) == "Uniquely mapped reads number":
                        uniquereads += int(get_val(line3))
                    if get_key(line3) == "Number of reads mapped to multiple loci":
                        multireads += int(get_val(line3))
                    if get_key(line3) == "Number of reads mapped to too many loci":
                        toomanyreads += int(get_val(line3))

    mismatch1 = 0
    mismatch2 = 0
    mismatch3 = 0
    if "1" in dict_unique_mismatch:
        mismatch1 += dict_unique_mismatch["1"]
    if "1" in dict_multi_mismatch:
        mismatch1 += dict_multi_mismatch["1"]
    if "2" in dict_unique_mismatch:
        mismatch2 += dict_unique_mismatch["2"]
    if "2" in dict_multi_mismatch:
        mismatch2 += dict_multi_mismatch["2"]
    if "3" in dict_unique_mismatch:
        mismatch3 += dict_unique_mismatch["3"]
    if "3" in dict_multi_mismatch:
        mismatch3 += dict_multi_mismatch["3"]
    output_file = "{}/{}_mapping_rate.txt".format(analysis_folder, library)
    with open(output_file, "w") as fout:
        fout.write("library\t{}\n".format(library))
        fout.write("total_reads\t{}\n".format(totalreads))
        fout.write("unique_aligned_reads\t{}\n".format(uniquereads))
        fout.write(
            "unique_aligned_ratio\t{}\n".format(
                "{0:.3g}".format(100 * uniquereads / totalreads)
            )
        )
        fout.write("multi_aligned_reads\t{}\n".format(multireads))
        fout.write(
            "multi_aligned_ratio\t{}\n".format(
                "{0:.3g}".format(100 * multireads / totalreads)
            )
        )
        fout.write("too_many_aligned_reads\t{}\n".format(toomanyreads))
        fout.write(
            "too_many_aligned_ratio\t{}\n".format(
                "{0:.3g}".format(100 * toomanyreads / totalreads)
            )
        )
        fout.write(
            "mismatch1_rate\t{}\n".format(
                "{0:.3g}".format(100 * mismatch1 / totalreads)
            )
        )
        fout.write(
            "mismatch2_rate\t{}\n".format(
                "{0:.3g}".format(100 * mismatch2 / totalreads)
            )
        )
        fout.write(
            "mismatch3_rate\t{}\n".format(
                "{0:.3g}".format(100 * mismatch3 / totalreads)
            )
        )

    if gen_updistance_plot:
        for i in range(len(lanes)):
            if libraries[i] != library:
                continue

            read1_file = "{}/{}.{}.read1.fastq".format(
                analysis_folder, library, lanes[i]
            )
            read2_file = "{}/{}.{}.read2.fastq".format(
                analysis_folder, library, lanes[i]
            )
            combined_bamfile = "{}/{}.{}.unmapped.bam".format(
                analysis_folder, library, lanes[i]
            )
            combined_baifile = "{}/{}.{}.unmapped.bai".format(
                analysis_folder, library, lanes[i]
            )

            commandStr = (
                f"java -Djava.io.tmpdir={tmpdir}"
                " -Dsamjdk.buffer_size=131072 -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50"
                " -XX:GCHeapFreeLimit=10 -Xmx8192m "
            )
            commandStr += (
                "-jar "
                + picard_folder
                + "/picard.jar MergeSamFiles TMP_DIR="
                + tmpdir
                + " CREATE_INDEX=true CREATE_MD5_FILE=false VALIDATION_STRINGENCY=SILENT "
            )
            commandStr += (
                "OUTPUT="
                + combined_bamfile
                + " SORT_ORDER=coordinate ASSUME_SORTED=true"
            )
            for lane_slice in slice_id[lanes[i]]:
                bamfile = "{}/{}.{}.{}.{}".format(
                    analysis_folder, flowcell_barcode, lanes[i], lane_slice, library
                )
                if barcodes[i]:
                    bamfile += "." + barcodes[i]
                bamfile += ".unmapped.bam"
                if not os.path.isfile(bamfile):
                    log.error(
                        f"{flowcell_barcode} - MergeSamFiles error: {bamfile} does not exist!"
                    )
                    raise Exception(f"{bamfile} does not exist!")
                commandStr += " INPUT=" + bamfile
            os.system(commandStr)

            # Convert bam to fastq
            commandStr = (
                f"java -Djava.io.tmpdir={tmpdir}"
                " -Xmx500m -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 "
            )
            commandStr += (
                "-jar "
                + picard_folder
                + "/picard.jar SamToFastq I="
                + combined_bamfile
                + " F="
                + read1_file
                + " F2="
                + read2_file
                + " VALIDATION_STRINGENCY=SILENT"
            )
            os.system(commandStr)

            if os.path.isfile(combined_bamfile):
                call(["rm", combined_bamfile])
            if os.path.isfile(combined_baifile):
                call(["rm", combined_baifile])
            if os.path.isfile(read2_file):
                call(["rm", read2_file])

            output_file = "{}/logs/run_analysis_UPdistance_{}_{}.log".format(
                output_folder, library, lanes[i]
            )
            submission_script = "{}/run_analysis_UPdistance.sh".format(scripts_folder)
            call_args = [
                "qsub",
                "-o",
                output_file,
                "-l",
                "h_vmem=35G",
                "-notify",
                "-l",
                "h_rt=10:0:0",
                "-j",
                "y",
                "-P",
                "macosko_lab",
                "-l",
                "os=RedHat7",
                submission_script,
                manifest_file,
                library,
                lanes[i],
                scripts_folder,
                output_folder,
                analysis_folder,
            ]
            call(call_args)

            break


if __name__ == "__main__":
    main()
