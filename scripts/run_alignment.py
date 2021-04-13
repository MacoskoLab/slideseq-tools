#!/usr/bin/python

# This script is to run the alignment steps

import csv
import os
import re
import sys
import traceback
from datetime import datetime
from subprocess import call

from slideseq.util import str2bool


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


# Write to log file
def write_log(log_file, flowcell_barcode, log_string):
    ...


def main():
    if len(sys.argv) != 6:
        print(
            "Please provide five arguments: manifest file, library ID, lane ID, slice ID and barcode!"
        )
        sys.exit()

    manifest_file = sys.argv[1]
    library = sys.argv[2]
    lane = sys.argv[3]
    lane_slice = sys.argv[4]
    barcode = sys.argv[5]

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
    STAR_folder = (
        options["STAR_folder"]
        if "STAR_folder" in options
        else "/broad/macosko/bin/dropseq-tools/3rdParty/STAR-2.5.2a"
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

    log_file = f"{output_folder}/logs/workflow.log"

    # Read info from metadata file
    bead_structure = ""
    reference = ""

    sequence = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
    base_quality = "10"
    email_address = ""
    experiment_date = ""
    with open(f"{output_folder}/parsed_metadata.txt", "r") as fin:
        reader = csv.reader(fin, delimiter="\t")
        rows = list(reader)
        row0 = rows[0]
        for i in range(1, len(rows)):
            row = rows[i]
            if row[row0.index("library")] == library:
                bead_structure = row[row0.index("bead_structure")]
                reference = row[row0.index("reference")]

                sequence = row[row0.index("start_sequence")]
                base_quality = row[row0.index("base_quality")]
                email_address = row[row0.index("email")]
                experiment_date = row[row0.index("date")]
                break

    reference_folder = reference[: reference.rfind("/")]
    referencePure = reference[reference.rfind("/") + 1 :]
    if referencePure.endswith(".gz"):
        referencePure = referencePure[: referencePure.rfind(".")]
    referencePure = referencePure[: referencePure.rfind(".")]
    genome_dir = f"{reference_folder}/STAR"
    intervals = f"{reference_folder}/{referencePure}.genes.intervals"
    annotations_file = f"{reference_folder}/{referencePure}.gtf"

    prefix_libraries = f"{library_folder}/{experiment_date}_{library}/{flowcell_barcode}.{lane}.{lane_slice}.{library}"
    if barcode:
        prefix_libraries += "." + barcode

    unmapped_bam = prefix_libraries + ".unmapped.bam"
    if not os.path.isfile(unmapped_bam):
        unmapped_bam1 = f"{output_folder}/{lane}/{lane_slice}/{library}/"
        if barcode:
            unmapped_bam1 += f"{barcode}/{flowcell_barcode}.{lane}.{lane_slice}.{library}.{barcode}.unmapped.bam"
        else:
            unmapped_bam1 += f"{flowcell_barcode}.{lane}.{lane_slice}.{library}.unmapped.bam"
        if os.path.isfile(unmapped_bam1):
            os.system("mv " + unmapped_bam1 + " " + unmapped_bam)

    bs_range1 = get_bead_structure_range(bead_structure, "C")
    bs_range2 = get_bead_structure_range(bead_structure, "M")

    folder_running = f"{output_folder}/status/running.alignment_{library}_{lane}_{lane_slice}_{barcode}"
    folder_finished = f"{output_folder}/status/finished.alignment_{library}_{lane}_{lane_slice}_{barcode}"
    folder_failed = f"{output_folder}/status/failed.alignment_{library}_{lane}_{lane_slice}_{barcode}"

    try:
        call(["mkdir", "-p", folder_running])

        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)

        # Tag bam with read sequence extended cellular
        commandStr = (
            dropseq_folder
            + "/TagBamWithReadSequenceExtended O="
            + prefix_libraries
            + ".unaligned_tagged_Cellular.bam COMPRESSION_LEVEL=0 TMP_DIR="
            + tmpdir
        )
        commandStr += (
            " SUMMARY="
            + prefix_libraries
            + ".unaligned_tagged_Cellular.bam_summary.txt BASE_RANGE="
            + bs_range1
            + " BASE_QUALITY="
            + base_quality
        )
        commandStr += (
            " BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 I="
            + unmapped_bam
            + " VALIDATION_STRINGENCY=SILENT"
        )
        write_log(
            log_file,
            flowcell_barcode,
            "TagBamWithReadSequenceExtended Cellular for "
            + library
            + " in Lane "
            + lane
            + " Command="
            + commandStr,
        )
        os.system(commandStr)
        write_log(
            log_file,
            flowcell_barcode,
            "TagBamWithReadSequenceExtended Cellular for "
            + library
            + " in Lane "
            + lane
            + " is done. ",
        )

        unaligned_cellular_file = "{}.unaligned_tagged_Cellular.bam".format(
            prefix_libraries
        )
        if not os.path.isfile(unaligned_cellular_file):
            write_log(
                log_file,
                flowcell_barcode,
                "TagBamWithReadSequenceExtended error: "
                + unaligned_cellular_file
                + " does not exist!",
            )
            raise Exception(
                "TagBamWithReadSequenceExtended error: "
                + unaligned_cellular_file
                + " does not exist!"
            )

        # Tag bam with read sequence extended molecular
        commandStr = (
            dropseq_folder
            + "/TagBamWithReadSequenceExtended O="
            + prefix_libraries
            + ".unaligned_tagged_Molecular.bam COMPRESSION_LEVEL=0 TMP_DIR="
            + tmpdir
        )
        commandStr += (
            " SUMMARY="
            + prefix_libraries
            + ".unaligned_tagged_Molecular.bam_summary.txt BASE_RANGE="
            + bs_range2
            + " BASE_QUALITY="
            + base_quality
        )
        commandStr += (
            " BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 I="
            + prefix_libraries
            + ".unaligned_tagged_Cellular.bam VALIDATION_STRINGENCY=SILENT"
        )
        write_log(
            log_file,
            flowcell_barcode,
            "TagBamWithReadSequenceExtended Molecular for "
            + library
            + " in Lane "
            + lane
            + " Command="
            + commandStr,
        )
        os.system(commandStr)
        write_log(
            log_file,
            flowcell_barcode,
            "TagBamWithReadSequenceExtended Molecular for "
            + library
            + " in Lane "
            + lane
            + " is done. ",
        )

        unaligned_molecular_file = "{}.unaligned_tagged_Molecular.bam".format(
            prefix_libraries
        )
        if not os.path.isfile(unaligned_molecular_file):
            write_log(
                log_file,
                flowcell_barcode,
                "TagBamWithReadSequenceExtended error: "
                + unaligned_molecular_file
                + " does not exist!",
            )
            raise Exception(
                "TagBamWithReadSequenceExtended error: "
                + unaligned_molecular_file
                + " does not exist!"
            )

        if os.path.isfile(unaligned_cellular_file):
            call(["rm", unaligned_cellular_file])

        # Filter low-quality reads
        commandStr = (
            dropseq_folder
            + "/FilterBam TAG_REJECT=XQ I="
            + prefix_libraries
            + ".unaligned_tagged_Molecular.bam "
        )
        commandStr += (
            "O="
            + prefix_libraries
            + ".unaligned.filtered.bam COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=SILENT TMP_DIR="
            + tmpdir
            + " OPTIONS_FILE=/broad/macosko/jilong/slideseq_pipeline/options.txt"
        )
        write_log(
            log_file,
            flowcell_barcode,
            "FilterBam for " + library + " in Lane " + lane + " Command=" + commandStr,
        )
        os.system(commandStr)
        write_log(
            log_file,
            flowcell_barcode,
            "FilterBam for " + library + " in Lane " + lane + " is done. ",
        )

        unaligned_filtered_file = "{}.unaligned.filtered.bam".format(prefix_libraries)
        if not os.path.isfile(unaligned_filtered_file):
            write_log(
                log_file,
                flowcell_barcode,
                "FilterBam error: " + unaligned_filtered_file + " does not exist!",
            )
            raise Exception(
                "FilterBam error: " + unaligned_filtered_file + " does not exist!"
            )

        if os.path.isfile(unaligned_molecular_file):
            call(["rm", unaligned_molecular_file])

        # Trim reads with starting sequence
        commandStr = (
            dropseq_folder
            + "/TrimStartingSequence INPUT="
            + prefix_libraries
            + ".unaligned.filtered.bam OUTPUT="
            + prefix_libraries
            + ".unaligned_trimstartingsequence.filtered.bam COMPRESSION_LEVEL=0 TMP_DIR="
            + tmpdir
            + " "
        )
        commandStr += (
            "OUTPUT_SUMMARY="
            + prefix_libraries
            + ".adapter_trimming_report.txt SEQUENCE="
            + sequence
            + " MISMATCHES=0 NUM_BASES=5 VALIDATION_STRINGENCY=SILENT"
        )
        write_log(
            log_file,
            flowcell_barcode,
            "TrimStartingSequence for "
            + library
            + " in Lane "
            + lane
            + " Command="
            + commandStr,
        )
        os.system(commandStr)
        write_log(
            log_file,
            flowcell_barcode,
            "TrimStartingSequence for " + library + " in Lane " + lane + " is done. ",
        )

        adapter_trim_file = "{}.unaligned_trimstartingsequence.filtered.bam".format(
            prefix_libraries
        )
        if not os.path.isfile(adapter_trim_file):
            write_log(
                log_file,
                flowcell_barcode,
                "TrimStartingSequence error: " + adapter_trim_file + " does not exist!",
            )
            raise Exception(
                "TrimStartingSequence error: " + adapter_trim_file + " does not exist!"
            )

        if os.path.isfile(unaligned_filtered_file):
            call(["rm", unaligned_filtered_file])

        # Adapter-aware poly A trimming
        commandStr = (
            dropseq_folder
            + "/PolyATrimmer I="
            + prefix_libraries
            + ".unaligned_trimstartingsequence.filtered.bam O="
            + prefix_libraries
            + ".unaligned_mc_tagged_polyA_filtered.bam TMP_DIR="
            + tmpdir
            + " "
        )
        commandStr += (
            "OUTPUT_SUMMARY="
            + prefix_libraries
            + ".polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6 VALIDATION_STRINGENCY=SILENT USE_NEW_TRIMMER=true"
        )
        write_log(
            log_file,
            flowcell_barcode,
            "PolyATrimmer for "
            + library
            + " in Lane "
            + lane
            + " Command="
            + commandStr,
        )
        os.system(commandStr)
        write_log(
            log_file,
            flowcell_barcode,
            "PolyATrimmer for " + library + " in Lane " + lane + " is done. ",
        )

        polyA_trim_file = "{}.unaligned_mc_tagged_polyA_filtered.bam".format(
            prefix_libraries
        )
        if not os.path.isfile(polyA_trim_file):
            write_log(
                log_file,
                flowcell_barcode,
                "PolyATrimmer error: " + polyA_trim_file + " does not exist!",
            )
            raise Exception(
                "PolyATrimmer error: " + polyA_trim_file + " does not exist!"
            )

        if os.path.isfile(adapter_trim_file):
            call(["rm", adapter_trim_file])

        # Convert bam to fastq
        commandStr = (
            "java -Djava.io.tmpdir="
            + tmpdir
            + " -Xmx500m -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 "
        )
        commandStr += (
            "-jar "
            + picard_folder
            + "/picard.jar SamToFastq I="
            + prefix_libraries
            + ".unaligned_mc_tagged_polyA_filtered.bam F="
            + prefix_libraries
            + ".fastq VALIDATION_STRINGENCY=SILENT"
        )
        write_log(
            log_file,
            flowcell_barcode,
            "SamToFastq for " + library + " in Lane " + lane + " Command=" + commandStr,
        )
        os.system(commandStr)
        write_log(
            log_file,
            flowcell_barcode,
            "SamToFastq for " + library + " in Lane " + lane + " is done. ",
        )

        fastq_file = "{}.fastq".format(prefix_libraries)
        if not os.path.isfile(fastq_file):
            write_log(
                log_file,
                flowcell_barcode,
                "SamToFastq error: " + fastq_file + " does not exist!",
            )
            raise Exception("SamToFastq error: " + fastq_file + " does not exist!")

        # Map reads to genome sequence using STAR
        commandStr = (
            STAR_folder
            + "/STAR --genomeDir "
            + genome_dir
            + " --readFilesIn "
            + prefix_libraries
            + ".fastq "
        )
        commandStr += (
            "--outFileNamePrefix "
            + prefix_libraries
            + ".star. --outStd Log --outSAMtype BAM Unsorted --outBAMcompression 0"
        )
        if is_NovaSeq or is_NovaSeq_S4:
            commandStr += " --limitOutSJcollapsed 5000000"
        write_log(
            log_file,
            flowcell_barcode,
            "Mapping using STAR for "
            + library
            + " in Lane "
            + lane
            + " Command="
            + commandStr,
        )
        os.system(commandStr)
        write_log(
            log_file,
            flowcell_barcode,
            "Mapping using STAR for " + library + " in Lane " + lane + " is done. ",
        )

        star_file = "{}.star.Aligned.out.bam".format(prefix_libraries)
        if not os.path.isfile(star_file):
            write_log(
                log_file,
                flowcell_barcode,
                "STAR error: " + star_file + " does not exist!",
            )
            raise Exception("STAR error: " + star_file + " does not exist!")

        if os.path.isfile(fastq_file):
            call(["rm", fastq_file])

        # Check alignments quality
        star_file2 = "{}.star.Aligned.out.sam".format(prefix_libraries)
        commandStr = "samtools view -h -o " + star_file2 + " " + star_file
        os.system(commandStr)
        commandStr = "{}/check_alignments_quality {}".format(scripts_folder, star_file2)
        write_log(
            log_file,
            flowcell_barcode,
            "Check alignments quality for "
            + library
            + " in Lane "
            + lane
            + " Command="
            + commandStr,
        )
        os.system(commandStr)
        write_log(
            log_file,
            flowcell_barcode,
            "Check alignments quality for "
            + library
            + " in Lane "
            + lane
            + " is done. ",
        )
        call(["rm", star_file2])

        # Sort aligned bam
        commandStr = (
            "java -Djava.io.tmpdir="
            + tmpdir
            + " -Xmx4000m -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 "
        )
        commandStr += (
            "-jar "
            + picard_folder
            + "/picard.jar SortSam I="
            + prefix_libraries
            + ".star.Aligned.out.bam "
        )
        commandStr += (
            "O="
            + prefix_libraries
            + ".aligned.sorted.bam SORT_ORDER=queryname VALIDATION_STRINGENCY=SILENT TMP_DIR="
            + tmpdir
        )
        write_log(
            log_file,
            flowcell_barcode,
            "SortSam for " + library + " in Lane " + lane + " Command=" + commandStr,
        )
        os.system(commandStr)
        write_log(
            log_file,
            flowcell_barcode,
            "SortSam for " + library + " in Lane " + lane + " is done. ",
        )

        sortsam_file = "{}.aligned.sorted.bam".format(prefix_libraries)
        if not os.path.isfile(sortsam_file):
            write_log(
                log_file,
                flowcell_barcode,
                "SortSam error: " + sortsam_file + " does not exist!",
            )
            raise Exception("SortSam error: " + sortsam_file + " does not exist!")

        if os.path.isfile(star_file):
            call(["rm", star_file])

        # Merge unmapped bam and aligned bam
        commandStr = (
            "java -Djava.io.tmpdir="
            + tmpdir
            + " -Xmx8192m -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 "
        )
        commandStr += (
            "-jar "
            + picard_folder
            + "/picard.jar MergeBamAlignment R="
            + reference
            + " UNMAPPED="
            + prefix_libraries
            + ".unaligned_mc_tagged_polyA_filtered.bam "
        )
        commandStr += (
            "ALIGNED="
            + prefix_libraries
            + ".aligned.sorted.bam O="
            + prefix_libraries
            + ".merged.bam COMPRESSION_LEVEL=0 INCLUDE_SECONDARY_ALIGNMENTS=false CLIP_ADAPTERS=false "
        )
        commandStr += "VALIDATION_STRINGENCY=SILENT TMP_DIR=" + tmpdir
        write_log(
            log_file,
            flowcell_barcode,
            "MergeBamAlignment for "
            + library
            + " in Lane "
            + lane
            + " Command="
            + commandStr,
        )
        os.system(commandStr)
        write_log(
            log_file,
            flowcell_barcode,
            "MergeBamAlignment for " + library + " in Lane " + lane + " is done. ",
        )

        mergedbam_file = "{}.merged.bam".format(prefix_libraries)
        if not os.path.isfile(mergedbam_file):
            write_log(
                log_file,
                flowcell_barcode,
                "MergeBamAlignment error: " + mergedbam_file + " does not exist!",
            )
            raise Exception(
                "MergeBamAlignment error: " + mergedbam_file + " does not exist!"
            )

        if os.path.isfile(polyA_trim_file):
            call(["rm", polyA_trim_file])
        if os.path.isfile(sortsam_file):
            call(["rm", sortsam_file])

        # Tag read with interval
        commandStr = (
            dropseq_folder
            + "/TagReadWithInterval I="
            + prefix_libraries
            + ".merged.bam O="
            + prefix_libraries
            + ".merged.TagReadWithInterval.bam "
        )
        commandStr += (
            "COMPRESSION_LEVEL=0 TMP_DIR="
            + tmpdir
            + " INTERVALS="
            + intervals
            + " TAG=XG VALIDATION_STRINGENCY=SILENT"
        )
        write_log(
            log_file,
            flowcell_barcode,
            "TagReadWithInterval for "
            + library
            + " in Lane "
            + lane
            + " Command="
            + commandStr,
        )
        os.system(commandStr)
        write_log(
            log_file,
            flowcell_barcode,
            "TagReadWithInterval for " + library + " in Lane " + lane + " is done. ",
        )

        merged_taginterval_file = "{}.merged.TagReadWithInterval.bam".format(
            prefix_libraries
        )
        if not os.path.isfile(merged_taginterval_file):
            write_log(
                log_file,
                flowcell_barcode,
                "TagReadWithInterval error: "
                + merged_taginterval_file
                + " does not exist!",
            )
            raise Exception(
                "TagReadWithInterval error: "
                + merged_taginterval_file
                + " does not exist!"
            )

        if os.path.isfile(mergedbam_file):
            call(["rm", mergedbam_file])

        # Tag read with gene function
        commandStr = (
            dropseq_folder
            + "/TagReadWithGeneFunction I="
            + prefix_libraries
            + ".merged.TagReadWithInterval.bam O="
            + prefix_libraries
            + ".star_gene_exon_tagged2.bam "
        )
        commandStr += (
            "ANNOTATIONS_FILE="
            + annotations_file
            + " TMP_DIR="
            + tmpdir
            + " VALIDATION_STRINGENCY=SILENT CREATE_INDEX=false"
        )
        write_log(
            log_file,
            flowcell_barcode,
            "TagReadWithGeneFunction for "
            + library
            + " in Lane "
            + lane
            + " Command="
            + commandStr,
        )
        os.system(commandStr)
        write_log(
            log_file,
            flowcell_barcode,
            "TagReadWithGeneFunction for "
            + library
            + " in Lane "
            + lane
            + " is done. ",
        )

        merged_taggenefunc_file = "{}.star_gene_exon_tagged2.bam".format(
            prefix_libraries
        )
        if not os.path.isfile(merged_taggenefunc_file):
            write_log(
                log_file,
                flowcell_barcode,
                "TagReadWithGeneFunction error: "
                + merged_taggenefunc_file
                + " does not exist!",
            )
            raise Exception(
                "TagReadWithGeneFunction error: "
                + merged_taggenefunc_file
                + " does not exist!"
            )

        if os.path.isfile(merged_taginterval_file):
            call(["rm", merged_taginterval_file])

        ToFolder = "{}/{}/{}/{}/".format(output_folder, lane, lane_slice, library)
        if barcode:
            ToFolder += barcode + "/"
        if os.path.isfile(prefix_libraries + ".star.Log.final.out"):
            call(["mv", prefix_libraries + ".star.Log.final.out", ToFolder])
        if os.path.isfile(prefix_libraries + ".star.Log.out"):
            call(["mv", prefix_libraries + ".star.Log.out", ToFolder])
        if os.path.isfile(prefix_libraries + ".star.Log.progress.out"):
            call(["mv", prefix_libraries + ".star.Log.progress.out", ToFolder])
        if os.path.isfile(prefix_libraries + ".star.SJ.out.tab"):
            call(["mv", prefix_libraries + ".star.SJ.out.tab", ToFolder])
        if os.path.isdir(prefix_libraries + ".star._STARtmp"):
            call(["mv", prefix_libraries + ".star._STARtmp", ToFolder])

        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
        print(dt_string)

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
                + " in lane "
                + lane
                + " slice "
                + lane_slice
                + " failed at the step of running alignment. Please check the log file for the issues. "
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
