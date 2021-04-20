#!/usr/bin/python

# This script is to run the alignment steps

import csv
import logging
import os
import re
import sys
from subprocess import call

from slideseq.util import str2bool

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
    if len(sys.argv) != 6:
        log.error(
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
        log.error("File {} does not exist. Exiting...".format(manifest_file))
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

    # Read info from metadata file
    bead_structure = ""
    reference = ""

    sequence = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
    base_quality = "10"
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

    unmapped_bam = f"{prefix_libraries}.unmapped.bam"
    if not os.path.isfile(unmapped_bam):
        unmapped_bam1 = f"{output_folder}/{lane}/{lane_slice}/{library}/"
        if barcode:
            unmapped_bam1 += f"{barcode}/{flowcell_barcode}.{lane}.{lane_slice}.{library}.{barcode}.unmapped.bam"
        else:
            unmapped_bam1 += (
                f"{flowcell_barcode}.{lane}.{lane_slice}.{library}.unmapped.bam"
            )
        if os.path.isfile(unmapped_bam1):
            os.system(f"mv {unmapped_bam1} {unmapped_bam}")

    bs_range1 = get_bead_structure_range(bead_structure, "C")
    bs_range2 = get_bead_structure_range(bead_structure, "M")

    # Tag bam with read sequence extended cellular
    commandStr = (
        f"{dropseq_folder}/TagBamWithReadSequenceExtended O={prefix_libraries}.unaligned_tagged_Cellular.bam"
        f" COMPRESSION_LEVEL=0 TMP_DIR={tmpdir} SUMMARY={prefix_libraries}.unaligned_tagged_Cellular.bam_summary.txt"
        f" BASE_RANGE={bs_range1} BASE_QUALITY={base_quality} BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC"
        f" NUM_BASES_BELOW_QUALITY=1 I={unmapped_bam} VALIDATION_STRINGENCY=SILENT"
    )
    log.info(
        f"{flowcell_barcode} TagBamWithReadSequenceExtended Cellular for {library} in lane {lane}"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} TagBamWithReadSequenceExtended Cellular for {library} in lane {lane} is done"
    )

    unaligned_cellular_file = "{}.unaligned_tagged_Cellular.bam".format(
        prefix_libraries
    )
    if not os.path.isfile(unaligned_cellular_file):
        log.error(
            f"{flowcell_barcode} - TagBamWithReadSequenceExtended error: {unaligned_cellular_file} does not exist"
        )
        raise Exception(
            f"TagBamWithReadSequenceExtended error: {unaligned_cellular_file} does not exist"
        )

    # Tag bam with read sequence extended molecular
    commandStr = (
        f"{dropseq_folder}/TagBamWithReadSequenceExtended O={prefix_libraries}.unaligned_tagged_Molecular.bam"
        f" COMPRESSION_LEVEL=0 TMP_DIR={tmpdir} SUMMARY={prefix_libraries}.unaligned_tagged_Molecular.bam_summary.txt"
        f" BASE_RANGE={bs_range2} BASE_QUALITY={base_quality} BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM"
        f" NUM_BASES_BELOW_QUALITY=1 I={prefix_libraries}.unaligned_tagged_Cellular.bam VALIDATION_STRINGENCY=SILENT"
    )
    log.info(
        f"{flowcell_barcode} - TagBamWithReadSequenceExtended Molecular for {library} in lane {lane}"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - TagBamWithReadSequenceExtended Molecular for {library} in lane {lane} done"
    )

    unaligned_molecular_file = f"{prefix_libraries}.unaligned_tagged_Molecular.bam"
    if not os.path.isfile(unaligned_molecular_file):
        log.error(
            f"TagBamWithReadSequenceExtended error: {unaligned_molecular_file} does not exist!"
        )
        raise Exception(
            f"TagBamWithReadSequenceExtended error: {unaligned_molecular_file} does not exist!"
        )

    if os.path.isfile(unaligned_cellular_file):
        call(["rm", unaligned_cellular_file])

    # Filter low-quality reads
    commandStr = (
        f"{dropseq_folder}/FilterBam TAG_REJECT=XQ I={prefix_libraries}.unaligned_tagged_Molecular.bam"
        f" O={prefix_libraries}.unaligned.filtered.bam COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=SILENT"
        f" TMP_DIR={tmpdir} OPTIONS_FILE=/broad/macosko/jilong/slideseq_pipeline/options.txt"
    )
    log.info(f"{flowcell_barcode} - FilterBam for {library} in Lane {lane}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - FilterBam for {library} in Lane {lane} is done")

    unaligned_filtered_file = f"{prefix_libraries}.unaligned.filtered.bam"
    if not os.path.isfile(unaligned_filtered_file):
        log.error(f"FilterBam error: {unaligned_filtered_file} does not exist!")
        raise Exception(f"FilterBam error: {unaligned_filtered_file} does not exist!")

    if os.path.isfile(unaligned_molecular_file):
        call(["rm", unaligned_molecular_file])

    # Trim reads with starting sequence
    commandStr = (
        f"{dropseq_folder}/TrimStartingSequence INPUT={prefix_libraries}.unaligned.filtered.bam"
        f" OUTPUT={prefix_libraries}.unaligned_trimstartingsequence.filtered.bam COMPRESSION_LEVEL=0 TMP_DIR={tmpdir}"
        f" OUTPUT_SUMMARY={prefix_libraries}.adapter_trimming_report.txt SEQUENCE={sequence}"
        f" MISMATCHES=0 NUM_BASES=5 VALIDATION_STRINGENCY=SILENT"
    )
    log.info(
        f"{flowcell_barcode} - TrimStartingSequence for {library} in Lane {lane}"
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - TrimStartingSequence for {library} in lane {lane} is done."
    )

    adapter_trim_file = f"{prefix_libraries}.unaligned_trimstartingsequence.filtered.bam"
    if not os.path.isfile(adapter_trim_file):
        log.error(
            f"{flowcell_barcode} - TrimStartingSequence error: {adapter_trim_file} does not exist!"
        )
        raise Exception(
            f"TrimStartingSequence error: {adapter_trim_file} does not exist!"
        )

    if os.path.isfile(unaligned_filtered_file):
        call(["rm", unaligned_filtered_file])

    # Adapter-aware poly A trimming
    commandStr = (
        f"{dropseq_folder}/PolyATrimmer I={prefix_libraries}.unaligned_trimstartingsequence.filtered.bam"
        f" O={prefix_libraries}.unaligned_mc_tagged_polyA_filtered.bam TMP_DIR={tmpdir}"
        f" OUTPUT_SUMMARY={prefix_libraries}.polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6"
        f" VALIDATION_STRINGENCY=SILENT USE_NEW_TRIMMER=true"
    )
    log.info(f"{flowcell_barcode} - PolyATrimmer for {library} in Lane {lane}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - PolyATrimmer for {library} in Lane {lane} is done.")

    polyA_trim_file = f"{prefix_libraries}.unaligned_mc_tagged_polyA_filtered.bam"
    if not os.path.isfile(polyA_trim_file):
        log.info(
            f"{flowcell_barcode} - PolyATrimmer error: {polyA_trim_file} does not exist!"
        )
        raise Exception(
            f"PolyATrimmer error: {polyA_trim_file} does not exist!",
        )

    if os.path.isfile(adapter_trim_file):
        call(["rm", adapter_trim_file])

    # Convert bam to fastq
    commandStr = (
        f"java -Djava.io.tmpdir={tmpdir} -Xmx500m -XX:+UseParallelGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=20"
        f" -XX:GCHeapFreeLimit=10 -jar {picard_folder}/picard.jar SamToFastq"
        f" I={prefix_libraries}.unaligned_mc_tagged_polyA_filtered.bam F={prefix_libraries}.fastq"
        f" VALIDATION_STRINGENCY=SILENT"
    )
    log.info(f"{flowcell_barcode} - SamToFastq for {library} in Lane {lane}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - SamToFastq for {library} in Lane {lane} is done.")

    fastq_file = f"{prefix_libraries}.fastq"
    if not os.path.isfile(fastq_file):
        log.error(
            f"{flowcell_barcode} - SamToFastq error: {fastq_file} does not exist!"
        )
        raise Exception(f"SamToFastq error: {fastq_file} does not exist!")

    # Map reads to genome sequence using STAR
    commandStr = (
        f"{STAR_folder}/STAR --genomeDir {genome_dir} --readFilesIn {prefix_libraries}.fastq"
        f" --outFileNamePrefix {prefix_libraries}.star. --outStd Log --outSAMtype BAM Unsorted --outBAMcompression 0"
        " --limitOutSJcollapsed 5000000" if is_NovaSeq else ""
    )
    log.info(f"{flowcell_barcode} - Mapping using STAR for {library} in Lane {lane}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - Mapping using STAR for {library} in lane {lane} is done."
    )

    star_file = f"{prefix_libraries}.star.Aligned.out.bam"
    if not os.path.isfile(star_file):
        log.error(f"{flowcell_barcode} - STAR error: {star_file} does not exist!")
        raise Exception(f"STAR error: {star_file} does not exist!")

    if os.path.isfile(fastq_file):
        call(["rm", fastq_file])

    # Check alignments quality
    star_file2 = f"{prefix_libraries}.star.Aligned.out.sam"
    commandStr = f"samtools view -h -o {star_file2} {star_file}"
    os.system(commandStr)
    commandStr = f"{scripts_folder}/check_alignments_quality {star_file2}"
    log.info(
        f"{flowcell_barcode} - Check alignments quality for {library} in Lane {lane} is done."
    )
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - Check alignments quality for"
        f" {library} in Lane {lane} is done."
    )
    call(["rm", star_file2])

    # Sort aligned bam
    commandStr = (
        f"java -Djava.io.tmpdir={tmpdir} -Xmx4000m -XX:+UseParallelGC"
        f" -XX:GCTimeLimit=20 -XX:GCHeapFreeLimit=10"
        f" -jar {picard_folder}/picard.jar SortSam"
        f" -I {prefix_libraries}.star.Aligned.out.bam"
        f" -O {prefix_libraries}.aligned.sorted.bam"
        f" --SORT_ORDER queryname --VALIDATION_STRINGENCY SILENT --TMP_DIR {tmpdir}"
    )
    log.info(f"{flowcell_barcode} - SortSam for {library} in lane {lane}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - SortSam for {library} in lane {lane} is done.")

    sortsam_file = f"{prefix_libraries}.aligned.sorted.bam"
    if not os.path.isfile(sortsam_file):
        log.error(f"{flowcell_barcode} - SortSam error: {sortsam_file} does not exist!")
        raise Exception(f"SortSam error: {sortsam_file} does not exist!")

    if os.path.isfile(star_file):
        call(["rm", star_file])

    # Merge unmapped bam and aligned bam
    commandStr = (
        f"java -Djava.io.tmpdir={tmpdir} -Xmx8192m -XX:+UseParallelGC -XX:GCTimeLimit=20"
        f" -XX:GCHeapFreeLimit=10 -jar {picard_folder}/picard.jar MergeBamAlignment -R {reference}"
        f" --UNMAPPED {prefix_libraries}.unaligned_mc_tagged_polyA_filtered.bam"
        f" --ALIGNED {prefix_libraries}.aligned.sorted.bam -O {prefix_libraries}.merged.bam"
        f" --COMPRESSION_LEVEL 0 --INCLUDE_SECONDARY_ALIGNMENTS false --CLIP_ADAPTERS false"
        f" --VALIDATION_STRINGENCY SILENT --TMP_DIR {tmpdir}"
    )
    log.info(f"{flowcell_barcode} - MergeBamAlignment for {library} in lane {lane}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(
        f"{flowcell_barcode} - MergeBamAlignment for {library} in Lane {lane} is done."
    )

    mergedbam_file = f"{prefix_libraries}.merged.bam"
    if not os.path.isfile(mergedbam_file):
        log.info(
            f"{flowcell_barcode} - MergeBamAlignment error: {mergedbam_file} does not exist!",
        )
        raise Exception(f"MergeBamAlignment error: {mergedbam_file} does not exist!")

    if os.path.isfile(polyA_trim_file):
        call(["rm", polyA_trim_file])
    if os.path.isfile(sortsam_file):
        call(["rm", sortsam_file])

    # Tag read with interval
    commandStr = (
        f"{dropseq_folder}/TagReadWithInterval I={prefix_libraries}.merged.bam"
        f" O={prefix_libraries}.merged.TagReadWithInterval.bam"
        f" COMPRESSION_LEVEL=0 TMP_DIR={tmpdir} INTERVALS={intervals} TAG=XG VALIDATION_STRINGENCY=SILENT"
    )
    log.info(f"{flowcell_barcode} - TagReadWithInterval for {library} in Lane {lane}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - TagReadWithInterval for {library} in Lane {lane} is done.")

    merged_taginterval_file = f"{prefix_libraries}.merged.TagReadWithInterval.bam"
    if not os.path.isfile(merged_taginterval_file):
        log.error(f"{flowcell_barcode} - TagReadWithInterval error: {merged_taginterval_file} does not exist!")
        raise Exception(
            f"TagReadWithInterval error: {merged_taginterval_file} does not exist!"
        )

    if os.path.isfile(mergedbam_file):
        call(["rm", mergedbam_file])

    # Tag read with gene function
    commandStr = (
        f"{dropseq_folder}/TagReadWithGeneFunction I={prefix_libraries}.merged.TagReadWithInterval.bam"
        f" O={prefix_libraries}.star_gene_exon_tagged2.bam"
        f" ANNOTATIONS_FILE={annotations_file} TMP_DIR={tmpdir} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=false"
    )
    log.info(f"{flowcell_barcode} - TagReadWithGeneFunction for {library} in lane {lane}")
    log.debug(f"Command = {commandStr}")
    os.system(commandStr)
    log.info(f"{flowcell_barcode} - TagReadWithGeneFunction for {library} in Lane {lane} is done.")

    merged_taggenefunc_file = f"{prefix_libraries}.star_gene_exon_tagged2.bam"
    if not os.path.isfile(merged_taggenefunc_file):
        log.error(
            f"{flowcell_barcode} - TagReadWithGeneFunction error: {merged_taggenefunc_file} does not exist!"
        )
        raise Exception(
            f"TagReadWithGeneFunction error: {merged_taggenefunc_file} does not exist!"
        )

    if os.path.isfile(merged_taginterval_file):
        call(["rm", merged_taginterval_file])

    ToFolder = f"{output_folder}/{lane}/{lane_slice}/{library}/"
    if barcode:
        ToFolder += f"{barcode}/"

    for file_name in (
        f"{prefix_libraries}.star.Log.final.out",
        f"{prefix_libraries}.star.Log.out",
        f"{prefix_libraries}.star.Log.progress.out",
        f"{prefix_libraries}.star.SJ.out.tab",
        f"{prefix_libraries}.star.Log.final.out",
    ):
        if os.path.isfile(file_name):
            call(["mv", file_name, ToFolder])


if __name__ == "__main__":
    main()
