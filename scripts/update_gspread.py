#!/usr/bin/python

# This script is to submit a request to the Slide-seq flowcell alignment pipeline

import argparse
import os
import re
import time
# silence warnings for pandas below
import warnings

import gspread
import numpy as np
from oauth2client.service_account import ServiceAccountCredentials

scope = [
    "https://spreadsheets.google.com/feeds",
    "https://www.googleapis.com/auth/drive",
]

json_file = "/broad/macosko/jilong/slideseq_pipeline/slideseq-6e435328493d.json"
credentials = ServiceAccountCredentials.from_json_keyfile_name(json_file, scope)

gc = gspread.authorize(credentials)

parser = argparse.ArgumentParser(
    description="Submit flowcells to the Slide-seq flowcell alignment pipeline."
)

parser.add_argument(
    "flowcells",
    metavar="flowcells",
    type=str,
    nargs="+",
    help="Names of flowcells to submit for processing (separated by spaces)",
)
parser.add_argument(
    "--spreadsheet",
    dest="spreadsheet",
    metavar="spreadsheet",
    default="1kwnKrkbl80LyE9lND0UZZJXipL4yfBbGjkTe6hcwJic",
    action="store",
    help="Optional spreadsheet key to open and populate (default is for Macosko Slide-seq Flowcell Alignment)",
)
parser.add_argument(
    "--worksheet_ind",
    dest="wks_ind",
    metavar="worksheet_ind",
    type=int,
    default=0,
    help="Which worksheet to open in spreadsheet (0-indexed) (default 0)",
)

args = parser.parse_args()

if os.path.isfile(args.flowcells[0]):
    with open(args.flowcells[0], "r") as name_file:
        lines = name_file.read().splitlines()
    flowcells = lines
else:
    flowcells = args.flowcells

print("Beginning update google spread for flowcells: {}".format(flowcells))

# convert non-decimals
non_decimal = re.compile(r"[^\d.]+")

wks = gc.open_by_key(args.spreadsheet).get_worksheet(args.wks_ind)

# IMPORTANT NOTE: row and columns are 1-indexed for range and other utilities
# these column indices are hardcoded for now based on spreadsheet organization
lib_i = 1
date_i = 4
reference_i = 14
locus_function_list_i = 16
base_quality_i = 18
min_transcripts_per_cell_i = 19
email_i = 20
puckcaller_i = 21

workflow_dir = "/broad/macosko/data/workflows/flowcell"
library_dir = "/broad/macosko/data/libraries"

for flowcell in flowcells:
    if "-" in flowcell and not args.diff_align:
        warnings.warn(
            (
                "Flowcell name follows syntax for different secondary alignment;"
                + " please set --diff_align flag before running."
            ).format(flowcell)
        )
        continue
    flow_cells = wks.findall(flowcell)
    if len(np.unique([cell.col for cell in flow_cells])) > 1:
        warnings.warn(
            (
                "Flowcell {} found in multiple columns;"
                + " not running until resolved."
            ).format(flowcell)
        )
        continue
    elif len(flow_cells) < 1:
        warnings.warn(
            "Flowcell {} not found in spreadsheet; please add to sheet.".format(
                flowcell
            )
        )
        continue
    flow_rows = [cell.row for cell in flow_cells]

    for row in flow_rows:
        library = wks.cell(row, 1).value
        alignment_folder = (
            library_dir + "/" + wks.cell(row, date_i).value + "_" + library
        )
        wks.update_acell("AA" + str(row), alignment_folder)

        file = "{}/{}.ReadQualityMetrics.txt".format(alignment_folder, library)
        total_reads = 0
        hq_reads = 0
        if os.path.isfile(file):
            mat = np.loadtxt(
                file,
                delimiter="\t",
                dtype="int",
                skiprows=3,
                max_rows=1,
                usecols=(1, 2, 3, 4),
            )
            wks.update_acell("AD" + str(row), str(mat[0]))
            wks.update_acell("AE" + str(row), str(mat[1]))
            wks.update_acell("AF" + str(row), str(mat[2]))
            total_reads = mat[0]
            hq_reads = mat[2]

        file = "{}/BeadBarcodes.txt".format(wks.cell(row, puckcaller_i).value)
        if os.path.isfile(file):
            l = 0
            with open(file, "r") as fin:
                for line in fin:
                    l += 1
            wks.update_acell("AB" + str(row), str(l))

        file = "{}/{}.fracIntronicExonic.txt".format(alignment_folder, library)
        if os.path.isfile(file):
            mat = np.loadtxt(
                file,
                delimiter="\t",
                dtype="float",
                skiprows=7,
                max_rows=1,
                usecols=(15, 16, 17, 18, 19),
            )
            mat[1] += mat[2]
            mat[2] = mat[1] + mat[3]
            # 'ribosomal','exonic','genic','intronic','intergenic'
            wks.update_acell("AG" + str(row), str(int(mat[4] * hq_reads)))
            wks.update_acell("AH" + str(row), str(int(mat[3] * hq_reads)))
            wks.update_acell("AI" + str(row), str(int(mat[1] * hq_reads)))
            wks.update_acell("AJ" + str(row), str(int(mat[2] * hq_reads)))
            wks.update_acell("AK" + str(row), str(int(mat[0] * hq_reads)))

        locus_function_list = wks.cell(row, locus_function_list_i).value
        lists = locus_function_list.split(",")
        reference = wks.cell(row, reference_i).value
        referencePure = reference[reference.rfind("/") + 1 :]
        if referencePure.endswith(".gz"):
            referencePure = referencePure[: referencePure.rfind(".")]
        referencePure = referencePure[: referencePure.rfind(".")]
        base_quality = wks.cell(row, base_quality_i).value
        min_transcripts_per_cell = wks.cell(row, min_transcripts_per_cell_i).value
        num_cells = ""
        fraction_reads = ""
        num_matched_cells = ""
        mean_reads = ""
        median_genes = ""
        median_umi = ""
        for li in lists:
            reference2 = referencePure + "." + li
            selected_cells = (
                "{}/{}/alignment/{}.{}_transcripts_mq_{}_selected_cells.txt".format(
                    alignment_folder,
                    reference2,
                    library,
                    min_transcripts_per_cell,
                    base_quality,
                )
            )
            if os.path.isfile(selected_cells):
                l = 0
                with open(selected_cells, "r") as fin:
                    for line in fin:
                        l += 1
                num_cells += str(l) + ";"

            file = "{}/{}/alignment/{}.ReadQualityMetrics.txt".format(
                alignment_folder, reference2, library
            )
            if os.path.isfile(file):
                mat = np.loadtxt(
                    file,
                    delimiter="\t",
                    dtype="int",
                    skiprows=3,
                    max_rows=1,
                    usecols=1,
                )
                fraction_reads += str(mat * 100 // total_reads) + "%;"

            dge_summary = "{}/{}/alignment/{}.digital_expression_summary.txt".format(
                alignment_folder, reference2, library
            )
            if os.path.isfile(dge_summary):
                dge_summary_reads = np.loadtxt(
                    dge_summary, delimiter="\t", dtype="int", skiprows=7, usecols=1
                )
                dge_summary_umi = np.loadtxt(
                    dge_summary, delimiter="\t", dtype="int", skiprows=7, usecols=2
                )
                dge_summary_gene = np.loadtxt(
                    dge_summary, delimiter="\t", dtype="int", skiprows=7, usecols=3
                )
                cou = len(dge_summary_reads)
                num_matched_cells += str(cou) + ";"
                mean_reads += str(sum(dge_summary_reads) // cou) + ";"
                median_genes += str(dge_summary_gene[cou // 2]) + ";"
                median_umi += str(dge_summary_umi[cou // 2]) + ";"
        wks.update_acell("AC" + str(row), num_cells[:-1])
        wks.update_acell("AL" + str(row), fraction_reads[:-1])
        wks.update_acell("AM" + str(row), num_matched_cells[:-1])
        wks.update_acell("AN" + str(row), mean_reads[:-1])
        wks.update_acell("AO" + str(row), median_genes[:-1])
        wks.update_acell("AP" + str(row), median_umi[:-1])

        time.sleep(30)
