# This script is to build genome reference

import os
import sys
from subprocess import call

if len(sys.argv) != 2:
    print("Please provide one argument: manifest file!")
    sys.exit()

manifest_file = sys.argv[1]

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

# Check if options exist
if "reference_fasta" not in options:
    print("reference_fasta is not specified in the manifest file. Exiting...")
    sys.exit()

if "gtf" not in options:
    print("gtf is not specified in the manifest file. Exiting...")
    sys.exit()

if "name" not in options:
    print("name is not specified in the manifest file. Exiting...")
    sys.exit()

if "output_folder" not in options:
    print("output_folder is not specified in the manifest file. Exiting...")
    sys.exit()

if "dropseq_folder" not in options:
    print("dropseq_folder is not specified in the manifest file. Exiting...")
    sys.exit()

if "picard_folder" not in options:
    print("picard_folder is not specified in the manifest file. Exiting...")
    sys.exit()

if "STAR_folder" not in options:
    print("STAR_folder is not specified in the manifest file. Exiting...")
    sys.exit()

if "bgzip_location" not in options:
    print("bgzip_location is not specified in the manifest file. Exiting...")
    sys.exit()

if "shell_script" not in options:
    print("shell_script is not specified in the manifest file. Exiting...")
    sys.exit()

reference_fasta = options["reference_fasta"]
gtf = options["gtf"]
name = options["name"]
output_folder = options["output_folder"]
dropseq_folder = options["dropseq_folder"]
picard_folder = options["picard_folder"]
STAR_folder = options["STAR_folder"]
bgzip_location = options["bgzip_location"]
submission_script = options["shell_script"]

MT_SEQUENCE = options["MT_SEQUENCE"] if "MT_SEQUENCE" in options else ""
filtered_gene_biotypes = (
    options["filtered_gene_biotypes"]
    if "filtered_gene_biotypes" in options
    else (
        "processed_pseudogene,unprocessed_pseudogene,transcribed_unprocessed_pseudogene,"
        "pseudogene,IG_V_pseudogene,transcribed_processed_pseudogene,TR_J_pseudogene,"
        "TR_V_pseudogene,unitary_pseudogene,polymorphic_pseudogene,IG_D_pseudogene,"
        "translated_processed_pseudogene,translated_unprocessed_pseudogene,IG_C_pseudogene"
    )
)

# Check if input folders and files exist
if not os.path.isdir(dropseq_folder):
    print("Folder {} does not exist. Exiting...".format(dropseq_folder))
    sys.exit()

if not os.path.isdir(picard_folder):
    print("Folder {} does not exist. Exiting...".format(picard_folder))
    sys.exit()

if not os.path.isdir(STAR_folder):
    print("Folder {} does not exist. Exiting...".format(STAR_folder))
    sys.exit()

if not os.path.isfile(bgzip_location):
    print("File {} does not exist. Exiting...".format(bgzip_location))
    sys.exit()

if not os.path.isfile(reference_fasta):
    print("File {} does not exist. Exiting...".format(reference_fasta))
    sys.exit()

if not os.path.isfile(gtf):
    print("File {} does not exist. Exiting...".format(gtf))
    sys.exit()

if not os.path.isfile(submission_script):
    print("File {} does not exist. Exiting...".format(submission_script))
    sys.exit()

print("Creating genome reference...")

filtered_gene_biotypes2 = ""
if len(filtered_gene_biotypes) > 0:
    types = filtered_gene_biotypes.split(",")
    for gene_type in types:
        filtered_gene_biotypes2 += "G={} ".format(gene_type)

# Create directories
if not os.path.isdir(output_folder):
    call(["mkdir", output_folder])
if not os.path.isdir(f"{output_folder}/STAR"):
    call(["mkdir", f"{output_folder}/STAR"])

sequence_dictionary = f"{output_folder}/{name}.dict"
if os.path.isdir(sequence_dictionary):
    call(["rm", "-r", sequence_dictionary])

output_file = f"{output_folder}/run.log"
call_args = [
    "qsub",
    "-o",
    output_file,
    submission_script,
    dropseq_folder,
    bgzip_location,
    output_folder,
    reference_fasta,
    gtf,
    name,
    filtered_gene_biotypes2,
    MT_SEQUENCE,
]
print(call_args)
call(call_args)
