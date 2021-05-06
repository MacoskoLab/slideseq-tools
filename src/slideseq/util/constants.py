from pathlib import Path

# number of times to try qsub before giving up
MAX_QSUB = 3

DROPSEQ_DIR = Path("/broad/macosko/bin/dropseq-tools-2.4.0")
PICARD = Path("/seq/software/picard-public/2.24.2/picard.jar")

REFERENCE_DIR = Path("/broad/macosko/reference")
# these biotypes are _removed_ from the reference GTF
FILTERED_BIOTYPES = [
    "processed_pseudogene",
    "unprocessed_pseudogene",
    "transcribed_unprocessed_pseudogene",
    "pseudogene",
    "IG_V_pseudogene",
    "transcribed_processed_pseudogene",
    "TR_J_pseudogene",
    "TR_V_pseudogene",
    "unitary_pseudogene",
    "polymorphic_pseudogene",
    "IG_D_pseudogene",
    "translated_processed_pseudogene",
    "translated_unprocessed_pseudogene",
    "IG_C_pseudogene",
]

WORKFLOW_DIR = Path("/broad/macosko/data/workflows/flowcell")
LIBRARY_DIR = Path("/broad/macosko/data/libraries")

PLATFORMS = {"MiniSeq", "NextSeq", "NovaSeq", "NovaSeqS4"}

# columns we need from the sequencing spreadsheet
METADATA_COLS = [
    "library",
    "date",
    "flowcell",
    "BCLPath",
    "IlluminaPlatform",
    "resubmit",
    "lane",
    "sample_barcode",
    "bead_structure",
    "estimated_num_cells",
    "estimated_num_beads",
    "reference",
    "run_barcodematching",
    "locus_function_list",
    "start_sequence",
    "base_quality",
    "min_transcripts_per_cell",
    "email",
    "puckcaller_path",
    "bead_type",
    "fdr_threshold",
    "gen_read1_plot",
    "gen_downsampling",
]

# for columns that pandas might not recognize automatically
METADATA_TYPES = {
    "date": str,
    "estimated_num_cells": int,
    "estimated_num_beads": int,
    "run_barcodematching": bool,
    "base_quality": int,
    "min_transcripts_per_cell": int,
    "fdr_threshold": float,
    "gen_read1_plot": bool,
    "gen_downsampling": bool,
}


# columns used from the output of GatherReadQualityMetrics
READ_QUALITY_METRICS = {
    "totalReads": "Total",
    "mappedReads": "Mapped",
    "hqMappedReads": "HQ",
    "hqMappedReadNoPCRDupes": "HQ No Dupes",
}


# columns used from the output of CollectRnaSeqMetrics
FRAC_INTRONIC_EXONIC = {
    "PCT_RIBOSOMAL_BASES": "ribosomal",
    "PCT_CODING_BASES": "coding_bases",
    "PCT_UTR_BASES": "utr_bases",
    "PCT_INTRONIC_BASES": "intronic",
    "PCT_INTERGENIC_BASES": "intergenic",
}
