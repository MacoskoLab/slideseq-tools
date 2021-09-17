# number of times to try qsub before giving up
MAX_QSUB = 3

# string used to mean "demux all lanes for this sample"
ALL_LANES = "{LANE}"

# string used to mean "no starting sequence to trim" (used for 10x runs)
NO_START_SEQUENCE = "no start sequence"

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


# columns we need from the sequencing spreadsheet
METADATA_COLS = [
    "library",
    "date",
    "flowcell",
    "run_name",
    "BCLPath",
    "lane",
    "sample_barcode",
    "bead_structure",
    "reference",
    "run_barcodematching",
    "locus_function_list",
    "start_sequence",
    "base_quality",
    "min_transcripts_per_cell",
    "email",
    "puckcaller_path",
    "bead_type",
    "gen_read1_plot",
    "gen_downsampling",
]


# a single library can span multiple rows in the sequencing spreadsheet,
# e.g. multiple flowcells, lanes, and sample barcodes. But the rest of the columns
# should be constant for consistent processing
VARIABLE_LIBRARY_COLS = ["bclpath", "flowcell", "lane", "sample_barcode"]


# for columns that pandas might not recognize automatically
METADATA_TYPES = {
    "date": str,
    "run_barcodematching": bool,
    "base_quality": int,
    "min_transcripts_per_cell": int,
    "gen_read1_plot": bool,
    "gen_downsampling": bool,
}


# columns used from the output of GatherReadQualityMetrics
READ_QUALITY_METRICS = {
    "totalReads": "Total",
    "mappedReads": "Mapped",
    "hqMappedReads": "HQ",
    "hqMappedReadsNoPCRDupes": "HQ No Dupes",
}


# columns used from the output of CollectRnaSeqMetrics
FRAC_INTRONIC_EXONIC = {
    "PCT_RIBOSOMAL_BASES": "ribosomal",
    "PCT_CODING_BASES": "coding",
    "PCT_UTR_BASES": "utr",
    "PCT_INTRONIC_BASES": "intronic",
    "PCT_INTERGENIC_BASES": "intergenic",
}
