from pathlib import Path

# number of times to try qsub before giving up
MAX_QSUB = 3

DROPSEQ_DIR = Path("/broad/macosko/bin/dropseq-tools")
PICARD = Path("/seq/software/picard-public/2.24.2/picard.jar")

WORKFLOW_DIR = Path("/broad/macosko/data/workflows/flowcell")
LIBRARY_DIR = Path("/broad/macosko/data/libraries")

PLATFORMS = ("MiniSeq", "NextSeq", "NovaSeq", "NovaSeqS4")

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
    "gen_updistance_plot",
    "gen_read1_plot",
    "gen_downsampling",
]
