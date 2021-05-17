import logging
import os
from pathlib import Path

import pandas as pd

from slideseq.util import dropseq_cmd, picard_cmd, run_command

log = logging.getLogger(__name__)


def downsample_dge(
    bam_file: Path, downsample_dir: Path, row: pd.Series, ratio: float, tmp_dir: Path
):
    downsampled_bam = bam_file.with_suffix(f".downsample_{ratio:.1f}.bam")
    digital_expression_summary = (
        downsample_dir / f"{row.library}_{ratio:.1f}.digital_expression_summary.txt"
    )

    # Downsample reads
    cmd = picard_cmd("DownsampleSam", tmp_dir)
    cmd.extend(["--INPUT", bam_file, "--OUTPUT", downsampled_bam, "-P", f"{ratio:.1f}"])
    run_command(cmd, "DownsampleSam", row.library)

    # output to /dev/null because we don't want to keep the DGE matrix
    cmd = dropseq_cmd("DigitalExpression", downsampled_bam, "/dev/null")
    cmd.extend(
        [
            "CELL_BARCODE_TAG=XC",
            f"SUMMARY={digital_expression_summary}",
            f"MIN_NUM_TRANSCRIPTS_PER_CELL={row.min_transcripts_per_cell}",
            f"READ_MQ={row.base_quality}",
            "OUTPUT_HEADER=false",
            f"UEI={row.library}",
        ]
    )
    if row.locus_function_list == "intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=null", "LOCUS_FUNCTION_LIST=INTRONIC"])
    elif row.locus_function_list == "exonic+intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=INTRONIC"])

    run_command(cmd, "DigitalExpression", row.library)

    log.debug(f"Finished with downsampling at ratio {ratio:.1f}")
    os.remove(downsampled_bam)

    return ratio, digital_expression_summary
