import logging
from pathlib import Path

import pandas as pd

from slideseq.util import dropseq_cmd, picard_cmd, start_popen

log = logging.getLogger(__name__)


def downsample_dge(
    bam_file: Path, downsample_dir: Path, row: pd.Series, ratio: float, tmp_dir: Path
):
    digital_expression_summary = (
        downsample_dir / f"{row.library}_{ratio}.digital_expression_summary.txt"
    )

    procs = []

    # Downsample reads
    cmd = picard_cmd("DownsampleSam", tmp_dir)
    cmd.extend(["--INPUT", bam_file, "--OUTPUT", "/dev/stdout", "-P", ratio])
    procs.append(start_popen(cmd, "DownsampleSam", row.library))

    # output to /dev/null because we don't want to keep the DGE matrix
    cmd = dropseq_cmd("DigitalExpression", "/dev/stdin", "/dev/null")
    cmd.extend(
        [
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

    procs.append(
        start_popen(cmd, "DigitalExpression", row.library, input_proc=procs[-1])
    )

    # close intermediate streams
    for p in procs[:-1]:
        p.stdout.close()

    # wait for final process to finish
    procs[-1].communicate()
    log.debug(f"Finished with downsampling at ratio {ratio}")
