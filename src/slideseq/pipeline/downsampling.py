import logging
from pathlib import Path

from slideseq.library import Library
from slideseq.util import dropseq_cmd, picard_cmd, run_command

log = logging.getLogger(__name__)


def downsample_dge(
    bam_file: Path, downsampled_bam: Path, library: Library, ratio: float, tmp_dir: Path
):
    digital_expression_summary = (
        library.downsample_dir / f"{library}_{ratio:.1f}.digital_expression_summary.txt"
    )

    # Downsample reads
    cmd = picard_cmd("DownsampleSam", tmp_dir)
    # true ratio is r / (r + 0.1) because we are using the previous output
    cmd.extend(
        [
            "--INPUT",
            bam_file,
            "--OUTPUT",
            downsampled_bam,
            "-P",
            f"{ratio / (ratio + 0.1):.3g}",
        ]
    )
    run_command(cmd, "DownsampleSam", library)

    # output to /dev/null because we don't want to keep the DGE matrix
    cmd = dropseq_cmd("DigitalExpression", downsampled_bam, "/dev/null", tmp_dir)
    cmd.extend(
        [
            "CELL_BARCODE_TAG=XC",
            f"SUMMARY={digital_expression_summary}",
            f"MIN_NUM_TRANSCRIPTS_PER_CELL={library.min_transcripts_per_cell}",
            f"READ_MQ={library.base_quality}",
            "OUTPUT_HEADER=false",
            f"UEI={library}",
        ]
    )
    if library.locus_function_list == "intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=null", "LOCUS_FUNCTION_LIST=INTRONIC"])
    elif library.locus_function_list == "exonic+intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=INTRONIC"])

    run_command(cmd, "DigitalExpression", library)

    log.debug(f"Finished with downsampling at ratio {ratio:.1f}")

    return ratio, digital_expression_summary
