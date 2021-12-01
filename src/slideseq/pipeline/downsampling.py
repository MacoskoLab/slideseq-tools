import logging
import os
from pathlib import Path

import click

from slideseq.config import Config, get_config
from slideseq.library import Library
from slideseq.metadata import Manifest
from slideseq.plot.plot_downsampling import plot_downsampling
from slideseq.util import give_group_access, rsync_to_google, run_command
from slideseq.util.logger import create_logger

log = logging.getLogger(__name__)


def downsample_dge(
    config: Config,
    bam_file: Path,
    downsampled_bam: Path,
    cell_tag: str,
    library: Library,
    ratio: float,
    tmp_dir: Path,
):
    digital_expression_summary = (
        library.downsample_dir / f"{library}_{ratio:.1f}.digital_expression_summary.txt"
    )

    # Downsample reads
    cmd = config.picard_cmd("DownsampleSam", tmp_dir)
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
    cmd = config.dropseq_cmd("DigitalExpression", downsampled_bam, "/dev/null", tmp_dir)
    cmd.extend(
        [
            f"CELL_BARCODE_TAG={cell_tag}",
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


@click.command("downsample_library")
@click.option(
    "-i",
    "--library-index",
    type=int,
    required=True,
    help="Which sample from the metadata to align",
)
@click.option(
    "--manifest-file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="YAML file containing the manifest",
)
@click.option("--debug", is_flag=True, help="Turn on debug logging")
@click.option("--log-file", type=click.Path(exists=False))
def main(
    library_index: int, manifest_file: str, debug: bool = False, log_file: str = None
):
    create_logger(debug=debug, log_file=log_file)
    config = get_config()

    log.debug(f"Reading manifest from {manifest_file}")
    manifest = Manifest.from_file(Path(manifest_file))

    # task array is 1-indexed
    library = manifest.get_library(library_index - 1)

    if not library.gen_downsampling:
        log.debug("Downsampling not requested, nothing to do")
        return

    log.info(f"Downsampling alignments for library {library.name}")
    library.downsample_dir.mkdir(exist_ok=True, parents=True)

    if library.run_barcodematching:
        # this will be much faster on matched bams
        downsample_input = library.matched
        downsample_tag = "XB"
    else:
        downsample_input = library.merged
        downsample_tag = "XC"

    downsample_output = []

    # Progressively downsample the BAM from largest to smallest
    input_bam = downsample_input.bam
    for n in range(9, 0, -1):
        ratio = n / 10

        downsampled_bam = downsample_input.downsampled_bam(ratio)
        downsample_output.append(
            downsample_dge(
                config=config,
                bam_file=input_bam,
                downsampled_bam=downsampled_bam,
                cell_tag=downsample_tag,
                library=library,
                ratio=ratio,
                tmp_dir=manifest.tmp_dir,
            )
        )

        if input_bam != downsample_input.bam:
            os.remove(input_bam)
        input_bam = downsampled_bam

    # remove final downsampled_bam
    if input_bam != downsample_input.bam:
        os.remove(input_bam)

    plot_downsampling(
        downsample_output,
        downsample_input.digital_expression_summary,
        downsample_input.downsampling_pdf,
    )

    log.debug("Setting group permissions")
    give_group_access(library.dir)
    if config.gs_path is not None:
        log.debug("Copying data to google storage")
        rsync_to_google(library.dir, config.gs_path / library.date_name)

    log.info(f"Downsampling for {library} complete")
