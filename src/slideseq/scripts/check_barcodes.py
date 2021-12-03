import gzip
import logging

import click

from slideseq.config import get_config
from slideseq.metadata import Manifest
from slideseq.util.logger import create_logger

log = logging.getLogger(__name__)


@click.command(name="check_barcodes", no_args_is_help=True)
@click.argument("run", metavar="RUN")
@click.option("--debug", is_flag=True, help="Turn on debug logging")
def main(run: str, debug: bool = False):
    """
    For a completed RUN, check for a barcode swap by counting the overlap of exact sequences
    """
    create_logger(debug=debug)
    config = get_config()

    run_dir = config.workflow_dir / run
    manifest_file = run_dir / "manifest.yaml"

    if not run_dir.exists():
        log.error(f"Could not find {run_dir}, has this run completed?")
        return
    elif not manifest_file.exists():
        log.error(f"Could not find {manifest_file}, has this run completed?")

    manifest = Manifest.from_file(manifest_file)

    seq_barcodes = dict()
    bead_barcodes = dict()

    log.info("Reading sequenced barcode and bead barcodes")
    for library in manifest.libraries:
        if library.merged.selected_cells.exists():
            with gzip.open(library.merged.selected_cells, "rt") as fh:
                seq_barcodes[library.name] = {line.strip() for line in fh}
        else:
            log.debug(f"No file {library.merged.selected_cells}, skipping")

        if library.bead_barcodes.exists():
            with library.bead_barcodes.open("r") as fh:
                bead_barcodes[library.name] = {
                    line.strip().replace(",", "") for line in fh
                }
        else:
            log.debug(f"No file {library.bead_barcodes}, skipping")

    log.debug(f"Read data for {len(seq_barcodes)} libraries")

    max_len = max(map(len, bead_barcodes)) + 1

    print(
        f"{'library':{max_len}s}",
        *(f"{n:{max_len}s}" for n in sorted(bead_barcodes)),
        sep="\t",
    )
    for lib_a in sorted(seq_barcodes):
        print(f"{lib_a:{max_len}s}", end="\t")
        for lib_b in sorted(bead_barcodes):
            print(
                f"{len(seq_barcodes[lib_a] & bead_barcodes[lib_b]):{max_len}d}",
                end="\t",
            )
        print()


if __name__ == "__main__":
    main()
