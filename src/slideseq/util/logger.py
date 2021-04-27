import logging
import pathlib
import sys

log = logging.getLogger(__name__)


def create_logger(
    debug: bool = False, dryrun: bool = False, log_file: pathlib.Path = None
):
    root_log = logging.getLogger()

    # google is noisy, turn up its logging level
    logging.getLogger("googleapiclient").setLevel(logging.WARNING)

    if debug:
        root_log.setLevel(logging.DEBUG)
    else:
        root_log.setLevel(logging.INFO)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)
    if dryrun:
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - DRYRUN - %(levelname)s - %(message)s"
        )
    else:
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
    stream_handler.setFormatter(formatter)

    root_log.addHandler(stream_handler)
    log.debug(msg="Added stream handler")

    if log_file is not None:
        log_handler = logging.FileHandler(log_file)
        log_handler.setLevel(logging.DEBUG)
        log_handler.setFormatter(formatter)
        root_log.addHandler(log_handler)
        log.debug(msg="Added file handler")

    if dryrun:
        log.info("DRY RUN ONLY -- No files will be written and no jobs run")
