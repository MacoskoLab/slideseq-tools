import logging
import pathlib
import sys


def create_logger(debug: bool = False, log_file: pathlib.Path = None):
    log = logging.getLogger()

    # google is noisy, turn up its logging level
    logging.getLogger("googleapiclient").setLevel(logging.WARNING)

    if debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    stream_handler.setFormatter(formatter)

    log.addHandler(stream_handler)
    log.debug(msg="Added stream handler")

    if log_file is not None:
        log_handler = logging.FileHandler(log_file)
        log_handler.setLevel(logging.DEBUG)
        log_handler.setFormatter(formatter)
        log.addHandler(log_handler)
        log.debug(msg="Added file handler")
