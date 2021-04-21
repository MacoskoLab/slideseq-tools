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
    stream_handler.setLevel(log_level)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    stream_handler.setFormatter(formatter)

    log.addHandler(stream_handler)
    log.info(msg="Added stream handler")

    log_handler = logging.FileHandler(log_file)
    log_handler.setLevel(log_level)
    log_handler.setFormatter(formatter)
    log.addHandler(log_handler)
