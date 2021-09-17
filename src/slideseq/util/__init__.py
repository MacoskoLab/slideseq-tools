from __future__ import annotations

import datetime
import logging
import os
import sys
from pathlib import Path
from subprocess import PIPE, Popen, run
from typing import Any

import slideseq.library as lib

log = logging.getLogger(__name__)


def get_env_name() -> str:
    # need to pass this to the qsub scripts so they can activate the right environment
    return sys.executable.split("/")[-3]


def give_group_access(path: Path):
    """Walks a directory and grants groups permissions.

    Sets rwx on drs and rw on dirs"""
    for dirpath, dirnames, filenames in os.walk(path):
        # set ug+rwx, o+r on directories
        os.chmod(dirpath, 0o774)
        for filename in filenames:
            # set ug+rw, o+r on files
            os.chmod(os.path.join(dirpath, filename), 0o664)


def rsync_to_google(path: Path, gs_path: Path):
    """
    Call gsutil rsync to upload a directory up to Google Storage, and sets Custom-Time
    metadata on BAM files for lifecycle rules.

    :param path: Local path to upload
    :param gs_path: Path to the destination, including bucket name but *not* 'gs://'
    """
    # -m for multithreading
    # -q to quiet output
    # -C to continue on errors
    # -e to ignore symlinks
    # -r to recurse into directories
    cmd = [
        "gsutil",
        "-m",
        "-q",
        "rsync",
        "-C",
        "-e",
        "-r",
        f"{path}",
        f"gs://{gs_path}",
    ]

    run_command(cmd, "gsutil rsync", path)

    date = datetime.datetime.now(tz=datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )

    cmd = [
        "gsutil",
        "-m",
        "-q",
        "setmeta",
        "-h",
        f"Custom-Time:{date}",
        f"gs://{gs_path}/**bam",
    ]

    run_command(cmd, "gsutil setmeta", path)


def qsub_args(
    log_file: Path = None, email: str = None, debug: bool = False, **kwargs: Any
) -> list[str]:
    """
    Returns command list starting with "qsub", adding configured options

    :param log_file: path to log file for output
    :param email: Email addresses to include with the -M option. If None, emails will
                  be sent to the submitting user
    :param debug: whether to turn on debug logging
    :param kwargs: additional keyword arguments are passed as environment variables.
                   Values will be converted to strings

    :return: configured list of options, including common args
    """

    # standard options, including resources requests, are specified in the shell scripts
    arg_list = ["qsub"]

    if log_file is not None:
        arg_list.extend(["-o", f"{log_file.absolute().resolve()}"])

    if email is not None:
        arg_list.extend(["-M", email])

    if debug:
        arg_list.extend(["-v", "DEBUG=--debug"])

    for name, value in kwargs.items():
        arg_list.extend(["-v", f"{name}={value}"])

    return arg_list


def run_command(cmd: list[Any], name: str, library: Any, lane: int = None):
    if lane is None:
        log.info(f"{name} for {library}")
    else:
        log.info(f"{name} for {library} in lane {lane}")

    # convert args to strings rather than relying on the caller
    cmd = [str(arg) for arg in cmd]
    log.debug(f"Command = {' '.join(cmd)}")

    proc = run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        log.error(f"Error running {name}:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info(f"{name} completed")


def start_popen(
    cmd: list[Any],
    name: str,
    library: lib.Library,
    lane: int = None,
    input_proc: Popen = None,
):
    if lane is None:
        log.info(f"{name} for {library}")
    else:
        log.info(f"{name} for {library} in lane {lane}")

    # convert args to strings rather than relying on the caller
    cmd = [str(arg) for arg in cmd]
    log.debug(f"Command = {' '.join(cmd)}")

    if input_proc is not None:
        return Popen(cmd, stdin=input_proc.stdout, stdout=PIPE)
    else:
        return Popen(cmd, stdout=PIPE)
