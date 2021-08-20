from __future__ import annotations

import datetime
import logging
import os
import sys
import xml.etree.ElementTree as et
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


def rsync_to_google(path: Path, gs_path: str):
    # -m for multithreading
    # -q to quiet output
    # -C to continue on errors
    # -e to ignore symlinks
    # -r to recurse into directories
    cmd = ["gsutil", "-m", "-q", "rsync", "-C", "-e", "-r", f"{path}", gs_path]

    proc = run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        log.error(f"Error running gsutil rsync:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info("gsutil rsync completed")

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
        f"{gs_path}/**bam",
    ]
    proc = run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        log.error(f"Error running gsutil setmeta:\n\t{proc.stderr}")
        sys.exit(1)
    else:
        log.info("gsutil setmeta completed")


def qsub_args(log_file: Path = None, debug: bool = False, **kwargs: Any) -> list[str]:
    """
    Returns command list starting with "qsub", adding configured options

    :param log_file: path to log file for output
    :param debug: whether to turn on debug logging
    :param kwargs: additional keyword arguments are passed as environment variables.
                   Values will be converted to strings

    :return: configured list of options, including common args
    """

    # standard options, including resources requests, are specified in the shell scripts
    arg_list = ["qsub"]

    if log_file is not None:
        arg_list.extend(["-o", f"{log_file.absolute().resolve()}"])

    if debug:
        arg_list.extend(["-v", "DEBUG=--debug"])

    for name, value in kwargs.items():
        arg_list.extend(["-v", f"{name}={value}"])

    return arg_list


def get_read_structure(run_info_file: Path) -> str:
    """
    Get read structure from RunInfo.xml. Assumes one index sequence only,
    will error otherwise

    :param run_info_file: path to RunInfo.xml in sequencing directory
    :return: Formatting string representing the read structure
    """
    # open the RunInfo.xml file and parse it with element tree
    with run_info_file.open() as f:
        run_info = et.parse(f)

    read_elems = run_info.findall("./Run/Reads/Read[@NumCycles][@Number]")
    read_elems.sort(key=lambda el: int(el.get("Number")))

    if len(read_elems) == 4:
        # two index reads. We will just ignore the second index
        log.warning(
            "This sequencing run has two index reads, we are ignoring the second one"
        )
        return "{}T{}B{}S{}T".format(*(el.get("NumCycles") for el in read_elems))
    elif len(read_elems) != 3:
        raise ValueError(f"Expected three reads, got {len(read_elems)}")

    return "{}T{}B{}T".format(*(el.get("NumCycles") for el in read_elems))


def get_lanes(run_info_file: Path) -> range:
    # open the RunInfo.xml file and parse it with element tree
    with run_info_file.open() as f:
        run_info = et.parse(f)

    lane_count = int(run_info.find("./Run/FlowcellLayout[@LaneCount]").get("LaneCount"))
    return range(1, lane_count + 1)


def get_flowcell(run_info_file: Path) -> str:
    # open the RunInfo.xml file and parse it with element tree
    with run_info_file.open() as f:
        run_info = et.parse(f)

    flowcell = run_info.find("./Run/Flowcell").text
    return flowcell


def run_command(cmd: list[Any], name: str, library: lib.Library, lane: int = None):
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
