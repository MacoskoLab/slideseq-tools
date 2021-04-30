import logging
import pathlib
import sys
import xml.etree.ElementTree as et
from subprocess import PIPE, Popen, run
from typing import Any, Union

import slideseq.util.constants as constants

log = logging.getLogger(__name__)


def get_env_name() -> str:
    # need to pass this to the qsub scripts so they can activate the right environment
    return sys.executable.split("/")[-3]


def qsub_args(
    log_file: pathlib.Path = None,
    debug: bool = False,
    **kwargs: Any,
) -> list[str]:
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


def dropseq_cmd(
    command: str,
    input_file: Union[pathlib.Path, str],
    output_file: Union[pathlib.Path, str],
    mem: str = "8g",
    compression: int = 0,
):
    """Return the beginning of a DropSeq command, with standard options

    :param command: name of the dropseq tool being invoked
    :param input_file: path to the input file
    :param output_file: path to the output file
    :param mem: memory for the heap. default is to share with other jobs
    :param compression: compression level for output. Use 0 for speed, 5 for storage
    """

    return [
        constants.DROPSEQ_DIR / command,
        "-m",
        mem,
        f"I={input_file}",
        f"O={output_file}",
        "VALIDATION_STRINGENCY=SILENT",
        f"COMPRESSION_LEVEL={compression}",
        "VERBOSITY=WARNING",
        "QUIET=true",
    ]


def picard_cmd(command: str, tmp_dir: pathlib.Path, mem: str = "62g"):
    """Return the beginning of a Picard command, with standard options

    :param command: name of the picard tool being invoked
    :param tmp_dir: Location of the tmp directory to use
    :param mem: Memory for the heap. Lower this for piped commands
    """
    return [
        "java",
        f"-Djava.io.tmp_dir={tmp_dir}",
        f"-Xms{mem}",
        f"-Xmx{mem}",
        "-XX:+UseParallelGC",
        "-XX:GCTimeLimit=20",
        "-XX:GCHeapFreeLimit=10",
        "-jar",
        constants.PICARD,
        command,
        "--TMP_DIR",
        tmp_dir,
        "--VALIDATION_STRINGENCY",
        "SILENT",
        "--VERBOSITY",
        "WARNING",
        "--QUIET",
        "true",
    ]


def get_read_structure(run_info_file: pathlib.Path) -> str:
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

    assert len(read_elems) == 3

    # I guess Picard wants this string format for some reason
    return "{}T{}B{}T".format(*(el.get("NumCycles") for el in read_elems))


def get_lanes(run_info_file: pathlib.Path) -> range:
    # open the RunInfo.xml file and parse it with element tree
    with run_info_file.open() as f:
        run_info = et.parse(f)

    lane_count = int(run_info.find("./Run/FlowcellLayout[@LaneCount]").get("LaneCount"))
    return range(1, lane_count + 1)


def run_command(cmd: list[Any], name: str, library: str, lane: int = None):
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
    library: str,
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
