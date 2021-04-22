#!/usr/bin/env python

import pathlib
import sys
import xml.etree.ElementTree as et
from typing import Any

import slideseq.util.constants as constants


def get_env_name() -> str:
    # need to pass this to the qsub scripts so they can activate the right environment
    return sys.executable.split("/")[-3]


def qsub_args(
    log_file: pathlib.Path = None,
    include_env: bool = False,
    **kwargs: Any,
) -> list[str]:
    """
    Returns command list starting with "qsub", adding configured options

    :param log_file: path to log file for output
    :param include_env: pass in the current conda environment as a variable
    :param kwargs: additional keyword arguments are passed as environment variables.
                   Values will be converted to strings

    :return: configured list of options, including common args
    """

    # standard options, including resources requests, are specified in the shell scripts
    arg_list = ["qsub"]

    if log_file is not None:
        arg_list.extend(["-o", f"{log_file.absolute().resolve()}"])

    if include_env:
        arg_list.extend(["-v", f"CONDA_ENV={get_env_name()}"])

    for name, value in kwargs.items():
        arg_list.extend(["-v", f"{name}='{value}'"])

    return arg_list


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


def get_tiles(run_info_file: pathlib.Path, lane: str) -> list[str]:
    # open the RunInfo.xml file and parse it with element tree
    with run_info_file.open() as f:
        run_info = et.parse(f)

    # extract all tile elements
    tile_elems = run_info.findall("./Run/FlowcellLayout/TileSet/Tiles/Tile")
    # convert to tuple of (lane, tile)
    tiles = [el.text.split("_") for el in tile_elems]
    # filter for desired lane, and sort
    tiles = sorted(tile for ln, tile in tiles if ln == lane)

    return tiles


def get_lanes(run_info_file: pathlib.Path) -> range:
    # open the RunInfo.xml file and parse it with element tree
    with run_info_file.open() as f:
        run_info = et.parse(f)

    lane_count = int(run_info.find('./Run/FlowcellLayout[@LaneCount]').get('LaneCount'))
    return range(1, lane_count + 1)


# Convert string to boolean
def str2bool(s):
    return s.lower() == "true"
