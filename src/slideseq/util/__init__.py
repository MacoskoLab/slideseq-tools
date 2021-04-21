#!/usr/bin/env python

import pathlib
import sys
import xml.etree.ElementTree as et

import slideseq.util.constants as constants


def get_env_name():
    # need to pass this to the qsub scripts so they activate the right thing
    return sys.executable.split("/")[-3]


def qsub_args(cpu: int = 1, mem: int = 1, hours: int = 1, minutes: int = 0) -> list[str]:
    """
    Returns a list of options to pass to qsub, including the constant options

    :param cpu: number of cpus to request
    :param mem: amount of memory (in gb) to request
    :param hours: expected hours needed for the job
    :param minutes: expected minutes needed for the job
    :return: configured list of options, including common args
    """

    arg_list = constants.QSUB_ARGS[:]
    if cpu > 1:
        arg_list.extend(["-pe", "smp", f"{cpu}", "-binding", f"linear:{cpu}"])

    arg_list.extend(["-l", f"h_vmem={mem}g", "-l", f"h_rt={hours}:{minutes}:0"])

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
    return range(lane_count)


# Convert string to boolean
def str2bool(s):
    return s.lower() == "true"
