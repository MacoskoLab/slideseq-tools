#!/usr/bin/env python

import pathlib

import xml.etree.ElementTree as et
from typing import List


# Get read structure from RunInfo.xml
def get_read_structure(run_info_file: pathlib.Path):
    # open the RunInfo.xml file and parse it with element tree
    with run_info_file.open() as f:
        run_info = et.parse(f)

    read_elems = run_info.findall("./Run/Reads/Read[@NumCycles][@Number]")
    read_elems.sort(key=lambda el: int(el.get("Number")))

    assert len(read_elems) == 3

    # I guess Picard wants this string format for some reason
    return "{}T{}B{}T".format(*(el.get("NumCycles") for el in read_elems))


def get_tiles(run_info_file: pathlib.Path, lane: str) -> List[str]:
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


# Convert string to boolean
def str2bool(s):
    return s.lower() == "true"
