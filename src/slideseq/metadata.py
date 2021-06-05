import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd
import yaml

import slideseq.util.constants as constants
from slideseq.library import Library, LibraryLane, Reference

log = logging.getLogger(__name__)


@dataclass
class Manifest:
    flowcell: str
    flowcell_dir: Path
    workflow_dir: Path
    library_dir: Path
    metadata_file: Path
    email_addresses: list[str]

    @staticmethod
    def from_file(input_file: Path):
        with input_file.open() as fh:
            data = yaml.safe_load(fh)

        log.debug(f"Read manifest file {input_file} for flowcell {data['flowcell']}")

        return Manifest(
            flowcell=data["flowcell"],
            flowcell_dir=Path(data["flowcell_dir"]),
            workflow_dir=Path(data["workflow_dir"]),
            library_dir=Path(data["library_dir"]),
            metadata_file=Path(data["metadata_file"]),
            email_addresses=data["email_addresses"],
        )

    def to_file(self, output_file: Path):
        data = {
            "flowcell": self.flowcell,
            "flowcell_dir": str(self.flowcell_dir),
            "workflow_dir": str(self.workflow_dir),
            "library_dir": str(self.library_dir),
            "metadata_file": str(self.metadata_file),
            "email_addresses": self.email_addresses,
        }

        with output_file.open("w") as out:
            yaml.safe_dump(data, stream=out)

        log.debug(f"Wrote manifest file {output_file} for flowcell {data['flowcell']}")

    def _get_library(self, library_index: int):
        metadata_df = pd.read_csv(self.metadata_file)

        library_name = sorted(set(metadata_df.library))[library_index]
        library_df = metadata_df.loc[metadata_df.library == library_name]

        # verify that all columns are consistent, except for lane
        validate_library_df(library_name, library_df)
        row = library_df.iloc[0]

        lanes = sorted(set(library_df.lane))

        return library_name, row, lanes, Reference(row.reference)

    def get_library(self, library_index: int) -> Library:
        return Library(*self._get_library(library_index), library_dir=self.library_dir)

    def get_library_lane(self, library_index: int, lane: int) -> Optional[LibraryLane]:
        library_name, row, lanes, reference = self._get_library(library_index)

        if lane not in lanes:
            log.warning("Library not present in this lane, nothing to do here")
            return None

        return LibraryLane(library_name, row, lanes, reference, self.library_dir, lane)

    @property
    def tmp_dir(self):
        return self.workflow_dir / "tmp"

    @property
    def log_dir(self):
        return self.workflow_dir / "logs"


def validate_flowcell_df(flowcell: str, flowcell_df: pd.DataFrame) -> bool:
    warning_logs = []

    if len(set(flowcell_df.library)) != len(flowcell_df.library):
        warning_logs.append(
            f"Flowcell {flowcell} does not have unique library names;"
            " please fill out naming metadata correctly or add suffixes."
        )

    if any(" " in name for name in flowcell_df.library):
        warning_logs.append(
            f"The 'library' column for {flowcell} contains spaces;"
            " please remove all spaces before running."
        )

    if flowcell_df.isna().values.any():
        warning_logs.append(
            f"Flowcell {flowcell} does not have complete submission metadata (orange"
            " and blue cols); please fill out before running."
        )

    bcl_path_set = set(flowcell_df.bclpath)
    if len(bcl_path_set) > 1:
        warning_logs.append(
            f"Flowcell {flowcell} has multiple BCLPaths associated with it;"
            " please correct to single path before running."
        )

    bcl_path = Path(bcl_path_set.pop())
    if not bcl_path.is_dir():
        warning_logs.append(
            f"Flowcell {flowcell} has incorrect BCLPath associated with it;"
            " please correct path before running."
        )

    run_info_file = bcl_path / "RunInfo.xml"
    if not run_info_file.exists():
        warning_logs.append(f"{run_info_file} for flowcell {flowcell} does not exist")

    # check if references exist
    if not all(os.path.isfile(build) for build in flowcell_df.reference):
        warning_logs.append(
            f"Reference for {flowcell} does not exist; please correct reference values"
            " before running.",
        )

    # for runs that need barcode matching, check for the existence of necessary files
    for _, row in flowcell_df.iterrows():
        if row.run_barcodematching:
            puck_dir = Path(row.puckcaller_path)
            if not (puck_dir / "BeadBarcodes.txt").exists():
                warning_logs.append(f"BeadBarcodes.txt not found in {puck_dir}")
            if not (puck_dir / "BeadLocations.txt").exists():
                warning_logs.append(f"BeadLocations.txt not found in {puck_dir}")

    if warning_logs:
        for msg in warning_logs:
            log.warning(msg)

        return False
    else:
        return True


def split_sample_lanes(flowcell_df: pd.DataFrame, lanes: range):
    new_rows = []

    for _, row in flowcell_df.iterrows():
        if row.lane == "{LANE}":
            row_lanes = lanes
        else:
            row_lanes = str(row.lane).split(",")

        for lane in row_lanes:
            row.lane = int(lane)
            new_rows.append(row.copy())

    return pd.DataFrame(new_rows, index=range(len(new_rows)))


def validate_library_df(library_name: str, library_df: pd.DataFrame):
    """Verify that all of the columns in the dataframe are constant, except
    for the lane which was expanded out earlier"""

    for col in constants.METADATA_COLS:
        if col.lower() == "lane":
            continue

        if len(set(library_df[col.lower()])) != 1:
            raise ValueError(
                f"Library {library_name} has multiple values in column {col}"
            )
