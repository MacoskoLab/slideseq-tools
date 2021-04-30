import logging
import os
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml

from slideseq.util import constants as constants

log = logging.getLogger(__name__)


@dataclass
class Manifest:
    flowcell: str
    flowcell_directory: Path
    output_directory: Path
    metadata_file: Path
    email_addresses: list[str]

    @staticmethod
    def from_file(input_file: Path):
        with input_file.open() as fh:
            data = yaml.safe_load(fh)

        log.debug(f"Read manifest file {input_file} for flowcell {data['flowcell']}")

        return Manifest(
            flowcell=data["flowcell"],
            flowcell_directory=Path(data["flowcell_directory"]),
            output_directory=Path(data["output_directory"]),
            metadata_file=Path(data["metadata_file"]),
            email_addresses=data["email_addresses"],
        )

    def to_file(self, output_file: Path):
        data = {
            "flowcell": self.flowcell,
            "flowcell_directory": str(self.flowcell_directory),
            "output_directory": str(self.output_directory),
            "metadata_file": str(self.metadata_file),
            "email_addresses": self.email_addresses,
        }

        with output_file.open("w") as out:
            yaml.safe_dump(data, stream=out)

        log.debug(f"Wrote manifest file {output_file} for flowcell {data['flowcell']}")

    @property
    def tmp_dir(self):
        return self.output_directory / "tmp"

    @property
    def log_dir(self):
        return self.output_directory / "logs"


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

    illumina_platform_set = set(flowcell_df.illuminaplatform)
    if len(illumina_platform_set) > 1:
        warning_logs.append(
            f"Flowcell {flowcell} has multiple Illumina platforms associated with it;"
            " please correct to single platform before running."
        )

    # check if platform supported
    if not illumina_platform_set & constants.PLATFORMS:
        warning_logs.append(
            f"Flowcell {flowcell} has unsupported platform; please correct to one of"
            f" {' '.join(constants.PLATFORMS)} before running."
        )

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
