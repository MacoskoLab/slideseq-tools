import datetime
import os
from pathlib import Path
from dataclasses import dataclass

import pandas as pd
import yaml
from slideseq.pipeline.submit_slideseq import log
from slideseq.util import constants as constants


@dataclass
class Manifest:
    flowcell_directory: Path
    output_directory: Path
    metadata_file: Path
    flowcell: str
    experiment_date: datetime.date
    is_novaseq: bool
    email_addresses: list[str]

    @staticmethod
    def from_file(input_file: Path):
        with input_file.open() as fh:
            data = yaml.safe_load_all(fh)

        return Manifest(
            flowcell_directory=Path(data["flowcell_directory"]),
            output_directory=Path(data["output_directory"]),
            metadata_file=Path(data["metadata_file"]),
            flowcell=data["flowcell"],
            experiment_date=datetime.date(data["experiment_date"]),
            is_novaseq=data["is_novaseq"],
            email_addresses=data["email_addresses"],
        )

    def to_file(self, output_file: Path):
        data = {
            "flowcell_directory": self.flowcell_directory,
            "output_directory": self.output_directory,
            "metadata_file": self.metadata_file,
            "flowcell": self.flowcell,
            "experiment_date": self.experiment_date,
            "is_novaseq": self.is_novaseq,
            "email_addresses": self.email_addresses,
        }

        with output_file.open("w") as out:
            yaml.safe_dump(data, stream=out)


@dataclass
class Metadata:
    ...


def validate_flowcell_df(flowcell: str, flowcell_df: pd.DataFrame) -> bool:
    warning_logs = []

    if len(set(flowcell_df.library)) != len(flowcell_df.library):
        warning_logs.append(
            f"Flowcell {flowcell} does not have unique library names;"
            f" please fill out naming metadata correctly or add suffixes."
        )

    if any(" " in name for name in flowcell_df.library):
        warning_logs.append(
            f"The 'library' column for {flowcell} contains spaces;"
            " please remove all spaces before running."
        )

    if any(flowcell_df.isna()):
        warning_logs.append(
            f"Flowcell {flowcell} does not have complete submission metadata (orange and blue cols);"
            " please fill out before running."
        )

    if len(set(flowcell_df.BCLPath)) > 1:
        warning_logs.append(
            f"Flowcell {flowcell} has multiple BCLPaths associated with it;"
            " please correct to single path before running."
        )

    if not os.path.isdir(flowcell_df.BCLPath.values[0]):
        warning_logs.append(
            f"Flowcell {flowcell} has incorrect BCLPath associated with it;"
            " please correct path before running."
        )

    if len(set(flowcell_df.IlluminaPlatform)) > 1:
        warning_logs.append(
            f"Flowcell {flowcell} has multiple Illumina platforms associated with it;"
            " please correct to single platform before running."
        )

    # check if platform supported
    if flowcell_df.IlluminaPlatform[0] not in constants.PLATFORMS:
        warning_logs.append(
            f"Flowcell {flowcell} has unsupported platform; please correct to one of"
            f" {' '.join(constants.PLATFORMS)} before running."
        )

    # check if references exist
    if not all(os.path.isfile(build) for build in flowcell_df.reference):
        warning_logs.append(
            f"Reference for {flowcell} does not exist; please correct reference values before running.",
        )

    if warning_logs:
        for msg in warning_logs:
            log.warning(msg)

        return False
    else:
        return True