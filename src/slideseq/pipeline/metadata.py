import datetime
import pathlib
import typing
from dataclasses import dataclass


@dataclass
class Manifest:
    flowcell_directory: pathlib.Path
    output_directory: pathlib.Path
    metadata_file: pathlib.Path
    flowcell: str
    experiment_date: datetime.date
    is_novaseq: bool
    email_addresses: typing.Sequence[str]


@dataclass
class Metadata:
    ...
