import importlib.resources
import logging
from dataclasses import dataclass
from pathlib import Path

import yaml

import slideseq

log = logging.getLogger(__name__)


@dataclass
class Config:
    picard: Path
    dropseq_dir: Path
    reference_dir: Path
    workflow_dir: Path
    library_dir: Path

    @staticmethod
    def from_file(input_file: Path):
        with input_file.open() as fh:
            data = yaml.safe_load(fh)

        log.debug(f"Read config file {input_file}")

        return Config(
            picard=Path(data["picard"]),
            dropseq_dir=Path(data["dropseq_dir"]),
            reference_dir=Path(data["reference_dir"]),
            workflow_dir=Path(data["workflow_dir"]),
            library_dir=Path(data["library_dir"]),
        )


def get_config() -> Config:
    with importlib.resources.path(slideseq, "config.yaml") as config_path:
        return Config.from_file(config_path)
