import importlib.resources
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

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
    gsecret_name: str
    gsheet_id: str
    worksheet: str
    gs_path: Optional[Path]

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
            gsecret_name=data["gsecret_name"],
            gsheet_id=data["gsheet_id"],
            worksheet=data["worksheet"],
            gs_path=Path(data["gs_path"]) if data["gs_path"] else None,
        )

    def dropseq_cmd(
        self,
        command: str,
        input_file: Union[Path, str],
        output_file: Union[Path, str],
        tmp_dir: Path,
        mem: str = "8g",
        compression: int = 0,
    ):
        """Return the beginning of a DropSeq command, with standard options

        :param command: name of the dropseq tool being invoked
        :param input_file: path to the input file
        :param output_file: path to the output file
        :param tmp_dir: Location of the tmp directory to use
        :param mem: memory for the heap. default is to share with other jobs
        :param compression: compression level for output. Use 0 for speed, 5 for storage
        """

        return [
            self.dropseq_dir / command,
            "-m",
            mem,
            f"I={input_file}",
            f"O={output_file}",
            f"TMP_DIR={tmp_dir}",
            "VALIDATION_STRINGENCY=SILENT",
            f"COMPRESSION_LEVEL={compression}",
            "VERBOSITY=WARNING",
            "QUIET=true",
        ]

    def picard_cmd(self, command: str, tmp_dir: Path, mem: str = "62g"):
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
            self.picard,
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


def get_config() -> Config:
    with importlib.resources.path(slideseq, "config.yaml") as config_path:
        return Config.from_file(config_path)
