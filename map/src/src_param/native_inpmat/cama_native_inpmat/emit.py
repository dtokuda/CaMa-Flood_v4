"""Emit runoff-input-specific shell configuration for CaMa-Flood."""

from __future__ import annotations

import json
from pathlib import Path
import shlex
from typing import Any

from .config import (
    configured_output_dir,
    configured_validation_forcing,
    configured_validation_reference_dir,
    resolve_input_files,
)
from .forcing import validate_forcing_collection
from .grid import read_rectilinear_grid
from .inventory import inspect_input


def calendar_to_lleapyr(calendar: str) -> str:
    if calendar in {"365_day", "noleap"}:
        return ".FALSE."
    if calendar in {"gregorian", "proleptic_gregorian", "standard"}:
        return ".TRUE."
    raise ValueError(f"unsupported calendar: {calendar}")


def _shell_assignment(name: str, value: Any) -> str:
    return f"{name}={shlex.quote(str(value))}"


def emit_input_config(config: dict[str, Any]) -> dict[str, Any]:
    output_dir = configured_output_dir(config)
    hash_decimals = config["grid"]["hash_decimals"]
    input_config = configured_validation_forcing(config)
    inventory = inspect_input(config, input_config)
    files = resolve_input_files(input_config)
    reference_file = files[0]
    variable = input_config["variable"]
    metadata = validate_forcing_collection(files, variable)
    grid = read_rectilinear_grid(reference_file, variable)
    grid_hash = grid.hash_at_precision(hash_decimals)
    grid_record = json.loads(
        (configured_validation_reference_dir(config) / "grid.json").read_text()
    )
    if inventory["grid_hash"] != grid_hash or grid_record["grid_hash"] != grid_hash:
        raise ValueError(
            f"{config['output']['dirname']}: validation forcing and mapping grids differ"
        )
    artifact_dir = output_dir
    values = {
        "OUTPUT_DIRNAME": config["output"]["dirname"],
        "RUNOFF_GRID_HASH": grid_hash,
        "INPUT_PATH": input_config["path"],
        "CAMA_CDIMINFO": artifact_dir / "diminfo.txt",
        "CAMA_CINPMAT": artifact_dir / "inpmat.bin",
        "CAMA_LINPCDF": ".TRUE.",
        "CAMA_LINPDAY": ".FALSE.",
        "CAMA_LINTERP": ".TRUE.",
        "CAMA_LITRPCDF": ".FALSE.",
        "CAMA_CVNTIME": "time",
        "CAMA_CVNROF": variable,
        "CAMA_LLEAPYR": calendar_to_lleapyr(str(metadata["calendar"])),
        "CAMA_SYEARIN": 0,
        "CAMA_SMONIN": 0,
        "CAMA_SDAYIN": 0,
        "CAMA_SHOURIN": 0,
        "CAMA_REFERENCE_CROFCDF": reference_file,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    env_path = output_dir / "cama-force.env"
    env_path.write_text(
        "\n".join(_shell_assignment(name, value) for name, value in values.items())
        + "\n"
    )
    report = {
        "config_version": 6,
        "status": "passed",
        "output": {"dirname": config["output"]["dirname"]},
        "grid_hash": grid_hash,
        "calendar": metadata["calendar"],
        "file_count": len(files),
        "env_file": str(env_path),
    }
    return report
