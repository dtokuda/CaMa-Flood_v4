"""Configuration loading and validation."""

from __future__ import annotations

import glob
import os
from pathlib import Path
from typing import Any

import yaml


ROOT_KEYS = {
    "version",
    "cama",
    "reference_mapping",
    "grid",
    "preprocess",
    "output",
    "input",
    "validation_forcing",
}
CAMA_KEYS = {"executable", "map_dir"}
REFERENCE_MAPPING_KEYS = {"diminfo_path", "inpmat_path", "input_grid"}
REFERENCE_GRID_KEYS = {
    "type",
    "longitude_order",
    "latitude_order",
}
GRID_KEYS = {"type", "coverage", "bounds_policy", "hash_decimals"}
PREPROCESS_KEYS = {
    "clip_negative",
    "use_areacella",
    "use_sftlf",
    "normalize_over_valid_area",
}
OUTPUT_KEYS = {"dirname"}
INPUT_KEYS = {"path", "variable"}
VALIDATION_FORCING_KEYS = {"use_input", "path", "variable"}


def _reject_unknown_keys(
    value: Any, allowed: set[str], section: str
) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{section} must be a mapping")
    unknown = set(value).difference(allowed)
    if unknown:
        raise ValueError(
            f"unknown {section} keys: {', '.join(sorted(unknown))}"
        )
    return value


def _resolve_input_entry(
    value: Any, config_dir: Path, section: str
) -> dict[str, Any]:
    entry = _reject_unknown_keys(value, INPUT_KEYS, section)
    missing = INPUT_KEYS.difference(entry)
    if missing:
        raise ValueError(f"{section} is missing keys: {', '.join(sorted(missing))}")
    entry["path"] = _resolve_pattern(str(entry["path"]), config_dir)
    variable = str(entry["variable"])
    if not variable:
        raise ValueError(f"{section}.variable must not be empty")
    entry["variable"] = variable
    return entry


def _resolve_path(value: str, config_dir: Path) -> str:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = config_dir / path
    return str(path.resolve())


def _resolve_pattern(value: str, config_dir: Path) -> str:
    expanded = os.path.expanduser(value)
    if not os.path.isabs(expanded):
        expanded = os.path.join(config_dir, expanded)
    return os.path.abspath(expanded)


def load_config(path: str | Path) -> dict[str, Any]:
    """Load the composed input-mapping configuration."""

    config_path = Path(path).expanduser().resolve()
    with config_path.open("r", encoding="utf-8") as stream:
        config = yaml.safe_load(stream)

    if not isinstance(config, dict):
        raise ValueError("configuration root must be a mapping")
    _reject_unknown_keys(config, ROOT_KEYS, "configuration")
    if config.get("version") != 6:
        raise ValueError("unsupported configuration version; expected version 6")

    for key in (
        "cama",
        "reference_mapping",
        "grid",
        "preprocess",
        "output",
        "input",
    ):
        if key not in config:
            raise ValueError(f"missing configuration key: {key}")

    config_dir = config_path.parent
    cama = _reject_unknown_keys(config["cama"], CAMA_KEYS, "cama")
    for key in ("executable", "map_dir"):
        if key not in cama:
            raise ValueError(f"missing cama configuration key: {key}")
        cama[key] = _resolve_path(str(cama[key]), config_dir)

    reference = _reject_unknown_keys(
        config["reference_mapping"],
        REFERENCE_MAPPING_KEYS,
        "reference_mapping",
    )
    for key in ("diminfo_path", "inpmat_path", "input_grid"):
        if key not in reference:
            raise ValueError(f"missing reference_mapping key: {key}")
    for key in ("diminfo_path", "inpmat_path"):
        reference[key] = _resolve_path(str(reference[key]), Path(cama["map_dir"]))
        if not Path(reference[key]).is_file():
            raise FileNotFoundError(f"reference {key} does not exist: {reference[key]}")
    reference_grid = _reject_unknown_keys(
        reference["input_grid"],
        REFERENCE_GRID_KEYS,
        "reference_mapping.input_grid",
    )
    if reference_grid.get("type") != "global_regular_latlon":
        raise ValueError("reference input_grid.type must be global_regular_latlon")
    if reference_grid.get("longitude_order") not in {
        "west_to_east",
        "east_to_west",
    }:
        raise ValueError("invalid reference longitude_order")
    if reference_grid.get("latitude_order") not in {
        "north_to_south",
        "south_to_north",
    }:
        raise ValueError("invalid reference latitude_order")

    grid = _reject_unknown_keys(config["grid"], GRID_KEYS, "grid")
    if grid.get("type") != "rectilinear":
        raise ValueError("only rectilinear source grids are supported")
    if grid.get("coverage") != "global":
        raise ValueError("grid.coverage must be global")
    if grid.get("bounds_policy") != "cf_or_derive":
        raise ValueError("bounds_policy must be cf_or_derive")
    decimals = int(grid.get("hash_decimals", 7))
    if not 4 <= decimals <= 12:
        raise ValueError("hash_decimals must be between 4 and 12")
    grid["hash_decimals"] = decimals

    preprocess = _reject_unknown_keys(
        config["preprocess"], PREPROCESS_KEYS, "preprocess"
    )
    required_false = (
        "use_areacella",
        "use_sftlf",
        "normalize_over_valid_area",
    )
    for key in required_false:
        if preprocess.get(key) is not False:
            raise ValueError(f"{key} must be false for this workflow")
    if preprocess.get("clip_negative") is not True:
        raise ValueError("clip_negative must be true")

    output = _reject_unknown_keys(config["output"], OUTPUT_KEYS, "output")
    if "dirname" not in output:
        raise ValueError("output is missing key: dirname")
    dirname = str(output["dirname"])
    if not dirname:
        raise ValueError("output.dirname must not be empty")
    if Path(dirname).name != dirname or dirname in {".", ".."}:
        raise ValueError("output.dirname must be a single directory name")
    output["dirname"] = dirname

    config["input"] = _resolve_input_entry(config["input"], config_dir, "input")

    if "validation_forcing" in config:
        validation = _reject_unknown_keys(
            config["validation_forcing"],
            VALIDATION_FORCING_KEYS,
            "validation_forcing",
        )
        use_input = validation.get("use_input")
        if use_input is True:
            if set(validation) != {"use_input"}:
                raise ValueError(
                    "validation_forcing.use_input cannot be combined with path or variable"
                )
        elif use_input is not None:
            raise ValueError("validation_forcing.use_input must be true when specified")
        else:
            config["validation_forcing"] = _resolve_input_entry(
                validation, config_dir, "validation_forcing"
            )
    return config


def configured_output_dir(config: dict[str, Any]) -> Path:
    """Return the configured human-readable output directory."""

    return (
        Path(config["cama"]["map_dir"])
        / "inpmat_comp-remap"
        / config["output"]["dirname"]
    )


def configured_validation_reference_dir(config: dict[str, Any]) -> Path:
    """Return the directory containing files needed only for validation."""

    return configured_output_dir(config) / "validation_reference"


def configured_validation_forcing(config: dict[str, Any]) -> dict[str, Any]:
    """Return the explicitly configured runoff used by forcing validations."""

    if "validation_forcing" not in config:
        raise ValueError(
            "this command requires validation_forcing; mapping generation itself "
            "only requires input"
        )
    value = config["validation_forcing"]
    return config["input"] if value.get("use_input") is True else value


def resolve_input_files(input_config: dict[str, Any]) -> list[Path]:
    """Resolve the configured input file or glob pattern."""

    pattern = input_config["path"]
    files = sorted(
        Path(value).resolve()
        for value in glob.glob(pattern, recursive=True)
        if Path(value).is_file()
    )
    if not files:
        raise FileNotFoundError(
            f"no files matched input.path={pattern}"
        )
    return files
