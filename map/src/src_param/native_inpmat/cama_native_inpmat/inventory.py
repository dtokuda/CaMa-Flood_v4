"""Inventory configured input datasets and their horizontal grids."""

from __future__ import annotations

from collections import Counter
import json
from typing import Any

import numpy as np
import xarray as xr

from .config import configured_validation_forcing, resolve_input_files
from .grid import grid_from_dataset


def _json_attribute(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, np.ndarray):
        return _json_attribute(value.tolist())
    if isinstance(value, np.generic):
        return _json_attribute(value.item())
    if isinstance(value, (list, tuple)):
        return [_json_attribute(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_attribute(item) for key, item in value.items()}
    if isinstance(value, float) and not np.isfinite(value):
        if np.isnan(value):
            return "NaN"
        return "Infinity" if value > 0.0 else "-Infinity"
    return value


def _attribute_key(value: Any) -> str:
    return json.dumps(_json_attribute(value), sort_keys=True, allow_nan=False)


def _only_value(values: Counter[str], dirname: str, label: str) -> Any:
    if len(values) != 1:
        raise ValueError(f"{dirname}: {label} changes across input files")
    return json.loads(next(iter(values)))


def inspect_input(
    config: dict[str, Any], input_config: dict[str, Any]
) -> dict[str, Any]:
    """Inspect every file selected by the generic input entry."""

    files = resolve_input_files(input_config)
    dirname = config["output"]["dirname"]
    variable_name = input_config["variable"]
    hash_decimals = config["grid"]["hash_decimals"]
    grid_summaries: dict[str, dict[str, Any]] = {}
    calendars: Counter[str] = Counter()
    time_units: Counter[str] = Counter()
    units: Counter[str] = Counter()
    dtypes: Counter[str] = Counter()
    dimensions: Counter[str] = Counter()
    packing: dict[str, Counter[str]] = {
        name: Counter()
        for name in ("_FillValue", "missing_value", "scale_factor", "add_offset")
    }
    time_steps = 0

    for path in files:
        with xr.open_dataset(
            path, decode_times=False, mask_and_scale=False
        ) as dataset:
            if variable_name not in dataset:
                raise ValueError(f"{path}: missing input variable {variable_name}")
            variable = dataset[variable_name]
            dimensions[_attribute_key(list(variable.dims))] += 1
            units[_attribute_key(variable.attrs.get("units"))] += 1
            dtypes[_attribute_key(str(variable.dtype))] += 1
            for name, values in packing.items():
                values[_attribute_key(variable.attrs.get(name))] += 1

            if "time" in dataset and "time" in variable.dims:
                calendar = dataset["time"].attrs.get("calendar", "standard")
                time_unit = dataset["time"].attrs.get("units")
                time_steps += int(dataset.sizes["time"])
            else:
                calendar = None
                time_unit = None
                time_steps += 1
            calendars[_attribute_key(calendar)] += 1
            time_units[_attribute_key(time_unit)] += 1

            grid = grid_from_dataset(dataset, variable_id=variable_name)
            grid_hash = grid.hash_at_precision(hash_decimals)
            grid_summaries.setdefault(grid_hash, grid.summary(hash_decimals))

    if len(grid_summaries) != 1:
        raise ValueError(
            f"{dirname}: input.path matched {len(grid_summaries)} indexed grids"
        )

    grid_hash, grid_summary = next(iter(grid_summaries.items()))
    return {
        "path": input_config["path"],
        "variable": variable_name,
        "file_count": len(files),
        "first_file": str(files[0]),
        "last_file": str(files[-1]),
        "time_steps": time_steps,
        "calendar": _only_value(calendars, dirname, "calendar"),
        "time_units": _only_value(time_units, dirname, "time units"),
        "units": _only_value(units, dirname, "input units"),
        "dtype": _only_value(dtypes, dirname, "input dtype"),
        "dimensions": _only_value(dimensions, dirname, "input dimensions"),
        "packing": {
            name: _only_value(values, dirname, name)
            for name, values in packing.items()
        },
        "grid_hash": grid_hash,
        "grid": grid_summary,
    }


def inspect_config(config: dict[str, Any]) -> dict[str, Any]:
    input_record = inspect_input(config, config["input"])
    result = {
        "inventory_version": 6,
        "status": "passed",
        "output": {"dirname": config["output"]["dirname"]},
        "unique_grid_count": 1,
        "unique_grid_hashes": [input_record["grid_hash"]],
        "input": input_record,
    }
    if "validation_forcing" in config:
        if config["validation_forcing"].get("use_input") is True:
            result["validation_forcing"] = {
                "use_input": True,
                "grid_hash": input_record["grid_hash"],
            }
            return result
        validation_record = inspect_input(config, configured_validation_forcing(config))
        if validation_record["grid_hash"] != input_record["grid_hash"]:
            raise ValueError(
                f"{config['output']['dirname']}: input and validation forcing grids differ"
            )
        result["validation_forcing"] = validation_record
    return result
