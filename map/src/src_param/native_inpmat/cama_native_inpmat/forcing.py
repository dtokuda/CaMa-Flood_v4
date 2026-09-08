"""Runoff preprocessing and real-field mapping validation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse
import xarray as xr

from .cama import read_binary_inpmat, read_diminfo
from .config import (
    configured_output_dir,
    configured_validation_forcing,
    configured_validation_reference_dir,
    resolve_input_files,
)
from .grid import read_rectilinear_grid
from .inventory import inspect_input
from .mapping import (
    grid_cell_area_factors,
    reference_grid_from_diminfo,
    relative_linf,
)


CAMA_RMIS_FLOAT32 = np.float32(1.0e20)


def preprocess_runoff(
    values: Any, *, missing_values: tuple[float, ...] = ()
) -> tuple[np.ndarray, dict[str, int | float | None]]:
    """Apply the same NaN-to-missing and negative clipping as CaMa forcing input."""

    if np.ma.isMaskedArray(values):
        array = np.asarray(values.filled(np.nan), dtype=np.float32)
    else:
        array = np.asarray(values, dtype=np.float32)
    array = array.copy()
    for marker in missing_values:
        array[array == np.float32(marker)] = np.nan
    missing = ~np.isfinite(array)
    negative = np.isfinite(array) & (array < 0.0)
    processed = np.nan_to_num(array, nan=0.0, posinf=0.0, neginf=0.0)
    processed = np.maximum(processed, np.float32(0.0))
    valid = array[~missing]
    return processed, {
        "missing_count": int(np.count_nonzero(missing)),
        "negative_count": int(np.count_nonzero(negative)),
        "valid_min": float(np.min(valid)) if valid.size else None,
        "valid_max": float(np.max(valid)) if valid.size else None,
        "processed_sum": float(np.sum(processed, dtype=np.float64)),
    }


def read_preprocessed_timestep(
    path: str | Path, variable_id: str, time_index: int
) -> tuple[np.ndarray, dict[str, int | float | None]]:
    with xr.open_dataset(path, decode_times=False, mask_and_scale=False) as dataset:
        variable = dataset[variable_id]
        values = variable.isel(time=time_index).values
        markers = tuple(
            float(variable.attrs[name])
            for name in ("_FillValue", "missing_value")
            if variable.attrs.get(name) is not None
        )
    return preprocess_runoff(values, missing_values=markers)


def forcing_metadata(path: str | Path, variable_id: str = "mrro") -> dict[str, Any]:
    """Inspect metadata that controls the direct CaMa NetCDF read."""

    with xr.open_dataset(
        path, decode_times=False, mask_and_scale=False
    ) as dataset:
        variable = dataset[variable_id]
        fill = variable.attrs.get("_FillValue")
        missing = variable.attrs.get("missing_value")
        if fill is not None and np.float32(fill) != CAMA_RMIS_FLOAT32:
            raise ValueError(f"{path}: _FillValue is not CaMa RMIS=1.E20 float32")
        if missing is not None and np.float32(missing) != CAMA_RMIS_FLOAT32:
            raise ValueError(f"{path}: missing_value is not CaMa RMIS=1.E20 float32")
        if variable.attrs.get("scale_factor") is not None:
            raise ValueError(f"{path}: packed scale_factor is not supported directly")
        if variable.attrs.get("add_offset") is not None:
            raise ValueError(f"{path}: packed add_offset is not supported directly")
        if variable.dims != ("time", "lat", "lon"):
            raise ValueError(f"{path}: expected {variable_id}(time, lat, lon)")
        return {
            "dtype": str(variable.dtype),
            "shape": [int(dataset.sizes[name]) for name in variable.dims],
            "units": variable.attrs.get("units"),
            "_FillValue": float(np.float32(fill)) if fill is not None else None,
            "missing_value": (
                float(np.float32(missing)) if missing is not None else None
            ),
            "scale_factor": None,
            "add_offset": None,
            "calendar": dataset["time"].attrs.get("calendar", "standard"),
            "time_units": dataset["time"].attrs.get("units"),
        }


def validate_forcing_collection(
    paths: list[Path], variable_id: str
) -> dict[str, Any]:
    """Check that every matched file is directly readable by the CaMa workflow."""

    records = [forcing_metadata(path, variable_id) for path in paths]
    return records[0]


def validate_real_forcing(config: dict[str, Any]) -> dict[str, Any]:
    """Compare explicit-reference and composed mappings on real runoff fields."""

    reference_config = config["reference_mapping"]
    base_diminfo = read_diminfo(reference_config["diminfo_path"])
    standard = read_binary_inpmat(
        reference_config["inpmat_path"], base_diminfo
    ).to_csr().astype(np.float32)
    output_dir = configured_output_dir(config)
    hash_decimals = config["grid"]["hash_decimals"]

    input_config = configured_validation_forcing(config)
    inventory = inspect_input(config, input_config)
    paths = resolve_input_files(input_config)
    path = paths[0]
    metadata = validate_forcing_collection(paths, input_config["variable"])
    grid = read_rectilinear_grid(path, input_config["variable"])
    grid_hash = grid.hash_at_precision(hash_decimals)
    artifact_dir = output_dir
    validation_reference_dir = configured_validation_reference_dir(config)
    grid_record = json.loads((validation_reference_dir / "grid.json").read_text())
    if inventory["grid_hash"] != grid_hash or grid_record["grid_hash"] != grid_hash:
        raise ValueError(
            f"{config['output']['dirname']}: validation forcing and mapping grids differ"
        )
    remap = sparse.load_npz(
        validation_reference_dir / "native_to_reference_flux.npz"
    ).tocsr()
    diminfo = read_diminfo(artifact_dir / "diminfo.txt")
    emitted = read_binary_inpmat(
        artifact_dir / "inpmat.bin", diminfo
    ).to_csr().astype(np.float32)
    target_grid = reference_grid_from_diminfo(
        base_diminfo, reference_config["input_grid"]
    )
    source_area = grid_cell_area_factors(grid)
    target_area = grid_cell_area_factors(target_grid)

    timestep_reports = []
    samples = ((paths[0], 0), (paths[-1], -1))
    for sample_path, time_index in samples:
        runoff, statistics = read_preprocessed_timestep(
            sample_path, input_config["variable"], time_index
        )
        native = runoff.ravel()
        reference_float64 = remap @ native.astype(np.float64)
        source_volume = float(np.dot(source_area, native.astype(np.float64)))
        target_volume = float(np.dot(target_area, reference_float64))
        mass_relative_error = abs(target_volume - source_volume) / max(
            abs(source_volume), np.finfo(float).tiny
        )
        reference_field = reference_float64.astype(np.float32)
        explicit_result = standard @ reference_field
        composed_result = emitted @ native
        explicit_sum = float(np.sum(explicit_result, dtype=np.float64))
        composed_sum = float(np.sum(composed_result, dtype=np.float64))
        mapped_total_relative_error = abs(composed_sum - explicit_sum) / max(
            abs(explicit_sum), np.finfo(float).tiny
        )
        timestep_reports.append(
            {
                "file": str(sample_path),
                "time_index": time_index,
                "input": statistics,
                "runtime_relative_linf": relative_linf(
                    composed_result, explicit_result
                ),
                "runtime_max_abs": float(
                    np.max(np.abs(composed_result - explicit_result), initial=0.0)
                ),
                "actual_native_to_reference_mass_relative_error": mass_relative_error,
                "global_mapped_sum_explicit": explicit_sum,
                "global_mapped_sum_composed": composed_sum,
                "mapped_total_relative_error": mapped_total_relative_error,
            }
        )

    thresholds = {
        "runtime_relative_linf": 1.0e-6,
        "actual_native_to_reference_mass_relative_error": 1.0e-12,
        "mapped_total_relative_error": 1.0e-6,
    }
    failures = {
        name: max(float(item[name]) for item in timestep_reports)
        for name, threshold in thresholds.items()
        if max(float(item[name]) for item in timestep_reports) > threshold
    }
    if failures:
        raise ValueError(
            f"{config['output']['dirname']}: real-field mapping errors {failures}"
        )
    report = {
        "validation_version": 5,
        "output": {"dirname": config["output"]["dirname"]},
        "first_file": str(paths[0]),
        "last_file": str(paths[-1]),
        "file_count": len(paths),
        "grid_hash": grid_hash,
        "forcing_metadata": metadata,
        "thresholds": thresholds,
        "timesteps": timestep_reports,
        "status": "passed",
    }
    return report


def forcing_start_date(
    path: str | Path, variable_id: str = "mrro"
) -> tuple[int, int, int]:
    """Read the first forcing date from the NetCDF time coordinate."""

    with xr.open_dataset(path, decode_times=True, mask_and_scale=False) as dataset:
        if variable_id not in dataset or "time" not in dataset[variable_id].dims:
            raise ValueError(f"{path}: {variable_id} has no time dimension")
        value = dataset["time"].values[0]
    if isinstance(value, np.datetime64):
        date_text = np.datetime_as_string(value, unit="D")
        year, month, day = date_text.split("-")
        return int(year), int(month), int(day)
    return int(value.year), int(value.month), int(value.day)
