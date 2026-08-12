"""Short CaMa-Flood integration test for explicit and composed mappings."""

from __future__ import annotations

from datetime import date, timedelta
import json
import hashlib
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any

import numpy as np
from scipy import sparse
import xarray as xr

from .cama import read_diminfo
from .config import (
    configured_output_dir,
    configured_validation_forcing,
    configured_validation_reference_dir,
    resolve_input_files,
)
from .emit import calendar_to_lleapyr
from .forcing import (
    forcing_metadata,
    forcing_start_date,
    read_preprocessed_timestep,
)
from .grid import read_rectilinear_grid
from .inventory import inspect_input
from .mapping import reference_grid_from_diminfo, relative_linf


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tool_source_sha256() -> str:
    """Hash Python source that defines generation and validation behavior."""

    root = Path(__file__).resolve().parents[1]
    paths = sorted((root / "cama_native_inpmat").glob("*.py"))
    paths.extend(sorted(root.glob("*.py")))
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path.relative_to(root)).encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def _read_native_fields(path: Path, variable_id: str, days: int) -> np.ndarray:
    with xr.open_dataset(path, decode_times=False, mask_and_scale=False) as dataset:
        if dataset.sizes["time"] < days:
            raise ValueError(f"{path}: fewer than {days} timesteps")
    fields = [
        read_preprocessed_timestep(path, variable_id, index)[0]
        for index in range(days)
    ]
    return np.stack(fields)


def _write_explicit_forcing(
    path: Path,
    native_fields: np.ndarray,
    remap: sparse.csr_matrix,
    target_grid,
    variable_name: str,
    units: str | None,
    start: tuple[int, int, int],
    calendar: str,
) -> np.ndarray:
    fields = np.stack(
        [
            (remap @ field.ravel().astype(np.float64))
            .astype(np.float32)
            .reshape(target_grid.shape)
            for field in native_fields
        ]
    )
    year, month, day = start
    dataset = xr.Dataset(
        data_vars={
            variable_name: (
                ("time", "lat", "lon"),
                fields,
                {"units": units or ""},
            )
        },
        coords={
            "time": (
                "time",
                np.arange(fields.shape[0], dtype=np.float64),
                {
                    "units": f"days since {year:04d}-{month:02d}-{day:02d} 00:00:00",
                    "calendar": calendar,
                },
            ),
            "lat": ("lat", target_grid.latitude),
            "lon": ("lon", target_grid.longitude),
        },
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_netcdf(
        path,
        encoding={
            variable_name: {"dtype": "float32", "_FillValue": np.float32(1.0e20)}
        },
    )
    return fields


def _namelist(
    *,
    map_dir: Path,
    diminfo: Path,
    inpmat: Path,
    forcing: Path,
    start: tuple[int, int, int],
    days: int,
    calendar: str,
    variable_name: str,
    forcing_start: tuple[int, int, int] | None = None,
    restart_path: Path | None = None,
) -> str:
    year, month, day = start
    end = date(year, month, day) + timedelta(days=days)
    lleapyr = calendar_to_lleapyr(calendar)
    restart_enabled = ".TRUE." if restart_path else ".FALSE."
    restart_value = str(restart_path) if restart_path else ""
    if forcing_start is None:
        forcing_time_namelist = "! forcing start time is read from NetCDF CF metadata"
    else:
        forcing_year, forcing_month, forcing_day = forcing_start
        forcing_time_namelist = f"""SYEARIN  = {forcing_year}
SMONIN   = {forcing_month}
SDAYIN   = {forcing_day}
SHOURIN  = 0"""
    return f"""&NRUNVER
LADPSTP  = .TRUE.
LPTHOUT  = .FALSE.
LDAMOUT  = .FALSE.
LRESTART = {restart_enabled}
LOUTPUT  = .TRUE.
LGRIDMAP = .TRUE.
LLEAPYR  = {lleapyr}
/
&NDIMTIME
CDIMINFO = "{diminfo}"
DT       = 86400
IFRQ_INP = 24
/
&NPARAM
PMANRIV  = 0.03D0
PMANFLD  = 0.10D0
PDSTMTH  = 10000.D0
PCADP    = 0.7
/
&NSIMTIME
SYEAR  = {year}
SMON   = {month}
SDAY   = {day}
SHOUR  = 0
EYEAR  = {end.year}
EMON   = {end.month}
EDAY   = {end.day}
EHOUR  = 0
/
&NMAP
LMAPCDF = .FALSE.
CNEXTXY = "{map_dir / 'nextxy.bin'}"
CGRAREA = "{map_dir / 'ctmare.bin'}"
CELEVTN = "{map_dir / 'elevtn.bin'}"
CNXTDST = "{map_dir / 'nxtdst.bin'}"
CRIVLEN = "{map_dir / 'rivlen.bin'}"
CFLDHGT = "{map_dir / 'fldhgt.bin'}"
CRIVWTH = "{map_dir / 'rivwth_gwdlr.bin'}"
CRIVHGT = "{map_dir / 'rivhgt.bin'}"
CRIVMAN = "{map_dir / 'rivman.bin'}"
CPTHOUT = "{map_dir / 'bifprm.txt'}"
/
&NRESTART
CRESTSTO = "{restart_value}"
CRESTDIR = "./"
CVNREST  = "restart"
LRESTCDF = .FALSE.
IFRQ_RST = 0
/
&NFORCE
LINPCDF  = .TRUE.
LINPDAY  = .FALSE.
LINTERP  = .TRUE.
LITRPCDF = .FALSE.
CINPMAT  = "{inpmat}"
CROFCDF  = "{forcing}"
CVNTIME  = "time"
CVNROF   = "{variable_name}"
{forcing_time_namelist}
/
&NOUTPUT
COUTDIR  = "./"
CVARSOUT = "runoff,rivsto,storge,rivout"
COUTTAG  = "_integration"
LOUTVEC  = .FALSE.
LOUTCDF  = .FALSE.
IFRQ_OUT = 24
/
"""


def _run_case(executable: Path, case_dir: Path, namelist: str) -> None:
    case_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "runoff_integration.bin",
        "rivsto_integration.bin",
        "storge_integration.bin",
        "rivout_integration.bin",
        "log_CaMa.txt",
        "stdout.log",
    ):
        path = case_dir / name
        if path.exists():
            path.unlink()
    (case_dir / "input_cmf.nam").write_text(namelist)
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    completed = subprocess.run(
        [str(executable)],
        cwd=case_dir,
        env=environment,
        text=True,
        capture_output=True,
        timeout=300,
        check=False,
    )
    log = completed.stdout + completed.stderr
    (case_dir / "stdout.log").write_text(log)
    if completed.returncode != 0:
        cama_log_path = case_dir / "log_CaMa.txt"
        cama_log = cama_log_path.read_text(errors="replace") if cama_log_path.is_file() else ""
        raise RuntimeError(
            f"CaMa failed in {case_dir} with code {completed.returncode}:\n"
            f"{(log + cama_log)[-8000:]}"
        )


def _read_output(path: Path, days: int, ny: int, nx: int) -> np.ndarray:
    values = np.fromfile(path, dtype="<f4")
    expected = days * ny * nx
    if values.size != expected:
        raise ValueError(f"{path}: found {values.size} values, expected {expected}")
    return values.reshape(days, ny, nx)


def _compare_outputs(explicit: np.ndarray, composed: np.ndarray) -> dict[str, float | int]:
    valid = np.isfinite(explicit) & np.isfinite(composed) & (explicit < 5.0e19) & (
        composed < 5.0e19
    )
    if not np.any(valid):
        raise ValueError("CaMa output comparison has no valid cells")
    difference = np.abs(composed[valid] - explicit[valid])
    absolute_expected = np.abs(explicit[valid])
    return {
        "valid_values": int(np.count_nonzero(valid)),
        "relative_linf": relative_linf(composed[valid], explicit[valid]),
        "relative_l1": float(
            np.sum(difference, dtype=np.float64)
            / max(np.sum(absolute_expected, dtype=np.float64), np.finfo(float).tiny)
        ),
        "max_abs": float(np.max(difference)),
    }


def _run_cama_integration_in_directory(
    config: dict[str, Any], integration_dir: Path, days: int
) -> dict[str, Any]:
    """Run CaMa twice and compare explicit-reference and composed paths."""

    input_config = configured_validation_forcing(config)
    dirname = config["output"]["dirname"]
    variable_name = input_config["variable"]
    forcing = resolve_input_files(input_config)[0]
    inventory = inspect_input(config, input_config)
    metadata = forcing_metadata(forcing, variable_name)
    calendar = str(metadata["calendar"])
    start = forcing_start_date(forcing, variable_name)
    source_grid = read_rectilinear_grid(forcing, variable_name)
    grid_hash = source_grid.hash_at_precision(config["grid"]["hash_decimals"])
    grid_record = json.loads(
        (configured_validation_reference_dir(config) / "grid.json").read_text()
    )
    if inventory["grid_hash"] != grid_hash or grid_record["grid_hash"] != grid_hash:
        raise ValueError(
            f"{dirname}: validation forcing and mapping grids differ"
        )

    cama_config = config["cama"]
    executable = Path(cama_config["executable"])
    if not executable.is_file():
        raise FileNotFoundError(f"CaMa executable not found: {executable}")
    map_dir = Path(cama_config["map_dir"])
    artifact_dir = configured_output_dir(config)
    explicit_forcing = integration_dir / "explicit_reference.nc"

    remap = sparse.load_npz(
        configured_validation_reference_dir(config) / "native_to_reference_flux.npz"
    ).tocsr()
    native_fields = _read_native_fields(forcing, variable_name, days)
    reference_config = config["reference_mapping"]
    base_diminfo = read_diminfo(reference_config["diminfo_path"])
    target_grid = reference_grid_from_diminfo(
        base_diminfo, reference_config["input_grid"]
    )
    _write_explicit_forcing(
        explicit_forcing,
        native_fields,
        remap,
        target_grid,
        variable_name,
        metadata["units"],
        start,
        calendar,
    )

    explicit_case = integration_dir / "explicit"
    composed_case = integration_dir / "composed"
    explicit_namelist = _namelist(
        map_dir=map_dir,
        diminfo=Path(reference_config["diminfo_path"]),
        inpmat=Path(reference_config["inpmat_path"]),
        forcing=explicit_forcing,
        start=start,
        days=days,
        calendar=calendar,
        variable_name=variable_name,
    )
    composed_namelist = _namelist(
        map_dir=map_dir,
        diminfo=artifact_dir / "diminfo.txt",
        inpmat=artifact_dir / "inpmat.bin",
        forcing=forcing,
        start=start,
        days=days,
        calendar=calendar,
        variable_name=variable_name,
    )
    _run_case(executable, explicit_case, explicit_namelist)
    _run_case(executable, composed_case, composed_namelist)

    comparisons: dict[str, Any] = {}
    outputs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for variable in ("runoff", "rivsto", "storge", "rivout"):
        explicit_output = _read_output(
            explicit_case / f"{variable}_integration.bin",
            days,
            base_diminfo.ny,
            base_diminfo.nx,
        )
        composed_output = _read_output(
            composed_case / f"{variable}_integration.bin",
            days,
            base_diminfo.ny,
            base_diminfo.nx,
        )
        outputs[variable] = (explicit_output, composed_output)
        comparisons[variable] = _compare_outputs(explicit_output, composed_output)

    if days >= 2:
        split_first = integration_dir / "restart_first_day"
        split_second = integration_dir / "restart_second_day"
        first_namelist = _namelist(
            map_dir=map_dir,
            diminfo=artifact_dir / "diminfo.txt",
            inpmat=artifact_dir / "inpmat.bin",
            forcing=forcing,
            start=start,
            days=1,
            calendar=calendar,
            variable_name=variable_name,
            forcing_start=start,
        )
        _run_case(executable, split_first, first_namelist)
        restart_date = date(*start) + timedelta(days=1)
        restart_path = split_first / (
            f"restart{restart_date.year:04d}{restart_date.month:02d}"
            f"{restart_date.day:02d}00.bin"
        )
        if not restart_path.is_file():
            raise FileNotFoundError(f"restart file was not produced: {restart_path}")
        second_start = (restart_date.year, restart_date.month, restart_date.day)
        second_namelist = _namelist(
            map_dir=map_dir,
            diminfo=artifact_dir / "diminfo.txt",
            inpmat=artifact_dir / "inpmat.bin",
            forcing=forcing,
            start=second_start,
            days=1,
            calendar=calendar,
            variable_name=variable_name,
            forcing_start=start,
            restart_path=restart_path,
        )
        _run_case(executable, split_second, second_namelist)
        for variable in ("runoff", "rivsto", "storge", "rivout"):
            restarted = _read_output(
                split_second / f"{variable}_integration.bin",
                1,
                base_diminfo.ny,
                base_diminfo.nx,
            )
            comparisons[f"restart_{variable}"] = _compare_outputs(
                outputs[variable][1][1:2], restarted
            )

    thresholds: dict[str, dict[str, float]] = {
        "runoff": {"relative_linf": 1.0e-5, "relative_l1": 1.0e-6},
        "rivsto": {"relative_linf": 1.0e-3, "relative_l1": 1.0e-5},
        "storge": {"relative_linf": 1.0e-3, "relative_l1": 1.0e-5},
        "rivout": {"relative_linf": 1.0e-3, "relative_l1": 1.0e-5},
    }
    if days >= 2:
        thresholds.update(
            {
                "restart_runoff": {
                    "relative_linf": 1.0e-7,
                    "relative_l1": 1.0e-7,
                },
                "restart_rivsto": {
                    "relative_linf": 1.0e-7,
                    "relative_l1": 1.0e-7,
                },
                "restart_storge": {
                    "relative_linf": 1.0e-7,
                    "relative_l1": 1.0e-7,
                },
                "restart_rivout": {
                    "relative_linf": 1.0e-7,
                    "relative_l1": 1.0e-7,
                },
            }
        )
    failures: dict[str, dict[str, float]] = {}
    for name, metric_thresholds in thresholds.items():
        failed_metrics = {
            metric: float(comparisons[name][metric])
            for metric, threshold in metric_thresholds.items()
            if float(comparisons[name][metric]) > threshold
        }
        if failed_metrics:
            failures[name] = failed_metrics
    if failures:
        raise ValueError(f"CaMa integration comparison failed: {failures}")

    report = {
        "integration_version": 6,
        "status": "passed",
        "output": {"dirname": dirname},
        "grid_hash": grid_hash,
        "days": days,
        "forcing": str(forcing),
        "executable": str(executable),
        "artifacts": {
            "forcing_sha256": _sha256(forcing),
            "executable_sha256": _sha256(executable),
            "reference_diminfo_sha256": _sha256(
                Path(reference_config["diminfo_path"])
            ),
            "reference_inpmat_sha256": _sha256(
                Path(reference_config["inpmat_path"])
            ),
            "composed_diminfo_sha256": _sha256(artifact_dir / "diminfo.txt"),
            "composed_inpmat_sha256": _sha256(artifact_dir / "inpmat.bin"),
            "tool_source_sha256": tool_source_sha256(),
        },
        "thresholds": thresholds,
        "comparisons": comparisons,
    }
    return report


def run_cama_integration(
    config: dict[str, Any], days: int = 2
) -> dict[str, Any]:
    """Run the end-to-end test in temporary storage and retain only its report."""

    with tempfile.TemporaryDirectory(prefix="cama-integration-") as directory:
        report = _run_cama_integration_in_directory(
            config, Path(directory), days
        )
    validation_reference_dir = configured_validation_reference_dir(config)
    validation_reference_dir.mkdir(parents=True, exist_ok=True)
    report_path = validation_reference_dir / "integration_validation.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report
