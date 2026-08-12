"""Generate source-grid-specific composed input matrices."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
from typing import Any

from scipy import sparse

from .cama import (
    read_binary_inpmat,
    read_diminfo,
    sha256_file,
    write_binary_inpmat,
    write_diminfo,
)
from .config import (
    configured_output_dir,
    configured_validation_reference_dir,
    resolve_input_files,
)
from .grid import RectilinearGrid, read_rectilinear_grid
from .mapping import (
    area_mean_composed_mapping,
    compose_mapping,
    conservative_remap_matrix,
    csr_to_cama_inpmat,
    normalize_mapping_rows,
    reference_grid_from_diminfo,
    validate_area_mean_mapping,
    validate_mapping,
)


MAPPING_THRESHOLDS = {
    "r_coverage_max_abs_error": 1.0e-12,
    "algebra_relative_linf": 1.0e-12,
    "constant_relative_linf": 1.0e-12,
    "global_mass_relative_error": 1.0e-12,
    "float32_emission_relative_linf": 1.0e-6,
    "preprocessed_float32_relative_linf": 1.0e-6,
}

AREA_MEAN_THRESHOLDS = {
    "reference_row_sum_max_abs_error": 1.0e-12,
    "composed_row_sum_max_abs_error": 1.0e-12,
    "constant_max_abs_error": 1.0e-12,
    "algebra_relative_linf": 1.0e-12,
    "convex_range_relative_violation": 1.0e-12,
    "float32_emission_relative_linf": 1.0e-6,
    "float32_constant_max_abs_error": 1.0e-6,
}


def _write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _temporary_path(path: Path) -> Path:
    """Return a deterministic sibling path that is never a production artifact."""

    if path.suffix == ".npz":
        return path.with_name(f".{path.stem}.tmp.npz")
    return path.with_name(f".{path.name}.tmp")


def _grid_provenance(grid: RectilinearGrid, hash_decimals: int) -> dict[str, Any]:
    result = grid.summary(hash_decimals)
    result.update(
        latitude=grid.latitude.tolist(),
        longitude=grid.longitude.tolist(),
        latitude_bounds=grid.latitude_bounds.tolist(),
        longitude_bounds=grid.longitude_bounds.tolist(),
    )
    return result


def _read_consistent_input_grid(
    input_config: dict[str, Any], hash_decimals: int
) -> tuple[RectilinearGrid, list[Path]]:
    """Read every matched file and require one indexed horizontal grid."""

    files = resolve_input_files(input_config)
    first_grid = read_rectilinear_grid(files[0], input_config["variable"])
    expected_hash = first_grid.hash_at_precision(hash_decimals)
    for path in files[1:]:
        candidate = read_rectilinear_grid(path, input_config["variable"])
        if candidate.hash_at_precision(hash_decimals) != expected_hash:
            raise ValueError(
                f"input.path matched inconsistent horizontal grids: "
                f"{files[0]} and {path}"
            )
    return first_grid, files


def _remove_obsolete_outputs(
    output_dir: Path, validation_reference_dir: Path
) -> list[str]:
    """Remove files from the former verbose layout after replacement succeeds."""

    obsolete = [
        output_dir / name
        for name in (
            "mapping_manifest.json",
            "input_inventory.json",
            "runoff_inventory.json",
            "forcing_validation.json",
            "input_config.json",
            "runoff_config.json",
            "production_validation.json",
            "native_to_reference_flux.npz",
            "native_to_cama_flux.npz",
            "grid.json",
            "validation.json",
            "cama-force.env",
        )
    ]
    obsolete.extend(
        [
            output_dir / "integration",
            validation_reference_dir / "native_to_reference_area_mean.npz",
            validation_reference_dir / "integration_validation.json",
        ]
    )
    removed: list[str] = []
    for path in obsolete:
        if path.is_symlink() or path.is_file():
            path.unlink()
            removed.append(str(path))
        elif path.is_dir():
            shutil.rmtree(path)
            removed.append(str(path))
    return removed


def generate_grid_artifacts(
    *,
    source_grid: RectilinearGrid,
    reference_to_cama: sparse.csr_matrix,
    reference_grid: RectilinearGrid,
    base_diminfo,
    output_dir: Path,
    validation_reference_dir: Path,
    aliases: list[dict[str, str]],
    reference_file: Path,
    hash_decimals: int,
) -> dict[str, Any]:
    native_to_reference = conservative_remap_matrix(source_grid, reference_grid)
    composed = compose_mapping(reference_to_cama, native_to_reference)
    emitted = csr_to_cama_inpmat(composed, base_diminfo, source_grid)
    validation = validate_mapping(
        reference_to_cama,
        native_to_reference,
        composed,
        emitted,
        source_grid,
        reference_grid,
    )

    thresholds = MAPPING_THRESHOLDS
    failures = {
        name: validation[name]
        for name, threshold in thresholds.items()
        if float(validation[name]) > threshold
    }
    if failures:
        raise ValueError(f"mapping validation failed: {failures}")

    reference_area_mean, composed_area_mean = area_mean_composed_mapping(
        reference_to_cama, native_to_reference
    )
    emitted_area_mean = csr_to_cama_inpmat(
        composed_area_mean,
        base_diminfo,
        source_grid,
        inpmat_name="inpmat_area_mean.bin",
    )
    area_mean_validation = validate_area_mean_mapping(
        reference_area_mean,
        native_to_reference,
        composed_area_mean,
        emitted_area_mean,
    )
    area_mean_failures = {
        name: area_mean_validation[name]
        for name, threshold in AREA_MEAN_THRESHOLDS.items()
        if float(area_mean_validation[name]) > threshold
    }
    if area_mean_failures:
        raise ValueError(f"area-mean mapping validation failed: {area_mean_failures}")

    output_dir.mkdir(parents=True, exist_ok=True)
    validation_reference_dir.mkdir(parents=True, exist_ok=True)
    r_path = validation_reference_dir / "native_to_reference_flux.npz"
    c_path = validation_reference_dir / "native_to_cama_flux.npz"
    binary_path = output_dir / "inpmat.bin"
    diminfo_path = output_dir / "diminfo.txt"
    grid_path = validation_reference_dir / "grid.json"
    validation_path = validation_reference_dir / "validation.json"
    c_area_mean_path = validation_reference_dir / "native_to_cama_area_mean.npz"
    binary_area_mean_path = output_dir / "inpmat_area_mean.bin"
    diminfo_area_mean_path = output_dir / "diminfo_area_mean.txt"
    validation_area_mean_path = validation_reference_dir / "validation_area_mean.json"
    final_paths = (
        r_path,
        c_path,
        binary_path,
        diminfo_path,
        c_area_mean_path,
        binary_area_mean_path,
        diminfo_area_mean_path,
        grid_path,
        validation_path,
        validation_area_mean_path,
    )
    temporary_paths = {path: _temporary_path(path) for path in final_paths}
    for path in temporary_paths.values():
        path.unlink(missing_ok=True)

    sparse.save_npz(temporary_paths[r_path], native_to_reference, compressed=True)
    sparse.save_npz(temporary_paths[c_path], composed, compressed=True)
    write_binary_inpmat(temporary_paths[binary_path], emitted)
    write_diminfo(temporary_paths[diminfo_path], emitted.diminfo)
    sparse.save_npz(
        temporary_paths[c_area_mean_path], composed_area_mean, compressed=True
    )
    write_binary_inpmat(
        temporary_paths[binary_area_mean_path], emitted_area_mean
    )
    write_diminfo(
        temporary_paths[diminfo_area_mean_path], emitted_area_mean.diminfo
    )

    grid_record = _grid_provenance(source_grid, hash_decimals)
    grid_record.update(
        reference_file=str(reference_file),
        aliases=aliases,
        reference_grid=reference_grid.summary(hash_decimals),
    )
    _write_json(temporary_paths[grid_path], grid_record)

    validation_record = {
        "validation_version": 2,
        "status": "passed",
        "thresholds": thresholds,
        "metrics": validation,
        "artifacts": {
            "inpmat_sha256": sha256_file(temporary_paths[binary_path]),
            "diminfo_sha256": sha256_file(temporary_paths[diminfo_path]),
            "r_sha256": sha256_file(temporary_paths[r_path]),
            "c_sha256": sha256_file(temporary_paths[c_path]),
            "grid_sha256": sha256_file(temporary_paths[grid_path]),
        },
    }
    _write_json(temporary_paths[validation_path], validation_record)
    area_mean_validation_record = {
        "validation_version": 1,
        "status": "passed",
        "assumptions": {
            "operation": "area_weighted_mean",
            "input_coverage": "complete_global_grid",
            "missing_values_supported": False,
            "weights_normalized_over": "configured_reference_inpmat_covered_area",
        },
        "thresholds": AREA_MEAN_THRESHOLDS,
        "metrics": area_mean_validation,
        "artifacts": {
            "inpmat_sha256": sha256_file(
                temporary_paths[binary_area_mean_path]
            ),
            "diminfo_sha256": sha256_file(
                temporary_paths[diminfo_area_mean_path]
            ),
            "c_sha256": sha256_file(temporary_paths[c_area_mean_path]),
            "grid_sha256": sha256_file(temporary_paths[grid_path]),
        },
    }
    _write_json(
        temporary_paths[validation_area_mean_path], area_mean_validation_record
    )
    for path in final_paths:
        os.replace(temporary_paths[path], path)
    return {
        "grid_hash": source_grid.hash_at_precision(hash_decimals),
        "output_dir": str(output_dir),
        "validation_reference_dir": str(validation_reference_dir),
        "aliases": aliases,
        "validation": validation_record,
        "area_mean_validation": area_mean_validation_record,
    }


def generate_configured_grids(config: dict[str, Any]) -> dict[str, Any]:
    reference_config = config["reference_mapping"]
    diminfo_path = Path(reference_config["diminfo_path"])
    inpmat_path = Path(reference_config["inpmat_path"])
    base_diminfo = read_diminfo(diminfo_path)
    reference_to_cama = read_binary_inpmat(inpmat_path, base_diminfo).to_csr()
    reference_grid = reference_grid_from_diminfo(
        base_diminfo, reference_config["input_grid"]
    )
    reference_cell_count = reference_grid.latitude.size * reference_grid.longitude.size
    if reference_to_cama.shape[1] != reference_cell_count:
        raise ValueError("reference inpmat and configured input grid are inconsistent")

    hash_decimals = config["grid"]["hash_decimals"]
    input_config = config["input"]
    grid, input_files = _read_consistent_input_grid(input_config, hash_decimals)
    reference_file = input_files[0]
    output_dir = configured_output_dir(config)
    validation_reference_dir = configured_validation_reference_dir(config)
    generated = generate_grid_artifacts(
        source_grid=grid,
        reference_to_cama=reference_to_cama,
        reference_grid=reference_grid,
        base_diminfo=base_diminfo,
        output_dir=output_dir,
        validation_reference_dir=validation_reference_dir,
        aliases=[
            {
                "output": {"dirname": config["output"]["dirname"]},
                "variable": input_config["variable"],
            }
        ],
        reference_file=reference_file,
        hash_decimals=hash_decimals,
    )
    removed_obsolete_outputs = _remove_obsolete_outputs(
        output_dir, validation_reference_dir
    )

    return {
        "generation_version": 7,
        "status": "passed",
        "output": {"dirname": config["output"]["dirname"]},
        "reference_diminfo": str(diminfo_path),
        "reference_inpmat": str(inpmat_path),
        "reference_inpmat_sha256": sha256_file(inpmat_path),
        "reference_grid": reference_grid.summary(hash_decimals),
        "grid": generated,
        "input_file_count": len(input_files),
        "removed_obsolete_outputs": removed_obsolete_outputs,
    }


def validate_configured_artifacts(config: dict[str, Any]) -> dict[str, Any]:
    """Reload all generated artifacts and repeat numerical/checksum validation."""

    reference_config = config["reference_mapping"]
    base_diminfo = read_diminfo(reference_config["diminfo_path"])
    reference_to_cama = read_binary_inpmat(
        reference_config["inpmat_path"], base_diminfo
    ).to_csr()
    reference_grid = reference_grid_from_diminfo(
        base_diminfo, reference_config["input_grid"]
    )
    hash_decimals = config["grid"]["hash_decimals"]
    artifact_dir = configured_output_dir(config)
    validation_reference_dir = configured_validation_reference_dir(config)

    input_config = config["input"]
    grid, input_files = _read_consistent_input_grid(input_config, hash_decimals)
    grid_hash = grid.hash_at_precision(hash_decimals)
    grid_record = json.loads((validation_reference_dir / "grid.json").read_text())
    stored = json.loads((validation_reference_dir / "validation.json").read_text())
    stored_area_mean = json.loads(
        (validation_reference_dir / "validation_area_mean.json").read_text()
    )
    if stored.get("validation_version") != 2 or stored.get("status") != "passed":
        raise ValueError(f"{artifact_dir}: validation.json is obsolete or incomplete")
    if stored.get("thresholds") != MAPPING_THRESHOLDS:
        raise ValueError(f"{artifact_dir}: mapping validation thresholds changed")
    if (
        stored_area_mean.get("validation_version") != 1
        or stored_area_mean.get("status") != "passed"
    ):
        raise ValueError(
            f"{artifact_dir}: validation_area_mean.json is obsolete or incomplete"
        )
    if stored_area_mean.get("thresholds") != AREA_MEAN_THRESHOLDS:
        raise ValueError(f"{artifact_dir}: area-mean validation thresholds changed")
    if grid_record["grid_hash"] != grid_hash:
        raise ValueError(f"{artifact_dir}: grid.json hash mismatch")

    native_to_reference = sparse.load_npz(
        validation_reference_dir / "native_to_reference_flux.npz"
    ).tocsr()
    composed = sparse.load_npz(
        validation_reference_dir / "native_to_cama_flux.npz"
    ).tocsr()
    diminfo = read_diminfo(artifact_dir / "diminfo.txt")
    emitted = read_binary_inpmat(artifact_dir / "inpmat.bin", diminfo)
    metrics = validate_mapping(
        reference_to_cama,
        native_to_reference,
        composed,
        emitted,
        grid,
        reference_grid,
    )
    failures = {
        name: metrics[name]
        for name, threshold in MAPPING_THRESHOLDS.items()
        if float(metrics[name]) > float(threshold)
    }
    if failures:
        raise ValueError(f"{artifact_dir}: validation failed: {failures}")

    actual_checksums = {
        "inpmat_sha256": sha256_file(artifact_dir / "inpmat.bin"),
        "diminfo_sha256": sha256_file(artifact_dir / "diminfo.txt"),
        "r_sha256": sha256_file(
            validation_reference_dir / "native_to_reference_flux.npz"
        ),
        "c_sha256": sha256_file(
            validation_reference_dir / "native_to_cama_flux.npz"
        ),
        "grid_sha256": sha256_file(validation_reference_dir / "grid.json"),
    }
    if actual_checksums != stored["artifacts"]:
        raise ValueError(f"{artifact_dir}: artifact checksum mismatch")
    composed_area_mean = sparse.load_npz(
        validation_reference_dir / "native_to_cama_area_mean.npz"
    ).tocsr()
    diminfo_area_mean = read_diminfo(artifact_dir / "diminfo_area_mean.txt")
    emitted_area_mean = read_binary_inpmat(
        artifact_dir / "inpmat_area_mean.bin", diminfo_area_mean
    )
    reference_area_mean = normalize_mapping_rows(reference_to_cama)
    area_mean_metrics = validate_area_mean_mapping(
        reference_area_mean,
        native_to_reference,
        composed_area_mean,
        emitted_area_mean,
    )
    area_mean_failures = {
        name: area_mean_metrics[name]
        for name, threshold in AREA_MEAN_THRESHOLDS.items()
        if float(area_mean_metrics[name]) > float(threshold)
    }
    if area_mean_failures:
        raise ValueError(
            f"{artifact_dir}: area-mean validation failed: {area_mean_failures}"
        )
    actual_area_mean_checksums = {
        "inpmat_sha256": sha256_file(artifact_dir / "inpmat_area_mean.bin"),
        "diminfo_sha256": sha256_file(artifact_dir / "diminfo_area_mean.txt"),
        "c_sha256": sha256_file(
            validation_reference_dir / "native_to_cama_area_mean.npz"
        ),
        "grid_sha256": sha256_file(validation_reference_dir / "grid.json"),
    }
    if actual_area_mean_checksums != stored_area_mean["artifacts"]:
        raise ValueError(f"{artifact_dir}: area-mean artifact checksum mismatch")
    report = {
        "grid_hash": grid_hash,
        "status": "passed",
        "metrics": metrics,
        "artifacts": actual_checksums,
    }
    area_mean_report = {
        "grid_hash": grid_hash,
        "status": "passed",
        "metrics": area_mean_metrics,
        "artifacts": actual_area_mean_checksums,
    }

    return {
        "validation_version": 4,
        "grid_count": 1,
        "input_file_count": len(input_files),
        "status": "passed",
        "grid": report,
        "area_mean": area_mean_report,
    }
