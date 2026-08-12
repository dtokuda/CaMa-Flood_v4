"""Aggregate all phase checks into a production-readiness manifest."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any

from .emit import emit_input_config
from .forcing import validate_real_forcing
from .generate import validate_configured_artifacts
from .inventory import inspect_config
from .config import configured_output_dir, configured_validation_reference_dir
from .integration import tool_source_sha256


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_revision(tool_root: Path) -> dict[str, Any]:
    repository = tool_root.parents[3]
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository,
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=repository,
        text=True,
        capture_output=True,
        check=True,
    ).stdout
    source_files = sorted(
        path
        for path in tool_root.rglob("*")
        if path.is_file()
        and "__pycache__" not in path.parts
        and path.suffix in {".py", ".yaml", ".md", ".txt"}
    )
    return {
        "git_revision": revision,
        "git_worktree_dirty": bool(status.strip()),
        "tool_files": {
            str(path.relative_to(tool_root)): _sha256(path) for path in source_files
        },
    }


def validate_production(config: dict[str, Any], tool_root: Path) -> dict[str, Any]:
    """Repeat all checks for the single configured input grid."""

    dirname = config["output"]["dirname"]
    input_dir = configured_output_dir(config)
    input_dir.mkdir(parents=True, exist_ok=True)

    inventory = inspect_config(config)
    mappings = validate_configured_artifacts(config)
    forcing = validate_real_forcing(config)
    emitted_config = emit_input_config(config)

    directory_names = {
        dirname,
        inventory["output"]["dirname"],
        forcing["output"]["dirname"],
        emitted_config["output"]["dirname"],
    }
    if len(directory_names) != 1:
        raise ValueError("output directory names differ across production artifacts")

    expected_hash = inventory["input"]["grid_hash"]
    if mappings["grid"]["grid_hash"] != expected_hash:
        raise ValueError("inventory and mapping grid hashes differ")

    integration_path = (
        configured_validation_reference_dir(config) / "integration_validation.json"
    )
    if not integration_path.is_file():
        raise FileNotFoundError(f"missing CaMa smoke report: {integration_path}")
    integration = json.loads(integration_path.read_text())
    if integration.get("integration_version") != 6:
        raise ValueError(f"{dirname}: CaMa integration report is obsolete")
    if integration.get("status") != "passed":
        raise ValueError(f"{dirname}: CaMa smoke test did not pass")
    if integration.get("output", {}).get("dirname") != dirname:
        raise ValueError(f"{dirname}: smoke-test output.dirname mismatch")
    if integration.get("grid_hash") != expected_hash:
        raise ValueError(f"{dirname}: smoke-test grid hash mismatch")
    if int(integration.get("days", 0)) < 2:
        raise ValueError(f"{dirname}: CaMa validation must include a restart run")
    required_comparisons = {
        "runoff",
        "rivsto",
        "storge",
        "rivout",
        "restart_runoff",
        "restart_rivsto",
        "restart_storge",
        "restart_rivout",
    }
    integration_thresholds = {
        "runoff": {"relative_linf": 1.0e-5, "relative_l1": 1.0e-6},
        "rivsto": {"relative_linf": 1.0e-3, "relative_l1": 1.0e-5},
        "storge": {"relative_linf": 1.0e-3, "relative_l1": 1.0e-5},
        "rivout": {"relative_linf": 1.0e-3, "relative_l1": 1.0e-5},
        "restart_runoff": {"relative_linf": 1.0e-7, "relative_l1": 1.0e-7},
        "restart_rivsto": {"relative_linf": 1.0e-7, "relative_l1": 1.0e-7},
        "restart_storge": {"relative_linf": 1.0e-7, "relative_l1": 1.0e-7},
        "restart_rivout": {"relative_linf": 1.0e-7, "relative_l1": 1.0e-7},
    }
    if integration.get("thresholds") != integration_thresholds:
        raise ValueError(f"{dirname}: integration validation thresholds changed")
    missing_comparisons = required_comparisons.difference(
        integration.get("comparisons", {})
    )
    if missing_comparisons:
        raise ValueError(
            f"{dirname}: integration report is missing comparisons: "
            f"{', '.join(sorted(missing_comparisons))}"
        )
    failed_comparisons: dict[str, dict[str, float]] = {}
    for name in required_comparisons:
        failed_metrics = {
            metric: float(integration["comparisons"][name][metric])
            for metric, threshold in integration_thresholds[name].items()
            if float(integration["comparisons"][name][metric]) > threshold
        }
        if failed_metrics:
            failed_comparisons[name] = failed_metrics
    if failed_comparisons:
        raise ValueError(
            f"{dirname}: stored integration comparisons failed: {failed_comparisons}"
        )

    validation_forcing = inventory.get("validation_forcing", inventory["input"])
    if validation_forcing.get("use_input") is True:
        validation_forcing = inventory["input"]
    forcing_path = Path(validation_forcing["first_file"])
    reference = config["reference_mapping"]
    expected_integration_artifacts = {
        "forcing_sha256": _sha256(forcing_path),
        "executable_sha256": _sha256(Path(config["cama"]["executable"])),
        "reference_diminfo_sha256": _sha256(Path(reference["diminfo_path"])),
        "reference_inpmat_sha256": _sha256(Path(reference["inpmat_path"])),
        "composed_diminfo_sha256": _sha256(input_dir / "diminfo.txt"),
        "composed_inpmat_sha256": _sha256(input_dir / "inpmat.bin"),
        "tool_source_sha256": tool_source_sha256(),
    }
    if integration.get("artifacts") != expected_integration_artifacts:
        raise ValueError(f"{dirname}: integration artifact checksum mismatch")

    result = {
        "production_validation_version": 7,
        "status": "passed",
        "output": {"dirname": dirname},
        "grid_hash": expected_hash,
        "file_count": inventory["input"]["file_count"],
        "source": _source_revision(tool_root),
        "artifacts": {
            "diminfo": str(input_dir / "diminfo.txt"),
            "inpmat": str(input_dir / "inpmat.bin"),
            "area_mean_diminfo": str(input_dir / "diminfo_area_mean.txt"),
            "area_mean_inpmat": str(input_dir / "inpmat_area_mean.bin"),
            "cama_environment": str(input_dir / "cama-force.env"),
            "validation_reference": str(
                configured_validation_reference_dir(config)
            ),
        },
        "inventory": inventory,
        "mapping_validation": mappings,
        "forcing_validation": forcing,
        "cama_configuration": emitted_config,
        "integration": {
            "days": integration["days"],
            "report": str(integration_path),
            "max_relative_linf": max(
                comparison["relative_linf"]
                for comparison in integration["comparisons"].values()
            ),
            "max_relative_l1": max(
                comparison["relative_l1"]
                for comparison in integration["comparisons"].values()
            ),
        },
    }
    return result
