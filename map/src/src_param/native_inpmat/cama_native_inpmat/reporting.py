"""Human-readable consolidated logging for command results."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

from .config import configured_output_dir


def validation_log_path(config: dict[str, Any]) -> Path:
    """Return the consolidated log path for one configured input grid."""

    return configured_output_dir(config) / "validation.log"


def write_validation_log(
    config: dict[str, Any],
    command: str,
    report: dict[str, Any],
    *,
    reset: bool = False,
) -> Path:
    """Append one timestamped command report, or start a new generation log."""

    path = validation_log_path(config)
    path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    rendered = json.dumps(report, indent=2, sort_keys=True, allow_nan=False)
    entry = f"[{timestamp}] {command}\n{rendered}\n\n"
    with path.open("w" if reset else "a", encoding="utf-8") as stream:
        stream.write(entry)
    return path
