#!/usr/bin/env python3
"""Reload and validate generated native-grid CaMa-Flood input matrices."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from cama_native_inpmat.config import load_config
from cama_native_inpmat.generate import validate_configured_artifacts
from cama_native_inpmat.reporting import write_validation_log


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "config" / "cmip6_mrro.yaml",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    report = validate_configured_artifacts(config)
    write_validation_log(config, "validate_native_inpmat.py", report)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
