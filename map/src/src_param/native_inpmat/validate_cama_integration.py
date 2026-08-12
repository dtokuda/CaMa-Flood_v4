#!/usr/bin/env python3
"""Validate two mapping paths with a short CaMa-Flood integration run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from cama_native_inpmat.config import load_config
from cama_native_inpmat.integration import run_cama_integration
from cama_native_inpmat.reporting import write_validation_log


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "config" / "cmip6_mrro.yaml",
    )
    parser.add_argument("--days", type=int, default=2)
    args = parser.parse_args()
    if args.days < 1:
        parser.error("--days must be positive")
    config = load_config(args.config)
    report = run_cama_integration(config, days=args.days)
    write_validation_log(config, "validate_cama_integration.py", report)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
