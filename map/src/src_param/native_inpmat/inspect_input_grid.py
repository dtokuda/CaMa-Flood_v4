#!/usr/bin/env python3
"""Inspect all input files selected by a YAML configuration."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from cama_native_inpmat.config import load_config
from cama_native_inpmat.inventory import inspect_config
from cama_native_inpmat.reporting import write_validation_log


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "config" / "cmip6_mrro.yaml",
        help="YAML file selecting the input grid",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="optional JSON output path; the report is always printed to stdout",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    report = inspect_config(config)
    rendered = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    write_validation_log(config, "inspect_input_grid.py", report)
    print(rendered)


if __name__ == "__main__":
    main()
