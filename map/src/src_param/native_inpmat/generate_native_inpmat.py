#!/usr/bin/env python3
"""Generate CaMa-format composed matrices for the configured source grid."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from cama_native_inpmat.config import load_config
from cama_native_inpmat.generate import generate_configured_grids
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
    report = generate_configured_grids(config)
    write_validation_log(config, "generate_native_inpmat.py", report, reset=True)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
