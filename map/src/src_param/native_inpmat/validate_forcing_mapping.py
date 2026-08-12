#!/usr/bin/env python3
"""Validate explicit and composed mappings on real runoff timesteps."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from cama_native_inpmat.config import load_config
from cama_native_inpmat.forcing import validate_real_forcing
from cama_native_inpmat.reporting import write_validation_log


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "config" / "cmip6_mrro.yaml",
    )
    args = parser.parse_args()
    config = load_config(args.config)
    report = validate_real_forcing(config)
    write_validation_log(config, "validate_forcing_mapping.py", report)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
