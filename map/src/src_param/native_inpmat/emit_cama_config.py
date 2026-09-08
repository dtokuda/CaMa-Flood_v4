#!/usr/bin/env python3
"""Emit a CaMa-Flood forcing environment file for the configured runoff."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from cama_native_inpmat.config import load_config
from cama_native_inpmat.emit import emit_input_config
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
    report = emit_input_config(config)
    write_validation_log(config, "emit_cama_config.py", report)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
