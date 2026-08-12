from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np
import xarray as xr

from cama_native_inpmat.inventory import inspect_input


class InventoryMetadataTests(unittest.TestCase):
    def test_arbitrary_runoff_file_does_not_require_cmip_attributes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "anything.nc"
            dataset = xr.Dataset(
                data_vars={
                    "runoff": (
                        ("time", "lat", "lon"),
                        np.ones((1, 2, 4), dtype=np.float32),
                        {"units": "mm day-1"},
                    )
                },
                coords={
                    "time": (
                        "time",
                        [0.0],
                        {
                            "units": "days since 2001-01-01",
                            "calendar": "standard",
                        },
                    ),
                    "lat": (
                        "lat",
                        [-45.0, 45.0],
                        {"units": "degrees_north"},
                    ),
                    "lon": (
                        "lon",
                        [45.0, 135.0, 225.0, 315.0],
                        {"units": "degrees_east"},
                    ),
                },
            )
            dataset.to_netcdf(path)
            report = inspect_input(
                {
                    "grid": {"hash_decimals": 7},
                    "output": {"dirname": "arbitrary-input"},
                },
                {
                    "path": str(path),
                    "variable": "runoff",
                },
            )
        self.assertEqual(report["path"], str(path))
        self.assertEqual(report["variable"], "runoff")
        self.assertEqual(report["units"], "mm day-1")
        self.assertEqual(report["dimensions"], ["time", "lat", "lon"])
        self.assertEqual(report["file_count"], 1)


if __name__ == "__main__":
    unittest.main()
