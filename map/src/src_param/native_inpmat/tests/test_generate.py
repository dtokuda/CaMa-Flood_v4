from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np
import xarray as xr

from cama_native_inpmat.generate import _read_consistent_input_grid


class GenerationInputTests(unittest.TestCase):
    @staticmethod
    def _write_grid(path: Path, longitude: np.ndarray) -> None:
        dataset = xr.Dataset(
            data_vars={
                "field": (
                    ("lat", "lon"),
                    np.zeros((2, longitude.size), dtype=np.float32),
                )
            },
            coords={
                "lat": (
                    "lat",
                    [-45.0, 45.0],
                    {"units": "degrees_north"},
                ),
                "lon": (
                    "lon",
                    longitude,
                    {"units": "degrees_east"},
                ),
            },
        )
        dataset.to_netcdf(path)

    def test_generation_rejects_mixed_grids_in_one_input_pattern(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_grid(
                root / "a.nc", np.array([45.0, 135.0, 225.0, 315.0])
            )
            self._write_grid(
                root / "b.nc", np.array([44.0, 134.0, 224.0, 314.0])
            )
            with self.assertRaisesRegex(ValueError, "inconsistent horizontal grids"):
                _read_consistent_input_grid(
                    {"path": str(root / "*.nc"), "variable": "field"}, 7
                )


if __name__ == "__main__":
    unittest.main()
