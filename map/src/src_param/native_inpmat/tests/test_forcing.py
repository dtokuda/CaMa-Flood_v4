from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np
import xarray as xr

from cama_native_inpmat.emit import calendar_to_lleapyr
from cama_native_inpmat.forcing import (
    preprocess_runoff,
    validate_forcing_collection,
)


class ForcingTests(unittest.TestCase):
    def test_preprocess_matches_cama_nan_and_negative_rules(self) -> None:
        values = np.array([[np.nan, -2.0, 0.0, 3.0]], dtype=np.float32)
        result, report = preprocess_runoff(values)
        np.testing.assert_array_equal(result, [[0.0, 0.0, 0.0, 3.0]])
        self.assertEqual(report["missing_count"], 1)
        self.assertEqual(report["negative_count"], 1)

    def test_raw_fill_value_is_converted_before_clipping(self) -> None:
        values = np.array([[np.float32(1.0e20), 2.0]], dtype=np.float32)
        result, report = preprocess_runoff(
            values, missing_values=(float(np.float32(1.0e20)),)
        )
        np.testing.assert_array_equal(result, [[0.0, 2.0]])
        self.assertEqual(report["missing_count"], 1)

    def test_calendar_mapping(self) -> None:
        self.assertEqual(calendar_to_lleapyr("365_day"), ".FALSE.")
        self.assertEqual(calendar_to_lleapyr("noleap"), ".FALSE.")
        self.assertEqual(calendar_to_lleapyr("gregorian"), ".TRUE.")
        self.assertEqual(calendar_to_lleapyr("proleptic_gregorian"), ".TRUE.")
        with self.assertRaisesRegex(ValueError, "unsupported calendar"):
            calendar_to_lleapyr("360_day")

    def test_every_forcing_file_is_checked_for_cama_missing_value(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for index, fill_value in enumerate((1.0e20, 9.0e19)):
                path = Path(directory) / f"forcing-{index}.nc"
                dataset = xr.Dataset(
                    data_vars={
                        "mrro": (
                            ("time", "lat", "lon"),
                            np.ones((1, 2, 2), dtype=np.float32),
                        )
                    },
                    coords={
                        "time": [0.0],
                        "lat": [-45.0, 45.0],
                        "lon": [90.0, 270.0],
                    },
                )
                dataset["time"].attrs.update(
                    units="days since 2000-01-01", calendar="standard"
                )
                dataset.to_netcdf(
                    path,
                    encoding={"mrro": {"_FillValue": np.float32(fill_value)}},
                )
                paths.append(path)
            with self.assertRaisesRegex(ValueError, "_FillValue is not CaMa"):
                validate_forcing_collection(paths, "mrro")


if __name__ == "__main__":
    unittest.main()
