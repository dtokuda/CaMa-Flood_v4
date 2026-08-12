from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np
import xarray as xr

from cama_native_inpmat.grid import (
    build_rectilinear_grid,
    read_rectilinear_grid,
)


class GridGeometryTests(unittest.TestCase):
    def test_regular_one_degree_grid(self) -> None:
        latitude = np.arange(-89.5, 90.0, 1.0)
        longitude = np.arange(0.5, 360.0, 1.0)
        grid = build_rectilinear_grid(latitude, longitude)
        self.assertEqual(grid.shape, (180, 360))
        np.testing.assert_allclose(grid.latitude_bounds[0], [-90.0, -89.0])
        np.testing.assert_allclose(grid.latitude_bounds[-1], [89.0, 90.0])
        self.assertAlmostEqual(
            float(np.sum(grid.longitude_bounds[:, 1] - grid.longitude_bounds[:, 0])),
            360.0,
        )

    def test_latitude_order_is_preserved_in_hash(self) -> None:
        latitude = np.array([-60.0, 0.0, 60.0])
        longitude = np.array([45.0, 135.0, 225.0, 315.0])
        ascending = build_rectilinear_grid(latitude, longitude)
        descending = build_rectilinear_grid(latitude[::-1], longitude)
        self.assertEqual(ascending.latitude_order, "south_to_north")
        self.assertEqual(descending.latitude_order, "north_to_south")
        self.assertNotEqual(ascending.grid_hash, descending.grid_hash)

    def test_periodic_longitude_conventions_cover_same_cells(self) -> None:
        latitude = np.array([-45.0, 45.0])
        zero_to_360 = build_rectilinear_grid(
            latitude, np.array([45.0, 135.0, 225.0, 315.0])
        )
        minus_180_to_180 = build_rectilinear_grid(
            latitude, np.array([-135.0, -45.0, 45.0, 135.0])
        )
        widths_a = np.sort(
            zero_to_360.longitude_bounds[:, 1]
            - zero_to_360.longitude_bounds[:, 0]
        )
        widths_b = np.sort(
            minus_180_to_180.longitude_bounds[:, 1]
            - minus_180_to_180.longitude_bounds[:, 0]
        )
        np.testing.assert_allclose(widths_a, widths_b)
        self.assertNotEqual(
            zero_to_360.grid_hash,
            minus_180_to_180.grid_hash,
            "source index order must remain part of the reusable mapping key",
        )

    def test_gaussian_bounds_are_recovered(self) -> None:
        nodes, _ = np.polynomial.legendre.leggauss(8)
        latitude = np.degrees(np.arcsin(nodes))
        longitude = np.arange(22.5, 360.0, 45.0)
        grid = build_rectilinear_grid(latitude, longitude)
        self.assertEqual(grid.latitude_bounds_source, "derived_gaussian")
        self.assertEqual(grid.latitude_bounds[0, 0], -90.0)
        self.assertEqual(grid.latitude_bounds[-1, 1], 90.0)
        sine_width = np.diff(np.sin(np.radians(grid.latitude_bounds)), axis=1)[:, 0]
        _, expected_weights = np.polynomial.legendre.leggauss(8)
        np.testing.assert_allclose(sine_width, expected_weights, atol=1.0e-14)

    def test_midpoint_bounds_include_poles(self) -> None:
        grid = build_rectilinear_grid(
            np.array([-90.0, -45.0, 0.0, 45.0, 90.0]),
            np.array([45.0, 135.0, 225.0, 315.0]),
        )
        self.assertEqual(grid.latitude_bounds_source, "derived_midpoint")
        np.testing.assert_allclose(grid.latitude_bounds[0], [-90.0, -67.5])
        np.testing.assert_allclose(grid.latitude_bounds[-1], [67.5, 90.0])

    def test_gap_in_cf_latitude_bounds_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "gap or overlap"):
            build_rectilinear_grid(
                np.array([-45.0, 45.0]),
                np.array([90.0, 270.0]),
                latitude_bounds=np.array([[-90.0, -1.0], [0.0, 90.0]]),
            )

    def test_arbitrary_longitude_permutation_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "monotonic traversal"):
            build_rectilinear_grid(
                np.array([-45.0, 45.0]),
                np.array([45.0, 225.0, 135.0, 315.0]),
            )

    def test_regional_grid_without_bounds_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "global grid"):
            build_rectilinear_grid(
                np.array([35.0, 45.0, 55.0]),
                np.array([0.0, 10.0, 20.0, 30.0]),
            )

    def test_descending_longitude_cf_bounds_are_supported(self) -> None:
        grid = build_rectilinear_grid(
            np.array([-45.0, 45.0]),
            np.array([315.0, 225.0, 135.0, 45.0]),
            longitude_bounds=np.array(
                [[360.0, 270.0], [270.0, 180.0], [180.0, 90.0], [90.0, 0.0]]
            ),
        )
        np.testing.assert_allclose(
            grid.longitude_bounds,
            [[270.0, 360.0], [180.0, 270.0], [90.0, 180.0], [0.0, 90.0]],
        )

    def test_tiny_cf_edge_mismatches_are_snapped(self) -> None:
        grid = build_rectilinear_grid(
            np.array([-45.0, 45.0]),
            np.array([45.0, 135.0, 225.0, 315.0]),
            latitude_bounds=np.array(
                [[-90.0, 3.0e-8], [-3.0e-8, 90.0]]
            ),
            longitude_bounds=np.array(
                [
                    [0.0, 90.0 + 3.0e-8],
                    [90.0 - 3.0e-8, 180.0],
                    [180.0, 270.0],
                    [270.0, 360.0],
                ]
            ),
        )
        self.assertIn("snapped_latitude_internal_edges", grid.repairs)
        self.assertIn("snapped_longitude_internal_edges", grid.repairs)
        self.assertEqual(grid.latitude_bounds[0, 1], grid.latitude_bounds[1, 0])
        self.assertEqual(grid.longitude_bounds[0, 1], grid.longitude_bounds[1, 0])

    def test_curvilinear_netcdf_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "curvilinear.nc"
            dataset = xr.Dataset(
                data_vars={
                    "mrro": (("time", "y", "x"), np.zeros((1, 2, 2))),
                    "lat": (("y", "x"), np.array([[-45.0, -45.0], [45.0, 45.0]])),
                    "lon": (("y", "x"), np.array([[0.0, 180.0], [0.0, 180.0]])),
                },
                coords={"time": [0], "y": [0, 1], "x": [0, 1]},
            )
            dataset.to_netcdf(path)
            with self.assertRaisesRegex(ValueError, "one-dimensional"):
                read_rectilinear_grid(path)

    def test_nondegree_coordinate_units_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "radians.nc"
            dataset = xr.Dataset(
                data_vars={"field": (("lat", "lon"), np.zeros((2, 4)))},
                coords={
                    "lat": (
                        "lat",
                        np.radians([-45.0, 45.0]),
                        {"units": "radian", "standard_name": "latitude"},
                    ),
                    "lon": (
                        "lon",
                        np.radians([45.0, 135.0, 225.0, 315.0]),
                        {"units": "radian", "standard_name": "longitude"},
                    ),
                },
            )
            dataset.to_netcdf(path)
            with self.assertRaisesRegex(ValueError, "units must be degrees"):
                read_rectilinear_grid(path, "field")


if __name__ == "__main__":
    unittest.main()
