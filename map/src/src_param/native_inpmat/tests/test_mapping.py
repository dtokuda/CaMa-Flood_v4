from __future__ import annotations

import unittest

import numpy as np
from scipy import sparse

from cama_native_inpmat.cama import CamaDimInfo
from cama_native_inpmat.grid import build_rectilinear_grid
from cama_native_inpmat.mapping import (
    area_mean_composed_mapping,
    compose_mapping,
    conservative_remap_matrix,
    csr_to_cama_inpmat,
    reference_grid_from_diminfo,
    grid_cell_area_factors,
    normalize_mapping_rows,
    standard_one_degree_grid,
    validate_area_mean_mapping,
)


class ConservativeMappingTests(unittest.TestCase):
    def test_half_degree_reference_grid_is_derived_from_diminfo(self) -> None:
        diminfo = CamaDimInfo(
            nx=4,
            ny=2,
            floodplain_layers=1,
            nxin=720,
            nyin=360,
            inpn=1,
            inpmat_name="reference-half-degree.bin",
            west=-180.0,
            east=180.0,
            north=90.0,
            south=-90.0,
        )
        grid = reference_grid_from_diminfo(
            diminfo,
            {
                "longitude_order": "west_to_east",
                "latitude_order": "north_to_south",
            },
        )
        self.assertEqual(grid.shape, (360, 720))
        self.assertAlmostEqual(grid.longitude[0], -179.75)
        self.assertAlmostEqual(grid.longitude[-1], 179.75)
        self.assertAlmostEqual(grid.latitude[0], 89.75)
        self.assertAlmostEqual(grid.latitude[-1], -89.75)

    def test_shifted_reference_grid_is_rejected(self) -> None:
        diminfo = CamaDimInfo(
            nx=4,
            ny=2,
            floodplain_layers=1,
            nxin=360,
            nyin=180,
            inpn=1,
            inpmat_name="shifted.bin",
            west=0.0,
            east=360.0,
            north=90.0,
            south=-90.0,
        )
        with self.assertRaisesRegex(ValueError, "unshifted global grid"):
            reference_grid_from_diminfo(
                diminfo,
                {
                    "longitude_order": "west_to_east",
                    "latitude_order": "north_to_south",
                },
            )

    def assert_coverage_and_mass(self, source) -> None:
        target = standard_one_degree_grid()
        matrix = conservative_remap_matrix(source, target)
        row_sums = np.asarray(matrix.sum(axis=1)).ravel()
        np.testing.assert_allclose(row_sums, 1.0, rtol=0.0, atol=1.0e-12)
        rng = np.random.default_rng(42)
        field = rng.random(source.latitude.size * source.longitude.size)
        source_volume = np.dot(grid_cell_area_factors(source), field)
        target_volume = np.dot(grid_cell_area_factors(target), matrix @ field)
        self.assertLess(abs(target_volume - source_volume) / source_volume, 1.0e-12)

    def test_standard_grid_is_identity(self) -> None:
        grid = standard_one_degree_grid()
        matrix = conservative_remap_matrix(grid, grid)
        self.assertEqual(matrix.shape, (360 * 180, 360 * 180))
        self.assertEqual(matrix.nnz, 360 * 180)
        np.testing.assert_allclose(matrix.diagonal(), 1.0, rtol=0.0, atol=1.0e-14)

    def test_two_degree_grid_conserves_area_integral(self) -> None:
        source = build_rectilinear_grid(
            np.arange(-89.0, 90.0, 2.0), np.arange(1.0, 360.0, 2.0)
        )
        self.assert_coverage_and_mass(source)

    def test_half_degree_grid_conserves_area_integral(self) -> None:
        source = build_rectilinear_grid(
            np.arange(-89.75, 90.0, 0.5), np.arange(0.25, 360.0, 0.5)
        )
        self.assert_coverage_and_mass(source)

    def test_offset_longitude_grid_conserves_area_integral(self) -> None:
        source = build_rectilinear_grid(
            np.arange(-89.0, 90.0, 2.0), np.arange(0.0, 360.0, 2.0)
        )
        self.assert_coverage_and_mass(source)

    def test_gaussian_grid_conserves_area_integral(self) -> None:
        nodes, _ = np.polynomial.legendre.leggauss(32)
        source = build_rectilinear_grid(
            np.degrees(np.arcsin(nodes)), np.arange(2.8125, 360.0, 5.625)
        )
        self.assert_coverage_and_mass(source)

    def test_dateline_indexing_maps_impulses_to_expected_longitudes(self) -> None:
        target = standard_one_degree_grid()
        source = build_rectilinear_grid(
            target.latitude,
            np.arange(0.5, 360.0, 1.0),
            latitude_bounds=target.latitude_bounds,
        )
        matrix = conservative_remap_matrix(source, target)
        equator_north_row = 89 * 360
        source_latitude_index = 89
        source_zero_index = source_latitude_index * 360
        source_180_index = source_latitude_index * 360 + 180
        self.assertEqual(matrix[equator_north_row + 180, source_zero_index], 1.0)
        self.assertEqual(matrix[equator_north_row, source_180_index], 1.0)

    def test_composition_and_float32_encoding(self) -> None:
        source = build_rectilinear_grid(
            np.arange(-89.0, 90.0, 2.0), np.arange(1.0, 360.0, 2.0)
        )
        remap = conservative_remap_matrix(source)
        rows = np.array([0, 0, 1, 2])
        columns = np.array([0, 1, 360, 360 * 180 - 1])
        areas = np.array([1.0, 2.0, 3.0, 4.0])
        cama = sparse.coo_matrix(
            (areas, (rows, columns)), shape=(4, 360 * 180)
        ).tocsr()
        composed = compose_mapping(cama, remap)
        diminfo = CamaDimInfo(
            nx=2,
            ny=2,
            floodplain_layers=1,
            nxin=360,
            nyin=180,
            inpn=1,
            inpmat_name="reference.bin",
            west=-180.0,
            east=180.0,
            north=90.0,
            south=-90.0,
        )
        encoded = csr_to_cama_inpmat(composed, diminfo, source)
        rng = np.random.default_rng(7)
        field = rng.random(source.latitude.size * source.longitude.size)
        np.testing.assert_allclose(
            encoded.to_csr() @ field,
            cama @ (remap @ field),
            rtol=1.0e-7,
            atol=1.0e-7,
        )

    def test_area_mean_mapping_is_normalized_and_convex(self) -> None:
        source = build_rectilinear_grid(
            np.arange(-89.0, 90.0, 2.0), np.arange(1.0, 360.0, 2.0)
        )
        remap = conservative_remap_matrix(source)
        rows = np.array([0, 0, 1, 2])
        columns = np.array([0, 1, 360, 360 * 180 - 1])
        areas = np.array([1.0, 3.0, 2.0, 4.0])
        cama = sparse.coo_matrix(
            (areas, (rows, columns)), shape=(4, 360 * 180)
        ).tocsr()
        reference_mean, composed = area_mean_composed_mapping(cama, remap)
        np.testing.assert_allclose(
            np.asarray(reference_mean.sum(axis=1)).ravel(),
            [1.0, 1.0, 1.0, 0.0],
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            np.asarray(composed.sum(axis=1)).ravel(),
            [1.0, 1.0, 1.0, 0.0],
            rtol=0.0,
            atol=1.0e-15,
        )

        diminfo = CamaDimInfo(
            nx=2,
            ny=2,
            floodplain_layers=1,
            nxin=360,
            nyin=180,
            inpn=1,
            inpmat_name="reference.bin",
            west=-180.0,
            east=180.0,
            north=90.0,
            south=-90.0,
        )
        encoded = csr_to_cama_inpmat(
            composed,
            diminfo,
            source,
            inpmat_name="inpmat_area_mean.bin",
        )
        metrics = validate_area_mean_mapping(
            reference_mean, remap, composed, encoded
        )
        self.assertLess(metrics["algebra_relative_linf"], 1.0e-12)
        self.assertLess(metrics["constant_max_abs_error"], 1.0e-12)
        self.assertEqual(metrics["convex_range_relative_violation"], 0.0)
        self.assertLess(metrics["float32_emission_relative_linf"], 1.0e-6)

    def test_area_mean_normalization_preserves_empty_rows(self) -> None:
        matrix = sparse.csr_matrix(
            np.array([[2.0, 6.0, 0.0], [0.0, 0.0, 0.0]])
        )
        normalized = normalize_mapping_rows(matrix)
        np.testing.assert_allclose(
            normalized.toarray(), [[0.25, 0.75, 0.0], [0, 0, 0]]
        )


if __name__ == "__main__":
    unittest.main()
