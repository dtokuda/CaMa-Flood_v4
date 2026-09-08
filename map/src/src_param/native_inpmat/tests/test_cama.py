from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

import numpy as np

from cama_native_inpmat.cama import (
    CamaDimInfo,
    CamaInpmat,
    read_binary_inpmat,
    read_diminfo,
    sha256_file,
    write_binary_inpmat,
    write_diminfo,
)


LOCAL_MAP = Path("/Users/dtokuda/work/data/cmf_v420_pkg/map/glb_15min")


class CamaInpmatTests(unittest.TestCase):
    def small_inpmat(self) -> CamaInpmat:
        diminfo = CamaDimInfo(
            nx=2,
            ny=2,
            floodplain_layers=1,
            nxin=3,
            nyin=2,
            inpn=2,
            inpmat_name="tiny.bin",
            west=-180.0,
            east=180.0,
            north=90.0,
            south=-90.0,
        )
        inpx = np.array(
            [[[1, 2], [3, 1]], [[0, 3], [0, 2]]], dtype=np.int32
        )
        inpy = np.array(
            [[[1, 1], [1, 2]], [[0, 2], [0, 2]]], dtype=np.int32
        )
        inpa = np.array(
            [[[1.0, 2.0], [3.0, 4.0]], [[0.0, 0.5], [0.0, 1.5]]],
            dtype=np.float32,
        )
        return CamaInpmat(diminfo=diminfo, inpx=inpx, inpy=inpy, inpa=inpa)

    def test_binary_round_trip_and_sparse_matrix(self) -> None:
        original = self.small_inpmat()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tiny.bin"
            write_binary_inpmat(path, original)
            recovered = read_binary_inpmat(path, original.diminfo)
        np.testing.assert_array_equal(recovered.inpx, original.inpx)
        np.testing.assert_array_equal(recovered.inpy, original.inpy)
        np.testing.assert_array_equal(recovered.inpa, original.inpa)
        matrix = recovered.to_csr()
        self.assertEqual(matrix.shape, (4, 6))
        self.assertEqual(matrix.nnz, 6)
        np.testing.assert_allclose(matrix @ np.ones(6), [1.0, 2.5, 3.0, 5.5])

    def test_diminfo_round_trip(self) -> None:
        original = self.small_inpmat().diminfo
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "diminfo.txt"
            write_diminfo(path, original)
            recovered = read_diminfo(path)
        self.assertEqual(recovered, original)

    def test_nonzero_index_in_padding_is_rejected(self) -> None:
        original = self.small_inpmat()
        bad_inpx = original.inpx.copy()
        bad_inpx[1, 0, 0] = 1
        invalid = replace(original, inpx=bad_inpx)
        with self.assertRaisesRegex(ValueError, "zero indices"):
            invalid.validate()

    def test_local_standard_one_degree_baseline(self) -> None:
        diminfo_path = LOCAL_MAP / "diminfo_test-1deg.txt"
        inpmat_path = LOCAL_MAP / "inpmat_test-1deg.bin"
        if not diminfo_path.exists() or not inpmat_path.exists():
            self.skipTest("local CaMa map package is unavailable")
        diminfo = read_diminfo(diminfo_path)
        self.assertEqual(
            (diminfo.nx, diminfo.ny, diminfo.nxin, diminfo.nyin, diminfo.inpn),
            (1440, 720, 360, 180, 7),
        )
        inpmat = read_binary_inpmat(inpmat_path, diminfo)
        summary = inpmat.summary()
        self.assertEqual(summary["active_cama_cells"], 252383)
        self.assertEqual(summary["nonzero_mappings"], 463046)
        self.assertEqual(summary["active_inpx_range"], [1, 360])
        self.assertEqual(summary["active_inpy_range"], [7, 146])
        self.assertEqual(
            sha256_file(inpmat_path),
            "202a78a451db0587ccca1b9283ee4703622e424b9eb583355a346a61aa0f1009",
        )
        matrix = inpmat.to_csr()
        self.assertEqual(matrix.shape, (1440 * 720, 360 * 180))
        self.assertEqual(matrix.nnz, 463046)
        with tempfile.TemporaryDirectory() as directory:
            round_trip = Path(directory) / "inpmat_test-1deg.bin"
            write_binary_inpmat(round_trip, inpmat)
            self.assertEqual(sha256_file(round_trip), sha256_file(inpmat_path))


if __name__ == "__main__":
    unittest.main()
