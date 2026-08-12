"""Readers and writers for the standard CaMa-Flood binary input matrix."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse


@dataclass(frozen=True)
class CamaDimInfo:
    nx: int
    ny: int
    floodplain_layers: int
    nxin: int
    nyin: int
    inpn: int
    inpmat_name: str
    west: float
    east: float
    north: float
    south: float


@dataclass(frozen=True)
class CamaInpmat:
    """CaMa input matrix in record order: (inpn, ny, nx)."""

    diminfo: CamaDimInfo
    inpx: np.ndarray
    inpy: np.ndarray
    inpa: np.ndarray

    def validate(self) -> None:
        expected = (self.diminfo.inpn, self.diminfo.ny, self.diminfo.nx)
        for name, values in (
            ("inpx", self.inpx),
            ("inpy", self.inpy),
            ("inpa", self.inpa),
        ):
            if values.shape != expected:
                raise ValueError(f"{name} has shape {values.shape}, expected {expected}")

        if not np.all(np.isfinite(self.inpa)):
            raise ValueError("inpa contains non-finite values")
        active = self.inpa > 0.0
        if np.any(self.inpa < 0.0):
            raise ValueError("inpa contains negative values")
        if np.any((self.inpx[active] < 1) | (self.inpx[active] > self.diminfo.nxin)):
            raise ValueError("active inpx is outside 1-based input-grid bounds")
        if np.any((self.inpy[active] < 1) | (self.inpy[active] > self.diminfo.nyin)):
            raise ValueError("active inpy is outside 1-based input-grid bounds")
        padding = ~active
        if np.any(self.inpx[padding] != 0) or np.any(self.inpy[padding] != 0):
            raise ValueError("inactive inpmat slots must have zero indices")

        seen_padding = np.maximum.accumulate(padding, axis=0)
        if np.any(active & seen_padding):
            raise ValueError("active entries occur after zero padding")

    def to_csr(self) -> sparse.csr_matrix:
        self.validate()
        active = self.inpa > 0.0
        ncell = self.diminfo.nx * self.diminfo.ny
        rows = np.broadcast_to(
            np.arange(ncell, dtype=np.int64).reshape(1, self.diminfo.ny, self.diminfo.nx),
            self.inpa.shape,
        )[active]
        columns = (
            (self.inpy[active].astype(np.int64) - 1) * self.diminfo.nxin
            + self.inpx[active].astype(np.int64)
            - 1
        )
        matrix = sparse.coo_matrix(
            (self.inpa[active].astype(np.float64), (rows, columns)),
            shape=(ncell, self.diminfo.nxin * self.diminfo.nyin),
        )
        return matrix.tocsr()

    def summary(self) -> dict[str, Any]:
        self.validate()
        active = self.inpa > 0.0
        active_cells = np.any(active, axis=0)
        result: dict[str, Any] = {
            "shape": [self.diminfo.nx, self.diminfo.ny, self.diminfo.inpn],
            "input_shape": [self.diminfo.nxin, self.diminfo.nyin],
            "active_cama_cells": int(np.count_nonzero(active_cells)),
            "nonzero_mappings": int(np.count_nonzero(active)),
            "inpa_sum_m2": float(np.sum(self.inpa, dtype=np.float64)),
        }
        if np.any(active):
            result.update(
                active_inpx_range=[
                    int(np.min(self.inpx[active])),
                    int(np.max(self.inpx[active])),
                ],
                active_inpy_range=[
                    int(np.min(self.inpy[active])),
                    int(np.max(self.inpy[active])),
                ],
                inpa_range_m2=[
                    float(np.min(self.inpa[active])),
                    float(np.max(self.inpa[active])),
                ],
            )
        return result


def _clean_diminfo_lines(path: Path) -> list[str]:
    lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        value = line.split("!!", 1)[0].strip()
        if value:
            lines.append(value)
    return lines


def read_diminfo(path: str | Path) -> CamaDimInfo:
    diminfo_path = Path(path)
    values = _clean_diminfo_lines(diminfo_path)
    if len(values) < 11:
        raise ValueError(f"{diminfo_path}: expected at least 11 data lines")
    try:
        integers = [int(values[index]) for index in range(6)]
        edges = [float(values[index]) for index in range(7, 11)]
    except ValueError as error:
        raise ValueError(f"{diminfo_path}: invalid numeric value") from error
    result = CamaDimInfo(
        nx=integers[0],
        ny=integers[1],
        floodplain_layers=integers[2],
        nxin=integers[3],
        nyin=integers[4],
        inpn=integers[5],
        inpmat_name=values[6],
        west=edges[0],
        east=edges[1],
        north=edges[2],
        south=edges[3],
    )
    if min(result.nx, result.ny, result.nxin, result.nyin, result.inpn) <= 0:
        raise ValueError(f"{diminfo_path}: dimensions must be positive")
    return result


def write_diminfo(path: str | Path, diminfo: CamaDimInfo) -> None:
    """Write the legacy text dimension file consumed by CaMa-Flood."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"{diminfo.nx:10d}     !! nXX",
        f"{diminfo.ny:10d}     !! nYY",
        f"{diminfo.floodplain_layers:10d}     !! floodplain layer",
        f"{diminfo.nxin:10d}     !! input nXX",
        f"{diminfo.nyin:10d}     !! input nYY",
        f"{diminfo.inpn:10d}     !! input num",
        diminfo.inpmat_name,
        f"{diminfo.west:12.3f}     !! west  edge",
        f"{diminfo.east:12.3f}     !! east  edge",
        f"{diminfo.north:12.3f}     !! north edge",
        f"{diminfo.south:12.3f}     !! south edge",
    ]
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_binary_inpmat(
    path: str | Path, diminfo: CamaDimInfo, *, endian: str = "<"
) -> CamaInpmat:
    """Read three groups of direct-access records written by Fortran."""

    inpmat_path = Path(path)
    cells_per_block = diminfo.nx * diminfo.ny * diminfo.inpn
    expected_bytes = cells_per_block * 4 * 3
    actual_bytes = inpmat_path.stat().st_size
    if actual_bytes != expected_bytes:
        raise ValueError(
            f"{inpmat_path}: size is {actual_bytes}, expected {expected_bytes} bytes"
        )

    raw = inpmat_path.read_bytes()
    int_dtype = np.dtype(f"{endian}i4")
    float_dtype = np.dtype(f"{endian}f4")
    shape = (diminfo.inpn, diminfo.ny, diminfo.nx)
    inpx = np.frombuffer(raw, dtype=int_dtype, count=cells_per_block).reshape(shape).copy()
    inpy = np.frombuffer(
        raw, dtype=int_dtype, count=cells_per_block, offset=cells_per_block * 4
    ).reshape(shape).copy()
    inpa = np.frombuffer(
        raw, dtype=float_dtype, count=cells_per_block, offset=cells_per_block * 8
    ).reshape(shape).copy()
    result = CamaInpmat(diminfo=diminfo, inpx=inpx, inpy=inpy, inpa=inpa)
    result.validate()
    return result


def write_binary_inpmat(
    path: str | Path, inpmat: CamaInpmat, *, endian: str = "<"
) -> None:
    """Write a CaMa binary input matrix in standard direct-record order."""

    inpmat.validate()
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as stream:
        stream.write(np.asarray(inpmat.inpx, dtype=f"{endian}i4").tobytes(order="C"))
        stream.write(np.asarray(inpmat.inpy, dtype=f"{endian}i4").tobytes(order="C"))
        stream.write(np.asarray(inpmat.inpa, dtype=f"{endian}f4").tobytes(order="C"))


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
