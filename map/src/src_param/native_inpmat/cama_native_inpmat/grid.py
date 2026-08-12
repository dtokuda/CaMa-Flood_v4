"""Rectilinear spherical-grid inspection and canonicalization."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr


GEOMETRY_TOLERANCE_DEGREES = 1.0e-7
DERIVED_LONGITUDE_SPACING_RTOL = 1.0e-4

_LATITUDE_DEGREE_UNITS = {
    "degree_north",
    "degrees_north",
    "degree_n",
    "degrees_n",
}
_LONGITUDE_DEGREE_UNITS = {
    "degree_east",
    "degrees_east",
    "degree_e",
    "degrees_e",
}


def _array_1d(values: Any, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size < 2:
        raise ValueError(f"{name} must be a one-dimensional array with at least 2 cells")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values")
    return array


def _bounds_array(values: Any, size: int, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.shape == (2, size):
        array = array.T
    if array.shape != (size, 2):
        raise ValueError(f"{name} must have shape ({size}, 2), got {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values")
    return array


def _strict_monotonic(values: np.ndarray, name: str) -> int:
    delta = np.diff(values)
    if np.all(delta > 0.0):
        return 1
    if np.all(delta < 0.0):
        return -1
    raise ValueError(f"{name} must be strictly monotonic")


def _normalized_longitudes(values: np.ndarray) -> np.ndarray:
    normalized = np.mod(values, 360.0)
    normalized[np.isclose(normalized, 360.0, atol=1.0e-12)] = 0.0
    return normalized


def _validate_periodic_longitude_order(longitude: np.ndarray) -> None:
    """Accept one cyclic eastward or westward traversal, not a permutation."""

    centers = _normalized_longitudes(longitude)
    cyclic_delta = np.diff(np.r_[centers, centers[0]])
    eastward = np.mod(cyclic_delta, 360.0)
    westward = np.mod(-cyclic_delta, 360.0)
    eastward_once = np.all(eastward > GEOMETRY_TOLERANCE_DEGREES) and np.isclose(
        np.sum(eastward), 360.0, rtol=0.0, atol=GEOMETRY_TOLERANCE_DEGREES
    )
    westward_once = np.all(westward > GEOMETRY_TOLERANCE_DEGREES) and np.isclose(
        np.sum(westward), 360.0, rtol=0.0, atol=GEOMETRY_TOLERANCE_DEGREES
    )
    if not (eastward_once or westward_once):
        raise ValueError("longitude must make one monotonic traversal of the globe")


def derive_longitude_bounds(longitude: np.ndarray) -> np.ndarray:
    """Derive periodic midpoint bounds while preserving source-cell indices."""

    centers = _normalized_longitudes(longitude)
    order = np.argsort(centers)
    sorted_centers = centers[order]
    gaps = np.diff(np.r_[sorted_centers, sorted_centers[0] + 360.0])
    if np.any(gaps <= GEOMETRY_TOLERANCE_DEGREES):
        raise ValueError("longitude centers are duplicated or not periodic")
    expected_gap = 360.0 / longitude.size
    if not np.allclose(
        gaps,
        expected_gap,
        rtol=DERIVED_LONGITUDE_SPACING_RTOL,
        atol=GEOMETRY_TOLERANCE_DEGREES,
    ):
        raise ValueError(
            "longitude bounds are absent and centers are not a uniformly spaced "
            "global grid"
        )

    previous = np.roll(gaps, 1)
    west = sorted_centers - 0.5 * previous
    east = sorted_centers + 0.5 * gaps
    sorted_bounds = np.column_stack((west, east))
    bounds = np.empty_like(sorted_bounds)
    bounds[order] = sorted_bounds
    return bounds


def _is_gaussian_latitude(latitude: np.ndarray, atol: float = 1.0e-6) -> bool:
    nodes, _ = np.polynomial.legendre.leggauss(latitude.size)
    expected = np.degrees(np.arcsin(nodes))
    candidate = latitude if latitude[0] < latitude[-1] else latitude[::-1]
    return bool(np.allclose(candidate, expected, rtol=0.0, atol=atol))


def derive_gaussian_latitude_bounds(latitude: np.ndarray) -> np.ndarray:
    """Recover Gaussian-cell bounds from Gauss-Legendre quadrature weights."""

    nodes, weights = np.polynomial.legendre.leggauss(latitude.size)
    expected = np.degrees(np.arcsin(nodes))
    ascending = latitude[0] < latitude[-1]
    candidate = latitude if ascending else latitude[::-1]
    if not np.allclose(candidate, expected, rtol=0.0, atol=1.0e-6):
        raise ValueError("latitude coordinates are not Gaussian nodes")

    sine_edges = np.r_[-1.0, -1.0 + np.cumsum(weights)]
    sine_edges[0] = -1.0
    sine_edges[-1] = 1.0
    edges = np.degrees(np.arcsin(np.clip(sine_edges, -1.0, 1.0)))
    bounds = np.column_stack((edges[:-1], edges[1:]))
    return bounds if ascending else bounds[::-1]


def derive_midpoint_latitude_bounds(latitude: np.ndarray) -> np.ndarray:
    """Derive midpoint bounds and close a global grid at both poles."""

    direction = _strict_monotonic(latitude, "latitude")
    ascending = latitude if direction > 0 else latitude[::-1]
    if ascending[0] < -90.0 - GEOMETRY_TOLERANCE_DEGREES:
        raise ValueError("latitude is south of -90 degrees")
    if ascending[-1] > 90.0 + GEOMETRY_TOLERANCE_DEGREES:
        raise ValueError("latitude is north of 90 degrees")
    south_gap = float(ascending[0] + 90.0)
    north_gap = float(90.0 - ascending[-1])
    south_is_global = south_gap <= GEOMETRY_TOLERANCE_DEGREES or np.isclose(
        south_gap,
        0.5 * float(ascending[1] - ascending[0]),
        rtol=0.25,
        atol=GEOMETRY_TOLERANCE_DEGREES,
    )
    north_is_global = north_gap <= GEOMETRY_TOLERANCE_DEGREES or np.isclose(
        north_gap,
        0.5 * float(ascending[-1] - ascending[-2]),
        rtol=0.25,
        atol=GEOMETRY_TOLERANCE_DEGREES,
    )
    if not south_is_global or not north_is_global:
        raise ValueError(
            "latitude bounds are absent and centers do not describe a global grid"
        )

    edges = np.empty(ascending.size + 1, dtype=np.float64)
    edges[0] = -90.0
    edges[-1] = 90.0
    edges[1:-1] = 0.5 * (ascending[:-1] + ascending[1:])
    bounds = np.column_stack((edges[:-1], edges[1:]))
    return bounds if direction > 0 else bounds[::-1]


def canonicalize_longitude_bounds(
    longitude: np.ndarray, bounds: np.ndarray
) -> np.ndarray:
    """Represent periodic intervals west-to-east, regardless of vertex order."""

    centers = _normalized_longitudes(longitude)
    result = np.empty_like(bounds)
    for index, ((first, second), center) in enumerate(zip(bounds, centers, strict=True)):
        candidates: list[tuple[float, float]] = []
        for west_raw, east_raw in ((first, second), (second, first)):
            raw_width = east_raw - west_raw
            width = float(np.mod(raw_width, 360.0))
            if np.isclose(width, 0.0, atol=GEOMETRY_TOLERANCE_DEGREES):
                if np.isclose(
                    abs(raw_width), 360.0, atol=GEOMETRY_TOLERANCE_DEGREES
                ):
                    width = 360.0
                else:
                    continue
            west = float(np.mod(west_raw, 360.0))
            center_distance = float(np.mod(center - west, 360.0))
            if center_distance <= width + GEOMETRY_TOLERANCE_DEGREES:
                candidates.append((width, west))
        if not candidates:
            raise ValueError(f"longitude center {index} is outside its bounds")
        width, west = min(candidates)
        result[index] = (west, west + width)
    return result


def canonicalize_latitude_bounds(
    latitude: np.ndarray, bounds: np.ndarray
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Order latitude pairs south-to-north and snap near-polar outer edges."""

    result = np.sort(bounds, axis=1)
    repairs: list[str] = []
    south_index = int(np.argmin(latitude))
    north_index = int(np.argmax(latitude))

    if abs(result[south_index, 0] + 90.0) <= GEOMETRY_TOLERANCE_DEGREES:
        if result[south_index, 0] != -90.0:
            repairs.append("snapped_south_pole")
        result[south_index, 0] = -90.0
    if abs(result[north_index, 1] - 90.0) <= GEOMETRY_TOLERANCE_DEGREES:
        if result[north_index, 1] != 90.0:
            repairs.append("snapped_north_pole")
        result[north_index, 1] = 90.0
    if np.any(latitude < result[:, 0] - GEOMETRY_TOLERANCE_DEGREES) or np.any(
        latitude > result[:, 1] + GEOMETRY_TOLERANCE_DEGREES
    ):
        raise ValueError("a latitude center is outside its bounds")

    order = np.argsort(latitude)
    ordered = result[order].copy()
    mismatch = ordered[1:, 0] - ordered[:-1, 1]
    if np.any(np.abs(mismatch) > GEOMETRY_TOLERANCE_DEGREES):
        raise ValueError("latitude bounds contain a gap or overlap")
    if np.any(mismatch != 0.0):
        shared = 0.5 * (ordered[1:, 0] + ordered[:-1, 1])
        ordered[:-1, 1] = shared
        ordered[1:, 0] = shared
        result[order] = ordered
        repairs.append("snapped_latitude_internal_edges")
    return result, tuple(repairs)


def snap_longitude_bounds(
    longitude: np.ndarray, bounds: np.ndarray
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Snap only sub-tolerance gaps/overlaps at periodic shared edges."""

    centers = _normalized_longitudes(longitude)
    order = np.argsort(centers)
    ordered = bounds[order].copy()
    unwrapped = np.empty_like(ordered)
    for position, ((west, east), center) in enumerate(
        zip(ordered, centers[order], strict=True)
    ):
        midpoint = 0.5 * (west + east)
        shift = 360.0 * np.round((center - midpoint) / 360.0)
        unwrapped[position] = (west + shift, east + shift)

    internal_mismatch = unwrapped[1:, 0] - unwrapped[:-1, 1]
    periodic_mismatch = unwrapped[0, 0] + 360.0 - unwrapped[-1, 1]
    mismatch = np.r_[internal_mismatch, periodic_mismatch]
    if np.any(np.abs(mismatch) > GEOMETRY_TOLERANCE_DEGREES):
        raise ValueError("longitude bounds contain a gap or overlap")
    if not np.any(mismatch != 0.0):
        return bounds, ()

    edges = np.empty(longitude.size + 1, dtype=np.float64)
    edges[0] = 0.5 * (unwrapped[0, 0] + unwrapped[-1, 1] - 360.0)
    edges[1:-1] = 0.5 * (unwrapped[:-1, 1] + unwrapped[1:, 0])
    edges[-1] = edges[0] + 360.0
    snapped_ordered = np.column_stack((edges[:-1], edges[1:]))
    snapped = np.empty_like(snapped_ordered)
    snapped[order] = snapped_ordered
    return canonicalize_longitude_bounds(longitude, snapped), (
        "snapped_longitude_internal_edges",
    )


def _validate_latitude_coverage(latitude: np.ndarray, bounds: np.ndarray) -> None:
    order = np.argsort(latitude)
    ordered = bounds[order]
    if abs(ordered[0, 0] + 90.0) > GEOMETRY_TOLERANCE_DEGREES:
        raise ValueError("latitude bounds do not reach the south pole")
    if abs(ordered[-1, 1] - 90.0) > GEOMETRY_TOLERANCE_DEGREES:
        raise ValueError("latitude bounds do not reach the north pole")
    mismatch = ordered[1:, 0] - ordered[:-1, 1]
    if np.any(np.abs(mismatch) > GEOMETRY_TOLERANCE_DEGREES):
        raise ValueError("latitude bounds contain a gap or overlap")


def _validate_longitude_coverage(bounds: np.ndarray) -> None:
    segments: list[tuple[float, float]] = []
    for west, east in bounds:
        if east <= 360.0 + GEOMETRY_TOLERANCE_DEGREES:
            segments.append((west, min(east, 360.0)))
        else:
            segments.append((west, 360.0))
            segments.append((0.0, east - 360.0))
    segments.sort()
    segments = [segment for segment in segments if segment[1] - segment[0] > 1.0e-12]
    if abs(segments[0][0]) > GEOMETRY_TOLERANCE_DEGREES:
        raise ValueError("longitude bounds do not start at the periodic seam")
    if abs(segments[-1][1] - 360.0) > GEOMETRY_TOLERANCE_DEGREES:
        raise ValueError("longitude bounds do not end at the periodic seam")
    for previous, current in zip(segments[:-1], segments[1:], strict=True):
        if abs(current[0] - previous[1]) > GEOMETRY_TOLERANCE_DEGREES:
            raise ValueError("longitude bounds contain a gap or overlap")


@dataclass(frozen=True)
class RectilinearGrid:
    """One indexed, global rectilinear source grid."""

    latitude: np.ndarray
    longitude: np.ndarray
    latitude_bounds: np.ndarray
    longitude_bounds: np.ndarray
    latitude_name: str = "lat"
    longitude_name: str = "lon"
    latitude_bounds_source: str = "derived_midpoint"
    longitude_bounds_source: str = "derived_midpoint"
    repairs: tuple[str, ...] = ()

    @property
    def shape(self) -> tuple[int, int]:
        return self.latitude.size, self.longitude.size

    @property
    def latitude_order(self) -> str:
        return "south_to_north" if self.latitude[0] < self.latitude[-1] else "north_to_south"

    @property
    def grid_hash(self) -> str:
        return self.hash_at_precision(7)

    def hash_at_precision(self, decimals: int) -> str:
        """Hash geometry and source index order after sub-metre quantization."""

        payload = {
            "hash_version": 1,
            "shape": list(self.shape),
            "latitude": np.round(self.latitude, decimals).tolist(),
            "longitude_mod_360": np.round(
                _normalized_longitudes(self.longitude), decimals
            ).tolist(),
            "latitude_bounds": np.round(self.latitude_bounds, decimals).tolist(),
            "longitude_bounds": np.round(self.longitude_bounds, decimals).tolist(),
        }
        encoded = json.dumps(
            payload, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def summary(self, hash_decimals: int = 7) -> dict[str, Any]:
        return {
            "shape": list(self.shape),
            "latitude_order": self.latitude_order,
            "latitude_range": [float(np.min(self.latitude)), float(np.max(self.latitude))],
            "longitude_range": [
                float(np.min(self.longitude)),
                float(np.max(self.longitude)),
            ],
            "latitude_bounds_source": self.latitude_bounds_source,
            "longitude_bounds_source": self.longitude_bounds_source,
            "repairs": list(self.repairs),
            "grid_hash": self.hash_at_precision(hash_decimals),
            "hash_decimals": hash_decimals,
        }


def build_rectilinear_grid(
    latitude: Any,
    longitude: Any,
    *,
    latitude_bounds: Any | None = None,
    longitude_bounds: Any | None = None,
    latitude_name: str = "lat",
    longitude_name: str = "lon",
    latitude_bounds_source: str | None = None,
    longitude_bounds_source: str | None = None,
) -> RectilinearGrid:
    """Build and validate a global rectilinear grid from coordinate arrays."""

    lat = _array_1d(latitude, latitude_name)
    lon = _array_1d(longitude, longitude_name)
    _strict_monotonic(lat, latitude_name)

    normalized_lon = _normalized_longitudes(lon)
    if np.unique(np.round(normalized_lon, 12)).size != lon.size:
        raise ValueError("longitude contains duplicate periodic coordinates")
    _validate_periodic_longitude_order(lon)

    repairs: list[str] = []
    if latitude_bounds is None:
        if _is_gaussian_latitude(lat):
            lat_bounds = derive_gaussian_latitude_bounds(lat)
            lat_source = "derived_gaussian"
        else:
            lat_bounds = derive_midpoint_latitude_bounds(lat)
            lat_source = "derived_midpoint"
    else:
        lat_bounds = _bounds_array(latitude_bounds, lat.size, "latitude bounds")
        lat_source = latitude_bounds_source or "cf"

    lat_bounds, lat_repairs = canonicalize_latitude_bounds(lat, lat_bounds)
    repairs.extend(lat_repairs)
    _validate_latitude_coverage(lat, lat_bounds)

    if longitude_bounds is None:
        lon_bounds = derive_longitude_bounds(lon)
        lon_source = "derived_periodic_midpoint"
    else:
        lon_bounds = _bounds_array(longitude_bounds, lon.size, "longitude bounds")
        lon_source = longitude_bounds_source or "cf"
    lon_bounds = canonicalize_longitude_bounds(lon, lon_bounds)
    lon_bounds, lon_repairs = snap_longitude_bounds(lon, lon_bounds)
    repairs.extend(lon_repairs)
    _validate_longitude_coverage(lon_bounds)

    return RectilinearGrid(
        latitude=lat,
        longitude=lon,
        latitude_bounds=lat_bounds,
        longitude_bounds=lon_bounds,
        latitude_name=latitude_name,
        longitude_name=longitude_name,
        latitude_bounds_source=lat_source,
        longitude_bounds_source=lon_source,
        repairs=tuple(repairs),
    )


def _coordinate_name(dataset: xr.Dataset, kind: str) -> str:
    preferred = "lat" if kind == "latitude" else "lon"
    if preferred in dataset.variables:
        return preferred
    standard_name = kind
    axis = "Y" if kind == "latitude" else "X"
    candidates = [
        name
        for name, variable in dataset.variables.items()
        if variable.attrs.get("standard_name") == standard_name
        or variable.attrs.get("axis") == axis
    ]
    if len(candidates) != 1:
        raise ValueError(f"could not uniquely identify {kind} coordinate")
    return candidates[0]


def _validate_coordinate_units(variable: xr.DataArray, kind: str) -> None:
    units = variable.attrs.get("units")
    if units is None:
        raise ValueError(f"{kind} coordinate is missing degree units")
    normalized = str(units).strip().lower().replace(" ", "_")
    allowed = (
        _LATITUDE_DEGREE_UNITS
        if kind == "latitude"
        else _LONGITUDE_DEGREE_UNITS
    )
    if normalized not in allowed:
        raise ValueError(
            f"{kind} coordinate units must be degrees, found {units!r}"
        )


def grid_from_dataset(dataset: xr.Dataset, variable_id: str = "mrro") -> RectilinearGrid:
    """Read the indexed horizontal grid without loading the runoff field."""

    if variable_id not in dataset:
        raise ValueError(f"variable not found: {variable_id}")
    latitude_name = _coordinate_name(dataset, "latitude")
    longitude_name = _coordinate_name(dataset, "longitude")
    latitude = dataset[latitude_name]
    longitude = dataset[longitude_name]
    if latitude.ndim != 1 or longitude.ndim != 1:
        raise ValueError("only one-dimensional rectilinear coordinates are supported")
    _validate_coordinate_units(latitude, "latitude")
    _validate_coordinate_units(longitude, "longitude")
    variable = dataset[variable_id]
    if latitude_name not in variable.dims or longitude_name not in variable.dims:
        raise ValueError(f"{variable_id} does not use the identified horizontal grid")

    lat_bounds_name = latitude.attrs.get("bounds")
    lon_bounds_name = longitude.attrs.get("bounds")
    lat_bounds = dataset[lat_bounds_name].values if lat_bounds_name else None
    lon_bounds = dataset[lon_bounds_name].values if lon_bounds_name else None

    return build_rectilinear_grid(
        latitude.values,
        longitude.values,
        latitude_bounds=lat_bounds,
        longitude_bounds=lon_bounds,
        latitude_name=latitude_name,
        longitude_name=longitude_name,
        latitude_bounds_source=f"cf:{lat_bounds_name}" if lat_bounds_name else None,
        longitude_bounds_source=f"cf:{lon_bounds_name}" if lon_bounds_name else None,
    )


def read_rectilinear_grid(
    path: str | Path, variable_id: str = "mrro"
) -> RectilinearGrid:
    """Open a NetCDF file and read only its horizontal-grid metadata."""

    with xr.open_dataset(
        Path(path), decode_times=False, mask_and_scale=False
    ) as dataset:
        return grid_from_dataset(dataset, variable_id=variable_id)
