"""Conservative spherical rectangle remapping and CaMa composition."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
from scipy import sparse

from .cama import CamaDimInfo, CamaInpmat
from .grid import RectilinearGrid, build_rectilinear_grid


def global_regular_grid(
    nlon: int,
    nlat: int,
    *,
    longitude_order: str = "west_to_east",
    latitude_order: str = "north_to_south",
) -> RectilinearGrid:
    """Construct an indexed, cell-centered global regular lat-lon grid."""

    if nlon < 1 or nlat < 1:
        raise ValueError("global regular grid dimensions must be positive")
    if longitude_order not in {"west_to_east", "east_to_west"}:
        raise ValueError(f"unsupported longitude order: {longitude_order}")
    if latitude_order not in {"north_to_south", "south_to_north"}:
        raise ValueError(f"unsupported latitude order: {latitude_order}")

    dx = 360.0 / nlon
    dy = 180.0 / nlat
    longitude = -180.0 + (np.arange(nlon, dtype=np.float64) + 0.5) * dx
    latitude = -90.0 + (np.arange(nlat, dtype=np.float64) + 0.5) * dy
    if longitude_order == "east_to_west":
        longitude = longitude[::-1]
    if latitude_order == "north_to_south":
        latitude = latitude[::-1]
    longitude_bounds = np.column_stack((longitude - 0.5 * dx, longitude + 0.5 * dx))
    latitude_bounds = np.column_stack((latitude - 0.5 * dy, latitude + 0.5 * dy))
    return build_rectilinear_grid(
        latitude,
        longitude,
        latitude_bounds=latitude_bounds,
        longitude_bounds=longitude_bounds,
        latitude_bounds_source="configured_global_regular_grid",
        longitude_bounds_source="configured_global_regular_grid",
    )


def reference_grid_from_diminfo(
    diminfo: CamaDimInfo, grid_config: dict[str, str]
) -> RectilinearGrid:
    """Build the indexed input grid expected by a configured reference inpmat."""

    expected_edges = (-180.0, 180.0, 90.0, -90.0)
    actual_edges = (diminfo.west, diminfo.east, diminfo.north, diminfo.south)
    if not np.allclose(actual_edges, expected_edges, rtol=0.0, atol=1.0e-7):
        raise ValueError(
            "reference diminfo must describe an unshifted global grid with "
            "edges west=-180, east=180, north=90, and south=-90 degrees"
        )

    return global_regular_grid(
        diminfo.nxin,
        diminfo.nyin,
        longitude_order=grid_config["longitude_order"],
        latitude_order=grid_config["latitude_order"],
    )


def standard_one_degree_grid() -> RectilinearGrid:
    """Return the conventional 360 x 180 CaMa reference input grid."""

    return global_regular_grid(360, 180)


def _ordinary_overlap_matrix(
    destination_bounds: np.ndarray,
    source_bounds: np.ndarray,
    *,
    transform=lambda values: values,
) -> sparse.csr_matrix:
    destination = transform(np.asarray(destination_bounds, dtype=np.float64))
    source = transform(np.asarray(source_bounds, dtype=np.float64))
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for destination_index, (destination_low, destination_high) in enumerate(destination):
        for source_index, (source_low, source_high) in enumerate(source):
            overlap = min(destination_high, source_high) - max(
                destination_low, source_low
            )
            if overlap > 0.0:
                rows.append(destination_index)
                columns.append(source_index)
                values.append(float(overlap))
    return sparse.coo_matrix(
        (values, (rows, columns)),
        shape=(destination.shape[0], source.shape[0]),
        dtype=np.float64,
    ).tocsr()


def _periodic_segments(bounds: np.ndarray) -> list[tuple[float, float, int]]:
    segments: list[tuple[float, float, int]] = []
    for index, (west, east) in enumerate(bounds):
        if east <= 360.0:
            segments.append((float(west), float(east), index))
        else:
            segments.append((float(west), 360.0, index))
            segments.append((0.0, float(east - 360.0), index))
    return [segment for segment in segments if segment[1] > segment[0]]


def _periodic_overlap_matrix(
    destination_bounds: np.ndarray, source_bounds: np.ndarray
) -> sparse.csr_matrix:
    destination_segments = _periodic_segments(destination_bounds)
    source_segments = _periodic_segments(source_bounds)
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for destination_low, destination_high, destination_index in destination_segments:
        for source_low, source_high, source_index in source_segments:
            overlap = min(destination_high, source_high) - max(
                destination_low, source_low
            )
            if overlap > 0.0:
                rows.append(destination_index)
                columns.append(source_index)
                values.append(float(overlap))
    matrix = sparse.coo_matrix(
        (values, (rows, columns)),
        shape=(destination_bounds.shape[0], source_bounds.shape[0]),
        dtype=np.float64,
    )
    matrix.sum_duplicates()
    return matrix.tocsr()


def conservative_remap_matrix(
    source: RectilinearGrid,
    destination: RectilinearGrid | None = None,
) -> sparse.csr_matrix:
    """Return R[dest, source] = spherical overlap / destination area."""

    target = destination or standard_one_degree_grid()
    sine = lambda values: np.sin(np.radians(values))
    latitude_overlap = _ordinary_overlap_matrix(
        target.latitude_bounds, source.latitude_bounds, transform=sine
    )
    longitude_overlap = _periodic_overlap_matrix(
        target.longitude_bounds, source.longitude_bounds
    )

    target_latitude_width = np.diff(sine(target.latitude_bounds), axis=1)[:, 0]
    target_longitude_width = np.diff(target.longitude_bounds, axis=1)[:, 0]
    if np.any(target_latitude_width <= 0.0) or np.any(target_longitude_width <= 0.0):
        raise ValueError("destination grid contains a non-positive cell area")

    latitude_weights = sparse.diags(1.0 / target_latitude_width) @ latitude_overlap
    longitude_weights = sparse.diags(1.0 / target_longitude_width) @ longitude_overlap
    result = sparse.kron(latitude_weights, longitude_weights, format="csr")
    result.sum_duplicates()
    result.eliminate_zeros()
    result.sort_indices()
    return result


def grid_cell_area_factors(grid: RectilinearGrid) -> np.ndarray:
    """Return cell areas without the common Earth-radius-squared factor."""

    latitude_factor = np.diff(
        np.sin(np.radians(grid.latitude_bounds)), axis=1
    )[:, 0]
    longitude_factor = np.radians(
        np.diff(grid.longitude_bounds, axis=1)[:, 0]
    )
    return np.outer(latitude_factor, longitude_factor).ravel()


def compose_mapping(
    reference_to_cama: sparse.spmatrix, native_to_reference: sparse.spmatrix
) -> sparse.csr_matrix:
    """Compose CaMa overlap areas with conservative reference-grid means."""

    if reference_to_cama.shape[1] != native_to_reference.shape[0]:
        raise ValueError(
            "mapping dimensions do not match: "
            f"{reference_to_cama.shape} @ {native_to_reference.shape}"
        )
    result = (reference_to_cama @ native_to_reference).tocsr()
    result.sum_duplicates()
    result.eliminate_zeros()
    result.sort_indices()
    if result.nnz and np.min(result.data) < 0.0:
        raise ValueError("composed mapping contains negative effective areas")
    return result


def normalize_mapping_rows(matrix: sparse.spmatrix) -> sparse.csr_matrix:
    """Normalize each nonempty nonnegative row to a dimensionless area mean."""

    result = matrix.tocsr(copy=True).astype(np.float64)
    result.sum_duplicates()
    result.eliminate_zeros()
    result.sort_indices()
    if result.nnz and (
        not np.all(np.isfinite(result.data)) or np.min(result.data) < 0.0
    ):
        raise ValueError("area-mean mapping contains invalid weights")
    row_sums = np.asarray(result.sum(axis=1)).ravel()
    inverse = np.zeros_like(row_sums)
    active = row_sums > 0.0
    inverse[active] = 1.0 / row_sums[active]
    result = (sparse.diags(inverse) @ result).tocsr()
    result.eliminate_zeros()
    result.sort_indices()
    return result


def area_mean_composed_mapping(
    reference_to_cama: sparse.spmatrix,
    native_to_reference: sparse.spmatrix,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """Return normalized reference-to-CaMa and composed area-mean operators."""

    reference_area_mean = normalize_mapping_rows(reference_to_cama)
    composed = compose_mapping(reference_area_mean, native_to_reference)
    return reference_area_mean, composed


def csr_to_cama_inpmat(
    matrix: sparse.spmatrix,
    base_diminfo: CamaDimInfo,
    source_grid: RectilinearGrid,
    *,
    inpmat_name: str = "inpmat.bin",
) -> CamaInpmat:
    """Encode a CaMa x native CSR matrix as padded inpx/inpy/inpa records."""

    csr = matrix.tocsr(copy=True)
    csr.sum_duplicates()
    csr.eliminate_zeros()
    csr.sort_indices()
    nlat, nlon = source_grid.shape
    expected_shape = (base_diminfo.nx * base_diminfo.ny, nlat * nlon)
    if csr.shape != expected_shape:
        raise ValueError(f"matrix has shape {csr.shape}, expected {expected_shape}")
    counts = np.diff(csr.indptr)
    inpn = int(np.max(counts, initial=0))
    if inpn <= 0:
        raise ValueError("composed mapping is empty")

    shape = (inpn, base_diminfo.ny, base_diminfo.nx)
    inpx = np.zeros(shape, dtype=np.int32)
    inpy = np.zeros(shape, dtype=np.int32)
    inpa = np.zeros(shape, dtype=np.float32)
    for row in np.flatnonzero(counts):
        start, end = csr.indptr[row], csr.indptr[row + 1]
        columns = csr.indices[start:end]
        values = csr.data[start:end]
        iy, ix = divmod(int(row), base_diminfo.nx)
        count = end - start
        inpx[:count, iy, ix] = columns % nlon + 1
        inpy[:count, iy, ix] = columns // nlon + 1
        inpa[:count, iy, ix] = values

    diminfo = replace(
        base_diminfo,
        nxin=nlon,
        nyin=nlat,
        inpn=inpn,
        inpmat_name=inpmat_name,
        west=-180.0,
        east=180.0,
        north=90.0,
        south=-90.0,
    )
    result = CamaInpmat(diminfo=diminfo, inpx=inpx, inpy=inpy, inpa=inpa)
    result.validate()
    return result


def relative_linf(actual: np.ndarray, expected: np.ndarray) -> float:
    difference = float(np.max(np.abs(actual - expected), initial=0.0))
    scale = max(float(np.max(np.abs(expected), initial=0.0)), np.finfo(float).tiny)
    return difference / scale


def validate_mapping(
    reference_to_cama: sparse.csr_matrix,
    native_to_reference: sparse.csr_matrix,
    composed: sparse.csr_matrix,
    emitted: CamaInpmat,
    source_grid: RectilinearGrid,
    destination_grid: RectilinearGrid | None = None,
    *,
    random_seed: int = 20260809,
) -> dict[str, float | int | list[int]]:
    """Run algebra, constant-field, coverage, and global-volume checks."""

    target = destination_grid or standard_one_degree_grid()
    row_sums = np.asarray(native_to_reference.sum(axis=1)).ravel()
    coverage_error = float(np.max(np.abs(row_sums - 1.0)))

    rng = np.random.default_rng(random_seed)
    field = rng.random(native_to_reference.shape[1])
    two_step = reference_to_cama @ (native_to_reference @ field)
    composed_result = composed @ field
    algebra_error = relative_linf(composed_result, two_step)

    constant_reference = np.asarray(reference_to_cama.sum(axis=1)).ravel()
    constant_composed = np.asarray(composed.sum(axis=1)).ravel()
    constant_error = relative_linf(constant_composed, constant_reference)

    source_area = grid_cell_area_factors(source_grid)
    target_area = grid_cell_area_factors(target)
    remapped = native_to_reference @ field
    source_volume = float(np.dot(source_area, field))
    target_volume = float(np.dot(target_area, remapped))
    mass_error = abs(target_volume - source_volume) / max(
        abs(source_volume), np.finfo(float).tiny
    )

    emitted_matrix = emitted.to_csr()
    emitted_error = relative_linf(emitted_matrix @ field, composed_result)

    stressed = rng.normal(size=native_to_reference.shape[1])
    stressed[::97] = np.nan
    preprocessed = np.maximum(np.nan_to_num(stressed, nan=0.0), 0.0)
    preprocessing_error = relative_linf(
        emitted_matrix @ preprocessed, composed @ preprocessed
    )

    return {
        "r_shape": list(native_to_reference.shape),
        "r_nnz": int(native_to_reference.nnz),
        "r_coverage_max_abs_error": coverage_error,
        "c_shape": list(composed.shape),
        "c_nnz": int(composed.nnz),
        "inpn": emitted.diminfo.inpn,
        "algebra_relative_linf": algebra_error,
        "constant_relative_linf": constant_error,
        "global_mass_relative_error": mass_error,
        "float32_emission_relative_linf": emitted_error,
        "preprocessed_float32_relative_linf": preprocessing_error,
    }


def validate_area_mean_mapping(
    reference_area_mean: sparse.csr_matrix,
    native_to_reference: sparse.csr_matrix,
    composed: sparse.csr_matrix,
    emitted: CamaInpmat,
    *,
    random_seed: int = 20260810,
) -> dict[str, float | int | list[int]]:
    """Validate normalization, composition, convexity, and binary encoding."""

    reference_row_sums = np.asarray(reference_area_mean.sum(axis=1)).ravel()
    composed_row_sums = np.asarray(composed.sum(axis=1)).ravel()
    active = reference_row_sums > 0.0
    if not np.any(active):
        raise ValueError("reference area-mean mapping has no active rows")
    reference_row_sum_error = float(
        np.max(np.abs(reference_row_sums[active] - 1.0), initial=0.0)
    )
    composed_row_sum_error = float(
        np.max(np.abs(composed_row_sums[active] - 1.0), initial=0.0)
    )

    rng = np.random.default_rng(random_seed)
    field = rng.normal(size=native_to_reference.shape[1])
    explicit = reference_area_mean @ (native_to_reference @ field)
    composed_result = composed @ field
    algebra_error = relative_linf(composed_result, explicit)

    constant = np.ones(native_to_reference.shape[1], dtype=np.float64)
    constant_result = composed @ constant
    constant_error = float(
        np.max(np.abs(constant_result[active] - 1.0), initial=0.0)
    )

    coo = composed.tocoo()
    row_min = np.full(composed.shape[0], np.inf, dtype=np.float64)
    row_max = np.full(composed.shape[0], -np.inf, dtype=np.float64)
    np.minimum.at(row_min, coo.row, field[coo.col])
    np.maximum.at(row_max, coo.row, field[coo.col])
    lower_violation = row_min[active] - composed_result[active]
    upper_violation = composed_result[active] - row_max[active]
    range_scale = max(float(np.ptp(field)), np.finfo(float).tiny)
    range_violation = max(
        float(np.max(lower_violation, initial=0.0)),
        float(np.max(upper_violation, initial=0.0)),
        0.0,
    ) / range_scale

    emitted_matrix = emitted.to_csr()
    emitted_error = relative_linf(emitted_matrix @ field, composed_result)
    emitted_row_sums = np.asarray(emitted_matrix.sum(axis=1)).ravel()
    emitted_constant_error = float(
        np.max(np.abs(emitted_row_sums[active] - 1.0), initial=0.0)
    )

    return {
        "reference_shape": list(reference_area_mean.shape),
        "reference_nnz": int(reference_area_mean.nnz),
        "reference_row_sum_max_abs_error": reference_row_sum_error,
        "c_shape": list(composed.shape),
        "c_nnz": int(composed.nnz),
        "inpn": emitted.diminfo.inpn,
        "composed_row_sum_max_abs_error": composed_row_sum_error,
        "constant_max_abs_error": constant_error,
        "algebra_relative_linf": algebra_error,
        "convex_range_relative_violation": range_violation,
        "float32_emission_relative_linf": emitted_error,
        "float32_constant_max_abs_error": emitted_constant_error,
    }
