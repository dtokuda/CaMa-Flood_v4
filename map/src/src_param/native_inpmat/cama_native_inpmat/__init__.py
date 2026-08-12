"""Utilities for native-grid CaMa-Flood runoff input matrices."""

from .cama import (
    CamaDimInfo,
    CamaInpmat,
    read_binary_inpmat,
    read_diminfo,
    write_binary_inpmat,
    write_diminfo,
)
from .grid import RectilinearGrid, build_rectilinear_grid, read_rectilinear_grid
from .mapping import (
    area_mean_composed_mapping,
    compose_mapping,
    conservative_remap_matrix,
    global_regular_grid,
    grid_cell_area_factors,
    normalize_mapping_rows,
    reference_grid_from_diminfo,
    standard_one_degree_grid,
)

__all__ = [
    "CamaDimInfo",
    "CamaInpmat",
    "RectilinearGrid",
    "build_rectilinear_grid",
    "area_mean_composed_mapping",
    "compose_mapping",
    "conservative_remap_matrix",
    "global_regular_grid",
    "grid_cell_area_factors",
    "normalize_mapping_rows",
    "read_binary_inpmat",
    "read_diminfo",
    "read_rectilinear_grid",
    "reference_grid_from_diminfo",
    "standard_one_degree_grid",
    "write_binary_inpmat",
    "write_diminfo",
]
