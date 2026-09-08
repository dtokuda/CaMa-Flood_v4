# Composed remapping inpmat

This directory contains tools for mapping runoff from a supported native
rectilinear latitude-longitude grid to CaMa-Flood unit catchments.

The scientific mapping is

```text
native runoff grid --R--> configured reference grid --M--> CaMa-Flood
```

`R` is a conservative native-to-reference remapping operator. `M` is an
existing CaMa-Flood `inpmat` for that reference grid. The generated production
matrix stores the composed operator

```text
C = M R
```

CaMa-Flood can therefore read the native runoff directly. No intermediate
runoff time series on the reference grid is written.

The same generation also writes a dimensionless area-mean operator for scalar
intensive variables for which an area-weighted mean is scientifically
appropriate. If `normalize_rows(M)` divides each active CaMa row by its covered
area, this operator is

```text
C_area_mean = normalize_rows(M) R
```

It produces unit-catchment area-weighted means rather than area-integrated
runoff. It is a geometry-only operator for complete fields; dynamic or
land-only missing-value renormalization is not supported.

## Tutorial

Run the commands below from this directory.

### Python environment

The required packages are listed in `requirements.txt`:

```text
numpy
scipy
xarray
netCDF4
cftime
PyYAML
```


### Configuration

One YAML file describes one source grid and one output directory. The current
example is `config/cmip6_mrro.yaml`.

```yaml
version: 6  # Configuration schema version understood by this tool.

cama:
  executable: /path/to/MAIN_cmf  # CaMa-Flood executable used only by the integration test.
  map_dir: /path/to/cama/map  # Base directory for reference files, generated files, and CaMa map parameters.

reference_mapping:
  diminfo_path: diminfo_reference.txt  # Path relative to cama.map_dir; NXIN and NYIN define the reference-grid size.
  inpmat_path: inpmat_reference.bin  # Existing operator M; path relative to cama.map_dir.
  input_grid:
    type: global_regular_latlon  # Interpret the reference input grid as a cell-centered regular global grid.
    longitude_order: west_to_east  # Longitude index increases from west to east.
    latitude_order: north_to_south  # Latitude index increases from north to south.

grid:
  type: rectilinear  # Source latitude and longitude coordinates must both be one-dimensional.
  coverage: global  # Only complete global source grids are supported.
  bounds_policy: cf_or_derive  # Use CF bounds when present; otherwise derive Gaussian or midpoint bounds.
  hash_decimals: 7  # Decimal precision used to hash coordinates and bounds for grid identity.

preprocess:
  clip_negative: true  # Replace negative runoff with zero before applying spatial weights.
  use_areacella: false  # Do not use externally supplied native-grid cell areas.
  use_sftlf: false  # Do not use native-grid land fractions.
  normalize_over_valid_area: false  # Do not renormalize weights over valid or land-only source cells.

output:
  dirname: example-grid  # Directory created below cama.map_dir/inpmat_comp-remap.

input:
  path: /data/example/runoff_*.nc  # NetCDF path or glob used to read the source grid.
  variable: runoff  # Variable carrying the required 1-D latitude and longitude coordinates.

validation_forcing:
  use_input: true  # Validate input.path and input.variable as CaMa runoff forcing.
```

Coordinates must use degree units. `lat` and `lon` are preferred; otherwise CF
`standard_name` or `axis` attributes must identify them uniquely. CF bounds are
optional and are derived when absent.

Mapping generation needs only `input`. For forcing validation, either set
`validation_forcing.use_input: true` or specify its `path` and `variable`; the
grid hashes must match. Direct CaMa input requires `(time, lat, lon)`, no
packing, and a float32 missing marker of `1e20`. Units are reported, not
validated or converted.

### Fixed assumptions

- Source cells are spherical latitude-longitude rectangles in degrees.
- Missing/ocean runoff and negative runoff are treated as zero.
- Time frequency, time-coordinate cadence, and gaps or overlaps between files
  are not validated in this tool version. Daily-forcing checks belong to a
  separate change; users must currently select appropriate daily runoff.

The configured reference grid must be a cell-centered, regular global grid
with edges at -180, 180, -90, and 90 degrees. Its dimensions are read from the
configured reference `diminfo`; 1-degree and 0.5-degree reference mappings are
therefore supported, but a shifted or regional reference grid is not.

Cell areas and overlaps are calculated analytically from the spherical
rectangle bounds; no separate `areacella` field is needed. The conservative
weight is `overlap area / full destination-cell area`. For
example, if valid runoff `q` covers 30% of a destination cell and the remaining
70% is missing, the remapped value is `0.3 q`, not `q`.

## Commands

### Generate the composed inpmat

Optionally inspect every file selected by `input.path`:

```sh
python inspect_input_grid.py --config config/cmip6_mrro.yaml
```

The inventory is printed and appended to `validation.log`. Use `--output PATH`
only when a separate JSON copy is explicitly needed.

Generate the runoff and area-mean composed matrices, their CaMa-format files,
and the references needed by later validations:

```sh
python generate_native_inpmat.py --config config/cmip6_mrro.yaml
```

Generation checks coverage, composition, constant fields, global area
integrals, area-mean normalization and convexity, float32 encoding, and CaMa
binary indexing before writing files. It starts a new `validation.log` and
invalidates results from an older generation.

Generate the shell settings used to connect the new files to a CaMa run:

```sh
python emit_cama_config.py --config config/cmip6_mrro.yaml
```

This writes `cama-force.env`. Simulation-specific settings such as runoff unit
conversion are intentionally not included.

### Validation

Validation consists of three stages, ordered from the mapping layer to the full
CaMa executable.

1. Reload and validate the generated mapping:

   ```sh
   python validate_native_inpmat.py --config config/cmip6_mrro.yaml
   ```

   This checks the sparse operators, checksums, CaMa binary encoding, and the
   mathematical equivalence `M(Rx) = (MR)x`. It does not run CaMa-Flood.

2. Validate real runoff values:

   ```sh
   python validate_forcing_mapping.py --config config/cmip6_mrro.yaml
   ```

   This inventories and checks every matched validation-forcing file, then
   compares explicit `M(Rx)` and binary-inpmat `Cx` for the first time step of
   the first file and the last time step of the last file. It also checks the
   area-integrated native-to-reference runoff volume. It does not run
   CaMa-Flood.

3. Run the end-to-end CaMa integration validation:

   ```sh
   python validate_cama_integration.py \
     --config config/cmip6_mrro.yaml \
     --days 2
   ```

   This runs CaMa with both the explicit reference-grid path and the composed
   native-grid path, then compares `runoff`, `rivsto`, `storge`, and `rivout`.
   With two or more days it also validates restart continuity. Production
   validation requires this restart result. Temporary forcing, namelists, logs,
   CaMa outputs, and restart files are removed afterward.

After all three validation stages pass, run the aggregate consistency check:

```sh
python validate_production.py --config config/cmip6_mrro.yaml
```

It reruns the inventory and the first two stages, reads the retained integration
result, verifies a common grid hash, and appends the final result to
`validation.log`.

### Reading validation results

Every command appends a timestamped, human-readable JSON block to
`validation.log`. A successful command reports `"status": "passed"`; a failed
threshold check exits with an error instead of reporting a pass. The compact
machine-readable results needed by later validation stages are retained in
`validation_reference/validation.json`,
`validation_reference/validation_area_mean.json`, and
`validation_reference/integration_validation.json`.

The main generation metrics in `validation.json` are:

| Metric | Interpretation |
| --- | --- |
| `r_coverage_max_abs_error` | Maximum absolute difference between an `R` row sum and one. It checks that the native grid completely covers every reference-grid cell. |
| `global_mass_relative_error` | Relative difference between the native and reference-grid area integrals for a test field. This is the direct native-to-reference runoff-conservation check under the full rectangular-cell assumption. |
| `algebra_relative_linf` | Relative L-infinity difference between explicit `M(Rx)` and composed `(MR)x` in float64. It checks that matrix composition does not change the mapping. |
| `constant_relative_linf` | Relative L-infinity difference between `M1` and `C1`. It checks the constant-field response against the configured reference `inpmat`. |
| `float32_emission_relative_linf` | Error introduced by writing `C` to the float32 CaMa `inpmat.bin` and reading it back. |
| `preprocessed_float32_relative_linf` | The same binary comparison after applying the runoff rules for missing and negative values. |
| `inpn`, `r_nnz`, `c_nnz`, `r_shape`, `c_shape` | Structural diagnostics rather than accuracy metrics. They describe the sparse matrices and emitted CaMa dimensions. |

Each accuracy metric must be no larger than its value in the `thresholds`
object. Values near machine precision are expected for the float64 algebra and
mass checks; errors involving the emitted CaMa matrix are larger because
`inpa` is stored as float32.

`validation_area_mean.json` contains the corresponding checks for the
intensive-variable operator:

| Metric | Interpretation |
| --- | --- |
| `reference_row_sum_max_abs_error` | Maximum difference between one and a nonempty row sum of normalized `M`. |
| `composed_row_sum_max_abs_error` | Maximum difference between one and a nonempty row sum of `C_area_mean`. |
| `constant_max_abs_error` | Maximum error when mapping a constant field of one. |
| `algebra_relative_linf` | Difference between the explicit and composed area-mean routes. |
| `convex_range_relative_violation` | Any output excursion outside the range of the source cells contributing to that output. Zero is expected. |
| `float32_emission_relative_linf` | Error introduced by the CaMa-format float32 binary. |
| `float32_constant_max_abs_error` | Constant-field error after float32 binary encoding. |

These checks assume every source cell contains a valid intensive value. A
single fixed matrix cannot renormalize a land-only or time-varying valid mask;
such input must not be mapped with `inpmat_area_mean.bin` without a separate
numerator/coverage workflow.

`validate_forcing_mapping.py` records the following for the first and last
samples of the matched real-input collection:

| Field | Interpretation |
| --- | --- |
| `runtime_relative_linf` | Relative L-infinity difference between the explicit reference-grid route and the generated binary composed route. It must be at most `1e-6`. |
| `runtime_max_abs` | Largest absolute mapped-value difference, reported as a scale-dependent diagnostic. |
| `actual_native_to_reference_mass_relative_error` | Area-integrated runoff difference between the actual native field and `Rx`. This must be at most `1e-12`. |
| `global_mapped_sum_explicit`, `global_mapped_sum_composed` | Sums of the two CaMa-mapped fields. Their agreement checks composition on real data; these are not an independent full-globe mass balance because the reference `M` defines the CaMa-covered area. |
| `mapped_total_relative_error` | Relative difference between those two mapped totals. This must be at most `1e-6`. |
| `missing_count`, `negative_count` | Counts of input values affected by runoff preprocessing. |

`integration_validation.json` compares the two routes after running CaMa. For
each of `runoff`, `rivsto`, `storge`, and `rivout`, both `relative_linf` and
`relative_l1` must not exceed their corresponding entries in `thresholds`;
`max_abs` provides the scale-dependent absolute difference. The direct runoff
comparison retains a strict L-infinity threshold. Routed states use both a
global L1 diagnostic and a less brittle L-infinity threshold because float32
forcing roundoff can be locally amplified by nonlinear and adaptive routing.
When at least two days are requested, the `restart_*` entries additionally
compare a continuous run with a restarted run using strict thresholds.
`valid_values` is the number of CaMa values included in each comparison.

Runoff conservation should therefore be assessed from
`global_mass_relative_error`, followed by `algebra_relative_linf` and the real
forcing comparison. The first metric tests the area-integrated
native-to-reference volume, while the latter two demonstrate that the composed
binary operator preserves the configured reference-`inpmat` route.

## Output files

Each input grid is written below

```text
<cama.map_dir>/inpmat_comp-remap/<output.dirname>/
```

The top level contains:

| File | Purpose |
| --- | --- |
| `inpmat.bin` | Standard CaMa-Flood binary `inpmat` containing `C = M R`. |
| `diminfo.txt` | Standard CaMa-Flood dimensions paired with `inpmat.bin`. |
| `inpmat_area_mean.bin` | Dimensionless normalized weights containing `C_area_mean`. |
| `diminfo_area_mean.txt` | CaMa-format dimensions paired with `inpmat_area_mean.bin`. |
| `cama-force.env` | Shell settings emitted for a CaMa run. |
| `validation.log` | Consolidated command summaries and validation results. |

Files that are reread only during validation are isolated below
`validation_reference/`:

| File | Purpose |
| --- | --- |
| `native_to_reference_flux.npz` | Conservative sparse operator `R`. |
| `native_to_cama_flux.npz` | Sparse composed operator `C` before CaMa encoding. |
| `native_to_cama_area_mean.npz` | Composed area-mean operator before CaMa encoding. |
| `grid.json` | Source coordinates, bounds, ordering, and grid hash. |
| `validation.json` | Generation thresholds, metrics, and artifact checksums. |
| `validation_area_mean.json` | Area-mean assumptions, thresholds, metrics, and checksums. |
| `integration_validation.json` | CaMa output and restart comparisons read by the final check. |

The current runoff configuration connects CaMa-Flood only to `inpmat.bin` and
`diminfo.txt`. The area-mean pair is not automatically connected to a CaMa
forcing variable and requires a consumer that interprets `inpa` as normalized,
dimensionless weights. The tools do not create separate manifest, inventory,
forcing-validation, configuration-summary, or production-validation JSON
reports; human-readable results are consolidated in `validation.log`.

## Developer guide

The main implementation files are:

```text
cama_native_inpmat/config.py       YAML validation and path resolution
cama_native_inpmat/grid.py         Rectilinear coordinates, bounds, and grid hash
cama_native_inpmat/mapping.py      Conservative R, composition C = M R, validation
cama_native_inpmat/cama.py         CaMa diminfo and binary inpmat I/O
cama_native_inpmat/generate.py     Production and validation-reference artifacts
cama_native_inpmat/forcing.py      Runoff preprocessing and real-field validation
cama_native_inpmat/integration.py  End-to-end CaMa and restart validation
cama_native_inpmat/emit.py         cama-force.env generation
cama_native_inpmat/production.py   Final aggregate consistency check
cama_native_inpmat/reporting.py    Consolidated validation.log writer
```

The command-line entry points are intentionally thin wrappers around these
modules. Unit and local-data regression tests are under `tests/`:

```sh
python -m unittest discover -s tests -v
```
