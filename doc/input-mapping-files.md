# Explicit mapping files for multi-variable inputs

Each `input_item` may select its own mapping with `diminfo_file` and
`inpmat_file`. Supply both or neither. Filenames and parent directories are
independent; relative paths are resolved from the model's working directory.
The fields apply to both NetCDF and binary inputs managed by `InputConf`.

```fortran
&input_item item='TAIR', fmt='nc', path='./forcing/tair.nc',
  diminfo_file='./mappings/atm_mean.info',
  inpmat_file='./mappings/atm_mean.bin' /
&input_nc item='TAIR', var_name='tair' /
```

Several experiments can share one directory, or share a diminfo file and use
different binary weight files. The cache key is the pair of canonical file
paths, not the directory or source array shape. Relative aliases and symbolic
links reuse a cache entry. Files are loaded once per pair and must not be
modified during a run. Different hard-link names are not deduplicated.

## Compatibility

With neither field set, the existing `intrp_map` / `nml_inpmat` selection and
prefix-based `.dim`, `.cfg`, `.bin` files remain in use. Explicit and legacy
inputs can coexist; explicit mappings never enter legacy automatic matching.
An explicit-only configuration does not need an `intrp_map` group. Its explicit
mappings also work when a legacy group has `LINTRP=.false.`. Legacy inputs keep
their existing `LINTRP` behavior.

`is_catm=.true.` inputs bypass mapping; specifying mapping files for them is an
error. This feature does not change the separate runoff input path:
`NFORCE/CINPMAT` and `NDIMTIME/CDIMINFO` retain their existing meanings.
It does not add NetCDF dimension slicing.

## Mapping format and validation

The new reader accepts the existing 11-line diminfo layout used by the map
tools: destination NX, NY, NLFP; input NX, NY; INPN; embedded binary filename;
west, east, north, south bounds. Numeric lines may retain their trailing
comments. The embedded filename is descriptive: **`inpmat_file` selects the
actual file**, even when its basename differs.

The binary contains three contiguous blocks of 32-bit values in Fortran array
order `(NX, NY, INPN)`: integer X indices, integer Y indices, then real weights.
It must have the native byte order used by the existing mapping reader.
The reader checks destination dimensions, source array dimensions, file size,
index bounds, finite non-negative weights, consistent zero padding, and
contiguous active slots. Entirely unused trailing slots are allowed. Mapping
uses the existing weighted-mean operation and missing-value handling; it does
not choose a different physical remapping operation from the directory name.

Diminfo bounds must be finite. Source coordinate order and grid identity are
not fully described by this format: shape validation cannot distinguish two
grids of the same size with different coordinates or ordering. The caller must
select weights generated for the exact source and destination grids. No new
metadata format or inference from filenames is introduced.

## Shell configuration

Run `bash gosh/etc/input-mapping-example.sh` to print two input definitions.
It only writes a namelist fragment to standard output, without running the
model or reading data. Set `INPMAT_DIR_ATM` / `INPMAT_DIR_LSM` for convenient
shared filenames, or override `DIMINFO_FILE_ATM`, `INPMAT_FILE_ATM`,
`DIMINFO_FILE_LSM`, `INPMAT_FILE_LSM` individually. The example variable names
and input paths must be adapted to the actual forcing data.

```sh
INPMAT_DIR_ATM=./mappings/atm \
DIMINFO_FILE_LSM=./mappings/experiment_b.info \
INPMAT_FILE_LSM=./mappings/experiment_b.bin \
bash gosh/etc/input-mapping-example.sh > input-mapping.nml
```

Append the resulting definitions to a complete run namelist, replacing any
existing definitions of the same input items. Shell grouping into ATM/LSM is
only a convenience; the Fortran interface has no such fixed groups.

## Tests

```sh
make -C src/mod test-input-mapping FCMP=gfortran
make -C src/mod test-input-mapping FCMP=gfortran INPUT_MAPPING_TEST_FLAGS=--netcdf
```

The runner requires Python 3 and gfortran. The optional NetCDF mode uses
`nf-config` from PATH. Compiler options may be supplied through FCMP. Both
single and double model precision are tested with runtime bounds checks.
Build products and synthetic fixtures are isolated in a temporary directory;
no external model data are needed. On macOS, configure the compiler's SDK
normally (for example with SDKROOT if required by the toolchain).
