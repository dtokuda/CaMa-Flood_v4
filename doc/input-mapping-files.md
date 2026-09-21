# Explicit mapping files for multi-variable inputs

Each `input_item` may select its own mapping with `diminfo_file` and
`inpmat_file`. Both are required for gridded inputs; omit both for
`is_catm=.true.` inputs. Filenames and parent directories are
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

## Compatibility scope

Runoff-only discharge calculations retain the existing `NFORCE/CINPMAT` and
`NDIMTIME/CDIMINFO` interfaces and readers. They do not use this multi-variable
mapping module, and do not require the new `input_item` fields. Standard
runoff-only run scripts and model routines are unchanged.

For multi-variable / heatlink inputs, the old `intrp_map` / `nml_inpmat`
registration, grid auto-selection and prefix-based `.dim/.cfg/.bin` reader
have been removed. Non-CaMa-grid inputs must specify both mapping files.
`is_catm=.true.` inputs bypass mapping and must omit both fields. An old
mapping directory must be replaced with appropriate diminfo/binary pairs;
simply renaming a `.dim` file does not convert its format.

The former heatlink shell variable `INPMAT_DIR` selected the old atmospheric
mapping directory only. It is removed; use `INPMAT_DIR_ATM` for the new
atmospheric mapping pairs. `INPMAT_DIR_LSM` in the example is a shell
convenience for the independent runoff `CDIMINFO` and `CINPMAT` paths.
This feature does not add NetCDF dimension slicing.

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

## Heatlink run script and namelist templates

The existing `gosh/etc/s01-simulation_heatlink.sh` expands `heat-link.nml`
and appends the atmospheric definitions from `atm_GSWP3.nml` before running
`src/MAIN_cmf`. There is no separate namelist-generation script to invoke.
Set forcing and mapping directories in the shell; filenames stay in the
namelists:

```fortran
&input_item item='TAIR', fmt='nc',
  path='@ATM_DIR@/GSWP3.BC.Tair.3hrMap.ILS.2000.nc',
  diminfo_file='@INPMAT_DIR_ATM@/diminfo.txt',
  inpmat_file='@INPMAT_DIR_ATM@/inpmat.bin' /
```

For example, from the repository root:

```sh
MAP_DIR=./map/glb_15min \
ATM_DIR=./forcing/gswp3 \
INPMAT_DIR_ATM=./mappings/atmosphere \
RUN_DIR=./out/heatlink-example \
bash gosh/etc/s01-simulation_heatlink.sh
```

`INPMAT_DIR_ATM` defaults to
`${MAP_DIR}/input_mappings/05deg_s-n_0e-360e/mean`. `ATM_DIR` is the forcing
data directory, not a mapping directory. `INPMAT_DIR_LSM` defaults to
`MAP_DIR`; the common template keeps the bundled binary-runoff filenames
`diminfo_test-1deg.txt` and `inpmat_test-1deg.bin` for `CDIMINFO` and `CINPMAT`.
`RUNOFF_DIR` still selects the binary-runoff data directory. The public example
retains the year-2000 annual run, forcing units, physics and output selections.

For several mapping pairs in the same directory, edit the individual
`diminfo_file` and `inpmat_file` filenames in the atmospheric template. You may
supply a different template through `NML_ATM` (or `NML_COMMON` for the common
settings). No per-file shell variables are required. Templates use
single-quoted paths; the script escapes embedded quotes and sed replacement
characters. Directory paths are resolved before changing to the run directory.
The script refuses a nonempty run directory and captures standard output,
standard error and elapsed-time reporting in `run_stdout.log` and
`run_stderr.log`. `OMP_NUM_THREADS`, if set, is inherited without imposing a
new thread-count default.

The example now uses explicit atmospheric mappings, so its old eight-entry
`intrp_map` / `nml_inpmat` list and the corresponding multi-variable Fortran
reader have been removed. The new mapping format requires the diminfo/binary pair described above; pointing at a legacy
`inpmat_01.dim/.cfg/.bin` directory does not convert it. Adapt template
filenames or generate the appropriate mapping pair before running.

## Tests

```sh
python3 gosh/etc/test/test_heatlink_script.py
make -C src/mod test-input-mapping FCMP=gfortran
make -C src/mod test-input-mapping FCMP=gfortran INPUT_MAPPING_TEST_FLAGS=--netcdf
```

The runner requires Python 3 and gfortran. The optional NetCDF mode uses
`nf-config` from PATH. Compiler options may be supplied through FCMP. Both
single and double model precision are tested with runtime bounds checks.
Build products and synthetic fixtures are isolated in a temporary directory;
no external model data are needed. On macOS, configure the compiler's SDK
normally (for example with SDKROOT if required by the toolchain).

The launcher tests use empty data directories and a fake executable to verify
namelist generation, log capture and failure handling without a model run.
