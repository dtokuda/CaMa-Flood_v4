# Water-temperature and ice heat-budget stability and diagnostics

## Stabilization

When nearly all water leaves a cell, dividing residual heat by residual water volume can produce extreme temperatures through roundoff. Reconstruct the temperature using nonnegative weights for retained and incoming water. Retain the original temperature when there is no inflow. This reconstruction does not modify hydraulic storage or discharge.

Retain the previous temperature when liquid volume is at or below the existing `STO_IGNORE` threshold. During dry conditions this is a remembered temperature, not a heat source for ice. After rewetting, determine the temperature from actual inflowing or melted water and its sensible heat.

Pass liquid sensible heat directly between freezing and melting operations, avoiding repeated conversion to temperature at very small volumes. Atmospheric heating can still melt ice in a dry cell. Ice transport also avoids dividing by liquid volumes at or below `STO_IGNORE`.

Initialize water sensible-heat and ice export limiters to one. Divide only when requested export exceeds availability. This avoids unnecessary quotient overflow for tiny positive transport without rounding all such fluxes to zero.

In detailed mode, record the difference between expected energy and the energy represented by the retained temperature as signed unapplied heat. Do not carry it forward, redistribute it, or reinject it into model state. An unchanged dry cell does not repeatedly accumulate the same previously unapplied heat.

## Monitoring settings

```fortran
&NHEATLINK
    LICE = .TRUE.
    LHEAT_DIAG = .FALSE.
    CHEAT_LOG = 'log_HEAT-LINK.txt'
/
```

`LHEAT_DIAG` defaults to `.FALSE.`. Set it to `.TRUE.` for the full heat-budget audit, including every internal advection step. `CHEAT_LOG` defaults to `log_HEAT-LINK.txt` in the run directory; a custom relative or absolute path is accepted. Its parent directory must exist. Each run replaces this log. An empty path, an inaccessible path or an already-open file (including the CaMa log) is an error. With `LHEATLINK = .FALSE.`, no heat log is opened.

| Check or operation | Always active | Additional work with `LHEAT_DIAG = .TRUE.` |
|---|---|---|
| Dry-state handling, mixing reconstruction and bounded export limiters | All physical safeguards | No change to physical state |
| Existing invalid-input and phase/storage consistency checks | Errors remain fatal | No change |
| End-of-outer-update water temperature | Finite-value check; wet/dry minima, maxima, maximum cell index and volume; warning above 350 K | Same checks |
| Calendar timestamps and heatlink errors | Separate heat log | Same log |
| Interval heat closure and independent extrema | Not evaluated | Advection, local and combined outer update: ΔE, Q, U, residuals, ratios and extrema |
| Cause-specific unapplied heat and cumulative closure | Not aggregated | Full `HEAT_RESIDUAL`, `HEAT_CLOSURE` and `HEAT_CONSERVATION` records |
| Transport/phase budget summaries | Not aggregated | Cell and domain maximum-error summaries |

The switch skips diagnostic snapshots, compensated domain reductions, cause scans and their log records; it is not just an output filter. The physical phase kernel still calculates its existing residual outputs for `RIVICE_*` output requests, independently of this switch. Minimal mode does **not** establish heat-budget closure: use detailed mode when auditing conservation. Neither mode redistributes or reinjects unapplied heat.

## Module responsibilities

- `heat_residual_mod` defines both types and operations for cause-specific unapplied-heat totals and the domain heat ledger accumulated since this run started. These are supporting diagnostics.
- `heat_step_monitor_mod` defines both types and operations for interval snapshots, compensated sums, state increments and independent extrema. Interval residuals avoid relying on subtraction of large annual totals.
- `heatlink_diagnostics_mod` connects those two diagnostic modules to canonical CaMa storage and read-only river temperature/ice arrays. It owns diagnostic state, assembles exchanges and writes the logs.
- `heatlink_log_mod` owns the separate log file and formats calendar markers supplied by the CaMa driver state.
- `heatlink_river_mod` retains physical state, forcing, solver order and output coupling, with calls to begin, record and finish the diagnostic intervals.

The first two modules are separated by diagnostic responsibility and time scale, not by separating type definitions from their use. Their numerical routines can be tested independently of the river driver.

## Heat-budget definitions

Represented energy consists of water sensible heat relative to the melting point and the latent energy of ice at the melting point. Ice-surface temperature is a massless diagnostic; it does not represent sensible heat stored inside ice.

- `delta[J]`: represented energy increment ΔE during the interval.
- `net[J]`: net external heat input Q during the same interval. Internal transport is excluded.
- `exchange_abs[J]`: absolute external heat exchange S during the interval. For local heating, sum the absolute **net input per cell**, not the absolute values of every radiation and turbulent flux component.
- `unapplied[J]`: signed unapplied heat U. Positive means unapplied heating; negative means unapplied cooling.
- `unapplied_abs[J]`: absolute unapplied heat without cancellation between its positive and negative contributions.
- `raw[J] = ΔE − Q`: residual before accounting for unapplied heat.
- `adjusted[J] = ΔE − Q + U`: residual after accounting for unapplied heat in the diagnostic ledger.
- `storage_scale[J]`: an energy-storage scale in which sensible heat and ice latent energy do not cancel.

Compute ΔE directly from cellwise volume, temperature and ice increments. Subtracting two domain energy totals is not the primary diagnostic. Use `JPRD` Neumaier compensated sums. The comparison fields `naive_delta[J]` and `naive_adjusted[J]` retain the large-total subtraction.

Record `raw/S`, `adjusted/S`, absolute unapplied heat/S and `adjusted/storage_scale`. Flag zero-denominator and unrepresentable ratios as invalid and exclude them from ratio extrema. A zero placeholder in such a field is not a valid zero ratio. Energy and ratio extrema are selected independently and may occur at different times.

## Reading the log

Heatlink-specific messages go to `CHEAT_LOG` (normally `log_HEAT-LINK.txt`), while CaMa and shared input/output messages remain in `log_CaMa.txt`. Detailed mode writes heat-budget column definitions during initialization. Minimal mode writes its monitoring setting explicitly and retains temperature records.

| Record | Meaning |
|---|---|
| `HEAT_TIME` | Model calendar marker: stage, CaMa step counter, YYYYMMDD and HHMM |
| `HEAT_STEP advection` | Domain heat budget for each internal advection step |
| `HEAT_STEP local` | Domain heat budget for each local heat update |
| `HEAT_STEP hour` | Combined advection and local budget over an outer update interval; one hour in the annual validation |
| `HEAT_STEP_COUNTS` | Number of intervals and invalid ratios for each process |
| `HEAT_STEP_EXTREME` | Minimum/maximum of each metric, interval index and the complete corresponding diagnostic record, written at shutdown |
| `HEAT_RESIDUAL` | Cause-specific unapplied heat: dry advection, advection reconstruction, dry local update, ice handling and the no-ice melting-point floor |
| `HEAT_CONSERVATION` | Domain budget since this run started; cumulative supporting information |
| `HEAT_CLOSURE` | Process closure information; overlaps the domain residual and must not be added to it |
| `THERMAL_WET` / `THERMAL_DRY` | Wet/dry temperature ranges and the cell index of the maximum temperature |

`HEAT_TIME BEGIN` records the current CaMa calendar before hydraulic advection. `HEAT_TIME LOCAL_END` identifies the target date/time of the local update; `HEAT_TIME END` follows the completed heat update and uses the same `JYYYYMMDD` and `JHHMM` as `CMF::DRV_ADVANCE END` in the CaMa log. The counter follows CaMa's own update convention. These markers include day/year rollover and bracket the internal-step records. Each completed outer update flushes the heat log.

In detailed mode, `elapsed_seconds` accumulates internal time steps from the start of this run and `dt_seconds` is the monitored interval length. These distinguish adaptive internal steps within the calendar markers; they are not wall-clock timestamps. Allow for accumulated internal-step roundoff when converting extrema times back to calendar dates. Extrema are temporal extrema of interval **domain** budgets, not spatial extrema of cellwise errors. Restart runs start new diagnostic accumulations.

For example, inspect interval records and final extrema with:

```sh
rg '^HEAT_STEP (advection|local|hour) ' /path/to/run/log_HEAT-LINK.txt
rg '^HEAT_STEP_(COUNTS|EXTREME) ' /path/to/run/log_HEAT-LINK.txt
rg '^(HEAT_RESIDUAL|HEAT_CONSERVATION|HEAT_CLOSURE|THERMAL_)' /path/to/run/log_HEAT-LINK.txt
```

Detailed monitoring evaluates and logs every internal transport interval, requiring extra memory, arithmetic and log I/O. It is disabled by default. It never corrects physical state using diagnostic values. Runtime comparisons should use identical compiler options, forcing, time intervals, output settings and OpenMP configuration; concurrent annual runs alone do not isolate monitoring cost.

With `LICE = .FALSE.`, unapplied cooling caused by the melting-point floor remains in `raw` and U. A small `adjusted` residual does not mean that cooling was physically applied. With `LICE = .TRUE.`, evaluate dry-state and other unapplied heat separately as well. These domain diagnostics do not rule out cancellation between cells or account for heat fluxes never computed because of dry-state checks.

## Measured monitoring cost

A matched short-run benchmark on Apple M4 Pro, GNU Fortran 13.2.0, double precision, `-O3`, NetCDF and eight OpenMP threads (`OMP_WAIT_POLICY=PASSIVE`, `OMP_DYNAMIC=FALSE`) used MIROC6 forcing for 1–3 January 2000. Each executable/ice setting had one warm-up, followed by three serial repetitions in rotating order. Median wall time includes startup, initialization and I/O.

| LICE | PR2 [s] | PR3 detail OFF [s] | Relative to PR2 | PR3 detail ON [s] | Relative to PR2 |
|---|---:|---:|---:|---:|---:|
| TRUE | 18.892 | 17.904 | -5.2% | 25.160 | +33.2% |
| FALSE | 9.717 | 9.379 | -3.5% | 14.120 | +45.3% |

OFF also omits the pre-existing transport/phase summary reductions that PR2 always performed. These figures compare complete run configurations, not just the new arithmetic. They do not predict annual-run overhead. All output/restart binaries are bit-identical between OFF, ON and the preceding stabilized PR3 for both ice settings; detailed diagnostic records are unchanged. The separate heat log is about 22 kB with detail OFF versus 565–608 kB with detail ON for these three-day runs.

## Regression tests

Use the existing `adm/Mkinclude` or make command-line settings for the compiler and NetCDF environment. Tests do not require personal configuration paths or real forcing data. Example build commands with heatlink enabled follow; use the same make-variable overrides for each invocation.

```sh
make -r -j1 -C src all EXT_LIBS='common mod phys io heatlink' DHEATLINK=-Dheatlink
make -r -j1 -C src/heatlink test DHEATLINK=-Dheatlink
make -r -j1 -C src/phys test DHEATLINK=-Dheatlink
```

The existing `test` targets build executables; run them separately from their corresponding directories. In particular, `test_heatlink_config` reads relative fixture paths from `src/heatlink`.

```sh
(cd src/heatlink && ./test_heatlink_config)
(cd src/heatlink && ./test_heatlink_log)
python3 src/heatlink/test/test_heatlink_log_errors.py
(cd src/heatlink && ./test_temperature_dry_state)
(cd src/heatlink && ./test_thermo_dry_exchange)
(cd src/heatlink && ./test_transport_limiter)
(cd src/heatlink && ./test_heat_step_monitor)
(cd src/phys && ./test_heat_budget)
```

Configuration/log tests also check default/reset behavior, custom filenames, file ownership, invalid settings, already-open/inaccessible paths and calendar marker formatting across a year boundary. The added numerical tests cover small reproductions of observed anomalous states, values below/at/above the dry threshold, drying, rewetting, refreezing, melting, signed unapplied heat, tiny-flux overflow, boundary heat exchange, cancellation-prone synthetic states, invalid ratios and independent extrema. Run them alongside the existing water/ice transport, boundary and phase-change tests. `test_river_water_advection_cold_inflow` is expected to stop because liquid inflow is below the melting point.

When switching precision, clean and rebuild all dependencies so that objects and modules with different kinds are not mixed. `DSINGLE=-DSinglePrec_CMF` selects single precision; `DSINGLE=` selects double precision. The existing single-precision ice-surface Newton convergence test failure is a separate issue; this stabilization does not change its convergence criterion.

A three-day single-precision integration with ice disabled also produces bit-identical output with detail OFF/ON. With ice enabled, both modes stop at the first local update on the existing ice-surface Newton nonconvergence (maximum residual about `1.1292e-3 W m-2`, tolerance `1e-6 W m-2`); the preceding PR3 executable reproduces the same 56,429 nonconverged cells. This monitoring change does not alter that solver or its tolerance.
