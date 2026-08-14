#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec "${script_dir}/run_cmip6_mrro_2001.sh" \
    CESM2 \
    /Users/dtokuda/work/data/cmip6/CMIP6/LS3MIP/NCAR/CESM2/land-hist/r1i1p1f1/day/mrro/gn/v20200124/mrro_day_CESM2_land-hist_r1i1p1f1_gn_19810102-20150101.nc \
    .FALSE.
