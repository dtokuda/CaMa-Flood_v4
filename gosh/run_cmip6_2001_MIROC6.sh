#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec "${script_dir}/run_cmip6_mrro_2001.sh" \
    MIROC6 \
    /Users/dtokuda/work/data/cmip6/CMIP6/LS3MIP/MIROC/MIROC6/land-hist/r1i1p1f1/day/mrro/gn/v20200501/mrro_day_MIROC6_land-hist_r1i1p1f1_gn_20000101-20091231.nc \
    .TRUE.
