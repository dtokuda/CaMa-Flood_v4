#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec "${script_dir}/run_cmip6_mrro_2001.sh" \
    CNRM-ESM2-1 \
    /Users/dtokuda/work/data/cmip6/CMIP6/LS3MIP/CNRM-CERFACS/CNRM-ESM2-1/land-hist/r1i1p1f2/day/mrro/gr/v20190820/mrro_day_CNRM-ESM2-1_land-hist_r1i1p1f2_gr_18500101-20141231.nc \
    .TRUE.
