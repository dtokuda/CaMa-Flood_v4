#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec "${script_dir}/run_cmip6_mrro_2001.sh" \
    IPSL-CM6A-LR \
    /Users/dtokuda/work/data/cmip6/CMIP6/LS3MIP/IPSL/IPSL-CM6A-LR/land-hist/r1i1p1f1/day/mrro/gr/v20190708/mrro_day_IPSL-CM6A-LR_land-hist_r1i1p1f1_gr_18500101-20141231.nc \
    .TRUE.
