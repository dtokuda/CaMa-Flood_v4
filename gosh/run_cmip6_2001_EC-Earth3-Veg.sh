#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec "${script_dir}/run_cmip6_mrro_2001.sh" \
    EC-Earth3-Veg \
    /Users/dtokuda/work/data/cmip6/CMIP6/LS3MIP/EC-Earth-Consortium/EC-Earth3-Veg/land-hist/r1i1p1f1/day/mrro/gr/v20201007/mrro_day_EC-Earth3-Veg_land-hist_r1i1p1f1_gr_20010101-20011231.nc \
    .TRUE.
