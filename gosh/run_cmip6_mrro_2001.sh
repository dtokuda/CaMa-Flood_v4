#!/usr/bin/env bash
# Run CaMa-Flood for calendar year 2001 with one CMIP6 LS3MIP mrro file.
# This helper is called by the model-specific scripts in this directory.

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 MODEL MRRO_FILE LLEAPYR" >&2
    exit 2
fi

model=$1
mrro_file=$2
lleapyr=$3

case "${lleapyr}" in
    .TRUE.|.FALSE.) ;;
    *)
        echo "LLEAPYR must be .TRUE. or .FALSE.: ${lleapyr}" >&2
        exit 2
        ;;
esac

cama_dir=/Users/dtokuda/work/model/CaMa-Flood_v4
map_dir=/Users/dtokuda/work/data/cmf_v420_pkg/map/glb_15min
program=${cama_dir}/src/MAIN_cmf
diminfo=${map_dir}/inpmat_comp-remap/${model}/diminfo.txt
inpmat=${map_dir}/inpmat_comp-remap/${model}/inpmat.bin
output_dir=${cama_dir}/out/${model}

for required_file in \
    "${program}" \
    "${mrro_file}" \
    "${diminfo}" \
    "${inpmat}" \
    "${map_dir}/nextxy.bin" \
    "${map_dir}/ctmare.bin" \
    "${map_dir}/elevtn.bin" \
    "${map_dir}/nxtdst.bin" \
    "${map_dir}/rivlen.bin" \
    "${map_dir}/fldhgt.bin" \
    "${map_dir}/rivwth_gwdlr.bin" \
    "${map_dir}/rivhgt.bin" \
    "${map_dir}/rivman.bin" \
    "${map_dir}/bifprm.txt"
do
    if [[ ! -r "${required_file}" ]]; then
        echo "Required file is not readable: ${required_file}" >&2
        exit 1
    fi
done

if [[ ! -x "${program}" ]]; then
    echo "CaMa-Flood executable is not executable: ${program}" >&2
    exit 1
fi

# Avoid silently mixing a new run with an earlier result. Move or remove the
# existing directory explicitly before rerunning a model.
if [[ -d "${output_dir}" ]] && [[ -n "$(find "${output_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    echo "Output directory is not empty: ${output_dir}" >&2
    exit 1
fi

mkdir -p "${output_dir}"

cat > "${output_dir}/input_cmf.nam" <<EOF
&NRUNVER
  LADPSTP  = .TRUE.
  LPTHOUT  = .FALSE.
  LRESTART = .FALSE.
  LOUTPUT  = .TRUE.
  LGRIDMAP = .TRUE.
  LLEAPYR  = ${lleapyr}
/

&NDIMTIME
  CDIMINFO = "${diminfo}"
  DT       = 86400
  IFRQ_INP = 24
/

&NPARAM
  PMANRIV = 0.03D0
  PMANFLD = 0.10D0
  PDSTMTH = 10000.D0
  PCADP   = 0.7D0
/

&NSIMTIME
  SYEAR = 2001
  SMON  = 1
  SDAY  = 1
  SHOUR = 0
  EYEAR = 2002
  EMON  = 1
  EDAY  = 1
  EHOUR = 0
/

&NMAP
  LMAPCDF = .FALSE.
  CNEXTXY = "${map_dir}/nextxy.bin"
  CGRAREA = "${map_dir}/ctmare.bin"
  CELEVTN = "${map_dir}/elevtn.bin"
  CNXTDST = "${map_dir}/nxtdst.bin"
  CRIVLEN = "${map_dir}/rivlen.bin"
  CFLDHGT = "${map_dir}/fldhgt.bin"
  CRIVWTH = "${map_dir}/rivwth_gwdlr.bin"
  CRIVHGT = "${map_dir}/rivhgt.bin"
  CRIVMAN = "${map_dir}/rivman.bin"
  CPTHOUT = "${map_dir}/bifprm.txt"
/

&NFORCE
  LINPCDF  = .TRUE.
  LINPDAY  = .FALSE.
  LINTERP  = .TRUE.
  LITRPCDF = .FALSE.
  CINPMAT  = "${inpmat}"
  CROFCDF  = "${mrro_file}"
  CVNTIME  = "time"
  CVNROF   = "mrro"
  DROFUNIT = 1000.D0
/

&NRESTART
  CRESTSTO = "unused"
  CRESTDIR = "./"
  CVNREST  = "restart_"
  LRESTCDF = .FALSE.
  LRESTDBL = .TRUE.
  IFRQ_RST = 0
/

&NOUTPUT
  COUTDIR  = "./"
  CVARSOUT = "outflw"
  COUTTAG  = "_2001"
  LOUTVEC  = .FALSE.
  LOUTCDF  = .FALSE.
  IFRQ_OUT = 24
/
EOF

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}

echo "Model:       ${model}"
echo "Runoff:      ${mrro_file}"
echo "Output:      ${output_dir}"
echo "OMP threads: ${OMP_NUM_THREADS}"

cd "${output_dir}"
"${program}" > run_stdout.log 2>&1

echo "Completed: ${model}"
echo "CaMa log:  ${output_dir}/log_CaMa.txt"
