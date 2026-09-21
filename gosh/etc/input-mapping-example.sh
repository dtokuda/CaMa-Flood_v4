#!/usr/bin/env bash
# Print a multi-input namelist fragment; do not run the model.
set -euo pipefail

INPMAT_DIR_ATM=${INPMAT_DIR_ATM:-./mappings/atm}
INPMAT_DIR_LSM=${INPMAT_DIR_LSM:-./mappings/lsm}
DIMINFO_NAME=${DIMINFO_NAME:-diminfo.txt}
INPMAT_NAME=${INPMAT_NAME:-inpmat.bin}
DIMINFO_FILE_ATM=${DIMINFO_FILE_ATM:-${INPMAT_DIR_ATM}/${DIMINFO_NAME}}
INPMAT_FILE_ATM=${INPMAT_FILE_ATM:-${INPMAT_DIR_ATM}/${INPMAT_NAME}}
DIMINFO_FILE_LSM=${DIMINFO_FILE_LSM:-${INPMAT_DIR_LSM}/${DIMINFO_NAME}}
INPMAT_FILE_LSM=${INPMAT_FILE_LSM:-${INPMAT_DIR_LSM}/${INPMAT_NAME}}
INPUT_FILE_ATM=${INPUT_FILE_ATM:-./forcing/tair.nc}
INPUT_FILE_LSM=${INPUT_FILE_LSM:-./forcing/runoff_temperature.nc}
INPUT_VAR_ATM=${INPUT_VAR_ATM:-tair}
INPUT_VAR_LSM=${INPUT_VAR_LSM:-runoff_temperature}

# Fortran character literals escape a quote by doubling it.
nml_quote() {
    local value=$1
    if [[ $value == *$'\n'* || $value == *$'\r'* ]]; then
        echo 'Namelist values must not contain newlines' >&2
        return 1
    fi
    value=${value//\'/\'\'}
    printf "'%s'" "$value"
}

emit_input() {
    local item=$1 path var diminfo inpmat
    path=$(nml_quote "$2")
    var=$(nml_quote "$3")
    diminfo=$(nml_quote "$4")
    inpmat=$(nml_quote "$5")
    printf "&input_item item='%s', fmt='nc', path=%s,\n  diminfo_file=%s, inpmat_file=%s /\n" \
        "$item" "$path" "$diminfo" "$inpmat"
    printf "&input_nc item='%s', var_name=%s /\n" "$item" "$var"
}

emit_input TAIR "$INPUT_FILE_ATM" "$INPUT_VAR_ATM" "$DIMINFO_FILE_ATM" "$INPMAT_FILE_ATM"
emit_input TROF "$INPUT_FILE_LSM" "$INPUT_VAR_LSM" "$DIMINFO_FILE_LSM" "$INPMAT_FILE_LSM"
