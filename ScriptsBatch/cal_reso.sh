#!/bin/bash
export SCRIPT_PATH=$(dirname "$0")/../src/compute_reso.py
echo "Running script at $SCRIPT_PATH"

export AN_RES_FILE="/home/wuct/MetaData/DATA/OO/apass2/resolution/corrected_evt_sels_pid/AnalysisResults.root"
export CENT_CLASS="k010 k020 k0100 k2050 k50100 k3050"
export VN_METHOD="sp"
export OUTPUT_DIR="/home/wuct/MetaData/DATA/OO/apass2/resolution/corrected_evt_sels_pid/"
export SUFFIX=""
export NWORKERS=8

function cal_reso() {
    local cent_class=$1
    local suffix=${SUFFIX}_${cent_class}
    local cmd="python3 $SCRIPT_PATH $AN_RES_FILE --centClass $cent_class --vn_method $VN_METHOD --outputdir $OUTPUT_DIR --suffix $suffix"
    echo "Executing: $cmd"
    eval "$cmd"
}
export -f cal_reso

echo $CENT_CLASS | tr ' ' '\n' | parallel -j $NWORKERS cal_reso {}