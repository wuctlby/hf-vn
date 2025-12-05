#!/bin/bash

INPUT_DIR=$1
WORKERS=$2

# Ensure INPUT_DIR is absolute
INPUT_DIR=$(realpath "$INPUT_DIR")

# Find all AO2D files
AODS=($(find "$INPUT_DIR" -name 'AO2D_*.root'))

if [ ${#AODS[@]} -eq 0 ]; then
    echo "No AO2D files found in $INPUT_DIR."
    exit 1
fi

mkdir -p "${INPUT_DIR}/compressed"
mkdir -p "${INPUT_DIR}/ao2ds_to_merge"

process_file() {
    local file=$1
    local input_dir=$2  # pass INPUT_DIR explicitly

    base_name=$(basename "$file" .root)

    # Write the merge list in absolute path
    merge_list="${input_dir}/ao2ds_to_merge/${base_name}.txt"
    echo "$file" > "$merge_list"

    # Output file with absolute path
    output_file="${input_dir}/compressed/compressed_${base_name}.root"

    echo "Merging $file into $output_file"
    o2-aod-merger --input "$merge_list" --output "$output_file" --max-size 1000000000000000 --verbosity 0
}

export -f process_file

# Pass INPUT_DIR explicitly to each parallel job

# Pipe files to parallel
find "$INPUT_DIR" -name 'AO2D_*.root' | parallel -j "$N_JOBS" process_file {} "$INPUT_DIR"

# parallel -j "$WORKERS" process_file {} "$INPUT_DIR" ::: "${AODS[@]}"



