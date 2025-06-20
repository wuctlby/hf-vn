#!/bin/bash

generate_cutset() {
    local yaml_file="$1"
    trial_num=$(basename "$yaml_file" | sed -E 's/config_trial_([0-9]+)\.yml/\1/')
    target_dir="$dir/trial_$trial_num"
    # echo "[$trial_num] Generating YAML cutsets for file $yaml_file into $target_dir"
    python3 /home/mdicosta/alice/hf-vn/src/make_cutsets_cfgs.py "$yaml_file" -o "$target_dir"
}

fit_trial() {
    local yaml_file="$1"
    trial_num=$(basename "$yaml_file" | sed -E 's/config_trial_([0-9]+)\.yml/\1/')
    target_dir="$dir/trial_$trial_num"
    
    proj_files=($(find "$target_dir/proj/" -maxdepth 1 -type f -name "proj_*.root"))
    for proj_file in "${proj_files[@]}"; do
        printf "\n[%d] Fitting %s with config %s, output in %s\n" "$trial_num" "$proj_file" "$yaml_file" "$target_dir"
        python3 /home/mdicosta/alice/hf-vn/src/get_vn_vs_mass.py "$yaml_file" "$proj_file" --batch --multitrial
    done
}

# Export function and necessary variables
export -f generate_cutset
export -f fit_trial
export dir  # make sure $dir is visible inside the function

export n_parallel=20
export config_modifies="/home/mdicosta/alice/hf-vn/syst/multitrial/systMultitrial.yml"
export config_default="/home/mdicosta/alice/hf-vn/configs/test_config.yml"
export output_dir="/home/mdicosta/alice/hf-vn/test_output/cutvar_test_combined"
export multitrial_cutsets_script="/home/mdicosta/alice/hf-vn/syst/make_cutsets_multitrial.py"

# Generate YAML file list with indices
echo "Generating yaml files for systematics ..."
python3 ./syst/multitrial/make_cutsets_multitrial.py $config_default -m $config_modifies -o $output_dir
yaml_files=($(find "$output_dir" -type f -name "config_var_*.yml" -print0 | xargs -0 -n1))
echo "Yaml files generated!"
echo "${yaml_files[@]}"

pt_dirs=($(ls -d "$output_dir"/syst/multitrial/pt_*))
for dir in "${pt_dirs[@]}"; do
    echo "Processing directory $dir"

    # Record start time for directory
    dir_start=$(date +%s)

    # Find YAML files in the current directory
    yaml_files=($(find "$dir" -maxdepth 2 -type f -name "config_trial_*.yml"))
    echo "Number of yaml files found: ${#yaml_files[@]}"

    # --- Generate YAML cutsets ---
    start=$(date +%s)
    parallel -j $n_parallel generate_cutset ::: "${yaml_files[@]}"
    end=$(date +%s)
    echo -e "\n\n   Cutsets generated in $((end - start)) seconds.\n\n"

    # --- Perform projections ---
    start=$(date +%s)
    python3 /home/mdicosta/alice/hf-vn/src/proj_thn.py $config_default --multitrial_folder "$dir"
    end=$(date +%s)
    echo -e "\n\n   Projections performed in $((end - start)) seconds.\n\n"

    # --- Perform fitting ---
    start=$(date +%s)
    parallel -j $n_parallel fit_trial ::: "${yaml_files[@]}"
    end=$(date +%s)
    echo -e "\n\n   Fits performed in $((end - start)) seconds.\n\n"

    # Total time for this pt_dir
    dir_end=$(date +%s)
    echo "✅ Total time for $dir: $((dir_end - dir_start)) seconds."
    echo "-----------------------------------------------"

done