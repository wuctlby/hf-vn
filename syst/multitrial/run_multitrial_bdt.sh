#!/bin/bash

generate_cutset() {
    local yaml_file="$1"
    local pt_dir="$2"

    # Safety checks
    if [[ ! -f "$yaml_file" ]]; then
        echo "ERROR: YAML file not found: $yaml_file"
        return 1
    fi
    if [[ -z "$pt_dir" ]]; then
        echo "ERROR: pt_dir is empty"
        return 1
    fi

    # Extract trial number from filename
    local trial_num
    trial_num=$(basename "$yaml_file" | sed -E 's/config_trial_([0-9]+)\.yml/\1/')

    # Target directory for this trial from directory of config_trial
    # yaml_file_dir=$(dirname "$yaml_file")
    # local target_dir="$target_base_dir/$trial_num"
    local target_dir=$(dirname "$yaml_file")

    # Create target directory
    mkdir -p "$target_dir"

    # Log start
    echo "[$trial_num] Generating YAML cutsets for file $yaml_file into $target_dir"

    # Run Python script
    if python3 "$path_to_src/src/make_cutsets_cfgs.py" "$yaml_file" -o "$target_dir"; then
        echo "[$trial_num] Cutsets generated successfully"
    else
        echo "[$trial_num] Failed to generate cutsets"
        return 1
    fi
}

perform_projections() {
    local yaml_file="$1"
    trial_num=$(basename "$yaml_file" | sed -E 's/config_trial_([0-9]+)\.yml/\1/')
    echo "[$trial_num] Projections for file $yaml_file"
    config_cutset="$(dirname "$yaml_file")/cutsets/cutset_00.yml"
    python3 $path_to_src/src/proj_thn.py $yaml_file --cutsetConfig $config_cutset
}

compute_efficiencies() {
    local yaml_file="$1"
    trial_num=$(basename "$yaml_file" | sed -E 's/config_trial_([0-9]+)\.yml/\1/')
    echo "[$trial_num] Computing efficiencies for file $yaml_file"
    proj_file="$(dirname "$yaml_file")/projs/proj_00.root"
    python3 $path_to_src/src/compute_efficiencies.py $yaml_file $proj_file
}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export ROOT_MAX_THREADS=1

# Export function and necessary variables
export -f generate_cutset
export -f perform_projections
export -f compute_efficiencies
export dir  # make sure $dir is visible inside the function

export n_parallel=50
export n_trials=2000
export path_to_src="/home/mdicosta/alice/hf-vn/"
export config_modifies_fit="/home/mdicosta/DFlowOO/SP_bdt_OO/020_correct_ev_sels/syst_multitrial_bdt_scan.yml"
export config_default="/home/mdicosta/DFlowOO/SP_bdt_OO/020_correct_ev_sels/config_020_remove_outliars.yml"
export output_dir="/home/mdicosta/DFlowOO/SP_bdt_OO/020_correct_ev_sels/cutvar_020_remove_outliars_combined"
export do_cms_fits=false
export do_compile_fitter=false
export do_cutset_generation=true
export do_projections=true
export do_efficiencies=true
export do_sim_fits=true
export produce_scan_plots=true
# export do_move_and_rename=false
# export do_df_frac=false
# export do_vn_vs_frac=false
# export produce_trials_plots=false

mkdir -p "$output_dir/syst/multitrial/bdt"

# Sanity check
if [ "$do_cms_fits" = true ] && [ "$do_sim_fits" = true ]; then
    echo "Error: Cannot perform both CMS fits and simulation fits in the same run. Please choose only one."
    exit 1
fi

if [ "$do_cutset_generation" = true ]; then
    # Generate YAML file list with indices
    python3 $path_to_src/syst/multitrial/make_configs_multitrial.py $config_default -m $config_modifies_fit -bm -o $output_dir > "$output_dir/syst/multitrial/bdt/log_make_cutsets.txt" 2>&1
    echo "Yaml files generated!"
fi

if [ "$do_compile_fitter" = true ]; then
    echo "Compiling InvMassFitter and VnVsMassFitter ..."
    rootcling -f $path_to_src/invmassfitter/vnfitter_dict.cxx -c /home/mdicosta/alice/hf-vn/invmassfitter/InvMassFitter.h /home/mdicosta/alice/hf-vn/invmassfitter/VnVsMassFitter.h /home/mdicosta/alice/hf-vn/invmassfitter/LinkDefVnFitter.h

    echo "Compiling fitter once for every trial ..."
    g++ -shared -fPIC `root-config --cflags --libs` \
        $path_to_src/invmassfitter/InvMassFitter.cxx $path_to_src/invmassfitter/VnVsMassFitter.cxx $path_to_src/invmassfitter/vnfitter_dict.cxx \
        -o $path_to_src/invmassfitter/libvnfitter.so
    echo "Compilation done!"
fi

pt_dirs=($(ls -d "$output_dir"/syst/multitrial/bdt/pt_*))

# Find YAML files, sort numerically by trial number, one per line
for dir in "${pt_dirs[@]}"; do
    find "$dir" -maxdepth 3 -type f -name "config_trial_*.yml" \
    | sort -V \
    > "$dir/multitrial_configs.txt"
done

# Loop over pt_dirs for cutset generation, projections, fitting with simultaneous fits
for dir in "${pt_dirs[@]}"; do
    echo -e "\nProcessing directory $dir"

    # Record start time for directory
    dir_start=$(date +%s)

    # Find YAML files in the current directory
    yaml_files=($(find "$dir" -maxdepth 3 -type f -name "config_trial_*.yml"))

    # Count them
    echo "----> Number of yaml files found in $dir: $(wc -l < "$dir/multitrial_configs.txt")"

    # --- Perform cutset generation ---
    if [ "$do_cutset_generation" = true ]; then
        log_file_cutset="$dir/log_cutset_generation.txt"
        start=$(date +%s)
        parallel -j "$n_parallel" generate_cutset ::: "${yaml_files[@]}" ::: "$dir" > "$log_file_cutset" 2>&1
        end=$(date +%s)
        echo -e "✅ Cutset generation done in $((end - start)) seconds. Starting projections ..."
    fi

    # --- Perform projections ---
    if [ "$do_projections" = true ]; then
        log_file_projections="$dir/log_projections.txt"
        start=$(date +%s)
        parallel -j "$n_parallel" perform_projections ::: "${yaml_files[@]}" > "$log_file_projections" 2>&1
        end=$(date +%s)
        echo -e "✅ Projections performed in $((end - start)) seconds. Starting fitting ..."
    fi

    # --- Compute efficiencies ---
    if [ "$do_efficiencies" = true ]; then
        log_file_efficiencies="$dir/log_efficiencies.txt"
        start=$(date +%s)
        parallel -j "$n_parallel" compute_efficiencies ::: "${yaml_files[@]}" > "$log_file_efficiencies" 2>&1
        end=$(date +%s)
        echo -e "✅ Efficiencies calculated in $((end - start)) seconds. Starting fitting ..."
    fi

    # --- Perform fitting ---
    if [ "$do_sim_fits" = true ]; then
        log_file_fits="$dir/log_simfits.txt"
        start=$(date +%s)
        python3 $path_to_src/syst/multitrial/run_fits.py "${yaml_files[@]}" --nproc $n_parallel > "$log_file_fits" 2>&1
        end=$(date +%s)
        echo -e "✅ Fitting done in $((end - start)) seconds."
    fi
done

pt_dirs=("$output_dir"/syst/multitrial/bdt/pt_*)

if [ "$do_cms_fits" = true ]; then
    log_file_fits="$dir/log_yieldfits.txt"
    start=$(date +%s)
    parallel -j "$n_parallel" \
        python3 "$path_to_src/src/get_vn_by_yield_extraction.py" "$config_default" \
        --multitrial \
        --multitrial_configs "{}/multitrial_configs.txt" \
        '>' "{}/log_yieldfits.txt" '2>&1' \
        ::: "${pt_dirs[@]}"

    end=$(date +%s)
    echo -e "----> CMS Fitting done in $((end - start)) seconds. Starting v2 vs frac extraction ..."
fi

# --- Produce final plots ---
if [ "$produce_scan_plots" = true ]; then
    echo "Producing final systematic plots ..."
    python3 $path_to_src/syst/multitrial/produce_bdt_multitrial_syst_plots.py "$config_default" "$output_dir" > "$output_dir/syst/multitrial/bdt/log_produce_scan_plots.txt" 2>&1
    echo "Final plots produced!"
fi

echo -e "\n--- Starting Vn vs fraction extraction ---"
# Run Vn vs fraction extraction if requested
if [ "$do_vn_vs_frac" = true ] || [ "$do_move_and_rename" = true ] || [ "$do_df_frac" = true ]; then
    for dir in "${pt_dirs[@]}"; do
        echo -e "\nProcessing directory $dir for Vn vs fraction extraction"
        log_file_vn_vs_frac="$dir/log_vn_vs_frac.txt"

        echo "Extracting Vn vs fraction for directory $dir"
        start=$(date +%s)

        # Build flags dynamically
        flags=()
        [ "$do_move_and_rename" = true ] && flags+=(--move_files)
        [ "$do_df_frac" = true ] && flags+=(--do_df_frac)

        [ "$do_vn_vs_frac" = true ] && flags+=(--do_vn_vs_frac)

        echo "DEBUG CMD:"
        echo python3 "$path_to_src/syst/multitrial/run_vn_vs_frac.py" \
            "$config_default" "$dir" "${flags[@]}"

        echo "Running with flags: ${flags[*]}"
        # Set niceness to 19
        nice -n 19 python3 "$path_to_src/syst/multitrial/run_vn_vs_frac.py" \
                           "$config_default" "$dir" --n_trials $n_trials "${flags[@]}" > "$log_file_vn_vs_frac" 2>&1

        end=$(date +%s)
        echo "✅ Vn vs fraction extraction done in $((end - start)) seconds."


        # {
        #     echo "Extracting Vn vs fraction for directory $dir"
        #     start=$(date +%s)

        #     # Build flags dynamically
        #     flags=()
        #     [ "$do_move_and_rename" = true ] && flags+=(--move_files)
        #     [ "$do_df_frac" = true ] && flags+=(--do_df_frac)

        #     [ "$do_vn_vs_frac" = true ] && flags+=(--do_vn_vs_frac)

        #     echo "DEBUG CMD:"
        #     echo python3 "$path_to_src/syst/multitrial/run_vn_vs_frac.py" \
        #         "$config_default" "$dir" "${flags[@]}"

        #     echo "Running with flags: ${flags[*]}"
        #     python3 "$path_to_src/syst/multitrial/run_vn_vs_frac.py" \
        #         "$config_default" "$dir" "${flags[@]}"

        #     end=$(date +%s)
        #     echo "✅ Vn vs fraction extraction done in $((end - start)) seconds."
        # } > "$log_file_vn_vs_frac" 2>&1
    done
fi


# --- Produce final plots ---
if [ "$produce_trials_plots" = true ]; then
    echo "Producing final systematic plots ..."
    python3 $path_to_src/syst/multitrial/produce_fit_multitrial_syst_plots.py "$config_default" "$output_dir" --multitrial_type bdt > "$output_dir/syst/multitrial/bdt/log_produce_trial_plots.txt" 2>&1
    echo "Final plots produced!"
fi
