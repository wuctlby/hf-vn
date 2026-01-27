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

    # Target directory for this trial
    local target_dir="$pt_dir/trials/$trial_num"

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

get_v2_vs_fracs() {
    local yaml_file="$1"
    local reference_frac_dir="$2"
    trial_num=$(basename "$yaml_file" | sed -E 's/config_trial_([0-9]+)\.yml/\1/')
    target_dir="$dir/trials/$trial_num"
    echo "[$trial_num] Computing v2 vs frac for file $yaml_file"
    python3 $path_to_src/src/get_v2_vs_frac.py "$yaml_file" "$target_dir/raw_yields/" "$reference_frac_dir" --multitrial
}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export ROOT_MAX_THREADS=1

# Export function and necessary variables
export -f generate_cutset
export -f get_v2_vs_fracs
export dir  # make sure $dir is visible inside the function

export n_parallel=30
export home_dir="/home/mdicosta/"
export path_to_src="$home_dir/alice/hf-vn/"
export config_modifies_fit="$home_dir/DFlowOO/SP_bdt_OO/020_correct_ev_sels/syst_multitrial_simfit.yml"
export config_default="$home_dir/DFlowOO/SP_bdt_OO/020_correct_ev_sels/config_020_remove_outliars.yml"
export output_dir="$home_dir/DFlowOO/SP_bdt_OO/020_correct_ev_sels/cutvar_020_remove_outliars_combined"
export do_cms_fits=false
export do_compile_fitter=false
export do_cutset_generation=true
export do_projections=true
export do_sim_fits=true
export do_v2_vs_frac=true
export produce_plots=true

mkdir -p "$output_dir/syst/multitrial/fit"

# Sanity check
if [ "$do_cms_fits" = true ] && [ "$do_sim_fits" = true ]; then
    echo "Error: Cannot perform both CMS fits and simulation fits in the same run. Please choose only one."
    exit 1
fi

if [ "$do_cutset_generation" = true ]; then
    # Generate YAML file list with indices
    python3 $path_to_src/syst/multitrial/make_configs_multitrial.py $config_default -m $config_modifies_fit -fm -o $output_dir > "$output_dir/syst/multitrial/fit/log_make_cutsets.txt" 2>&1
    echo "Yaml files generated!"
fi

if [ "$do_compile_fitter" = true ]; then
    echo "Compiling InvMassFitter and VnVsMassFitter ..."
    rootcling -f $path_to_src/invmassfitter/vnfitter_dict.cxx -c $home_dir/alice/hf-vn/invmassfitter/InvMassFitter.h $home_dir/alice/hf-vn/invmassfitter/VnVsMassFitter.h $home_dir/alice/hf-vn/invmassfitter/LinkDefVnFitter.h

    echo "Compiling fitter once for every trial ..."
    g++ -shared -fPIC `root-config --cflags --libs` \
        $path_to_src/invmassfitter/InvMassFitter.cxx $path_to_src/invmassfitter/VnVsMassFitter.cxx $path_to_src/invmassfitter/vnfitter_dict.cxx \
        -o $path_to_src/invmassfitter/libvnfitter.so
    echo "Compilation done!"
fi

pt_dirs=($(ls -d "$output_dir"/syst/multitrial/fit/pt_*))

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

    # --- Perform projections ---
    if [ "$do_projections" = true ]; then
        log_file_projections="$dir/log_projections.txt"
        start=$(date +%s)
        python3 $path_to_src/src/proj_thn.py $config_default --multitrial_folder "$dir" --multitrial_workers $n_parallel > "$log_file_projections" 2>&1
        end=$(date +%s)
        echo -e "✅ Projections performed in $((end - start)) seconds. Starting fitting ..."
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

pt_dirs=("$output_dir"/syst/multitrial/fit/pt_*)

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

# --- Obtain v2 vs frac ---
if [ "$do_v2_vs_frac" = true ]; then
    for dir in "${pt_dirs[@]}"; do
        # Record start time for directory
        dir_start=$(date +%s)
        log_file_v2_vs_frac="$dir/log_v2_vs_frac.txt"
        yaml_files=($(find "$dir" -maxdepth 3 -type f -name "config_trial_*.yml"))
        reference_frac_dir="$output_dir/frac/"
        start=$(date +%s)
        parallel -j $n_parallel get_v2_vs_fracs ::: "${yaml_files[@]}" ::: "$reference_frac_dir" > "$log_file_v2_vs_frac" 2>&1
        # Total time for this pt_dir
        dir_end=$(date +%s)
        echo "-----------------------------------------------"
        echo -e "✅ V2 vs frac for $dir extracted in $((dir_end - dir_start)) seconds.\n"
    done
fi

# --- Evaluate systematics and produce final plots ---
export max_chi2=10
export min_signif=5
export max_signif=1000
export force_prompt_enhanced=false
if [ "$produce_plots" = true ]; then
    echo "Producing final systematic plots ..."
    mkdir -p "$output_dir/syst/multitrial/fit/summary"
    cmd=(
        python3 "$path_to_src/syst/multitrial/produce_fit_multitrial_syst_plots.py"
        "$config_default"
        "$output_dir"
        --multitrial_type fit
        --max_chi2 "$max_chi2"
        --min_signif "$min_signif"
        --max_signif "$max_signif"
    )

    if [ "$force_prompt_enhanced" = "true" ]; then
        cmd+=(--force_prompt_enhanced)
    fi

    "${cmd[@]}" > "$output_dir/syst/multitrial/fit/summary/log_produce_plots.txt" 2>&1
    echo "Final plots produced!"
fi
