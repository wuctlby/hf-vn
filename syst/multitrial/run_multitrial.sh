#!/bin/bash

generate_cutset() {
    local yaml_file="$1"
    trial_num=$(basename "$yaml_file" | sed -E 's/config_trial_([0-9]+)\.yml/\1/')
    target_dir="$dir/trial_$trial_num"
    echo "[$trial_num] Generating YAML cutsets for file $yaml_file into $target_dir"
    python3 /home/mdicosta/alice/hf-vn/src/make_cutsets_cfgs.py "$yaml_file" -o "$target_dir"
}

get_v2_vs_fracs() {
    local yaml_file="$1"
    local reference_frac_dir="$2"
    trial_num=$(basename "$yaml_file" | sed -E 's/config_trial_([0-9]+)\.yml/\1/')
    target_dir="$dir/trial_$trial_num"
    echo "[$trial_num] Computing v2 vs frac for file $yaml_file"
    python3 /home/mdicosta/alice/hf-vn/src/get_v2_vs_frac.py "$yaml_file" "$target_dir/raw_yields/" "$reference_frac_dir" --multitrial
}

# Export function and necessary variables
export -f generate_cutset
export -f get_v2_vs_fracs
export dir  # make sure $dir is visible inside the function

export n_parallel=20
export config_modifies="/home/mdicosta/FlowDplus/TestRCT/systMultitrial.yml"
export config_default="/home/mdicosta/FlowDplus/TestRCT/config.yml"
export output_dir="/home/mdicosta/FlowDplus/TestRCT/cutvar_RCTFlagsOn_combined"
export do_cutset_generation=false
export do_compile_fitter=false
export do_projections=false
export do_fitting=false
export do_v2_vs_frac=false
export produce_plots=true

mkdir -p "$output_dir/syst/multitrial"


if [ "$do_cutset_generation" = true ]; then
    # Generate YAML file list with indices
    echo "Generating yaml files for systematics ..."
    python3 /home/mdicosta/alice/hf-vn/syst/multitrial/make_cutsets_multitrial.py $config_default -m $config_modifies -o $output_dir > "$output_dir/syst/multitrial/log_make_cutsets.txt" 2>&1
    echo "Yaml files generated!"
fi

if [ "$do_compile_fitter" = true ]; then
    echo "Compiling InvMassFitter and VnVsMassFitter ..."
    rootcling -f ./invmassfitter/vnfitter_dict.cxx -c ./invmassfitter/InvMassFitter.h ./invmassfitter/VnVsMassFitter.h ./syst/multitrial/LinkDefVnFitter.h

    g++ -shared -fPIC `root-config --cflags --libs` \
        ./invmassfitter/InvMassFitter.cxx ./invmassfitter/VnVsMassFitter.cxx ./invmassfitter/vnfitter_dict.cxx \
        -o ./invmassfitter/libvnfitter.so
    echo "Compilation done!"
fi

pt_dirs=($(ls -d "$output_dir"/syst/multitrial/pt_*))
for dir in "${pt_dirs[@]}"; do
    log_file="$dir/log.txt"
    echo -e "\nProcessing directory $dir"

    # Record start time for directory
    dir_start=$(date +%s)

    # Find YAML files in the current directory
    yaml_files=($(find "$dir" -maxdepth 2 -type f -name "config_trial_*.yml"))
    echo "Number of yaml files found: ${#yaml_files[@]}"

    # --- Generate YAML cutsets ---
    if [ "$do_cutset_generation" = true ]; then
        start=$(date +%s)
        parallel -j $n_parallel generate_cutset ::: "${yaml_files[@]}" > "$log_file" 2>&1
        end=$(date +%s)
        echo -e "----> Cutsets generated in $((end - start)) seconds. Starting projections ..."
    fi

    # --- Perform projections ---
    if [ "$do_projections" = true ]; then
        start=$(date +%s)
        python3 /home/mdicosta/alice/hf-vn/src/proj_thn.py $config_default --multitrial_folder "$dir" > "$log_file" 2>&1
        end=$(date +%s)
        echo -e "----> Projections performed in $((end - start)) seconds. Starting fitting ..."
    fi

    # --- Perform fitting ---
    if [ "$do_fitting" = true ]; then
        start=$(date +%s)
        python3 /home/mdicosta/alice/hf-vn/syst/multitrial/run_fits.py "${yaml_files[@]}" --nproc $n_parallel > "$log_file" 2>&1
        end=$(date +%s)
        echo -e "----> Fitting done in $((end - start)) seconds. Starting v2 vs frac extraction ..."
    fi

    # --- Obtain v2 vs frac ---
    if [ "$do_v2_vs_frac" = true ]; then
        reference_frac_dir="$output_dir/frac/"
        start=$(date +%s)
        parallel -j $n_parallel get_v2_vs_fracs ::: "${yaml_files[@]}" ::: "$reference_frac_dir" > "$log_file" 2>&1
        end=$(date +%s)
        echo -e "----> V2 vs frac extracted in $((end - start)) seconds."
    fi

    # Total time for this pt_dir
    dir_end=$(date +%s)
    echo -e "✅ Total time for $dir: $((dir_end - dir_start)) seconds.\n"
    echo "-----------------------------------------------"

done

# --- Produce final plots ---
if [ "$produce_plots" = true ]; then
    echo "Producing final systematic plots ..."
    python3 /home/mdicosta/alice/hf-vn/syst/multitrial/produce_multitrial_syst_plots.py "$output_dir" > "$output_dir/syst/multitrial/log_produce_plots.txt" 2>&1
    echo "Final plots produced!"
fi
