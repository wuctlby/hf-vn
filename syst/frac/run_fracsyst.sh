#!/bin/bash
# Bash script to compute resolution systematics in flow analysis

########################################
# Define paths and configuration
########################################

# Default cut variation path (uncorrected)
export cutvar_default="/home/spolitan/alice/ds_v2_oo/output/181225/cutvar_20_correlated/cutVar"

# Array of systematic variation paths (e.g., for "random" and "even" variations)
cutvar_systvariation_paths=(
  "/home/spolitan/alice/ds_v2_oo/output/181225/cutvar_20_correlated/cutVar_random"
  "/home/spolitan/alice/ds_v2_oo/output/181225/cutvar_20_correlated/cutVar_even"
  "/home/spolitan/alice/ds_v2_oo/output/181225/cutvar_20_correlated/cutVar_odd"
  "/home/spolitan/alice/ds_v2_oo/output/181225/cutvar_20_correlated/cutVar_1in4"
  "/home/spolitan/alice/ds_v2_oo/output/181225/cutvar_20_correlated/cutVar_1in3"
  "/home/spolitan/alice/ds_v2_oo/output/181225/cutvar_20_correlated/cutVar_minus3low"
  "/home/spolitan/alice/ds_v2_oo/output/181225/cutvar_20_correlated/cutVar_minus3high"
  "/home/spolitan/alice/ds_v2_oo/output/181225/cutvar_20_correlated/cutVar_minus3"
)

# Other important paths
export baseOutputDir="/home/spolitan/alice/ds_v2_oo/output/181225/"
export corrPath=$baseOutputDir"cutvar_20_correlated/"
export combinedPath=$baseOutputDir"cutvar_20_combined/"
export output_dir="$baseOutputDir/syst/frac/"
export suffix="test"
export config="/home/spolitan/alice/hf-vn/configs/test_config_ds_small.yml"
export n_parallel=2  # Number of parallel jobs

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

########################################
# Parallel function to run analysis for each configuration
########################################

parallel_func() {

    ########################################
    # Build and execute Python commands
    ########################################

    cutVarSystPath=$1
    # Check if the directory exists
    if [ ! -d "$cutVarSystPath" ]; then
        echo "Error: Directory $cutVarSystPath does not exist."
        return 1
    fi
    echo "Processing systematic variation path: $cutVarSystPath"
    
    
    # Command to compute data driven fraction
    echo "Computing data driven fraction for cut variation at: $cutVarSystPath"
    local cmd=(python3 "./src/data_driven_fraction.py" "$cutVarSystPath/cutVar.root" "$combinedPath/effs" -b -o "$cutVarSystPath")
    echo "Executing: ${cmd[*]}"
    "${cmd[@]}"

    # Command to compute prompt FD v2
    local cmd=(python3 "./src/get_v2_vs_frac.py" "$config" "$combinedPath/raw_yields" "$cutVarSystPath/frac" -b -o "$cutVarSystPath")
    echo "Executing: ${cmd[*]}"
    "${cmd[@]}"

    echo "Finished processing systematic variation path: $cutVarSystPath"
}
export -f parallel_func

########################################
# Run the analysis for each configuration in parallel
########################################

# Change directory to the hf-vn folder (adjust as necessary)
cd ../../ || { echo "Failed to change directory to hf-vn"; exit 1; }

echo "Running parallel jobs with n_parallel: $n_parallel"

# Use GNU Parallel to run the function for each systematic variation path in the array
parallel -j "$n_parallel" parallel_func ::: "${cutvar_systvariation_paths[@]}"

echo "All systematic variations processed."

########################################
# Compute systematics
########################################
echo "Computing systematics for data driven fraction ..."
cd ./syst/frac/ || { echo "Failed to change directory to syst/frac"; exit 1; }

# Create list of v2 results
all_v2_files=()

# Add reference as first entry
all_v2_files+=("$baseOutputDir/v2/v2VsFrac.root")

for cutVarSystPath in "${cutvar_systvariation_paths[@]}"; do
    v2_file="$cutVarSystPath/v2/v2VsFrac.root"
    if [ -f "$v2_file" ]; then
        all_v2_files+=("$v2_file")
    else
        echo "Warning: v2 results file not found at $v2_file"
    fi
done

echo "All v2 files to be processed for systematics:"
for file in "${all_v2_files[@]}"; do
    echo " - $file"
done
# Produce systematic plots 
echo "python3 produce_frac_syst_plots.py ${all_v2_files[@]} -s $suffix -o $output_dir"
python3 produce_frac_syst_plots.py ${all_v2_files[@]} -s $suffix -o $output_dir
echo "Systematic plots produced in $output_dir"