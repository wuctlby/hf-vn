#!/usr/bin/env python3

import argparse
import os
from multiprocessing import Pool
script_dir = os.path.dirname(os.path.realpath(__file__))
os.sys.path.append(os.path.join(script_dir, '../..', 'src'))
from get_vn_vs_mass import get_vn_vs_mass
from ROOT import gSystem

# Load your precompiled library once
script_dir = os.path.dirname(os.path.realpath(__file__))
gSystem.Load(f"{script_dir}/../../invmassfitter/libvnfitter.so")

def fit_proj(args):
    yaml_file, proj_file = args
    print(f"Fitting projection file: {proj_file} with config: {yaml_file}")
    get_vn_vs_mass(yaml_file, proj_file, batch=True, isMultitrial=True)

def main():
    parser = argparse.ArgumentParser(description="Run vn vs mass for multiple YAMLs")
    parser.add_argument('yaml_files', nargs='+', help='List of YAML config files')
    parser.add_argument('--nproc', type=int, default=8, help='Number of parallel processes')
    args = parser.parse_args()

    # Prepare all (yaml, proj_file) pairs
    all_tasks = []
    for yaml_file in args.yaml_files:
        trial_dir = os.path.dirname(yaml_file)
        proj_dir = os.path.join(trial_dir, "projs")
        # check if proj_dir exists
        if not os.path.exists(proj_dir):
            print(f"Projection directory {proj_dir} does not exist. Skipping.")
            continue
        proj_files = sorted([os.path.join(proj_dir, f)
                             for f in os.listdir(proj_dir) if f.startswith("proj_") and f.endswith(".root")])
        all_tasks.extend([(yaml_file, pf) for pf in proj_files])

    # Parallel execution
    with Pool(args.nproc) as pool:
        pool.map(fit_proj, all_tasks)

if __name__ == "__main__":
    main()
