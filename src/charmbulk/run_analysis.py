#!/usr/bin/env python3
"""
run_analysis.py — Orchestrator for the charm bulk workflow.

Adapted for the hf-vn-dev/dev-v0/ repository structure.

Reads a flow configuration YAML and runs the requested nodes:

  1. preprocess              — Pre-process large ROOT files (split by pT, skim)
  2. make_yaml               — Generate cutset YAML files
  3. proj_data / proj_mc     — Projection of data / MC distributions
  4. efficiencies            — Compute efficiencies from MC projections
  5. mass_fit / get_vn_vs_mass   — Fit invariant mass (+ vn) distributions
  6. get_vn_yield_extraction — Yield-extraction-based vn method  (via flarefly)
  7. cut_variation           — Corrected yields via minimisation
  8. data_driven_fraction    — Prompt/FD fractions
  9. get_v2_vs_frac          — Prompt/non-prompt v2 via linear extrapolation

Usage:
    python3 run_analysis.py config.yml [--workers N]
"""
import os
import sys
import yaml
import argparse
import time
import concurrent.futures
import subprocess

# ── Paths relative to the dev-v0 repo root ──────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

# Point to the repo's utility modules
sys.path.insert(0, os.path.join(REPO_ROOT, "utils"))
from utils import logger

# The python executable to use (with ROOT + any dependencies)
PYTHON = os.environ.get("WF_PYTHON", "conda run -n alice python3")

SRC_DIR = os.path.join(REPO_ROOT, "src")


def get_paths():
    """Map node names to their scripts in the dev-v0 repository."""
    return {
        "preprocess": os.path.join(SRC_DIR, "pre_process.py"),
        "make_yaml": os.path.join(SRC_DIR, "make_cutsets_cfgs.py"),
        "proj_thn": os.path.join(SRC_DIR, "proj_thn.py"),
        "mass_fit": os.path.join(SRC_DIR, "mass_fit.py"),
        "efficiencies": os.path.join(SRC_DIR, "compute_efficiencies.py"),
        "get_vn_vs_mass": os.path.join(SRC_DIR, "get_vn_vs_mass.py"),
        "get_vn_yield_extraction": os.path.join(SRC_DIR, "get_vn_by_yield_extraction.py"),
        "cut_variation": os.path.join(SRC_DIR, "cut_variation.py"),
        "data_driven_fraction": os.path.join(SRC_DIR, "data_driven_fraction.py"),
        "get_v2_vs_frac": os.path.join(SRC_DIR, "get_v2_vs_frac.py"),
    }


def run_cmd(cmd, timeout=600):
    """Run a shell command with logging."""
    logger(f"Running: {cmd}", "COMMAND")
    result = subprocess.run(cmd, shell=True, timeout=timeout, text=True, stderr=subprocess.STDOUT)
    if result.stdout:
        for line in result.stdout.strip().split("\n"):
            logger(line, "OUTPUT")
    if result.returncode != 0:
        logger(f"Command failed (rc={result.returncode})", "ERROR")
        if result.stderr:
            for line in result.stderr.strip().split("\n")[-20:]:
                logger(f"  {line}", "ERROR")
        raise RuntimeError(f"Command failed: {cmd}")
    return result


def main():
    parser = argparse.ArgumentParser(description="Charm bulk workflow")
    parser.add_argument("flow_config", help="Flow configuration YAML file")
    parser.add_argument("--workers", "-w", type=int, default=1,
                        help="Number of workers for parallel steps")
    parser.add_argument("--test", "-t", action="store_true",
                        help="Test mode (3 pT bins)")
    args = parser.parse_args()

    start_time = time.time()

    # ── Load config ────────────────────────────────────────────────────
    with open(args.flow_config, "r") as f:
        config = yaml.safe_load(f)

    operations = config.get("operations", {})
    base_outdir = config.get("outdir", ".")
    suffix = config.get("suffix", "default")

    os.makedirs(base_outdir, exist_ok=True)

    # Copy config
    cfg_copy_dir = os.path.join(base_outdir, "config_flow")
    os.makedirs(cfg_copy_dir, exist_ok=True)
    cfg_copy = os.path.join(cfg_copy_dir, os.path.basename(args.flow_config))
    os.system(f"cp {args.flow_config} {cfg_copy}")

    # For testing: modify to use only first 3 pT bins
    if args.test:
        ptbins = config.get("ptbins", [])
        if len(ptbins) > 4:
            logger(f"TEST MODE: reducing from {len(ptbins)-1} to 3 pT bins", "WARNING")
            config["ptbins"] = ptbins[:4]
            test_cfg_path = os.path.join(
                base_outdir, "config_flow", f"test_{os.path.basename(args.flow_config)}")
            with open(test_cfg_path, "w") as f:
                yaml.dump(config, f, default_flow_style=True)
            flow_config = test_cfg_path
            logger(f"Using test config: {test_cfg_path}", "INFO")
        else:
            flow_config = args.flow_config
    else:
        flow_config = args.flow_config

    paths = get_paths()
    n_workers = args.workers

    # ════════════════════════════════════════════════════════════════════
    # Node 0: preprocess
    # ════════════════════════════════════════════════════════════════════
    if operations.get("preprocess", False):
        logger("Step 0: Pre-processing...", "INFO")
        cmd = f"{PYTHON} {paths['preprocess']} {flow_config} -w {n_workers}"
        run_cmd(cmd)
    else:
        logger("Step 0: preprocess — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Node 1: make_yaml
    # ════════════════════════════════════════════════════════════════════
    if operations.get("make_yaml", False):
        logger("Step 1: Generating cutset YAML files...", "INFO")
        outdir = os.path.join(base_outdir, f"cutvar_{suffix}_combined")
        cmd = f"{PYTHON} {paths['make_yaml']} {flow_config} -o {outdir} -c"
        run_cmd(cmd)
    else:
        logger("Step 1: make_yaml — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Node 2: proj_data / proj_mc  (projections over all cutsets)
    # ════════════════════════════════════════════════════════════════════
    if operations.get("proj_data", False) or operations.get("proj_mc", False):
        logger("Step 2: Projecting distributions over all cutsets...", "INFO")

        # Find cutset directory (generated by make_yaml in step 1)
        cutset_dir = os.path.join(base_outdir, f"cutvar_{suffix}_combined", "cutsets")
        if not os.path.isdir(cutset_dir):
            cutset_dir = os.path.join(base_outdir, "cutsets")
        cutset_files = sorted(
            [f for f in os.listdir(cutset_dir) if f.endswith(".yml")]
        ) if os.path.isdir(cutset_dir) else []
        if not cutset_files:
            logger(f"No cutset files found in {cutset_dir}", "ERROR")
            sys.exit(1)

        proj_out = os.path.join(base_outdir, "projs")
        os.makedirs(proj_out, exist_ok=True)
        m_cutsets = len(cutset_files)
        logger(f"Found {m_cutsets} cutsets, projecting each...", "INFO")

        def run_projection(i):
            cutset_config = os.path.join(cutset_dir, f"cutset_{i:02d}.yml")
            cmd = (f"{PYTHON} {paths['proj_thn']} {flow_config} "
                   f"--cutsetConfig {cutset_config} -o {proj_out}")
            return run_cmd(cmd)

        if n_workers > 1 and m_cutsets > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
                list(executor.map(run_projection, range(m_cutsets)))
        else:
            for i in range(m_cutsets):
                run_projection(i)
    else:
        logger("Step 2: projections — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Node 3: efficiencies  (compute from MC projections)
    # ════════════════════════════════════════════════════════════════════
    if operations.get("efficiencies", False):
        logger("Step 3: Computing efficiencies...", "INFO")

        proj_dir = config.get("inputs", {}).get("projs", "") or os.path.join(base_outdir, "projs")
        eff_dir = os.path.join(base_outdir, "effs")
        os.makedirs(eff_dir, exist_ok=True)

        proj_files = sorted(
            f for f in os.listdir(proj_dir)
            if f.startswith("proj_") and f.endswith(".root"))

        def run_eff(proj_file):
            proj_path = os.path.join(proj_dir, proj_file)
            cmd = f"{PYTHON} {paths['efficiencies']} {flow_config} {proj_path} -b"
            return run_cmd(cmd)

        if n_workers > 1 and len(proj_files) > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
                list(executor.map(run_eff, proj_files))
        else:
            for pf in sorted(proj_files):
                run_eff(pf)
    else:
        logger("Step 3: efficiencies — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Node 4: mass_fit  (flarefly-based mass-only fit → raw yields)
    #
    #  Uses mass_fit.py which relies on flarefly's F2MassFitter.
    #  Config expects a "v2extraction" or "fitConfig" section with:
    #    MassFitRanges, SgnFunc, BkgFunc, Rebin, Sigma, etc.
    # ════════════════════════════════════════════════════════════════════
    if operations.get("mass_fit", False):
        logger("Step 4: Fitting invariant mass distributions (flarefly)...", "INFO")

        custom_projs = config.get("inputs", {}).get("projs", "")
        proj_dir = custom_projs if custom_projs else os.path.join(base_outdir, "projs")

        ry_dir = os.path.join(base_outdir, "raw_yields")
        os.makedirs(ry_dir, exist_ok=True)

        proj_files = sorted(
            f for f in os.listdir(proj_dir)
            if f.startswith("proj_") and f.endswith(".root"))

        def run_mass_fit(proj_file):
            proj_path = os.path.join(proj_dir, proj_file)
            cmd = f"CUDA_VISIBLE_DEVICES=-1 {PYTHON} {paths['mass_fit']} {flow_config} {proj_path} -o {ry_dir} -b"
            return run_cmd(cmd)

        if n_workers > 1 and len(proj_files) > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
                list(executor.map(run_mass_fit, proj_files))
        else:
            for pf in sorted(proj_files):
                run_mass_fit(pf)
    else:
        ry_dir = os.path.join(base_outdir, "raw_yields")
        if config.get("inputs", {}).get("raw_yields", ""):
            ry_dir = config["inputs"]["raw_yields"]
        logger("Step 4: mass_fit — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Node 4b: get_vn_vs_mass  (simultaneous mass + vn fit)
    #
    #  Uses get_vn_vs_mass.py which performs a simultaneous mass + vn fit
    #  using VnVsMassFitter.
    #
    #  Config expects a "v2extraction" or "fitConfig" section with:
    #    MassFitRanges, SgnFunc, BkgFunc, BkgFuncVn, Rebin, etc.
    # ════════════════════════════════════════════════════════════════════
    if operations.get("get_vn_vs_mass", False):
        logger("Step 4b: Fitting mass + vn simultaneously...", "INFO")

        custom_projs = config.get("inputs", {}).get("projs", "")
        proj_dir = custom_projs if custom_projs else os.path.join(base_outdir, "projs")

        ry_dir = os.path.join(base_outdir, "raw_yields")
        os.makedirs(ry_dir, exist_ok=True)

        proj_files = sorted(
            f for f in os.listdir(proj_dir)
            if f.startswith("proj_") and f.endswith(".root"))

        def run_vn_vs_mass(proj_file):
            proj_path = os.path.join(proj_dir, proj_file)
            cmd = f"{PYTHON} {paths['get_vn_vs_mass']} {flow_config} {proj_path} -b"
            return run_cmd(cmd)

        if n_workers > 1 and len(proj_files) > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
                list(executor.map(run_vn_vs_mass, proj_files))
        else:
            for pf in sorted(proj_files):
                run_vn_vs_mass(pf)
    else:
        logger("Step 4b: get_vn_vs_mass — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Node 5: get_vn_yield_extraction  (yield-extraction-based vn method)
    # ════════════════════════════════════════════════════════════════════
    if operations.get("get_vn_yield_extraction", False):
        logger("Step 5: Extracting vn via yield extraction (flarefly)...", "INFO")

        custom_projs = config.get("inputs", {}).get("projs", "")
        proj_dir = custom_projs if custom_projs else os.path.join(base_outdir, "projs")

        cmd = f"{PYTHON} {paths['get_vn_yield_extraction']} {flow_config} -b"
        run_cmd(cmd)
    else:
        logger("Step 5: get_vn_yield_extraction — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Node 6: cut_variation
    # ════════════════════════════════════════════════════════════════════
    if operations.get("cut_variation", False):
        logger("Step 6: Cut variation...", "INFO")

        custom_ry = config.get("inputs", {}).get("raw_yields", "")
        ry_path = custom_ry if custom_ry else ry_dir

        custom_eff = config.get("inputs", {}).get("effs", "")
        eff_path = custom_eff if custom_eff else os.path.join(base_outdir, "effs")

        cmd = f"{PYTHON} {paths['cut_variation']} {flow_config} {ry_path} {eff_path} -b"
        run_cmd(cmd)
    else:
        logger("Step 6: cut_variation — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Node 7: data_driven_fraction
    # ════════════════════════════════════════════════════════════════════
    if operations.get("data_driven_fraction", False):
        logger("Step 7: Data-driven fraction...", "INFO")

        custom_eff = config.get("inputs", {}).get("effs", "")
        eff_path = custom_eff if custom_eff else os.path.join(base_outdir, "effs")

        cutvar_path = os.path.join(base_outdir, "cutVar", "cutVar.root")
        if not os.path.exists(cutvar_path):
            logger(f"cutVar.root not found at {cutvar_path}, using default", "WARNING")

        cmd = f"{PYTHON} {paths['data_driven_fraction']} {cutvar_path} {eff_path} -o {base_outdir} -b"
        run_cmd(cmd)
    else:
        logger("Step 7: data_driven_fraction — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Node 8: get_v2_vs_frac  (prompt/non-prompt v2 via linear extrapolation)
    # ════════════════════════════════════════════════════════════════════
    if operations.get("get_v2_vs_frac", False):
        logger("Step 8: v2 vs fd fraction extrapolation...", "INFO")

        frac_dir = os.path.join(base_outdir, "frac")
        v2_dir = os.path.join(base_outdir, "v2")
        os.makedirs(v2_dir, exist_ok=True)

        cmd = f"{PYTHON} {paths['get_v2_vs_frac']} {flow_config} {ry_dir} {frac_dir} -b"
        run_cmd(cmd)
    else:
        logger("Step 8: get_v2_vs_frac — SKIPPED", "WARNING")

    elapsed = time.time() - start_time
    logger(f"Analysis completed in {elapsed:.1f} seconds", "INFO")


if __name__ == "__main__":
    main()
