#!/usr/bin/env python3
"""
run_analysis.py — Orchestrator for the charm bulk workflow (flook charmbulk).

Adapted from the top-level run_analysis.py, but:
  - Only correlated cuts (no uncorrelated branch) — fully removed.
  - Uses proj_data_charmbulk projections from proj_XX.root.
  - Splits mass fit, efficiencies, and cut variation by (side, MeanPtBin).

Expected workflow (controlled by YAML operations):
    preprocess → make_yaml → proj_data → proj_mc → efficiencies → mass_fit → do_cut_variation

Output directory structure:
    projs/                                                   — projection files (proj_XX.root)
    effs/{side}/{meanpt_label}/eff_XX.root                   — efficiency per (side, meanPt, cutset)
    raw_yields/{side}/{meanpt_label}/raw_yields_XX.root      — raw yields per (side, meanPt, cutset)
    cutVar/{side}/{meanpt_label}/cutVar/cutVar.root          — cut variation results per (side, meanPt)

Usage:
    python3 run_analysis.py config_flow.yml [--workers N] [--test]
"""
import os
import sys
import yaml
import argparse
import time
import concurrent.futures
import subprocess
import shutil

# ── Paths relative to the dev-v0 repo root ──────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
from charmbulk_utils import get_paths
SCRIPTS = get_paths()
print(f"{SCRIPTS['Utils']}")
sys.path.append(SCRIPTS["Utils"])
from utils import logger

# Use current python3 if inside alice env, otherwise plain python3
if os.environ.get("CONDA_DEFAULT_ENV") == "alice":
    # PYTHON = "conda run -n alice python3"
    PYTHON = "python3"
else:
    PYTHON = "python3"

def run_cmd(cmd, timeout=600):
    """Run a shell command with logging."""
    logger(f"{cmd}", "COMMAND")
    result = subprocess.run(cmd, shell=True, text=True, stderr=subprocess.STDOUT)

    if result.stdout:
        for line in result.stdout.strip().split("\n"):
            logger(line, "OUTPUT")

    if result.returncode != 0:
        logger(f"Command failed (rc={result.returncode})", "ERROR")
        if result.stdout:
            logger("--- Last 25 lines of output before failure ---", "ERROR")
            for line in result.stdout.strip().split("\n")[-25:]:
                logger(f"  {line}", "ERROR")
        raise RuntimeError(f"Command failed: {cmd}")

    return result

def get_cutset_count(cutset_dir):
    """Return number of cutset YAML files in *cutset_dir*."""
    if not os.path.isdir(cutset_dir):
        return 0
    return len([f for f in os.listdir(cutset_dir)
                if f.endswith(".yml") and f.startswith("cutset_")])

def discover_mean_pt_labels(proj_dir):
    """
    Discover MeanPtBin labels from the first projection file available.
    Returns sorted list like ['pt_0_71', 'pt_71_74', ..., 'pt_79_97'].
    Returns empty list if unable to discover.
    """
    if not os.path.isdir(proj_dir):
        return []
    proj_files = sorted(f for f in os.listdir(proj_dir)
                        if f.startswith("proj_") and f.endswith(".root"))
    if not proj_files:
        return []
    try:
        import ROOT
        f = ROOT.TFile.Open(os.path.join(proj_dir, proj_files[0]))
        pt_labels = sorted(k.GetName() for k in f.GetListOfKeys()
                          if k.GetClassName() == "TDirectoryFile" and
                          k.GetName().startswith("pt_"))
        if not pt_labels:
            f.Close()
            return []
        pt_dir = f.Get(pt_labels[0])
        if not pt_dir:
            f.Close()
            return []
        labels = sorted(k.GetName() for k in pt_dir.GetListOfKeys()
                       if k.GetClassName() == "TDirectoryFile")
        f.Close()
        pt_labels = [l for l in labels if l.startswith("pt_")]
        return pt_labels
    except Exception:
        return []


def get_mean_pt_labels_from_config(config, proj_dir):
    """
    Get MeanPtBin labels from config MeanPtBins, or discover from projection file.
    """
    mean_pt_bins = config.get("projections", {}).get("proj_data", {}).get("MeanPtBins", [])
    if mean_pt_bins and len(mean_pt_bins) >= 2:
        labels = []
        for lo, hi in zip(mean_pt_bins[:-1], mean_pt_bins[1:]):
            labels.append(f"pt_{int(lo*100)}_{int(hi*100)}")
        return labels
    return discover_mean_pt_labels(proj_dir)


# ═══════════════════════════════════════════════════════════════════════
# Steps of the workflow
# ═══════════════════════════════════════════════════════════════════════
def step_preprocess(config_path, n_workers):
    """Step 0: Pre-process large ROOT files."""
    logger("Preprocess will be performed", level="INFO")
    cmd = f"{PYTHON} {SCRIPTS['Preprocess']} {config_path} -w {n_workers}"
    run_cmd(cmd)


def step_make_yaml(config_path, outdir):
    """Step 1: Generate cutset YAML files (correlated only)."""
    cmd = f"{PYTHON} {SCRIPTS['YamlCuts']} {config_path} -o {outdir} -c"
    run_cmd(cmd)


def step_projections(config_path, outdir, n_workers, mCutsets):
    """Step 2: Project distributions for all cutsets."""
    proj_out = os.path.join(outdir, "projs")
    os.makedirs(proj_out, exist_ok=True)
    cutset_dir = os.path.join(outdir, "cutsets")

    def run_projections(i):
        iCutSets = f"{i:02d}"
        logger(f"Performing projections for cutset {iCutSets}...", level="INFO")
        cutset_config = os.path.join(cutset_dir, f"cutset_{i:02d}.yml")
        cmd = (f"{PYTHON} {SCRIPTS['Projections']} {config_path} "
               f"--cutsetConfig {cutset_config} -o {proj_out}")
        return run_cmd(cmd)

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        list(executor.map(run_projections, range(mCutsets)))

def step_efficiencies(config_path, cutvar_dir, n_workers, m_cutsets):
    """Step 3: Compute efficiencies for each cutset."""
    eff_out = os.path.join(cutvar_dir, "effs")
    os.makedirs(eff_out, exist_ok=True)
    cutset_dir = os.path.join(cutvar_dir, "cutsets")

    def run_efficiency(i):
        iCutSet = f"{i:02d}"
        logger(f"Computing efficiencies for cutset {iCutSet}...", level="INFO")
        proj_cutset = os.path.join(cutvar_dir, "projs", f"proj_{iCutSet}.root")
        cmd = f"{PYTHON} {SCRIPTS['Efficiencies']} {config_path} {proj_cutset} -b --mode charm_bulk"
        return run_cmd(cmd)
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        list(executor.map(run_efficiency, range(m_cutsets)))

def step_mass_fit(config_path, cutvar_dir, n_workers, m_cutsets):
    """Step 4: Fit invariant mass distributions for each cutset."""
    ry_out = os.path.join(cutvar_dir, "raw_yields")
    os.makedirs(ry_out, exist_ok=True)
    cutset_dir = os.path.join(cutvar_dir, "cutsets")

    # fix sigma to the sigma of the first cutset to stabilize fits across cutsets

    def run_mass_fit(i):
        iCutSet = f"{i:02d}"
        logger(f"Fitting mass distributions for cutset {iCutSet}...", level="INFO")
        proj_cutset = os.path.join(cutvar_dir, "projs", f"proj_{iCutSet}.root")
        if i == 0:
            load_sigma=""
        else:
            load_sigma="--load-sigma"
        cmd = f"{PYTHON} {SCRIPTS['MassFit']} {config_path} {proj_cutset} -b {load_sigma}"
        return run_cmd(cmd)

    run_mass_fit(0)  # Run first cutset to get sigma values

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        list(executor.map(run_mass_fit, range(1, m_cutsets)))


def step_cut_variation(config_path, cutvar_dir, n_workers):
    """
    Step 5: Cut variation for each (side, MeanPtBin) combination.

    cut_variation.py uses str.replace('eff', 'cutset') on the full eff file path
    to locate cutset YAMLs, and saves cutVar/ at os.path.dirname(eff_path)/cutVar/.

    To keep per-MeanPtBin results separate, we create isolated wrapper directories
    with symlinks to the actual eff and cutset files. After cut_variation runs,
    cutVar/cutVar.root appears at the expected per-MeanPtBin location.
    """
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    proj_dir = os.path.join(cutvar_dir, "projs")
    mean_pt_labels = get_mean_pt_labels_from_config(config, proj_dir)
    if not mean_pt_labels:
        logger("No MeanPtBin labels — cannot run cut variation", "FATAL")

    ry_base = os.path.join(cutvar_dir, "raw_yields")
    eff_base = os.path.join(cutvar_dir, "effs")
    cutset_base = os.path.join(cutvar_dir, "cutsets")
    cutvar_base = os.path.join(cutvar_dir, "cutVar")
    sides = ["a_side", "b_side"]
    paths = get_paths()

    def run_one(side_meanpt):
        side, mpl = side_meanpt
        ry_path = os.path.join(ry_base, side, mpl)
        eff_path_actual = os.path.join(eff_base, side, mpl)

        if not os.path.isdir(ry_path):
            logger(f"  Raw yields dir not found: {ry_path}, skipping", "WARNING")
            return
        if not os.path.isdir(eff_path_actual):
            logger(f"  Eff dir not found: {eff_path_actual}, skipping", "WARNING")
            return

        # # Create isolated wrapper dir
        # wrapper_dir = os.path.join(cutvar_base, side, mpl, "_inputs")
        # os.makedirs(wrapper_dir, exist_ok=True)

        # # Symlink eff files
        # for fname in os.listdir(eff_path_actual):
        #     if fname.endswith(".root"):
        #         src = os.path.relpath(os.path.join(eff_path_actual, fname), wrapper_dir)
        #         dst = os.path.join(wrapper_dir, fname)
        #         if not os.path.exists(dst):
        #             os.symlink(src, dst)

        # # Symlink cutset files (named to match str.replace('eff', 'cutset'))
        # for fname in os.listdir(cutset_base):
        #     if fname.startswith("cutset_") and fname.endswith(".yml"):
        #         src = os.path.relpath(os.path.join(cutset_base, fname), wrapper_dir)
        #         # The str.replace('eff','cutset') maps eff_00.root → cutset_00.yml
        #         # Our cutset files are already named cutset_00.yml
        #         dst = os.path.join(wrapper_dir, fname)
        #         if not os.path.exists(dst):
        #             os.symlink(src, dst)

        cmd = f"{PYTHON} {paths['CutVariation']} {config_path} {ry_path} {eff_path_actual} -b"
        try:
            run_cmd(cmd)
            cv_result = os.path.join(cutvar_base, side, mpl, "cutVar", "cutVar.root")
            if os.path.exists(cv_result):
                logger(f"  Cut variation done for {side}/{mpl}", "INFO")
            else:
                logger(f"  cutVar.root not found at {cv_result}", "WARNING")
        except RuntimeError as e:
            logger(f"  Cut variation failed for {side}/{mpl}: {e}", "ERROR")

        # # Clean up
        # shutil.rmtree(wrapper_dir, ignore_errors=True)

    combinations = [(s, mpl) for s in sides for mpl in mean_pt_labels]

    if n_workers > 1 and len(combinations) > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
            list(executor.map(run_one, combinations))
    else:
        for s, mpl in combinations:
            run_one((s, mpl))


# ═══════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Charm bulk workflow")
    parser.add_argument("config_path", help="Flow configuration YAML file")
    parser.add_argument("--workers", "-w", type=int, default=1,
                        help="Number of workers for parallel steps")
    parser.add_argument("--test", "-t", action="store_true",
                        help="Test mode (3 pT bins)")
    args = parser.parse_args()

    start_time = time.time()

    # ── Load config ────────────────────────────────────────────────────
    with open(args.config_path, "r") as f:
        config = yaml.safe_load(f)

    operations = config.get("operations", {})
    base_outdir = config.get("outdir", ".")
    suffix = config.get("suffix", "")
    n_workers = args.workers
    if args.test:
        suffix += "_test"
        logger("TEST MODE ENABLED: using only first 3 pT bins and modified output dirs", "WARNING")
    os.makedirs(base_outdir, exist_ok=True)

    # Copy config
    nfile = 0
    cfg_copy_dir = os.path.join(base_outdir, f"cutvar_{suffix}", "config")
    os.makedirs(cfg_copy_dir, exist_ok=True)
    while os.path.exists(f'{cfg_copy_dir}/{os.path.splitext(os.path.basename(args.config_path))[0]}_{config["suffix"]}_{nfile}.yml'):
        nfile = nfile + 1
    os.system(f'cp {args.config_path} {cfg_copy_dir}/{os.path.splitext(os.path.basename(args.config_path))[0]}_{config["suffix"]}_{nfile}.yml')

    # For testing: modify to use only first 3 pT bins
    if args.test:
        ptbins = config.get("ptbins", [])
        if len(ptbins) > 4:
            logger(f"TEST MODE: reducing from {len(ptbins)-1} to 3 pT bins", "WARNING")
            config["ptbins"] = ptbins[:4]
            test_cfg_path = os.path.join(cfg_copy_dir, f"cutvar_{suffix}_config.yml")
            with open(test_cfg_path, "w") as f:
                yaml.dump(config, f, default_flow_style=True)
            config_path = test_cfg_path
            logger(f"TEST MODE: using test config: {test_cfg_path}", "INFO")
        else:
            config_path = args.config_path
    else:
        config_path = args.config_path

    # ── Build output directory name ─────────────────────────────────────
    outdir = os.path.join(base_outdir, f"cutvar_{suffix}")
    os.makedirs(outdir, exist_ok=True)

    # ════════════════════════════════════════════════════════════════════
    # Step 0: preprocess
    # ════════════════════════════════════════════════════════════════════
    if operations.get("preprocess", False):
        logger("Step 0: Pre-processing...", "INFO")
        step_preprocess(config_path, n_workers)
    else:
        logger("Step 0: Pre-processing — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Step 1: make_yaml (correlated only)
    # ════════════════════════════════════════════════════════════════════
    if operations.get("make_yaml", False):
        logger("Step 1: Generating cutset YAML files (correlated cuts)...", "INFO")
        step_make_yaml(config_path, outdir)
    else:
        logger("Step 1: Generating cutset YAML files — SKIPPED", "WARNING")

    # ── Determine number of cutsets ────────────────────────────────────
    mCutSets = len([f for f in os.listdir(f"{outdir}/cutsets") if os.path.isfile(os.path.join(f"{outdir}/cutsets", f))])
    logger(f"Number of cutsets: {mCutSets}", "INFO")

    # ════════════════════════════════════════════════════════════════════
    # Step 2: proj_data / proj_mc
    # ═══════════════════════════════════════════════════════════════════════
    if operations.get("proj_data", False) or operations.get("proj_mc", False):
        logger("Step 2: Projecting distributions over all cutsets...", "INFO")
        step_projections(config_path, outdir, n_workers, mCutSets)
    else:
        logger("Step 2: projections — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Step 3: efficiencies (a_side / b_side separated)
    # ════════════════════════════════════════════════════════════════════
    if operations.get("efficiencies", False):
        logger("Step 3: Computing efficiencies (a_side / b_side)...", "INFO")
        
        step_efficiencies(config_path, outdir, n_workers, mCutSets)
    else:
        logger("Step 3: efficiencies — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Step 4: mass_fit (flarefly-based mass-only fit → raw yields)
    # ════════════════════════════════════════════════════════════════════
    if operations.get("mass_fit", False):
        logger("Step 4: Fitting invariant mass distributions...", "INFO")
        step_mass_fit(config_path, outdir, n_workers, mCutSets)
    else:
        logger("Step 4: mass_fit — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════════════════════════
    # Step 5: cut_variation
    # ════════════════════════════════════════════════════════════════════
    if operations.get("cut_variation", False):
        logger("Step 5: Cut variation for each (side, mean pT)...", "INFO")
        step_cut_variation(config_path, outdir, n_workers)
    else:
        logger("Step 5: cut_variation — SKIPPED", "WARNING")

    elapsed = time.time() - start_time
    logger(f"Analysis completed in {elapsed:.1f} seconds", "INFO")


if __name__ == "__main__":
    main()
