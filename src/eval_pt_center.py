import argparse
import sys
import os
import re
import array
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # pylint: disable=wrong-import-position
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
import itertools
sys.path.append("./flareflyfitter/")
from raw_yield_fitter import RawYieldFitter
import yaml
import ROOT
from ROOT import TFile
import awkward as ak
sys.path.append(f"{os.path.dirname(os.path.abspath(__file__))}/../utils")
from utils import logger
from matplotlib import gridspec
import uproot
from multiprocessing import Pool, cpu_count
ROOT.gROOT.SetBatch(True)

def fit_control_var(df, i_bin, pt_min, pt_max, cfg_full, cfg_fit, output_dir, sel_string):
    print("\n\n")
    part_name = cfg_full['Dmeson']
    fitter = RawYieldFitter(part_name, os.path.basename(output_dir), True)
    fitter.set_rebin(cfg_fit.get('rebin', 1))
    fitter.set_fit_range(cfg_fit['MassFitRanges'][i_bin][0], cfg_fit['MassFitRanges'][i_bin][1])
    fitter.set_data_to_fit_df(df)

    bkg_funcs = cfg_fit['BkgFunc'][i_bin] if isinstance(cfg_fit['BkgFunc'], list) else [cfg_fit['BkgFunc']]
    for bkg_func in bkg_funcs:
        fitter.add_bkg_func(bkg_func, "Comb. bkg")

    sgn_funcs = cfg_fit['SgnFunc'][i_bin] if isinstance(cfg_fit['SgnFunc'], list) else [cfg_fit['SgnFunc']]
    for i_sgn, sgn_func in enumerate(sgn_funcs):
        fitter.add_sgn_func(sgn_func, f"signal_{i_sgn}")

    # Add correlated background if specified
    if cfg_full.get('corr_bkgs'):
        fitter.add_corr_bkgs(cfg_full['corr_bkgs'], sel_string, pt_min, pt_max)

    fitter.setup()
    fit_res = fitter.fit()

    sgn_sweights = fitter.get_sweights_sgn(len(sgn_funcs)-1)

    os.makedirs(output_dir, exist_ok=True)
    with open(f"{output_dir}/../fits_status.txt", "a") as f:
        computed_sweights = True if sgn_sweights is not None else False
        f.write(
                f"{output_dir}: fit_res.valid -> {fit_res.valid}, "
                f"fit_res.status -> {fit_res.status}, "
                f"fit_res.converged -> {fit_res.converged}, "
                f"sweights computed -> {computed_sweights} \n"
                )

    fig, _ = fitter.plot_fit(False, True, loc=["lower left", "upper left"]) # (log, show_extra_info)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    fig.savefig(
        os.path.join(
            output_dir,
            f'fits/fM_fit.png'
        ),
        dpi=300, bbox_inches="tight"
    )

    return sgn_sweights

def load_aod_file(i_file, aod_file):
    logger(f"Loading AOD file [{i_file}]: {aod_file}", "INFO")
    branches = ["fPt", "fM", "fMlScore0", "fMlScore1", "fScalarProd", "fCent"]
    try:
        f = uproot.open(aod_file)
        o2_keys = [k for k in f.keys() if "O2" in k]
        dfs = []
        for key in o2_keys:
            arr = f[key].arrays(branches, library="ak")
            arr = ak.zip({b: arr[b] for b in arr.fields})  # row-wise structure
            df_tree = pd.DataFrame({b: arr[b].to_numpy() for b in arr.fields})
            dfs.append(df_tree)

        # concatenate all O2 trees in this file
        if dfs:
            df_file = pd.concat(dfs, ignore_index=True)
        else:
            df_file = pd.DataFrame()

        logger(f"[{i_file}] Loaded {aod_file} with {len(df_file)} total entries from {len(o2_keys)} O2 trees", "INFO")
        return df_file
    except Exception as e:
        logger(f"Failed to load {aod_file}: {e}", "ERROR")
        return pd.DataFrame()  # return empty DataFrame on failure

def load_input_df(input_aod_cfg, workers):

    # Input files
    if isinstance(input_aod_cfg, str) and input_aod_cfg.endswith('.root'):
        logger("Single AOD file detected.", "INFO")
        input_aods = [input_aod_cfg]
    elif isinstance(input_aod_cfg, list):
        logger("List of AOD files detected.", "INFO")
        input_aods = input_aod_cfg
    elif isinstance(input_aod_cfg, str) and os.path.isdir(input_aod_cfg):
        logger("Directory of AOD files detected.", "INFO")
        input_aods = [os.path.join(input_aod_cfg, f) for f in os.listdir(input_aod_cfg) if f.endswith('.root') and 'AO2D' in f]
    else:
        logger("Invalid input_aod configuration.", "ERROR")
        sys.exit(1)
    input_aods = sorted(input_aods, key=lambda x: int(__import__('re').search(r'AO2D_(\d+)', x).group(1)))
    print(f"Total AOD files to load: {len(input_aods)}")
    with Pool(processes=workers) as pool:
        dfs = pool.starmap(load_aod_file, enumerate(input_aods))

    # Concatenate
    df = pd.concat(dfs, ignore_index=True)
    logger(f"Total entries: {len(df)}", "INFO")
    return df

def eval_pt_center(cfg_file_name, workers=1):
    # Read the configuration file
    with open(cfg_file_name, 'r') as cfg_file:
        cfg = yaml.safe_load(cfg_file)

    # Retrieve cutsets configs
    cutsets_dir = os.path.join(cfg['outdir'], f"cutvar_{cfg['suffix']}_combined/cutsets")
    cutset_files = [os.path.join(cutsets_dir, f) for f in os.listdir(cutsets_dir) if f.endswith('.yml')]
    cutset_files.sort(key=lambda x: int(re.search(r'(\d+)', os.path.basename(x)).group(1)))
    print(f"Found {cutset_files} cutset files in {cutsets_dir}")

    df = load_input_df(cfg["pt_centering"]["input_aod"], workers=workers)
    infer_vars = df.columns.tolist()
    infer_vars.remove('fM')
    infer_vars.remove('fMlScore0')
    infer_vars.remove('fMlScore1')

    # Loop over cutset configs
    s_weights = {}
    cfg_fit = cfg["v2extraction"]
    for cutset_file in cutset_files:
        with open(cutset_file, 'r') as cs_file:
            cutset_cfg = yaml.safe_load(cs_file)
        out_file_path = cutset_file.replace('cutset', 'ptcenter').replace('.yml', '.root')
        os.makedirs(os.path.dirname(out_file_path), exist_ok=True)
        out_file = TFile.Open(out_file_path, "recreate")

        histos_avgs = {}
        for var in infer_vars:
            histos_avgs[f"h_{var}_sgn"] = ROOT.TH1F(f"h_{var}_sgn", f"h_{var}_sgn", len(cutset_cfg["Pt"]["min"]), array.array('d', cutset_cfg["Pt"]["min"] + [cutset_cfg["Pt"]["max"][-1]]))
            histos_avgs[f"h_{var}_bkg"] = ROOT.TH1F(f"h_{var}_bkg", f"h_{var}_bkg", len(cutset_cfg["Pt"]["min"]), array.array('d', cutset_cfg["Pt"]["min"] + [cutset_cfg["Pt"]["max"][-1]]))

        # Loop over pt bins
        for i_bin, (pt_min, pt_max, score_bkg_min, score_bkg_max, score_fd_min, score_fd_max) in \
            enumerate(zip(cutset_cfg["Pt"]["min"], cutset_cfg["Pt"]["max"],
                          cutset_cfg["ScoreBkg"]["min"], cutset_cfg["ScoreBkg"]["max"],
                          cutset_cfg["ScoreFD"]["min"], cutset_cfg["ScoreFD"]["max"])):

            out_dir_pt = f"{out_file_path.replace('.root', '')}/pt_{int(pt_min*10)}_{int(pt_max*10)}"
            logger(f"Processing pt bin: {pt_min} - {pt_max} of cutset file: {cutset_file}, output directory: {out_dir_pt}", level="INFO")

            os.makedirs(out_dir_pt, exist_ok=True)
            os.makedirs(out_dir_pt + "/fits", exist_ok=True)
            os.makedirs(out_dir_pt + "/vars", exist_ok=True)

            # Query the dataframe
            mass_min, mass_max = cfg_fit["MassFitRanges"][i_bin]
            sel_string = f"{pt_min} <= fPt < {pt_max} and " \
                         f"{score_bkg_min} <= fMlScore0 < {score_bkg_max} and " \
                         f"{score_fd_min} <= fMlScore1 < {score_fd_max} and " \
                         f"{mass_min} <= fM <= {mass_max}"
            sel_string_corr_bkgs = f"fMlScore0 < {score_bkg_max} && " \
                                   f"fMlScore0 >= {score_bkg_min} && " \
                                   f"fMlScore1 < {score_fd_max} && " \
                                   f"fMlScore1 >= {score_fd_min} && " \
                                   f"fM >= {mass_min} && fM < {mass_max}"
            sel_df = df.query(sel_string).reset_index(drop=True)
            fig, ax = plt.subplots(1, 1, figsize=(12, 10))
            sel_df['fM'].hist(bins=100, alpha=0.5, range=(mass_min, mass_max))
            ax.set_xlabel('fM')
            ax.set_ylabel('Counts')
            if not os.path.exists(out_dir_pt):
                os.makedirs(out_dir_pt)
            fig.savefig(
                os.path.join(
                        out_dir_pt,
                        f'fits/fM_raw.png'
                ),
                dpi=300, bbox_inches="tight"
            )

            s_weights = fit_control_var(sel_df, i_bin, pt_min, pt_max, cfg, cfg_fit, out_dir_pt, sel_string_corr_bkgs)
            for var in infer_vars:
                print(f"    Drawing {var}")

                # Create figure with two subplots (distros and ratio)
                fig, ax = plt.subplots(figsize=(8, 8))

                bins = 200
                var_range = (min(sel_df[var]), max(sel_df[var]))
                sgn_vals, bin_edges = np.histogram(sel_df[var], bins=bins, range=var_range, weights=s_weights, density=True)
                bkg_vals, _ = np.histogram(sel_df[var], bins=bins, range=var_range, weights=(1 - s_weights), density=True)
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

                # Plot signal and background
                ax.hist(bin_centers, bins=bins, weights=sgn_vals, label="Signal", color="#1f77b4", alpha=0.5, histtype='step', log=True)
                ax.hist(bin_centers, bins=bins, weights=bkg_vals, label="Bkg", color="#ff7f0e", alpha=0.5, histtype='step', log=True)

                ax.set_ylabel("Entries")
                ax.set_xlabel(var)
                ax.set_title(var)
                ax.legend()

                # Save figure
                fig.tight_layout()
                fig.savefig(
                    os.path.join(out_dir_pt, f'vars/{var}.png'),
                    dpi=300, bbox_inches="tight"
                )
                plt.close(fig)
                
                # Fill histograms for averages
                histos_avgs[f"h_{var}_sgn"].SetBinContent(i_bin + 1, np.average(sgn_vals))
                histos_avgs[f"h_{var}_sgn"].SetBinError(i_bin + 1, np.std(sgn_vals))
                histos_avgs[f"h_{var}_bkg"].SetBinContent(i_bin + 1, np.average(bkg_vals))
                histos_avgs[f"h_{var}_bkg"].SetBinError(i_bin + 1, np.std(bkg_vals))
        
        # Write histograms to output ROOT file
        out_file.cd()
        for hist in histos_avgs.values():
            hist.Write()
        out_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate pt-centering with sPlot and FlareFly')
    parser.add_argument('config_file', help='Path to the input configuration file')
    parser.add_argument("--workers", "-w", type=int, default=1, help="number of workers")
    args = parser.parse_args()

    eval_pt_center(args.config_file, args.workers)
