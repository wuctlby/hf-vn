import argparse
import sys
import os
import re
import array
import time
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # pylint: disable=wrong-import-position
import pandas as pd
import numpy as np
import seaborn as sns
import awkward as ak
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
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
import tensorflow as tf
tf.config.threading.set_intra_op_parallelism_threads(20)
tf.config.threading.set_inter_op_parallelism_threads(20)

def load_aod_file(aod_file, has_sp_cent, num_workers=16, chunk_size=1_000_000):
    start_total = time.time()
    branches = ["fPt", "fM", "fMlScore0", "fMlScore1", "fScalarProd", "fCent"] \
               if has_sp_cent else ["fPt", "fM", "fMlScore0", "fMlScore1"]
    key = "ptcentersp" if has_sp_cent else "ptcenter"

    # Open with parallel decompression
    t0 = time.time()
    f = uproot.open(aod_file, num_workers=8)
    entries = f[key].num_entries
    logger(f"Opened file in {time.time()-t0:.2f}s using {num_workers} workers, total entries: {entries}", level="INFO")

    # Iterate in chunks
    dfs = []
    for df_chunk in f[key].iterate(filter_name=branches, step_size=chunk_size, library="np"):
        dfs.append(pd.DataFrame(df_chunk))

    # Concatenate
    t2 = time.time()
    df = pd.concat(dfs, ignore_index=True)
    logger(f"TOTAL time: {time.time()-start_total:.2f}s", level="INFO")
    return df

def eval_pt_center(cfg_file_name, minimizer, workers=1):
    # Read the configuration file
    with open(cfg_file_name, 'r') as cfg_file:
        cfg = yaml.safe_load(cfg_file)
    print(tf.config.threading.get_intra_op_parallelism_threads())

    # Retrieve cutsets configs
    cutsets_dir = os.path.join(cfg['outdir'], f"cutvar_{cfg['suffix']}_combined/cutsets")
    cutset_files = [os.path.join(cutsets_dir, f) for f in os.listdir(cutsets_dir) if f.endswith('.yml')]
    cutset_files.sort(key=lambda x: int(re.search(r'(\d+)', os.path.basename(x)).group(1)))

    infer_vars = ['fPt', 'fScalarProd', 'fCent']

    # Loop over cutset configs
    s_weights = {}
    cfg_fit = cfg["v2extraction"]

    for i_bin, (pt_min, pt_max) in enumerate(zip(cfg["ptbins"][:-1], cfg["ptbins"][1:])):
        logger(f"Pt bin {i_bin}: {pt_min} - {pt_max}", level="INFO")

        # Load input
        pt_str = f"pt_{int(pt_min*10)}_{int(pt_max*10)}"
        has_sp_cent = True
        prep_dir = f"{cfg['outdir']}/preprocess/{pt_str}/TreesPtCenterSp"
        if not os.path.exists(prep_dir):
            prep_dir = f"{cfg['outdir']}/preprocess/{pt_str}/TreesPtCenter"
            has_sp_cent = False
            logger(f"Using tree without SP and centrality!", level="WARNING")
        df = load_aod_file(f"{prep_dir}/AO2D_{pt_str}.root", has_sp_cent)

        out_dir_pt = f"{cfg['outdir']}/cutvar_{cfg['suffix']}_combined/ptcenter_{minimizer}/{pt_str}"
        os.makedirs(out_dir_pt, exist_ok=True)
        os.makedirs(out_dir_pt + "/fits", exist_ok=True)
        os.makedirs(out_dir_pt + "/vars", exist_ok=True)

        out_file = TFile.Open(f"{out_dir_pt}/pt_center.root", "recreate")
        histos_avgs = {}
        for var in infer_vars:
            histos_avgs[f"h_{var}_sgn"] = ROOT.TH1F(f"h_{var}_sgn", f"h_{var}_sgn", len(cutset_files)-1, array.array('d', [i for i in range(len(cutset_files))]))
            histos_avgs[f"h_{var}_bkg"] = ROOT.TH1F(f"h_{var}_bkg", f"h_{var}_bkg", len(cutset_files)-1, array.array('d', [i for i in range(len(cutset_files))]))

        sgn_funcs = {} # More info for signal functions, a dictionary is better
        sgn_funcs[cfg_fit['SgnFuncLabel']] = {
            'func': cfg_fit['SgnFunc'][i_bin] if isinstance(cfg_fit['SgnFunc'], list) else cfg_fit['SgnFunc'],
            'part': cfg['Dmeson']
        }

        histos_avgs[f"h_ry_{cfg_fit['SgnFuncLabel']}"] = ROOT.TH1F(f"h_ry_{cfg_fit['SgnFuncLabel']}",
                                                                   f"h_ry_{cfg_fit['SgnFuncLabel']}",
                                                                   len(cutset_files)-1,
                                                                   array.array('d', [i for i in range(len(cutset_files))]))
        print(f"Adding signal function: {sgn_funcs[cfg_fit['SgnFuncLabel']]}, {cfg_fit['SgnFuncLabel']} ... ")
        if cfg_fit.get('InclSecPeak'):
            print("Including secondary peak signal function ... ")
            include_sec_peak = cfg_fit['InclSecPeak'][i_bin] if isinstance(cfg_fit['InclSecPeak'], list) else cfg_fit['InclSecPeak']
            print(f"include_sec_peak = {include_sec_peak}")
            if include_sec_peak:
                print(f"Adding secondary peak signal function: {cfg_fit['SgnFuncSecPeak'][i_bin]} ... ")
                sgn_funcs[cfg_fit['SgnFuncSecPeakLabel']] = {
                    'func': cfg_fit['SgnFuncSecPeak'][i_bin] if isinstance(cfg_fit['SgnFuncSecPeak'], list) else cfg_fit['SgnFuncSecPeak'],
                    'part': 'Dplus' if cfg['Dmeson'] == 'Ds' else 'Dstar',
                }
                histos_avgs[f"h_ry_{cfg_fit['SgnFuncSecPeakLabel']}"] = ROOT.TH1F(f"h_ry_{cfg_fit['SgnFuncSecPeakLabel']}",
                                                                                  f"h_ry_{cfg_fit['SgnFuncSecPeakLabel']}",
                                                                                  len(cutset_files)-1, 
                                                                                  array.array('d', [i for i in range(len(cutset_files))]))

        print(f"\n\n sgn_funcs = {sgn_funcs}\n\n")

        # Initialize fitter
        fitter = RawYieldFitter(cfg['Dmeson'], pt_min, pt_max, pt_str, minimizer)
        fitter.set_fit_range(cfg_fit['MassFitRanges'][i_bin][0], cfg_fit['MassFitRanges'][i_bin][1])

        for i_cutset, cutset_file in enumerate(cutset_files):
            with open(cutset_file, 'r') as cs_file:
                cutset_cfg = yaml.safe_load(cs_file)
            logger(f"Processing cutset file {cutset_file} ... ", level="INFO")
            cutset_suffix = os.path.basename(cutset_file).replace('.yml', '').split('_')[-1]

            # Setup fitter
            fitter.add_bkg_func(cfg_fit['BkgFunc'][i_bin] if isinstance(cfg_fit['BkgFunc'], list) else cfg_fit['BkgFunc'], "Comb. bkg")
            for i_sgn, (label, sgn_func) in enumerate(sgn_funcs.items()):
                print(f"Adding signal function: {sgn_func}, {label} ... ")
                fitter.add_sgn_func(sgn_func['func'], label, sgn_func['part'])

            score_bkg_min = cutset_cfg["ScoreBkg"]["min"][i_bin]
            score_bkg_max = cutset_cfg["ScoreBkg"]["max"][i_bin]
            score_fd_min = cutset_cfg["ScoreFD"]["min"][i_bin]
            score_fd_max = cutset_cfg["ScoreFD"]["max"][i_bin]

            # Query the dataframe
            mass_min, mass_max = cfg_fit["MassFitRanges"][i_bin]
            sel_string = f"fPt >= {pt_min} and fPt < {pt_max} and " \
                         f"fMlScore0 >= {score_bkg_min} and fMlScore0 < {score_bkg_max} and " \
                         f"fMlScore1 >= {score_fd_min} and fMlScore1 < {score_fd_max} and " \
                         f"fM >= {mass_min} and fM < {mass_max}"
            sel_df = df.query(sel_string).reset_index(drop=True)
            fitter.set_name(f"{pt_str}_{cutset_suffix}")
            fitter.set_data_to_fit_df(sel_df)
            if cfg_fit.get('Rebin'):
                fitter.set_rebin(cfg_fit['Rebin'][i_bin]) if isinstance(cfg_fit['Rebin'], list) else fitter.set_rebin(cfg_fit['Rebin'])

            # Add correlated background if specified
            if cfg.get('corr_bkgs'):
                fitter.add_corr_bkgs(cfg['corr_bkgs'], sel_string.replace(' and ', ' && '), pt_min, pt_max)

            fig, ax = plt.subplots(1, 1, figsize=(12, 10))
            sel_df['fM'].hist(bins=100, alpha=0.5, range=(mass_min, mass_max))
            ax.set_xlabel('fM')
            ax.set_ylabel('Counts')
            fig.savefig(f"{out_dir_pt}/fits/fM_raw_{cutset_suffix}.pdf", dpi=300, bbox_inches="tight")

            fitter.setup()

            if cfg_fit.get('InitPars'):
                fitter.set_fit_pars(cfg_fit['InitPars'], pt_min, pt_max)

            # Prefit the MC prompt enhanced cut to fix the tails, binned fit
            if cfg_fit.get('FixSgnFromMC'):
                fitter.set_fix_sgn_to_mc_prefit(True)
                if i_cutset == 0:
                    fitter.prefit_mc(f"{cfg['outdir']}/corr_bkgs/templs_{pt_str}.root")
                    fitter.plot_mc_prefit(False, True, loc=["lower left", "upper left"],
                                          path=f"{out_dir_pt}/", out_file=out_file)
                    fitter.plot_raw_residuals_mc_prefit(path=f"{out_dir_pt}/fM_mc_prefit_residuals_{cutset_suffix}.pdf")

            status, converged = fitter.fit()

            fitter.plot_fit(False, True, loc=["lower left", "upper left"], \
                            path=f"{out_dir_pt}/fits/fM_fit_{cutset_suffix}.pdf",
                            out_file=out_file) # (log, show_extra_info)

            fit_info, _, _, _, _ = fitter.get_fit_info()
            for label in sgn_funcs.keys():
                histos_avgs[f"h_ry_{label}"].SetBinContent(i_cutset + 1, fit_info[label]["ry"])
                histos_avgs[f"h_ry_{label}"].SetBinError(i_cutset + 1, fit_info[label]["ry_unc"])

            if minimizer != "flarefly":
                logger("Skipping sWeights computation: not using flarefly minimizer", level="WARNING")
                continue

            s_weights_sgn = fitter.get_sweights_sgn(cfg_fit['SgnFuncLabel'])
            s_weights_sec_peak = fitter.get_sweights_sgn(cfg_fit['SgnFuncSecPeakLabel']) if cfg_fit.get('InclSecPeak') else None
            with open(f"{out_dir_pt}/fits_status.txt", "a") as f:
                computed_sweights = True if s_weights_sgn is not None else False
                f.write(
                        f"{fitter.get_name()}: "
                        f"fit_res.status -> {status}, "
                        f"fit_res.converged -> {converged}, "
                        f"sweights computed -> {computed_sweights} \n"
                        )

            for var in infer_vars:
                print(f"    Drawing {var}")

                # Create figure with two subplots (distros and ratio)
                fig, ax = plt.subplots(figsize=(8, 8))

                bins = 200
                var_range = (min(sel_df[var]), max(sel_df[var]))
                sgn_vals, bin_edges = np.histogram(sel_df[var], bins=bins, range=var_range, weights=s_weights_sgn, density=True)
                if s_weights_sec_peak is not None:
                    bkg_sweights = np.ones(len(s_weights_sgn)) - np.asarray(s_weights_sgn) - np.asarray(s_weights_sec_peak)
                else:
                    bkg_sweights = np.ones(len(s_weights_sgn)) - np.asarray(s_weights_sgn)
                bkg_vals, _ = np.histogram(sel_df[var], bins=bins, range=var_range, weights=bkg_sweights, density=True)
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
                    os.path.join(out_dir_pt, f'vars/{var}_{cutset_suffix}.pdf'),
                    dpi=300, bbox_inches="tight"
                )
                plt.close(fig)

                # Fill histograms for averages
                histos_avgs[f"h_{var}_sgn"].SetBinContent(i_cutset + 1, np.average(sel_df[var], weights=s_weights_sgn))
                histos_avgs[f"h_{var}_sgn"].SetBinError(i_cutset + 1, np.std(sgn_vals))
                histos_avgs[f"h_{var}_bkg"].SetBinContent(i_cutset + 1, np.average(sel_df[var], weights=bkg_sweights))
                histos_avgs[f"h_{var}_bkg"].SetBinError(i_cutset + 1, np.std(bkg_vals))

            # Reset fitter for new cutset: correlated bkg fracs will change
            fitter.reset()

        # Compute average pt
        histos_avgs["h_avg_pt"] = ROOT.TH1F("h_avg_pt", "h_avg_pt", 1, 0, 1)
        avg_pt = 0
        for i_bin in range(len(cutset_files)-1):
            avg_pt += (histos_avgs["h_fPt_sgn"].GetBinContent(i_bin + 1) * \
                       histos_avgs[f"h_ry_{cfg_fit['SgnFuncLabel']}"].GetBinContent(i_bin + 1)) / \
                       histos_avgs[f"h_ry_{cfg_fit['SgnFuncLabel']}"].Integral()

        histos_avgs["h_avg_pt"].SetBinContent(1, avg_pt)

        # Write histograms to output ROOT file
        out_file.cd()
        for hist in histos_avgs.values():
            hist.Write()
        out_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate pt-centering with sPlot and flarefly/roofit')
    parser.add_argument('config_file', help='Path to the input configuration file')
    parser.add_argument("--minimizer", "-m", type=str, default="flarefly", help="minimizer to use")
    parser.add_argument("--workers", "-w", type=int, default=1, help="number of workers")
    args = parser.parse_args()

    eval_pt_center(args.config_file, args.minimizer, args.workers)
