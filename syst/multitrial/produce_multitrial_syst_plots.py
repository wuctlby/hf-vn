import sys
import os
import numpy as np
import argparse
import re
import glob
import array
import pandas as pd
import yaml
import ast
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kError  # Only show errors and above
from ROOT import TFile, TCanvas, TH1F, TGraphAsymmErrors, TLegend, kOrange, kAzure, kBlack
script_dir = os.path.dirname(os.path.realpath(__file__))
os.sys.path.append(os.path.join(script_dir, '../..', 'utils'))
from utils import logger
from StyleFormatter import SetGlobalStyle, SetObjectStyle
from PyPDF2 import PdfMerger
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

SetGlobalStyle(titleoffsety=1.1, maxdigits=3, topmargin=0.1, bottommargin=0.4, leftmargin=0.3, rightmargin=0.15,
               labelsizey=0.04, setoptstat=0, setopttitle=0, setdecimals=True,titleoffsetx=0.74)

def plot_populations(trials_df, unique_cols, pt_min, pt_max, output_dir):
    """
    Produce publication-quality 2-panel figures:
      - Left: individual trials with jitter
      - Right: binned distributions + means

    One PDF per cutset and per Prompt / NonPrompt.
    Each PDF contains one page per column in unique_cols.
    """

    jitter_val = 0.45

    if "VnVsMassBins" in unique_cols:
        trials_df["VnVsMassBinWidths (MeV/$c^2$)"] = trials_df["VnVsMassBins"].apply(
            lambda s: int(round((ast.literal_eval(s)[1] - ast.literal_eval(s)[0]) * 1000))
        )

        trials_df = trials_df.drop(columns=["VnVsMassBins"])
        unique_cols[unique_cols.index("VnVsMassBins")] = "VnVsMassBinWidths (MeV/$c^2$)"

    # Make grouping columns hashable
    for col in unique_cols:
        trials_df[col] = trials_df[col].apply(
            lambda x: str(x) if isinstance(x, (list, tuple)) else x
        )

    # Global style
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 15,
        "legend.fontsize": 12,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13
    })

    colors = plt.cm.tab10.colors
    cutsets = sorted(trials_df["Cutset"].unique())

    # ============================================================
    #  CUTSET PDFs
    # ============================================================
    for cutset in cutsets:
        df_cut = trials_df[trials_df["Cutset"] == cutset]
        df_trials = df_cut[df_cut["TrialIdx"] != "Reference"]
        df_ref = df_cut[df_cut["TrialIdx"] == "Reference"]

        pdf_path = f"{output_dir}/V2_cutset_{cutset}.pdf"

        with PdfPages(pdf_path) as pdf:
            for col in unique_cols:
                
                groups = sorted(df_trials[col].unique())
                x_ticks_labels = groups
                # Setup labels for mass range
                if col == "MassFitRanges":
                    x_ticks_labels = [g.replace(",", "-").replace('[', '').replace(']', '') for g in groups]

                x_pos = np.arange(len(groups))

                fig = plt.figure(figsize=(11, 5.5))
                gs = fig.add_gridspec(1, 2, width_ratios=[1, 1])

                ax_scatter = fig.add_subplot(gs[0, 0])
                ax_hist = fig.add_subplot(gs[0, 1])

                # ------------------ TITLE ------------------
                fig.suptitle(
                    rf"{col} variations, cutset {cutset}, "
                    rf"${pt_min} < p_T < {pt_max}$",
                    y=0.96
                )

                # ================== SCATTER ==================
                for i, g in enumerate(groups):
                    df_g = df_trials[df_trials[col] == g]
                    jitter = np.random.uniform(-jitter_val, jitter_val, len(df_g))

                    ax_scatter.errorbar(x_pos[i] + jitter,
                        df_g["V2"], yerr=df_g["V2Unc"],
                        fmt="o", color=colors[i % len(colors)],
                        alpha=0.85, markersize=6, zorder=3, label=x_ticks_labels[i]
                    )

                ax_scatter.set_xticks([])
                ax_scatter.set_xlim(-0.6, len(groups) - 0.4)
                ax_scatter.set_xlabel(col, labelpad=12)
                ax_scatter.set_ylabel(r"$v_2$")
                ax_scatter.grid(True, linestyle="--", alpha=0.25)
                ax_scatter.set_box_aspect(1)

                # ---------- SINGLE REFERENCE (ON TOP) ----------
                if not df_ref.empty:
                    ref_val = df_ref["V2"].iloc[0]
                    ref_unc = df_ref["V2Unc"].iloc[0]

                    xmin, xmax = ax_scatter.get_xlim()
                    ax_scatter.fill_between(
                        [xmin, xmax], ref_val - ref_unc, ref_val + ref_unc,
                        color="black", alpha=0.12, zorder=4
)
                    ax_scatter.axhline(
                        ref_val, color="black",
                        linestyle="--", linewidth=2.3, zorder=5
                    )

                # ================== HISTOGRAM ==================
                vmin = df_trials["V2"].min() - 0.02
                vmax = df_trials["V2"].max() + 0.02
                bins = np.linspace(vmin, vmax, 120)

                for i, g in enumerate(groups):
                    df_g = df_trials[df_trials[col] == g]

                    ax_hist.hist(
                        df_g["V2"], bins=bins, color=colors[i % len(colors)],
                        alpha=0.55, edgecolor="black", linewidth=0.8
                    )

                    ax_hist.axvline(
                        df_g["V2"].mean(), color=colors[i % len(colors)],
                        linestyle="--", linewidth=2, zorder=3
                    )

                ax_hist.set_xlabel(r"$v_2$")
                ax_hist.set_ylabel("Counts")
                ax_hist.grid(True, linestyle="--", alpha=0.25)
                ax_hist.set_box_aspect(1)

                # ---------- SINGLE HISTO REFERENCE ----------
                if not df_ref.empty:
                    ymin, ymax = ax_hist.get_ylim()
                    ax_hist.fill_betweenx(
                        [ymin, ymax], ref_val - ref_unc, ref_val + ref_unc,
                        color="black", alpha=0.12, zorder=4
                    )
                    ax_hist.axvline(
                        ref_val, color="black", linestyle="--",
                        linewidth=2.3, zorder=5
                    )

                # Legend only once
                ax_scatter.legend(ncol=2, frameon=False)

                # ---------- LAYOUT CONTROL ----------
                fig.subplots_adjust(top=0.90, bottom=0.12, left=0.08, right=0.98, wspace=0.25)

                pdf.savefig(fig)
                plt.close(fig)

        print(f"✔ Saved {pdf_path}")

    # ============================================================
    #  PROMPT / NON-PROMPT PDFs
    # ============================================================
    for var, label in [
        ("PromptV2", r"Prompt $v_2$"),
        ("NonPromptV2", r"Non-Prompt $v_2$")
    ]:
        df_trials = trials_df[trials_df["TrialIdx"] != "Reference"]
        df_ref = trials_df[trials_df["TrialIdx"] == "Reference"]

        pdf_path = f"{output_dir}/{var}.pdf"

        with PdfPages(pdf_path) as pdf:
            for col in unique_cols:
                groups = sorted(df_trials[col].unique())
                x_pos = np.arange(len(groups))

                fig = plt.figure(figsize=(11, 5.5))
                gs = fig.add_gridspec(1, 2)

                ax_scatter = fig.add_subplot(gs[0, 0])
                ax_hist = fig.add_subplot(gs[0, 1])

                fig.suptitle(
                    rf"{col} variations, {label}, "
                    rf"${pt_min} < p_T < {pt_max}$",
                    y=0.96
                )

                for i, g in enumerate(groups):
                    df_g = df_trials[df_trials[col] == g]
                    jitter = np.random.uniform(-jitter_val, jitter_val, len(df_g))

                    ax_scatter.errorbar(
                        x_pos[i] + jitter,
                        df_g[var], yerr=df_g[f"{var}Unc"],
                        fmt="o", color=colors[i % len(colors)],
                        alpha=0.85, markersize=6, zorder=3, label=g
                    )

                ax_scatter.set_xticks([])
                ax_scatter.set_xlim(-0.6, len(groups) - 0.4)
                ax_scatter.set_xlabel(col, labelpad=12)
                ax_scatter.set_ylabel(label)
                ax_scatter.grid(True, linestyle="--", alpha=0.25)
                ax_scatter.set_box_aspect(1)

                # Global reference
                ref_val = df_ref[var].mean()
                ref_unc = df_ref[f"{var}Unc"].mean()

                xmin, xmax = ax_scatter.get_xlim()
                ax_scatter.fill_between(
                    [xmin, xmax], ref_val - ref_unc, ref_val + ref_unc,
                    color="black", alpha=0.12, zorder=4
                )
                ax_scatter.axhline(
                    ref_val, color="black", linestyle="--",
                    linewidth=2.3, zorder=5
                )

                bins = np.linspace(df_trials[var].min() - 0.02, df_trials[var].max() + 0.02, 120)

                for i, g in enumerate(groups):
                    df_g = df_trials[df_trials[col] == g]
                    ax_hist.hist(
                        df_g[var], bins=bins, color=colors[i % len(colors)],
                        alpha=0.55, edgecolor="black", linewidth=0.8
                    )
                    ax_hist.axvline(
                        df_g[var].mean(), color=colors[i % len(colors)],
                        linestyle="--", linewidth=2
                    )

                ymin, ymax = ax_hist.get_ylim()
                ax_hist.fill_betweenx(
                    [ymin, ymax], ref_val - ref_unc, ref_val + ref_unc,
                    color="black", alpha=0.12, zorder=4
                )
                ax_hist.axvline(
                    ref_val, color="black", linestyle="--",
                    linewidth=2.3, zorder=5
                )

                ax_hist.set_xlabel(r"$v_2$")
                ax_hist.set_ylabel("Counts")
                ax_hist.grid(True, linestyle="--", alpha=0.25)
                ax_hist.set_box_aspect(1)

                ax_scatter.legend(ncol=2, frameon=False)

                fig.subplots_adjust(top=0.90, bottom=0.12, left=0.08, right=0.98, wspace=0.25)

                pdf.savefig(fig)
                plt.close(fig)

        print(f"✔ Saved {pdf_path}")


def plot_variations_multitrial(cfg_ref, v2_vs_frac_ref, ry_cutsets_ref, mult_dir, out_dir, syst_uncs):

    print(f"PRODUCING VARIATIONS PLOTS FOR: {mult_dir}")
    # Retrieve the final prompt and FD v2 from reference results for
    # the selected pt bin and the v2 for the single cutsets
    pt_bin_ref = None
    _, pt_min_times_10, pt_max_times_10 = mult_dir.split('_')
    pt_min = float(pt_min_times_10) / 10.
    pt_max = float(pt_max_times_10) / 10.
    pt_center = (pt_min + pt_max) / 2.
    reference_results = {key: None for key in ry_cutsets_ref}
    for cutset in ry_cutsets_ref.keys():
        h_vn_cutset = ry_cutsets_ref[cutset].Get('hVnSimFit')
        h_signif_cutset = ry_cutsets_ref[cutset].Get('hRawYieldsSignificanceSimFit')
        if pt_bin_ref is None:
            for i_bin in range(1, h_vn_cutset.GetNbinsX()+1):
                pt_low_edge = h_vn_cutset.GetBinLowEdge(i_bin)
                pt_up_edge = pt_low_edge + h_vn_cutset.GetBinWidth(i_bin)
                if pt_low_edge <= pt_center < pt_up_edge:
                    pt_bin_ref = i_bin
                    break
        reference_results[cutset] = {}
        reference_results[cutset]['V2']  = h_vn_cutset.GetBinContent(pt_bin_ref)
        reference_results[cutset]['V2Unc'] = h_vn_cutset.GetBinError(pt_bin_ref)
        reference_results[cutset]['Significance'] = h_signif_cutset.GetBinContent(pt_bin_ref)
        reference_results[cutset]['SignificanceUnc'] = h_signif_cutset.GetBinError(pt_bin_ref)
        reference_results[cutset]['Chi2'] = ry_cutsets_ref[cutset].Get('hRedChi2SimFit').GetBinContent(pt_bin_ref)
        reference_results[cutset]['Chi2Unc'] = ry_cutsets_ref[cutset].Get('hRedChi2SimFit').GetBinError(pt_bin_ref)

    trials_dirs = glob.glob(f"{args.results_dir}/syst/multitrial/{mult_dir}/trials/*")

    ry_trials_results = {key: None for key in ry_cutsets_ref}
    for key in ry_trials_results:
        ry_trials_results[key] = {
            'Significances': [],
            'SignificancesUncs': [],
            'Chi2s': [],
            'Chi2sUncs': [],
            'V2s': [],
            'V2sUncs': [],
            'BkgFuncs': [],
            'BkgFuncsVn': [],
            'SgnFuncs': [],
            'Rebins': [],
            'VnVsMassBins': [],
            'MassFitRanges': [],
            'TrialIdxs': [],
            'PromptV2': [],
            'PromptV2Unc': [],
            'NonPromptV2': [],
            'NonPromptV2Unc': []
        }

    v2_prompt_ref = v2_vs_frac_ref.Get('hV2VsPtPrompt').GetBinContent(pt_bin_ref)
    v2_prompt_ref_unc = v2_vs_frac_ref.Get('hV2VsPtPrompt').GetBinError(pt_bin_ref)
    v2_non_prompt_ref = v2_vs_frac_ref.Get('hV2VsPtFD').GetBinContent(pt_bin_ref)
    v2_non_prompt_ref_unc = v2_vs_frac_ref.Get('hV2VsPtFD').GetBinError(pt_bin_ref)
    v2_prompt_trials, v2_prompt_trials_uncs, v2_non_prompt_trials, v2_non_prompt_trials_uncs = [], [], [], []
    for trial_dir in trials_dirs:
        trial_idx = os.path.split(trial_dir)[-1]
        try:
            v2_vs_frac_trial = TFile.Open(f"{trial_dir}/v2/v2VsFrac.root", "read")
            v2_prompt_trials.append(v2_vs_frac_trial.Get('hV2VsPtPrompt').GetBinContent(1))
            v2_prompt_trials_uncs.append(v2_vs_frac_trial.Get('hV2VsPtPrompt').GetBinError(1))
            v2_non_prompt_trials.append(v2_vs_frac_trial.Get('hV2VsPtFD').GetBinContent(1))
            v2_non_prompt_trials_uncs.append(v2_vs_frac_trial.Get('hV2VsPtFD').GetBinError(1))

            ry_files = [f for f in os.listdir(f"{trial_dir}/raw_yields/") if re.match(r"raw_yields_\d+\.root$", f)]
            cutsets = sorted((re.search(r"raw_yields_(\d+)\.root", f).group(1) for f in ry_files), key=int)
            ry_cutsets_dict_trial = {
                cutset: TFile.Open(os.path.join(f"{trial_dir}/raw_yields/", f"raw_yields_{cutset}.root"), "READ")
                for cutset in cutsets
            }

            # Retrieve fit configuration
            trial_number = os.path.split(trial_dir)[-1]
            print(f"Opening {trial_dir}/config_trial_{trial_number}.yml")
            with open(f"{trial_dir}/config_trial_{trial_number}.yml", 'r') as CfgFlow:
                cfg_trial = yaml.safe_load(CfgFlow)
                trial_bkg_func = cfg_trial['v2extraction']['BkgFunc'][0] if isinstance(cfg_trial['v2extraction']['BkgFunc'], list) else cfg_trial['v2extraction']['BkgFunc']
                trial_bkg_func_vn = cfg_trial['v2extraction']['BkgFuncVn'][0] if isinstance(cfg_trial['v2extraction']['BkgFuncVn'], list) else cfg_trial['v2extraction']['BkgFuncVn']
                trial_sgn_func = cfg_trial['v2extraction']['SgnFunc'][0] if isinstance(cfg_trial['v2extraction']['SgnFunc'], list) else cfg_trial['v2extraction']['SgnFunc']
                trial_rebins = cfg_trial['v2extraction']['Rebin'][0]
                trial_vn_vs_mass_bins = cfg_trial['projections']['inv_mass_bins'][0]
                trial_mass_fit_range = cfg_trial['v2extraction']['MassFitRanges'][0]

            for cutset, file in ry_cutsets_dict_trial.items():
                ry_trials_results[cutset]['Significances'].append(file.Get('hRawYieldsSignificanceSimFit').GetBinContent(1))
                ry_trials_results[cutset]['SignificancesUncs'].append(file.Get('hRawYieldsSignificanceSimFit').GetBinError(1))
                ry_trials_results[cutset]['Chi2s'].append(file.Get('hRedChi2SimFit').GetBinContent(1))
                ry_trials_results[cutset]['Chi2sUncs'].append(file.Get('hRedChi2SimFit').GetBinError(1))
                ry_trials_results[cutset]['V2s'].append(file.Get('hVnSimFit').GetBinContent(1))
                ry_trials_results[cutset]['V2sUncs'].append(file.Get('hVnSimFit').GetBinError(1))
                ry_trials_results[cutset]['BkgFuncs'].append(trial_bkg_func)
                ry_trials_results[cutset]['BkgFuncsVn'].append(trial_bkg_func_vn)
                ry_trials_results[cutset]['SgnFuncs'].append(trial_sgn_func)
                ry_trials_results[cutset]['Rebins'].append(trial_rebins)
                ry_trials_results[cutset]['VnVsMassBins'].append(trial_vn_vs_mass_bins)
                ry_trials_results[cutset]['MassFitRanges'].append(trial_mass_fit_range)
                ry_trials_results[cutset]['TrialIdxs'].append(trial_idx)
                ry_trials_results[cutset]['PromptV2'].append(v2_prompt_trials[-1])
                ry_trials_results[cutset]['PromptV2Unc'].append(v2_prompt_trials_uncs[-1])
                ry_trials_results[cutset]['NonPromptV2'].append(v2_non_prompt_trials[-1])
                ry_trials_results[cutset]['NonPromptV2Unc'].append(v2_non_prompt_trials_uncs[-1])

        except Exception as e:
            logger(f"Error opening files for trial_dir {trial_dir}: {e}", "WARNING")
            continue

    # Flatten ry_trials_results into a DataFrame
    trials_df = []
    for cutset, data in ry_trials_results.items():
        n_trials = len(data['V2s'])
        for i in range(n_trials):
            trials_df.append({
                "Cutset": cutset,
                "TrialIdx": data['TrialIdxs'][i],
                "V2": data['V2s'][i],
                "V2Unc": data['V2sUncs'][i],
                "Significance": data['Significances'][i],
                "SignificanceUnc": data['SignificancesUncs'][i],
                "Chi2": data['Chi2s'][i],
                "Chi2Unc": data['Chi2sUncs'][i],
                "BkgFunc": data['BkgFuncs'][i],
                "BkgFuncVn": data['BkgFuncsVn'][i],
                "SgnFunc": data['SgnFuncs'][i],
                "Rebin": data['Rebins'][i],
                "VnVsMassBins": str(data['VnVsMassBins'][i]),  # convert lists to string
                "MassFitRanges": data['MassFitRanges'][i],
                "PromptV2": data['PromptV2'][i],
                "PromptV2Unc": data['PromptV2Unc'][i],
                "NonPromptV2": data['NonPromptV2'][i],
                "NonPromptV2Unc": data['NonPromptV2Unc'][i]
            })

    # Insert the reference results as trial_index -1
    for cutset, ref_data in reference_results.items():
        trials_df.append({
            "Cutset": cutset,
            "TrialIdx": "Reference",
            "V2": ref_data['V2'],
            "V2Unc": ref_data['V2Unc'],
            "Significance": ref_data['Significance'],
            "SignificanceUnc": ref_data['SignificanceUnc'],
            "Chi2": ref_data['Chi2'],
            "Chi2Unc": ref_data['Chi2Unc'],
            "BkgFunc": cfg_ref['v2extraction']['BkgFunc'][pt_bin_ref-1] if isinstance(cfg_ref['v2extraction']['BkgFunc'], list) else cfg_ref['v2extraction']['BkgFunc'],
            "BkgFuncVn": cfg_ref['v2extraction']['BkgFuncVn'][pt_bin_ref-1] if isinstance(cfg_ref['v2extraction']['BkgFuncVn'], list) else cfg_ref['v2extraction']['BkgFuncVn'],
            "SgnFunc": cfg_ref['v2extraction']['SgnFunc'][pt_bin_ref-1] if isinstance(cfg_ref['v2extraction']['BkgFunc'], list) else cfg_ref['v2extraction']['SgnFunc'],
            "Rebin": cfg_ref['v2extraction']['Rebin'][pt_bin_ref-1] if isinstance(cfg_ref['v2extraction']['Rebin'], list) else cfg_ref['v2extraction']['Rebin'],
            "VnVsMassBins": str(cfg_ref['projections']['inv_mass_bins'][pt_bin_ref-1]),  # convert lists to string
            "MassFitRanges": cfg_ref['v2extraction']['MassFitRanges'][pt_bin_ref-1],
            "PromptV2": v2_prompt_ref,
            "PromptV2Unc": v2_prompt_ref_unc,
            "NonPromptV2": v2_non_prompt_ref,
            "NonPromptV2Unc": v2_non_prompt_ref_unc
        })

    trials_df = pd.DataFrame(trials_df)
    trials_df.to_parquet(f"{out_dir}/ry_trials_results_{mult_dir}.parquet", index=False)

    # plot_populations(trials_df, ['BkgFunc', 'SgnFunc', 'BkgFuncVn', 'Rebin'], pt_min, pt_max,
    plot_populations(trials_df, ['BkgFunc', 'SgnFunc', 'BkgFuncVn', 'Rebin', 'MassFitRanges', 'VnVsMassBins'],
                     pt_min, pt_max, output_dir=f"{out_dir}/pt_{pt_min_times_10}_{pt_max_times_10}")

def compute_syst_multitrial(v2_vs_frac_ref, ry_cutsets_ref, mult_dir, out_dir):

    # Retrieve the final prompt and FD v2 from reference results for
    # the selected pt bin and the v2 for the single cutsets
    pt_bin_ref = None
    _, pt_min_times_10, pt_max_times_10 = mult_dir.split('_')
    pt_min = float(pt_min_times_10) / 10.
    pt_max = float(pt_max_times_10) / 10.
    pt_center = (pt_min + pt_max) / 2.
    reference_results = {key: None for key in ry_cutsets_ref}
    for cutset in ry_cutsets_ref.keys():
        h_vn_cutset = ry_cutsets_ref[cutset].Get('hVnSimFit')
        if pt_bin_ref is None:
            for i_bin in range(1, h_vn_cutset.GetNbinsX()+1):
                pt_low_edge = h_vn_cutset.GetBinLowEdge(i_bin)
                pt_up_edge = pt_low_edge + h_vn_cutset.GetBinWidth(i_bin)
                if pt_low_edge <= pt_center < pt_up_edge:
                    pt_bin_ref = i_bin
                    break
        reference_results[cutset] = {}
        reference_results[cutset]['V2']  = h_vn_cutset.GetBinContent(pt_bin_ref)
        reference_results[cutset]['V2Unc'] = h_vn_cutset.GetBinError(pt_bin_ref)

    trials_dirs = glob.glob(f"{args.results_dir}/syst/multitrial/{mult_dir}/trials/*")

    ry_trials_results = {key: None for key in ry_cutsets_ref}
    for key in ry_trials_results:
        ry_trials_results[key] = {
            'Significances': [],
            'SignificancesUncs': [],
            'Chi2s': [],
            'Chi2sUncs': [],
            'V2s': [],
            'V2sUncs': []
        }

    v2_prompt_trials, v2_prompt_trials_uncs, v2_non_prompt_trials, v2_non_prompt_trials_uncs = [], [], [], []
    for trial_dir in trials_dirs:
        try:
            v2_vs_frac_trial = TFile.Open(f"{trial_dir}/v2/v2VsFrac.root", "read")
            v2_prompt_trials.append(v2_vs_frac_trial.Get('hV2VsPtPrompt').GetBinContent(1))
            v2_prompt_trials_uncs.append(v2_vs_frac_trial.Get('hV2VsPtPrompt').GetBinError(1))
            v2_non_prompt_trials.append(v2_vs_frac_trial.Get('hV2VsPtFD').GetBinContent(1))
            v2_non_prompt_trials_uncs.append(v2_vs_frac_trial.Get('hV2VsPtFD').GetBinError(1))
            ry_files = [f for f in os.listdir(f"{trial_dir}/raw_yields/") if re.match(r"raw_yields_\d+\.root$", f)]
            cutsets = sorted((re.search(r"raw_yields_(\d+)\.root", f).group(1) for f in ry_files), key=int)
            ry_cutsets_dict_trial = {
                cutset: TFile.Open(os.path.join(f"{trial_dir}/raw_yields/", f"raw_yields_{cutset}.root"), "READ")
                for cutset in cutsets
            }

            for cutset, file in ry_cutsets_dict_trial.items():
                ry_trials_results[cutset]['Significances'].append(file.Get('hRawYieldsSignificanceSimFit').GetBinContent(1))
                ry_trials_results[cutset]['SignificancesUncs'].append(file.Get('hRawYieldsSignificanceSimFit').GetBinError(1))
                ry_trials_results[cutset]['Chi2s'].append(file.Get('hRedChi2SimFit').GetBinContent(1))
                ry_trials_results[cutset]['Chi2sUncs'].append(file.Get('hRedChi2SimFit').GetBinError(1))
                ry_trials_results[cutset]['V2s'].append(file.Get('hVnSimFit').GetBinContent(1))
                ry_trials_results[cutset]['V2sUncs'].append(file.Get('hVnSimFit').GetBinError(1))

        except Exception as e:
            logger(f"Error opening files for trial_dir {trial_dir}: {e}", "WARNING")
            continue

    syst_uncs = {}
    for cutset in ry_cutsets_ref.keys():
        # Some pt bins don't have all trials and repeat the last one in the reference file,
        # while they are not generated in the multitrial as only single pt-bins are processed
        if ry_trials_results[cutset]['V2s'] == []:
            logger(f"No results for cutset {cutset}, skipping.", "WARNING")
            continue

        leg_header = f'Cutset {cutset}, {pt_min} < #it{{p}}_{{T}} < {pt_max} GeV/#it{{c}}'
        syst_uncs[cutset] = compute_systematics_cutset(f"{out_dir}/{mult_dir}", cutset, leg_header, \
                                   reference_results[cutset]['V2'], reference_results[cutset]['V2Unc'], \
                                   ry_trials_results[cutset]['V2s'], ry_trials_results[cutset]['V2sUncs'], \
                                   ry_trials_results[cutset]['Significances'], ry_trials_results[cutset]['SignificancesUncs'], \
                                   ry_trials_results[cutset]['Chi2s'], ry_trials_results[cutset]['Chi2sUncs'])

    # Prompt v2 systematic uncertainty wrt reference
    v2_prompt_ref = v2_vs_frac_ref.Get('hV2VsPtPrompt').GetBinContent(pt_bin_ref)
    v2_prompt_ref_unc = v2_vs_frac_ref.Get('hV2VsPtPrompt').GetBinError(pt_bin_ref)
    syst_uncs['PromptV2'] = compute_systematics_fin_val(f"{out_dir}/{mult_dir}/SystV2Prompt.root", f"Prompt #it{{v}}_{{2}}, {pt_min} < #it{{p}}_{{T}} < {pt_max} GeV/#it{{c}}", \
                                   v2_prompt_ref, v2_prompt_ref_unc, \
                                   v2_prompt_trials, v2_prompt_trials_uncs)

    # FD v2 systematic uncertainty wrt reference
    v2_non_prompt_ref = v2_vs_frac_ref.Get('hV2VsPtFD').GetBinContent(pt_bin_ref)
    v2_non_prompt_ref_unc = v2_vs_frac_ref.Get('hV2VsPtFD').GetBinError(pt_bin_ref)
    syst_uncs['NonPromptV2'] = compute_systematics_fin_val(f"{out_dir}/{mult_dir}/SystV2FD.root", f"Non-prompt #it{{v}}_{{2}}, {pt_min} < #it{{p}}_{{T}} < {pt_max} GeV/#it{{c}}", \
                                   v2_non_prompt_ref, v2_non_prompt_ref_unc, \
                                   v2_non_prompt_trials, v2_non_prompt_trials_uncs)

    return syst_uncs

def compute_systematics_fin_val(out_file_path, leg_header, vn_ref, vn_ref_unc, vn_trials, vn_trials_unc):

    h_vn_vs_trial = TH1F('h_vn_vs_trial', 'h_vn_vs_trial;trial;#it{v}_{n}', len(vn_trials), 0, len(vn_trials)+1)
    h_syst = TH1F('h_syst', 'h_syst;#it{v}_{n}(trial) - #it{v}_{n}(ref.);Counts', 200, -0.05, 0.05)
    canvas = TCanvas(f'c_syst_multitrial', f'c_syst_multitrial', 1600, 800)

    canvas.cd().Divide(2, 1)

    canvas.SetLeftMargin(0)
    canvas.SetRightMargin(0)
    canvas.SetTopMargin(0)
    canvas.SetBottomMargin(0)

    # Restrict the panels to fit the axis labels on both canvases
    canvas.cd(1).SetLeftMargin(0.16)
    canvas.cd(1).SetRightMargin(0.12)
    canvas.cd(1).SetTopMargin(0.05)
    canvas.cd(1).SetBottomMargin(0.12)
    canvas.cd(2).SetLeftMargin(0.16)
    canvas.cd(2).SetRightMargin(0.12)
    canvas.cd(2).SetTopMargin(0.05)
    canvas.cd(2).SetBottomMargin(0.12)

    SetObjectStyle(h_vn_vs_trial, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    SetObjectStyle(h_syst, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)

    for itrial, (vn, vn_unc) in enumerate(zip(vn_trials, vn_trials_unc)):    # loop over trials
        h_vn_vs_trial.SetBinContent(itrial+1, vn)
        h_vn_vs_trial.SetBinError(itrial+1, vn_unc)
        h_syst.Fill(vn-vn_ref)

    # Compute systematic uncertainty
    rms = h_syst.GetRMS()
    mean = h_syst.GetMean()
    syst_unc = np.sqrt(rms**2 + mean**2)
    max_syst = h_syst.GetMaximum()

    g_syst = TGraphAsymmErrors()
    g_syst.SetPoint(0, 0, max_syst*0.5)
    g_syst.SetPointError(0, syst_unc, syst_unc,
                            max_syst*0.5, max_syst*0.5)
    SetObjectStyle(g_syst, markerstyle=20, markercolor=kOrange+2,
                   markersize=1, linecolor=kOrange+2,
                      linewidth=2, fillcolor=kOrange+2, fillstyle=3153,
                      fillalpha=0.5, linestyle=9)

    # Pad 1: vn vs trial
    canvas.cd(1).SetGrid()
    # Define reference line
    g_ref = TGraphAsymmErrors()
    g_ref.SetPoint(0, 0, vn_ref)
    g_ref.SetPointError(0, 0, 0, vn_ref_unc, vn_ref_unc)
    g_ref.SetPoint(1, len(vn_trials), vn_ref)
    g_ref.SetPointError(1, 0, 0, vn_ref_unc, vn_ref_unc)
    SetObjectStyle(g_ref, markerstyle=20, markercolor=kAzure+2,
                   markersize=0, linecolor=kAzure+2,
                   linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
    h_vn_vs_trial.Draw('same')
    h_vn_vs_trial.GetXaxis().SetTitleOffset(1.05)
    if not "Non-prompt" in leg_header:
        h_vn_vs_trial.GetYaxis().SetTitleOffset(1.6)
    g_ref.Draw('c3 same')
    canvas.cd(2)

    # Legend
    leg = TLegend(0.25, 0.65, 0.85, 0.85)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetHeader(leg_header)
    # leg.AddEntry(hvn, f'vn vs trial ({hvn.GetEntries()} trials)', 'p')
    leg.AddEntry(g_syst, f'#sqrt{{shift^{{2}} + rms^{{2}}}} = {syst_unc:.3f}', 'f')
    leg.AddEntry(g_ref, 'vn stat. unc.', 'f')

    # Define reference vertical line at 1
    g_ref_two = g_ref.Clone()
    g_ref_two.SetPoint(0, 0, max_syst*0.5)
    g_ref_two.SetPointError(0, vn_ref_unc, vn_ref_unc, max_syst*0.5, max_syst*0.5)
    SetObjectStyle(g_ref_two, markerstyle=20, markercolor=kAzure+2,
                   markersize=1, linecolor=kAzure+2,
                   linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
    h_syst.GetYaxis().SetRangeUser(0, max_syst*1.8)
    h_syst.GetXaxis().SetRangeUser(h_syst.GetXaxis().GetXmin(), h_syst.GetXaxis().GetXmax())
    h_syst.GetXaxis().SetTitleOffset(1.05)
    # h_syst.GetYaxis().SetTitleOffset(0.6)
    h_syst.Draw('same')
    g_syst.Draw('2')
    g_ref_two.Draw('2')
    g_syst.Draw('2')
    leg.Draw()

    os.makedirs(os.path.dirname(out_file_path), exist_ok=True)
    canvas.SaveAs(f'{out_file_path.replace('.root', '.pdf')}')

    out_file = TFile(out_file_path, 'recreate')
    h_syst.Write()
    logger(f'Saved systematic results to {out_file_path}', "INFO")
    out_file.Close()

    return syst_unc

def compute_systematics_cutset(out_dir, suffix, leg_header, vn_ref, vn_ref_unc, vn_trials, vn_trials_unc, signif_trials, signif_trials_unc, chi2_trials, chi2_trials_unc):

    h_vn_vs_trial = TH1F('h_vn_vs_trial', 'h_vn_vs_trial;trial;#it{v}_{n}', len(vn_trials), 0, len(vn_trials)+1)
    h_signif_vs_trial = TH1F('h_signif_vs_trial', 'h_signif_vs_trial;trial;significance', len(vn_trials), 0, len(vn_trials)+1)
    h_chi2_vs_trial = TH1F('h_chi2_vs_trial', 'h_chi2_vs_trial;trial;#chi^{2}', len(vn_trials), 0, len(vn_trials)+1)
    h_syst = TH1F('h_syst', 'h_syst;#it{v}_{n}(trial) - #it{v}_{n}(ref.);Counts', 200, -0.05, 0.05)

    # Set fonts, sizes, ticks before drawing
    h_vn_vs_trial.GetXaxis().SetTitleSize(0.05)
    h_vn_vs_trial.GetXaxis().SetTitleOffset(0.85)
    h_signif_vs_trial.GetXaxis().SetTitleSize(0.05)
    h_signif_vs_trial.GetXaxis().SetTitleOffset(0.85)
    h_chi2_vs_trial.GetXaxis().SetTitleSize(0.05)
    h_chi2_vs_trial.GetXaxis().SetTitleOffset(0.85)
    h_syst.GetXaxis().SetTitleSize(0.05)
    h_syst.GetXaxis().SetTitleOffset(0.85)

    h_vn_vs_trial.GetYaxis().SetTitleSize(0.05)
    h_vn_vs_trial.GetYaxis().SetTitleOffset(1.01)
    h_signif_vs_trial.GetYaxis().SetTitleSize(0.05)
    h_signif_vs_trial.GetYaxis().SetTitleOffset(1.01)
    h_chi2_vs_trial.GetYaxis().SetTitleSize(0.05)
    h_chi2_vs_trial.GetYaxis().SetTitleOffset(1.01)
    h_syst.GetYaxis().SetTitleSize(0.05)
    h_syst.GetYaxis().SetTitleOffset(1.01)

    canvas = TCanvas(f'c_syst_multitrial', f'c_syst_multitrial', 800, 800)
    canvas.cd().Divide(2, 2)
    canvas.cd(1).SetLeftMargin(0.12)
    canvas.cd(1).SetTopMargin(0.12)
    canvas.cd(2).SetRightMargin(0.12)
    canvas.cd(2).SetTopMargin(0.12)
    canvas.cd(3).SetLeftMargin(0.12)
    canvas.cd(3).SetBottomMargin(0.12)
    canvas.cd(4).SetRightMargin(0.12)
    canvas.cd(4).SetBottomMargin(0.12)

    SetObjectStyle(h_vn_vs_trial, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    SetObjectStyle(h_signif_vs_trial, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    SetObjectStyle(h_chi2_vs_trial, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    SetObjectStyle(h_syst, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)

    for itrial, (vn, vn_unc, chi2, chi2_unc, signif, signif_unc) in \
        enumerate(zip(vn_trials, vn_trials_unc, chi2_trials, \
                      chi2_trials_unc, signif_trials, signif_trials_unc)):    # loop over trials

        h_chi2_vs_trial.SetBinContent(itrial+1, chi2)
        h_chi2_vs_trial.SetBinError(itrial+1, chi2_unc)
        h_signif_vs_trial.SetBinContent(itrial+1, signif)
        h_signif_vs_trial.SetBinError(itrial+1, signif_unc)
        h_vn_vs_trial.SetBinContent(itrial+1, vn)
        h_vn_vs_trial.SetBinError(itrial+1, vn_unc)
        h_syst.Fill(vn - vn_ref)

    # Compute systematic uncertainty
    rms = h_syst.GetRMS()
    mean = h_syst.GetMean()
    syst_unc = np.sqrt(rms**2 + mean**2)
    max_syst = h_syst.GetMaximum()
    g_syst = TGraphAsymmErrors()
    g_syst.SetPoint(0, 0, max_syst*0.5)
    g_syst.SetPointError(0, syst_unc, syst_unc,
                            max_syst*0.5, max_syst*0.5)
    SetObjectStyle(g_syst, markerstyle=20, markercolor=kOrange+2,
                   markersize=1, linecolor=kOrange+2,
                      linewidth=2, fillcolor=kOrange+2, fillstyle=3153,
                      fillalpha=0.5, linestyle=9)

    # Pad 1: vn vs trial
    canvas.cd(1).SetGrid()
    # Define reference line
    g_ref = TGraphAsymmErrors()
    g_ref.SetPoint(0, 0, vn_ref)
    g_ref.SetPointError(0, 0, 0, vn_ref_unc, vn_ref_unc)
    g_ref.SetPoint(1, len(vn_trials), vn_ref)
    g_ref.SetPointError(1, 0, 0, vn_ref_unc, vn_ref_unc)
    SetObjectStyle(g_ref, markerstyle=20, markercolor=kAzure+2,
                   markersize=0, linecolor=kAzure+2,
                   linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
    h_vn_vs_trial.Draw('same')
    g_ref.Draw('c3 same')
    canvas.cd(2)

    # Legend
    leg = TLegend(0.25, 0.65, 0.85, 0.85)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetHeader(leg_header)
    # leg.AddEntry(hvn, f'vn vs trial ({hvn.GetEntries()} trials)', 'p')
    leg.AddEntry(g_syst, f'#sqrt{{shift^{{2}} + rms^{{2}}}} = {syst_unc:.3f}', 'f')
    leg.AddEntry(g_ref, 'vn stat. unc.', 'f')
    # Define reference vertical line at 1
    g_ref_two = g_ref.Clone()
    g_ref_two.SetPoint(0, 0, max_syst*0.5)
    g_ref_two.SetPointError(0, vn_ref_unc, vn_ref_unc, max_syst*0.5, max_syst*0.5)
    SetObjectStyle(g_ref_two, markerstyle=20, markercolor=kAzure+2,
                   markersize=1, linecolor=kAzure+2,
                   linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
    h_syst.GetYaxis().SetRangeUser(0, max_syst*1.8)
    h_syst.GetXaxis().SetRangeUser(h_syst.GetXaxis().GetXmin(), h_syst.GetXaxis().GetXmax())
    h_syst.Draw('same')
    g_syst.Draw('2')
    g_ref_two.Draw('2')
    g_syst.Draw('2')
    leg.Draw()
    # Pad 3: chi2 vs trial
    canvas.cd(3)
    h_chi2_vs_trial.Draw('same')
    canvas.Update()
    # Pad 4: significance vs trial
    canvas.cd(4)
    h_signif_vs_trial.Draw('same')
    canvas.Update()

    os.makedirs(out_dir, exist_ok=True)
    canvas.SaveAs(f'{out_dir}/SystRy_{suffix}.pdf')

    out_file_name = os.path.join(out_dir, f'SystRy_{suffix}.root')
    out_file = TFile(out_file_name, 'recreate')
    h_syst.Write()
    logger(f'Saved systematic results to {out_file_name}', "INFO")
    out_file.Close()

    return syst_unc

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('default_cfg', metavar='text', default='path to the reference yaml config')
    parser.add_argument('results_dir', metavar='text', default='path to the directory containing the .root files from multitrial')
    args = parser.parse_args()

    # Open reference config
    with open(args.default_cfg, 'r') as CfgFlow:
        cfg_ref = yaml.safe_load(CfgFlow)

    # Retrieve results of reference cfg
    v2_vs_frac_ref = TFile.Open(f"{args.results_dir}/v2/v2VsFrac.root", "read")
    ry_files = [f for f in os.listdir(f"{args.results_dir}/raw_yields/") if re.match(r"raw_yields_\d+\.root$", f)]
    cutsets = sorted((re.search(r"raw_yields_(\d+)\.root", f).group(1) for f in ry_files), key=int)

    ry_cutsets_dict = {
        cutset: TFile.Open(os.path.join(f"{args.results_dir}/raw_yields/", f"raw_yields_{cutset}.root"), "READ")
        for cutset in cutsets
    }

    # List all pt-bins for which multitrial results are available
    multitrial_pt_dirs = [pt_dir for pt_dir in os.listdir(f"{args.results_dir}/syst/multitrial/") \
                          if os.path.isdir(os.path.join(args.results_dir, 'syst/multitrial', pt_dir)) \
                          and pt_dir.startswith('pt_')]

    syst_uncs_all = {}
    pt_bins = []
    for mult_dir in multitrial_pt_dirs:
        out_dir = os.path.join(args.results_dir, 'syst/multitrial', f"summary")
        syst_uncs_all[mult_dir] = compute_syst_multitrial(v2_vs_frac_ref, ry_cutsets_dict, mult_dir, out_dir)
        print(f"STARTING TO PLOT VARIATIONS FOR {mult_dir} ...")
        plot_variations_multitrial(cfg_ref, v2_vs_frac_ref, ry_cutsets_dict, mult_dir, out_dir, syst_uncs_all[mult_dir])
        _, pt_min_times_10, pt_max_times_10 = mult_dir.split('_')
        pt_min = float(pt_min_times_10) / 10.
        pt_max = float(pt_max_times_10) / 10.
        if pt_min not in pt_bins:
            pt_bins.append(pt_min)
        if pt_max not in pt_bins:
            pt_bins.append(pt_max)

    pt_bins = sorted(pt_bins)
    # Produce plots with all pt-bins summary for prompt and non-prompt v2
    h_prompt_syst = TH1F('h_prompt_syst', 'h_prompt_syst;#it{p}_{T} (GeV/#it{c});Fit syst. unc. on prompt #it{v}_{2}', len(pt_bins)-1, array.array('d', pt_bins))
    h_non_prompt_syst = TH1F('h_non_prompt_syst', 'h_non_prompt_syst;#it{p}_{T} (GeV/#it{c});Fit syst. unc. on non-prompt #it{v}_{2}', len(pt_bins)-1, array.array('d', pt_bins))
    for i_pt, mult_dir in enumerate(sorted(multitrial_pt_dirs)):
        logger(f"Processing final summary for {mult_dir} ...", "INFO")
        h_prompt_syst.SetBinContent(i_pt+1, syst_uncs_all[mult_dir]['PromptV2'])
        h_non_prompt_syst.SetBinContent(i_pt+1, syst_uncs_all[mult_dir]['NonPromptV2'])

    # Draw canvas for prompt
    canvas_prompt = TCanvas('c_prompt_syst', 'c_prompt_syst', 600, 600)
    canvas_prompt.SetTopMargin(0.12)
    canvas_prompt.SetLeftMargin(0.15)
    canvas_prompt.SetRightMargin(0.12)
    canvas_prompt.SetBottomMargin(0.15)
    SetObjectStyle(h_prompt_syst, markerstyle=20, markercolor=kOrange+2,
                   markersize=1.,linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.7)
    SetObjectStyle(h_prompt_syst, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    h_prompt_syst.Draw('hist')
    h_prompt_syst.GetXaxis().SetTitleOffset(1.07)
    canvas_prompt.SaveAs(f'{args.results_dir}/syst/multitrial/summary/SystV2Prompt_AllPtBins.pdf')

    # Draw canvas for non-prompt
    canvas_non_prompt = TCanvas('c_non_prompt_syst', 'c_non_prompt_syst', 600, 600)
    canvas_non_prompt.SetTopMargin(0.12)
    canvas_non_prompt.SetLeftMargin(0.15)
    canvas_non_prompt.SetRightMargin(0.12)
    canvas_non_prompt.SetBottomMargin(0.15)
    SetObjectStyle(h_non_prompt_syst, markerstyle=20, markercolor=kOrange+2,
                   markersize=1.,linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.7)
    SetObjectStyle(h_non_prompt_syst, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    h_non_prompt_syst.Draw('hist')
    h_non_prompt_syst.GetXaxis().SetTitleOffset(1.07)
    canvas_non_prompt.SaveAs(f'{args.results_dir}/syst/multitrial/summary/SystV2FD_AllPtBins.pdf')

    print("Combining all SystV2Prompt and SystV2FD plots into multi-page pdfs...")
    # Create a multi-page pdf with all SystV2Prompt and SystV2FD plots
    all_prompt_files = sorted(glob.glob(f'{args.results_dir}/syst/multitrial/summary/**/SystV2Prompt.pdf', recursive=True))
    all_non_prompt_files = sorted(glob.glob(f'{args.results_dir}/syst/multitrial/summary/**/SystV2FD.pdf', recursive=True))

    # --- Merge Prompt ---
    if all_prompt_files:
        merger = PdfMerger()
        for pdf in all_prompt_files:
            merger.append(pdf)
        out_prompt = f"{args.results_dir}/syst/multitrial/summary/TrialsSystV2Prompt.pdf"
        merger.write(out_prompt)
        merger.close()

    # --- Merge Non-prompt ---
    if all_non_prompt_files:
        merger = PdfMerger()
        for pdf in all_non_prompt_files:
            merger.append(pdf)
        out_fd = f"{args.results_dir}/syst/multitrial/summary/TrialsSystV2FD.pdf"
        merger.write(out_fd)
        merger.close()
