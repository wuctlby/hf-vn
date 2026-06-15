#!/usr/bin/env python3
"""
mass_fit.py — Fit invariant mass distributions from projection files
using flarefly's F2MassFitter.

Adapted for hf-vn-dev/dev-v0/ repository structure.

Takes projection ROOT file → fits hMassData per pT bin → outputs raw yields histograms.

Usage:
    CUDA_VISIBLE_DEVICES=-1 python3 mass_fit.py flow_config.yml proj_file.root -o output_dir
"""
import argparse
import os
os.environ['ZFIT_DISABLE_TF_WARNINGS'] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
import sys
import yaml
import numpy as np
import ROOT
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from pdf_merger import collage_pdf_pages_to_single, collage_pdfs_by_page

# ── Repo paths ──────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
from charmbulk_utils import get_pt_dependent_param

# Add utility modules from dev-v0
sys.path.insert(0, os.path.join(REPO_ROOT, "utils"))
from utils import logger

from flarefly import F2MassFitter
from flarefly.data_handler import DataHandler

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# Map from our naming to flarefly naming
SGN_MAP = {
    "kGaus": "gaussian",
    "k2Gaus": "doublegaus",
    "k2GausSigmaRatioPar": "doublegaus",
    "crystalball": "crystalball",
    "doublecb": "doublecb",
    "gausexptail": "gausexptail",
    "genergausexptail": "genergausexptail",
    "genergaussian": "genergaussian",
    "genercrystalball": "genercrystalball",
    "bifurgaus": "bifurgaus",
}

BKG_MAP = {
    "kExpo": "expo",
    "kLin": "chebpol",
    "kPol2": "chebpol",
    "kPol3": "chebpol",
}

PARTICLE_LABELS = {
    "D0": r"$\mathrm{D}^{0}$",
    "Dzero": r"$\mathrm{D}^{0}$",
    "Dplus": r"$\mathrm{D}^{+}$",
    "Ds": r"$\mathrm{D}_{s}^{+}$",
    "Lc": r"$\Lambda_{c}^{+}$",
}

RAW_PART_LABELS = {
    D: label.replace('$', '') for D, label in PARTICLE_LABELS.items()
}

MASS_LABELS = {
    "D0": r"$M(\mathrm{K^- \pi^+})\ \mathrm{(GeV/}c^2)$",
    "Dzero": r"$M(\mathrm{K^- \pi^+})\ \mathrm{(GeV/}c^2)$",
    "Dplus": r"$M(\mathrm{\pi^+ K^- \pi^+})\ \mathrm{(GeV/}c^2)$",
    "Ds": r"$M(\mathrm{K^+ K^- \pi^+})\ \mathrm{(GeV/}c^2)$",
    "Lc": r"$M(\mathrm{p K^+ \pi^-})\ \mathrm{(GeV/}c^2)$",
}

DECAY_CHANNELS = {
    "D0": [r"$\mathrm{D}^{0} \rightarrow \mathrm{K}^{-} \pi^{+}$"],
    "Dzero": [r"$\mathrm{D}^{0} \rightarrow \mathrm{K}^{-} \pi^{+}$"],
    "Dplus": [r"$\mathrm{D}^{+} \rightarrow \mathrm{K}^{-} \pi^{+} \pi^{+}$"],
    "Ds": [r"$\mathrm{D}_{s}^{+} \rightarrow \mathrm{K}^{+} \mathrm{K}^{-} \pi^{+}$"],
    "Lc": [r"$\Lambda_{c}^{+} \rightarrow \mathrm{p} \mathrm{K}^{-} \pi^{+}$"],
}


def get_bkg_config(bkg_str):
    """Return (bkg_name, bkg_pdf_arg: [str]) for given bkg_str"""
    if bkg_str == "kExpo":
        return "expo", ["expo"]
    elif bkg_str == "kLin":
        return "chebpol", ["chebpol1"]
    elif bkg_str == "kPol2":
        return "chebpol", ["chebpol2"]
    elif bkg_str == "kPol3":
        return "chebpol", ["chebpol3"]
    else:
        logger(f"Unknown bkg {bkg_str}, falling back to expo", "WARNING")
        return "expo", ["expo"]


def try_fit(h_mass, mass_min, mass_max, rebin, sgn_name, bkg_name,
            bkg_pdf_arg, sig0, fix_sigma=None, Dmeson="D0"):
    """Attempt a mass fit; returns (fitter, True) or (None, False).
    
    If fix_sigma is not None, the sigma parameter is fixed to that value
    via set_signal_initpar(..., fix=True)."""
    h = h_mass.Clone()
    h.SetDirectory(0)
    h.GetXaxis().SetRangeUser(mass_min, mass_max)
    if rebin > 1:
        h.Rebin(rebin)

    dh = DataHandler(h, var_name="mass", limits=[mass_min, mass_max])
    fitter = F2MassFitter(
        dh,
        name_signal_pdf=[sgn_name],
        name_background_pdf=bkg_pdf_arg,
        label_signal_pdf=DECAY_CHANNELS.get(Dmeson, ["Signal"]),
        label_bkg_pdf=["Background"],
    )

    # Set initial parameters
    fitter.set_signal_initpar(0, "mu", 1.86)
    if fix_sigma is not None:
        fitter.set_signal_initpar(0, "sigma", fix_sigma, fix=True)
    else:
        fitter.set_signal_initpar(0, "sigma", sig0)
    fitter.set_background_initpar(0, "lamb", -0.5)

    fitter.mass_zfit(do_prefit=False)
    return fitter, True

def run_mass_fit(fit_config, proj_file, batch=True):
    """Fit invariant mass distributions with following configuration:
        - D meson: ["Dmeson"] (e.g. "D0", "Dplus", "Ds", "Lc")
        - mass path: the mass histogram path (e.g. "a_side/pt_79_97/hMassData")
        - Signal function: ["SgnFunc"] (e.g. "kGaus")
        - Background function: ["fitConfig"]["BkgFunc"] (e.g. "kExpo")
        - Mass fit range: ["MassFitRanges"] (e.g. [[1.7, 2.05]] or [[1.7, 2.05], [1.6, 2.1], ...] per pt bin)
        - Initial sigma: ["Sigma"] (e.g. 0.02 or [0.02, 0.018, ...] per pt bin)
        - Rebin factor: ["fitConfig"]["Rebin"] (e.g. 1 or [1, 2, 4, ...] per pt bin)
        Outputs histograms of raw yields, mean, sigma, chi2/ndf, significance, and S/B vs pT to {output_dir}/raw_yields_XX.root
    """
    if batch:
        ROOT.gROOT.SetBatch(True)

    if isinstance(fit_config, dict):
        fit_cfg = fit_config
    else:
        with open(fit_config, "r") as f:
            fit_cfg = yaml.safe_load(f)

    output_subdir = fit_cfg["output_subdir"]
    ptbins = fit_cfg["ptbins"]
    ptmins = ptbins[:-1]
    ptmaxs = ptbins[1:]
    nPtBins = len(ptmins)
    Dmeson = fit_cfg["Dmeson"]

    # sgn_name = SGN_MAP.get(sgn_func, "gaussian")
    # bkg_name, bkg_pdf_arg = get_bkg_config(bkg_func)

    # Open projection file
    infile = ROOT.TFile.Open(proj_file)
    if not infile or not infile.IsOpen():
        logger(f"Cannot open {proj_file}", "ERROR")
        sys.exit(1)

    # Create output ROOT and PDF files
    proj_dir = os.path.dirname(proj_file)
    ry_basename = os.path.basename(proj_file).replace("proj_", "raw_yields_")
    ry_out = proj_dir.replace("projs", "raw_yields")
    ry_out_path = os.path.join(ry_out, f"{output_subdir}")
    os.makedirs(ry_out_path, exist_ok=True)
    out_path = os.path.join(ry_out_path, ry_basename)
    pdf_path = out_path.replace(".root", ".pdf")
    pdf = PdfPages(pdf_path)
    residual_pdf_path = os.path.join(os.path.dirname(pdf_path), os.path.basename(pdf_path).replace("raw_yields", "residuals"))
    residual_pdf = PdfPages(residual_pdf_path)

    # Create output histograms
    pt_arr = np.array(ptbins, dtype="d")
    hRawYields = ROOT.TH1F("hRawYields", ";p_{T} (GeV/c);raw yield", nPtBins, pt_arr)
    hMean = ROOT.TH1F("hMean", ";p_{T} (GeV/c);mean", nPtBins, pt_arr)
    hSigma = ROOT.TH1F("hSigma", ";p_{T} (GeV/c);#sigma", nPtBins, pt_arr)
    hChi2 = ROOT.TH1F("hRedChi2", ";p_{T} (GeV/c);#chi^{2}/ndf", nPtBins, pt_arr)
    hSignif = ROOT.TH1F("hSignificance", ";p_{T} (GeV/c);significance", nPtBins, pt_arr)
    hSoverB = ROOT.TH1F("hSoverB", ";p_{T} (GeV/c);S/B", nPtBins, pt_arr)

    for iPt, (pt_min, pt_max) in enumerate(zip(ptmins, ptmaxs)):
        pt_label = f"pt_{int(pt_min*10)}_{int(pt_max*10)}"
        logger(f"Fitting {pt_label}: {pt_min:.1f}–{pt_max:.1f} GeV/c")

        # fit parameters for this pt bin (use last value if list is shorter than nPtBins)
        sgn_func = fit_cfg["sgn_func"][min(iPt, len(fit_cfg["sgn_func"]) - 1)]
        bkg_func = fit_cfg["bkg_func"][min(iPt, len(fit_cfg["bkg_func"]) - 1)]
        sgn_name = SGN_MAP.get(sgn_func, "gaussian")
        bkg_name, bkg_pdf_arg = get_bkg_config(bkg_func)
        mass_min, mass_max = fit_cfg["mass_fit_range"][min(iPt, len(fit_cfg["mass_fit_range"]) - 1)]
        rebin = fit_cfg["rebin"][min(iPt, len(fit_cfg["rebin"]) - 1)]
        sigma_init = fit_cfg["sigma_init"][min(iPt, len(fit_cfg["sigma_init"]) - 1)]
        fix_sigma = fit_cfg["fix_sigma"][min(iPt, len(fit_cfg["fix_sigma"]) - 1)]
        fix_val = sigma_init if fix_sigma else None

        # Get mass histogram
        h_mass = infile.Get(f"{pt_label}/{fit_cfg['mass_path']}")
        if not h_mass:
            logger(f"  No hMassData for {pt_label}, skipping", "WARNING")
            continue
        h_mass.SetDirectory(0)

        # Fit
        try:
            fitter, ok = try_fit(h_mass, mass_min, mass_max, rebin,
                                 sgn_name, bkg_name, bkg_pdf_arg, sigma_init, fix_sigma=fix_val,
                                 Dmeson=Dmeson)
            if not ok:
                continue
        except Exception as e:
            logger(f"  Fit failed: {e}", "WARNING")
            continue

        # Extract results
        try:
            # Raw yield: tuple (value, error)
            ry_raw = fitter.get_raw_yield(0)
            ry, ry_err = float(ry_raw[0]), float(ry_raw[1])

            # Signal pars: list of dicts
            sgn_pars = fitter.get_signal_pars()[0]
            sgn_pars_uncs = fitter.get_signal_pars_uncs()[0]

            mean  = float(sgn_pars["mu"])
            mean_unc  = float(sgn_pars_uncs["mu"])
            sigma = float(sgn_pars["sigma"])
            sigma_unc = float(sgn_pars_uncs["sigma"])
            chi2      = fitter.get_chi2()
            ndf       = fitter.get_ndf()
            chi2_ndf  = float(chi2) / float(ndf)
            bkg, bkg_err = fitter.get_background(0, nsigma=3)  # background in ±3σ
            bkg, bkg_err = float(bkg), float(bkg_err)

            # Significance and S/B: tuples (value, error)
            signif, signif_err = float(fitter.get_significance(0)[0]), float(fitter.get_significance(0)[1])
            sob, sob_err = float(fitter.get_signal_over_background(0)[0]), float(fitter.get_signal_over_background(0)[1])

        except Exception as e:
            logger(f"  Error extracting results: {e}", "WARNING")
            continue

        bin_idx = iPt + 1
        hRawYields.SetBinContent(bin_idx, ry)
        hRawYields.SetBinError(bin_idx, ry_err)
        hMean.SetBinContent(bin_idx, mean)
        hMean.SetBinError(bin_idx, mean_unc)
        hSigma.SetBinContent(bin_idx, sigma)
        hSigma.SetBinError(bin_idx, sigma_unc)
        hChi2.SetBinContent(bin_idx, chi2_ndf)
        hSignif.SetBinContent(bin_idx, signif)
        hSignif.SetBinError(bin_idx, signif_err)
        hSoverB.SetBinContent(bin_idx, sob)
        hSoverB.SetBinError(bin_idx, sob_err)

        # Draw fit result and save to PDF
        try:
            # ----- mass fit plot -----
            fig, axs = fitter.plot_mass_fit(
                show_extra_info=False,
                figsize=(10, 8),
                style='ATLAS',
                legend_loc='upper right',
                extra_info_loc=['lower right', 'upper left'], # chi2/ndf in lower right, significance and S/B in upper left
                extra_info_massnsigma=3, # the n sigma range to calculate S, B, S/B, significance in the plot's extra info box
                axis_title= MASS_LABELS.get(Dmeson, r"$\mathrm{(GeV/}c^2)$"),
            )
            axs.text(
                0.04, 0.97,
                fr"${pt_min:.1f} < p_{{\mathrm{{T}}}}^{{{RAW_PART_LABELS.get(Dmeson)}}} < {pt_max:.1f} \ \mathrm{{GeV/}}c$",
                transform=axs.transAxes,
                ha='left', va='top',
                fontsize=16,
                # bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.95)
            )
            fit_info_text = (
                    fr"$\mu = {mean*1000:.1f} \pm {mean_unc*1000:.1f}$ MeV/$c^2$" "\n"
                    fr"$\sigma = {sigma*1000:.1f} \pm {sigma_unc*1000:.1f}$ MeV/$c^2$" "\n"
                    fr"$S = {ry:.0f} \pm {ry_err:.0f}$" "\n"
                    fr"$B(3\sigma) = {bkg:.0f} \pm {bkg_err:.0f}$" "\n"
                    fr"$S/B(3\sigma) = {sob:.2f} \pm {sob_err:.2f}$" "\n"
                    fr"Signif.$(3\sigma) = {signif:.1f} \pm {signif_err:.1f}$"
                )
            axs.text(
                0.04, 0.92, fit_info_text,
                transform=axs.transAxes, ha='left', va='top',
                fontsize=16,
            )
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            
            # ----- residual plot -----
            # plot_std_residuals: (data - total_fit) / sqrt(data)
            # plot_raw_residuals: data - fitted background + fitted signal
            fig_residual, ax_residual = fitter.plot_raw_residuals(
                figsize=(10, 8),
                style='ATLAS',
            )
            axs.text(
                0.04, 0.97,
                fr"${pt_min:.1f} < p_{{\mathrm{{T}}}}^{{{RAW_PART_LABELS.get(Dmeson)}}} < {pt_max:.1f} \ \mathrm{{GeV/}}c$",
                transform=axs.transAxes,
                ha='left', va='top',
                fontsize=16,
                # bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.95)
            )
            fig_residual.tight_layout()
            residual_pdf.savefig(fig_residual)
            plt.close(fig_residual)
        except Exception as e:
            logger(f"  Plot failed for {pt_label}/{output_subdir}: {e}", "WARNING")

    infile.Close()
    pdf.close()
    residual_pdf.close()
    infile.Close()

    # Save output
    outfile = ROOT.TFile(out_path, "RECREATE")

    hRawYields.Write()
    hMean.Write()
    hSigma.Write()
    hChi2.Write()
    hSignif.Write()
    hSoverB.Write()

    outfile.Close()
    logger(f"Results saved to {out_path}")
    import pathlib
    pdf_path_obj = pathlib.Path(pdf_path)
    output_pdf = pdf_path_obj.parent / "all_pt_bins" / pdf_path_obj.name
    output_png = output_pdf.with_suffix(".png")
    merged_pdf_path, merged_png_path = collage_pdf_pages_to_single(pdf_path, str(output_pdf), output_png=str(output_png), render_zoom=3.0)

    residual_pdf_path_obj = pathlib.Path(residual_pdf_path)
    output_residual_pdf = residual_pdf_path_obj.parent / "all_pt_bins" / residual_pdf_path_obj.name
    output_residual_png = output_residual_pdf.with_suffix(".png")
    merged_residual_pdf_path, merged_residual_png_path = collage_pdf_pages_to_single(residual_pdf_path, str(output_residual_pdf), output_png=str(output_residual_png), render_zoom=3.0)

    return out_path
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mass fit from projection files")
    parser.add_argument("config", help="Charm-bulk configuration YAML")
    parser.add_argument("proj_file", help="Projection ROOT file")
    parser.add_argument("-b", "--batch", action="store_true", default=True, help="Batch mode")
    parser.add_argument("--load-sigma", action="store_true", help="Load sigma values from first fit instead of config")
    args = parser.parse_args()
    output_dir = os.path.join(os.path.dirname(args.proj_file), "raw_yields")
    
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)
    Dmeson = config["Dmeson"]
    ptbins = config["ptbins"]
    n_ptbins = len(ptbins) - 1
    mean_ptbins = config["projections"]["proj_data"]["MeanPtBins"]

    # build the fit config
    fit_config = config['fitConfig']
    task_config = []
    sgn_funcs       = get_pt_dependent_param(fit_config["SgnFunc"], n_ptbins)
    bkg_funcs       = get_pt_dependent_param(fit_config["BkgFunc"], n_ptbins)
    mass_fit_ranges = get_pt_dependent_param(fit_config["MassFitRanges"], n_ptbins, isList=True)
    sigma_inits     = get_pt_dependent_param(fit_config["Sigma"], n_ptbins)
    rebins          = get_pt_dependent_param(fit_config["Rebin"], n_ptbins)
    fix_sigmas      = get_pt_dependent_param(False, n_ptbins)

    out_paths = {}
    for iMeanPt, (mean_pt_min, mean_pt_max) in enumerate(zip(mean_ptbins[:-1], mean_ptbins[1:])):
        mean_pt_label = f"pt_{int(mean_pt_min*100)}_{int(mean_pt_max*100)}"
        for side in ["a_side", "b_side"]:
            if args.load_sigma:
                first_cutset_path = os.path.join(os.path.dirname(os.path.dirname(args.proj_file)), "raw_yields", side, mean_pt_label, "raw_yields_00.root")
                with ROOT.TFile.Open(first_cutset_path) as f:
                    if f and f.IsOpen() and f.Get("hSigma"):
                        sigma_inits = [f.Get("hSigma").GetBinContent(i+1) for i in range(n_ptbins)]
                        logger(f"Loaded sigma values from {first_cutset_path}: {sigma_inits}")
                        fix_sigmas = get_pt_dependent_param(True, n_ptbins)  # fix sigma for all pt bins
                    else:
                        logger(f"Cannot load sigma values from {first_cutset_path}, file or hSigma histogram not found. Falling back to config values.", "WARNING")
            task_cfg = {
                "Dmeson":           Dmeson,
                "output_subdir":    f"{side}/{mean_pt_label}",
                "mass_path":        f"{side}/{mean_pt_label}/hMassData",
                "sgn_func":         sgn_funcs,
                "bkg_func":         bkg_funcs,
                "mass_fit_range":   mass_fit_ranges,
                "sigma_init":       sigma_inits,
                "rebin":            rebins,
                "ptbins":           ptbins,
                "fix_sigma":        fix_sigmas,
            }
            out_path = run_mass_fit(task_cfg, args.proj_file, batch=args.batch)
            out_paths[f"{side}/{mean_pt_label}"] = out_path


    # Run fits
    #TODO: process all projection files in parallel, take the sigma from the results of the first cut
    # sigma_results = []
    # if args.extract_sigma:
    #     for task in task_config:
