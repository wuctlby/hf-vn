#!/usr/bin/env python3
"""Rewrite of FitCorrel.C that runs the DhCorrelationFitter workflow via PyROOT and a YAML config."""

import argparse
import math
import os
import shutil
from array import array
from ctypes import c_double, c_int

import numpy as np
import ROOT
import yaml
from ROOT import gSystem

script_dir = os.path.dirname(os.path.realpath(__file__))
gSystem.CompileMacro(os.path.join(script_dir, "DhCorrelationFitter.cxx"), "kO")
from ROOT import DhCorrelationFitter


def set_canvas_style():
    style = ROOT.gStyle
    style.SetOptStat(0)
    style.SetPadLeftMargin(0.2)
    style.SetPadRightMargin(0.005)
    style.SetPadBottomMargin(0.2)
    style.SetFrameLineWidth(2)
    style.SetLineWidth(2)
    style.SetCanvasDefH(1126)
    style.SetCanvasDefW(1840)


def set_th1_style(
    histo,
    h_title,
    h_xaxis_title,
    h_yaxis_title,
    marker_style,
    marker_color,
    marker_size,
    line_color,
    line_width,
    h_title_xaxis_offset=1.3,
    h_title_yaxis_offset=1.3,
    h_title_xaxis_size=0.045,
    h_title_yaxis_size=0.045,
    h_label_xaxis_size=0.045,
    h_label_yaxis_size=0.045,
    center_xaxis_title=False,
    center_yaxis_title=False,
):
    histo.SetTitle(h_title)
    histo.GetXaxis().SetTitle(h_xaxis_title)
    histo.GetYaxis().SetTitle(h_yaxis_title)
    histo.SetMarkerStyle(marker_style)
    histo.SetMarkerColor(marker_color)
    histo.SetMarkerSize(marker_size)
    histo.SetLineColor(line_color)
    histo.SetLineWidth(line_width)
    histo.GetXaxis().SetTitleOffset(h_title_xaxis_offset)
    histo.GetYaxis().SetTitleOffset(h_title_yaxis_offset)
    histo.GetXaxis().SetTitleSize(h_title_xaxis_size)
    histo.GetYaxis().SetTitleSize(h_title_yaxis_size)
    histo.GetXaxis().SetLabelSize(h_label_xaxis_size)
    histo.GetYaxis().SetLabelSize(h_label_yaxis_size)
    histo.GetXaxis().CenterTitle(center_xaxis_title)
    histo.GetYaxis().CenterTitle(center_yaxis_title)


def ensure_clean_dir(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


def format_bin_dir(prefix, low, high, scale):
    return f"{prefix}_{int(low * scale):.0f}_{int(high * scale):.0f}"


def ensure_per_bin(value, n_bins, name, default):
    if value is None:
        return [default] * n_bins
    if isinstance(value, list):
        if len(value) < n_bins:
            raise ValueError(f"{name} must provide at least {n_bins} entries, got {len(value)}")
        return value[:n_bins]
    return [value] * n_bins


def fit_correl(cfg_path):
    ROOT.gROOT.SetBatch(True)
    set_canvas_style()

    with open(cfg_path, "r") as cfg_file:
        config = yaml.safe_load(cfg_file)

    suffix = config.get("suffix")
    outdir = config.get("outdir")
    if not suffix or not outdir:
        raise ValueError("Config must define both 'suffix' and 'outdir'.")

    code_name = config.get("CodeName") or suffix
    out_root_dir = f"Output_CorrelationFitting_{code_name}_Root"
    out_png_dir = f"Output_CorrelationFitting_{code_name}_png"
    ensure_clean_dir(out_root_dir)
    ensure_clean_dir(out_png_dir)

    outdir_full = os.path.join(outdir, f"CorrelExtract_{suffix}")
    input_file = os.path.join(outdir_full, "CorrelationsResults", "CorrelationsResults.root")
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Correlation results file not found: {input_file}")

    pt_bins_cand = config.get("ptBinsCand")
    pt_bins_had = config.get("ptBinsHad")
    inv_mass_bins = config.get("invMassBins")
    if not (pt_bins_cand and pt_bins_had and inv_mass_bins):
        raise ValueError("Config must provide 'ptBinsCand', 'ptBinsHad' and 'invMassBins'.")

    pt_bins_cand = [float(x) for x in pt_bins_cand]
    pt_bins_had = [float(x) for x in pt_bins_had]
    inv_mass_bins = [[float(edge) for edge in edges] for edges in inv_mass_bins]

    n_bins_pt_cand = len(pt_bins_cand) - 1
    n_bins_pt_had = len(pt_bins_had) - 1
    if len(inv_mass_bins) != n_bins_pt_cand:
        raise ValueError("Number of inv. mass bin definitions must match number of candidate pt bins.")

    n_bins_inv_mass = len(inv_mass_bins[0]) - 1
    if n_bins_inv_mass <= 0:
        raise ValueError("Each inverse mass list must contain at least two entries.")
    if any(len(edges) != len(inv_mass_bins[0]) for edges in inv_mass_bins):
        raise ValueError("All inv mass configurations must hold the same number of edges.")

    bins_pt_cand_array = array("d", pt_bins_cand)

    method = config.get("method", "MassBinning").strip()
    method_lower = method.lower()
    if method_lower not in ("massbinning", "deltaphibinning"):
        raise ValueError("method must be either 'MassBinning' or 'DeltaPhiBinning'.")

    delta_phi_bins = config.get("deltaPhiBins")
    if method_lower == "deltaphibinning":
        if not delta_phi_bins:
            delta_phi_bins = np.linspace(-1.571, 4.712, 24).tolist()
        delta_phi_bins = [float(x) for x in delta_phi_bins]
        # Extraction writes only the first delta phi slice into CorrPhiDs_FinalPlots
        delta_phi_dir = format_bin_dir("DeltaPhiBin", delta_phi_bins[0], delta_phi_bins[1], 1000)
    else:
        delta_phi_dir = None

    fit_functions = config.get("FitFunction")
    fit_functions = fit_functions or config.get("fitConfig", {}).get("FitFunction")
    fit_functions = ensure_per_bin(fit_functions, n_bins_pt_cand, "FitFunction", 1)

    fix_baseline = config.get("FixBaseline", 0)
    fix_mean = config.get("FixMean", 0)
    n_baseline_points = config.get("nBaselinePoints", 0)
    points_for_baseline = config.get("binsForBaseline", [])
    if n_baseline_points and len(points_for_baseline) != n_baseline_points:
        raise ValueError("'binsForBaseline' length must match 'nBaselinePoints'.")
    points_for_baseline = [int(x) for x in points_for_baseline]

    par_vals = config.get("parVals", [])
    par_low_bounds = config.get("parLowBounds", [])
    par_upper_bounds = config.get("parUpperBounds", [])
    if par_vals and not (len(par_vals) == len(par_low_bounds) == len(par_upper_bounds)):
        raise ValueError("parVals, parLowBounds and parUpperBounds must share the same length.")
    par_vals = [float(x) for x in par_vals]
    par_low_bounds = [float(x) for x in par_low_bounds]
    par_upper_bounds = [float(x) for x in par_upper_bounds]

    is_reflected = config.get("IsReflected")
    if is_reflected is None:
        is_reflected = config.get("IsRiflected", False)
    shift_base_up = config.get("ShiftBaseUp", False)
    shift_base_down = config.get("ShiftBaseDown", False)

    do_correlation = config.get("doCorrelation", False)
    remove_ns_peak_low_pt = config.get("removeNSPeakLowPt", False)

    print("===========================")
    print("Input variables from config")
    for idx, func in enumerate(fit_functions):
        print(f"iPt = {idx + 1}  FitFunction = {func}")
    print(f"FixBaseline = {fix_baseline}")
    print(f"FixMean = {fix_mean}")
    print("===========================\n")

    f_min = -0.5 * math.pi
    f_max = 1.5 * math.pi

    pt_cand_pairs = list(zip(pt_bins_cand[:-1], pt_bins_cand[1:]))
    pt_had_pairs = list(zip(pt_bins_had[:-1], pt_bins_had[1:]))
    inv_mass_edges = inv_mass_bins[0]

    in_file = ROOT.TFile.Open(input_file)
    if not in_file or in_file.IsZombie():
        raise FileNotFoundError(f"Could not open {input_file}")

    h_baseline = [[None] * n_bins_inv_mass for _ in range(n_bins_pt_had)]
    h_ns_yield = [[None] * n_bins_inv_mass for _ in range(n_bins_pt_had)]
    h_ns_sigma = [[None] * n_bins_inv_mass for _ in range(n_bins_pt_had)]
    h_as_yield = [[None] * n_bins_inv_mass for _ in range(n_bins_pt_had)]
    h_as_sigma = [[None] * n_bins_inv_mass for _ in range(n_bins_pt_had)]
    h_beta = [[None] * n_bins_inv_mass for _ in range(n_bins_pt_had)]
    h_ns_bin_count = [[None] * n_bins_inv_mass for _ in range(n_bins_pt_had)]
    h_as_bin_count = [[None] * n_bins_inv_mass for _ in range(n_bins_pt_had)]
    h_v2_delta = [[None] * n_bins_inv_mass for _ in range(n_bins_pt_had)]

    for i_pt_had in range(n_bins_pt_had):
        for i_mass in range(n_bins_inv_mass):
            tag = f"PtBinAssoc{i_pt_had + 1}_InvMassBin{i_mass + 1}"
            edges = array("d", pt_bins_cand)
            h_baseline[i_pt_had][i_mass] = ROOT.TH1D(f"hBaselin_{tag}", "", n_bins_pt_cand, edges)
            h_ns_yield[i_pt_had][i_mass] = ROOT.TH1D(f"hNSYield_{tag}", "", n_bins_pt_cand, array("d", pt_bins_cand))
            h_ns_sigma[i_pt_had][i_mass] = ROOT.TH1D(f"hNSSigma_{tag}", "", n_bins_pt_cand, array("d", pt_bins_cand))
            h_as_yield[i_pt_had][i_mass] = ROOT.TH1D(f"hASYield_{tag}", "", n_bins_pt_cand, array("d", pt_bins_cand))
            h_as_sigma[i_pt_had][i_mass] = ROOT.TH1D(f"hASSigma_{tag}", "", n_bins_pt_cand, array("d", pt_bins_cand))
            h_beta[i_pt_had][i_mass] = ROOT.TH1D(f"hBeta_{tag}", "", n_bins_pt_cand, array("d", pt_bins_cand))
            h_ns_bin_count[i_pt_had][i_mass] = ROOT.TH1D(f"hNSYieldBinCount_{tag}", "", n_bins_pt_cand, array("d", pt_bins_cand))
            h_as_bin_count[i_pt_had][i_mass] = ROOT.TH1D(f"hASYieldBinCount_{tag}", "", n_bins_pt_cand, array("d", pt_bins_cand))
            h_v2_delta[i_pt_had][i_mass] = ROOT.TH1D(f"hv2Delta_{tag}", "", n_bins_pt_cand, array("d", pt_bins_cand))

    final_hist_file = ROOT.TFile(os.path.join(out_root_dir, f"CorrPhiDs_FinalPlots.root"), "RECREATE")

    for i_mass in range(n_bins_inv_mass):
        inv_min = inv_mass_edges[i_mass]
        inv_max = inv_mass_edges[i_mass + 1]
        print(f"[INFO] InvMass: {inv_min} - {inv_max}")
        for i_pt_had, (pt_had_min, pt_had_max) in enumerate(pt_had_pairs):
            print(f"[INFO] PtHad: {pt_had_min} - {pt_had_max}")
            canvas_name = f"CanvasCorrPhi_PtBinAssoc{i_pt_had + 1}_InvMassBin{i_mass + 1}"
            canvas_title = f"CorrPhiDs_PtBinAssoc{i_pt_had + 1}_InvMassBin{i_mass + 1}"
            canvas = ROOT.TCanvas(canvas_name, canvas_title)
            if n_bins_pt_cand <= 4:
                canvas.Divide(2, 2)
            elif n_bins_pt_cand <= 6:
                canvas.Divide(3, 2)
            else:
                pads_x = min(4, n_bins_pt_cand)
                pads_y = math.ceil(n_bins_pt_cand / pads_x)
                canvas.Divide(pads_x, pads_y)

            for i_pt_cand, (pt_cand_min, pt_cand_max) in enumerate(pt_cand_pairs):
                print(f"[INFO] PtCand: {pt_cand_min} - {pt_cand_max}")
                pt_cand_dir = format_bin_dir("PtCandBin", pt_cand_min, pt_cand_max, 10)
                pt_had_dir = format_bin_dir("PtHadBin", pt_had_min, pt_had_max, 10)
                if method_lower == "massbinning":
                    mass_dir = format_bin_dir("InvMassBin", inv_min, inv_max, 1000)
                    hist_path = f"{pt_cand_dir}/{pt_had_dir}/{mass_dir}/hCorrectedCorrHisto"
                else:
                    hist_path = f"{pt_cand_dir}/{pt_had_dir}/{delta_phi_dir}/hCorrectedCorrHisto"

                h_corr_phi = in_file.Get(hist_path)
                if not h_corr_phi:
                    raise RuntimeError(f"Missing histogram {hist_path}")

                canvas.cd(i_pt_cand + 1)
                canvas.SetTickx()
                canvas.SetTicky()

                set_th1_style(
                    h_corr_phi,
                    "",
                    "#Delta#phi [rad]",
                    "#frac{dN^{assoc}}{d#Delta#phi} [rad^{-1}]",
                    ROOT.kFullCircle,
                    ROOT.kRed + 1,
                    1.4,
                    ROOT.kRed + 1,
                    3,
                )
                h_corr_phi.SetStats(0)
                h_corr_phi.SetMinimum(0)

                corr_fitter = DhCorrelationFitter(h_corr_phi, f_min, f_max)
                corr_fitter.SetHistoIsReflected(False)
                corr_fitter.SetFixBaseline(fix_baseline)
                corr_fitter.SetBaselineUpOrDown(bool(shift_base_up), bool(shift_base_down))
                if n_baseline_points > 0:
                    baseline_array = (c_int * n_baseline_points)(*points_for_baseline)
                    corr_fitter.SetPointsForBaseline(n_baseline_points, baseline_array)
                corr_fitter.SetReflectedCorrHisto(bool(is_reflected))
                corr_fitter.SetFixMean(fix_mean)
                corr_fitter.SetPtRanges(pt_cand_min, pt_cand_max, pt_had_min, pt_had_max)
                if par_vals:
                    npars = len(par_vals)
                    vals = (c_double * npars)(*par_vals)
                    lows = (c_double * npars)(*par_low_bounds)
                    ups = (c_double * npars)(*par_upper_bounds)
                    corr_fitter.SetExternalValsAndBounds(npars, vals, lows, ups)
                corr_fitter.SetFuncType(int(fit_functions[i_pt_cand]))
                corr_fitter.Fitting(True, True)

                pttext = ROOT.TPaveText(0.15, 0.9, 0.85, 0.95, "NDC")
                pttext.SetFillStyle(0)
                pttext.SetBorderSize(0)
                pttext.AddText(
                    0.0,
                    0.8,
                    f"{pt_cand_min:.0f} < p_{{T}}^{{D_{{s}}}} < {pt_cand_max:.0f} GeV/c, p_{{T}}^{{assoc}} > {pt_had_min:.1f} GeV/c",
                )

                bin_index = i_pt_cand + 1
                if do_correlation:
                    h_baseline[i_pt_had][i_mass].SetBinContent(bin_index, corr_fitter.GetPedestal())
                    h_baseline[i_pt_had][i_mass].SetBinError(bin_index, corr_fitter.GetPedestalError())

                    if i_pt_cand == 0 and remove_ns_peak_low_pt:
                        for hist in (h_ns_yield, h_ns_sigma, h_beta):
                            hist[i_pt_had][i_mass].SetBinContent(bin_index, -1)
                            hist[i_pt_had][i_mass].SetBinError(bin_index, 0)
                    else:
                        h_ns_yield[i_pt_had][i_mass].SetBinContent(bin_index, corr_fitter.GetNSYield())
                        h_ns_yield[i_pt_had][i_mass].SetBinError(bin_index, corr_fitter.GetNSYieldError())
                        if fit_functions[i_pt_cand] not in (5, 6):
                            h_ns_sigma[i_pt_had][i_mass].SetBinContent(bin_index, corr_fitter.GetNSSigma())
                            h_ns_sigma[i_pt_had][i_mass].SetBinError(bin_index, corr_fitter.GetNSSigmaError())
                        else:
                            ns_sigma = math.sqrt(1.0 / corr_fitter.GetNSSigma())
                            err_rel = corr_fitter.GetNSSigmaError() / corr_fitter.GetNSSigma() / 2.0
                            h_ns_sigma[i_pt_had][i_mass].SetBinContent(bin_index, ns_sigma)
                            h_ns_sigma[i_pt_had][i_mass].SetBinError(bin_index, err_rel * ns_sigma)
                        h_as_yield[i_pt_had][i_mass].SetBinContent(bin_index, corr_fitter.GetASYield())
                        h_as_yield[i_pt_had][i_mass].SetBinError(bin_index, corr_fitter.GetASYieldError())
                        h_as_sigma[i_pt_had][i_mass].SetBinContent(bin_index, corr_fitter.GetASSigma())
                        h_as_sigma[i_pt_had][i_mass].SetBinError(bin_index, corr_fitter.GetASSigmaError())
                        if fit_functions[i_pt_cand] == 4:
                            h_beta[i_pt_had][i_mass].SetBinContent(bin_index, corr_fitter.GetBeta())
                            h_beta[i_pt_had][i_mass].SetBinError(bin_index, corr_fitter.GetBetaError())

                    h_ns_bin_count[i_pt_had][i_mass].SetBinContent(bin_index, corr_fitter.GetBinCountingNSYield())
                    h_ns_bin_count[i_pt_had][i_mass].SetBinError(bin_index, corr_fitter.GetBinCountingNSYieldErr())
                    h_as_bin_count[i_pt_had][i_mass].SetBinContent(bin_index, corr_fitter.GetBinCountingASYield())
                    h_as_bin_count[i_pt_had][i_mass].SetBinError(bin_index, corr_fitter.GetBinCountingASYieldErr())
                else:
                    h_v2_delta[i_pt_had][i_mass].SetBinContent(bin_index, corr_fitter.Getv2Delta())
                    h_v2_delta[i_pt_had][i_mass].SetBinError(bin_index, corr_fitter.Getv2DeltaError())

                h_corr_phi.Draw("same")
                pttext.Draw("same")

            png_name = os.path.join(out_png_dir, f"CorrPhiDs_PtBinAssoc{i_pt_had + 1}_InvMassBin{i_mass + 1}.png")
            root_name = os.path.join(out_root_dir, f"CorrPhiDs_PtBinAssoc{i_pt_had + 1}_InvMassBin{i_mass + 1}.root")
            canvas.SaveAs(png_name)
            canvas.SaveAs(root_name)
            canvas.Close()

    final_hist_file.cd()
    for i_mass in range(n_bins_inv_mass):
        for i_pt_had in range(n_bins_pt_had):
            if do_correlation:
                h_baseline[i_pt_had][i_mass].Write()
                h_ns_yield[i_pt_had][i_mass].Write()
                h_ns_sigma[i_pt_had][i_mass].Write()
                h_as_yield[i_pt_had][i_mass].Write()
                h_as_sigma[i_pt_had][i_mass].Write()
                h_beta[i_pt_had][i_mass].Write()
                h_ns_bin_count[i_pt_had][i_mass].Write()
                h_as_bin_count[i_pt_had][i_mass].Write()
            else:
                h_v2_delta[i_pt_had][i_mass].Write()
    final_hist_file.Close()
    in_file.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fit Dh azimuthal correlations")
    parser.add_argument(
        "config",
        nargs="?",
        default="config_CorrAnalysis_v2_010_negDeta.yaml",
        help="Path to the YAML configuration file",
    )
    args = parser.parse_args()
    fit_correl(args.config)
