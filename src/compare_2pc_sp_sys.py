#!/usr/bin/env python3
"""
Compare 2pc prompt v2 (with systematic uncertainties) against Biao's sp result.

Reads the output of compute_prompt_v2_unfold_sys.py and plots:
  - 2pc central values with statistical error bars (line + markers)
  - systematic uncertainty band (shaded region)
  - total systematic box
  - sp result for comparison

Adapted from compare.ipynb (cells 1-2).
"""

import os
import ROOT

ROOT.gROOT.SetBatch(True)


def read_hist_bins(h):
    """Return (centers, contents, errors, half_widths) from a TH1."""
    n = h.GetNbinsX()
    xs, ys, eys, exs = [], [], [], []
    for i in range(1, n + 1):
        xs.append(h.GetBinCenter(i))
        exs.append(0.5 * h.GetBinWidth(i))
        ys.append(h.GetBinContent(i))
        eys.append(h.GetBinError(i))
    return xs, ys, eys, exs


import ctypes


def read_graph_points(g):
    """Return (xs, ys, exs, eys) from a TGraph*."""
    n = g.GetN()
    xs, ys, exs, eys = [], [], [], []
    for i in range(n):
        x = ctypes.c_double(0.0)
        y = ctypes.c_double(0.0)
        g.GetPoint(i, x, y)
        xs.append(float(x.value))
        ys.append(float(y.value))
        exs.append(g.GetErrorX(i))
        eys.append(g.GetErrorY(i))
    return xs, ys, eys, exs


def read_graph_asymm_points(g):
    """Return (xs, ys, exs, eylows, eyhighs) from a TGraphAsymmErrors."""
    n = g.GetN()
    xs, ys, exs, eylows, eyhighs = [], [], [], [], []
    for i in range(n):
        x = ctypes.c_double(0.0)
        y = ctypes.c_double(0.0)
        g.GetPoint(i, x, y)
        xs.append(float(x.value))
        ys.append(float(y.value))
        exs.append(g.GetErrorX(i))
        eylows.append(g.GetErrorYlow(i))
        eyhighs.append(g.GetErrorYhigh(i))
    return xs, ys, exs, eylows, eyhighs


def draw_sys_band(xs, central_ys, sys_lows, sys_highs, exs, color=ROOT.kOrange + 1, alpha=0.35):
    """Draw a shaded systematic uncertainty band around central values (asymmetric)."""
    n = len(xs)
    band = ROOT.TGraphAsymmErrors(n)
    band.SetName("band_sys")
    for i in range(n):
        band.SetPoint(i, xs[i], central_ys[i])
        band.SetPointError(i, exs[i], exs[i], sys_lows[i], sys_highs[i])
    band.SetFillColorAlpha(color, alpha)
    band.SetLineWidth(0)
    band.Draw("E2 SAME")
    return band


def draw_tot_sys_box(xs, central_ys, tot_lows, tot_highs, exs, color=ROOT.kGray, alpha=0.25):
    """Draw total systematic uncertainty as boxes for each pT bin (asymmetric)."""
    boxes = []
    for i in range(len(xs)):
        xlo = xs[i] - exs[i]
        xhi = xs[i] + exs[i]
        ylo = central_ys[i] - tot_lows[i]
        yhi = central_ys[i] + tot_highs[i]
        box = ROOT.TBox(xlo, ylo, xhi, yhi)
        box.SetFillColorAlpha(color, alpha)
        box.SetLineWidth(0)
        box.Draw("SAME")
        boxes.append(box)
    return boxes


def compare_2pc_sp(
    sys_file, sys_hname, sys_gname_sys, sys_gname_tot,
    sp_file, sp_hname,
    output_dir, suffix,
    legend_2pc="2pc prompt", legend_sp="Biao sp prompt",
    ymin=-0.01, ymax=0.20,
):
    """
    Main comparison function.

    Parameters
    ----------
    sys_file : str
        ROOT file from compute_prompt_v2_unfold_sys.py
    sys_hname : str
        Name of the central TH1 with stat errors (e.g. "hV2Prompt")
    sys_gname_sys : str
        Name of the TGraphErrors with sys uncertainty (e.g. "gSysUnc")
    sys_gname_tot : str
        Name of the TGraphErrors with total uncertainty (e.g. "gTotUnc")
    sp_file : str
        ROOT file with Biao's sp result
    sp_hname : str
        Name of the sp TH1 (e.g. "hV2VsPtPrompt")
    output_dir : str
        Directory to save output plots
    suffix : str
        Suffix for output filenames
    """
    # ── Colors & styles ────────────────────────────────────────────
    color_2pc = ROOT.kRed + 1
    color_sp = ROOT.kBlue + 1
    color_sys_band = ROOT.kOrange + 1
    color_tot_box = ROOT.kGray + 1

    # ── Read 2pc sys result ────────────────────────────────────────
    f2pc = ROOT.TFile.Open(sys_file, "READ")
    if not f2pc or f2pc.IsZombie():
        raise RuntimeError(f"Cannot open {sys_file}")

    h_central = f2pc.Get(sys_hname)
    h_central.Scale(0.07/0.0625) 
    print("Using the 0.07/0.0625")
    input("Enter to continue")
    g_sys_asym = f2pc.Get(sys_gname_sys)
    g_tot = f2pc.Get(sys_gname_tot)  # optional — computed from scratch below

    if not h_central:
        raise RuntimeError(f"Missing {sys_hname} in {sys_file}")
    if not g_sys_asym:
        raise RuntimeError(f"Missing {sys_gname_sys} in {sys_file}")
    # gTotUnc is optional — tot_lows/tot_highs computed from sys+stat

    h_central.SetDirectory(0)

    xs, central_ys, stat_eys, exs = read_hist_bins(h_central)
    _, _, _, sys_lows, sys_highs = read_graph_asymm_points(g_sys_asym)
    _, tot_ys, _, _ = read_graph_points(g_tot)
    # asymmetric total: sqrt(stat² + sys_low²), sqrt(stat² + sys_high²)
    tot_lows = [math.sqrt(stat_eys[i] ** 2 + sys_lows[i] ** 2) for i in range(len(xs))]
    tot_highs = [math.sqrt(stat_eys[i] ** 2 + sys_highs[i] ** 2) for i in range(len(xs))]

    # ── Read sp result ──────────────────────────────────────────────
    fsp = ROOT.TFile.Open(sp_file, "READ")
    if not fsp or fsp.IsZombie():
        raise RuntimeError(f"Cannot open {sp_file}")

    h_sp = fsp.Get(sp_hname)
    if not h_sp:
        raise RuntimeError(f"Missing {sp_hname} in {sp_file}")
    h_sp.SetDirectory(0)

    sp_xs, sp_ys, sp_eys, sp_exs = read_hist_bins(h_sp)

    # ── Build TGraphErrors for 2pc (stat only) ─────────────────────
    n_2pc = len(xs)
    g_2pc_stat = ROOT.TGraphErrors(n_2pc)
    g_2pc_stat.SetName("g2pc_stat")
    for i in range(n_2pc):
        g_2pc_stat.SetPoint(i, xs[i], central_ys[i])
        g_2pc_stat.SetPointError(i, exs[i], stat_eys[i])
    g_2pc_stat.SetLineColor(color_2pc)
    g_2pc_stat.SetLineWidth(3)
    g_2pc_stat.SetMarkerColor(color_2pc)
    g_2pc_stat.SetMarkerStyle(20)
    g_2pc_stat.SetMarkerSize(1.8)

    # ── Build TGraphErrors for sp ───────────────────────────────────
    n_sp = len(sp_xs)
    g_sp = ROOT.TGraphErrors(n_sp)
    g_sp.SetName("g_sp")
    for i in range(n_sp):
        g_sp.SetPoint(i, sp_xs[i], sp_ys[i])
        g_sp.SetPointError(i, sp_exs[i], sp_eys[i])
    g_sp.SetLineColor(color_sp)
    g_sp.SetLineWidth(3)
    g_sp.SetMarkerColor(color_sp)
    g_sp.SetMarkerStyle(21)
    g_sp.SetMarkerSize(1.8)

    # ── Canvas 1: Main comparison ───────────────────────────────────
    c = ROOT.TCanvas("c_compare", "2pc vs SP comparison", 1600, 1200)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.10)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.04)

    frame = c.DrawFrame(0, ymin, 12, ymax*1.1)
    frame.GetYaxis().SetTitle("D^{0} v_{2}")
    frame.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
    frame.GetYaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitleSize(0.045)
    frame.GetXaxis().SetLabelSize(0.04)

    # Draw sys band FIRST (behind data points)
    band_sys = draw_sys_band(xs, central_ys, sys_lows, sys_highs, exs, color=color_sys_band, alpha=0.25)

    # Draw tot_sys box (asymmetric)
    boxes_tot = draw_tot_sys_box(xs, central_ys, tot_lows, tot_highs, exs, color=color_tot_box, alpha=0.15)

    # Draw 2pc stat
    g_2pc_stat.Draw("P E1 SAME")

    # Draw sp
    g_sp.Draw("P E1 SAME")

    # ── Legend ──────────────────────────────────────────────────────
    x_legend = 0.18
    y_legend = 0.65
    legend = ROOT.TLegend(x_legend, y_legend, x_legend + 0.35, y_legend + 0.22)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    legend.AddEntry(g_2pc_stat, legend_2pc, "lp")
    legend.AddEntry(g_sp, legend_sp, "lp")
    legend.AddEntry(band_sys, "sys. unc. (ratio scan)", "f")
    legend.AddEntry(boxes_tot[0] if boxes_tot else None, "tot. unc. (stat#oplussys)", "f")
    legend.Draw()

    # ALICE label
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    # latex.DrawLatex(0.18, 0.88, "This work, OO #sqrt{#it{s}_{NN}} = 5.36 TeV")
    # latex.DrawLatex(0.18, 0.83, "Prompt D^{0}, 0-20% c.c.")

    c.Update()

    # ── Canvas 2: Ratio plot ────────────────────────────────────────
    c_ratio = ROOT.TCanvas("c_ratio", "Ratio to SP", 1600, 1200)
    c_ratio.SetLeftMargin(0.12)
    c_ratio.SetBottomMargin(0.10)
    c_ratio.SetRightMargin(0.04)
    c_ratio.SetTopMargin(0.04)

    x_min = min(min(xs), min(sp_xs)) - 0.5
    x_max = max(max(xs), max(sp_xs)) + 0.5

    frame_ratio = c_ratio.DrawFrame(x_min, 0, 8, 2)
    frame_ratio.GetYaxis().SetTitle("Ratio to SP")
    frame_ratio.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
    frame_ratio.GetYaxis().SetTitleSize(0.045)
    frame_ratio.GetYaxis().SetLabelSize(0.04)
    frame_ratio.GetXaxis().SetTitleSize(0.045)
    frame_ratio.GetXaxis().SetLabelSize(0.04)

    # Build ratio graph: 2pc / sp (bin-by-bin for overlapping pT range)
    n_ratio = min(n_2pc, n_sp)
    g_ratio = ROOT.TGraphErrors(n_ratio)
    g_ratio.SetName("g_ratio")
    g_ratio.SetLineColor(color_2pc)
    g_ratio.SetLineWidth(3)
    g_ratio.SetMarkerColor(color_2pc)
    g_ratio.SetMarkerStyle(20)
    g_ratio.SetMarkerSize(1.8)

    g_ratio_sys = ROOT.TGraphAsymmErrors(n_ratio)
    g_ratio_sys.SetName("g_ratio_sys")
    g_ratio_sys.SetFillColorAlpha(color_sys_band, 0.25)
    g_ratio_sys.SetLineWidth(0)
    
    g_ratio_tot = ROOT.TGraphAsymmErrors(n_ratio)
    g_ratio_tot.SetName("g_ratio_tot")
    g_ratio_tot.SetFillColorAlpha(color_tot_box, 0.15)

    for i in range(n_ratio):
        num = central_ys[i]
        den = sp_ys[i]
        if den != 0:
            ratio_val = num / den
            g_ratio.SetPoint(i, xs[i], ratio_val)
            # stat error propagation
            ratio_stat = ratio_val * math.sqrt(
                (stat_eys[i] / num) ** 2 + (sp_eys[i] / den) ** 2
            ) if num != 0 else 0
            g_ratio.SetPointError(i, exs[i], ratio_stat)
            # # sys error propagation (asymmetric)
            # ratio_sys_low = ratio_val * (sys_lows[i] / num) if num != 0 else 0
            # ratio_sys_high = ratio_val * (sys_highs[i] / num) if num != 0 else 0
            # g_ratio_sys.SetPoint(i, xs[i], ratio_val)
            # g_ratio_sys.SetPointError(i, exs[i], exs[i], ratio_sys_low, ratio_sys_high)
            # # total error propagation (asymmetric)
            # ratio_tot_low = ratio_val * (tot_lows[i] / num) if num != 0 else 0
            # ratio_tot_high = ratio_val * (tot_highs[i] / num) if num != 0 else 0
            # g_ratio_tot.SetPoint(i, xs[i], ratio_val)
            # g_ratio_tot.SetPointError(i, exs[i], exs[i], ratio_tot_low, ratio_tot_high)
        else:
            g_ratio.SetPoint(i, xs[i], 0)
            g_ratio.SetPointError(i, exs[i], 0)
            g_ratio_sys.SetPoint(i, xs[i], 0)
            g_ratio_sys.SetPointError(i, exs[i], exs[i], 0, 0)

    # g_ratio_sys.Draw("E2 SAME")
    g_ratio.Draw("P E1 SAME")
    # g_ratio_tot.Draw("E2 SAME")

    # Unity line
    line = ROOT.TLine(x_min, 1, 8, 1)
    line.SetLineStyle(7)
    line.SetLineWidth(2)
    line.SetLineColor(ROOT.kBlack)
    line.Draw("SAME")

    # Ratio legend
    leg_ratio = ROOT.TLegend(0.50, 0.18, 0.88, 0.36)
    leg_ratio.SetBorderSize(0)
    leg_ratio.SetFillStyle(0)
    leg_ratio.SetTextFont(42)
    leg_ratio.SetTextSize(0.035)
    leg_ratio.AddEntry(g_ratio, "2pc / SP (stat only)", "lp")
    # leg_ratio.AddEntry(g_ratio_sys, "sys. unc. (ratio scan)", "f")
    # leg_ratio.AddEntry(g_ratio_tot, "total unc.", "f")
    leg_ratio.Draw()

    c_ratio.Update()

    # ── Save outputs ─────────────────────────────────────────────────
    os.makedirs(output_dir, exist_ok=True)

    png_main = os.path.join(output_dir, f"compare_{suffix}.png")
    png_ratio = os.path.join(output_dir, f"compare_ratio_{suffix}.png")
    root_out = os.path.join(output_dir, f"compare_{suffix}.root")

    c.SaveAs(png_main)
    c_ratio.SaveAs(png_ratio)
    print(f"[OK] Main plot: {png_main}")
    print(f"[OK] Ratio plot: {png_ratio}")

    # Save ROOT file
    fout = ROOT.TFile.Open(root_out, "RECREATE")
    c.Write("c_compare")
    c_ratio.Write("c_ratio")
    g_2pc_stat.Write()
    g_sp.Write()
    g_ratio.Write()
    fout.Close()
    print(f"[OK] ROOT file: {root_out}")

    f2pc.Close()
    fsp.Close()


# ── Math import (needed for ratio error propagation) ─────────────────
import math


# ── Main ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Compare 2pc vs SP with sys uncertainties")
    ap.add_argument("--sys-file", type=str,
                    default="/home/wuct/ALICE/reps/hf-vn-dev/dev/src/v2_prompt_ratio_0d2_1d3.root",
                    help="ROOT file from compute_prompt_v2_unfold_sys.py")
    ap.add_argument("--sys-hname", type=str, default="hV2Prompt",
                    help="Name of central TH1 in sys file")
    ap.add_argument("--sys-gname-sys", type=str, default="gSysUncAsym",
                    help="Name of sys uncertainty TGraphAsymmErrors")
    ap.add_argument("--sys-gname-tot", type=str, default="gTotUnc",
                    help="Name of total uncertainty TGraphErrors")
    ap.add_argument("--sp-file", type=str,
                    default="/home/wuct/MetaData/DATA/OO/apass2/corr/results/fthlook_Biao/large/old/k020/v2VsFracD0020Biao.root",
                    help="ROOT file with Biao's sp result")
    ap.add_argument("--sp-hname", type=str, default="hV2VsPtPrompt",
                    help="Name of sp TH1")
    ap.add_argument("--output-dir", type=str,
                    default="/home/wuct/MetaData/DATA/OO/apass2/corr/results/fifth/k020_gausPer/etaVariation/comparison/prompt_sys_new_v2hh",
                    help="Output directory for plots")
    ap.add_argument("--suffix", type=str, default="0d2_prompt_sys",
                    help="Suffix for output filenames")
    ap.add_argument("--legend-2pc", type=str,
                    default="2pc prompt 0.2<|#Delta#eta|<1.3 0-20% (r=0.5)")
    ap.add_argument("--legend-sp", type=str,
                    default="Biao sp prompt 0-20%")
    ap.add_argument("--ymin", type=float, default=-0.1)
    ap.add_argument("--ymax", type=float, default=0.25)

    args = ap.parse_args()

    compare_2pc_sp(
        sys_file=args.sys_file,
        sys_hname=args.sys_hname,
        sys_gname_sys=args.sys_gname_sys,
        sys_gname_tot=args.sys_gname_tot,
        sp_file=args.sp_file,
        sp_hname=args.sp_hname,
        output_dir=args.output_dir,
        suffix=args.suffix,
        legend_2pc=args.legend_2pc,
        legend_sp=args.legend_sp,
        ymin=args.ymin,
        ymax=args.ymax,
    )
