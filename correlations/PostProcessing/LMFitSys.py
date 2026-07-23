#!/usr/bin/env python3
"""
LM fit systematic uncertainty — vary LM template parameters via
multivariate Gaussian sampling, re-fit to get v2 spread.

Workflow per (PtCand, PtHad):
  1. Read data histogram from CorrPhiD0 canvas ROOT file
  2. Read LM template (fLMOutput) + covariance + parameters from hLMtemplate ROOT
  3. Sample LM parameters from multivariate Gaussian (N samples)
  4. Build varied GausPeriodic LM templates
  5. Fit with DhCorrelationFitter -> collect v2_delta per sample
  6. Plot v2 vs pT with systematic band, and probability vs v2 per pT
"""

import argparse
import gc
import math
import os
import sys
from array import array
from ctypes import c_double, c_int

import numpy as np
from scipy.stats import multivariate_normal
import ROOT
from alive_progress import alive_bar

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gStyle.SetOptStat(0)  # No stat boxes
ROOT.gStyle.SetPalette(ROOT.kRainBow)  # kRainBow = 55

# ── Compile DhCorrelationFitter ────────────────────────────────────────
ROOT.gSystem.AddIncludePath("-I/home/wuct/Software/miniforge3/envs/alice/include")
_fitter_cxx = os.path.join(os.path.dirname(__file__), "DhCorrelationFitter.cxx")
ROOT.gSystem.CompileMacro(_fitter_cxx, "kO+")  # force recompile with ACLiC
from ROOT import DhCorrelationFitter

# ── Paths (from FitCorrel config) ──────────────────────────────────────
BASE = "/home/wuct/MetaData/DATA/OO/apass2/corr/results/fifth/k020_gausPer/etaVariation"
EXTRACT = f"{BASE}/CorrelExtract_0d2_1d3_AppDeltaPhi"
OUT_ROOT = f"{EXTRACT}/CorrelationFitResults/Output_CorrelationFitting_Root"
OUT_PNG  = f"{EXTRACT}/CorrelationFitResults/outputPathOutput_CorrelationFitting_png"
OUT_LM   = f"{EXTRACT}/CorrelationFitResults/LMSysResults"
os.makedirs(OUT_LM, exist_ok=True)

# Config values (mirror FitCorrel config)
DMESON = "D0"
PT_CAND = [0., 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6., 7., 8., 12.]
PT_HAD  = [0.2, 3.0]
MASS_RANGE = (1.72, 2.02)
TEMP_FUNC = 4  # GausPeriodic
FIX_LM_FACTOR = True
WITH_PED_LM = False
FIX_BASELINE = -4
FIT_FUNC_TYPE = 10  # kTemplateFit
F_MIN = -0.5 * math.pi
F_MAX = 1.5 * math.pi
N_SAMPLES = 100  # Default, override with --n-samples
RY_TRIGGER = None  # Will be read from file or set from config
LM_RY_TRIGGER = None
RY_ERR = 0.
LM_RY_ERR = 0.

_GAUS_PER_FORMULA = (
    "[0]"
    "+[1]/(TMath::Sqrt(2*TMath::Pi())*[2])"
    "*TMath::Exp(-x*x/(2*[2]*[2]))"
    "+[1]/(TMath::Sqrt(2*TMath::Pi())*[2])"
    "*TMath::Exp(-(x-2*TMath::Pi())*(x-2*TMath::Pi())/(2*[2]*[2]))"
    "+[1]/(TMath::Sqrt(2*TMath::Pi())*[2])"
    "*TMath::Exp(-(x+2*TMath::Pi())*(x+2*TMath::Pi())/(2*[2]*[2]))"
    "+[3]/(TMath::Sqrt(2*TMath::Pi())*[4])"
    "*TMath::Exp(-(x-TMath::Pi())*(x-TMath::Pi())/(2*[4]*[4]))"
    "+[3]/(TMath::Sqrt(2*TMath::Pi())*[4])"
    "*TMath::Exp(-(x-3*TMath::Pi())*(x-3*TMath::Pi())/(2*[4]*[4]))"
    "+[3]/(TMath::Sqrt(2*TMath::Pi())*[4])"
    "*TMath::Exp(-(x+TMath::Pi())*(x+TMath::Pi())/(2*[4]*[4]))"
)

def _hsv_to_rgb(h, s, v):
    """Convert HSV (0..1) to RGB (0..1)."""
    if s == 0.0:
        return v, v, v
    i = int(h * 6.0)
    f = h * 6.0 - i
    p = v * (1.0 - s)
    q = v * (1.0 - s * f)
    t = v * (1.0 - s * (1.0 - f))
    i %= 6
    if i == 0: return v, t, p
    if i == 1: return q, v, p
    if i == 2: return p, v, t
    if i == 3: return p, q, v
    if i == 4: return t, p, v
    return v, p, q

def get_custom_palette(name="rainbow", n_colors=100):
    """
    Generate a list of ROOT TColor indices for a specified palette.
    Available names: 'rainbow', 'coolwarm', 'viridis', 'heat', 'gray'
    """
    colors = []
    for i in range(n_colors):
        # Fractional position along the palette [0.0, 1.0]
        frac = float(i) / max(n_colors - 1, 1)

        if name == "rainbow":
            # Red -> Green -> Blue -> Violet
            h = frac
            s = 0.8
            v = 1.0

        elif name == "coolwarm":
            # Diverging: Blue -> White/Gray -> Red
            h = 0.66 * (1.0 - frac)     # 0.66 is blue, 0.0 is red
            s = abs(frac - 0.5) * 2.0   # Saturation dips to 0 in the middle
            v = 0.95                    # Keep overall brightness high

        elif name == "viridis":
            # Sequential: Dark Purple -> Teal -> Green -> Yellow
            h = 0.8 - frac * 0.65       # Shift hue from purple to yellow
            s = 0.4 + 0.6 * frac        # Increase saturation
            v = 0.3 + 0.7 * frac        # Dark to light

        elif name == "heat":
            # Sequential: Black -> Red -> Orange -> Yellow -> White
            h = 0.16 * frac             # Shift hue slightly towards yellow
            s = 1.0 - max(0.0, (frac - 0.5) * 2.0) # Desaturate at the high end
            v = min(1.0, frac * 2.0)    # Ramp up brightness quickly

        elif name == "gray":
            # Grayscale: Black -> White
            h, s = 0.0, 0.0
            v = frac

        else:
            # Fallback to simple rainbow
            h, s, v = frac, 0.8, 1.0

        r, g, b = _hsv_to_rgb(h, s, v)
        c_idx = ROOT.TColor.GetColor(int(r * 255), int(g * 255), int(b * 255))
        colors.append(c_idx)

    return colors

def build_lm_tf1(par, name="fLMTemplate"):
    """
    Build a TF1 GausPeriodic with given parameters (5-element array).
    Bypasses the internal LM template fit in DhCorrelationFitter.
    """
    func = ROOT.TF1(name, _GAUS_PER_FORMULA, F_MIN, F_MAX)
    for i in range(5):
        func.SetParameter(i, par[i])
    func.SetNpx(3 * 100)  # 3x oversampling for smoothness
    ROOT.SetOwnership(func, False)  # Python must not delete TF1; fitter holds raw pointer
    return func

def extract_lm_data_histo(lm_root_path):
    f = ROOT.TFile.Open(lm_root_path)
    if not f or f.IsZombie(): return None
    for k in f.GetListOfKeys():
        obj = f.Get(k.GetName())
        if obj and obj.InheritsFrom("TCanvas"):
            for ip in range(obj.GetListOfPrimitives().GetSize()):
                p = obj.GetListOfPrimitives().At(ip)
                if p and "LM" in p.GetName():
                    h = p.Clone("hLMData"); h.SetDirectory(0)
                    f.Close(); return h
            break
    f.Close(); return None

def extract_data_histo(canvas_root_path):
    """Extract the data histogram (hPairsYields) from a canvas ROOT file."""
    f = ROOT.TFile.Open(canvas_root_path)
    if not f or f.IsZombie():
        print(f"[ERROR] Cannot open {canvas_root_path}")
        return None, None

    keys = [k.GetName() for k in f.GetListOfKeys()]
    canvas = None
    for k in keys:
        obj = f.Get(k)
        if obj and obj.InheritsFrom("TCanvas"):
            canvas = obj
            break

    if not canvas:
        print(f"[ERROR] No TCanvas found in {canvas_root_path}")
        f.Close()
        return None, None

    h_data = None
    prims = canvas.GetListOfPrimitives()
    for ip in range(prims.GetSize()):
        obj = prims.At(ip)
        if obj and obj.InheritsFrom("TH1"):
            name = obj.GetName()
            if "h_corr" in name.lower() or "corr" in name.lower():
                h_data = obj.Clone("hData")
                h_data.SetDirectory(0)  # Safe PyROOT memory detachment
                break
                
    if not h_data:
        # Fallback: take first TH1
        for ip in range(prims.GetSize()):
            obj = prims.At(ip)
            if obj and obj.InheritsFrom("TH1") and "hFit" in obj.GetName():
                h_data = obj.Clone("hData")
                h_data.SetDirectory(0)
                break

    f.Close()
    if not h_data:
        print(f"[ERROR] No data histogram found in canvas")
    return h_data, canvas

def extract_lm_covariance(lm_root_path, pt_cand_idx, pt_had_idx, mass_label):
    """Extract LM fit covariance matrix and parameter means from hLMtemplate ROOT."""
    f = ROOT.TFile.Open(lm_root_path)
    if not f or f.IsZombie():
        print(f"[ERROR] Cannot open {lm_root_path}")
        return None, None, None

    cov_name = f"hLMCovMatrix_PtCand{pt_cand_idx}_PtHad{pt_had_idx}_InvMassBin{mass_label}"
    h_cov = f.Get(cov_name)
    if not h_cov:
        print(f"[ERROR] Cov TH2D '{cov_name}' not found in {lm_root_path}")
        f.Close()
        return None, None, None

    npar = h_cov.GetNbinsX()
    mean_vec = np.zeros(npar)
    cov_mat = np.zeros((npar, npar))

    for i in range(npar):
        for j in range(npar):
            cov_mat[i, j] = h_cov.GetBinContent(i + 1, j + 1)

    par_name = f"hLMParValues_PtCand{pt_cand_idx}_PtHad{pt_had_idx}_InvMassBin{mass_label}"
    h_par = f.Get(par_name)
    if h_par:
        for i in range(npar):
            mean_vec[i] = h_par.GetBinContent(i + 1)
    else:
        print(f"[WARNING] Par TH1D '{par_name}' not found, using covariance diag only")

    f.Close()
    return mean_vec, cov_mat, npar

def read_lm_factor(canvas_root_path, pt_cand_idx, pt_had_idx, mass_label):
    """Read F (LM Scale) from hParValues TH1D in the canvas ROOT file."""
    f = ROOT.TFile.Open(canvas_root_path)
    if not f or f.IsZombie():
        return None, None
    hname = f"hParValues_PtCand{pt_cand_idx}_PtHad{pt_had_idx}_InvMassBin{mass_label}"
    h = f.Get(hname)
    if not h:
        for k in f.GetListOfKeys():
            kn = k.GetName()
            if "hParValues" in kn and f"PtCand{pt_cand_idx}" in kn:
                h = f.Get(kn)
                break
    if not h:
        f.Close()
        return None, None
        
    f_val = h.GetBinContent(1)
    f_err = h.GetBinError(1)
    f.Close()
    return f_val, f_err

def fit_with_lm_template(h_data, par, ry_val=1.0, ry_err=0.0, lm_ry_val=1.0, lm_ry_err=0.0, tag=""):
    """Fit data with externally-set LM template function.
    `tag` is a unique suffix to avoid ROOT object name collisions in loops."""
    h_fit = ROOT.TH1F(f"hFit_{tag}", "", h_data.GetNbinsX(),
                       h_data.GetXaxis().GetXmin(),
                       h_data.GetXaxis().GetXmax())
    for ib in range(1, h_data.GetNbinsX() + 1):
        h_fit.SetBinContent(ib, h_data.GetBinContent(ib))
        h_fit.SetBinError(ib, h_data.GetBinError(ib))
    ROOT.SetOwnership(h_fit, False)

    c = ROOT.TCanvas(f"c_fit_{tag}", "", 1840, 1126)
    c.SetBottomMargin(0.08); c.SetLeftMargin(0.12)
    c.SetRightMargin(0.02); c.SetTopMargin(0.1)
    c.cd()

    fitter = DhCorrelationFitter(h_fit, F_MIN, F_MAX)
    ROOT.SetOwnership(fitter, True)

    fitter.SetHistoIsReflected(False)
    fitter.SetWithPedLM(WITH_PED_LM)
    fitter.SetFixLMFactor(FIX_LM_FACTOR)
    fitter.SetTempFunc(TEMP_FUNC)
    fitter.SetFixBaseline(FIX_BASELINE)
    fitter.SetBaselineUpOrDown(False, False)
    fitter.SetReflectedCorrHisto(True)
    fitter.SetFixMean(0)
    fitter.SetRyTrigger(ry_val, ry_err)
    fitter.SetLMRyTrigger(lm_ry_val, lm_ry_err)
    # TH1D (fTempHisto) + params (fixed, skip LM fit in BuildLMOutput)
    h_lm = ROOT.TH1D(f"hLM_{tag}", "", h_data.GetNbinsX(),
                      h_data.GetXaxis().GetXmin(), h_data.GetXaxis().GetXmax())
    f_tmp = build_lm_tf1(par, f"_tmp_{tag}")
    for ib in range(1, h_data.GetNbinsX() + 1):
        h_lm.SetBinContent(ib, f_tmp.Eval(h_lm.GetBinCenter(ib)))
    h_lm.SetDirectory(0)
    fitter.SetLMTemplate(h_lm)
    c_arr = (c_double * len(par))(*par)
    fitter.SetLMTemplateParams(len(par), c_arr)
    fitter.SetFuncType(DhCorrelationFitter.FunctionType(FIT_FUNC_TYPE))

    try:
        fitter.Fitting(True, False)
    except Exception as e:
        import traceback
        print(f"[WARNING] Fit failed: {e}")
        traceback.print_exc()
        c.Close(); del fitter, h_fit
        return None, None, None, None

    h_fit.GetYaxis().SetRangeUser(h_fit.GetMinimum() * 0.96, h_fit.GetMaximum() * 1.05)
    c.Update()

    h_data.SetStats(0); h_data.SetMinimum(0)
    h_data.SetMarkerStyle(ROOT.kFullCircle)
    h_data.SetMarkerColor(ROOT.kRed + 1)
    h_data.SetMarkerSize(1.4)
    h_data.SetLineColor(ROOT.kRed + 1)
    h_data.SetLineWidth(3)
    h_data.Draw("same")

    v2 = fitter.Getv2Delta()
    v2_err = fitter.Getv2DeltaError()

    return v2, v2_err, c, fitter

def save_trial(trial_dir, canvas, fitter, f_fit, f_lm,
               pt_cand_idx, pt_had_idx, mass_label, trial_id,
               pt_cand_min, pt_cand_max, pt_had_min, pt_had_max,
               mass_min, mass_max):
    """Save a single trial result: canvas+params to trial_dir/."""
    os.makedirs(trial_dir, exist_ok=True)

    # Write v2_result.txt FIRST — canvas.SaveAs may modify fitter state
    with open(os.path.join(trial_dir, "v2_result.txt"), "w") as fp:
        fp.write(f"v2_delta = {fitter.Getv2Delta():.8f} +/- {fitter.Getv2DeltaError():.8f}\n")
        if f_fit and f_fit.GetNpar() > 3:
            fp.write(f"v3_delta = {f_fit.GetParameter(3):.8f} +/- {f_fit.GetParError(3):.8f}\n")

    canvas.SaveAs(os.path.join(trial_dir, "CorrFit.png"))
    canvas.SaveAs(os.path.join(trial_dir, "CorrFit.root"))

    if f_fit:
        npar = f_fit.GetNpar()
        with open(os.path.join(trial_dir, "parFit.txt"), "w") as fp:
            fp.write(f"# Trial {trial_id}  PtCand [{pt_cand_min:.1f},{pt_cand_max:.1f}]"
                     f"  PtHad [{pt_had_min:.1f},{pt_had_max:.1f}]\n")
            fp.write(f"{'Param':<6} {'Name':<20} {'Value':<18} {'Error':<18}\n")
            fp.write("-" * 60 + "\n")
            for ip in range(npar):
                fp.write(f"{ip:<6} {f_fit.GetParName(ip):<20} "
                         f"{f_fit.GetParameter(ip):<18.8f} {f_fit.GetParError(ip):<18.8f}\n")

    if f_lm:
        npar_lm = f_lm.GetNpar()
        with open(os.path.join(trial_dir, "parLM.txt"), "w") as fp:
            fp.write(f"# Trial {trial_id}  LM tempFunc={TEMP_FUNC}\n")
            fp.write(f"{'Param':<6} {'Name':<20} {'Value':<18} {'Error':<18}\n")
            fp.write("-" * 60 + "\n")
            for ip in range(npar_lm):
                fp.write(f"{ip:<6} {f_lm.GetParName(ip):<20} "
                         f"{f_lm.GetParameter(ip):<18.8f} {f_lm.GetParError(ip):<18.8f}\n")


def _read_v2_from_txt(path):
    if not os.path.exists(path):
        return None, None
    with open(path, "r") as f:
        for line in f:
            if line.startswith("v2_delta"):
                parts = line.strip().split()
                return float(parts[2]), float(parts[4])
    return None, None


def _read_parlm_from_txt(path):
    if not os.path.exists(path):
        return None
    pars = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("---") or line.startswith("Param"):
                continue
            parts = line.rsplit(maxsplit=2)
            if len(parts) >= 3:
                pars.append(float(parts[-2]))
    return pars if len(pars) == 5 else None


def run_systematics(n_samples=100, load_from_exist=False, pt_only=0):
    """Run full systematic uncertainty estimation for all PtCand bins."""
    N_DIGITS = int(math.ceil(math.log10(n_samples + 1)))

    pt_cand_bins = PT_CAND
    n_pt_cand = len(pt_cand_bins) - 1
    pt_had_idx = 1
    mass_label = 1

    # ── Results storage ──────────────────────────────────────────────
    all_v2_central = np.zeros(n_pt_cand)
    all_v2_err     = np.zeros(n_pt_cand)
    all_prob_central = np.zeros(n_pt_cand)
    all_v2_samples = {}
    all_probs      = {}
    all_fit_funcs  = {}
    all_h_data     = {}
    all_par_samples = {}
    all_lm_data    = {}
    all_central_pars = {}  # <--- [新增] 用于存储每个 pT bin 的中心参数

    rainbow_colors = get_custom_palette("rainbow", n_colors=n_samples)

    for i_pt in range(n_pt_cand):
        pt_cand_idx = i_pt + 1
        if pt_only > 0 and pt_cand_idx != pt_only:
            continue
        pt_min = pt_cand_bins[i_pt]
        pt_max = pt_cand_bins[i_pt + 1]
        print(f"\n{'='*60}")
        print(f"[PtCand {pt_cand_idx}/{n_pt_cand}]  [{pt_min}, {pt_max}] GeV/c")

        canvas_path = os.path.join(OUT_ROOT, f"CorrPhi{DMESON}_PtBinCand{pt_cand_idx}_PtBinAssoc{pt_had_idx}_InvMassBin{mass_label}.root")
        h_data, _ = extract_data_histo(canvas_path)
        if not h_data:
                print(f"  [SKIP] No data histogram found")
                continue
        print(f"  Data: {h_data.GetNbinsX()} bins in [{h_data.GetXaxis().GetXmin():.2f}, {h_data.GetXaxis().GetXmax():.2f}]")
        all_h_data[i_pt] = h_data

        if load_from_exist:
            pt_dir = os.path.join(OUT_LM, f"PtCand_{pt_cand_idx}")
            if not os.path.isdir(pt_dir):
                print(f"  PtCand {pt_cand_idx}: no directory, skip")
                continue

            lm_path = os.path.join(OUT_ROOT, f"hLMtemplate_PtBinCand{pt_cand_idx}_PtBinAssoc{pt_had_idx}.root")
            mean_vec, cov_mat, npar = extract_lm_covariance(lm_path, pt_cand_idx, pt_had_idx, mass_label)
            if mean_vec is None:
                print(f"  [SKIP] No LM covariance found")
                continue

            # Read central from trial_0000
            central_path = os.path.join(pt_dir, f"trial_{0:0{N_DIGITS}d}", "v2_result.txt")
            v2_cent, v2_cent_err = _read_v2_from_txt(central_path)
            central_par = _read_parlm_from_txt(os.path.join(pt_dir, f"trial_{0:0{N_DIGITS}d}", "parLM.txt"))
            if central_par is None:
                central_par = mean_vec.copy()
                print(f"  [WARNING] Cannot read central parLM, using mean_vec")
            else:
                # Convert list to numpy array for consistency
                central_par = np.array(central_par)
            print(f"  Central LM par = {central_par}")
            lm_func_central = build_lm_tf1(central_par, "lm_central")
            all_v2_central[i_pt] = v2_cent if v2_cent is not None else 0.
            all_v2_err[i_pt] = v2_cent_err if v2_cent_err is not None else 0.

            dist = multivariate_normal(mean=mean_vec, cov=cov_mat, allow_singular=True)

            v2_list, prob_list, lm_funcs, par_list, pdf_list = [], [], [], [], []
            with alive_bar(n_samples + 1, title=f"Loading PtCand {pt_cand_idx} trials") as bar:
                for itrial in range(0, n_samples + 1):
                    trial_path = os.path.join(pt_dir, f"trial_{itrial:0{N_DIGITS}d}")
                    if not os.path.isdir(trial_path):
                        break
                    v2_val, _ = _read_v2_from_txt(os.path.join(trial_path, "v2_result.txt"))
                    if v2_val is None:
                        continue
                    par_lm = _read_parlm_from_txt(os.path.join(trial_path, "parLM.txt"))
                    if par_lm is None:
                        continue
                    v2_list.append(v2_val)
                    pdf_list.append(dist.pdf(par_lm))
                    par_list.append(par_lm)
                    lm_funcs.append(build_lm_tf1(par_lm, f"lm_{itrial}"))
                    bar()  # Finish progress bar
            # for MC sampling, rvs
            prob_list = pdf_list
            # # for Uniform Grid Scan
            # # Normalize probabilities
            # pdf_sum = np.sum(pdf_list) if pdf_list else 1.0
            # prob_list = [p / pdf_sum for p in pdf_list]
            all_prob_central[i_pt] = prob_list[0] if prob_list else 1.0
            all_v2_samples[i_pt] = np.array(v2_list)
            all_probs[i_pt] = np.array(prob_list)
            all_fit_funcs[i_pt] = (lm_func_central, lm_funcs)
            all_par_samples[i_pt] = par_list
            all_lm_data[i_pt] = extract_lm_data_histo(lm_path)
            all_central_pars[i_pt] = central_par
            print(f"  PtCand {pt_cand_idx}: central v2={all_v2_central[i_pt]:.6f}, loaded {len(v2_list)} trials, "
                  f"lm_data={'OK' if all_lm_data[i_pt] else 'MISSING'}, "
                  f"lm_central={'OK' if lm_func_central else 'MISSING'}")
            continue

        f_lm, f_lm_err = read_lm_factor(canvas_path, pt_cand_idx, pt_had_idx, mass_label)
        if f_lm is not None:
                ry_val, lm_ry_val, ry_err, lm_ry_err = f_lm, 1.0, 0.0, 0.0
                print(f"  LM Factor F = {f_lm:.6f}  (ry={ry_val:.1f}, lm_ry={lm_ry_val:.1f})")
        else:
                print(f"  [WARNING] LM Factor not found, using ry=lm_ry=1")
                ry_val, lm_ry_val, ry_err, lm_ry_err = 1.0, 1.0, 0.0, 0.0

        lm_path = os.path.join(OUT_ROOT, f"hLMtemplate_PtBinCand{pt_cand_idx}_PtBinAssoc{pt_had_idx}.root")
        mean_vec, cov_mat, npar = extract_lm_covariance(lm_path, pt_cand_idx, pt_had_idx, mass_label)
        if mean_vec is None:
                print(f"  [SKIP] No LM covariance found")
                continue
        print(f"  LM params: {npar} ({mean_vec})")

        dist = multivariate_normal(mean=mean_vec, cov=cov_mat, allow_singular=True)
        samples = dist.rvs(size=n_samples)
        pdf_vals = dist.pdf(samples)
        probs = pdf_vals # mc sampling
        # probs = pdf_vals / np.sum(pdf_vals)
        central_par = mean_vec.copy()
        print(f"  Central LM par = {central_par}")
        all_prob_central[i_pt] = dist.pdf(central_par)

        v2_central, v2_central_err, c_cent, fitter_cent = fit_with_lm_template(
                h_data, central_par, ry_val=ry_val, ry_err=ry_err,
                lm_ry_val=lm_ry_val, lm_ry_err=lm_ry_err, tag="central"
        )
        # Fix zero check
        all_v2_central[i_pt] = v2_central if v2_central is not None else 0.
        all_v2_err[i_pt] = v2_central_err if v2_central_err is not None else 0.
        print(f"  Central fit: v2 = {all_v2_central[i_pt]:.6f} ± {all_v2_err[i_pt]:.6f}")

        pt_dir = os.path.join(OUT_LM, f"PtCand_{pt_cand_idx}")
        trial_central = os.path.join(pt_dir, f"trial_{0:0{N_DIGITS}d}")
        lm_func_central = build_lm_tf1(central_par, "lm_central")
        if c_cent and fitter_cent:
                f_fit_cent = fitter_cent.GetFitFunction()
                f_lm_cent = fitter_cent.GetLMTemplateFunc()
                save_trial(trial_central, c_cent, fitter_cent, f_fit_cent, f_lm_cent,
                           pt_cand_idx, pt_had_idx, mass_label, "central",
                           pt_min, pt_max, PT_HAD[0], PT_HAD[1],
                           MASS_RANGE[0], MASS_RANGE[1])
                # Overwrite parLM.txt with original central par
                lm_names = ["Baseline", "Amplitude NS", "Sigma NS", "Amplitude AS", "Sigma AS"]
                with open(os.path.join(trial_central, "parLM.txt"), "w") as fp:
                    fp.write(f"# Trial central par  LM tempFunc={TEMP_FUNC}\n")
                    fp.write(f"{'Param':<6} {'Name':<20} {'Value':<18} {'Error':<18}\n")
                    fp.write("-"*60+"\n")
                    for ip in range(len(central_par)):
                        fp.write(f"{ip:<6} {lm_names[ip]:<20} {central_par[ip]:<18.8f} {0.0:<18.8f}\n")
                c_cent.Close(); del c_cent, fitter_cent
        gc.collect()

        v2_list, prob_list, fit_funcs, par_list_fit = [], [], [], []
        n_converged = 0
        with alive_bar(n_samples, title=f"Fitting PtCand {pt_cand_idx} trials") as bar:
            for isamp in range(n_samples):
                par = samples[isamp]
                if par[2] <= 0.05 or par[4] <= 0.05 or par[1] <= 0 or par[3] <= 0:
                    continue

                v2, _, c_var, fitter_var = fit_with_lm_template(
                    h_data, par, ry_val=ry_val, ry_err=ry_err,
                    lm_ry_val=lm_ry_val, lm_ry_err=lm_ry_err, tag=f"s{isamp}"
                )
                if v2 is not None:
                    v2_list.append(v2)
                    prob_list.append(probs[isamp])
                    n_converged += 1
                    fit_funcs.append(build_lm_tf1(par, f"lm_{isamp}"))
                    par_list_fit.append(par.copy())

                    trial_dir = os.path.join(pt_dir, f"trial_{n_converged:0{N_DIGITS}d}")
                    f_fit_var = fitter_var.GetFitFunction()
                    f_lm_var = fitter_var.GetLMTemplateFunc()
                    save_trial(trial_dir, c_var, fitter_var, f_fit_var, f_lm_var,
                                pt_cand_idx, pt_had_idx, mass_label, f"{n_converged}",
                                pt_min, pt_max, PT_HAD[0], PT_HAD[1],
                                MASS_RANGE[0], MASS_RANGE[1])
                    # Overwrite parLM.txt with original sampled par
                    with open(os.path.join(trial_dir, "parLM.txt"), "w") as fp:
                        fp.write(f"# Trial sampled par  LM tempFunc={TEMP_FUNC}\n")
                        fp.write(f"{'Param':<6} {'Name':<20} {'Value':<18} {'Error':<18}\n")
                        fp.write("-"*60+"\n")
                        _lm_names = ["Baseline", "Amplitude NS", "Sigma NS", "Amplitude AS", "Sigma AS"]
                        for ip in range(len(par)):
                            fp.write(f"{ip:<6} {_lm_names[ip]:<20} {par[ip]:<18.8f} {0.0:<18.8f}\n")
                    c_var.Close(); del c_var, fitter_var

                if isamp % 20 == 19:
                    ROOT.gROOT.GetListOfCanvases().Clear()
                    ROOT.gROOT.GetListOfFunctions().Clear()
                    ROOT.gDirectory = ROOT.nullptr
                    gc.collect()
                bar()

        all_v2_samples[i_pt] = np.array(v2_list)
        all_probs[i_pt] = np.array(prob_list)
        all_fit_funcs[i_pt] = (lm_func_central, fit_funcs)
        all_par_samples[i_pt] = par_list_fit
        all_lm_data[i_pt] = extract_lm_data_histo(lm_path)
        all_central_pars[i_pt] = central_par
        print(f"  Varied fits: {n_converged}/{n_samples} converged")
        if len(v2_list) > 0:
                print(f"  v2 range: [{np.min(v2_list):.6f}, {np.max(v2_list):.6f}]")
                print(f"  v2 mean ± RMS: {np.mean(v2_list):.6f} ± {np.std(v2_list):.6f}")

    # ═══════════════════════════════════════════════════════════════════
    # ── [OPTIMIZED] 新增通用散点图绘制函数，支持按概率热力图上色 ──
    # ═══════════════════════════════════════════════════════════════════
    def draw_colored_scatter(x_vals, p_vals, title, xtitle, ytitle, is_par_plot=False):
        """使用 TMultiGraph 和颜色分组实现逐点变色，高概率点置于顶层绘制。"""
        mg = ROOT.TMultiGraph()
        mg.SetTitle(f"{title};{xtitle};{ytitle}")
        
        p_min, p_max = np.min(p_vals), np.max(p_vals)
        n_colors = ROOT.TColor.GetNumberOfColors()
        
        # 按照颜色索引创建字典，将相同颜色的点归为一个 TGraph（极大提高 ROOT 绘制性能）
        graphs = {}
        
        # 获取按概率排序的索引，确保低概率的先画（底层），高概率的后画（顶层）
        sorted_indices = np.argsort(p_vals)
        
        for idx in sorted_indices:
            x, p = x_vals[idx], p_vals[idx]
            
            # 将概率映射到 0~1 的比例
            frac = (p - p_min) / (p_max - p_min) if p_max > p_min else 0.5
            c_idx = int(frac * (n_colors - 1))
            
            if c_idx not in graphs:
                graphs[c_idx] = ROOT.TGraph()
                graphs[c_idx].SetMarkerStyle(ROOT.kFullCircle)
                # 参数图点稍大，v2 散点稍小
                graphs[c_idx].SetMarkerSize(0.8 if is_par_plot else 0.6) 
                
                # 从 ROOT 色盘中获取当前概率对应的颜色
                color = ROOT.TColor.GetColorPalette(c_idx)
                alpha_color = ROOT.TColor.GetColorTransparent(color, 0.6)
                graphs[c_idx].SetMarkerColor(alpha_color)
                
            n_pts = graphs[c_idx].GetN()
            graphs[c_idx].SetPoint(n_pts, x, p)
            
        # 将字典按键值(颜色索引)排序加入，顺理成章实现高概率(暖色)在上的覆盖效果
        for c_idx in sorted(graphs.keys()):
            mg.Add(graphs[c_idx], "P")
            
        return mg, graphs  # 返回 graphs 防止 PyROOT 的垃圾回收机制清除底层对象


# ═══════════════════════════════════════════════════════════════════
    # ── 6-7. PLOT ────────────────────────────────────────────────────
    # ═══════════════════════════════════════════════════════════════════

    # ═════════════════════════════════════════════════════════════════════
    # ── Probability plots, Overlays, and Parameter plots ──
    # ═════════════════════════════════════════════════════════════════════
    f_out = ROOT.TFile.Open(os.path.join(OUT_LM, "probV2_PtCand.root"), "RECREATE")
    n_colors = ROOT.TColor.GetNumberOfColors()
    sigmas = []
    mus = []

    with alive_bar(n_pt_cand, title="Plotting Probability vs v2") as bar:
        for i_pt in range(n_pt_cand):
            if i_pt not in all_v2_samples or len(all_v2_samples[i_pt]) < 5: continue
            
            v2_arr = all_v2_samples[i_pt]
            prob_arr = all_probs[i_pt]
            p_min, p_max = np.min(prob_arr), np.max(prob_arr)

            d = f_out.GetDirectory(f"PtCand_{i_pt + 1}") or f_out.mkdir(f"PtCand_{i_pt + 1}")
            d.cd()

            v2_mean, v2_std = np.mean(v2_arr), np.std(v2_arr)
            v2_cent = all_v2_central[i_pt]
            prob_cent = all_prob_central[i_pt] if i_pt < len(all_prob_central) else p_max

            # ── 2. Probability vs v2 Scatter Plot (Colored 2D Graph) ──
            c2 = ROOT.TCanvas(f"c_pV2_PtCand{i_pt + 1}", f"Prob vs v2", 1000, 800)
            c2.SetLeftMargin(0.15); c2.SetRightMargin(0.05); c2.SetBottomMargin(0.12)

            mg_pv2, _h_refs_1 = draw_colored_scatter(
                v2_arr, prob_arr, 
                f"Probability vs v_{{2}}^{{#Delta}} (PtCand {PT_CAND[i_pt]}-{PT_CAND[i_pt+1]});v_{{2}}^{{#Delta}};Probability density",
                "v_{2}^{#Delta}", "Probability density", is_par_plot=False
            )
            
            mg_pv2.Draw("A")
            ROOT.gPad.Update()
            
            y_min = p_min
            y_max = np.max([p_max, prob_cent]) * 1.1 
            y_bottom = y_min * 0.9 if y_min > 0 else 0
            mg_pv2.GetYaxis().SetRangeUser(y_bottom, y_max)
            
            # 1 sigma 辅助线
            line_mean = ROOT.TLine(v2_mean, y_bottom, v2_mean, y_max * 0.9)
            line_mean.SetLineColor(ROOT.kGreen + 1); line_mean.SetLineStyle(2); line_mean.Draw()
            line_std_low = ROOT.TLine(v2_mean - v2_std, y_bottom, v2_mean - v2_std, y_max * 0.9)
            line_std_low.SetLineColor(ROOT.kGreen + 1); line_std_low.SetLineStyle(3); line_std_low.Draw()
            line_std_high = ROOT.TLine(v2_mean + v2_std, y_bottom, v2_mean + v2_std, y_max * 0.9)
            line_std_high.SetLineColor(ROOT.kGreen + 1); line_std_high.SetLineStyle(3); line_std_high.Draw()

            # 中心值标注 (黑色虚线 + 红色五角星)
            line = ROOT.TLine(v2_cent, y_bottom, v2_cent, prob_cent)
            line.SetLineColor(ROOT.kBlack); line.SetLineStyle(5); line.Draw()
            
            gr_cent_v2 = ROOT.TGraph(1)
            gr_cent_v2.SetPoint(0, v2_cent, prob_cent)
            gr_cent_v2.SetMarkerStyle(ROOT.kFullStar)
            gr_cent_v2.SetMarkerSize(2.0) 
            gr_cent_v2.SetMarkerColor(ROOT.kRed + 1)
            gr_cent_v2.Draw("P SAME")

            c2.SaveAs(os.path.join(OUT_LM, f"probV2_PtCand{i_pt + 1}.png"))
            c2.Write(); c2.Close()

            # ── 3. [恢复] v2 概率密度直方图 (1D TH1D) ──
            c_h_pv2 = ROOT.TCanvas(f"c_hProbV2_PtCand{i_pt + 1}", f"1D Probability vs v2", 1000, 800)
            c_h_pv2.SetLeftMargin(0.15); c_h_pv2.SetRightMargin(0.05); c_h_pv2.SetBottomMargin(0.12)
            
            v2_min_arr, v2_max_arr = np.min(v2_arr), np.max(v2_arr)
            if v2_min_arr == v2_max_arr:
                pad = max(abs(v2_min_arr) * 0.1, 1e-6)
                v2_min_arr, v2_max_arr = v2_min_arr - pad, v2_max_arr + pad

            # concentrate bins around the central value for better visualization
            diffs = np.abs(v2_arr - v2_cent)
            half_width = np.max(diffs)
            
            # 增加 5% 留白，防止最边缘的点落在边界上导致 bin 计算误差
            half_width *= 1.05
            
            v2_min_sym = v2_cent - half_width
            v2_max_sym = v2_cent + half_width
            
            # 2. 设定奇数个 bin (例如 25 个)，确保中心 bin 的中心严格对齐 v2_cent
            n_bins = 50 
            
            h_pv2 = ROOT.TH1D(f"hProbV2_PtCand{i_pt + 1}", 
                            f"Distribution of v_{{2}}^{{#Delta}} (PtCand {PT_CAND[i_pt]}-{PT_CAND[i_pt + 1]});v_{{2}}^{{#Delta}};Counts",
                            n_bins, v2_min_sym, v2_max_sym)
            for v in v2_arr:
                h_pv2.Fill(v)
            
            # h_pv2.SetLineColor(ROOT.kBlue + 1)
            # h_pv2.SetLineWidth(2)
            h_pv2.SetFillColorAlpha(ROOT.kBlue - 7, 0.5)
            h_pv2.Draw("HIST")
            
            
            # 直方图上的均值红虚线
            l_mean_hist = ROOT.TLine(v2_mean, 0, v2_mean, h_pv2.GetMaximum() * 1.05)
            l_mean_hist.SetLineColor(ROOT.kRed + 1)
            l_mean_hist.SetLineStyle(2)
            l_mean_hist.SetLineWidth(2)
            l_mean_hist.Draw("SAME")
            # 1σ 范围虚线
            l_std_low_hist = ROOT.TLine(v2_mean - v2_std, 0, v2_mean - v2_std, h_pv2.GetMaximum() * 1.05)
            l_std_low_hist.SetLineColor(ROOT.kRed + 1)
            l_std_low_hist.SetLineStyle(5)  # dashed line style
            l_std_low_hist.SetLineWidth(2)
            l_std_low_hist.Draw("SAME")
            l_std_high_hist = ROOT.TLine(v2_mean + v2_std, 0, v2_mean + v2_std, h_pv2.GetMaximum() * 1.05)
            l_std_high_hist.SetLineColor(ROOT.kRed + 1)
            l_std_high_hist.SetLineStyle(5)  # dashed line style
            l_std_high_hist.SetLineWidth(2)
            l_std_high_hist.Draw("SAME")
            # 3σ 范围虚线
            l_std_low = ROOT.TLine(v2_mean - v2_std, 0, v2_mean - v2_std, h_pv2.GetMaximum() * 1.05)
            l_std_low.SetLineColor(ROOT.kRed + 1)
            l_std_low.SetLineStyle(3)
            l_std_low.SetLineWidth(2)
            l_std_low.Draw("SAME")
            l_std_high = ROOT.TLine(v2_mean + v2_std, 0, v2_mean + v2_std, h_pv2.GetMaximum() * 1.05)
            l_std_high.SetLineColor(ROOT.kRed + 1)
            l_std_high.SetLineStyle(3)
            l_std_high.SetLineWidth(2)
            l_std_high.Draw("SAME")

            line_cent = ROOT.TLine(v2_cent, 0, v2_cent, h_pv2.GetMaximum() * 1.1)
            line_cent.SetLineColor(ROOT.kGreen + 1)
            line_cent.SetLineWidth(2)
            line_cent.SetLineStyle(2)
            line_cent.Draw("SAME")

            text_mean_cent = ROOT.TLatex()
            text_mean_cent.SetNDC()
            text_mean_cent.SetTextSize(0.035)
            text_mean_cent.SetTextColor(ROOT.kRed + 1)
            text_mean_cent.DrawLatex(0.15, 0.85, f"Mean = {v2_mean:.6f} #pm {v2_std:.6f}")
            text_mean_cent.SetTextColor(ROOT.kGreen + 1)
            text_mean_cent.DrawLatex(0.15, 0.80, f"Central = {v2_cent:.6f}")
            
            fit_res = h_pv2.Fit("gaus", "LSQ+")
            fit_func = h_pv2.GetFunction("gaus")

            if fit_func:
                mu = fit_func.GetParameter(1)
                sigma = fit_func.GetParameter(2)
                mus.append(mu)
                sigmas.append(sigma)
                text_fit = ROOT.TLatex()
                text_fit.SetNDC()
                text_fit.SetTextSize(0.035)
                text_fit.SetTextColor(ROOT.kMagenta + 1) # 把字体也改成和线一样的颜色，更直观
                text_fit.DrawLatex(0.20, 0.75, f"#mu = {mu:.6f}")
                text_fit.DrawLatex(0.20, 0.70, f"#sigma = {sigma:.6f}")
                
                cloned_fit_func = fit_func.Clone(f"fit_gaus_PtCand{i_pt + 1}")
                cloned_fit_func.SetLineColor(ROOT.kMagenta + 1)
                cloned_fit_func.SetLineWidth(2)
                # 显式地将克隆的函数画在当前画布上
                cloned_fit_func.Draw("SAME")
            c_h_pv2.Update()
            h_pv2.Write()
            c_h_pv2.SaveAs(os.path.join(OUT_LM, f"hProbV2_PtCand{i_pt + 1}.png"))

            c_h_pv2.Write(f"c_hProbV2_PtCand{i_pt + 1}")
            c_h_pv2.Close()

            # ── 4. LM Fit Functions Overlay ──
            lm_central, lm_funcs = all_fit_funcs.get(i_pt, (None, []))
            h_lm_data = all_lm_data.get(i_pt)
            if h_lm_data and lm_central and lm_funcs:
                c_lm = ROOT.TCanvas(f"c_LMFits_PtCand{i_pt + 1}", "LM fits", 1400, 900)
                c_lm.SetLeftMargin(0.12); c_lm.SetRightMargin(0.03)
                c_lm.SetBottomMargin(0.10); c_lm.SetTopMargin(0.05)
                
                h_lm_data.SetStats(0)
                h_lm_data.SetMarkerStyle(ROOT.kFullCircle)
                h_lm_data.SetMarkerSize(1.2)
                h_lm_data.SetMarkerColor(ROOT.kBlack)
                h_lm_data.SetLineColor(ROOT.kBlack)
                h_lm_data.Draw("E")
                
                sorted_indices = np.argsort(prob_arr)
                total_draw = len(lm_funcs)
                base_alpha = 0.05 if total_draw > 500 else (0.15 if total_draw > 100 else 0.35)
                
                for idx in sorted_indices:
                    f = lm_funcs[idx]
                    p = prob_arr[idx]
                    frac = (p - p_min) / (p_max - p_min) if p_max > p_min else 0.5
                    c_idx = int(frac * (n_colors - 1))
                    
                    color = ROOT.TColor.GetColorPalette(c_idx)
                    alpha_idx = ROOT.TColor.GetColorTransparent(color, base_alpha * (0.5 + 0.5 * frac)) 
                    f.SetLineColor(alpha_idx)
                    f.SetLineWidth(1)
                    f.Draw("same")

                lm_central.SetLineColor(ROOT.kBlack)
                lm_central.SetLineWidth(3)
                lm_central.Draw("same")  
                h_lm_data.Draw("same E") 
                
                c_lm.Update()
                out_lm_png = os.path.join(OUT_LM, "LMSys", f"LMFits_PtCand{i_pt + 1}.png")
                os.makedirs(os.path.dirname(out_lm_png), exist_ok=True)
                
                lm_central.Write(f"lm_central_PtCand{i_pt + 1}")
                for i, f in enumerate(lm_funcs[:50]):
                    f.Write(f"lm_func_{i}_PtCand{i_pt + 1}")
                c_lm.SaveAs(out_lm_png)
                c_lm.Write(f"c_LMFits_PtCand{i_pt + 1}")
                c_lm.Close()
            bar()  # Finish progress bar
        f_out.cd()
    f_out.Close()

    print(f"\n[PLOT] Summaries and Systematics...")
    pt_centers = np.array([0.5 * (PT_CAND[i] + PT_CAND[i + 1]) for i in range(n_pt_cand)])

    # ── 1. v2 vs pT Systematic Band ──
    c1 = ROOT.TCanvas("c_v2_vs_pt", "v2 vs pT with LM sys", 1400, 900)
    c1.SetLeftMargin(0.15); c1.SetRightMargin(0.05); c1.SetBottomMargin(0.12)
    gr_cent = ROOT.TGraphErrors(n_pt_cand)
    gr_cent.SetName("gV2Central")
    gr_cent.SetTitle(";p_{T}^{D} (GeV/c);v_{2}^{#Delta}")
    gr_cent.SetMarkerStyle(ROOT.kFullCircle)
    gr_cent.SetMarkerColor(ROOT.kBlack)
    gr_cent.SetLineColor(ROOT.kBlack)

    gr_band = ROOT.TGraphAsymmErrors(n_pt_cand)
    gr_band.SetName("gV2SysBand")
    gr_band.SetFillColorAlpha(ROOT.kBlue - 7, 0.35)
    # gr_band.SetLineColor(ROOT.kBlue + 2)

    for i_pt in range(n_pt_cand):
        pt_c = pt_centers[i_pt]
        gr_cent.SetPoint(i_pt, pt_c, all_v2_central[i_pt])
        gr_cent.SetPointError(i_pt, 0.5 * (PT_CAND[i_pt + 1] - PT_CAND[i_pt]), all_v2_err[i_pt])

        if i_pt in all_v2_samples and len(all_v2_samples[i_pt]) > 0:
            v2_arr = all_v2_samples[i_pt]
            v2_mean, v2_std = np.mean(v2_arr), np.std(v2_arr)
            gr_band.SetPoint(i_pt, pt_c, v2_mean)
            gr_band.SetPointError(i_pt,
                    0.3 * (PT_CAND[i_pt + 1] - PT_CAND[i_pt]),
                    0.3 * (PT_CAND[i_pt + 1] - PT_CAND[i_pt]),
                    v2_std, v2_std)
        else:
            gr_band.SetPoint(i_pt, pt_c, all_v2_central[i_pt])
            gr_band.SetPointError(i_pt, 0., 0., 0., 0.)
    
    # band from fit
    gr_band_fit = ROOT.TGraphAsymmErrors(n_pt_cand)
    gr_band_fit.SetName("gV2SysBandFit")
    gr_band_fit.SetFillColorAlpha(ROOT.kRed - 7, 0.35)
    # gr_band_fit.SetLineColor(ROOT.kRed + 2)
    for i_pt in range(n_pt_cand):
        pt_c = pt_centers[i_pt]
        mu = mus[i_pt]
        sigma = sigmas[i_pt]
        gr_band_fit.SetPoint(i_pt, pt_c, mu)
        gr_band_fit.SetPointError(i_pt,
                0.3 * (PT_CAND[i_pt + 1] - PT_CAND[i_pt]),
                0.3 * (PT_CAND[i_pt + 1] - PT_CAND[i_pt]),
                sigma, sigma)

    gr_band.Draw("A2")
    gr_band.GetYaxis().SetRangeUser(0., max(all_v2_central.max() + 0.01, 0.01))
    gr_band_fit.Draw("2 SAME")
    gr_cent.Draw("P SAME")

    leg = ROOT.TLegend(0.65, 0.72, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0)
    leg.AddEntry(gr_cent, "Central fit", "lep")
    leg.AddEntry(gr_band, "Mean #pm1# std sigma band", "f")
    leg.AddEntry(gr_band_fit, "Gaussian #pm1#sigma band", "f")
    leg.Draw()

    c1.SaveAs(os.path.join(OUT_LM, "v2_vs_pT_LMSys.png"))
    out_file = ROOT.TFile.Open(os.path.join(OUT_LM, "v2_vs_pT_LMSys.root"), "RECREATE")
    gr_cent.Write(); gr_band.Write(); gr_band_fit.Write(); c1.Write()
    out_file.Close()

    # ── 5. Parameter-Probability Scatter Plots (Colored) ──
    print(f"\n[INFO] Parameter-probability scatter plots...")
    par_names = ["Baseline", "Amplitude NS", "Sigma NS", "Amplitude AS", "Sigma AS"]
    
    for i_pt in range(n_pt_cand):
        if i_pt not in all_par_samples or not all_par_samples[i_pt]: continue
        
        pars = all_par_samples[i_pt]
        probs_arr = all_probs.get(i_pt, [])
        p_max = np.max(probs_arr) if len(probs_arr) > 0 else 1.0
        prob_cent = all_prob_central[i_pt] if i_pt < len(all_prob_central) else p_max
        
        if i_pt in all_central_pars and all_central_pars[i_pt] is not None:
            curr_central_par = all_central_pars[i_pt]
        else:
            curr_central_par = None

        for ipar in range(5):
            c_pp = ROOT.TCanvas(f"c_ParProb_PtCand{i_pt+1}_par{ipar}", f"{par_names[ipar]} vs Prob", 800, 600)
            c_pp.SetLeftMargin(0.15); c_pp.SetRightMargin(0.05); c_pp.SetBottomMargin(0.12)
            
            x_vals = np.array([p[ipar] for p in pars])
            x_min, x_max = np.min(x_vals), np.max(x_vals)
            
            if x_min == x_max:
                pad = max(abs(x_min) * 0.1, 1e-6)
                x_min, x_max = x_min - pad, x_max + pad
                
            mg_pp, _h_refs_2 = draw_colored_scatter(
                x_vals, probs_arr,
                f"{par_names[ipar]} vs Prob (PtCand {PT_CAND[i_pt]}-{PT_CAND[i_pt+1]});{par_names[ipar]};Probability density",
                par_names[ipar], "Probability density", is_par_plot=True
            )
            
            mg_pp.Draw("A")
            
            margin = (x_max - x_min) * 0.1
            mg_pp.GetXaxis().SetLimits(x_min - margin, x_max + margin)
            
            y_max_ax = np.max([p_max, prob_cent]) * 1.1
            mg_pp.GetYaxis().SetRangeUser(0, y_max_ax)
            ROOT.gPad.Update()
            
            cent_val = curr_central_par[ipar] if curr_central_par is not None else np.mean(x_vals)
            
            line = ROOT.TLine(cent_val, 0, cent_val, prob_cent)
            line.SetLineColor(ROOT.kBlack); line.SetLineStyle(5); line.Draw()
            
            gr_cent_p = ROOT.TGraph(1)
            gr_cent_p.SetPoint(0, cent_val, prob_cent)
            gr_cent_p.SetMarkerStyle(ROOT.kFullStar)
            gr_cent_p.SetMarkerSize(2.0)
            gr_cent_p.SetMarkerColor(ROOT.kRed + 1)
            gr_cent_p.Draw("P SAME")

            c_pp.Update()
            
            out_png = os.path.join(OUT_LM, "LMSys", f"ParProb_{par_names[ipar].replace(' ','_')}_PtCand{i_pt+1}.png")
            os.makedirs(os.path.dirname(out_png), exist_ok=True)
            c_pp.SaveAs(out_png)
            c_pp.Close()

    print(f"\n[DONE] Outputs in: {OUT_LM}")
    print(f"  v2_vs_pT_LMSys.png  — v2 vs pT with systematics band")
    print(f"  probV2_PtCand*.png  — probability vs v2 per pT bin (Colored & Sorted)")
    print(f"  hProbV2_PtCand*.png — 1D counts distribution of v2")


# ═══════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LM fit systematic uncertainty via parameter sampling")
    parser.add_argument("--n-samples", type=int, default=N_SAMPLES, help="Number of samples")
    parser.add_argument("--pt-only", type=int, default=0, help="Process only one PtCand bin")
    parser.add_argument("--load-from-exist", '-l', action="store_true", help="Load v2 from exist")
    args = parser.parse_args()
    run_systematics(n_samples=args.n_samples, load_from_exist=args.load_from_exist, pt_only=args.pt_only)