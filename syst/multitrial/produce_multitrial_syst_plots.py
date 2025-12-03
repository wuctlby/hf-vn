import sys
import os
import numpy as np
import argparse
import re
import glob
import array
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kError  # Only show errors and above
from ROOT import TFile, TCanvas, TH1F, TGraphAsymmErrors, TLegend, kOrange, kAzure, kBlack
script_dir = os.path.dirname(os.path.realpath(__file__))
os.sys.path.append(os.path.join(script_dir, '../..', 'utils'))
from utils import logger
from StyleFormatter import SetGlobalStyle, SetObjectStyle

SetGlobalStyle(titleoffsety=1.1, maxdigits=3, topmargin=0.1, bottommargin=0.4, leftmargin=0.3, rightmargin=0.15,
               labelsizey=0.04, setoptstat=0, setopttitle=0, setdecimals=True,titleoffsetx=0.74)

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

    trials_dirs = glob.glob(f"{args.results_dir}/syst/multitrial/{mult_dir}/trial_*")

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
    syst_uncs['PromptV2'] = compute_systematics_fin_val(f"{out_dir}/{mult_dir}/SystV2Prompt.root", "Prompt #it{v}_{2}", \
                                   v2_prompt_ref, v2_prompt_ref_unc, \
                                   v2_prompt_trials, v2_prompt_trials_uncs)

    # FD v2 systematic uncertainty wrt reference
    v2_non_prompt_ref = v2_vs_frac_ref.Get('hV2VsPtFD').GetBinContent(pt_bin_ref)
    v2_non_prompt_ref_unc = v2_vs_frac_ref.Get('hV2VsPtFD').GetBinError(pt_bin_ref)
    syst_uncs['NonPromptV2'] = compute_systematics_fin_val(f"{out_dir}/{mult_dir}/SystV2FD.root", "Non-prompt #it{v}_{2}", \
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
    parser.add_argument('results_dir', metavar='text', default='path to the directory containing the .root files from multitrial')
    args = parser.parse_args()

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
        _, pt_min_times_10, pt_max_times_10 = mult_dir.split('_')
        pt_min = float(pt_min_times_10) / 10.
        pt_max = float(pt_max_times_10) / 10.
        if pt_min not in pt_bins:
            pt_bins.append(pt_min)
        if pt_max not in pt_bins:
            pt_bins.append(pt_max)

    pt_bins = sorted(pt_bins)
    # Produce plots with all pt-bins summary for prompt and non-prompt v2
    h_prompt_syst = TH1F('h_prompt_syst', 'h_prompt_syst;#it{p}_{T} (GeV/#it{c});syst (%)', len(pt_bins)-1, array.array('d', pt_bins))
    h_non_prompt_syst = TH1F('h_non_prompt_syst', 'h_non_prompt_syst;#it{p}_{T} (GeV/#it{c});syst (%)', len(pt_bins)-1, array.array('d', pt_bins))
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
