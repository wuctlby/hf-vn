import argparse
import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import yaml
import ROOT
import array
os.environ["CUDA_VISIBLE_DEVICES"] = "" # pylint: disable=wrong-import-position
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
from ROOT import TFile, TH1F, TGraphAsymmErrors, kBlack, kFullCircle
script_dir = os.path.dirname(os.path.realpath(__file__))
os.sys.path.append(os.path.join(script_dir, '..', 'utils'))
from utils import logger
from load_utils import load_aod_file
from StyleFormatter import SetGlobalStyle, SetObjectStyle
from matplotlib import gridspec
ROOT.gROOT.SetBatch(True)
import multiprocessing as mp
import sys
sys.path.append("./flareflyfitter/")
from raw_yield_fitter import RawYieldFitter
mp.set_start_method("spawn", force=True)
msg_service = ROOT.RooMsgService.instance()

def process_pt_bin(i_cutset, i_pt, config, summary, rebin_factor, sp_edges, pt_label, data, reso, outdir, sel_string):

    # Initialize the fitter for sp-integrated yield extraction
    fit_cfg = config['v2extraction']
    fitter = RawYieldFitter(config['Dmeson'], config['ptbins'][i_pt], config['ptbins'][i_pt + 1],
                            f"sp_integrated_{pt_label}_fit", config['V2ExtractionByYield']['Minimizer'])
    fitter.set_rebin(fit_cfg['Rebin'][i_pt] if fit_cfg.get('Rebin') else 1)
    fitter.set_fit_range(fit_cfg['MassFitRanges'][i_pt][0], fit_cfg['MassFitRanges'][i_pt][1])
    if isinstance(data, ROOT.TH2):
        data_mass_int = data.ProjectionX()
        data_mass_int = data
        fitter.set_data_to_fit_hist(data_mass_int)
    else:
        data = data.query(sel_string.replace(" && ", " and "))
        fitter.set_data_to_fit_df(data, 'fM')

    # Add model components
    fitter.add_bkg_func(fit_cfg['BkgFunc'][i_pt] if isinstance(fit_cfg['BkgFunc'], list) else fit_cfg['BkgFunc'], "Comb. bkg")
    sgn_funcs = {} # More info for signal functions, a dictionary is better
    sgn_funcs[fit_cfg['SgnFuncLabel']] = {
        'func': fit_cfg['SgnFunc'][i_pt] if isinstance(fit_cfg['SgnFunc'], list) else fit_cfg['SgnFunc'],
        'part': config['Dmeson']
    }
    if fit_cfg.get('InclSecPeak'):
        include_sec_peak = fit_cfg['InclSecPeak'][i_pt] if isinstance(fit_cfg['InclSecPeak'], list) else fit_cfg['InclSecPeak']
        if include_sec_peak:
            sgn_funcs[fit_cfg['SgnFuncSecPeakLabel']] = {
                'func': fit_cfg['SgnFuncSecPeak'][i_pt] if isinstance(fit_cfg['SgnFuncSecPeak'], list) else fit_cfg['SgnFuncSecPeak'],
                'part': 'Dplus' if config['Dmeson'] == 'Ds' else 'Dstar',
            }
    for i_sgn, (label, sgn_func) in enumerate(sgn_funcs.items()):
        fitter.add_sgn_func(sgn_func['func'], label, sgn_func['part'])

    # Add correlated background if specified
    if config.get('corr_bkgs'):
        fitter.add_corr_bkgs(config['corr_bkgs'], sel_string, config['ptbins'][i_pt], config['ptbins'][i_pt + 1])

    fitter.setup()
    if fit_cfg.get('InitPars'):
        fitter.set_fit_pars(fit_cfg['InitPars'], pt_min, pt_max)

    # Prefit the MC prompt enhanced cut to fix the tails, binned fit
    if fit_cfg.get('FixSgnFromMC'):
        fitter.set_fix_sgn_to_mc_prefit(True)
        fitter.prefit_mc(f"{config['outdir']}/corrbkgs/templs_{pt_label}.root")
        fitter.plot_mc_prefit(False, True, loc=["lower left", "upper left"],
                              path=f"{outdir}/")
        fitter.plot_raw_residuals_mc_prefit(path=f"{outdir}/fM_mc_prefit_residuals_{i_cutset}_{pt_label}.pdf")

    if config['V2ExtractionByYield'].get('FixParsToIntFit'):  # Fix mean to results of sp-integrated fit
        logger(f"Will fix all parameters of signal functions to sp-integrated fit result", "WARNING")
        fitter.fix_sgn_pars_to_first_fit()
    status, converged = fitter.fit()
    fig_int_fit_path = f"{outdir}/fit_int_{pt_label}.pdf"
    fitter.plot_fit(False, True, loc=["lower left", "upper left"],  # (log, show_extra_info)
                    path=fig_int_fit_path)

    fit_info, sgn_pars, sgn_pars_uncs, _, _ = fitter.get_fit_info()
    label = list(sgn_funcs.keys())[0]       # Take only the first signal function, which is the peak of interest
    raw_yield, raw_yield_unc = fit_info[label]["ry"], fit_info[label]["ry_unc"]
    mean_int, mean_unc = sgn_pars[f"mu_{label}"], sgn_pars_uncs[f"mu_{label}"]
    sigma_int, sigma_unc = sgn_pars[f"sigma_{label}"], sgn_pars_uncs[f"sigma_{label}"]

    summary['hRawYieldsSimFit'].SetBinContent(i_pt + 1, raw_yield)
    summary['hRawYieldsSimFit'].SetBinError(i_pt + 1, raw_yield_unc)
    summary['hMeanSimFit'].SetBinContent(i_pt + 1, mean_int)
    summary['hMeanSimFit'].SetBinError(i_pt + 1, mean_unc)
    summary['hSigmaSimFit'].SetBinContent(i_pt + 1, sigma_int)
    summary['hSigmaSimFit'].SetBinError(i_pt + 1, sigma_unc)

    # Extract raw yields in each sp bin
    stats = [
        'means', 'means_unc',
        'sigmas', 'sigmas_unc',
        'sp_ry', 'sp_ry_unc',
        'weight_av_unc', 'weight_av_unc_first_term', 'weight_av_unc_second_term',
    ]
    hist_stats, vals_stats = {}, {}
    for stat in stats:
        hist_stats[stat] = ROOT.TH1F(f"hist_{stat}", f"hist_{stat}", len(sp_edges['bins']) - 1, array.array('f', sp_edges['sp_mins'] + [sp_edges['sp_maxs'][-1]]))
        hist_stats[stat].SetDirectory(0)
        vals_stats[stat] = [0] * (len(sp_edges['bins']) - 1)
    vals_stats['sp_center'] = [0] * (len(sp_edges['bins']) - 1)
    if isinstance(data, ROOT.TH2):
        hist_stats['sp_histos'] = [None] * (len(sp_edges['bins']) - 1)

    for isp, (sp_left_bin, sp_right_bin, sp_min, sp_max) in enumerate(zip(sp_edges['bins'][:-1], sp_edges['bins'][1:], \
                                                                          sp_edges['sp_mins'][:-1], sp_edges['sp_maxs'][:-1])):
        sp_center = (sp_min + sp_max) / 2
        sp_label = f"sp_bins_{sp_left_bin}_{sp_right_bin-1}_range_{sp_min:.2f}_{sp_max:.2f}"
        logger(f"Processing sp_label {sp_label}", "INFO")

        # For each sp bin, update fitter changing the histogram to fit, name and rebin
        if isinstance(data, ROOT.TH2):
            data.GetYaxis().SetRange(sp_left_bin, sp_right_bin)
            data_sp = data.ProjectionX(f"h_mass_{sp_label}")
            data_sp.SetDirectory(0)
            fitter.set_rebin(rebin_factor)
            fitter.set_data_to_fit_hist(data_sp)
        else:
            data_sp = data.query(f"fScalarProd >= {sp_min} and fScalarProd < {sp_max}")
            fitter.set_data_to_fit_df(data_sp, 'fM')

        fitter.set_name(sp_label)
        fitter.setup()
        if isinstance(data, ROOT.TH2):
            hist_stats['sp_histos'][isp] = histo_sp
        vals_stats['sp_center'][isp] = sp_center
        try:
            status, converged = fitter.fit()

            fig_sp_pt_path = f"{outdir}/{pt_label}/{sp_label}.pdf"
            fitter.plot_fit(False, True, loc=["lower left", "upper left"], # (log, show_extra_info)
                            path=fig_sp_pt_path)
            fit_info, sgn_pars, sgn_pars_uncs, _, _ = fitter.get_fit_info()
            label = list(sgn_funcs.keys())[0]
            vals_stats['sp_ry'][isp], vals_stats['sp_ry_unc'][isp] = fit_info[label]["ry"], fit_info[label]["ry_unc"]
            vals_stats['means'][isp], vals_stats['means_unc'][isp] = sgn_pars[f"mu_{label}"], sgn_pars_uncs[f"mu_{label}"]
            vals_stats['sigmas'][isp], vals_stats['sigmas_unc'][isp] = sgn_pars[f"sigma_{label}"], sgn_pars_uncs[f"sigma_{label}"]
        except Exception as e:
            logger(f"Error fitting sp {sp_min:.2f} - {sp_max:.2f}: {e}", "ERROR")
            vals_stats['sp_ry'][isp], vals_stats['sp_ry_unc'][isp] = 0, 0
            vals_stats['means'][isp], vals_stats['means_unc'][isp] = 0, 0
            vals_stats['sigmas'][isp], vals_stats['sigmas_unc'][isp] = 0, 0

    # Compute weighted average and fill histograms
    sum_yields = np.sum(vals_stats['sp_ry'])
    sum_yields_unc = np.sqrt(np.sum(np.array(vals_stats['sp_ry_unc'])**2))
    summary["hSummedSpYields"].SetBinContent(i_pt + 1, sum_yields)
    summary["hSummedSpYields"].SetBinError(i_pt + 1, sum_yields_unc)
    sum_weighted = np.sum(np.array(vals_stats['sp_ry']) * (np.array(vals_stats['sp_center']) / reso))
    weighted_avg = sum_weighted / sum_yields
    weighted_avg_unc = 0
    for i_sp, (i_v2, i_ry, i_ry_unc) in enumerate(zip(vals_stats['sp_center'], vals_stats['sp_ry'], vals_stats['sp_ry_unc'])):
        d_yields_dRyi = (i_v2/reso) / sum_yields
        d_v2_dRyi = (sum_weighted) / (sum_yields**2)
        weighted_avg_unc = weighted_avg_unc + ((d_yields_dRyi - d_v2_dRyi)**2)*i_ry_unc**2

        vals_stats['weight_av_unc_first_term'][i_sp] = (d_yields_dRyi)**2 * i_ry_unc**2
        hist_stats['weight_av_unc_first_term'].SetBinContent(i_sp + 1, vals_stats['weight_av_unc_first_term'][i_sp])
        vals_stats['weight_av_unc_second_term'][i_sp] = (d_v2_dRyi)**2 * i_ry_unc**2
        hist_stats['weight_av_unc_second_term'].SetBinContent(i_sp + 1, vals_stats['weight_av_unc_second_term'][i_sp])
        hist_stats['weight_av_unc'].SetBinContent(i_sp + 1, np.sqrt(vals_stats['weight_av_unc_first_term'][i_sp] + \
                                                                    vals_stats['weight_av_unc_second_term'][i_sp]))

        for stat in stats:
            if 'unc' in stat:
                continue
            hist_stats[stat].SetBinContent(i_sp + 1, vals_stats[stat][i_sp])
            hist_stats[stat].SetBinError(i_sp + 1, vals_stats[stat + '_unc'][i_sp])
            hist_stats[stat + '_unc'].SetBinContent(i_sp + 1, vals_stats[stat + '_unc'][i_sp])

    weighted_avg_unc = np.sqrt(weighted_avg_unc)

    summary['hVnSimFit'].SetBinContent(i_pt + 1, weighted_avg)
    summary['hVnSimFit'].SetBinError(i_pt + 1, weighted_avg_unc)
    pt_min, pt_max = config['ptbins'][i_pt], config['ptbins'][i_pt+1]
    summary['gVnSimFit'].SetPoint(i_pt, (pt_min+pt_max)/2, weighted_avg)
    summary['gVnSimFit'].SetPointError(i_pt, (pt_max-pt_min)/2, (pt_max-pt_min)/2, weighted_avg_unc, weighted_avg_unc)
    summary['gVnUnc'].SetPoint(i_pt, (pt_min+pt_max)/2, weighted_avg_unc)
    summary['gVnUnc'].SetPointError(i_pt, (pt_max-pt_min)/2, (pt_max-pt_min)/2, 1.e-20, 1.e-20)

    return hist_stats, weighted_avg, weighted_avg_unc

def get_sp_bin_edges(axis_sp, cfg, i_cut, i_pt):

    sp_abs_val_max = cfg['SpRanges'][i_cut]
    if sp_abs_val_max <= 0:
        raise ValueError(f"Invalid sp_abs_val_max: {sp_abs_val_max}. It must be positive.")
        sys.exit(1)

    SpWindowsNbins = cfg['SpWindowsNbins'][i_cut][i_pt]
    sp_bins = axis_sp.GetNbins()

    # Find first bin
    for i in range(1, sp_bins + 1):
        if axis_sp.GetBinLowEdge(i)-0.0001 > -sp_abs_val_max:
            logger(f"First bin found at i={i}, low edge={axis_sp.GetBinLowEdge(i)}", "WARNING")
            first_bin = i-1
            break

    # Find last bin
    for i in range(1, sp_bins + 1):
        if axis_sp.GetBinLowEdge(i)-0.0001 > sp_abs_val_max:
            logger(f"Last bin found at i={i}, low edge={axis_sp.GetBinLowEdge(i)}", "WARNING")
            last_bin = i-1
            break
    
    edges = {}
    edges['bins'] = list(range(first_bin, last_bin, SpWindowsNbins))
    edges['sp_mins'] = [axis_sp.GetBinLowEdge(b) for b in edges['bins']]
    edges['sp_maxs'] = [axis_sp.GetBinUpEdge(b + SpWindowsNbins - 1) for b in edges['bins']]
    return edges

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('input_config', metavar='text', default='config.yml')
    parser.add_argument('infile', metavar='text', default='proj_XX.root')
    parser.add_argument('--batch', '-b', help='suppress video output', action='store_true')
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)

    with open(args.input_config, 'r') as CfgFlow:
        cfg_flow = yaml.safe_load(CfgFlow)

    # Set outfile name
    out_file_name = os.path.join(os.path.dirname(os.path.dirname(args.infile)),
                                 'raw_yields',
                                 os.path.basename(args.infile).replace('proj', 'raw_yields'))
    os.makedirs(os.path.dirname(out_file_name), exist_ok=True)

    pt_bins = cfg_flow['ptbins']
    hists_summary = {
        'hRawYieldsSimFit': TH1F("hRawYieldsSimFit", "hRawYieldsSimFit", len(pt_bins) - 1, np.array(pt_bins)),
        'hSummedSpYields': TH1F("hSummedSpYields", "hSummedSpYields", len(pt_bins) - 1, np.array(pt_bins)),
        'hMeanSimFit': TH1F("hMeanSimFit", "hMeanSimFit", len(pt_bins) - 1, np.array(pt_bins)),
        'hSigmaSimFit': TH1F("hSigmaSimFit", "hSigmaSimFit", len(pt_bins) - 1, np.array(pt_bins)),
        'hVnSimFit': TH1F("hVnSimFit", "hVnSimFit", len(pt_bins) - 1, np.array(pt_bins)),
        'gVnSimFit': TGraphAsymmErrors(1),
        'gVnUnc': TGraphAsymmErrors(1),
    }
    hists_summary['gVnSimFit'].SetName("gVnSimFit")
    hists_summary['gVnUnc'].SetName("gVnUnc")

    pt_bins = cfg_flow['ptbins']
    
    # Retrieve number of cutset
    _, cutset_str = os.path.splitext(os.path.basename(args.infile))[0].split('_')
    i_cutset = int(cutset_str)
    with open((args.infile).replace('proj', 'cutset').replace('.root', '.yml'), 'r') as CfgFlow:
        cfg_cutset = yaml.safe_load(CfgFlow)

    rebin_factors = cfg_flow['V2ExtractionByYield']['RebinCutsets'][i_cutset] if cfg_flow['V2ExtractionByYield'].get('RebinCutsets') else [1]*(len(pt_bins)-1)

    proj_file = TFile.Open(args.infile, "READ")
    h_resolution = proj_file.Get("hResolution")
    reso = h_resolution.GetBinContent(1)
    outfile = TFile.Open(out_file_name, 'RECREATE')
    for i_pt, (pt_min, pt_max, mass_range, rebin_factor) in enumerate(zip(pt_bins[:-1], pt_bins[1:], \
                                                                          cfg_flow['v2extraction']['MassFitRanges'], \
                                                                          rebin_factors)):
        logger(f"\nProcessing pt bin {i_pt+1}/{len(pt_bins)-1}: {pt_min} - {pt_max} GeV/c", "INFO")
        pt_label = f"pt_{int(pt_min*10)}_{int(pt_max*10)}"
        hist_mass_sp = proj_file.Get(f"{pt_label}/hMassSpData")
        hist_mass_sp.SetDirectory(0)
        sp_edges = get_sp_bin_edges(hist_mass_sp.GetYaxis(), cfg_flow['V2ExtractionByYield'], i_cutset, i_pt)
        if cfg_flow['V2ExtractionByYield'].get('UseTree'):
            data = load_aod_file(f"{cfg_flow['outdir']}/preprocess/{pt_label}/TreesPtCenterSp/AO2D_{pt_label}.root", True)
        else:
            data = hist_mass_sp

        sel_string_cutset = f"fMlScore0 < {cfg_cutset['ScoreBkg']['max'][i_pt]} && " \
                            f"fMlScore0 >= {cfg_cutset['ScoreBkg']['min'][i_pt]} && " \
                            f"fMlScore1 < {cfg_cutset['ScoreFD']['max'][i_pt]} && " \
                            f"fMlScore1 >= {cfg_cutset['ScoreFD']['min'][i_pt]} && " \
                            f"fM >= {mass_range[0]} && fM < {mass_range[1]}"
        hists_stats, weighted_avg, weighted_avg_unc = process_pt_bin(i_cutset, i_pt, cfg_flow, \
                                                                     hists_summary, \
                                                                     rebin_factor, sp_edges, \
                                                                     pt_label, data, reso, \
                                                                     f"{os.path.dirname(out_file_name)}/scan_{cutset_str}/", \
                                                                     sel_string_cutset)

        outfile.mkdir(f"{pt_label}/sp_bins")
        for hist_name, hist in hists_stats.items():
            if 'sp_histos' in hist_name:
                outfile.cd(f"{pt_label}/sp_bins")
                for hist_sp in hist:
                    hist_sp.Write()
            else:
                outfile.cd(pt_label)
                hist.Write(hist_name)

    proj_file.Close()

    outfile.cd()
    for hist_name, hist in hists_summary.items():
        SetObjectStyle(hist, color=kBlack, markerstyle=kFullCircle)
        hist.Write(hist_name)
    outfile.Close()
    logger(f"\n\nProcessed {args.infile} and config file {args.input_config} and saved results to {out_file_name}", "INFO")
