import argparse
import os
import re
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import yaml
import ROOT
import array
os.environ["CUDA_VISIBLE_DEVICES"] = "" # pylint: disable=wrong-import-position
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
from ROOT import TFile, TH1F, TGraphAsymmErrors, kBlack, kFullCircle
script_dir = os.path.dirname(os.path.realpath(__file__))
os.sys.path.append(os.path.join(script_dir, '..', 'utils'))
from utils import logger, make_dir_root_file
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

def sp_scan_pt_bin(fitter, data_cutset, sgn_func_label, pt_label, reso, sp_intervals, outdir):

    vals_stats = {}

    # Extract raw yields in each sp bin
    stats = [
        'Means', 'MeansUnc',
        'Sigmas', 'SigmasUnc',
        'SpRy', 'SpRyUnc',
        'Chi2', 'Chi2OverNdf',
        'VnSimFitUnc', 'SpErrorContrib',
        'SpFirstTermErrorContrib', 'SpSecondTermErrorContrib',
    ]
    for stat in stats:
        vals_stats[stat] = [0] * (len(sp_intervals) - 1)
    vals_stats['SpCenter'] = [0] * (len(sp_intervals) - 1)

    for isp, (sp_min, sp_max) in enumerate(zip(sp_intervals[:-1], sp_intervals[1:])):
        sp_center = (sp_min + sp_max) / 2
        sp_label = f"sp_range_{sp_min:.2f}_{sp_max:.2f}"
        logger(f"Processing sp_label {sp_label}", "INFO")

        # For each sp bin, update fitter changing the histogram to fit, name and rebin
        data_cutset_sp = data_cutset.query(f"fScalarProd >= {sp_min} and fScalarProd < {sp_max}")
        if len(data_cutset_sp) == 0:
            logger(f"No entries found in sp bin {sp_label}. Setting all values to 0.", "WARNING")
            vals_stats['SpRy'][isp], vals_stats['SpRyUnc'][isp] = 0, 0
            vals_stats['Means'][isp], vals_stats['MeansUnc'][isp] = 0, 0
            vals_stats['Sigmas'][isp], vals_stats['SigmasUnc'][isp] = 0, 0
            vals_stats['Chi2'][isp], vals_stats['Chi2OverNdf'][isp] = -1, -1
            continue
        fitter.set_data_to_fit_df(data_cutset_sp, 'fM')

        fitter.set_name(sp_label)
        fitter.setup()
        vals_stats['SpCenter'][isp] = sp_center
        try:
            status, converged = fitter.fit(verbose=False)

            fig_sp_pt_path = f"{outdir}/{pt_label}/{sp_label}.pdf"
            fitter.plot_fit(False, True, loc=["lower left", "upper left"], # (log, show_extra_info)
                            path=fig_sp_pt_path)
            fit_info, sgn_pars, sgn_pars_uncs, _, _ = fitter.get_fit_info()
            vals_stats['SpRy'][isp], vals_stats['SpRyUnc'][isp] = fit_info[sgn_func_label]["ry"], fit_info[sgn_func_label]["ry_unc"]
            vals_stats['Means'][isp], vals_stats['MeansUnc'][isp] = sgn_pars[f"mu_{sgn_func_label}"], sgn_pars_uncs[f"mu_{sgn_func_label}"]
            vals_stats['Sigmas'][isp], vals_stats['SigmasUnc'][isp] = sgn_pars[f"sigma_{sgn_func_label}"], sgn_pars_uncs[f"sigma_{sgn_func_label}"]
            vals_stats['Chi2'][isp], vals_stats['Chi2OverNdf'][isp] = fit_info["chi2"], fit_info["chi2_over_ndf"]
        except Exception as e:
            logger(f"Error fitting sp {sp_min:.2f} - {sp_max:.2f}: {e}", "WARNING")
            logger(f"Setting all values to 0 for this sp bin.", "WARNING")
            vals_stats['SpRy'][isp], vals_stats['SpRyUnc'][isp] = 0, 0
            logger("sp set", "WARNING")
            vals_stats['Means'][isp], vals_stats['MeansUnc'][isp] = 0, 0
            logger("means set", "WARNING")
            vals_stats['Sigmas'][isp], vals_stats['SigmasUnc'][isp] = 0, 0
            logger("sigmas set", "WARNING")
            vals_stats['Chi2'][isp], vals_stats['Chi2OverNdf'][isp] = -1, -1
            logger("chi2 set", "WARNING")

    # Compute weighted averages and fill histograms
    vals_stats['SummedSpYields'] = np.sum(vals_stats['SpRy'])
    vals_stats['SummedSpYieldsUnc'] = np.sqrt(np.sum(np.array(vals_stats['SpRyUnc'])**2))
    vals_stats['WeightedSum'] = np.sum(np.array(vals_stats['SpRy']) * (np.array(vals_stats['SpCenter']) / reso))
    vals_stats['VnSimFit'] = vals_stats['WeightedSum'] / vals_stats['SummedSpYields']
    weighted_avg_unc = 0
    for i_sp, (i_v2, i_ry, i_ry_unc) in enumerate(zip(vals_stats['SpCenter'], vals_stats['SpRy'], vals_stats['SpRyUnc'])):
        d_yields_dRyi = (i_v2/reso) / vals_stats['SummedSpYields']
        d_v2_dRyi = (vals_stats['WeightedSum']) / (vals_stats['SummedSpYields']**2)
        weighted_avg_unc = weighted_avg_unc + ((d_yields_dRyi - d_v2_dRyi)**2)*i_ry_unc**2
        vals_stats['SpFirstTermErrorContrib'][i_sp] = (d_yields_dRyi)**2 * i_ry_unc**2
        vals_stats['SpSecondTermErrorContrib'][i_sp] = (d_v2_dRyi)**2 * i_ry_unc**2
        vals_stats['SpErrorContrib'][i_sp] = np.sqrt(abs(vals_stats['SpFirstTermErrorContrib'][i_sp] - \
                                                         vals_stats['SpSecondTermErrorContrib'][i_sp]))
    vals_stats['VnSimFitUnc'] = np.sqrt(weighted_avg_unc)

    return vals_stats

def run_pt_bin_worker(cutset_cfgs, i_pt, pt_min, pt_max, cfg_flow, reso):
    pt_label = f"pt_{int(pt_min*10)}_{int(pt_max*10)}"
    data = load_aod_file(f"{cfg_flow['outdir']}/preprocess/{pt_label}/TreesPtCenterSp/AO2D_{pt_label}.root", True)

    # Initialize the fitter for sp-integrated yield extraction
    fit_cfg = cfg_flow['v2extraction']
    part = cfg_flow['Dmeson']
    fitter = RawYieldFitter(part, pt_min, pt_max, f"sp_integrated_{pt_label}_fit", fit_cfg['Minimizer'])
    fitter.set_fit_range(fit_cfg['MassFitRanges'][i_pt][0], fit_cfg['MassFitRanges'][i_pt][1])

    # Add model components
    fitter.add_bkg_func(fit_cfg['BkgFunc'][i_pt] if isinstance(fit_cfg['BkgFunc'], list) else fit_cfg['BkgFunc'], "Comb. bkg")
    sgn_funcs = {} # More info for signal functions, a dictionary is better
    sgn_funcs[fit_cfg['SgnFuncLabel']] = {
        'func': fit_cfg['SgnFunc'][i_pt] if isinstance(fit_cfg['SgnFunc'], list) else fit_cfg['SgnFunc'],
        'part': part
    }
    if fit_cfg.get('InclSecPeak'):
        include_sec_peak = fit_cfg['InclSecPeak'][i_pt] if isinstance(fit_cfg['InclSecPeak'], list) else fit_cfg['InclSecPeak']
        if include_sec_peak:
            sgn_funcs[fit_cfg['SgnFuncSecPeakLabel']] = {
                'func': fit_cfg['SgnFuncSecPeak'][i_pt] if isinstance(fit_cfg['SgnFuncSecPeak'], list) else fit_cfg['SgnFuncSecPeak'],
                'part': 'Dplus' if part == 'Ds' else 'Dstar',
            }
    for i_sgn, (label, sgn_func) in enumerate(sgn_funcs.items()):
        fitter.add_sgn_func(sgn_func['func'], label, sgn_func['part'])

    stats = {}
    for i_cutset, cutset_file in enumerate(cutset_cfgs):
        with open(cutset_file, 'r') as CfgCutset:
            cfg_cutset = yaml.safe_load(CfgCutset)
        outdir = os.path.dirname(cutset_file).replace('cutset', "raw_yield")

        sel_string_cutset = (
            f"fMlScore0 < {cfg_cutset['ScoreBkg']['max'][i_pt]} && "
            f"fMlScore0 >= {cfg_cutset['ScoreBkg']['min'][i_pt]} && "
            f"fMlScore1 < {cfg_cutset['ScoreFD']['max'][i_pt]} && "
            f"fMlScore1 >= {cfg_cutset['ScoreFD']['min'][i_pt]} && "
            f"fM >= {cfg_flow['v2extraction']['MassFitRanges'][i_pt][0]} && "
            f"fM < {cfg_flow['v2extraction']['MassFitRanges'][i_pt][1]}"
        )

        data_cutset = data.query(sel_string_cutset.replace(" && ", " and "))
        fitter.set_data_to_fit_df(data_cutset, 'fM')

        # Add correlated background if specified
        if cfg_flow.get('corr_bkgs'):
            fitter.add_corr_bkgs(cfg_flow['corr_bkgs'], sel_string, pt_min, pt_max)

        fitter.setup()
        if fit_cfg.get('InitPars'):
            fitter.set_fit_pars(fit_cfg['InitPars'], pt_min, pt_max)

        # Prefit the MC prompt enhanced cut to fix the tails, binned fit
        if fit_cfg.get('FixSgnFromMC'):
            fitter.set_fix_sgn_to_mc_prefit(True)
            fitter.prefit_mc(f"{cfg_flow['outdir']}/corrbkgs/templs_{pt_label}.root")
            fitter.plot_mc_prefit(False, True, loc=["lower left", "upper left"], path=outdir)
            fitter.plot_raw_residuals_mc_prefit(path=f"{outdir}/fM_mc_prefit_residuals_{pt_label}.pdf")

        if fit_cfg.get('FixParsToIntFit'):  # Fix mean to results of sp-integrated fit
            logger(f"Will fix all parameters of signal functions to sp-integrated fit result", "WARNING")
            fitter.fix_sgn_pars_to_first_fit()
        status, converged = fitter.fit(verbose=True)
        fig_int_fit_path = f"{outdir}/fit_int_{pt_label}.pdf"
        fitter.plot_fit(False, True, loc=["lower left", "upper left"],  # (log, show_extra_info)
                        path=fig_int_fit_path)

        fit_info_pt_int, sgn_pars_pt_int, sgn_pars_uncs_pt_int, _, _ = fitter.get_fit_info()
        label = list(sgn_funcs.keys())[0]       # Take only the first signal function, which is the peak of interest
        sgn_func_label = fit_cfg['SgnFuncLabel']

        step = fit_cfg['SpWindowWidth'][i_cutset][i_pt]
        sp_abs_val_max = fit_cfg['SpRanges'][i_pt]
        sp_intervals = np.arange(-sp_abs_val_max, sp_abs_val_max + 0.5 * step, step).tolist()
        stats[i_cutset] = {}
        stats[i_cutset] = sp_scan_pt_bin(fitter, data_cutset, sgn_func_label, pt_label, reso, sp_intervals,
                                         f"{outdir}/scan_{i_cutset:02d}")
        stats[i_cutset]['SpIntervals'] = sp_intervals
        stats[i_cutset]['RawYieldsSimFit'] = fit_info_pt_int[sgn_func_label]["ry"]
        stats[i_cutset]['RawYieldsSimFitUnc'] = fit_info_pt_int[sgn_func_label]["ry_unc"]
        stats[i_cutset]['MeanSimFit'] = sgn_pars_pt_int[f"mu_{sgn_func_label}"]
        stats[i_cutset]['MeanSimFitUnc'] = sgn_pars_uncs_pt_int[f"mu_{sgn_func_label}"]
        stats[i_cutset]['SigmaSimFit'] = sgn_pars_pt_int[f"sigma_{sgn_func_label}"]
        stats[i_cutset]['SigmaSimFitUnc'] = sgn_pars_uncs_pt_int[f"sigma_{sgn_func_label}"]
        # Clear data_cutset to save memory
        del data_cutset

    # IMPORTANT: return ONLY python objects
    return pt_label, i_pt, stats

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('input_config', metavar='text', default='config.yml')
    parser.add_argument('infile', metavar='text', default='proj_XX.root')
    parser.add_argument('--batch', '-b', help='suppress video output', action='store_true')
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    sys.exit(0)
    with open(args.input_config, 'r') as CfgFlow:
        cfg_flow = yaml.safe_load(CfgFlow)

    # Retrieve cutsets configs
    try:
        cutsets_dir = os.path.join(cfg_flow['outdir'], f"cutvar_{cfg_flow['suffix']}_combined/cutsets")
        cutset_files = [os.path.join(cutsets_dir, f) for f in os.listdir(cutsets_dir) if f.endswith('.yml')]
        out_dir_type = "combined"
    except Exception as e:
        logger(f"Could not find combined cutsets, trying correlated cutsets ... ", level="WARNING")
        cutsets_dir = os.path.join(cfg_flow['outdir'], f"cutvar_{cfg_flow['suffix']}_correlated/cutsets")
        cutset_files = [os.path.join(cutsets_dir, f) for f in os.listdir(cutsets_dir) if f.endswith('.yml')]
        out_dir_type = "correlated"
    cutset_files.sort(key=lambda x: int(re.search(r'(\d+)', os.path.basename(x)).group(1)))

    logger(f"Using cutsets {cutset_files}\n", "INFO")

    # Read the Resolution and MassSp histograms for all pt bins
    proj_file = TFile.Open(args.infile, "READ")
    h_resolution = proj_file.Get("hResolution")
    reso = h_resolution.GetBinContent(1)
    proj_file.Close()

    tasks = []
    vals_stats = {}
    pt_bins = cfg_flow['ptbins']
    for i_pt, (pt_min, pt_max) in enumerate(zip(pt_bins[:-1], pt_bins[1:])):
        tasks.append((cutset_files, i_pt, pt_min, pt_max, cfg_flow, reso))
    with ProcessPoolExecutor(max_workers=cfg_flow['v2extraction'].get('Workers', 1)) as executor:
        results = [executor.submit(run_pt_bin_worker, *task) for task in tasks]
        for task in as_completed(results):
            pt_label, i_pt, stats = task.result()
            vals_stats[pt_label] = (i_pt, stats)

    # Sort vals_stats by i_pt
    vals_stats = dict(sorted(vals_stats.items(), key=lambda item: item[1][0]))

    out_files, summaries = {}, {}
    for i_cutset, cutset in enumerate(cutset_files):
        os.makedirs(os.path.dirname(cutset).replace('cutset', 'raw_yield'), exist_ok=True)
        out_files[i_cutset] = cutset.replace('cutset', 'raw_yield').replace('.yml', '.root')
        summaries[i_cutset] = {
            'hRawYieldsSimFit': TH1F("hRawYieldsSimFit", "hRawYieldsSimFit", len(pt_bins) - 1, np.array(pt_bins)),
            'hSummedSpYields': TH1F("hSummedSpYields", "hSummedSpYields", len(pt_bins) - 1, np.array(pt_bins)),
            'hMeanSimFit': TH1F("hMeanSimFit", "hMeanSimFit", len(pt_bins) - 1, np.array(pt_bins)),
            'hSigmaSimFit': TH1F("hSigmaSimFit", "hSigmaSimFit", len(pt_bins) - 1, np.array(pt_bins)),
            'hVnSimFit': TH1F("hVnSimFit", "hVnSimFit", len(pt_bins) - 1, np.array(pt_bins)),
            'hVnSimFitUnc': TH1F("hVnSimFitUnc", "hVnSimFitUnc", len(pt_bins) - 1, np.array(pt_bins)),
            'hWeightedSum': TH1F("hWeightedSum", "hWeightedSum", len(pt_bins) - 1, np.array(pt_bins)),
            'gVnSimFit': TGraphAsymmErrors(1),
            'gVnUnc': TGraphAsymmErrors(1),
        }
        summaries[i_cutset]['gVnSimFit'].SetName("gVnSimFit")
        summaries[i_cutset]['gVnUnc'].SetName("gVnUnc")

    for pt_label, (i_pt, pt_stats) in vals_stats.items():
        pt_min = pt_bins[i_pt]
        pt_max = pt_bins[i_pt + 1]

        for cutset, cutset_vals in pt_stats.items():
            for var, vals in cutset_vals.items():
                if isinstance(vals, list):
                    summaries[cutset][f"{pt_label}/{var}"] = TH1F(f"hist_{var}", f"hist_{var}", len(cutset_vals['SpIntervals']) - 1, array.array('f', cutset_vals['SpIntervals']))
                    summaries[cutset][f"{pt_label}/{var}"].SetDirectory(0)
                    for i_val, val in enumerate(vals):
                        summaries[cutset][f"{pt_label}/{var}"].SetBinContent(i_val + 1, val)
                        if 'Unc' in var:
                            summaries[cutset][f"{pt_label}/{var.replace('Unc', '')}"].SetBinError(i_val + 1, val)
                else:
                    if "Unc" in var:
                        continue
                    summaries[cutset][f"h{var}"].SetBinContent(i_pt + 1, vals)
                    try:
                        summaries[cutset][f"h{var}"].SetBinError(i_pt + 1, cutset_vals[f"{var}Unc"])
                    except KeyError:
                        pass
                    if "VnSimFit" in var:
                        summaries[cutset]['gVnSimFit'].SetPoint(i_pt, (pt_min+pt_max)/2, vals)
                        summaries[cutset]['gVnSimFit'].SetPointError(i_pt, (pt_max-pt_min)/2, (pt_max-pt_min)/2, cutset_vals[f"{var}Unc"], cutset_vals[f"{var}Unc"])
                        summaries[cutset]['gVnUnc'].SetPoint(i_pt, (pt_min+pt_max)/2, cutset_vals[f"{var}Unc"])
                        summaries[cutset]['gVnUnc'].SetPointError(i_pt, (pt_max-pt_min)/2, (pt_max-pt_min)/2, 1.e-20, 1.e-20)

    for i_cutset, cutset in enumerate(cutset_files):
        out_file_name = cutset.replace('cutsets', 'raw_yields').replace('cutset', 'raw_yields').replace('.yml', '.root')
        outfile = TFile.Open(out_file_name, 'RECREATE')
        for pt_label in vals_stats.keys():
            make_dir_root_file(pt_label, outfile, verbose=False)
        for hist_name, hist in summaries[i_cutset].items():
            SetObjectStyle(hist, color=kBlack, markerstyle=kFullCircle)
            if "/" in hist_name:
                pt_dir, hist_name = hist_name.split("/")
                outfile.cd(pt_dir)
                hist.SetName(hist_name)
                hist.Write(hist_name)
            else:
                outfile.cd()
                hist.Write(hist_name)
        outfile.Close()

    logger(f"Processed {args.infile} and config file {args.input_config} and saved results to {out_file_name}\n\n", "INFO")
