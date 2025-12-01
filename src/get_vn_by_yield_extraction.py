import argparse
import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import yaml
import ROOT
import array
os.environ["CUDA_VISIBLE_DEVICES"] = "" # pylint: disable=wrong-import-position
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
from ROOT import TFile
from matplotlib import gridspec
ROOT.gROOT.SetBatch(True)
import multiprocessing as mp
import sys
sys.path.append("./flareflyfitter/")
from raw_yield_fitter import RawYieldFitter
mp.set_start_method("spawn", force=True)
msg_service = ROOT.RooMsgService.instance()

def process_pt_bin(i_pt, config, histos_summary, rebin_factor, sp_edges, pt_label, hist_mass_sp_int, reso, outdir):
    hist_mass_int = hist_mass_sp_int.ProjectionX()
    fit_cfg = config['v2extraction']

    # Initialize the fitter for sp-integrated yield extraction
    fitter = RawYieldFitter(config['Dmeson'], f"sp_integrated_{pt_label}_fit", config['V2ExtractionByYield']['UseFlareFly'])
    fitter.set_rebin(config.get('rebin', 1))
    fitter.set_fit_range(fit_cfg['MassFitRanges'][i_pt][0], fit_cfg['MassFitRanges'][i_pt][1])
    fitter.set_data_to_fit_hist(hist_mass_int)

    bkg_funcs = fit_cfg['BkgFunc'][i_pt] if isinstance(fit_cfg['BkgFunc'], list) else [fit_cfg['BkgFunc']]
    for bkg_func in bkg_funcs:
        fitter.add_bkg_func(bkg_func, "Comb. bkg")

    sgn_funcs = fit_cfg['SgnFunc'][i_pt] if isinstance(fit_cfg['SgnFunc'], list) else [fit_cfg['SgnFunc']]
    for i_sgn, sgn_func in enumerate(sgn_funcs):
        fitter.add_sgn_func(sgn_func, f"signal_{i_sgn}")

    fitter.setup()
    # fitter.init_corr_bkgs(config, i_pt)
    fitter.fit()
    fig, _ = fitter.plot_fit(False, True, loc=["lower left", "upper left"]) # (log, show_extra_info)
    fig_int_fit_path = f"{outdir}/fit_int_{pt_label}.png"
    os.makedirs(os.path.dirname(fig_int_fit_path), exist_ok=True)
    fig.savefig(fig_int_fit_path, dpi=300, bbox_inches="tight")

    raw_yield, raw_yield_unc = fitter.get_fitter().get_raw_yield()
    mean_int, mean_unc = fitter.get_fitter().get_mass()
    sigma_int, sigma_unc = fitter.get_fitter().get_sigma()

    histos_summary['integrated_yields'].SetBinContent(i_pt + 1, raw_yield)

    # Extract raw yields in each sp bin in parallel
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
    hist_stats['sp_histos'] = [None] * (len(sp_edges['bins']) - 1)

    for isp, (sp_left_bin, sp_right_bin) in enumerate(zip(sp_edges['bins'][:-1], sp_edges['bins'][1:])):
        sp_min = hist_mass_sp_int.GetYaxis().GetBinLowEdge(sp_left_bin)
        sp_max = hist_mass_sp_int.GetYaxis().GetBinUpEdge(sp_right_bin-1)
        sp_center = (sp_min + sp_max) / 2
        sp_label = f"sp_bins_{sp_left_bin}_{sp_right_bin-1}_range_{sp_min:.2f}_{sp_max:.2f}"
        print(f"Processing sp_label {sp_label}")

        hist_mass_sp_int.GetYaxis().SetRange(sp_left_bin, sp_right_bin)
        histo_sp = hist_mass_sp_int.ProjectionX(f"h_mass_{sp_label}")
        histo_sp.SetDirectory(0)

        # For each sp bin, update fitter changing the histogram to fit and name
        fitter.set_data_to_fit_hist(histo_sp)
        fitter.set_name(sp_label)
        fitter.setup()
        if config['V2ExtractionByYield'].get('FixMeanToInt'):  # Fix mean to results of sp-integrated fit
            print("Fixing mean to sp-integrated fit result:", mean_int)
            fitter.fix_sgn_par(len(sgn_funcs)-1, "mu", mean_int)
        if config['V2ExtractionByYield'].get('FixSigmaToInt'): # Fix mean to results of sp-integrated fit
            print("Fixing sigma to sp-integrated fit result:", sigma_int)
            fitter.fix_sgn_par(len(sgn_funcs)-1, "sigma", sigma_int)
        hist_stats['sp_histos'][isp] = histo_sp
        vals_stats['sp_center'][isp] = sp_center
        try:
            fitter.fit()

            fig, _ = fitter.plot_fit(False, True, loc=["lower left", "upper left"]) # (log, show_extra_info)
            fig_sp_pt_path = f"{outdir}/{pt_label}/{sp_label}.png"
            os.makedirs(os.path.dirname(fig_sp_pt_path), exist_ok=True)
            fig.savefig(fig_sp_pt_path, dpi=300, bbox_inches="tight")

            vals_stats['sp_ry'][isp], vals_stats['sp_ry_unc'][isp] = fitter.get_fitter().get_raw_yield()
            vals_stats['means'][isp], vals_stats['means_unc'][isp] = fitter.get_fitter().get_mass()
            vals_stats['sigmas'][isp], vals_stats['sigmas_unc'][isp] = fitter.get_fitter().get_sigma()
        except Exception as e:
            print(f"Error fitting sp {sp_min:.2f} - {sp_max:.2f}: {e}")
            vals_stats['sp_ry'][isp], vals_stats['sp_ry_unc'][isp] = 0, 0
            vals_stats['means'][isp], vals_stats['means_unc'][isp] = 0, 0
            vals_stats['sigmas'][isp], vals_stats['sigmas_unc'][isp] = 0, 0

    # Compute weighted average and fill histograms
    sum_yields = np.sum(vals_stats['sp_ry'])
    histos_summary["summed_sp_yields"].SetBinContent(i_pt + 1, sum_yields)
    sum_weighted = np.sum(np.array(vals_stats['sp_ry']) * (np.array(vals_stats['sp_center']) / reso))
    weighted_avg = sum_weighted / sum_yields
    weighted_avg_unc = 0
    for i_sp, (i_v2, i_ry_unc) in enumerate(zip(vals_stats['sp_center'], vals_stats['sp_ry_unc'])):
        d_yields_dRyi = ((i_v2/reso) / sum_yields)
        d_v2_dRyi = (sum_weighted) / (sum_yields**2)
        weighted_avg_unc = weighted_avg_unc + ((d_yields_dRyi - d_v2_dRyi)**2)*i_ry_unc**2

        vals_stats['weight_av_unc_first_term'][i_sp] = ((i_v2/reso) / sum_yields)**2 * i_ry_unc**2
        hist_stats['weight_av_unc_first_term'].SetBinContent(i_sp + 1, vals_stats['weight_av_unc_first_term'][i_sp])
        vals_stats['weight_av_unc_second_term'][i_sp] = ((sum_weighted) / (sum_yields**2))**2 * i_ry_unc**2
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
    return hist_stats, weighted_avg, weighted_avg_unc

def get_sp_bin_edges(axis_sp, cfg, i_cut, i_pt):

    sp_abs_val_max = cfg['sp_ranges'][i_cut]
    if sp_abs_val_max <= 0:
        raise ValueError(f"Invalid sp_abs_val_max: {sp_abs_val_max}. It must be positive.")
        sys.exit(1)

    sp_windows_nbins = cfg['sp_windows_nbins'][i_cut][i_pt]
    sp_bins = axis_sp.GetNbins()

    # Find first bin
    for i in range(1, sp_bins + 1):
        if axis_sp.GetBinLowEdge(i)-0.0001 > -sp_abs_val_max:
            print(f"First bin found at i={i}, low edge={axis_sp.GetBinLowEdge(i)}")
            first_bin = i-1
            break

    # Find last bin
    for i in range(1, sp_bins + 1):
        if axis_sp.GetBinLowEdge(i)-0.0001 > sp_abs_val_max:
            print(f"Last bin found at i={i}, low edge={axis_sp.GetBinLowEdge(i)}")
            last_bin = i-1
            break
    
    edges = {}
    edges['bins'] = list(range(first_bin, last_bin, sp_windows_nbins))
    edges['sp_mins'] = [axis_sp.GetBinLowEdge(b) for b in edges['bins']]
    edges['sp_maxs'] = [axis_sp.GetBinUpEdge(b + sp_windows_nbins - 1) for b in edges['bins']]
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
    out_file_name = os.path.join(args.infile.replace('proj', 'raw_yield'))
    os.makedirs(os.path.dirname(out_file_name), exist_ok=True)

    pt_bins = cfg_flow['ptbins']
    histos_summary_yields = {
        'integrated_yields': ROOT.TH1F("hist_integrated_yields", "hist_integrated_yields", len(pt_bins) - 1, np.array(pt_bins)),
        'summed_sp_yields': ROOT.TH1F("hist_summed_sp_yields", "hist_summed_sp_yields", len(pt_bins) - 1, np.array(pt_bins))
    }

    pt_bins = cfg_flow['ptbins']
    
    # Retrieve number of cutset
    _, cutset_str = os.path.splitext(os.path.basename(args.infile))[0].split('_')
    i_cutset = int(cutset_str)
    rebin_factors = cfg_flow['V2ExtractionByYield']['rebin_cutsets'][i_cutset]

    proj_file = TFile.Open(args.infile, "READ")
    h_resolution = proj_file.Get("hResolution")
    reso = h_resolution.GetBinContent(1)
    outfile = TFile.Open(out_file_name, 'RECREATE')
    hist_v2_values = ROOT.TH1F("hist_v2_values", "hist_v2_values", len(pt_bins) - 1, np.array(pt_bins))
    for i_pt, (pt_min, pt_max, mass_range, rebin_factor) in enumerate(zip(pt_bins[:-1], pt_bins[1:], cfg_flow['MassFitRanges'], rebin_factors)):
        if i_pt != 3:
            continue
        print(f"\nProcessing pt bin {i_pt+1}/{len(pt_bins)-1}: {pt_min} - {pt_max} GeV/c")
        pt_label = f"pt_{int(pt_min*10)}_{int(pt_max*10)}"
        hist_mass_sp_int = proj_file.Get(f"{pt_label}/hMassSpData")
        hist_mass_sp_int.SetDirectory(0)
        sp_edges = get_sp_bin_edges(hist_mass_sp_int.GetYaxis(), cfg_flow['V2ExtractionByYield'], i_cutset, i_pt)
        hists_stats, weighted_avg, weighted_avg_unc = process_pt_bin(i_pt, cfg_flow, \
                                                                     histos_summary_yields, \
                                                                     rebin_factor, sp_edges, \
                                                                     pt_label, hist_mass_sp_int, reso, \
                                                                     f"{os.path.dirname(out_file_name)}/scan_{cutset_str}/")
        print(f"Weighted average v2: {weighted_avg} ± {weighted_avg_unc}")
        hist_v2_values.SetBinContent(i_pt + 1, weighted_avg)
        hist_v2_values.SetBinError(i_pt + 1, weighted_avg_unc)

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
    for hist_name, hist in histos_summary_yields.items():
        hist.Write(hist_name)
    hist_v2_values.Write("hist_v2_values")
    outfile.Close()
    print(f"\n\nProcessed {args.infile} and config file {args.input_config} and saved results to {out_file_name}")
