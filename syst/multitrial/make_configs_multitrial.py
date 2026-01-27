import yaml

import copy
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import ROOT
ROOT.DisableImplicitMT()

import numpy as np
import argparse
import itertools
import sys
sys.path.append(f"{os.path.dirname(os.path.abspath(__file__))}/../../utils")
from utils import logger

def get_reference_config_pt_bin(config_flow, iPtBin):
    # Deepcopy to ensure that we are working with a copy and not modifying the original config_flow
    pt_bin_config = copy.deepcopy(config_flow)

    # Separate bin-specific keys (those with lists as values) from common settings
    for key, value in pt_bin_config.items():
        if isinstance(value, list):
            if key == 'ptbins':
                pt_bin_config[key] = [value[iPtBin], value[iPtBin+1]]
            else:
                pt_bin_config[key] = value[iPtBin]

    # # Modify cut_variation values but leave the original config_flow intact
    pt_bin_config["cut_variation"]["corr_bdt_cut"]["bkg_max"] = [pt_bin_config["cut_variation"]["corr_bdt_cut"]["bkg_max"][iPtBin]]
    pt_bin_config["cut_variation"]["corr_bdt_cut"]["sig"]["max"] = [pt_bin_config["cut_variation"]["corr_bdt_cut"]["sig"]["max"][iPtBin]]
    pt_bin_config["cut_variation"]["corr_bdt_cut"]["sig"]["min"] = [pt_bin_config["cut_variation"]["corr_bdt_cut"]["sig"]["min"][iPtBin]]
    pt_bin_config["cut_variation"]["corr_bdt_cut"]["sig"]["step"] = [pt_bin_config["cut_variation"]["corr_bdt_cut"]["sig"]["step"][iPtBin]]
    pt_bin_config["cut_variation"]["uncorr_bdt_cut"]["bkg_max"] = [pt_bin_config["cut_variation"]["uncorr_bdt_cut"]["bkg_max"][iPtBin]]
    pt_bin_config["cut_variation"]["uncorr_bdt_cut"]["sig"] = [pt_bin_config["cut_variation"]["uncorr_bdt_cut"]["sig"][iPtBin]]


    pt_bin_config['v2extraction']['MassFitRanges'] = [pt_bin_config['v2extraction']['MassFitRanges'][iPtBin]]
    pt_bin_config['projections']['inv_mass_bins'] = [pt_bin_config['projections']['inv_mass_bins'][iPtBin]]
    keys = ["Sigma", "Rebin", "BkgFunc", "SgnFunc", "BkgFuncVn"]
    v2extr = pt_bin_config["v2extraction"]
    for key in keys:
        v2extr[key] = [v2extr[key][iPtBin]] if isinstance(v2extr[key], list) else [v2extr[key]]
    if "SpRanges" in pt_bin_config["v2extraction"]:
        pt_bin_config["v2extraction"]["SpRanges"] = pt_bin_config["v2extraction"]["SpRanges"][iPtBin]
    if "NSigma4SB" in pt_bin_config["v2extraction"]:
        pt_bin_config["v2extraction"]["NSigma4SB"] = pt_bin_config["v2extraction"]["NSigma4SB"][iPtBin]

    return pt_bin_config

def produce_trials_fit_configs(cfg_flow, ptmin, ptmax, ref_config_ptbin, setting, multitrial_pt_dir):

    varied_configs = {}
    for var in setting["multitrial"]:
        print(f"Processing variation for {var}")
        if var in ref_config_ptbin:
            varied_configs[var] = list(set(ref_config_ptbin[var] + setting["multitrial"][var]))
        if var in ref_config_ptbin["projections"]:
            varied_configs[var] = list(set(ref_config_ptbin["projections"][var] + setting["multitrial"][var]))
        if var in ref_config_ptbin["v2extraction"]:
            varied_configs[var] = list(set(ref_config_ptbin["v2extraction"][var] + setting["multitrial"][var]))
        else:
            varied_configs[var] = setting["multitrial"][var]
        varied_configs[var] = np.unique(varied_configs[var]).tolist()
        print(f"  Values: {varied_configs[var]}")

    # Perform itertools.product to get all combinations
    keys = list(varied_configs.keys())
    values = list(varied_configs.values())
    combinations = list(itertools.product(*values))
    config_variants = [dict(zip(keys, combination)) for combination in combinations] 

    # Print results
    for idx, variant in enumerate(config_variants):
        cfg_variant = copy.deepcopy(ref_config_ptbin)
        cfg_variant["iTrial"] = idx

        cfg_variant["outdir"] = f"{multitrial_pt_dir}/trials/{idx}/"
        if setting['MaxChi2'] > 10:
            logger(f"MaxChi2 is set to {setting['MaxChi2']} for pt bin {ptmin}-{ptmax}, please consider a tight selection!", "WARNING")
        cfg_variant["MaxChi2"] = setting['MaxChi2']
        if setting['MinSignificance'] < 5:
            logger(f"MinSignificance is set to {setting['MinSignificance']}, please consider a higher lower limit!", "WARNING")
        cfg_variant["MinSignificance"] = setting['MinSignificance']
        if setting['MaxSignificance'] > 1000:
            logger(f"MaxSignificance is set to {setting['MaxSignificance']}, please consider a lower upper limit!", "WARNING")
        cfg_variant["MaxSignificance"] = setting['MaxSignificance']
        for varied_var in variant:
            if ref_config_ptbin.get(varied_var):
                if isinstance(ref_config_ptbin[varied_var], list):
                    cfg_variant[varied_var] = [variant[varied_var]]
                else:
                    cfg_variant[varied_var] = variant[varied_var]

        if variant.get("MassMin") and variant.get("MassMax"):
            cfg_variant["v2extraction"]["MassFitRanges"] = [[variant["MassMin"], variant["MassMax"]]]
            if variant.get("SpWindowWidth"):
                cfg_variant["v2extraction"]["SpWindowWidth"] = [variant["SpWindowWidth"]]
            if variant.get("InvMassBinSteps"): # Not needed for CMS method, only for simfit
                cfg_variant["projections"]["inv_mass_bins"] = [np.arange(variant['MassMin'],
                                                                         variant['MassMax'] + variant['InvMassBinSteps'],
                                                                         variant['InvMassBinSteps']).tolist()]
            # else the region used for the prefit of the bkg function would be empty, leading to a crash
            if cfg_flow['Dmeson'] == 'Dplus' and (variant["MassMin"] > 1.75 or variant["MassMax"] < 1.95):
                cfg_variant["NSigma4SB"] = [2]
                cfg_variant["v2extraction"]["Sigma"] = [0.01]
        if variant.get("BkgFunc"):
            cfg_variant["v2extraction"]["BkgFunc"] = variant["BkgFunc"]
        if variant.get("SgnFunc"):
            cfg_variant["v2extraction"]["SgnFunc"] = variant["SgnFunc"]
        if variant.get("BkgFuncVn"):
            cfg_variant["v2extraction"]["BkgFuncVn"] = variant["BkgFuncVn"]
        if variant.get("SpRanges"):
            cfg_variant["v2extraction"]["SpRanges"] = [variant["SpRanges"]]
        if variant.get("Rebin"):
            cfg_variant["v2extraction"]["Rebin"] = variant["Rebin"]
        if variant.get("NSigma4SB"):
            cfg_variant["v2extraction"]["NSigma4SB"] = [variant["NSigma4SB"]]

        os.makedirs(os.path.join(f"{multitrial_pt_dir}/trials/{idx}"), exist_ok=True)
        output_file = os.path.join(f"{multitrial_pt_dir}/trials/{idx}", f"config_trial_{idx}.yml")

        # No preprocess and operations are controlled in the bash script
        cfg_variant.pop('operations', None)
        cfg_variant.pop('preprocess', None)

        # Cut variation is taken from reference results
        cfg_variant.pop('minimisation', None)

        with open(output_file, 'w') as out_file:
            yaml.dump(cfg_variant, out_file, default_flow_style=False, sort_keys=False)

def produce_trials_bdt_configs(cfg_flow, ref_config_ptbin, setting, multitrial_pt_dir):

    varied_configs = {}
    for i_cutset, cutset_cfg in enumerate(setting["multitrial"]["cutsets"]):
        if cutset_cfg.get("BkgScoreStep"):
            bdt_bkg_cuts = np.arange(cutset_cfg["BkgScoreRange"][0], cutset_cfg["BkgScoreRange"][1] + \
                                    cutset_cfg["BkgScoreStep"], cutset_cfg["BkgScoreStep"]).tolist()
        else:
            bdt_bkg_cuts = cutset_cfg["BkgScoreVals"]
        for i_bkg_cut, bkg_cut in enumerate(bdt_bkg_cuts):
            bkg_cut = round(bkg_cut, 3)  # To avoid floating point precision issues
            # Have equal strings, so pad 0.0x to 0.0x0 to match, e.g., 0.05 and 0.052 number of digits
            bkg_cut_str = f"{bkg_cut:.3f}"
            out_dir = f"{multitrial_pt_dir}/trials_cutset_{i_cutset}/bkg_{bkg_cut_str}/"
            os.makedirs(out_dir, exist_ok=True)

            cfg_variant = copy.deepcopy(ref_config_ptbin)
            cfg_variant["iTrial"] = i_bkg_cut
            cfg_variant["cut_variation"]["uncorr_bdt_cut"]["bkg_max"] = [[bkg_cut]]
            cfg_variant["cut_variation"]["uncorr_bdt_cut"]["sig"] = [[ref_config_ptbin["cut_variation"]["uncorr_bdt_cut"]["sig"][0][i_cutset],
                                                                      ref_config_ptbin["cut_variation"]["uncorr_bdt_cut"]["sig"][0][i_cutset+1]]]
            cfg_variant["outdir"] = out_dir

            # Setup operations
            cfg_variant["operations"]['preprocess'] = False
            cfg_variant["operations"]['make_yaml'] = True
            cfg_variant["operations"]['proj_data'] = True
            cfg_variant["operations"]['proj_mc'] = True
            cfg_variant["operations"]['efficiencies'] = True
            cfg_variant["operations"]['get_vn_vs_mass'] = True
            cfg_variant["operations"]['do_cut_variation'] = False
            cfg_variant["operations"]['data_driven_fraction'] = False
            cfg_variant["operations"]['get_v2_vs_frac'] = False

            # Pick up pre-processed files
            cfg_variant["outdirPrep"] = ref_config_ptbin["outdir"]
            with open(f"{out_dir}/config_trial_{i_bkg_cut}.yml", 'w') as out_file:
                yaml.dump(cfg_variant, out_file, default_flow_style=False, sort_keys=False)

def produce_trials_cfgs(config_flow, config_mod, output_dir, multitrial_type='fit'):

    with open(config_flow, 'r') as CfgFlow:
        cfg_flow = yaml.safe_load(CfgFlow)
    with open(config_mod, 'r') as CfgMod:
        cfg_mod = yaml.safe_load(CfgMod)

    multitrial_dir = f"{output_dir}/syst/multitrial/{multitrial_type}/"
    os.makedirs(multitrial_dir, exist_ok=True)
    os.makedirs(f"{multitrial_dir}/config_history/", exist_ok=True)

    with open(os.path.join(f"{multitrial_dir}/config_history/", f"config_reference.yml"), 'w') as out_file:
        yaml.dump(cfg_flow, out_file, default_flow_style=False)
    with open(os.path.join(f"{multitrial_dir}/config_history/", f"config_modifications.yml"), 'w') as out_file:
        yaml.dump(cfg_mod, out_file, default_flow_style=False)

    for _, setting in enumerate(cfg_mod['ptbins']):
        for ptmin, ptmax in setting['ranges']:
            logger(f"Generating multitrial configs for pt bin {ptmin}-{ptmax} GeV/c", "INFO")

            pt_bin_index = cfg_flow['ptbins'].index(ptmin)
            ref_config_ptbin = get_reference_config_pt_bin(cfg_flow, pt_bin_index)
            multitrial_pt_dir = f"{multitrial_dir}/pt_{int(ptmin*10)}_{int(ptmax*10)}/"
            
            os.makedirs(multitrial_pt_dir, exist_ok=True)
            output_file = os.path.join(multitrial_pt_dir, f"config_reference.yml")
            with open(output_file, 'w') as out_file:
                yaml.dump(ref_config_ptbin, out_file, default_flow_style=False)
                
            if multitrial_type == 'fit':
                produce_trials_fit_configs(cfg_flow, ptmin, ptmax, ref_config_ptbin, setting, multitrial_pt_dir)
            elif multitrial_type == 'bdt':
                produce_trials_bdt_configs(cfg_flow, ref_config_ptbin, setting, multitrial_pt_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('input_config', metavar='text', default='config_Ds_Fit.yml')
    parser.add_argument('--fit_multitrial', "-fm", action='store_true', default=False)
    parser.add_argument('--bdt_multitrial', "-bm", action='store_true', default=False)
    parser.add_argument('--modifications_config', "-m", metavar='text', default='')
    parser.add_argument("--outputdir", "-o", metavar="text", default=".", help="output directory")
    args = parser.parse_args()

    produce_trials_cfgs(args.input_config,
                        args.modifications_config,
                        args.outputdir,
                        'fit' if args.fit_multitrial else 'bdt')
