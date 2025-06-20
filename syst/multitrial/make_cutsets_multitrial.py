import yaml
import ROOT
import copy
import os
import numpy as np
import argparse
import itertools
import sys
sys.path.append(f"{os.path.dirname(os.path.abspath(__file__))}/../../utils")
from utils import check_dir, logger

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

    pt_bin_config['simfit']['Sigma'] = [pt_bin_config['simfit']['Sigma'][iPtBin]]

    return pt_bin_config

def modify_yaml_bdt(config_flow, config_mod, output_dir):
    
    with open(config_flow, 'r') as CfgFlow:
        cfg_flow = yaml.safe_load(CfgFlow)
    with open(config_mod, 'r') as CfgMod:
        cfg_mod = yaml.safe_load(CfgMod)

    multitrial_dir = f"{output_dir}/syst/multitrial/"
    check_dir(multitrial_dir)
    os.makedirs(multitrial_dir, exist_ok=True)
    os.makedirs(f"{multitrial_dir}/config_history/", exist_ok=True)

    with open(os.path.join(f"{multitrial_dir}/config_history/", f"config_reference.yml"), 'w') as out_file:
        yaml.dump(cfg_flow, out_file, default_flow_style=False)
    with open(os.path.join(f"{multitrial_dir}/config_history/", f"config_modifications.yml"), 'w') as out_file:
        yaml.dump(cfg_mod, out_file, default_flow_style=False)

    for _, pt_bin in enumerate(cfg_mod['ptbins']):
        ptmin = pt_bin['range'][0]
        ptmax = pt_bin['range'][1]
        pt_bin_index = cfg_flow['ptbins'].index(ptmin)
        ref_config_ptbin = get_reference_config_pt_bin(cfg_flow, pt_bin_index)
        multitrial_pt_dir = f"{multitrial_dir}/pt_{int(ptmin*10)}_{int(ptmax*10)}/"
        
        os.makedirs(multitrial_pt_dir, exist_ok=True)
        output_file = os.path.join(multitrial_pt_dir, f"config_reference.yml")
        with open(output_file, 'w') as out_file:
            yaml.dump(ref_config_ptbin, out_file, default_flow_style=False)

        varied_configs = {}
        for var in pt_bin["multitrial"]:
            if var in ref_config_ptbin:
                varied_configs[var] = list(set(ref_config_ptbin[var] + pt_bin["multitrial"][var]))
            else:
                varied_configs[var] = pt_bin["multitrial"][var]

        # Perform itertools.product to get all combinations
        keys = list(varied_configs.keys())
        values = list(varied_configs.values())
        combinations = list(itertools.product(*values))
        config_variants = [dict(zip(keys, combination)) for combination in combinations] 

        # Print results
        for idx, variant in enumerate(config_variants):
            cfg_variant = copy.deepcopy(ref_config_ptbin)
            cfg_variant["iTrial"] = idx

            cfg_variant["outdir"] = f"{multitrial_dir}trials/"
            if pt_bin['MaxChi2'] > 10:
                logger(f"MaxChi2 is set to {pt_bin['MaxChi2']} for pt bin {ptmin}-{ptmax}, please consider a tight selection!", "WARNING")
            cfg_variant["MaxChi2"] = pt_bin['MaxChi2']
            if pt_bin['MinSignificance'] < 5:
                logger(f"MinSignificance is set to {pt_bin['MinSignificance']}, please consider a higher lower limit!", "WARNING")
            cfg_variant["MinSignificance"] = pt_bin['MinSignificance']
            if pt_bin['MaxSignificance'] > 300:
                logger(f"MaxSignificance is set to {pt_bin['MaxSignificance']}, please consider a lower upper limit!", "WARNING")
            cfg_variant["MaxSignificance"] = pt_bin['MaxSignificance']
            for varied_var in variant:
                if ref_config_ptbin.get(varied_var):
                    if isinstance(ref_config_ptbin[varied_var], list):
                        cfg_variant[varied_var] = [variant[varied_var]]
                    else:
                        cfg_variant[varied_var] = variant[varied_var]

            cfg_variant["simfit"]["MassFitRanges"] = [[variant["MassMin"], variant["MassMax"]]]
            cfg_variant["simfit"]["BkgFunc"] = [variant["BkgFunc"]]
            cfg_variant["simfit"]["SgnFunc"] = [variant["SgnFunc"]]
            cfg_variant["simfit"]["BkgFuncVn"] = [variant["BkgFuncVn"]]
            cfg_variant["projections"]["inv_mass_bins"] = [np.arange(variant['MassMin'], variant['MassMax'] + variant['inv_mass_bins_steps'], variant['inv_mass_bins_steps']).tolist()]
            if cfg_flow['Dmeson'] == 'Dplus' and (variant["MassMin"] > 1.75 or variant["MassMax"] < 1.95):
                # else the region used for the prefit of the bkg function would be empty, leading to a crash
                cfg_variant["NSigma4SB"] = [2]

            os.makedirs(os.path.join(f"{multitrial_pt_dir}/trial_{idx}"), exist_ok=True)
            output_file = os.path.join(f"{multitrial_pt_dir}/trial_{idx}", f"config_trial_{idx}.yml")
            
            # No preprocess and operations are controlled in the bash script
            cfg_variant.pop('operations', None)
            cfg_variant.pop('preprocess', None)

            with open(output_file, 'w') as out_file:
                yaml.dump(cfg_variant, out_file, default_flow_style=False, sort_keys=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('input_config', metavar='text', default='config_Ds_Fit.yml')
    parser.add_argument('--modifications_config', "-m", metavar='text', default='')
    parser.add_argument("--outputdir", "-o", metavar="text", default=".", help="output directory")
    args = parser.parse_args()

    modify_yaml_bdt(args.input_config,
                    args.modifications_config,
                    args.outputdir)
