import yaml, os, sys
from itertools import product
import pathlib as PATH
from alive_progress import alive_bar
import argparse

PROJECT_ROOT = os.path.abspath(os.path.join(os.getcwd(), '..', '..'))
print(f"[INFO] Project root: {PROJECT_ROOT}")

def run_command(cmd, log, label):
    log.parent.mkdir(parents=True, exist_ok=True)
    print(f"[{label}] {cmd} > {log} 2>&1")
    out = os.system(f"{cmd} > {log} 2>&1")
    if out != 0:
        print(f"[{label}] WARNING: Command exited with code {out}")

def task_modify_config(config_path, inner_edge, outer_edge, task_LM=None):
    with open(config_path) as f:
        config = yaml.safe_load(f)
    config['deltaEtaBins'] = [
        [-float(outer_edge), -float(inner_edge)],
        [float(inner_edge), float(outer_edge)]
    ]
    if config.get("method") == "MassBinning":
        config['rebinDeltaPhi'] = 4
    if task_LM is not None:
        config['pathFileSE'] = task_LM['pathFileSE']
        config['pathFileME'] = task_LM['pathFileME']
        config['pathFileMass'] = task_LM['pathFileMass']
        config['outdir'] = task_LM['outdir']
        config['nDeltaPhiBins'] = task_LM['nDeltaPhiBins']
        if config.get("method") == "MassBinning":
            config['rebinDeltaPhi'] = 4

    # =======+++++++=======+++++++=======+++++++======= !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    final_suffix = f'AppMass' # the fianl suffix is here /_\|/_\|/_\|/_/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|/_\|
    # =======+++++++=======+++++++=======+++++++======= !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    config['suffix'] = f'{inner_edge.replace(".", "d")}_{outer_edge.replace(".", "d")}_{final_suffix}'
    output_dir = PATH.Path(config['outdir']) / f"CorrelExtract_{config['suffix']}"
    output_dir.mkdir(parents=True, exist_ok=True)
    out_config_path = output_dir / config_path.name.replace('.yaml', f'_{config["suffix"]}.yaml')
    with open(out_config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    print(f"Generated config file: {out_config_path}")
    return out_config_path

def task_extract_correl(config_path):
    run_command(f"cd {PROJECT_ROOT}/correlations/PostProcessing && python3 {PROJECT_ROOT}/correlations/PostProcessing/ExtractOutputCorrel.py {config_path}",
                config_path.parent / f"log_extract_{config_path.stem}.txt", "Extract")

def task_fit_mass(config_path):
    run_command(f"cd {PROJECT_ROOT}/src && python3 {PROJECT_ROOT}/src/ry_interface.py {config_path}",
                config_path.parent / f"log_fit_mass_{config_path.stem}.txt", "FitMass")

def task_fit_correl(config_path):
    run_command(f"cd {PROJECT_ROOT}/correlations/PostProcessing && python3 {PROJECT_ROOT}/correlations/PostProcessing/FitCorrel.py {config_path}",
                config_path.parent / f"log_fit_correl_{config_path.stem}.txt", "FitCorrel")

def task_extract_ry_trigger(config_path):
    run_command(f"cd {PROJECT_ROOT}/correlations/PostProcessing && python3 {PROJECT_ROOT}/correlations/PostProcessing/construct_v2_mass.py extract-ry-trigger {config_path}",
                config_path.parent / f"log_ry_trigger_{config_path.stem}.txt", "RyTrig")

def task_construct_v2_mass(config_path):
    run_command(f"cd {PROJECT_ROOT}/correlations/PostProcessing && python3 {PROJECT_ROOT}/correlations/PostProcessing/construct_v2_mass.py build-mass-v2 {config_path}",
                config_path.parent / f"log_construct_v2_mass_{config_path.stem}.txt", "MassV2")

def task_simfit(config_path):
    with open(config_path) as f:
        config = yaml.safe_load(f)
    if config.get("method") != "MassBinning":
        print("[SimFit] Skipping since method is not MassBinning")
        return
    outdir = PATH.Path(config['outdir'])
    suffix = config['suffix']
    mass_v2_dir = outdir / f"CorrelExtract_{suffix}" / "CorrelationFitResults" / "MassVsV2"
    pt_bins_cand = [float(x) for x in config['ptBinsCand']]
    pt_bins_had = [float(x) for x in config['ptBinsHad']]
    fit_config = config.get('fitConfig', {})
    Dmeson = config.get("Dmeson", "Dzero")
    get_vn_script = f"{PROJECT_ROOT}/src/get_vn_vs_mass.py"
    if not os.path.exists(get_vn_script):
        print("[SimFit] WARNING: get_vn_vs_mass.py not found, skipping SimFit")
        return
    simfit_config = {
        "Dmeson": Dmeson,
        "ptbins": pt_bins_cand,
        "centrality": config.get("centrality", "k020"),
        "v2extraction": {
            "SgnFunc": fit_config.get("SgnFunc", "kGaus"),
            "BkgFunc": fit_config.get("BkgFunc", "kExpo"),
            "BkgFuncVn": fit_config.get("BkgFuncVn", "kLin"),
            "Sigma": fit_config.get("Sigma", [0.02] * (len(pt_bins_cand)-1)),
            "MassFitRanges": fit_config.get("MassFitRanges", [[pt_bins_cand[0], pt_bins_cand[-1]]] * (len(pt_bins_cand)-1)),
            "Rebin": fit_config.get("Rebin", 4),
            "FixSigma": fit_config.get("FixSigma", 0),
            "FixMean": fit_config.get("FixMean", 0),
            "InclSecPeak": fit_config.get("InclSecPeak", 0),
            "enableRef": fit_config.get("enableRef", False),
            "ReflFunc": fit_config.get("ReflFunc", "2gaus")
        }
    }
    for i_pt_had in range(len(pt_bins_had)-1):
        pt_had_label = f"PtAssoc{int(pt_bins_had[i_pt_had]*10):02d}to{int(pt_bins_had[i_pt_had+1]*10):02d}"
        mass_v2_file = mass_v2_dir / f"InvMassVsV2_{pt_had_label}.root"
        if not mass_v2_file.exists():
            print(f"[SimFit] WARNING: {mass_v2_file} not found, skipping {pt_had_label}")
            continue
        simfit_config_file = mass_v2_dir / f"simfit_config_{pt_had_label}.yaml"
        with open(simfit_config_file, 'w') as f:
            yaml.dump(simfit_config, f, default_flow_style=False)
        log_file = outdir / f"CorrelExtract_{suffix}" / f"log_simfit_{pt_had_label}.txt"
        print(f"[SimFit] Running SimFit for {pt_had_label}...")
        os.system(f"python3 {get_vn_script} {simfit_config_file} {mass_v2_file} -b > {log_file} 2>&1")

import ROOT
def produce_final_results(config_path):
    with open(config_path) as f:
        config = yaml.safe_load(f)
    if config.get("method") != "DeltaPhiBinning":
        print("[FinalResults] Skipping since method is not DeltaPhiBinning")
        return
    outdir = PATH.Path(config['outdir'])
    suffix = config['suffix']
    input_file = outdir / f"CorrelExtract_{suffix}" / "CorrelationFitResults" / "Output_CorrelationFitting_Root" / "CorrPhiD0_FinalPlots.root"
    output_file = outdir / f"CorrelExtract_{suffix}" / "final_results.root"
    print(f"[FinalResults] Processing {input_file} -> {output_file}")
    if not input_file.exists():
        print(f"[FinalResults] WARNING: Input file {input_file} does not exist, skipping")
        return
    f = ROOT.TFile.Open(str(input_file))
    for k in f.GetListOfKeys():
        print(f"[FinalResults] Processing histogram: {k.GetName()}")
        h = f.Get(k.GetName())
        h.SetDirectory(0)
        h.Scale(1/0.07)
        o = ROOT.TFile.Open(str(output_file), "UPDATE")
        o.cd()
        h.Write(k.GetName(), ROOT.TObject.kOverwrite)
        o.Close()
    f.Close()

def main(config_path=None, inner_edges=None, outer_edges=None, special_cases=None):
    if config_path is None:
        config_path = f"{PROJECT_ROOT}/correlations/PostProcessing/config_CorrAnalysis_v2_010_negDeta.yaml"
    if inner_edges is None:
        inner_edges = ["0.8"]
    if outer_edges is None:
        outer_edges = ["1.3"]
    if special_cases is None:
        special_cases = []
    # ===== MAIN LOOP =====
    config = config_path
    # list_inner_edges = ["0.", "0.2", "0.4", "0.6"]
    list_inner_edges = ["0.8"]
    list_outer_edges = ["1.3"]
    special_cases = []

    with open(config, 'r') as f:
        base_config = yaml.safe_load(f)
    method = base_config.get("method", "DeltaPhiBinning")
    print(f"[INFO] Using method: {method}")

    all_pairs = list(product(list_inner_edges, list_outer_edges)) + special_cases
    with alive_bar(len(all_pairs), title="Processing deltaEta configs", length=30) as bar:
        for inner_edge, outer_edge in all_pairs:
            if float(outer_edge) <= float(inner_edge) + 0.1:
                print(f"[INFO] Skipping invalid config: inner_edge={inner_edge}, outer_edge={outer_edge}")
                continue
            print(f"\n{'='*60}\n  inner={inner_edge}, outer={outer_edge}\n{'='*60}")
            out_config_path = task_modify_config(PATH.Path(config), inner_edge, outer_edge)
            
            # === Step 1: Extract correlations ===
            print(">>> Step 1: Extract")
            task_extract_correl(out_config_path)

            if method == "MassBinning":
                print(">>> Step 2a: Extract ry_trigger")
                task_extract_ry_trigger(out_config_path)
            if method == "DeltaPhiBinning":
                print(">>> Step 2: Fit mass")
                task_fit_mass(out_config_path)
            
            ## === LM template ===
            if 'task_LM' in base_config and base_config['task_LM'].get('do', False):
                out_config_path_lm = task_modify_config(PATH.Path(config), inner_edge, outer_edge, base_config['task_LM'])
                print(">>> LM Template: Extract correlations")
                task_extract_correl(out_config_path_lm)
                if method == "MassBinning":
                    print(">>> LM Template: Extract ry_trigger")
                    task_extract_ry_trigger(out_config_path_lm)
                if method == "DeltaPhiBinning":
                    print(">>> LM Template: Fit mass")
                    task_fit_mass(out_config_path_lm)

            #== Step 3: Fit correlations ===
            print(">>> Step 3: Fit correlations")
            task_fit_correl(out_config_path)
            
            # === MassBinning-only steps ===
            if method == "MassBinning":
                print(">>> Step 4: Build v2-vs-mass")
                task_construct_v2_mass(out_config_path)
                print(">>> Step 5: Simultaneous fit")
                task_simfit(out_config_path)
            
            # === DeltaPhiBinning-only steps ===
            if method == "DeltaPhiBinning":
                print(">>> Final: Produce final results")
                produce_final_results(out_config_path)
            bar()

    print("\n[INFO] All done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the full analysis pipeline")
    parser.add_argument("--config", default=f"config_CorrAnalysis_v2_010_negDeta.yaml",
                        help="Path to the base config file (default: config_CorrAnalysis_v2_010_negDeta.yaml)")
    parser.add_argument("--inner", nargs="+", default=["0.8"],
                        help="Inner delta-eta edges (default: 0.8)")
    parser.add_argument("--outer", nargs="+", default=["1.3"],
                        help="Outer delta-eta edges (default: 1.3)")
    parser.add_argument("--special", nargs=2, action="append", default=[],
                        metavar=("INNER", "OUTER"),
                        help="Special (inner, outer) pairs, e.g. --special 0. 1.3")
    args = parser.parse_args()

    main(args.config, args.inner, args.outer, [(s[0], s[1]) for s in args.special])