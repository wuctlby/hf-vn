import os
import sys
import argparse
import yaml
import concurrent.futures
import time
import subprocess
# from concurrent.futures import ProcessPoolExecutor
work_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f"{work_dir}/utils")
from utils import check_dir, logger

paths = {
	"Preprocess": os.path.join(work_dir, "./src/pre_process.py"),
	"YamlCuts": os.path.join(work_dir, "./src/make_cutsets_cfgs.py"),
	"Projections": os.path.join(work_dir, "./src/proj_thn.py"),
	"Efficiencies": os.path.join(work_dir, "./src/compute_efficiencies.py"),
	"GetVnVsMass": os.path.join(work_dir, "./src/get_vn_vs_mass.py"),
	"CutVariation": os.path.join(work_dir, "./src/cut_variation.py"),
	"DataDrivenFraction": os.path.join(work_dir, "./src/data_driven_fraction.py"),
	"GetV2VsFrac": os.path.join(work_dir, "./src/get_v2_vs_frac.py"),
}

def make_yaml(flow_config, outdir, correlated=False):
	logger("YAML file will be created", level="INFO")
	check_dir(f"{outdir}/cutsets")

	method = "--correlated" if correlated else ""
	cmd = (
		f'python3 {paths["YamlCuts"]} {flow_config} -o {outdir} {method}'
	)

	logger(f"{cmd}", level="COMMAND")
	os.system(cmd)

#! I don't see any reason to use bool for correlated and combined here, this is just a flag to decide whether to add a suffix
def project(flow_config, outdir, nworkers, mCutSets):
	logger("Projections will be performed", level="INFO")
	os.makedirs(f"{outdir}/proj", exist_ok=True)
	logger(f"Output directory for projections: {outdir}/proj", level="INFO")

	def run_projections(i, mCutSets):
		"""Run sparse projection for a given cutset index."""
		iCutSets = f"{i:02d}"
		print(f"\033[32mProcessing cutset {iCutSets}...\033[0m")

		config_cutset = f"{outdir}/cutsets/cutset_{iCutSets}.yml"
		cmd = (
			f"python3 {paths['Projections']} {flow_config} -cc {config_cutset} --mCutSets {mCutSets}"
		)
		print(f"\033[32m{cmd}\033[0m")
		os.system(cmd)

	with concurrent.futures.ThreadPoolExecutor(max_workers=nworkers) as executor:
		results_proj = list(executor.map(run_projections, range(mCutSets), [mCutSets] * mCutSets))

def efficiencies(flow_config, outdir, nworkers, mCutSets):
	print("\033[32mINFO: Efficiencies will be computed\033[0m")
	check_dir(f"{outdir}/eff")

	def run_efficiency(i):
		"""Run efficiency computation for a given cutset index."""
		iCutSet = f"{i:02d}"
		print(f"\033[32mProcessing cutset {iCutSet}...\033[0m")

		proj_cutset = f"{outdir}/proj/proj_{iCutSet}.root"
		if not os.path.exists(proj_cutset):
			logger(f"Projection file {proj_cutset} does not exist, skipping cutset {iCutSet}", level="ERROR")
			return
		cmd = (
			f"python3 {paths['Efficiencies']} {flow_config} {proj_cutset} -b"
		)
		print(f"\033[32m{cmd}\033[0m")
		os.system(cmd)

	with concurrent.futures.ThreadPoolExecutor(max_workers=nworkers) as executor:
		results_eff = list(executor.map(run_efficiency, range(mCutSets)))

def get_vn_vs_mass(flow_config, outdir, nworkers, mCutSets, batch=False):
	print("\033[32mINFO: Fit v2 vs mass will be performed\033[0m")
	check_dir(f"{outdir}/raw_yields")

	def run_fit(i):
		"""Run simultaneous fit for a given cutset index."""
		iCutSets = f"{i:02d}"
		print(f"\033[32mProcessing cutset {iCutSets}...\033[0m")

		proj_cutset = f"{outdir}/proj/proj_{iCutSets}.root"
		cmd = (
			f"python3 {paths['GetVnVsMass']} {flow_config} {proj_cutset}"
		)
		if batch:
			cmd += " --batch"
		print(f"\033[32m{cmd}\033[0m")
		os.system(cmd)

	with concurrent.futures.ThreadPoolExecutor(max_workers=nworkers) as executor:
		results_fit = list(executor.map(run_fit, range(mCutSets)))

def cut_variaion(flow_config, outdir, correlated, combined=False, operations=None, batch=False):
	check_dir(f"{outdir}/cutVar")
	if correlated:
		logger("Cut variation will be performed", level="INFO")

		ry_path = f"{outdir}/raw_yields"
		eff_path = f"{outdir}/eff"
		method = "--correlated" if correlated else ""

		cmd = (
			f"python3 {paths['CutVariation']} {flow_config} {ry_path} {eff_path} {method}"
		)
		if batch:
			cmd += " --batch"
		logger(f"{cmd}", level="COMMAND")
		os.system(cmd)

	if combined:
		outdirCorr = os.path.join(os.path.dirname(outdir), os.path.basename(outdir).replace("_combined", "_correlated"))
		if os.path.exists(f"{outdirCorr}/cutVar/cutVar.root"):
			logger("Cut variation with correlated cuts found!", level="INFO")
			cmd = (
				f"cp {outdirCorr}/cutVar/cutVar.root {outdir}/cutVar/cutVar.root"
			)
			logger(f"{cmd}", level="COMMAND")
			os.system(cmd)
		else:
			logger(f"Cut variation with correlated cuts not found in {outdirCorr}, running correlated cut variation", level="WARNING")

			logger("Running analysis with correlated cut variation", level="WARNING")
			logger("========================================================================================", level="WARNING")
			import copy
			operationsCorr = copy.deepcopy(operations)
			for key in operationsCorr:
				if key in ['do_cut_variation', 'make_yaml', 'proj_mc', 'proj_data', 'efficiencies', 'get_vn_vs_mass']:
					operationsCorr[key] = True
				else:
					operationsCorr[key] = False
			os.makedirs(os.path.join(outdir, "config_flow"), exist_ok=True)
			flow_configCorr = os.path.join(outdir, "config_flow", f"corr_{os.path.basename(flow_config)}")
			with open(flow_config, 'r') as cfgFlow:
				config = yaml.safe_load(cfgFlow)
				config['operations'] = operationsCorr
				config['outdir'] = outdirCorr #! make sure all output directories are set from function arguments not from config file
				with open(flow_configCorr, 'w') as cfgFlow:
					yaml.dump(config, cfgFlow, default_flow_style=True)
			logger(f"Using temporary flow config: {flow_configCorr}", level="WARNING")
			run_correlated_cut_variation(flow_configCorr, operationsCorr, nworkers, outdirCorr)
			logger("========================================================================================", level="WARNING")

			logger("Copying cut variation with correlated cuts to combined cut variation", level="WARNING")
			cmd = (
				f"cp {outdirCorr}/cutVar/cutVar.root {outdir}/cutVar/cutVar.root"
			)
			logger(f"{cmd}", level="COMMAND")
			os.system(cmd)

def data_driven_fraction(outdir, syst_frac, batch=False):
	logger("Data driven fraction will be performed", level="INFO")
	check_dir(f"{outdir}/frac")
 
	cutvar_file = f"{outdir}/cutVar/cutVar.root"
	eff_path = f"{outdir}/eff"

	cmd = (
		f"python3 {paths['DataDrivenFraction']} {cutvar_file} {eff_path}"
	)
	if syst_frac:
		cmd += " --syst_frac"
	if batch:
		cmd += " --batch"
	logger(f"{cmd}", level="COMMAND")
	os.system(cmd)
 
def get_v2_vs_frac(flow_config, outdir, correlated=False, batch=False):
	logger("Fit v2 vs fd fraction will be performed", level="INFO")
	check_dir(f"{outdir}/v2")

	ry_path = f"{outdir}/raw_yields"
	frac_path = f"{outdir}/frac"
	cmd = (
		f"python3 {paths['GetV2VsFrac']} {flow_config} {ry_path} {frac_path}"
	)
	if correlated:
		cmd += " --correlated"
	if batch:
		cmd += " --batch"
	logger(f"{cmd}", level="COMMAND")
	os.system(cmd)

def run_correlated_cut_variation(flow_config, operations, nworkers, outdir):

	#___________________________________________________________________________________________________________________________
	# make yaml file
	if operations.get('make_yaml', False):
		make_yaml(flow_config, outdir, correlated=True)
	else:
		print("\033[33mWARNING: Make yaml will not be performed\033[0m")

	mCutSets = len([f for f in os.listdir(f"{outdir}/cutsets") if os.path.isfile(os.path.join(f"{outdir}/cutsets", f))])
	print(f"mCutSets: {mCutSets}")

	#___________________________________________________________________________________________________________________________
	# Projection for MC and apply the ptweights
	if operations.get('proj_mc', False) or operations.get('proj_data', False):
		project(flow_config, outdir, nworkers, mCutSets)
	else:
		print("\033[33mWARNING: Projections will not be performed\033[0m")

	#___________________________________________________________________________________________________________________________
	# Efficiencies
	if operations.get('efficiencies', False):
		efficiencies(flow_config, outdir, nworkers, mCutSets)
	else:
		print("\033[33mWARNING: Efficiencies will not be computed\033[0m")

	#___________________________________________________________________________________________________________________________
	# Simultaneous fit
	if operations.get('get_vn_vs_mass', False):
		get_vn_vs_mass(flow_config, outdir, nworkers, mCutSets, True)
	else:
		print("\033[33mWARNING: Fit v2 vs mass will not be performed\033[0m")

	#___________________________________________________________________________________________________________________________
	# Cut variation
	if operations.get('do_cut_variation'):
		cut_variaion(flow_config, outdir, correlated=True, combined=False, operations=operations, batch=True)
	else:
		logger("Cut variation will not be performed", level="WARNING")
  
	#___________________________________________________________________________________________________________________________
 	# Data driven fraction
	if operations.get('data_driven_fraction', False):
		data_driven_fraction(outdir, syst_frac=False, batch=True)
	else:
		logger("Data driven fraction will not be performed", level="WARNING")

	#___________________________________________________________________________________________________________________________
 	# linear fit of vn vs fd fraction
	if operations.get('get_v2_vs_frac'):
		get_v2_vs_frac(flow_config, outdir, correlated=True, batch=True)
	else:
		logger("Fit v2 vs fd fraction will not be performed", level="WARNING")

def run_combined_cut_variation(flow_config, operations, nworkers, outdir):

	#___________________________________________________________________________________________________________________________
	# make yaml file
	if operations.get('make_yaml', False):
		make_yaml(flow_config, outdir, correlated=False)
	else:
		print("\033[33mWARNING: Make yaml will not be performed\033[0m")

	mCutSets = len([f for f in os.listdir(f"{outdir}/cutsets") if os.path.isfile(os.path.join(f"{outdir}/cutsets", f))])
	print(f"mCutSets: {mCutSets}")

	#___________________________________________________________________________________________________________________________
	# Projection for MC and apply the ptweights
	if operations.get('proj_mc', False) or operations.get('proj_data', False):
		project(flow_config, outdir, nworkers, mCutSets)
	else:
		print("\033[33mWARNING: Projections will not be performed\033[0m")

	#___________________________________________________________________________________________________________________________
	# Efficiencies
	if operations.get('efficiencies', False):
		efficiencies(flow_config, outdir, nworkers, mCutSets)
	else:
		print("\033[33mWARNING: Efficiencies will not be computed\033[0m")

	#___________________________________________________________________________________________________________________________
	# Simultaneous fit
	if operations.get('get_vn_vs_mass', False):
		get_vn_vs_mass(flow_config, outdir, nworkers, mCutSets, True)
	else:
		print("\033[33mWARNING: Fit v2 vs mass will not be performed\033[0m")

	#___________________________________________________________________________________________________________________________
	# Cut variation
	if operations.get('do_cut_variation', False):
		cut_variaion(flow_config, outdir, correlated=False, combined=True, operations=operations, batch=True)
	else:
		logger("Cut variation will not be performed", level="WARNING")
  
	#___________________________________________________________________________________________________________________________
	# Data driven fraction
	if operations.get('data_driven_fraction', False):
		data_driven_fraction(outdir, syst_frac=True, batch=True)
	else:
		logger("Data driven fraction will not be performed", level="WARNING")
  
	#___________________________________________________________________________________________________________________________
	# linear fit of vn vs fd fraction
	if operations.get('get_v2_vs_frac'):
		get_v2_vs_frac(flow_config, outdir, correlated=False, batch=True)
	else:
		logger("Fit v2 vs fd fraction will not be performed", level="WARNING")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Arguments')
	parser.add_argument('flow_config', metavar='text', default='config_flow_d0.yml', help='configuration file')
	parser.add_argument("--correlated", "-corr", action="store_true", help="perform correlated analysis")
	parser.add_argument("--combined", "-comb", action="store_true", help="perform combined analysis")
	args = parser.parse_args()

	start_time = time.time()

	# Load and copy the configuration file
	with open(args.flow_config, 'r') as cfgFlow:
		config = yaml.safe_load(cfgFlow)

	operations = config['operations']
	nworkers = config['nworkers']
	if args.correlated:
		outdir = f"{config['outdir']}/cutvar_{config['suffix']}" + "_correlated"
	else:
		outdir = f"{config['outdir']}/cutvar_{config['suffix']}" + "_combined"
	os.system(f"mkdir -p {outdir}")

	# copy the configuration file
	nfile = 0
	os.makedirs(f'{outdir}/config_flow', exist_ok=True)
	while os.path.exists(f'{outdir}/config_flow/{os.path.splitext(os.path.basename(args.flow_config))[0]}_{config["suffix"]}_{nfile}.yml'):
		nfile = nfile + 1
	os.system(f'cp {args.flow_config} {outdir}/config_flow/{os.path.splitext(os.path.basename(args.flow_config))[0]}_{config["suffix"]}_{nfile}.yml')

	if operations.get('preprocess_data') or operations.get('preprocess_mc'):
		print("\033[32mINFO: Preprocess will be performed\033[0m")
		os.system(f"python3 {paths['Preprocess']} {args.flow_config}")
	else:
		print("\033[33mWARNING: Preprocess will not be performed\033[0m")

	if not args.correlated and not args.combined:
		print("\033[33mWARNING: No cut variation will be performed\033[0m")
		sys.exit(0)
	if args.correlated:
		run_correlated_cut_variation(args.flow_config, operations, nworkers, outdir)
	if args.combined:
		run_combined_cut_variation(args.flow_config, operations, nworkers, outdir)

	end_time = time.time()
	execution_time = end_time - start_time
	print(f"\033[34mTotal execution time: {execution_time:.2f} seconds\033[0m")
