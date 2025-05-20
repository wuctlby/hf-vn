import os
import sys
import argparse
import yaml
import concurrent.futures
import time
import subprocess
from concurrent.futures import ProcessPoolExecutor
work_dir = os.path.dirname(os.path.realpath(__file__))
from utils import check_dir

paths = {
	"YamlCuts": os.path.join(work_dir, "./make_cutsets_cfgs.py"),
}

def make_yaml(config, outdir, correlated=False, combined=False):
	print("\033[32mINFO: Make yaml will be performed\033[0m")
	check_dir(f"{outdir}/cutsets")
	print(f"\033[32mpython3 {paths['YamlCuts']} {config} -o {outdir} --correlated\033[0m")
	if correlated:
		os.system(f"python3 {paths['YamlCuts']} {config} -o {outdir} --correlated")
	if combined:
		os.system(f"python3 {paths['YamlCuts']} {config} -o {outdir}")

def run_correlated_cut_variation(config, operations, nworkers, outdir):

	#___________________________________________________________________________________________________________________________
	# make yaml file
	make_yaml(config, outdir, correlated=True) if operations['make_yaml'] else print("\033[33mWARNING: Make yaml will not be performed\033[0m")

def run_combined_cut_variation(config, operations, nworkers, outdir):

	#___________________________________________________________________________________________________________________________
	# make yaml file
	make_yaml(config, outdir, combined=True) if operations['make_yaml'] else print("\033[33mWARNING: Make yaml will not be performed\033[0m")

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
	while os.path.exists(f'{outdir}/config_flow/{os.path.splitext(os.path.basename(args.flow_config))[0]}_{config['suffix']}_{nfile}.yml'):
		nfile = nfile + 1
	os.system(f'cp {args.flow_config} {outdir}/config_flow/{os.path.splitext(os.path.basename(args.flow_config))[0]}_{config['suffix']}_{nfile}.yml')

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
