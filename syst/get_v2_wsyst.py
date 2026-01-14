import ROOT
import yaml
import sys
import os
import numpy as np
import argparse
from ROOT import TFile, TGraphAsymmErrors, kRed, kOrange, kAzure, kBlack, kGreen, kOpenCircle
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(script_dir, '..', 'utils'))
from utils import get_particle_info, logger
from StyleFormatter import SetObjectStyle

def th1_to_tgraph(hist):
    """
    Convert a TH1 histogram into a TGraph.
    
    Parameters:
        hist (ROOT.TH1): Input histogram
    
    Returns:
        ROOT.TGraphAsymmErrors: Converted graph
    """
    n_points = hist.GetNbinsX()
    graph = TGraphAsymmErrors(n_points)
    
    for i in range(1, n_points + 1):
        x = hist.GetBinCenter(i)
        y = hist.GetBinContent(i)
        graph.SetPoint(i - 1, x, y)
        graph.SetPointError(i - 1, hist.GetBinWidth(i)/2, hist.GetBinWidth(i)/2, hist.GetBinError(i), hist.GetBinError(i))
    
    return graph

def assign_syst(histo, syst, name, isrelative=False):
    """
    Assign systematic uncertainties to a TGraph from a TH1 histogram.
    
    Parameters:
        histo (ROOT.TH1): Input histogram
        syst (list): List of systematic uncertainties
        name (string): Name of the output graph
        isrelative (bool): If True, treat uncertainties as relative
    
    Returns:
        ROOT.TGraphAsymmErrors: Graph with systematic uncertainties applied
    """
    gisto = th1_to_tgraph(histo)
    gisto.SetName(name)
    for ipoint in range(histo.GetNbinsX()):
        unc = syst[ipoint]
        if isrelative:
            print(f'WARNING: {name} considered as relative uncertainty!')
            unc *= gisto.GetPointY(ipoint)
        gisto.SetPointError(ipoint, 0.3, 0.3, unc, unc)

    return gisto

def compute_vn_with_syst(in_file_name, config, dmeson, origin,  centrality, output, suffix):
    """
    Compute D v2 with systematic uncertainties.
    
    Parameters:
        in_file_name (str): Input ROOT file path
        config (str): YAML configuration file path
        dmeson (str): D meson specie
        origin: (str): D meson origin
        centrality (str): Centrality class
        output (str): Output directory path
        suffix (str): Suffix for output files
    """
    with open(config, 'r') as yml_cfg_file:
        config = yaml.load(yml_cfg_file, yaml.FullLoader)

    infile_vn = TFile.Open(in_file_name)
    if not infile_vn or not infile_vn.IsOpen():
        print(f"ERROR: Unable to open file {in_file_name}")
        return
    
    hvn_stat = infile_vn.Get('hV2VsPtPrompt') if origin == 'prompt' else infile_vn.Get('hV2VsPtFD')
    if not hvn_stat:
        print("ERROR: hV2VsPtFD does not exist! Exiting.")
        return
    
    gvn_stat = th1_to_tgraph(hvn_stat)
    gvn_stat.SetName("gvn_prompt_stat") if origin == 'prompt' else gvn_stat.SetName("gvn_np_stat")
    SetObjectStyle(gvn_stat, color=kBlack, markerstyle=kOpenCircle)
    
    fd_syst = config[dmeson][centrality]['fd_syst']
    fit_syst = config[dmeson][centrality]['fit_syst']
    reso_syst = config[dmeson][centrality]['reso_syst']
    pt_bins = gvn_stat.GetN()
    
    if not isinstance(fd_syst, list):
        fd_syst = [fd_syst] * pt_bins
    if not isinstance(fit_syst, list):
        fit_syst = [fit_syst] * pt_bins
    if not isinstance(reso_syst, list):
        reso_syst = [reso_syst] * pt_bins
    
    if len(fd_syst) != pt_bins or len(fit_syst) != pt_bins or len(reso_syst) != pt_bins:
        logger("Error: Inconsistent number of bins! Exiting.", level='ERROR')
    
    abs_reso_syst = [reso_syst[ipt] * hvn_stat.GetBinContent(ipt + 1) for ipt in range(pt_bins)]
    tot_syst = [np.sqrt(fd**2 + fit**2 + reso**2) for fd, fit, reso in zip(fd_syst, fit_syst, abs_reso_syst)]
    
    gvn_fdsyst = assign_syst(hvn_stat, fd_syst, 'fd_syst')
    gvn_fitsyst = assign_syst(hvn_stat, fit_syst, 'fit_syst')
    gvn_resosyst = assign_syst(hvn_stat, reso_syst, 'reso_syst', isrelative=True)
    gvn_totsyst = assign_syst(hvn_stat, tot_syst, 'tot_syst')
    
    SetObjectStyle(gvn_fdsyst, markerstyle=kOpenCircle, fillalpha=0.2, color=kAzure+6)
    SetObjectStyle(gvn_fitsyst, markerstyle=kOpenCircle, fillalpha=0.2, color=kGreen+1)
    SetObjectStyle(gvn_resosyst, markerstyle=kOpenCircle, fillalpha=0.2, color=kRed-4)
    SetObjectStyle(gvn_totsyst, markerstyle=kOpenCircle, fillalpha=0.2, color=kBlack)
    
    output_dir = f'{output}/v2_wsyst/'
    os.makedirs(output_dir, exist_ok=True)
    
    outFile = TFile(f'{output_dir}/v2_prompt_wsyst_{suffix}.root', "recreate") if origin == 'prompt' else TFile(f'{output_dir}/v2_np_wsyst_{suffix}.root', "recreate")
    gvn_fdsyst.Write()
    gvn_fitsyst.Write()
    gvn_resosyst.Write()
    gvn_totsyst.Write()
    gvn_stat.Write()
    outFile.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute prompt v2 with systematic uncertainties.')
    parser.add_argument('infile', metavar='input_file', help='Path to input ROOT file')
    parser.add_argument('config_syst', metavar='config_file', help='Systematic summary YAML file')
    parser.add_argument('--dmeson', '-d', default='D0', help='D meson type')
    parser.add_argument('--origin', '-or', default='prompt', help='D meson origin (prompt/np)')
    parser.add_argument('--centrality', '-c', default='k3050', help='Centrality class')
    parser.add_argument('--suffix', '-s', default='', help='Suffix for output files')
    parser.add_argument('--outputdir', '-o', default='.', help='Output directory')
    
    args = parser.parse_args()
    
    compute_vn_with_syst(args.infile, args.config_syst, args.dmeson, args.origin, args.centrality, args.outputdir, args.suffix)