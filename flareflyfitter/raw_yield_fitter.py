import argparse
import yaml
import os
import sys
import ROOT
from ROOT import TFile, TH1, TH1D, TH1F # pylint: disable=import-error,no-name-in-module
script_dir = os.path.dirname(os.path.realpath(__file__))
os.sys.path.append(os.path.join(script_dir, '..', 'utils'))
from correlated_bkgs import get_corr_bkg
from utils import logger, get_centrality_bins
import zfit
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
import uproot
import copy
msg_service = ROOT.RooMsgService.instance()
msg_service.setGlobalKillBelow(ROOT.RooFit.FATAL)  # Only show FATAL errors (you can also use ROOT.RooFit.ERROR or INFO)
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # pylint: disable=wrong-import-position

class RawYieldFitter:
    """ 
    Fitter of invariant mass spectra to extract raw yields using the FlareFly package
    """

    def __init__(self, particle, label, minimize_flarefly=True):
        print("\nInitializing RawYieldFitter")
        if minimize_flarefly:
            self.minimize_flarefly = True
            self.maximize_roofit = False
        else:
            self.minimize_flarefly = False
            self.maximize_roofit = True
        self.rebin = None
        self.fit_range_min = None
        self.fit_range_max = None
        self.sgn_funcs = None
        self.sgn_funcs_labels = None
        self.sgn_main_label = None
        self.bkg_funcs = None
        self.bkg_funcs_labels = None
        self.corr_bkgs = None
        self.corr_bkgs_labels = None
        self.corr_bkgs_histos = None
        self.corr_bkgs_trees = None
        self.fitter = None
        self.data = None
        self.hist = None
        self.fit_result = None
        self.particle = None
        self.x_axis_label = None
        self.fit_name = f"{particle}_{label}"
        self.label = label
        self.particle_pdg = None
        self.set_particle(particle)

    def set_particle(self, particle_name):
        print(f"\nSetting particle to fit: {particle_name}")
        self.particle = particle_name
        if particle_name == "Dplus":
            self.x_axis_label = r"$M(\mathrm{\pi^+ K^- \pi^+})\ \mathrm{(GeV/}c^2)$"
            self.sgn_main_label = "DplusToPiKPi"
            self.particle_pdg = 411
        elif particle_name == "D0":
            self.x_axis_label = r"$M(\mathrm{K^- \pi^+})\ \mathrm{(GeV/}c^2)$"
            self.sgn_main_label = "D0ToKPi"
            self.particle_pdg = 421
        elif particle_name == "Ds":
            self.x_axis_label = r"$M(\mathrm{K^+ K^- \pi^+})\ \mathrm{(GeV/}c^2)$"
            self.sgn_main_label = "DsToKKPi"
            self.particle_pdg = 431
        else:
            print(f"Particle {particle_name} not recognized!")
            sys.exit(1)

    def set_name(self, fit_name):
        print(f"\nSetting fit name: {fit_name}")
        self.fit_name = fit_name

    def set_fit_range(self, fit_range_min, fit_range_max):
        print(f"\nSetting fit range: {fit_range_min} - {fit_range_max} GeV/c")
        self.fit_range_min = fit_range_min
        self.fit_range_max = fit_range_max

    def set_rebin(self, rebin):
        print(f"\nSetting rebin: {rebin}")
        self.rebin = rebin

    def set_data_to_fit_hist(self, data):
        print(f"\nSetting data to fit from histogram with limits {self.fit_range_min} - {self.fit_range_max} GeV/c"
              f" for fitter with name {self.fit_name}")
        self.hist = data
        self.data = DataHandler(data, limits=[self.fit_range_min, self.fit_range_max], rebin=self.rebin)
    
    def set_data_to_fit_df(self, data_df, var_name='fM'):
        print(f"\nSetting data to fit from dataframe with variable {var_name} and limits {self.fit_range_min} - {self.fit_range_max} GeV/c")
        self.data = DataHandler(data_df, var_name=var_name, limits=[self.fit_range_min, self.fit_range_max])

    def set_data_to_fit_tree(self, data_tree, var_name):
        print(f"\nSetting data to fit from tree with variable {var_name} and limits {self.fit_range_min} - {self.fit_range_max} GeV/c")
        self.data = DataHandler(data_tree, var_name=var_name, limits=[self.fit_range_min, self.fit_range_max])

    def add_sgn_func(self, sgn_func, label):
        print(f"\nAdding signal function with label {label}")
        if self.sgn_funcs is None:
            self.sgn_funcs = []
            self.sgn_funcs_labels = []
        if sgn_func.startswith("k"):
            if sgn_func == "kGaus":
                sgn_func = "gaussian"
            else:
                print(f"Signal function {sgn_func} not recognized!")
                sys.exit(1)
        self.sgn_funcs.append(sgn_func)
        self.sgn_funcs_labels.append(label)

    def add_bkg_func(self, bkg_func, label):
        print(f"\nAdding background function with label {label}, self.bkg_funcs: {self.bkg_funcs}, self.bkg_funcs_labels: {self.bkg_funcs_labels}")
        if self.bkg_funcs is None:
            self.bkg_funcs = []
            self.bkg_funcs_labels = []
        if bkg_func.startswith("k"):
            if bkg_func == "kExpo":
                bkg_func = "expo"
            elif bkg_func == "kLin":
                bkg_func = "chebpol1"
            elif bkg_func == "kPol2":
                bkg_func = "chebpol2"
            else:
                print(f"Background function {bkg_func} not recognized!")
                sys.exit(1)
        self.bkg_funcs.append(bkg_func)
        self.bkg_funcs_labels.append(label)

    def add_corr_bkg(self, corr_bkg_distro, corr_bkg_label, var_name='fM'):
        print(f"\nAdding correlated background with label {corr_bkg_label}")
        if self.corr_bkgs is None:
            self.corr_bkgs = []
            self.corr_bkgs_labels = []
            self.corr_bkgs_histos = []
            self.corr_bkgs_trees = []
        if isinstance(corr_bkg_distro, TH1):
            self.corr_bkgs.append("hist")
            self.corr_bkgs_histos.append(DataHandler(corr_bkg_distro, limits=[self.fit_range_min, self.fit_range_max], rebin=self.rebin))
        else:
            self.corr_bkgs.append("kde_grid")
            df = corr_bkg_distro.arrays(library="pd")
            self.corr_bkgs_histos.append(DataHandler(df, var_name=var_name, limits=[self.fit_range_min, self.fit_range_max]))

        self.corr_bkgs_labels.append(corr_bkg_label)

    def fix_sgn_par(self, sgn_func_idx, par_name, par_val):
        print(f"\nFixing signal parameter {par_name} of function index {sgn_func_idx} to value {par_val}")
        self.fitter.set_signal_initpar(sgn_func_idx, par_name, par_val, fix=True)

    def fix_bkg_par(self, bkg_func_idx, par_name, par_val):
        print(f"\nFixing background parameter {par_name} of function index {bkg_func_idx} to value {par_val}")
        self.fitter.set_background_initpar(bkg_func_idx, par_name, par_val, fix=True)

    def set_sgn_par(self, sgn_func_idx, par_name, par_val, lims):
        print(f"\nSetting signal parameter {par_name} of function index {sgn_func_idx} to value {par_val}")
        self.fitter.set_signal_initpar(sgn_func_idx, par_name, par_val, limits=lims)

    def set_bkg_par(self, bkg_func_idx, par_name, par_val, lims):
        print(f"\nSetting background parameter {par_name} of function index {bkg_func_idx} to value {par_val}")
        self.fitter.set_background_initpar(bkg_func_idx, par_name, par_val, limits=lims)

    def fix_frac_to_sgn_pdf(self, pdf_idx, anchor_pdf_idx, frac):
        print(f"\nFixing background fraction of function index {pdf_idx} to signal PDF {anchor_pdf_idx} fraction {frac}")
        self.fitter.fix_bkg_frac_to_signal_pdf(pdf_idx, anchor_pdf_idx, frac)

    def fix_frac_to_bkg_pdf(self, pdf_idx, anchor_pdf_idx, frac):
        print(f"\nFixing background fraction of function index {pdf_idx} to background PDF {anchor_pdf_idx} fraction {frac}")
        self.fitter.fix_bkg_frac_to_bkg_pdf(pdf_idx, anchor_pdf_idx, frac)

    def set_frac_to_sgn_pdf(self, pdf_idx, anchor_pdf_idx, frac):
        print(f"\nSetting background fraction of function index {pdf_idx} to signal PDF {anchor_pdf_idx} fraction {frac}")
        self.fitter.set_bkg_frac_to_signal_pdf(pdf_idx, anchor_pdf_idx, frac)

    def set_frac_to_bkg_pdf(self, pdf_idx, anchor_pdf_idx, frac):
        print(f"\nSetting background fraction of function index {pdf_idx} to background PDF {anchor_pdf_idx} fraction {frac}")
        self.fitter.set_bkg_frac_to_bkg_pdf(pdf_idx, anchor_pdf_idx, frac)

    def setup(self):
        if self.minimize_flarefly:
            self.setup_flarefly()
        else:
            self.setup_roofit()

    def fit(self):
        if self.minimize_flarefly:
            return self.perform_fit_flarefly()
        else:
            return self.perform_fit_roofit()

    def setup_flarefly(self):
        self.sgn_funcs_labels[-1] = self.sgn_main_label  # Set last signal function label as main particle label
        
        self.fitter = F2MassFitter(self.data,
                                   label_signal_pdf=self.sgn_funcs_labels,
                                   name_signal_pdf=self.sgn_funcs,
                                   name_background_pdf=self.corr_bkgs + self.bkg_funcs if self.corr_bkgs is not None else self.bkg_funcs,
                                   label_bkg_pdf=self.corr_bkgs_labels + self.bkg_funcs_labels if self.corr_bkgs_labels is not None else self.bkg_funcs_labels,
                                   name=self.fit_name)

        # Set particle mass and sigma initial par for the main signal
        # function here, so they can be overridden later if needed
        self.fitter.set_particle_mass(len(self.sgn_funcs)-1, pdg_id=self.particle_pdg)
        self.fitter.set_signal_initpar(len(self.sgn_funcs)-1, "sigma", 0.015, limits=[0.005, 0.05])
        n_bkg_funcs = len(self.corr_bkgs) + len(self.bkg_funcs) if self.corr_bkgs is not None else len(self.bkg_funcs)
        self.fitter.set_background_initpar(n_bkg_funcs-1, "c0", 1000.0, limits=[0.0, 1e6])         # Resonable value for c0
        self.fitter.set_background_initpar(n_bkg_funcs-1, "c1", 0.0, limits=[-1000.0, 1000.0])     # Resonable value for c1
        self.fitter.set_background_initpar(n_bkg_funcs-1, "c2", 0.0, limits=[-1000.0, 1000.0])     # Resonable value for c2

        # Initialize corr bkgs templs
        if self.corr_bkgs is not None:
            for idx, corr_bkg in enumerate(self.corr_bkgs):
                self.fitter.set_background_template(idx, self.corr_bkgs_histos[idx]) if self.data.get_is_binned() else \
                self.fitter.set_background_kde(idx, self.corr_bkgs_trees[idx])

    def perform_fit_flarefly(self):
        print(f"\nPerforming fit on data with fit range {self.fit_range_min} - {self.fit_range_max} GeV/c")
        self.fit_result = self.fitter.mass_zfit()
        return self.fit_result

    def plot_fit(self, logy, show_extra_info, loc=None):
        fig, axs = self.fitter.plot_mass_fit(style="ATLAS",
                                             figsize=(8, 8),
                                             axis_title=self.x_axis_label,
                                             show_extra_info=show_extra_info,
                                             logy=logy,
                                             extra_info_loc=loc if loc is not None else ["lower right", "lower left"]
                                            )
        return fig, axs

    def plot_raw_residuals(self):
        fig_res = self.fitter.plot_raw_residuals(style="ATLAS",
                                                 figsize=(8, 8),
                                                 axis_title=self.x_axis_label)
        return fig_res

    def plot_std_residuals(self):
        fig_pulls = self.fitter.plot_std_residuals(style="ATLAS",
                                                   figsize=(8, 8),
                                                   axis_title=self.x_axis_label)
        return fig_pulls
    
    def get_sweights_sgn(self, func_idx):
        print(f"\nGetting sWeights for signal function index {func_idx}")
        return self.fitter.get_sweights()[f"signal{func_idx}"]
    
    def get_sweights_bkg(self, func_idx):
        print(f"\nGetting sWeights for background function index {func_idx}")
        return self.fitter.get_sweights()[f"bkg{func_idx}"]
    
    def get_fitter(self):
        return self.fitter
    
    def get_sgn_labels(self):
        return self.sgn_funcs_labels
    
    def get_bkg_labels(self):
        return self.bkg_funcs_labels

    def get_fit_info(self):

        fit_info = {}
        sgn_func_idx = len(self.sgn_funcs) - 1
        try:
            fit_info['ry'] = self.fitter.get_raw_yield(sgn_func_idx)[0]
            fit_info['ry_unc'] = self.fitter.get_raw_yield(sgn_func_idx)[1]
        except Exception as e:
            logger(f"Could not get raw yield: {e}", "WARNING")
            fit_info['ry'] = -1.
            fit_info['ry_unc'] = -1.
        try:
            fit_info['ry_bin_counting'] = self.fitter.get_raw_yield_bincounting(sgn_func_idx, nsigma=5)[0]
            fit_info['ry_bin_counting_unc'] = self.fitter.get_raw_yield_bincounting(sgn_func_idx, nsigma=5)[1]
        except Exception as e:
            logger(f"Could not get raw yield from bin counting: {e}", "WARNING")
            fit_info['ry_bin_counting'] = -1.
            fit_info['ry_bin_counting_unc'] = -1.
        try:
            fit_info['signif'] = self.fitter.get_significance(sgn_func_idx)[0]
            fit_info['signif_unc'] = self.fitter.get_significance(sgn_func_idx)[1]
        except Exception as e:
            logger(f"Could not get significance: {e}", "WARNING")
            fit_info['signif'] = -1.
            fit_info['signif_unc'] = -1.
        try:
            fit_info['s_over_b'] = self.fitter.get_signal_over_background(sgn_func_idx)[0]
            fit_info['s_over_b_unc'] = self.fitter.get_signal_over_background(sgn_func_idx)[1]
        except Exception as e:
            logger(f"Could not get signal over background: {e}", "WARNING")
            fit_info['s_over_b'] = -1.
            fit_info['s_over_b_unc'] = -1.
        try:
            fit_info['chi2_fits'] = float(self.fitter.get_chi2())
            fit_info['chi2_over_ndf_fits'] = float(self.fitter.get_chi2())/self.fitter.get_ndf()
        except Exception as e:
            logger(f"Could not get chi2: {e}", "WARNING")
            fit_info['chi2_fits'] = -1.
            fit_info['chi2_over_ndf_fits'] = -1.

        signal_pars = self.fitter.get_signal_pars()
        signal_pars_uncs = self.fitter.get_signal_pars_uncs()
        bkg_pars = self.fitter.get_bkg_pars()
        bkg_pars_uncs = self.fitter.get_bkg_pars_uncs()

        return fit_info, signal_pars, signal_pars_uncs, bkg_pars, bkg_pars_uncs

    # def fit_histo_roofit(histo, bkgfuncs, sgnfuncs, self.fit_range_min, self.fit_range_max, outfilename, mean_int=None, sigma_int=None):
    def perform_fit_roofit():

        mass = ROOT.RooRealVar("mass", "Invariant Mass", self.fit_range_min, self.fit_range_max)
        mass.setRange("fitRange", self.fit_range_min, self.fit_range_max)
        vnvsmass = ROOT.RooRealVar("vnvsmass", "Vn Vs Mass", self.fit_range_min, self.fit_range_max)
        vnvsmass.setRange("fitRange", self.fit_range_min, self.fit_range_max)
        data_hist = ROOT.RooDataHist("data_hist", "Dataset from histogram", ROOT.RooArgList(mass), self.hist)

        # # Define Double-Sided Crystal Ball (DSCB) components for the signal
        # if mean_int is not None and sigma_int is not None:
        #     sgn_mean = ROOT.RooRealVar("sgn_mean", "Signal Mean", mean_int, 1.85, 1.88)
        #     sgn_width = ROOT.RooRealVar("sgn_width", "Signal Width", sigma_int, 0.01, 0.02)
        #     print(f"Using mean_int = {mean_int}, sigma_int = {sigma_int}")
        #     sgn_mean.setConstant(True)
        #     sgn_width.setConstant(True)
        # else:
        #     sgn_mean = ROOT.RooRealVar("sgn_mean", "Signal Mean", 1.86, 1.85, 1.88)
        #     sgn_width = ROOT.RooRealVar("sgn_width", "Signal Width", 0.01, 0.005, 0.02)
        #     print("Using default mean = 1.86, sigma = 0.01")

        gaussian = ROOT.RooGaussian("gaussian", "gaussian", mass, sgn_mean, sgn_width)

        # Define Background Polynomial component (pol is implememted as 1 + ...)
        c1_bkg = ROOT.RooRealVar("c1_bkg", "c1_bkg", -1, -3, 3)
        c2_bkg = ROOT.RooRealVar("c2_bkg", "c2_bkg", 0, -1, 1)
        c3_bkg = ROOT.RooRealVar("c3_bkg", "c3_bkg", 0, -1, 1)
        background = ROOT.RooPolynomial("background", "Polynomial Background", mass, ROOT.RooArgList(c1_bkg, c2_bkg, c3_bkg))

        # Create the RooAddPdf with components and their fractions
        n_signal = ROOT.RooRealVar("n_signal", "Number of signal events", 5000, 0, 1e6)
        n_background = ROOT.RooRealVar("n_background", "Number of background events", 10000, 0, 1e6)

        # Construct the extended model using RooAddPdf
        total_pdf_mass = ROOT.RooAddPdf("total_pdf_mass", "Signal + Background",
                                        ROOT.RooArgList(gaussian, background),
                                        ROOT.RooArgList(n_signal, n_background))

        # Perform the extended maximum likelihood fit
        fit_result = total_pdf_mass.fitTo(data_hist, ROOT.RooFit.Save(),
                                          ROOT.RooFit.Extended(True),
                                          ROOT.RooFit.Range("fitRange"),
                                          ROOT.RooFit.PrintLevel(-1),       # Suppresses messages
                                          ROOT.RooFit.Verbose(False),       # Turns off verbose mode
                                          ROOT.RooFit.Warnings(False),      # Suppresses warnings
                                          ROOT.RooFit.Timer(False)          # Disables timing info
                                        )
        fit_result.printMultiline(ROOT.std.cout, 3, True)

        # Create frame for plotting
        frame = mass.frame()
        data_hist.plotOn(frame, ROOT.RooFit.Range("fitRange"))  # Plot within fit range
        total_pdf_mass.plotOn(frame, ROOT.RooFit.Range("fitRange"))  # Plot the full model within fit range
        total_pdf_mass.plotOn(frame, ROOT.RooFit.Components("signal"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Range("fitRange"))  # Gaussian
        total_pdf_mass.plotOn(frame, ROOT.RooFit.Components("background"), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Range("fitRange"))  # Background

        # Draw the plot
        canvas = ROOT.TCanvas("canvas", "Fit Result", 800, 600)
        frame.Draw()
        os.makedirs(os.path.dirname(outfilename), exist_ok=True)
        canvas.SaveAs(outfilename)

        return (n_signal.getVal(), n_signal.getError()), (sgn_mean.getVal(), sgn_mean.getError()), (sgn_width.getVal(), sgn_width.getError())