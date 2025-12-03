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
from correlated_bkgs import get_corr_bkg
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
        logger("\nInitializing RawYieldFitter")
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
        self.templs = None
        self.sgn_templ_frac = None
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
        logger(f"\nSetting particle to fit: {particle_name}")
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
            logger(f"Particle {particle_name} not recognized!")
            sys.exit(1)

    def set_name(self, fit_name):
        logger(f"\nSetting fit name: {fit_name}")
        self.fit_name = fit_name

    def set_fit_range(self, fit_range_min, fit_range_max):
        logger(f"\nSetting fit range: {fit_range_min} - {fit_range_max} GeV/c")
        self.fit_range_min = fit_range_min
        self.fit_range_max = fit_range_max

    def set_rebin(self, rebin):
        logger(f"\nSetting rebin: {rebin}")
        self.rebin = rebin

    def set_data_to_fit_hist(self, data):
        logger(f"\nSetting data to fit from histogram with limits {self.fit_range_min} - {self.fit_range_max} GeV/c"
              f" for fitter with name {self.fit_name}")
        self.hist = data
        self.data = DataHandler(data, limits=[self.fit_range_min, self.fit_range_max], rebin=self.rebin)

    def set_data_to_fit_df(self, data_df, var_name='fM'):
        logger(f"\nSetting data to fit from dataframe with variable {var_name} and limits {self.fit_range_min} - {self.fit_range_max} GeV/c")
        self.data = DataHandler(data_df, var_name=var_name, limits=[self.fit_range_min, self.fit_range_max])

    def set_data_to_fit_tree(self, data_tree, var_name):
        logger(f"\nSetting data to fit from tree with variable {var_name} and limits {self.fit_range_min} - {self.fit_range_max} GeV/c")
        self.data = DataHandler(data_tree, var_name=var_name, limits=[self.fit_range_min, self.fit_range_max])

    def add_sgn_func(self, sgn_func, label):
        logger(f"\nAdding signal function with label {label}")
        if self.sgn_funcs is None:
            self.sgn_funcs = []
            self.sgn_funcs_labels = []
        if sgn_func.startswith("k"):
            if sgn_func == "kGaus":
                sgn_func = "gaussian"
            else:
                logger(f"Signal function {sgn_func} not recognized!")
                sys.exit(1)
        self.sgn_funcs.append(sgn_func)
        self.sgn_funcs_labels.append(label)

    def add_bkg_func(self, bkg_func, label):
        logger(f"\nAdding background function with label {label}, self.bkg_funcs: {self.bkg_funcs}, self.bkg_funcs_labels: {self.bkg_funcs_labels}")
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
                logger(f"Background function {bkg_func} not recognized!")
                sys.exit(1)
        self.bkg_funcs.append(bkg_func)
        self.bkg_funcs_labels.append(label)

    def add_corr_bkgs(self, cfg, sel_string, pt_min, pt_max, var_name='fM'):

        # find correlated bkg cocktail associated to this pt-bin
        cocktail_cfg = None
        for cocktail in cfg['cocktails']:
            for pt_range in cocktail['pt_ranges']:
                pt_center = (pt_range[0] + pt_range[1]) / 2
                if (pt_center >= pt_min and pt_center < pt_max):
                    cocktail_cfg = cocktail
                    break
            if cocktail_cfg is not None:
                break

        if cocktail_cfg is None:
            logger(f"No correlated background cocktail found for pt range {pt_min} - {pt_max} GeV/c", level="WARNING")
            return

        if self.templs is None:
            self.templs = {}

        pt_label = f"pt_{int(pt_min*10)}_{int(pt_max*10)}"
        corr_bkg_file = TFile.Open(f"{cfg['input_files']}_{pt_label}.root", "READ")

        logger(f"Adding correlated backgrounds to the fitter for pt range {pt_min} - {pt_max} GeV/c", level="INFO")
        sgn_fin_state = cfg['sgn_fin_state']
        count_sgn_templs, count_bkg_templs = 0, 0
        for chn in cocktail_cfg['channels']:
            name = chn['name']
            templ, frac = get_corr_bkg(corr_bkg_file, name, sel_string, pt_label, cfg['templ_type'], cfg['output_type'])
            if frac < 1e-10:
                logger(f"Skipping correlated bkg source {name} with negligible fraction {frac}", level="WARNING")
                continue

            self.templs[name] = {}
            self.templs[name]['frac'] = frac
            logger(f"\nSetting correlated bkg source {chn}")
            if chn.get('sgn_func'):
                self.templs[name]['type'] = 'sgn'
                # Allow for modeling signal with templates
                if chn['sgn_func'] == 'templ':
                    # set_data_to_fit_hist or set_data_to_fit_df must be called before to add_corr_bkgs
                    self.templs[name]['func'] = 'hist' if self.data.get_is_binned() else 'kde_grid'
                else:
                    self.templs[name]['func'] = chn['sgn_func']
                self.templs[name]['idx'] = count_sgn_templs
                count_sgn_templs += 1
            elif chn.get('bkg_func'):
                self.templs[name]['type'] = 'bkg'
                # set_data_to_fit_hist or set_data_to_fit_df must be called before to add_corr_bkgs
                self.templs[name]['func'] = 'hist' if self.data.get_is_binned() else 'kde_grid'
                self.templs[name]['idx'] = count_bkg_templs
                count_bkg_templs += 1
            else:
                logger(f"Correlated bkg source {chn} without 'sgn_func' or 'bkg_func' key", level="ERROR")
                continue

            logger(f"templ: {templ}, type: {type(templ)}")
            self.templs[name]['data'] = templ
            if chn.get('fix_to'):
                self.templs[name]['anchor_fix'] = chn.get('fix_to', None)
            elif chn.get('init_to'):
                self.templs[name]['anchor_init'] = chn.get('init_to', None)
            else:
                logger(f"Correlated bkg source {chn} without 'fix_to' or 'init_to' key", level="WARNING")

        _, self.sgn_templ_frac = get_corr_bkg(corr_bkg_file, cfg['sgn_fin_state'], sel_string, pt_label, cfg['templ_type'], cfg['output_type'])
        corr_bkg_file.Close()

    def fix_sgn_par(self, sgn_func_idx, par_name, par_val):
        logger(f"\nFixing signal parameter {par_name} of function index {sgn_func_idx} to value {par_val}")
        self.fitter.set_signal_initpar(sgn_func_idx, par_name, par_val, fix=True)

    def fix_bkg_par(self, bkg_func_idx, par_name, par_val):
        logger(f"\nFixing background parameter {par_name} of function index {bkg_func_idx} to value {par_val}")
        self.fitter.set_background_initpar(bkg_func_idx, par_name, par_val, fix=True)

    def set_sgn_par(self, sgn_func_idx, par_name, par_val, lims):
        logger(f"\nSetting signal parameter {par_name} of function index {sgn_func_idx} to value {par_val}")
        self.fitter.set_signal_initpar(sgn_func_idx, par_name, par_val, limits=lims)

    def set_bkg_par(self, bkg_func_idx, par_name, par_val, lims):
        logger(f"\nSetting background parameter {par_name} of function index {bkg_func_idx} to value {par_val}")
        self.fitter.set_background_initpar(bkg_func_idx, par_name, par_val, limits=lims)

    def set_pdf_frac(self, pdf_idx, frac, pdf_type):
        logger(f"\nSetting fraction of function index {pdf_idx} and type {pdf_type} to {frac}")
        if pdf_type == 'sgn':
            self.fitter.set_signal_initpar(pdf_idx, "frac", frac, limits=[0., 1.])
        else:
            self.fitter.set_background_initpar(pdf_idx, "frac", frac, limits=[0., 1.])

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

        templs_sgn, templs_bkg = {}, {}
        if self.templs is not None:
            for name, templ in list(self.templs.items())[::-1]:   # reverse loop over dict, else the indices get messed up
                if templ['type'] == 'sgn':
                    # Insert label at the beginning
                    self.sgn_funcs.insert(0, templ['func'])
                    self.sgn_funcs_labels.insert(0, name)
                else:
                    self.bkg_funcs.insert(0, templ['func'])
                    self.bkg_funcs_labels.insert(0, name)

        logger(f"\nSetting up FlareFly F2MassFitter with signal functions {self.sgn_funcs} and background functions {self.bkg_funcs}")
        logger(f"sgn_funcs_labels: {self.sgn_funcs_labels}, bkg_funcs_labels: {self.bkg_funcs_labels}")
        self.fitter = F2MassFitter(self.data, name=self.fit_name,
                                   label_signal_pdf=self.sgn_funcs_labels, name_signal_pdf=self.sgn_funcs,
                                   name_background_pdf=self.bkg_funcs, label_bkg_pdf=self.bkg_funcs_labels)
                                #    extended=self.data.get_is_binned())

        # Set particle mass and sigma initial par for the main signal
        # function here, so they can be overridden later if needed
        self.fitter.set_particle_mass(len(self.sgn_funcs)-1, pdg_id=self.particle_pdg)
        self.fitter.set_signal_initpar(len(self.sgn_funcs)-1, "sigma", 0.015, limits=[0.005, 0.05])
        self.fitter.set_background_initpar(len(self.bkg_funcs)-1, "c0", 1000.0, limits=[0.0, 1e6])         # Resonable value for c0
        self.fitter.set_background_initpar(len(self.bkg_funcs)-1, "c1", 0.0, limits=[-1000.0, 1000.0])     # Resonable value for c1
        self.fitter.set_background_initpar(len(self.bkg_funcs)-1, "c2", 0.0, limits=[-1000.0, 1000.0])     # Resonable value for c2

        if self.templs is None:
            return
        # Setup templates data handlers and fractions
        for name, templ in self.templs.items():
            logger(f"\nSetting up template for correlated source {name} with template {templ}")
            # Set data handler
            if self.data.get_is_binned():
                data_hdl = DataHandler(templ['data'], limits=(self.fit_range_min, self.fit_range_max), \
                                       rebin=self.rebin)
                logger(f"Setting binned template histogram for source {name}")
                if templ['type'] == 'sgn':
                    logger(f"Setting signal template for source {name}, idx {templ['idx']}")
                    self.fitter.set_signal_template(templ['idx'], data_hdl)
                else:
                    logger(f"Setting background template for source {name}, idx {templ['idx']}")
                    self.fitter.set_background_template(templ['idx'], data_hdl)
            else:
                logger(f"Setting unbinned template KDE for source {name}")
                data_hdl = DataHandler(templ['data'], limits=(self.fit_range_min, self.fit_range_max), \
                                       nbins=100, var_name="fM")
                logger(f"data_hdl: {data_hdl}")
                if templ['type'] == 'sgn':
                    logger(f"Setting signal KDE for source {name}, idx {templ['idx']}")
                    self.fitter.set_signal_kde(templ['idx'], data_hdl)
                else:
                    logger(f"Setting background KDE for source {name}, idx {templ['idx']}")
                    self.fitter.set_background_kde(templ['idx'], data_hdl)
            # Set fraction
            if templ.get('anchor_fix'):
                anchor_mode = 'anchor_fix'
            elif templ.get('anchor_init'):
                anchor_mode = 'anchor_init'
            else:
                logger(f"Correlated bkg source {name} without 'fix_to' or 'init_to' key", level="WARNING")
                continue
            anchor_func = templ[anchor_mode]
            if anchor_func == "signal":
                anchor_pdf_idx = len(self.sgn_funcs)-1
                anchor_pdf_frac = self.sgn_templ_frac
                frac = templ['frac'] / anchor_pdf_frac
                logger(f"frac of chn {name} wrt anchor {templ[anchor_mode]}: {templ['frac']} / {anchor_pdf_frac} = {frac}", "INFO")
                if anchor_mode == 'anchor_fix':
                    self.fitter.fix_bkg_frac_to_signal_pdf(templ['idx'], anchor_pdf_idx, frac)
                else:
                    self.set_pdf_frac(templ['idx'], frac, templ['type'])
                continue

            anchor_pdf_idx = self.templs[anchor_func]['idx']
            anchor_pdf_frac = self.templs[anchor_func]['frac']
            frac = templ['frac'] / anchor_pdf_frac
            logger(f"frac of chn {name} wrt anchor {templ[anchor_mode]}: {templ['frac']} / 67{anchor_pdf_frac} = {frac}", "INFO")

            if anchor_mode == 'anchor_fix':
                if templs[anchor_func]['type'] == 'sgn' and templ['type'] == 'bkg':
                    self.fitter.fix_bkg_frac_to_signal_pdf(templ['idx'], anchor_pdf_idx, frac)
                elif templs[anchor_func]['type'] == 'bkg' and templ['type'] == 'bkg':
                    self.fitter.fix_bkg_frac_to_bkg_pdf(templ['idx'], anchor_pdf_idx, frac)
                elif templs[anchor_func]['type'] == 'sgn' and templ['type'] == 'sgn':
                    self.fitter.fix_signal_frac_to_signal_pdf(templ['idx'], anchor_pdf_idx, frac)
                else:
                    self.fitter.fix_bkg_frac_to_bkg_pdf(templ['idx'], anchor_pdf_idx, frac)
            else:
                self.set_pdf_frac(templ['idx'], frac, templ['type'])

    def perform_fit_flarefly(self):
        logger(f"\nPerforming fit on data with fit range {self.fit_range_min} - {self.fit_range_max} GeV/c")
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
        # Count corr bkgs before main signal
        for name, templ in self.templs.items():
            if templ['type'] == 'sgn':
                func_idx += 1
        logger(f"\nGetting sWeights for signal function index {func_idx}")
        return self.fitter.get_sweights()[f"signal{func_idx}"]

    def get_sweights_bkg(self, func_idx):
        # Count corr bkgs before main signal
        for name, templ in self.templs.items():
            if templ['type'] == 'bkg':
                func_idx += 1
        logger(f"\nGetting sWeights for background function index {func_idx}")
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

    # # def fit_histo_roofit(histo, bkgfuncs, sgnfuncs, self.fit_range_min, self.fit_range_max, outfilename, mean_int=None, sigma_int=None):
    # def perform_fit_roofit():

    #     mass = ROOT.RooRealVar("mass", "Invariant Mass", self.fit_range_min, self.fit_range_max)
    #     mass.setRange("fitRange", self.fit_range_min, self.fit_range_max)
    #     vnvsmass = ROOT.RooRealVar("vnvsmass", "Vn Vs Mass", self.fit_range_min, self.fit_range_max)
    #     vnvsmass.setRange("fitRange", self.fit_range_min, self.fit_range_max)
    #     data_hist = ROOT.RooDataHist("data_hist", "Dataset from histogram", ROOT.RooArgList(mass), self.hist)

    #     # # Define Double-Sided Crystal Ball (DSCB) components for the signal
    #     # if mean_int is not None and sigma_int is not None:
    #     #     sgn_mean = ROOT.RooRealVar("sgn_mean", "Signal Mean", mean_int, 1.85, 1.88)
    #     #     sgn_width = ROOT.RooRealVar("sgn_width", "Signal Width", sigma_int, 0.01, 0.02)
    #     #     logger(f"Using mean_int = {mean_int}, sigma_int = {sigma_int}")
    #     #     sgn_mean.setConstant(True)
    #     #     sgn_width.setConstant(True)
    #     # else:
    #     #     sgn_mean = ROOT.RooRealVar("sgn_mean", "Signal Mean", 1.86, 1.85, 1.88)
    #     #     sgn_width = ROOT.RooRealVar("sgn_width", "Signal Width", 0.01, 0.005, 0.02)
    #     #     logger("Using default mean = 1.86, sigma = 0.01")

    #     gaussian = ROOT.RooGaussian("gaussian", "gaussian", mass, sgn_mean, sgn_width)

    #     # Define Background Polynomial component (pol is implememted as 1 + ...)
    #     c1_bkg = ROOT.RooRealVar("c1_bkg", "c1_bkg", -1, -3, 3)
    #     c2_bkg = ROOT.RooRealVar("c2_bkg", "c2_bkg", 0, -1, 1)
    #     c3_bkg = ROOT.RooRealVar("c3_bkg", "c3_bkg", 0, -1, 1)
    #     background = ROOT.RooPolynomial("background", "Polynomial Background", mass, ROOT.RooArgList(c1_bkg, c2_bkg, c3_bkg))

    #     # Create the RooAddPdf with components and their fractions
    #     n_signal = ROOT.RooRealVar("n_signal", "Number of signal events", 5000, 0, 1e6)
    #     n_background = ROOT.RooRealVar("n_background", "Number of background events", 10000, 0, 1e6)

    #     # Construct the extended model using RooAddPdf
    #     total_pdf_mass = ROOT.RooAddPdf("total_pdf_mass", "Signal + Background",
    #                                     ROOT.RooArgList(gaussian, background),
    #                                     ROOT.RooArgList(n_signal, n_background))

    #     # Perform the extended maximum likelihood fit
    #     fit_result = total_pdf_mass.fitTo(data_hist, ROOT.RooFit.Save(),
    #                                       ROOT.RooFit.Extended(True),
    #                                       ROOT.RooFit.Range("fitRange"),
    #                                       ROOT.RooFit.PrintLevel(-1),       # Suppresses messages
    #                                       ROOT.RooFit.Verbose(False),       # Turns off verbose mode
    #                                       ROOT.RooFit.Warnings(False),      # Suppresses warnings
    #                                       ROOT.RooFit.Timer(False)          # Disables timing info
    #                                     )
    #     fit_result.loggerMultiline(ROOT.std.cout, 3, True)

    #     # Create frame for plotting
    #     frame = mass.frame()
    #     data_hist.plotOn(frame, ROOT.RooFit.Range("fitRange"))  # Plot within fit range
    #     total_pdf_mass.plotOn(frame, ROOT.RooFit.Range("fitRange"))  # Plot the full model within fit range
    #     total_pdf_mass.plotOn(frame, ROOT.RooFit.Components("signal"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Range("fitRange"))  # Gaussian
    #     total_pdf_mass.plotOn(frame, ROOT.RooFit.Components("background"), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Range("fitRange"))  # Background

    #     # Draw the plot
    #     canvas = ROOT.TCanvas("canvas", "Fit Result", 800, 600)
    #     frame.Draw()
    #     os.makedirs(os.path.dirname(outfilename), exist_ok=True)
    #     canvas.SaveAs(outfilename)

    #     return (n_signal.getVal(), n_signal.getError()), (sgn_mean.getVal(), sgn_mean.getError()), (sgn_width.getVal(), sgn_width.getError())