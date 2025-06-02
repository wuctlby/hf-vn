'''
Script for fitting D+, D0 and Ds+ invariant-mass spectra
run: python get_vn_vs_mass.py fitConfigFileName.yml centClass inputFileName.root outFileName.root
            [--refFileName][--isMC][--batch]
'''
import re
import sys
import argparse
import ctypes
import numpy as np
import yaml
from ROOT import TLatex, TFile, TCanvas, TLegend, TH1D, TH1F, TDatabasePDG, TGraphAsymmErrors # pylint: disable=import-error,no-name-in-module
from ROOT import gROOT, gPad, gInterpreter, kBlack, kRed, kAzure, kCyan, kBlue, kGray, kOrange, kGreen, kMagenta, kFullCircle, kFullSquare, kOpenCircle # pylint: disable=import-error,no-name-in-module
import os
work_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(work_dir, '../')))
gInterpreter.ProcessLine(f'#include "{work_dir}/invmassfitter/InvMassFitter.cxx"')
gInterpreter.ProcessLine(f'#include "{work_dir}/invmassfitter/VnVsMassFitter.cxx"')
from ROOT import InvMassFitter, VnVsMassFitter
from utils.flow_analysis_utils import get_centrality_bins, get_vnfitter_results, get_refl_histo, get_particle_info # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.histo_operator import rebin_histo
from utils.check import check_ptbinned_pars

def get_vn_vs_mass(flow_config, centClass, inFileName,
                   outputdir, suffix, batch):


    with open(flow_config, 'r') as f:
        cfg = yaml.safe_load(f)

    # extract the cutset number from the suffix
    cut_var_suffix = re.search(rf"_(\d{2})", suffix)
    cut_var_suffix = cut_var_suffix.group(1) if cut_var_suffix else None

    gROOT.SetBatch(batch)
    SetGlobalStyle(padleftmargin=0.14, padbottommargin=0.12, padtopmargin=0.12, opttitle=1)
    cent, centMinMax = get_centrality_bins(centClass)

    # fixed parameters
    harmonic = 2
    vn_method = 'sp'

    # read global configuration______________________________________________________________________
    ptMins = cfg['ptbins'][:-1]
    ptMaxs = cfg['ptbins'][1:]
    ptLims = list(ptMins)
    nPtBins = len(ptMins)
    ptLims.append(ptMaxs[-1])
    ptBinsArr = np.asarray(ptLims, 'd')
    ptTit = '#it{p}_{T} (GeV/#it{c})'
    fixSigma, fixMean, rebins = check_ptbinned_pars(nPtBins,
        cfg['FixSigma'], cfg['FixMean'], cfg['Rebin'])
    sigma = check_ptbinned_pars(nPtBins, cfg.get('Sigma', 0.01))[0]
    fixSigmaFromFile = cfg.get('FixSigmaFromFile', [])
    massMins = [range[0] for range in cfg['MassFitRanges']]
    massMaxs = [range[1] for range in cfg['MassFitRanges']]
    particleName = cfg['Dmeson']
    
    # optional parameters
    useRefl = cfg.get('UseRefl', False)
    reflFile = cfg.get('ReflFile', '')
    reflFuncStr = cfg.get('ReflFunc', '2Gaus')

    inclSecPeak = check_ptbinned_pars(nPtBins, cfg.get('IncludeSecondPeak', 0))[0]
    if particleName == 'Dplus':
        fixVnSecPeakToSgn = cfg.get('FixVnSecPeakToSgn', True)
    elif particleName == 'Ds':
        fixVnSecPeakToSgn = cfg.get('FixVnSecPeakToSgn', False)
    sigmaRatioFile = cfg.get('SigmaRatioFile', '')
    sigmaSecPeak = cfg.get('SigmaSecPeak')
    
    drawVnComps = cfg.get('DrawVnComps', False)

    # corr bkg
    prefitMC = cfg.get("PrefitMC", False)
    includeTempls = cfg.get('IncludeTempls', False)
    templsNames = cfg.get('TemplsNames', [])
    if includeTempls:
        anchor_templs_mode = cfg['AnchorTemplsMode']
        fix_vn_templ_to_sgn = cfg.get('FixVnTemplToSgn', False)
        init_weights = cfg['InitWeights']
        min_weights = cfg['MinWeights']
        max_weights = cfg['MaxWeights']
        vn_init_weights = cfg['VnInitWeights']
        vn_min_weights = cfg['VnMinWeights']
        vn_max_weights = cfg['VnMaxWeights']
    initFitPars = cfg.get('InitFitPars', [])
    drawSingleTempls = cfg.get('DrawSingleTempls', False)
    drawDynRange = cfg.get('DrawDynRange', True)

    # read fit configuration
    SgnFuncStr, BkgFuncStr, BkgFuncVnStr = check_ptbinned_pars(nPtBins,
        cfg['SgnFunc'], cfg['BkgFunc'], cfg['BkgFuncVn'])
    #________________________________________________________________________________________________

    # sanity check of fit configuration
    SgnFunc, BkgFunc, BkgFuncVn, degPol = [], [], [], []
    for iPt, (bkgStr, sgnStr, bkgVnStr) in enumerate(zip(BkgFuncStr, SgnFuncStr, BkgFuncVnStr)):
        degPol.append(-1)
        if bkgStr == 'kExpo':
            BkgFunc.append(InvMassFitter.kExpo)
        elif bkgStr == 'kLin':
            BkgFunc.append(InvMassFitter.kLin)
        elif bkgStr == 'kPol2':
            BkgFunc.append(InvMassFitter.kPol2)
        elif bkgStr == 'kPol3':
            BkgFunc.append(6)
            degPol[-1] = 3
        elif bkgStr == 'kPol4':
            BkgFunc.append(6)
            degPol[-1] = 4
            if nPtBins > 1 and inclSecPeak[iPt] == 1:
                print('ERROR: Pol3 and Pol4 fits work only with one bin if you have the secondary peak! Exit!')
                sys.exit()
        elif bkgStr == 'kPow':
            BkgFunc.append(InvMassFitter.kPow)
        elif bkgStr == 'kPowEx':
            BkgFunc.append(InvMassFitter.kPowEx)
        else:
            print('ERROR: only kExpo, kLin, kPol2, kPol3, kPol4, kPow, and kPowEx background functions supported! Exit')
            sys.exit()
        if bkgVnStr == 'kExpo':
            BkgFuncVn.append(InvMassFitter.kExpo)
        elif bkgVnStr == 'kLin':
            BkgFuncVn.append(InvMassFitter.kLin)
        elif bkgVnStr == 'kPol2':
            BkgFuncVn.append(InvMassFitter.kPol2)
        else:
            print('ERROR: only kExpo, kLin, and kPol2 background functions supported for vn! Exit')
            sys.exit()
        if sgnStr == 'kGaus':
            SgnFunc.append(InvMassFitter.kGaus)
        elif sgnStr == 'k2Gaus':
            SgnFunc.append(InvMassFitter.k2Gaus)
        elif sgnStr == 'kDoubleCBAsymm':
            SgnFunc.append(InvMassFitter.kDoubleCBAsymm)
        elif sgnStr == 'kDoubleCBSymm':
            SgnFunc.append(InvMassFitter.kDoubleCBSymm)
        elif sgnStr == 'k2GausSigmaRatioPar':
            SgnFunc.append(InvMassFitter.k2GausSigmaRatioPar)
        else:
            print('ERROR: only kGaus, k2Gaus, kDoubleCBAsymm, kDoubleCBSymm and k2GausSigmaRatioPar signal functions supported! Exit!')
            sys.exit()

    # Retrieve histogram to fix signal
    if prefitMC:
        mcFitFile = TFile.Open(f"{outputdir}/Prefit_mc_prompt_enhanced.root", 'r')
        mcFitFile.cd(f"{sgnStr}/")
        print(f"mcFitFile: {mcFitFile.GetName()}")
        print(f"sgnStr: {sgnStr}")
        directory = mcFitFile.GetDirectory(f"{sgnStr}/")

        # List all histogram names in the directory
        histosPars = []
        for key in directory.GetListOfKeys():
            obj = key.ReadObj()
            if obj.InheritsFrom("TH1") and not obj.GetName() == "hChi2":
                histosPars.append(obj)

        print("Histograms found:", histosPars)

    # set particle configuration
    if particleName == 'Dzero':
        _, massAxisTit, decay, massForFit = get_particle_info(particleName)
    if particleName == 'Ds':
        _, massAxisTit, decay, massForFit = get_particle_info(particleName)
        massDplus = TDatabasePDG.Instance().GetParticle(411).Mass()
    if particleName == 'Dplus':
        _, massAxisTit, decay, massForFit = get_particle_info(particleName)
        massDstar = TDatabasePDG.Instance().GetParticle(413).Mass()
    else:
        _, massAxisTit, decay, massForFit = get_particle_info(particleName)

    # load histos
    infile = TFile.Open(inFileName)
    if not infile or not infile.IsOpen():
        print(f'ERROR: file "{inFileName}" cannot be opened! Exit!')
        sys.exit()
    hMass, hMassForFit, hVn, hVnForFit, hPulls, hPullsPrefit, hRel = [], [], [], [], [], [], []
    fTotFuncMass, fTotFuncVn, fSgnFuncMass, fBkgFuncMass, fMassBkgRflFunc, fMassSecPeakFunc, fBkgFuncVn, fVnSecPeakFunc, fVnCompFuncts = [], [], [], [], [], [], [], [], []
    hMCSgn, hMCRefl = [], []
    
    fMassTemplFuncts = [[None]*len(templsNames) for _ in range(nPtBins)] if includeTempls and (particleName == 'Dplus' or particleName == 'Ds') else []
    fMassTemplTotFuncts = [None]*nPtBins if includeTempls and (particleName == 'Dplus' or particleName == 'Ds') else []
    fVnCompFuncts = [] # [[None]*len(fitConfig['TemplsNames']) for _ in range(nPtBins)] if fitConfig.get('DrawVnComps') else []

    hist_reso = infile.Get('hist_reso')
    hist_reso.SetDirectory(0)
    reso = hist_reso.GetBinContent(1)
    for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
        print(f'loading: cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_vn_{vn_method}_pt{ptMin}_{ptMax}')
        hMass.append(infile.Get(f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_mass_cent{cent}_pt{ptMin}_{ptMax}'))
        hVn.append(infile.Get(f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_vn_{vn_method}_pt{ptMin}_{ptMax}'))
        hVn[iPt].SetDirectory(0)
        hMass[iPt].SetDirectory(0)
        SetObjectStyle(hMass[iPt], color=kBlack, markerstyle=kFullCircle)
        SetObjectStyle(hVn[iPt], color=kBlack, markerstyle=kFullCircle)
    infile.Close()

    hSigmaToFix = None
    if sigmaRatioFile:
        # load sigma of first gaussian
        infileSigma = TFile.Open(sigmaRatioFile)
        if not infileSigma:
            print(f'ERROR: file "{infileSigma}" cannot be opened! Exit!')
            sys.exit()
        hSigmaToFix = infileSigma.Get('hRawYieldsSigma')
        hSigmaToFix.SetDirectory(0)
        if hSigmaToFix.GetNbinsX() != nPtBins:
            print('WARNING: Different number of bins for this analysis and histo for fix sigma')
        infileSigma.Close()
        # load sigma of second gaussian
        infileSigma2 = TFile.Open(sigmaRatioFile)
        if not infileSigma2:
            print(f'ERROR: file "{infileSigma2}" cannot be opened! Exit!')
            sys.exit()
        hSigmaToFix2 = infileSigma2.Get('hRawYieldsSigma2')
        hSigmaToFix2.SetDirectory(0)
        if hSigmaToFix2.GetNbinsX() != nPtBins:
            print('WARNING: Different number of bins for this analysis and histo for fix sigma')
        infileSigma2.Close()

    # check reflections
    if useRefl and particleName == 'Dzero':
        if reflFile == '':
            reflFile = inFileName
            useRefl, hMCSgn, hMCRefl = get_refl_histo(reflFile, centMinMax, ptMins, ptMaxs)
        else:
            useRefl, hMCSgn, hMCRefl = get_refl_histo(reflFile, centMinMax, ptMins, ptMaxs)
    else:
        useRefl = False

    # create histos for fit results
    hSigmaSimFit = TH1D('hSigmaSimFit', f';{ptTit};#sigma', nPtBins, ptBinsArr)
    hMeanSimFit = TH1D('hMeanSimFit', f';{ptTit};mean', nPtBins, ptBinsArr)
    hMeanSecPeakFitMass = TH1D('hMeanSecondPeakFitMass', f';{ptTit};mean second peak mass fit', nPtBins, ptBinsArr)
    hMeanSecPeakFitVn = TH1D('hMeanSecondPeakFitVn', f';{ptTit};mean second peak vn fit', nPtBins, ptBinsArr)
    hTemplOverSgn = TH1D('hTemplOverSgn', f';{ptTit};Templ / Sgn', nPtBins, ptBinsArr)
    hSigmaSecPeakFitMass = TH1D('hSigmaSecondPeakFitMass',
                                f';{ptTit};width second peak mass fit', nPtBins, ptBinsArr)
    hSigmaSecPeakFitVn = TH1D('hSigmaSecondPeakFitVn', f';{ptTit};width second peak vn fit', nPtBins, ptBinsArr)
    hRawYieldsSimFit = TH1D('hRawYieldsSimFit', f';{ptTit};raw yield', nPtBins, ptBinsArr)
    hRawYieldsTrueSimFit = TH1D('hRawYieldsTrueSimFit', f';{ptTit};raw yield true', nPtBins, ptBinsArr)
    hRawYieldsSecPeakSimFit = TH1D('hRawYieldsSecondPeakSimFit',
                                    f';{ptTit};raw yield second peak', nPtBins, ptBinsArr)
    hRawYieldsSignificanceSimFit = TH1D('hRawYieldsSignificanceSimFit',
                                        f';{ptTit};significance', nPtBins, ptBinsArr)
    hRawYieldsSoverBSimFit = TH1D('hRawYieldsSoverBSimFit', f';{ptTit};S/B', nPtBins, ptBinsArr)
    hRedChi2SimFit = TH1D('hRedChi2SimFit', f';{ptTit};#chi^{{2}}/#it{{ndf}}', nPtBins, ptBinsArr)
    hProbSimFit = TH1D('hProbSimFit', f';{ptTit};prob', nPtBins, ptBinsArr)
    hRedChi2SBVnPrefit = TH1D('hRedChi2SBVnPrefit', f';{ptTit};#chi^{{2}}/#it{{ndf}}', nPtBins, ptBinsArr)
    hProbSBVnPrefit = TH1D('hProbSBVnPrefit', f';{ptTit};prob', nPtBins, ptBinsArr)
    hvnSimFit = TH1D('hvnSimFit',f';{ptTit};V2 ({vn_method})', nPtBins, ptBinsArr)

    SetObjectStyle(hSigmaSimFit, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hMeanSimFit, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hMeanSecPeakFitMass, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hSigmaSecPeakFitMass, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hMeanSecPeakFitVn, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hSigmaSecPeakFitVn, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hRawYieldsSimFit, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hRawYieldsTrueSimFit, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hRawYieldsSecPeakSimFit, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hRawYieldsSignificanceSimFit, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hRawYieldsSoverBSimFit, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hRedChi2SimFit, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hProbSimFit, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hRedChi2SBVnPrefit, color=kRed, markerstyle=kFullSquare)
    SetObjectStyle(hProbSBVnPrefit, color=kRed, markerstyle=kFullSquare)
    SetObjectStyle(hvnSimFit, color=kBlack, markerstyle=kFullCircle)

    gvnSimFit = TGraphAsymmErrors(1)
    gvnSimFit.SetName('gvnSimFit')
    gvnSimFitSecPeak = TGraphAsymmErrors(1)
    gvnSimFitSecPeak.SetName('gvnSimFitSecPeak')
    gvnUnc = TGraphAsymmErrors(1)
    gvnUnc.SetName('gvnUnc')
    gvnUncSecPeak = TGraphAsymmErrors(1)
    gvnUncSecPeak.SetName('gvnUncSecPeak')
    gvnTempls = []
    gvnTemplsUncs = []
    if includeTempls:
        for iTempl in templsNames:
            gvnTempl = TGraphAsymmErrors(1)
            gvnTempl.SetName('gvnTemplUnc')
            gvnTempls.append(gvnTempl)
            gvnTemplUnc = TGraphAsymmErrors(1)
            gvnTemplUnc.SetName('gvnUncSecPeak')
            gvnTemplsUncs.append(gvnTemplUnc)
    SetObjectStyle(gvnSimFit, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(gvnSimFitSecPeak, color=kRed, markerstyle=kOpenCircle)
    SetObjectStyle(gvnUnc, color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(gvnUncSecPeak, color=kBlack, markerstyle=kOpenCircle)

    # create canvases
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)

    cSimFit = []
    # TODO: add residual plots, cResiduals
    cInvMassPrefits = []
    for i in range(nPtBins):
        ptLow = ptMins[i]
        ptHigh = ptMaxs[i]
        cSimFit.append(TCanvas(f'cSimFit_Pt{ptLow}_{ptHigh}', f'cSimFit_Pt{ptLow}_{ptHigh}', 400, 900))
        cInvMassPrefits.append(TCanvas(f'cMassPrefit_Pt{ptLow}_{ptHigh}', f'cMassPrefit_Pt{ptLow}_{ptHigh}', 400, 900))
        cSimFit[-1].Divide(1, 2)
    canvVn = TCanvas('cVn', 'cVn', 900, 900)
    canvVnUnc = TCanvas('canvVnUnc', 'canvVnUnc', 900, 900)

    #_____________________________________________________
    # Vn estimation with Scalar Product      
    vnFitter = []
    for iPt, (hM, hV, ptMin, ptMax, reb, sgnEnum, bkgEnum, bkgVnEnum, secPeak, massMin, massMax) in enumerate(
            zip(hMass, hVn, ptMins, ptMaxs, rebins, SgnFunc, BkgFunc, BkgFuncVn, inclSecPeak, massMins, massMaxs)):
        iCanv = iPt
        hMassForFit.append(TH1F())
        hVnForFit.append(TH1F())
        rebin_histo(hM, reb).Copy(hMassForFit[iPt]) #to cast TH1D to TH1F
        hMassForFit[iPt].SetDirectory(0)
        xbins = np.asarray(hV.GetXaxis().GetXbins())
        hDummy = TH1F('hDummy', '', len(xbins)-1, xbins)
        for iBin in range(1, hV.GetNbinsX()+1):
            hDummy.SetBinContent(iBin, hV.GetBinContent(iBin))
            hDummy.SetBinError(iBin, hV.GetBinError(iBin))
        hVnForFit[iPt] = hDummy
        hVnForFit[iPt].SetDirectory(0)
        hVnForFit[iPt].GetXaxis().SetTitle(massAxisTit)
        hVnForFit[iPt].GetYaxis().SetTitle(f'#it{{v}}{harmonic}')
        binWidth = hMassForFit[iPt].GetBinWidth(1)
        if cut_var_suffix is not None:
            hMassForFit[iPt].SetTitle((f'{ptMin:0.1f} < #it{{p}}_{{T}} < {ptMax:0.1f} GeV/#it{{c}}, cutset {cut_var_suffix};{massAxisTit};'
                                    f'Counts per {binWidth*1000:.0f} MeV/#it{{c}}^{{2}}'))
        else:
            hMassForFit[iPt].SetTitle((f'{ptMin:0.1f} < #it{{p}}_{{T}} < {ptMax:0.1f} GeV/#it{{c}};{massAxisTit};'
                                    f'Counts per {binWidth*1000:.0f} MeV/#it{{c}}^{{2}}'))
        hMassForFit[iPt].SetName(f'MassForFit{iPt}')
        SetObjectStyle(hMassForFit[iPt], color=kBlack, markerstyle=kFullCircle, markersize=1)
        SetObjectStyle(hVnForFit[iPt], color=kBlack, markerstyle=kFullCircle, markersize=0.8)

        print(f'Fitting {ptMin} - {ptMax} GeV/c')
        vnFitter.append(VnVsMassFitter(hMassForFit[iPt], hVnForFit[iPt],
                                            massMin, massMax, bkgEnum, sgnEnum, bkgVnEnum))
        vnFitter[iPt].SetHarmonic(harmonic)

        #_____________________________________________________
        # set the parameters for the fit
        # Mean
        vnFitter[iPt].SetInitialGaussianMean(massForFit, 1)
        if fixMean[iPt]:
            vnFitter[iPt].FixMeanFromMassFit()
        # Sigma
        if fixSigma[iPt]:
            if fixSigmaFromFile != []:
                sigmaFile = TFile.Open(fixSigmaFromFile[0], 'r')
                # get the sigma histo from config file
                hSigmaFromFile = sigmaFile.Get(fixSigmaFromFile[1])
                hSigmaFromFile.SetDirectory(0)
                sigmaBin = hSigmaFromFile.FindBin((ptMin+ptMax)/2)
                if hSigmaFromFile.GetBinLowEdge(sigmaBin) != ptMin:
                    bin_low_edge = hSigmaFromFile.GetBinLowEdge(sigmaBin)
                    bin_up_edge = bin_low_edge + hSigmaFromFile.GetBinWidth(sigmaBin)
                    if not (np.isclose(bin_low_edge, ptMin) and np.isclose(bin_up_edge, ptMax)):
                        print(f'ERROR: bin edges do not match! Expected [{ptMin}, {ptMax}], got [{bin_low_edge}, {bin_up_edge}]. Cannot load sigma from {fixSigmaFromFile[0]}! Exit!')
                        sys.exit()
                print(f'Fixing sigma from file {fixSigmaFromFile}: {hSigmaFromFile.GetBinContent(sigmaBin)}')
                vnFitter[iPt].SetInitialGaussianSigma(hSigmaFromFile.GetBinContent(sigmaBin), 2)
            else:
                vnFitter[iPt].SetInitialGaussianSigma(sigma[iPt], 2)
        else:
            vnFitter[iPt].SetInitialGaussianSigma(sigma[iPt], 1)
        # nSigma4SB
        if 'NSigma4SB' in cfg:
            nSigma4SB = check_ptbinned_pars(nPtBins, cfg['NSigma4SB'])[0]
            print(f'NSigma4SB = {nSigma4SB[iPt]}')
            vnFitter[iPt].SetNSigmaForVnSB(nSigma4SB[iPt])
        # Second peak (Ds specific)
        if secPeak and particleName == 'Ds':
            vnFitter[iPt].IncludeSecondGausPeak(massDplus, False, sigmaSecPeak[iPt], False, 1, fixVnSecPeakToSgn)
            if fixSigma[iPt]:
                vnFitter[iPt].SetInitialGaussianSigma2Gaus(sigmaSecPeak[iPt], 2)
        # Second peak (Dplus specific)
        if secPeak and particleName == 'Dplus':
            vnFitter[iPt].IncludeSecondGausPeak(massDstar, True, sigmaSecPeak[iPt], True, 1, fixVnSecPeakToSgn)
            if fixSigma[iPt]:
                vnFitter[iPt].SetInitialGaussianSigma2Gaus(sigmaSecPeak[iPt], 2)

        # Reflections for D0
        if useRefl:
            SoverR = (hMCRefl[iPt].Integral(hMCRefl[iPt].FindBin(massMin*1.0001),hMCRefl[iPt].FindBin(massMax*0.9999)))/(
                hMCSgn[iPt].Integral(hMCSgn[iPt].FindBin(massMin*1.0001),hMCSgn[iPt].FindBin(massMax*0.9999)))
            vnFitter[iPt].SetTemplateReflections(hMCRefl[iPt],reflFuncStr,massMin,massMax)
            vnFitter[iPt].SetFixReflOverS(SoverR)
            vnFitter[iPt].SetReflVnOption(0) # kSameVnSignal

        useTemplates = False
        if includeTempls:
            useTemplates = True
            weightsFile = TFile.Open(f"{cfg['out_dir']}/cutvar_{cfg['suffix']}/ry/weights.root", 'r')
            Templs, TemplsRelWeights = [], []
            for name in templsNames:
                Templs.append(weightsFile.Get(f"cutset_{cut_var_suffix}/{name}/hMassRaw"))
                if anchor_templs_mode == 2:    # Templates anchored to signal
                    TemplsRelWeights.append(weightsFile.Get(f"cutset_{cut_var_suffix}/{name}/hRelWeightToSgn").GetBinContent(1))
                else:
                    TemplsRelWeights.append(1.)
                    print("Relative template weighting still to be implemented!")
            print(f"Templs: {Templs}")
            vnFitter[iPt].SetTemplatesHisto(anchor_templs_mode, TemplsRelWeights, templsNames, Templs,
                                            init_weights[iPt] if anchor_templs_mode != 2 else [],
                                            min_weights[iPt] if anchor_templs_mode != 2 else [],
                                            max_weights[iPt] if anchor_templs_mode != 2 else [],
                                            vn_init_weights[iPt] if not fix_vn_templ_to_sgn else [],
                                            vn_min_weights[iPt] if not fix_vn_templ_to_sgn else [],
                                            vn_max_weights[iPt] if not fix_vn_templ_to_sgn else [],
                                            fix_vn_templ_to_sgn)
            print("Histo templates set!")

        initPars = []
        if initFitPars and initFitPars[iPt] is not None:
            initPars = initFitPars[iPt]
        if prefitMC:
            print(f"Fixing signal parameters from MC")
            sgnParsFromHisto = [[histoPar.GetName(), 
                                    histoPar.GetBinContent(iPt+1), # getBinContent is 1-based
                                    0 if histoPar.GetName() not in ["Mean", "Sigma"] else histoPar.GetBinContent([iPt]+1)-10, 
                                    -1 if histoPar.GetName() not in ["Mean", "Sigma"] else histoPar.GetBinContent([iPt]+1)+10] 
                                for histoPar in histosPars] 
            print(f"sgnParsFromHisto: {sgnParsFromHisto}")
            initPars.extend(sgnParsFromHisto)
        if len(initPars) > 0:
            vnFitter[iPt].SetInitPars(initPars)
            
        # collect fit results
        vnFitter[iPt].SimultaneousFit()
        vnResults = get_vnfitter_results(vnFitter[iPt], secPeak, useRefl, useTemplates, drawVnComps)
        hPulls.append(vnResults['pulls'])
        fTotFuncMass.append(vnResults['fTotFuncMass'])
        fTotFuncVn.append(vnResults['fTotFuncVn'])
        
        fSgnFuncMass.append(vnResults['fSgnFuncMass'])
        fBkgFuncMass.append(vnResults['fBkgFuncMass'])
        fBkgFuncVn.append(vnResults['fBkgFuncVn'])
        if secPeak:
            fMassSecPeakFunc.append(vnResults['fMassSecPeakFunc'])
            fVnSecPeakFunc.append(vnResults['fVnSecPeakFunct'])
        if useTemplates:
            fMassTemplFuncts[iPt] = vnResults['fMassTemplFuncts']
            fMassTemplTotFuncts[iPt] = vnResults['fMassTemplTotFunc']
        if drawVnComps:
            fVnCompFuncts.append(vnResults['fVnCompsFuncts'])

        if useRefl:
            fMassBkgRflFunc.append(vnResults['fMassBkgRflFunc'])
            hRel.append(vnResults['fMassRflFunc'])

        hSigmaSimFit.SetBinContent(iPt+1, vnResults['sigma'])
        hSigmaSimFit.SetBinError(iPt+1, vnResults['sigmaUnc'])
        hMeanSimFit.SetBinContent(iPt+1, vnResults['mean'])
        hMeanSimFit.SetBinError(iPt+1, vnResults['meanUnc'])
        hRedChi2SimFit.SetBinContent(iPt+1, vnResults['chi2'])
        hRedChi2SimFit.SetBinError(iPt+1, 1.e-20)
        hProbSimFit.SetBinContent(iPt+1, vnResults['prob'])
        hProbSimFit.SetBinError(iPt+1, 1.e-20)
        hRawYieldsSimFit.SetBinContent(iPt+1, vnResults['ry'])
        hRawYieldsSimFit.SetBinError(iPt+1, vnResults['ryUnc'])
        hRawYieldsTrueSimFit.SetBinContent(iPt+1, vnResults['ryTrue'])
        hRawYieldsTrueSimFit.SetBinError(iPt+1, vnResults['ryTrueUnc'])
        hRawYieldsSignificanceSimFit.SetBinContent(iPt+1, vnResults['signif'])
        hRawYieldsSignificanceSimFit.SetBinError(iPt+1, vnResults['signifUnc'])
        hvnSimFit.SetBinContent(iPt+1, vnResults['vn'])
        hvnSimFit.SetBinError(iPt+1, vnResults['vnUnc'])
        gvnSimFit.SetPoint(iPt, (ptMin+ptMax)/2, vnResults['vn'])
        gvnSimFit.SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2, vnResults['vnUnc'], vnResults['vnUnc'])
        gvnUnc.SetPoint(iPt, (ptMin+ptMax)/2, vnResults['vnUnc'])
        gvnUnc.SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2, 1.e-20, 1.e-20)

        
        if secPeak:
            hMeanSecPeakFitMass.SetBinContent(iPt+1, vnResults['secPeakMeanMass'])
            hMeanSecPeakFitMass.SetBinError(iPt+1, vnResults['secPeakMeanMassUnc'])
            hSigmaSecPeakFitMass.SetBinContent(iPt+1, vnResults['secPeakSigmaMass'])
            hSigmaSecPeakFitMass.SetBinError(iPt+1, vnResults['secPeakSigmaMassUnc'])
            hMeanSecPeakFitVn.SetBinContent(iPt+1, vnResults['secPeakMeanVn'])
            hMeanSecPeakFitVn.SetBinError(iPt+1, vnResults['secPeakMeanVnUnc'])
            hSigmaSecPeakFitVn.SetBinContent(iPt+1, vnResults['secPeakSigmaVn'])
            hSigmaSecPeakFitVn.SetBinError(iPt+1, vnResults['secPeakSigmaVnUnc'])
            gvnSimFitSecPeak.SetPoint(iPt, (ptMin+ptMax)/2, vnResults['vnSecPeak'])
            gvnSimFitSecPeak.SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2,
                                            vnResults['vnSecPeakUnc'],
                                            vnResults['vnSecPeakUnc'])
            gvnUncSecPeak.SetPoint(iPt, (ptMin+ptMax)/2, vnResults['vnSecPeakUnc'])
            gvnUncSecPeak.SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2, 1.e-20, 1.e-20)
        
        if useTemplates:
            print(f"Filling TemplOverSgn: {vnFitter[iPt].GetTemplOverSig()}")
            hTemplOverSgn.SetBinContent(iPt+1, vnFitter[iPt].GetTemplOverSig())
            for iTempl, (templVn, templVnUnc) in enumerate(zip(vnResults["vnTemplates"], vnResults["vnTemplatesUncs"])):
                gvnTempls[iTempl].SetPoint(iPt, (ptMin+ptMax)/2, templVn)
                gvnTempls[iTempl].SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2, 1.e-20, 1.e-20)
                gvnTemplsUncs[iTempl].SetPoint(iPt, (ptMin+ptMax)/2, templVnUnc)
                gvnTemplsUncs[iTempl].SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2, 1.e-20, 1.e-20)

        if vnResults['vn'] != 0:
            cSimFit[iPt].cd(1)
            hMassForFit[iPt].GetYaxis().SetRangeUser(0.2*hMassForFit[iPt].GetMinimum(),
                                                        1.6*hMassForFit[iPt].GetMaximum())
            hMassForFit[iPt].GetYaxis().SetMaxDigits(3)
            hMassForFit[iPt].GetXaxis().SetRangeUser(massMin, massMax)
            hMassForFit[iPt].Draw('E')
            SetObjectStyle(fBkgFuncMass[iPt], color=kOrange+1, linestyle=9, linewidth=2)
            SetObjectStyle(fTotFuncMass[iPt], color=kAzure+4, linewidth=3)
            SetObjectStyle(fSgnFuncMass[iPt], fillcolor=kAzure+4, fillstyle=3245, linewidth=0)
            if useRefl:
                SetObjectStyle(hRel[iPt], fillcolor=kGreen+1, fillstyle=3254, linewidth=0)
                SetObjectStyle(fMassBkgRflFunc[iPt], color=kRed+1, linestyle=7, linewidth=2)
            fSgnFuncMass[iPt].Draw('fc same')
            fBkgFuncMass[iPt].Draw('same')
            fTotFuncMass[iPt].Draw('same')
            if useRefl:
                fMassBkgRflFunc[iPt].Draw('same')
                hRel[iPt].Draw('same')
            if secPeak:
                SetObjectStyle(fMassSecPeakFunc[-1], fillcolor=kGreen+1, fillstyle=3254, linewidth=0)
                fMassSecPeakFunc[-1].Draw('same')
            latex.DrawLatex(0.18, 0.80, f'#mu = {vnResults["mean"]:.3f} #pm {vnResults["meanUnc"]:.3f} GeV/c^{2}')
            latex.DrawLatex(0.18, 0.75, f'#sigma = {vnResults["sigma"]:.3f} #pm {vnResults["sigmaUnc"]:.3f} GeV/c^{2}')
            latex.DrawLatex(0.18, 0.70, f'S = {vnResults["ry"]:.0f} #pm {vnResults["ryUnc"]:.0f}')
            latex.DrawLatex(0.18, 0.65, f'S/B (3#sigma) = {vnResults["ry"]/vnResults["bkg"]:.2f}')
            latex.DrawLatex(0.18, 0.60, f'Signif. (3#sigma) = {round(vnResults["signif"], 2)}')
            if useRefl:
                latex.DrawLatex(0.18, 0.20, f'RoverS = {SoverR:.2f}')
            if useTemplates:
                SetObjectStyle(fMassTemplTotFuncts[iPt], color=kRed, linewidth=2)
                fMassTemplTotFuncts[iPt].Draw('same')
                cSimFit[iCanv].Modified()
                cSimFit[iCanv].Update()
                if drawSingleTempls:
                    for iMassTemplFunct, massTemplFunct in enumerate(fMassTemplFuncts[iPt]):
                        SetObjectStyle(massTemplFunct, color=kMagenta+2+iMassTemplFunct*2, linewidth=3)
                        massTemplFunct.Draw('same')
                        cSimFit[iCanv].Modified()
                        cSimFit[iCanv].Update()
            cSimFit[iPt].cd(2)
            
            if drawDynRange:
                if hVnForFit[iPt].GetMinimum() > 0:
                    minCounts = hVnForFit[iPt].GetMinimum() - hVnForFit[iPt].GetBinError(hVnForFit[iPt].GetMinimumBin())
                else:
                    minCounts = hVnForFit[iPt].GetMinimum() + hVnForFit[iPt].GetBinError(hVnForFit[iPt].GetMinimumBin())
                maxCounts = hVnForFit[iPt].GetMaximum() + hVnForFit[iPt].GetBinError(hVnForFit[iPt].GetMaximumBin())
                vnMaxRange = maxCounts + 0.15 * maxCounts
                vnMinRange = minCounts - 0.15 * minCounts if minCounts > 0 else minCounts + 0.15 * minCounts
                hVnForFit[iPt].GetYaxis().SetRangeUser(vnMinRange, vnMaxRange)
            else:
                hVnForFit[iPt].GetYaxis().SetRangeUser(-0.2, 0.4)
            hVnForFit[iPt].GetYaxis().SetTitle(f'#it{{v}}_{{{harmonic}}} ({vn_method})')
            hVnForFit[iPt].GetXaxis().SetRangeUser(massMin, massMax)
            hVnForFit[iPt].Draw('E')
            SetObjectStyle(fBkgFuncVn[iPt], color=kOrange+1, linestyle=7, linewidth=2)
            SetObjectStyle(fTotFuncVn[iPt], color=kAzure+4, linewidth=3)
            fBkgFuncVn[iPt].Draw('same')
            fTotFuncVn[iPt].Draw('same')
            latex.DrawLatex(0.18, 0.23, f'#it{{R}}_{{{harmonic}}} = {reso:.3f}')
            latex.DrawLatex(0.18, 0.18, f'#chi^{{2}}/ndf = {vnResults["chi2"]:.2f}')
            latex.DrawLatex(0.18, 0.80,
                            f'#it{{v}}{harmonic}({particleName}) = {vnResults["vn"]:.3f} #pm {vnResults["vnUnc"]:.3f}')
            if secPeak and particleName == "Ds":
                latex.DrawLatex(0.18, 0.75,
                                f'#it{{v}}{harmonic}(D^{{+}}) = {vnResults["vnSecPeak"]:.3f} #pm {vnResults["vnSecPeakUnc"]:.3f}')
            if secPeak and particleName == "Dplus":
                latex.DrawLatex(0.18, 0.75,
                                f'#it{{v}}{harmonic}(D^{{*}}) = {vnResults["vnSecPeak"]:.3f} #pm {vnResults["vnSecPeakUnc"]:.3f}')
            if useTemplates:
                for iVnTempl, (vnCoeff, vnCoeffUnc) in enumerate(zip(vnResults["vnTemplates"], vnResults["vnTemplatesUncs"])):
                    latex.DrawLatex(0.18, 0.70-iVnTempl*0.05,
                                f'#it{{v}}{harmonic}(Templ{iVnTempl}) = {vnCoeff:.3f} #pm {vnCoeffUnc:.3f}')
                    cSimFit[iCanv].Modified()
                    cSimFit[iCanv].Update()
            if drawVnComps:
                legVnCompn = TLegend(0.72, 0.15, 0.9, 0.35)
                legVnCompn.SetBorderSize(0)
                legVnCompn.SetFillStyle(0)
                legVnCompn.SetTextSize(0.03)
                legVnCompn.AddEntry(fBkgFuncVn[iPt], f'#it{{v}}{harmonic} Bkg Func.', 'l')
                legVnCompn.AddEntry(fTotFuncVn[iPt], f'#it{{v}}{harmonic} Tot Func.', 'l')
                print(f"\n\n")
                print(f"fVnCompFuncts: {fVnCompFuncts}")
                print(f"\n\n")
                SetObjectStyle(fVnCompFuncts[iPt]['vnSgn'], fillcolor=kAzure+4, fillstyle=3245, linewidth=0)
                legVnCompn.AddEntry(fVnCompFuncts[iPt]['vnSgn'], f"Signal #it{{v}}{harmonic}", 'f')
                SetObjectStyle(fVnCompFuncts[iPt]['vnBkg'], color=kOrange+1, linestyle=1, linewidth=2)
                legVnCompn.AddEntry(fVnCompFuncts[iPt]['vnBkg'], f"Bkg #it{{v}}{harmonic}", 'l')
                if secPeak:
                    SetObjectStyle(fVnCompFuncts[iPt]['vnSecPeak'], fillcolor=kGreen+1, fillstyle=3254, linewidth=0)
                    legVnCompn.AddEntry(fVnCompFuncts[iPt]['vnSecPeak'], f"Second peak #it{{v}}{harmonic}", 'f')
                if useTemplates:
                    for iTempl in range(len(fVnCompFuncts[iPt])-2-secPeak):
                        SetObjectStyle(fVnCompFuncts[iPt][f'vnTempl{iTempl}'], color=kMagenta+2+iTempl*2, linewidth=3)
                        legVnCompn.AddEntry(fVnCompFuncts[iPt][f'vnTempl{iTempl}'], f"Templ{iTempl} #it{{v}}{harmonic}", 'l')
                for vnCompKey, vnCompFunct in fVnCompFuncts[iPt].items():
                    vnCompFunct.Draw('same')
                    cSimFit[iCanv].Modified()
                    cSimFit[iCanv].Update()
                legVnCompn.Draw()

            cSimFit[iCanv].Modified()
            cSimFit[iCanv].Update()

        invMassPrefit = vnFitter[iPt].GetMassPrefitObject()
        hPullsPrefit.append(invMassPrefit.GetPullDistribution())
        histoMassPrefit = invMassPrefit.GetHistoClone()
        totFuncMassPrefit = invMassPrefit.GetMassFunc()
        bkgFuncMassPrefit = invMassPrefit.GetBackgroundRecalcFunc()
        sgnFuncMassPrefit = invMassPrefit.GetSignalFunc()
        cInvMassPrefits[iPt] = TCanvas(f"cMass_{ptMin*10:.0f}_{ptMax*10:.0f}", f"Mass Fit {ptMin}-{ptMax} GeV/c", 800, 600)
        histoMassPrefit.SetStats(0)
        histoMassPrefit.Draw("E")
        bkgFuncMassPrefit.SetLineColor(kGreen+2)
        bkgFuncMassPrefit.SetLineWidth(2)
        bkgFuncMassPrefit.SetLineWidth(3)
        bkgFuncMassPrefit.Draw("same")
        sgnFuncMassPrefit.SetLineColor(kBlue)
        sgnFuncMassPrefit.SetLineWidth(2)
        sgnFuncMassPrefit.SetLineWidth(3)
        sgnFuncMassPrefit.Draw("same")
        if useTemplates:
            templFuncMassPrefit = invMassPrefit.GetTemplFunc()
            print(f"invMassPrefit.GetTemplOverSig(): {invMassPrefit.GetTemplOverSig()}")
            print(f"vnFitter[iPt].GetTemplOverSig(): {vnFitter[iPt].GetTemplOverSig()}")
            templFuncMassPrefit.SetLineColor(kMagenta)
            templFuncMassPrefit.SetLineWidth(2)
            templFuncMassPrefit.SetLineWidth(3)
            templFuncMassPrefit.Draw("same")

        totFuncMassPrefit.SetLineColor(kRed)
        totFuncMassPrefit.SetLineWidth(2)
        totFuncMassPrefit.SetLineWidth(3)
        totFuncMassPrefit.Draw("same")
    
    canvVn.cd().SetLogx()
    hframe = canvVn.DrawFrame(0.5, -0.5, gvnSimFit.GetXaxis().GetXmax()+0.5, 0.5,
                              f';#it{{p}}_{{T}} (GeV/c); v_{{{harmonic}}} ({vn_method})')
    hframe.GetYaxis().SetDecimals()
    hframe.GetXaxis().SetNdivisions(504)
    hframe.GetXaxis().SetMoreLogLabels()
    gPad.SetGridy()
    gvnSimFit.Draw('same pez')
    if secPeak:
        gvnSimFitSecPeak.Draw('pez same')
    latex.DrawLatex(0.2, 0.2, f'#it{{R}}_{{{harmonic}}} = {reso:.3f}')
    latex.DrawLatexNDC(0.20, 0.80, 'This work')
    latex.DrawLatexNDC(0.20, 0.75, f'Pb#minusPb #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV ({centMinMax[0]}#minus{centMinMax[1]}%)')
    latex.DrawLatexNDC(0.20, 0.70, decay)
    canvVn.Modified()
    canvVn.Update()
    canvVnUnc.cd()
    gvnUnc.Draw('apez same')
    if secPeak:
        gvnUncSecPeak.Draw('pez same')
    canvVnUnc.Modified()
    canvVnUnc.Update()
    if not batch:
        input('Press Enter to continue...')

    #save output histos
    print(f'Saving output to {outputdir}')
    for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
        if iPt == 0:
            suffix_pdf = '('
        elif iPt == nPtBins-1:
            suffix_pdf = ')'
        else:
            suffix_pdf = ''
        if nPtBins==1:
            cSimFit[iPt].SaveAs(f'{outputdir}/SimFit{suffix}_{particleName}.pdf')
        else:
            cSimFit[iPt].SaveAs(f'{outputdir}/SimFit{suffix}_{particleName}.pdf{suffix_pdf}')

    outfile_name = f'{outputdir}/raw_yields{suffix}.root'
    outFile = TFile(outfile_name, 'recreate')
    for canv in cSimFit:
        canv.Write()
    for canvPrefit in cInvMassPrefits:
        canvPrefit.Write()
    for hist in hMass:
        hist.Write('hist_mass')
    for hist in hPulls:
        hist.Write('hist_pulls')
    for hist in hPullsPrefit:
        hist.Write('hist_pulls_prefit')
    for hist in hVn:
        hist.Write('hist_vn')
    for ipt, (ptmin, ptmax) in enumerate(zip(ptMins, ptMaxs)):
        try:
            fTotFuncMass[ipt].Write(f'fTotFuncMass_pt{ptmin*10:.0f}_{ptmax*10:.0f}')
            fTotFuncVn[ipt].Write(f'fTotFuncVn_pt{ptmin*10:.0f}_{ptmax*10:.0f}')
            fSgnFuncMass[ipt].Write(f'fSgnFuncMass_pt{ptmin*10:.0f}_{ptmax*10:.0f}')
            fBkgFuncMass[ipt].Write(f'fBkgFuncMass_pt{ptmin*10:.0f}_{ptmax*10:.0f}')
            fBkgFuncVn[ipt].Write(f'fBkgFuncVn_pt{ptmin*10:.0f}_{ptmax*10:.0f}')
        except Exception:
            print(f'WARNING: Fit failed for pt {ptmin*10:.0f}-{ptmax*10:.0f}')
            
            
    hSigmaSimFit.Write()
    hMeanSimFit.Write()
    hMeanSecPeakFitMass.Write()
    hMeanSecPeakFitVn.Write()
    hTemplOverSgn.Write()
    hSigmaSecPeakFitMass.Write()
    hSigmaSecPeakFitVn.Write()
    hRawYieldsSimFit.Write()
    hRawYieldsTrueSimFit.Write()
    hRawYieldsSecPeakSimFit.Write()
    hRawYieldsSignificanceSimFit.Write()
    hRawYieldsSoverBSimFit.Write()
    hRedChi2SimFit.Write()
    hProbSimFit.Write()
    hRedChi2SBVnPrefit.Write()
    hProbSBVnPrefit.Write()
    hvnSimFit.Write()
   
    gvnSimFit.Write()
    gvnUnc.Write()
    if secPeak:
        gvnSimFitSecPeak.Write()
        gvnUncSecPeak.Write()
    if includeTempls:
        for iTempl in range(len(templsNames)):
            gvnTempls[iTempl].Write()
            gvnTemplsUncs[iTempl].Write()
    hist_reso.Write()

    outFile.Close()

    # Debug: Check ROOT objects before exit
    print("\n=== DEBUG: ROOT Object Cleanup Checks ===")
    try:
        import ROOT
        print(f"Number of active canvases: {len(ROOT.gROOT.GetListOfCanvases())}")
        print(f"Memory usage: {ROOT.gSystem.GetMemoryInfo()}")
        
        # Explicit cleanup
        print("Performing explicit cleanup...")
        for c in cSimFit + cInvMassPrefits + [canvVn, canvVnUnc]:
            if c and c.IsOnHeap():
                c.Close()
        
        # Clear lists
        cSimFit.clear()
        cInvMassPrefits.clear()
        
        print("Cleanup complete. Exiting.")
    except Exception as e:
        print(f"Error during cleanup: {str(e)}")
    
    if not batch:
        input('Press enter to exit')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('flow_config', metavar='text', default='config_flow.yml')
    parser.add_argument('centClass', metavar='text', default='')
    parser.add_argument('inFileName', metavar='text', default='')
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument('--batch', help='suppress video output', action='store_true')
    args = parser.parse_args()

    get_vn_vs_mass(
        args.flow_config,
        args.centClass,
        args.inFileName,
        args.outputdir,
        args.suffix,
        args.batch
    )