'''
Analysis utilities for flow analysis
'''
import sys
import ctypes
from .check import check_file_exists
from ROOT import TFile, TDatabasePDG, TH1

def get_centrality_bins(centrality):
    '''
    Get centrality bins

    Input:
        - centrality:
            str, centrality class (e.g. 'k3050')

    Output:
        - cent_bins:
            list of floats, centrality bins
        - cent_label:
            str, centrality label
    '''
    if centrality == 'k010':
        return '0_10', [0, 10]
    if centrality == 'k020':
        return '0_20', [0, 20]
    if centrality == 'k2030':
        return '20_30', [20, 30]
    elif centrality == 'k3040':
        return '30_40', [30, 40]
    elif centrality == 'k3050':
        return '30_50', [30, 50]
    elif centrality == 'k4050':
        return '40_50', [40, 50]
    elif centrality == 'k2060':
        return '20_60', [20, 60]
    elif centrality == 'k4060':
        return '40_60', [40, 60]
    elif centrality == 'k6070':
        return '60_70', [60, 70]
    elif centrality == 'k6080':
        return '60_80', [60, 80]
    elif centrality == 'k7080':
        return '70_80', [70, 80]
    elif centrality == 'k0100':
        return '0_100', [0, 100]
    elif centrality == 'k5080':
        return '50_80', [50, 80]
    else:
        print(f"ERROR: cent class \'{centrality}\' is not supported! Exit")
    sys.exit()
    
def get_vnfitter_results(vnFitter, secPeak, useRefl, useTempl, DrawVnComps):
    '''
    Get vn fitter results:
    0: BkgInt
    1: BkgSlope
    2: SgnInt
    3: Mean
    4: Sigma
    5: SecPeakInt
    6: SecPeakMean
    7: SecPeakSigma
    8: ConstVnBkg
    9: SlopeVnBkg
    10: v2Sgn
    11: v2SecPeak
    12: reflection

    Input:
        - vnfitter:
            VnVsMassFitter, vn fitter object
        - secPeak:
            bool, if True, save secondary peak results
        - useRefl:
            bool, if True, save the results with reflection
        - useTempl:
            bool, if True, save the results with templates for corelated background
        - DrawVnComps:
            bool, if True, save the vn components functions

    Output:
        - vn_results:
            dict, dictionary with vn results
            vn: vn value
            vnUnc: uncertainty of vn value
            mean: mean value
            meanUnc: uncertainty of mean value
            sigma: sigma value
            sigmaUnc: uncertainty of sigma value
            ry: raw yield
            ryUnc: uncertainty of raw yield
            ryTrue: true raw yield
            ryTrueUnc: uncertainty of true raw yield
            signif: significance
            signifUnc: uncertainty of significance
            chi2: reduced chi2
            prob: fit probability
            fTotFuncMass: total fit function for mass
            fTotFuncVn: total fit function for vn
            secPeakMeanMass: secondary peak mean mass
            secPeakMeanMassUnc: uncertainty of secondary peak mean mass
            secPeakSigmaMass: secondary peak sigma mass
            secPeakSigmaMassUnc: uncertainty of secondary peak sigma mass
            secPeakMeanVn: secondary peak mean vn
            secPeakMeanVnUnc: uncertainty of secondary peak mean vn
            secPeakSigmaVn: secondary peak sigma vn
            secPeakSigmaVnUnc: uncertainty of secondary peak sigma vn
            vnSecPeak: vn secondary peak
            vnSecPeakUnc: uncertainty of vn secondary peak
            fMassRflFunc: mass reflection function
            fMassBkgRflFunc: mass background reflection function
            fMassSecPeakFunc: mass secondary peak function
            fVnSecPeakFunct: vn secondary peak function
            fVnCompsFuncts: dictionary with vn components functions
            fMassTemplFuncts: dictionary with mass template functions
            vnTemplates: list of vn templates
            vnTemplatesUncs: list of vn templates uncertainties
    '''
    vn_results = {}
    vn_results['vn'] = vnFitter.GetVn()
    vn_results['vnUnc'] = vnFitter.GetVnUncertainty()
    vn_results['mean'] = vnFitter.GetMean()
    vn_results['meanUnc'] = vnFitter.GetMeanUncertainty()
    vn_results['sigma'] = vnFitter.GetSigma()
    vn_results['sigmaUnc'] = vnFitter.GetSigmaUncertainty()
    vn_results['ry'] = vnFitter.GetRawYield()
    vn_results['ryUnc'] = vnFitter.GetRawYieldUncertainty()
    vn_results['chi2'] = vnFitter.GetReducedChiSquare()
    vn_results['prob'] = vnFitter.GetFitProbability()
    vn_results['fTotFuncMass'] = vnFitter.GetMassTotFitFunc()
    vn_results['fTotFuncVn'] = vnFitter.GetVnVsMassTotFitFunc()
    vn_results['fBkgFuncMass'] = vnFitter.GetMassBkgFitFunc()
    vn_results['fBkgFuncVn'] = vnFitter.GetVnVsMassBkgFitFunc()
    vn_results['fSgnFuncMass'] = vnFitter.GetMassSignalFitFunc()
    vn_results['pulls'] = vnFitter.GetPullDistribution()

    if DrawVnComps:
        vn_results['fVnCompsFuncts'] = {}
        vnComps = vnFitter.GetVnCompsFuncts()
        vn_results['fVnCompsFuncts']['vnSgn'] = vnComps[0]
        vn_results['fVnCompsFuncts']['vnBkg'] = vnComps[1]
        if useTempl:
            for iTempl in range(len(vn_results['fMassTemplFuncts'])):
                vn_results['fVnCompsFuncts'][f'vnTempl{iTempl}'] = vnComps[2+secPeak+iTempl]
    if secPeak and DrawVnComps:
        vn_results['fVnCompsFuncts']['vnSecPeak'] = vnComps[2]
    vn_results['fMassTemplFuncts'] = [] 
    if useTempl:
        vn_results['fMassTemplTotFunc'] = vnFitter.GetMassTemplFitFunc()
        vn_results['fMassTemplFuncts'] = vnFitter.GetMassTemplFuncts()
    
    bkg, bkgUnc = ctypes.c_double(), ctypes.c_double()
    vnFitter.Background(3, bkg, bkgUnc)
    vn_results['bkg'] = bkg.value
    vn_results['bkgUnc'] = bkgUnc.value
    sgn, sgnUnc = ctypes.c_double(), ctypes.c_double()
    vnFitter.Signal(3, sgn, sgnUnc)
    vn_results['ryTrue'] = sgn.value
    vn_results['ryTrueUnc'] = sgnUnc.value
    signif, signifUnc = ctypes.c_double(), ctypes.c_double()
    vnFitter.Significance(3, signif, signifUnc)
    vn_results['signif'] = signif.value
    vn_results['signifUnc'] = signifUnc.value

    massSgnPars = vnFitter.GetNMassSgnPars()
    massBkgPars = vnFitter.GetNMassBkgPars()
    massSecPeakPars = vnFitter.GetNMassSecPeakPars()
    massReflPars = vnFitter.GetNMassReflPars()
    massTemplPars = len(vn_results['fMassTemplFuncts'])
    totMassPars = massSgnPars + massBkgPars + massSecPeakPars +  massReflPars + massTemplPars
    vnSgnPars = vnFitter.GetNVnSgnPars()
    vnBkgPars = vnFitter.GetNVnBkgPars()

    if secPeak:
        vn_results['fMassSecPeakFunc'] = vnFitter.GetMassSecPeakFunc()
        vn_results['fVnSecPeakFunct'] = vnFitter.GetVnSecPeakFunc()
        vn_results['secPeakMeanMass'] = vn_results['fTotFuncMass'].GetParameter(vn_results['fTotFuncMass'].GetParName(massSgnPars + massBkgPars + 1))
        vn_results['secPeakMeanMassUnc'] = vn_results['fTotFuncMass'].GetParError(massSgnPars + massBkgPars + 1)
        vn_results['secPeakSigmaMass'] = vn_results['fTotFuncMass'].GetParameter(vn_results['fTotFuncMass'].GetParName(massSgnPars + massBkgPars + 2))
        vn_results['secPeakSigmaMassUnc'] = vn_results['fTotFuncMass'].GetParError(massSgnPars + massBkgPars + 2)
        vn_results['secPeakMeanVn'] = vn_results['fTotFuncVn'].GetParameter(vn_results['fTotFuncVn'].GetParName(totMassPars + vnSgnPars + vnBkgPars + 1))
        vn_results['secPeakMeanVnUnc'] = vn_results['fTotFuncVn'].GetParError(vnSgnPars + vnBkgPars + 1)
        vn_results['secPeakSigmaVn'] = vn_results['fTotFuncVn'].GetParameter(vn_results['fTotFuncVn'].GetParName(totMassPars + vnSgnPars + vnBkgPars + 2))
        vn_results['secPeakSigmaVnUnc'] = vn_results['fTotFuncVn'].GetParError(vnSgnPars + vnBkgPars + 2)
        vn_results['vnSecPeak'] = vn_results['fTotFuncVn'].GetParameter(vn_results['fTotFuncVn'].GetParName(totMassPars + vnSgnPars + vnBkgPars))
        vn_results['vnSecPeakUnc'] = vn_results['fTotFuncVn'].GetParError(totMassPars + vnSgnPars + vnBkgPars)

    if useRefl:
        vn_results['fMassRflFunc'] = vnFitter.GetMassRflFunc()
        vn_results['fMassBkgRflFunc'] = vnFitter.GetMassBkgRflFunc()
    
    if useTempl:
        vn_results['vnTemplates'] = list(vnFitter.GetVnTemplates())
        vn_results['vnTemplatesUncs'] = list(vnFitter.GetVnTemplatesUncertainties())

    return vn_results

def get_refl_histo(reflFile, centMinMax, ptMins, ptMaxs):
    '''
    Method that loads MC histograms for the reflections of D0

    Input:
        - reflFile:
           TFile, ROOT file, include reflections of D0
        - centMinMax:
            list, min and max centrality
        - ptMins:
            list, min pt bins
        - ptMaxs:
            list, max pt bins
    
    Output:
        - useRefl:
            bool, if True, MC histograms for the reflections of D0 exists
        - hMCSgn:
            lsit, signal histograms of D0
        - hMCRefl:
            list, reflection histograms of D0
    '''
    hMCSgn, hMCRefl = [], []
    if not check_file_exists(reflFile):
        print(f'Error: reflection file {reflFile} does not exist! Turning off reflections usage')
        return False
    
    reflFile = TFile(reflFile, 'READ')
    for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
        dirName = f'cent_bins{centMinMax[0]}_{centMinMax[1]}/pt_bins{ptMin}_{ptMax}'
        if not reflFile.GetDirectory(dirName):
            print(f'In directory {dirName}, no histograms found! Turning off reflections usage')
            return False

        hMCSgn.append(reflFile.Get(f'{dirName}/hFDMass'))
        hMCSgn[iPt].Add(reflFile.Get(f'{dirName}/hPromptMass'), 1)
        if not isinstance(hMCSgn[iPt], TH1) or hMCSgn[iPt] == None:
            print(f'In directory {dirName}, hMCSgnMass_{ptMin}_{ptMax} not found! Turning off reflections usage')
            return False
        hMCSgn[iPt].SetName(f'histSgn_{iPt}')
        hMCSgn[iPt].SetDirectory(0)

        hMCRefl.append(reflFile.Get(f'{dirName}/hReflMass'))
        if not isinstance(hMCRefl[iPt], TH1) or hMCRefl[iPt] == None:
            print(f'In directory {dirName}, hReflMass_{ptMin}_{ptMax} not found! Turning off reflections usage')
            return False
        hMCRefl[iPt].SetName(f'histRfl_{iPt}')
        hMCRefl[iPt].SetDirectory(0)

        if hMCRefl[iPt].Integral() <= 0 or hMCSgn[iPt].Integral() <= 0:
            print(f'WARNING: Empty reflection template for pt bin {ptMin}-{ptMax} GeV/c')
            return False

    reflFile.Close()

    return True, hMCSgn, hMCRefl

def get_particle_info(particleName):
    '''
    Get particle information

    Input:
        - particleName: 
            the name of the particle

    Output:
        - particleTit: 
            the title of the particle
        - massAxisTit: 
            the title of the mass axis
        - decay: 
            the decay of the particle
        - massForFit: 
            float, the mass of the particle
    '''

    if particleName == 'Dplus':
        particleTit = 'D^{+}'
        massAxisTit = '#it{M}(K#pi#pi) (GeV/#it{c}^{2})'
        massForFit = TDatabasePDG.Instance().GetParticle(411).Mass()
        decay = 'D^{+} #rightarrow K^{#minus}#pi^{+}#pi^{+}'
    elif particleName == 'Ds':
        particleTit = 'D_{s}^{+}'
        massAxisTit = '#it{M}(KK#pi) (GeV/#it{c}^{2})'
        decay = 'D_{s}^{+} #rightarrow #phi#pi^{+} #rightarrow K^{+}K^{#minus}#pi^{+}'
        massForFit = TDatabasePDG.Instance().GetParticle(431).Mass()
    elif particleName == 'LctopKpi':
        particleTit = '#Lambda_{c}^{+}'
        massAxisTit = '#it{M}(pK#pi) (GeV/#it{c}^{2})'
        decay = '#Lambda_{c}^{+} #rightarrow pK^{#minus}#pi^{+}'
        massForFit = TDatabasePDG.Instance().GetParticle(4122).Mass()
    elif particleName == 'LctopK0s':
        massAxisTit = '#it{M}(pK^{0}_{s}) (GeV/#it{c}^{2})'
        decay = '#Lambda_{c}^{+} #rightarrow pK^{0}_{s}'
        massForFit = 2.25 # please calfully check the mass of Lc->pK0s, it is constant
        # massForFit = ROOT.TDatabasePDG.Instance().GetParticle(4122).Mass()
    elif particleName == 'Dstar':
        particleTit = 'D^{*+}'
        massAxisTit = '#it{M}(K#pi#pi) - #it{M}(K#pi) (GeV/#it{c}^{2})'
        decay = 'D^{*+} #rightarrow D^{0}#pi^{+} #rightarrow K^{#minus}#pi^{+}#pi^{+}'
        massForFit = TDatabasePDG.Instance().GetParticle(413).Mass() - TDatabasePDG.Instance().GetParticle(421).Mass()
    elif particleName == 'Dzero':
        particleTit = 'D^{0}'
        massAxisTit = '#it{M}(K#pi) (GeV/#it{c}^{2})'
        decay = 'D^{0} #rightarrow K^{#minus}#pi^{+}'
        massForFit = TDatabasePDG.Instance().GetParticle(421).Mass()
    else:
        print(f'ERROR: the particle "{particleName}" is not supported! Choose between Dzero, Dplus, Ds, Dstar, and Lc. Exit!')
        sys.exit()

    return particleTit, massAxisTit, decay, massForFit