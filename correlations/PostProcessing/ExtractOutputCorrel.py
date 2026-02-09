import ROOT
import json
import os
import sys
import yaml
import argparse
from ROOT import gSystem
script_dir = os.path.dirname(os.path.realpath(__file__))
gSystem.CompileMacro(f"{script_dir}/DhCorrelationExtraction.cxx", "kO")
from ROOT import DhCorrelationExtraction as CorrelExtractor
from itertools import product
from concurrent.futures import ProcessPoolExecutor, as_completed
from gc import collect
import numpy as np

def SetCanvasStyle():
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadBottomMargin(0.15)
    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetCanvasDefH(1126)
    ROOT.gStyle.SetCanvasDefW(1840)

def get_pt_dependent_param(param, nPtBins, isList=False):
    ''' handle parameters that can be either a single value or a list of values for each pt bin.
    If isList is True, the parameter is expected to be a single value to be repeated for all pt bins.
    If isList is False, the parameter can be a list of values or a single value to be repeated.
    '''
    if isList:
        parameter = [param] * nPtBins
    else:
        if isinstance(param, list):
            if len(param) < nPtBins:
                print(f"[ERROR] Length of parameter list {len(param)} does not match number of pt bins {nPtBins}")
                sys.exit(1)
        else:
            parameter = [param] * nPtBins
    return parameter

def process_correlation_task(tempExtractor, task):
    extractor = CorrelExtractor.CreateCopy(tempExtractor)
    ROOT.SetOwnership(extractor, True)

    iPtCand = task["iPtCand"]
    ptMin = task["ptMin"]
    ptMax = task["ptMax"]
    rebinDeltaEta = task["rebinDeltaEta"]
    rebinDeltaPhi = task["rebinDeltaPhi"]
    iPtHad = task["iPtHad"]
    ptHadMin = task["ptHadMin"]
    ptHadMax = task["ptHadMax"]
    iMass = task["iMass"]
    invMassMin = task["invMassMin"]
    invMassMax = task["invMassMax"]

    results = {}
    results["task"] = task

    extractor.SetRebin2DcorrelHisto(rebinDeltaEta, rebinDeltaPhi)
    extractor.SetCandAndHadBins((ptMin, ptMax), (ptHadMin, ptHadMax))
    extractor.SetInvMassBins((invMassMin, invMassMax))
    if task["method"] == "DeltaPhiBinning":
        extractor.SetMethod(CorrelExtractor.kDeltaPhiBinning)
        deltaPhiMin = task["deltaPhiMin"]
        deltaPhiMax = task["deltaPhiMax"]
        extractor.SetDeltaPhiBins((deltaPhiMin, deltaPhiMax))
    extractor.ExtractCorrelations()

    results['hCorrectedCorrHisto'] = extractor.GetCorrectedCorrHisto()
    results['hCorrectedCorrHisto_BaselineSubtr'] = extractor.GetCorrectedCorrHisto_BaselineSubtr()
    results['hCorrectedCorrHisto_Reflected'] = extractor.GetCorrectedCorrHisto_Reflected()
    results['hCorrectedCorrHisto_Reflected_BaselineSubtr'] = extractor.GetCorrectedCorrHisto_Reflected_BaselineSubtr()

    results['hCorrel_SE_2D'] = extractor.GetCorrel_SE_2D()
    results['hCorrel_ME_2D'] = extractor.GetCorrel_ME_2D()
    results['hCorrectedCorrel_2D'] = extractor.GetCorrectedCorrel_2D()
    results['hNonNormalizedCorrHisto'] = extractor.GetNonNormalizedCorrHisto()

    results['hVecOriginalCorrel_SE_2D'] = extractor.GetVecOriginalCorrel_SE_2D()
    results['hVecOriginalCorrel_ME_2D'] = extractor.GetVecOriginalCorrel_ME_2D()
    results['hVecCorrel_SE_2D'] = extractor.GetVecCorrel_SE_2D()
    results['hVecCorrel_ME_2D'] = extractor.GetVecCorrel_ME_2D()
    results['hVecCorrel_ME_norm_2D'] = extractor.GetVecCorrel_ME_norm_2D()
    results['hVecCorrectedCorrel_2D'] = extractor.GetVecCorrectedCorrel_2D()

#   std::vector<TH1D*> GetVecCorrectedRatioVsDeltaEta() { return fVecCorrectedRatioVsDeltaEta; }
#   std::vector<std::vector<TH1D*>> GetVecVecCorrectedMassPairsVsDeltaEta() { return fVecVecCorrectedMassPairsVsDeltaEta; }
    if task["method"] == "DeltaPhiBinning":
        results['hCorrectedMassPairs'] = extractor.GetCorrectedMassPairs()
        results['hVecVecCorrectedMassPairsVsDeltaEta'] = extractor.GetVecVecCorrectedMassPairsVsDeltaEta()
        results['hVecCorrectedRatioVsDeltaEta'] = extractor.GetVecCorrectedRatioVsDeltaEta()
        results['hCorrectedRatioVsDeltaEta'] = extractor.GetCorrectedRatioVsDeltaEta()
        results['hVecMassPairsVsDeltaEta_2D'] = extractor.GetVecMassPairsVsDeltaEta_2D()

    del extractor
    collect()

    return results

def ExtractOutputCorrel(cfgFile):
    SetCanvasStyle()

    with open(cfgFile, 'r') as f:
        config = yaml.safe_load(f)

    #  Input files
    pathFileSE = config["pathFileSE"]
    pathFileME = config["pathFileME"]
    print(f"[INFO] Using SE file: {pathFileSE}")
    # pathFileMass = config["pathFileMass"]
    pathFileSecPart = config.get("pathFileSecPart", "")

    # General info for corelation extraction
    Dmeson = config["Dmeson"]
    deltaEtaBins = config["deltaEtaBins"]
    if deltaEtaBins[0][1] > deltaEtaBins[1][0]:
        print(f"[ERROR] The deltaEta bins for correlation extraction overlap: {deltaEtaBins}")
        sys.exit(1)
    valDeltaPhiMEnorm = config.get("ValDeltaPhiMEnorm", 0)
    valDeltaEtaMEnorm = config.get("ValDeltaEtaMEnorm", 0)
    nPools = config.get("nPools", 10)
    doPoolByPool = config.get("doPoolByPool", False)
    method = config.get("method", "MassBinning")

    # Binning operations
    ptBinsCand = config["ptBinsCand"]
    ptBinsHad = config["ptBinsHad"]
    nPtBinsCand = len(ptBinsCand) - 1
    invMassBins = config["invMassBins"]
    deltaPhiBins = config.get("deltaPhiBins", list(np.linspace(-1.571, 4.712, 24)) if method == "DeltaPhiBinning" else [-9,9])  # default 64 bins from -pi/2 to 3pi/2
    rebinsDeltaEta = get_pt_dependent_param(config.get("rebinsDeltaEta", 1), len(ptBinsCand)-1, isList=False)
    rebinsDeltaPhi = get_pt_dependent_param(config.get("rebinsDeltaPhi", 1), len(ptBinsCand)-1, isList=False)

    # Optional settings
    doSecPartContamination = config.get("doSecPartContamination", False)
    doRebinSecPart = config.get("doRebinSecPart", False)
    debug = config.get("debug", 0)

    # Define the default extractor
    tempExtractor = CorrelExtractor.CreateDefault() # aotumatic settings for directories and sparse names
    tempExtractor.SetInputFilenameSE(pathFileSE)
    tempExtractor.SetInputFilenameME(pathFileME)
    tempExtractor.SetInputFilenameMass(pathFileSE) # mass sparse is in the same file as SE
    if Dmeson == "D0" or Dmeson == "Dzero":
        tempExtractor.SetDmesonSpecie(0)
    elif Dmeson == "Dplus":
        tempExtractor.SetDmesonSpecie(1)
    elif Dmeson == "Ds":
        tempExtractor.SetDmesonSpecie(2)
        if pathFileSecPart != "" and doSecPartContamination:
            tempExtractor.SetInputFilenameSecPart(pathFileSecPart)
    else:
        print(f"[ERROR] Unknown D meson specie {Dmeson}")
        sys.exit(1)
    tempExtractor.SetPoolSettings(nPools, doPoolByPool)
    if method == "MassBinning":
        tempExtractor.SetMethod(CorrelExtractor.kMassBinning)
    elif method == "DeltaPhiBinning":
        tempExtractor.SetMethod(CorrelExtractor.kDeltaPhiBinning)
        if any(len (invMassBins[i]) > 2 for i in range(len(invMassBins))):
            print(f"Using DeltaPhiBinning method, but invMassBins has more than 2 edges, replace with 2 edges for min and max")
            invMassBins = [[invMassBins[i][0], invMassBins[i][-1]] for i in range(len(invMassBins))]
    else:
        print(f"[ERROR] Unknown extraction method {method}")
        sys.exit(1)
    tempExtractor.SetBinDeltaEtaLeft(deltaEtaBins[0][0], deltaEtaBins[0][1])
    tempExtractor.SetBinDeltaEtaRight(deltaEtaBins[1][0], deltaEtaBins[1][1])
    tempExtractor.SetBinDeltaPhiEtaForMEnorm(valDeltaPhiMEnorm, valDeltaEtaMEnorm)
    tempExtractor.SetSecPartContamination(doSecPartContamination)
    tempExtractor.SetdoRebinSecondaryPart(doRebinSecPart)
    tempExtractor.SetDebugLevel(debug)
    

    outdir = config["outdir"]
    suffix = config["suffix"]
    outdirFull = os.path.join(outdir, f"CorrelExtract_{suffix}")
    if not os.path.exists(outdirFull):
        os.makedirs(outdirFull)
    outdirMass = os.path.join(outdir, "InvMass")
    if not os.path.exists(outdirMass):
        os.makedirs(outdirMass)
    
    # mass vs pt
    tempExtractor.ProjMassVsPt()
    hMassVsPt = tempExtractor.GetMassVsPtHist2D()
    outMassVsPtFile = ROOT.TFile(os.path.join(outdirMass, f"InvMassVsPt.root"), "RECREATE")
    hMassVsPt.Write()
    outMassVsPtFile.Close()

    # extract for all pt bins, associated hadron pt bins, and inv. mass bins
    tasks = []
    for iPtCand, (ptMin, ptMax, rebinDeltaEta, rebinDeltaPhi) in enumerate(zip(ptBinsCand[:-1], ptBinsCand[1:], rebinsDeltaEta, rebinsDeltaPhi)):
        for iPtHad, (ptHadMin, ptHadMax) in enumerate(zip(ptBinsHad[:-1], ptBinsHad[1:])):
            invMassBinsPtCand = invMassBins[iPtCand]
            for iMass, (invMassMin, invMassMax) in enumerate(zip(invMassBinsPtCand[:-1], invMassBinsPtCand[1:])):
                if method == "MassBinning":
                    task = {
                        "iPtCand": iPtCand, "ptMin": ptMin, "ptMax": ptMax,
                        "rebinDeltaEta": rebinDeltaEta, "rebinDeltaPhi": rebinDeltaPhi,
                        "iPtHad": iPtHad, "ptHadMin": ptHadMin, "ptHadMax": ptHadMax,
                        "iMass": iMass, "invMassMin": invMassMin, "invMassMax": invMassMax,
                        "method": method
                    }
                    tasks.append(task)
                elif method == "DeltaPhiBinning":
                    for iDeltaPhi, (deltaPhiMin, deltaPhiMax) in enumerate(zip(deltaPhiBins[:-1], deltaPhiBins[1:])):
                        task = {
                            "iPtCand": iPtCand, "ptMin": ptMin, "ptMax": ptMax,
                            "rebinDeltaEta": rebinDeltaEta, "rebinDeltaPhi": rebinDeltaPhi,
                            "iPtHad": iPtHad, "ptHadMin": ptHadMin, "ptHadMax": ptHadMax,
                            "iMass": iMass, "invMassMin": invMassMin, "invMassMax": invMassMax,
                            "iDeltaPhi": iDeltaPhi, "deltaPhiMin": deltaPhiMin, "deltaPhiMax": deltaPhiMax,
                            "method": method
                        }
                        tasks.append(task)

    all_results = []
    with ProcessPoolExecutor(max_workers=24) as executor:
        future_to_task = {executor.submit(process_correlation_task, tempExtractor, task): task for task in tasks}
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                result = future.result()
                all_results.append(result)
            except Exception as exc:
                print(f"[ERROR] Task {task} generated an exception: {exc}")

    if method == "MassBinning":
        all_results.sort(key=lambda x: (x['task']['iMass'], x['task']['iPtHad'], x['task']['iPtCand']))
    elif method == "DeltaPhiBinning":
        all_results.sort(key=lambda x: (x['task']['iMass'], x['task']['iPtHad'], x['task']['iPtCand'], x['task']['iDeltaPhi']))

    # Save outputs
    outdirCorrelation = os.path.join(outdirFull, "CorrelationsResults")
    if not os.path.exists(outdirCorrelation):
        os.makedirs(outdirCorrelation)
    outOriginalHistFile = ROOT.TFile(os.path.join(outdirCorrelation, "CorrelationsResults_Original_2D.root"), "RECREATE")
    outHistFile = ROOT.TFile(os.path.join(outdirCorrelation, "CorrelationsResults.root"), "RECREATE")
    outHistFile_BaselineSubtr = ROOT.TFile(os.path.join(outdirCorrelation, "CorrelationsResults_BaselineSubtr.root"), "RECREATE")
    outHistFile_Reflected = ROOT.TFile(os.path.join(outdirCorrelation, "CorrelationsResults_Reflected.root"), "RECREATE")
    outHistFile_Reflected_BaselineSubtr = ROOT.TFile(os.path.join(outdirCorrelation, "CorrelationsResults_Reflected_BaselineSubtr.root"), "RECREATE")
    outFiles = [outOriginalHistFile, outHistFile, outHistFile_BaselineSubtr, outHistFile_Reflected, outHistFile_Reflected_BaselineSubtr]
    for results in all_results:
        subOutdirPtCand = f"PtCandBin_{int(results['task']['ptMin']*10):.0f}_{int(results['task']['ptMax']*10):.0f}"
        subOutdirPtHad = f"PtHadBin_{int(results['task']['ptHadMin']*10):.0f}_{int(results['task']['ptHadMax']*10):.0f}"
        subOutdirInvMass = f"InvMassBin_{int(results['task']['invMassMin']*1000):.0f}_{int(results['task']['invMassMax']*1000):.0f}"
        subOutdirDeltaPhi = f"DeltaPhiBin_{int(results['task']['deltaPhiMin']*1000):.0f}_{int(results['task']['deltaPhiMax']*1000):.0f}"
        for file in outFiles:
            if not file.GetDirectory(subOutdirPtCand):
                file.mkdir(subOutdirPtCand)
            file.cd(subOutdirPtCand)
            if not ROOT.gDirectory.GetDirectory(subOutdirPtHad):
                ROOT.gDirectory.mkdir(subOutdirPtHad)
            ROOT.gDirectory.cd(subOutdirPtHad)
            if not ROOT.gDirectory.GetDirectory(subOutdirInvMass if method == "MassBinning" else subOutdirDeltaPhi):
                ROOT.gDirectory.mkdir(subOutdirInvMass if method == "MassBinning" else subOutdirDeltaPhi)
        # Write histograms
        subOutdir = os.path.join(subOutdirPtCand, subOutdirPtHad, subOutdirInvMass if method == "MassBinning" else subOutdirDeltaPhi)

        # Write main histograms
        # correlationFlow = ROOT.TCanvas("correctionFlow", f"Correlation workflow", 1800, 1200)
        # correlationFlow.Divide(2,3)
        # correlationFlow.cd(1)
        # results['hCorrel_SE_2D'].Draw("lego 2")
        # correlationFlow.cd(2)
        # results['hCorrel_ME_2D'].Draw("lego 2")
        # correlationFlow.cd(3)
        # results['hCorrectedCorrel_2D'].Draw("lego 2")
        # correlationFlow.cd(4)
        # results['hCorrectedCorrHisto'].Draw("hist")
        # correlationFlow.cd(5)
        # results['hNonNormalizedCorrHisto'].Draw("hist")
        if method == "DeltaPhiBinning":
            # correlationFlow.cd(6)
            # results['hCorrectedMassPairs'].Draw("hist")
            # correlationFlowMass = ROOT.TCanvas("massDistribution", f"Corrected mass distribution", 1200, 800)
            # correlationFlowMass.Divide(2,3)
            # correlationFlowMass.cd(1)
            # results['hVecMassPairsVsDeltaEta_2D'][0].Draw("lego 2")
            # correlationFlowMass.cd(2)
            # results['hVecVecCorrectedMassPairsVsDeltaEta'][0][5].Draw("hist")
            # correlationFlowMass.cd(3)
            # results['hVecCorrectedRatioVsDeltaEta'][0].Draw("hist")
            # correlationFlowMass.cd(4)
            # results['hCorrectedRatioVsDeltaEta'].Draw("hist")
            # correlationFlowMass.cd(5)
            # results['hCorrectedMassPairs'].Draw("hist")
            # correlationFlowMass.Update()
            if results['task']['iDeltaPhi'] == 0:
                outHistFile.cd(f"{subOutdirPtCand}/{subOutdirPtHad}")
                results['hCorrectedCorrHisto'].Write()
                results['hCorrel_SE_2D'].Write()
                results['hCorrel_ME_2D'].Write()
                results['hCorrectedCorrel_2D'].Write()
                results['hNonNormalizedCorrHisto'].Write()
            # correlationFlow.Write()
            outHistFile.cd(subOutdir)
            results['hCorrectedMassPairs'].Write()
            results['hCorrectedRatioVsDeltaEta'].Write()
            # correlationFlowMass.Write()
        else:
            outHistFile.cd(subOutdir)
            results['hCorrectedCorrHisto'].Write()
            results['hCorrel_SE_2D'].Write()
            results['hCorrel_ME_2D'].Write()
            results['hCorrectedCorrel_2D'].Write()
            results['hNonNormalizedCorrHisto'].Write()
        # correlationFlow.Update()
        # correlationFlow.Write()

        outOriginalHistFile.cd(subOutdir)
        for i, (hVecOriCorrel_SE_2D, hVecOriCorrel_ME_2D, hVecCorrel_SE_2D, hVecCorrel_ME_2D, hVecCorrel_ME_norm_2D, hVecCorrectedCorrel_2D) in enumerate(zip(
            results['hVecOriginalCorrel_SE_2D'],
            results['hVecOriginalCorrel_ME_2D'],
            results['hVecCorrel_SE_2D'],
            results['hVecCorrel_ME_2D'],
            results['hVecCorrel_ME_norm_2D'],
            results['hVecCorrectedCorrel_2D']
        )):
            hVecOriCorrel_SE_2D.Write(f'hOriginalCorrel_SE_2D_{i}')
            hVecOriCorrel_ME_2D.Write(f'hOriginalCorrel_ME_2D_{i}')
            hVecCorrel_SE_2D.Write(f'hCorrel_SE_2D_{i}')
            hVecCorrel_ME_2D.Write(f'hCorrel_ME_2D_{i}')
            hVecCorrel_ME_norm_2D.Write(f'hCorrel_ME_norm_2D_{i}')
            hVecCorrectedCorrel_2D.Write(f'hCorrectedCorrel_2D_{i}')
            results['hVecMassPairsVsDeltaEta_2D'][i].Write(f'hMassPairsVsDeltaEta_2D_{i}')
            if method == "DeltaPhiBinning":
                for hVecCorrectedMassPairsVsDeltaEta in results['hVecVecCorrectedMassPairsVsDeltaEta']:
                    for hCorrectedMassPairsVsDeltaEta in hVecCorrectedMassPairsVsDeltaEta:
                        hCorrectedMassPairsVsDeltaEta.Write()
                for hCorrectedRatioVsDeltaEta in results['hVecCorrectedRatioVsDeltaEta']:
                    hCorrectedRatioVsDeltaEta.Write()

        outHistFile_BaselineSubtr.cd(subOutdir)
        results['hCorrectedCorrHisto_BaselineSubtr'].Write()

        outHistFile_Reflected.cd(subOutdir)
        results['hCorrectedCorrHisto_Reflected'].Write()

        outHistFile_Reflected_BaselineSubtr.cd(subOutdir)
        results['hCorrectedCorrHisto_Reflected_BaselineSubtr'].Write()
    outOriginalHistFile.Close()
    outHistFile.Close()
    outHistFile_BaselineSubtr.Close()
    outHistFile_Reflected.Close()
    outHistFile_Reflected_BaselineSubtr.Close()

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_dir, '../../', 'utils'))
from utils import check_dir

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract output correlations")
    parser.add_argument("config", nargs="?", default="config_CorrAnalysis_v2_010_negDeta.yaml", help="Configuration file")
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    suffix = config["suffix"]
    outdir = config["outdir"]

    check_dir(os.path.join(outdir, f"CorrelExtract_{suffix}"))
    check_dir(os.path.join(outdir, f"InvMass"))

    ExtractOutputCorrel(args.config)