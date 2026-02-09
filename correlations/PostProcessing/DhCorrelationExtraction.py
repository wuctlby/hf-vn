import ROOT
import math
import sys

class DhCorrelationExtraction:
    # Enums
    kD0toKpi = 0
    kDplusKpipi = 1
    kDsToKKPi = 2
    kDStarD0pi = 3

    kSE = 0
    kME = 1

    kPrompt = 0
    kFD = 1

    kPrimaryPart = 0
    kAllPart = 1

    def __init__(self):
        self.fFileMass = None
        self.fFileSE = None
        self.fFileME = None
        self.fFileFDTemplate = None
        self.fFileFDPromptFrac = None
        self.fFileSecPart = None
        self.fDirMass = None
        self.fDirSE = None
        self.fDirME = None
        self.fDirSecPart = None
        self.fFilePromptMc = None
        self.fFileNonPromptMc = None
        
        self.fCorrHisto_2D_SE = None
        self.fCorrHisto_2D_ME = None
        self.fCorrectedCorrHisto_2D = None
        self.fMassVsPt = None
        self.fCorrectedCorrHisto = None
        self.fCorrectedCorrHisto_BaselineSubtr = None
        self.fCorrectedCorrHisto_Reflected = None
        self.fCorrectedCorrHisto_Reflected_BaselineSubtr = None
        
        self.fDmesonSpecies = self.kDsToKKPi
        self.fDmesonLabel = "Ds"
        self.fNpools = 9
        self.fDeltaEtaMin = -1.
        self.fDeltaEtaMax = 1.
        self.fCorrectPoolsSeparately = True
        self.fSubtractSoftPiME = False
        
        self.fFileNameSE = ""
        self.fFileNameME = ""
        self.fFileSecPartName = ""
        self.fFilePromptMcRecName = ""
        self.fFileNonPromptMcRecName = ""
        self.fDirNameSE = ""
        self.fDirNameME = ""
        self.fDirSecPartName = ""
        self.fMassHistoNameSgn = ""
        self.fMassHistoNameBkg = ""
        self.fMassHistoNameSBs = ""
        self.fSECorrelHistoName = ""
        self.fMECorrelHistoName = ""
        self.fSECorrelSignalRegionName = ""
        self.fSECorrelSidebandsName = ""
        self.fSECorrelSidebandLeftName = ""
        self.fSECorrelSidebandRightName = ""
        self.fMECorrelSignalRegionName = ""
        self.fMECorrelSidebandsName = ""
        self.fMECorrelSidebandLeftName = ""
        self.fMECorrelSidebandRightName = ""
        self.fFileFDTemplateName = ""
        self.fFileFDPromptFracName = ""
        self.fHistoFDTemplatePromptName = ""
        self.fHistoFDTemplateNonPromptName = ""
        self.fHistoFDPromptFracName = ""
        self.fHistoPrimaryPartName = ""
        self.fHistoAllPartName = ""
        
        self.fBkgScaleFactor = 1.
        self.fSgnYieldNorm = 1.
        self.fBkgYield = 1.
        self.fValPhiMEnorm = 1.
        self.fValEtaMEnorm = 1.
        
        self.fRebinAngCorr = False
        self.fRebinFDCorr = False
        self.fRebinSecPart = False
        self.fSidebandDivided = False
        self.fUseSidebLeft = False
        self.fUseSidebRight = False
        
        self.fRebinAxisDeltaEta = 1
        self.fRebinAxisDeltaPhi = 1
        self.fBinPtCand = 0
        self.fBinPtHad = 0
        self.fBinInvMass = 0
        self.fDebug = 0
        self.fFDsubtraction = 0
        self.fSecPartContamination = 0
        self.fCorrBiasBtoD = 0

    def SetDmesonSpecie(self, k):
        if k < 0 or k > 3:
            print("[ERROR] D meson specie not correctly set!")
            return False
        elif k == 0:
            self.fDmesonLabel = "Dzero"
        elif k == 1:
            self.fDmesonLabel = "Dplus"
        elif k == 2:
            self.fDmesonLabel = "Ds"
        else:
            self.fDmesonLabel = "Dstar"
        self.fDmesonSpecies = k
        return True

    # Setters
    def SetInputFilenameMass(self, filenameMass): self.fFileNameMass = filenameMass
    def SetInputFilenameSE(self, filenameSE): self.fFileNameSE = filenameSE
    def SetInputFilenameME(self, filenameME): self.fFileNameME = filenameME
    def SetInputFilenameSecPart(self, filenameSecPart): self.fFileSecPartName = filenameSecPart
    def SetInputFilenameBiasBtoD(self, filenamePromptMcRec, filenameNonPromptMcRec):
        self.fFilePromptMcRecName = filenamePromptMcRec
        self.fFileNonPromptMcRecName = filenameNonPromptMcRec
    def SetDirNameSE(self, dirNameSE): self.fDirNameSE = dirNameSE
    def SetDirNameME(self, dirNameME): self.fDirNameME = dirNameME
    def SetDirNameSecPart(self, dirNameSecPart): self.fDirSecPartName = dirNameSecPart
    def SetMassHistoNameSgn(self, massHistoNameSgn): self.fMassHistoNameSgn = massHistoNameSgn
    def SetMassHistoNameBkg(self, massHistoNameBkg): self.fMassHistoNameBkg = massHistoNameBkg
    def SetMassHistoNameSBs(self, massHistoNameSBs): self.fMassHistoNameSBs = massHistoNameSBs
    def SetCorrelHistoSE(self, correlNameSE): self.fSECorrelHistoName = correlNameSE
    def SetCorrelHistoME(self, correlNameME): self.fMECorrelHistoName = correlNameME
    def SetSECorrelHistoSignalName(self, correlNameSigSE): self.fSECorrelSignalRegionName = correlNameSigSE
    def SetSECorrelHistoSidebandName(self, correlNameSbSE): self.fSECorrelSidebandsName = correlNameSbSE
    def SetSECorrelHistoSidebandLeftName(self, correlNameSbSE): self.fSECorrelSidebandLeftName = correlNameSbSE
    def SetSECorrelHistoSidebandRightName(self, correlNameSbSE): self.fSECorrelSidebandRightName = correlNameSbSE
    def SetMECorrelHistoSignalName(self, correlNameSigME): self.fMECorrelSignalRegionName = correlNameSigME
    def SetMECorrelHistoSidebandName(self, correlNameSbME): self.fMECorrelSidebandsName = correlNameSbME
    def SetMECorrelHistoSidebandLeftName(self, correlNameSbME): self.fMECorrelSidebandLeftName = correlNameSbME
    def SetMECorrelHistoSidebandRightName(self, correlNameSbME): self.fMECorrelSidebandRightName = correlNameSbME
    def SetHistoSecPartName(self, histoPrimaryPartName, histoAllPartName):
        self.fHistoPrimaryPartName = histoPrimaryPartName
        self.fHistoAllPartName = histoAllPartName
    def SetInputFilenameFDTemplate(self, filenameFDTemplate): self.fFileFDTemplateName = filenameFDTemplate
    def SetInputFilenameFDPromptFrac(self, filenameFDPromptFrac): self.fFileFDPromptFracName = filenameFDPromptFrac
    def SetInputHistoNameFDTemplatePrompt(self, hNameFDTemplatePrompt): self.fHistoFDTemplatePromptName = hNameFDTemplatePrompt
    def SetInputHistoNameFDTemplateNonPrompt(self, hNameFDTemplateNonPrompt): self.fHistoFDTemplateNonPromptName = hNameFDTemplateNonPrompt
    def SetInputHistoNameFDPromptFrac(self, hNameFDPromptFrac): self.fHistoFDPromptFracName = hNameFDPromptFrac
    def SetNpools(self, npools): self.fNpools = npools
    def SetCorrectPoolsSeparately(self, usePools): self.fCorrectPoolsSeparately = usePools
    def SetDeltaEtaRange(self, etaLow=-1., etaHigh=1):
        self.fDeltaEtaMin = etaLow
        self.fDeltaEtaMax = etaHigh
    def SetSubtractSoftPiInMEdistr(self, subtractSoftPiME): self.fSubtractSoftPiME = subtractSoftPiME
    def SetBkgScaleFactor(self, scaleFactor): self.fBkgScaleFactor = scaleFactor
    def SetSignalYieldforNorm(self, sgnYield): self.fSgnYieldNorm = sgnYield
    def SetBkgYield(self, bkgYield): self.fBkgYield = bkgYield
    def SetSBYield(self, sbYield): self.fSBYield = sbYield
    def SetRebinOptions(self, rebinAngCorr, rebinFDCorr, rebinSecPart):
        self.fRebinAngCorr = rebinAngCorr
        self.fRebinFDCorr = rebinFDCorr
        self.fRebinSecPart = rebinSecPart
    def SetRebin2DcorrelHisto(self, rebinAxisDeltaEta, rebinAxisDeltaPhi):
        self.fRebinAxisDeltaEta = rebinAxisDeltaEta
        self.fRebinAxisDeltaPhi = rebinAxisDeltaPhi
    def SetBinDeltaPhiEtaForMEnorm(self, valPhiMEnorm, valEtaMEnorm):
        self.fValPhiMEnorm = valPhiMEnorm
        self.fValEtaMEnorm = valEtaMEnorm
    def SetDebugLevel(self, debug): self.fDebug = debug
    def SetFDSubtraction(self, fdSubtraction): self.fFDsubtraction = fdSubtraction
    def SetSecPartContamination(self, secPartContamination): self.fSecPartContamination = secPartContamination
    def SetCorrBiasBtoD(self, corrBiasBtoD): self.fCorrBiasBtoD = corrBiasBtoD
    def SetDividedSidebands(self, dividedSidebands, useSidebLeft, useSidebRight):
        self.fSidebandDivided = dividedSidebands
        self.fUseSidebLeft = useSidebLeft
        self.fUseSidebRight = useSidebRight
    def SetBinCandAndHad(self, binPtCand, binPtHad, binInvMass):
        self.fBinPtCand = binPtCand
        self.fBinPtHad = binPtHad
        self.fBinInvMass = binInvMass

    # Getters
    def GetCorrHisto2D_SE(self): return self.fCorrHisto_2D_SE
    def GetCorrHisto2D_ME(self): return self.fCorrHisto_2D_ME
    def GetInvMassVsPtHisto(self): return self.fMassVsPt
    def GetCorrectedCorrHisto2D(self): return self.fCorrectedCorrHisto_2D
    def GetCorrectedCorrHisto(self): return self.fCorrectedCorrHisto
    def GetCorrectedCorrHisto_BaselineSubtr(self): return self.fCorrectedCorrHisto_BaselineSubtr
    def GetCorrectedCorrHisto_Reflected(self): return self.fCorrectedCorrHisto_Reflected
    def GetCorrectedCorrHisto_Reflected_BaselineSubtr(self): return self.fCorrectedCorrHisto_Reflected_BaselineSubtr

    def ReadInputSEandME(self):
        self.fFileSE = ROOT.TFile.Open(self.fFileNameSE)
        if not self.fFileSE:
            print(f"[ERROR] File {self.fFileNameSE} cannot be opened! check your file path!")
            return False

        self.fFileME = ROOT.TFile.Open(self.fFileNameME)
        if not self.fFileME:
            print(f"[ERROR] File {self.fFileNameME} cannot be opened! check your file path!")
            return False

        self.fDirSE = self.fFileSE.Get(self.fDirNameSE)
        self.fDirME = self.fFileME.Get(self.fDirNameME)

        print("=====================")
        print("Read inputs SE and ME")
        print(f"TFile SE    = {self.fFileNameSE}")
        print(f"TFile ME    = {self.fFileNameME}")
        print(f"TDir SE    = {self.fDirNameSE}")
        print(f"TDir ME    = {self.fDirNameME}")
        print("=====================")
        print(" ")
        return True

    def ReadInputInvMass(self):
        self.fFileMass = ROOT.TFile.Open(self.fFileNameMass)
        if not self.fFileMass:
            print(f"[ERROR] File {self.fFileNameMass} cannot be opened! check your file path!")
            return False

        print("=====================")
        print("Read inputs inv. mass")
        print(f"TFile Mass    = {self.fFileNameMass}")
        print(f"Histo Mass name    = {self.fMassHistoNameSgn}")
        print("=====================")
        print(" ")
        return True

    def ReadInputFDSubtr(self):
        self.fFileFDTemplate = ROOT.TFile.Open(self.fFileFDTemplateName)
        self.fFileFDPromptFrac = ROOT.TFile.Open(self.fFileFDPromptFracName)
        if not self.fFileFDTemplate:
            print(f"[ERROR] File {self.fFileFDTemplateName} cannot be opened! check your file path!")
            return False
        if not self.fFileFDPromptFrac:
            print(f"[ERROR] File {self.fFileFDPromptFracName} cannot be opened! check your file path!")
            return False

        print("=====================")
        print("Read inputs FD template")
        print(f"TFile FD template    = {self.fFileFDTemplateName}")
        print(f"TFile FD Prompt Frac    = {self.fFileFDPromptFracName}")
        print(f"Histo FD template Prompt    = {self.fHistoFDTemplatePromptName}")
        print(f"Histo FD template Non Prompt     = {self.fHistoFDTemplateNonPromptName}")
        print(f"Histo FD Prompt Frac     = {self.fHistoFDPromptFracName}")
        print("=====================")
        print(" ")
        return True

    def ReadInputSecondaryPartContamination(self):
        self.fFileSecPart = ROOT.TFile.Open(self.fFileSecPartName)
        if not self.fFileSecPart:
            print(f"[ERROR] File {self.fFileSecPartName} cannot be opened! check your file path!")
            return False

        self.fDirSecPart = self.fFileSecPart.Get(self.fDirSecPartName)

        if not self.fDirSecPart:
            print(f"[ERROR] Directory {self.fDirSecPartName} cannot be opened! check your file path!")
            return False

        print("=====================")
        print("Read inputs SE and ME")
        print(f"TFile Sec. part.    = {self.fFileSecPartName}")
        print(f"TDir Sec. part.    = {self.fDirSecPartName}")
        print("=====================")
        print(" ")
        return True

    def GetCorrelHisto(self, SEorME, pool, PtCandMin, PtCandMax, PtHadMin, PtHadMax, InvMassMin, InvMassMax):
        h2D = ROOT.TH2D()
        hSparse = None
        if SEorME == self.kSE:
            hSparse = self.fDirSE.Get(self.fSECorrelHistoName)
        else:
            hSparse = self.fDirME.Get(self.fMECorrelHistoName)

        if not hSparse:
            print("[ERROR] hSparse is null! Check that the object name exists in the directory and the file is open.")
            raise RuntimeError("hSparse is null")

        binExtPtCandMin = hSparse.GetAxis(1).FindBin(PtCandMin + 0.0001)
        binExtPtCandMax = hSparse.GetAxis(1).FindBin(PtCandMax - 0.0001)
        binExtPtHadMin = hSparse.GetAxis(2).FindBin(PtHadMin + 0.0001)
        binExtPtHadMax = hSparse.GetAxis(2).FindBin(PtHadMax - 0.0001)
        binExtInvMassMin = hSparse.GetAxis(5).FindBin(InvMassMin + 0.0001)
        binExtInvMassMax = hSparse.GetAxis(5).FindBin(InvMassMax - 0.0001)
        
        binExtPoolMin = 0
        binExtPoolMax = 0
        if self.fCorrectPoolsSeparately:
            binExtPoolMin = hSparse.GetAxis(0).FindBin(pool + 0.01)
            binExtPoolMax = hSparse.GetAxis(0).FindBin(pool + 0.99)
        else:
            binExtPoolMin = 1
            binExtPoolMax = hSparse.GetAxis(0).GetNbins()
            
        binExtEtaMin = hSparse.GetAxis(3).FindBin(self.fDeltaEtaMin + 0.0001)
        binExtEtaMax = hSparse.GetAxis(3).FindBin(self.fDeltaEtaMax - 0.0001)
        if binExtEtaMax > hSparse.GetAxis(3).GetNbins():
            binExtEtaMax = hSparse.GetAxis(3).GetNbins()
        if binExtEtaMin < 1:
            binExtEtaMin = 1
            
        hSparse.GetAxis(0).SetRange(binExtPoolMin, binExtPoolMax)
        hSparse.GetAxis(1).SetRange(binExtPtCandMin, binExtPtCandMax)
        hSparse.GetAxis(2).SetRange(binExtPtHadMin, binExtPtHadMax)
        hSparse.GetAxis(3).SetRange(binExtEtaMin, binExtEtaMax)
        hSparse.GetAxis(5).SetRange(binExtInvMassMin, binExtInvMassMax)
        
        h2D = hSparse.Projection(4, 3) # axis4: deltaPhi, axis3: deltaEta
        if SEorME == self.kSE:
            h2D.SetName(f"hCorr_SE_2D_PtCandBin{binExtPtCandMin}_PtHadBin{binExtPtHadMin}_InvMassBin{binExtInvMassMin}_iPool{pool}")
        else:
            h2D.SetName(f"hCorr_ME_2D_PtCandBin{binExtPtCandMin}_PtHadBin{binExtPtHadMin}_InvMassBin{binExtInvMassMin}_iPool{pool}")
            
        return h2D

    def NormalizeMEplot(self, histoME, histoMEsoftPi):
        bin0phi = histoME.GetYaxis().FindBin(0.)
        bin0eta = histoME.GetXaxis().FindBin(0.)

        factorNorm = 0
        for i_n in range(-1, 1):
            factorNorm += histoME.GetBinContent(bin0eta, bin0phi + i_n)
        for i_n in range(-1, 1):
            factorNorm += histoME.GetBinContent(bin0eta - 1, bin0phi + i_n)
        factorNorm /= 4.

        if factorNorm == 0:
            bin0phi = histoME.GetYaxis().FindBin(self.fValPhiMEnorm)
            bin0eta = histoME.GetXaxis().FindBin(self.fValEtaMEnorm)
            for i_n in range(-1, 1):
                factorNorm += histoME.GetBinContent(bin0eta, bin0phi + i_n)
            for i_n in range(-1, 1):
                factorNorm += histoME.GetBinContent(bin0eta - 1, bin0phi + i_n)
            factorNorm /= 4.

        print(f"bin 0 phi: {bin0phi}")
        print(f"bin 0 eta: {bin0eta}")
        print(f"Factor norm. ME: {factorNorm}")
        print(f"Bin content (0,0) ME: {histoME.GetBinContent(bin0eta, bin0phi)}")

        if self.fSubtractSoftPiME and histoMEsoftPi:
            histoME.Add(histoMEsoftPi, -1)

        if factorNorm != 0:
            histoME.Scale(1. / factorNorm)

    def GetFDTemplateHisto(self, PromptOrFD, PtCandMin, PtCandMax, PtHadMin, PtHadMax):
        h2D = None
        if PromptOrFD == self.kPrompt:
            h2D = self.fFileFDTemplate.Get(f"{self.fHistoFDTemplatePromptName}{PtCandMin:.0f}_{PtCandMax:.0f}_ptassoc{PtHadMin:.0f}_{PtHadMax:.0f}")
        else:
            h2D = self.fFileFDTemplate.Get(f"{self.fHistoFDTemplateNonPromptName}{PtCandMin:.0f}_{PtCandMax:.0f}_ptassoc{PtHadMin:.0f}_{PtHadMax:.0f}")

        if not h2D:
             print(f"[ERROR] FD Template histogram not found for {PtCandMin}-{PtCandMax}, {PtHadMin}-{PtHadMax}")
             return ROOT.TH2D()

        binExtEtaMin = h2D.GetXaxis().FindBin(self.fDeltaEtaMin + 0.000001)
        binExtEtaMax = h2D.GetXaxis().FindBin(self.fDeltaEtaMax - 0.000001)
        if binExtEtaMax > h2D.GetXaxis().GetNbins():
            binExtEtaMax = h2D.GetXaxis().GetNbins()
        if binExtEtaMin < 1:
            binExtEtaMin = 1

        h2D.GetXaxis().SetRange(binExtEtaMin, binExtEtaMax)
        if PromptOrFD == self.kPrompt:
            h2D.SetName(f"hFDTemplatePrompt_2D_PtCand{PtCandMin:.0f}to{PtCandMax:.0f}_PtHad{PtHadMin:.0f}to{PtHadMax:.0f}")
        else:
            h2D.SetName(f"hFDTemplateNonPrompt_2D_PtCand{PtCandMin:.0f}to{PtCandMax:.0f}_PtHad{PtHadMin:.0f}to{PtHadMax:.0f}")
        h2D.GetYaxis().SetTitle("#Delta#phi (rad)")
        h2D.GetXaxis().SetTitle("#Delta#eta")

        return h2D

    def GetFDPromptFrac(self, PtCandMin, PtCandMax, PtHadMin, PtHadMax):
        h1D = self.fFileFDPromptFrac.Get(self.fHistoFDPromptFracName)
        binPtCandMin = h1D.GetXaxis().FindBin(PtCandMin + 0.01)
        binPtCandMax = h1D.GetXaxis().FindBin(PtCandMax - 0.01)
        PromptFraction = 0.
        if binPtCandMin == binPtCandMax:
            PromptFraction = h1D.GetBinContent(binPtCandMin)
        else:
            print("[ERROR] Different bin obtained from PtCandMin and PtCandMax")
            return 0.
        return PromptFraction

    def GetCorrelHistoSecondaryPart(self, PartType, PtCandMin, PtCandMax, PtHadMin, PtHadMax):
        hSparse = None
        if PartType == self.kPrimaryPart:
            hSparse = self.fDirSecPart.Get(self.fHistoPrimaryPartName)
        else:
            hSparse = self.fDirSecPart.Get(self.fHistoAllPartName)
            
        binExtPtCandMin = hSparse.GetAxis(2).FindBin(PtCandMin + 0.01)
        binExtPtCandMax = hSparse.GetAxis(2).FindBin(PtCandMax - 0.01)
        binExtPtHadMin = hSparse.GetAxis(3).FindBin(PtHadMin + 0.01)
        binExtPtHadMax = hSparse.GetAxis(3).FindBin(PtHadMax - 0.01)
        
        binExtPoolMin = 0
        binExtPoolMax = 0
        if PartType == self.kAllPart:
            binExtPoolMin = 1
            binExtPoolMax = hSparse.GetAxis(4).GetNbins()
            
        binExtEtaMin = hSparse.GetAxis(1).FindBin(self.fDeltaEtaMin + 0.0001)
        binExtEtaMax = hSparse.GetAxis(1).FindBin(self.fDeltaEtaMax - 0.0001)
        if binExtEtaMax > hSparse.GetAxis(1).GetNbins():
            binExtEtaMax = hSparse.GetAxis(1).GetNbins()
        if binExtEtaMin < 1:
            binExtEtaMin = 1
            
        hSparse.GetAxis(1).SetRange(binExtEtaMin, binExtEtaMax)
        hSparse.GetAxis(2).SetRange(binExtPtCandMin, binExtPtCandMax)
        hSparse.GetAxis(3).SetRange(binExtPtHadMin, binExtPtHadMax)
        if PartType == self.kAllPart:
            hSparse.GetAxis(4).SetRange(binExtPoolMin, binExtPoolMax)
            
        h1D = hSparse.Projection(0)
        if PartType == self.kPrimaryPart:
            h1D.SetName(f"hPrimaryPartCorr_PtD{PtCandMin:.0f}to{PtCandMax:.0f}_PtHad{PtHadMin:.0f}to{PtHadMax:.0f}")
        else:
            h1D.SetName(f"hAllPartCorr_PtD{PtCandMin:.0f}to{PtCandMax:.0f}_PtHad{PtHadMin:.0f}to{PtHadMax:.0f}")
            
        return h1D

    def CalculateBaseline(self, histo, totalRange, reflected):
        baseline = 0.
        nBinsPhi = histo.GetNbinsX()
        binPhiHalf = nBinsPhi // 2
        binPhiHalfMinus1 = nBinsPhi // 2 - 1
        binPhiHalfPlus1 = nBinsPhi // 2 + 1
        binPhiHalfPlus2 = nBinsPhi // 2 + 2

        # Helper to get content and error squared
        def get_val_err2(bin_idx):
            c = histo.GetBinContent(bin_idx)
            e = histo.GetBinError(bin_idx)
            if e == 0: return 0, 1e-9 # Avoid div by zero
            return c, e*e

        if totalRange:
            if nBinsPhi >= 32:
                bins = [1, 2, binPhiHalfMinus1, binPhiHalf, binPhiHalfPlus1, binPhiHalfPlus2, nBinsPhi-1, nBinsPhi]
            else:
                bins = [1, binPhiHalf, binPhiHalfPlus1, nBinsPhi]
        else:
            if reflected:
                bins = [binPhiHalfMinus1, binPhiHalf, binPhiHalfPlus1, binPhiHalfPlus2]
            else:
                if nBinsPhi >= 32:
                    bins = [binPhiHalfMinus1, binPhiHalf, binPhiHalfPlus1, binPhiHalfPlus2]
                else:
                    bins = [binPhiHalf, binPhiHalfPlus1]

        numerator = 0.
        denominator = 0.
        
        for b in bins:
            val, err2 = get_val_err2(b)
            weight = 1. / err2
            numerator += val * weight
            denominator += weight
            
        if denominator > 0:
            baseline = numerator / denominator
            
        return baseline

    def CalculateBaselineError(self, histo, totalRange, reflected):
        errBaseline = 0.
        nBinsPhi = histo.GetNbinsX()
        binPhiHalf = nBinsPhi // 2
        binPhiHalfMinus1 = nBinsPhi // 2 - 1
        binPhiHalfPlus1 = nBinsPhi // 2 + 1
        binPhiHalfPlus2 = nBinsPhi // 2 + 2

        def get_err2(bin_idx):
            e = histo.GetBinError(bin_idx)
            if e == 0: return 1e-9
            return e*e

        if totalRange:
            if nBinsPhi >= 32:
                bins = [1, 2, binPhiHalfMinus1, binPhiHalf, binPhiHalfPlus1, binPhiHalfPlus2, nBinsPhi-1, nBinsPhi]
            else:
                bins = [1, binPhiHalf, binPhiHalfPlus1, nBinsPhi]
        else:
            if reflected:
                bins = [binPhiHalfMinus1, binPhiHalf, binPhiHalfPlus1, binPhiHalfPlus2]
            else:
                if nBinsPhi >= 32:
                    bins = [binPhiHalfMinus1, binPhiHalf, binPhiHalfPlus1, binPhiHalfPlus2]
                else:
                    bins = [binPhiHalf, binPhiHalfPlus1]

        sum_inv_err2 = 0.
        for b in bins:
            sum_inv_err2 += 1. / get_err2(b)
            
        if sum_inv_err2 > 0:
            errBaseline = 1. / math.sqrt(sum_inv_err2)
            
        return errBaseline

    def ReflectCorrHistogram(self, histo):
        nBinsPhi = histo.GetNbinsX()
        nBinsPhiRefl = nBinsPhi // 2
        bin0Phi = nBinsPhi // 4 + 1
        binPiPhi = 3 * nBinsPhi // 4

        h1D = ROOT.TH1D("h1D_Reflected", "", nBinsPhiRefl, 0., math.pi)
        h1D.Sumw2()

        for iBin in range(nBinsPhiRefl // 2):
            val1 = histo.GetBinContent(bin0Phi - iBin - 1)
            val2 = histo.GetBinContent(bin0Phi + iBin)
            reflectedContent = (val1 + val2) / 2.
            
            err1 = histo.GetBinError(bin0Phi - iBin - 1)
            err2 = histo.GetBinError(bin0Phi + iBin)
            reflectedContentError = 0.5 * math.sqrt(err1**2 + err2**2)
            
            h1D.SetBinContent(iBin + 1, reflectedContent)
            h1D.SetBinError(iBin + 1, reflectedContentError)

        for iBin in range(nBinsPhiRefl // 2, nBinsPhiRefl):
            val1 = histo.GetBinContent(bin0Phi + iBin)
            # binPiPhi + 2 * bin0Phi - iBin - 2  <-- from C++
            idx2 = binPiPhi + 2 * bin0Phi - iBin - 2
            val2 = histo.GetBinContent(idx2)
            reflectedContent = (val1 + val2) / 2.
            
            err1 = histo.GetBinError(bin0Phi + iBin)
            err2 = histo.GetBinError(idx2)
            reflectedContentError = 0.5 * math.sqrt(err1**2 + err2**2)
            
            h1D.SetBinContent(iBin + 1, reflectedContent)
            h1D.SetBinError(iBin + 1, reflectedContentError)
            
        return h1D

    def SetTH1HistoStyle(self, histo, hTitle, hXaxisTitle, hYaxisTitle, markerStyle, markerColor, markerSize, lineColor, lineWidth, hTitleXaxisOffset=1.0, hTitleYaxisOffset=1.0, hTitleXaxisSize=0.05, hTitleYaxisSize=0.05, hLabelXaxisSize=0.05, hLabelYaxisSize=0.05, centerXaxisTitle=False, centerYaxisTitle=False):
        histo.SetTitle(hTitle)
        histo.GetXaxis().SetTitle(hXaxisTitle)
        histo.GetYaxis().SetTitle(hYaxisTitle)
        histo.SetMarkerStyle(markerStyle)
        histo.SetMarkerColor(markerColor)
        histo.SetMarkerSize(markerSize)
        histo.SetLineColor(lineColor)
        histo.SetLineWidth(lineWidth)
        histo.GetXaxis().SetTitleOffset(hTitleXaxisOffset)
        histo.GetYaxis().SetTitleOffset(hTitleYaxisOffset)
        histo.GetXaxis().SetTitleSize(hTitleXaxisSize)
        histo.GetYaxis().SetTitleSize(hTitleYaxisSize)
        histo.GetXaxis().SetLabelSize(hLabelXaxisSize)
        histo.GetYaxis().SetLabelSize(hLabelYaxisSize)
        histo.GetXaxis().CenterTitle(centerXaxisTitle)
        histo.GetYaxis().CenterTitle(centerYaxisTitle)

    def EvaluateMCClosModulations(self, PtCandMin, PtCandMax, PtHadMin, PtHadMax):
        # Placeholder for EvaluateMCClosModulations as it requires specific MC files and logic
        # that might need more context or files not fully described.
        # Implementing basic structure based on C++ code.
        
        hModul = ROOT.TH1D()
        
        self.fFilePromptMc = ROOT.TFile.Open(self.fFilePromptMcRecName)
        self.fFileNonPromptMc = ROOT.TFile.Open(self.fFileNonPromptMcRecName)
        
        if not self.fFilePromptMc:
            print("[ERROR] File prompt MC rec cannot be opened!")
        if not self.fFileNonPromptMc:
            print("[ERROR] File non-prompt MC rec cannot be opened!")
            
        hRecPrompt = self.fFilePromptMc.Get(f"h1D_Rec_iPtD{self.fBinPtCand}_iPtAssoc{self.fBinPtHad}")
        hRecNonPrompt = self.fFileNonPromptMc.Get(f"h1D_Rec_iPtD{self.fBinPtCand}_iPtAssoc{self.fBinPtHad}")
        hGenPrompt = self.fFilePromptMc.Get(f"h1D_Gen_iPtD{self.fBinPtCand}_iPtAssoc{self.fBinPtHad}")
        hGenNonPrompt = self.fFileNonPromptMc.Get(f"h1D_Gen_iPtD{self.fBinPtCand}_iPtAssoc{self.fBinPtHad}")
        
        if not hRecPrompt or not hRecNonPrompt or not hGenPrompt or not hGenNonPrompt:
             print("[ERROR] MC histograms not found")
             return hModul

        hRecNonPrompt.Sumw2()
        hGenPrompt.Sumw2()
        hGenNonPrompt.Sumw2()
        
        hRatioNonPrompt = hRecNonPrompt.Clone("hRatioNonPrompt")
        hRatioNonPrompt.Sumw2()
        hRatioNonPrompt.Divide(hRecNonPrompt, hGenNonPrompt, 1., 1., "B")
        hModul = hRatioNonPrompt.Clone("hModul")
        
        funFit = ROOT.TF1("funFit", "[0]", math.pi * 3. / 8., math.pi * 3 / 2)
        hRatioNonPrompt.Fit(funFit, "R")
        fitVal = funFit.GetParameter(0)
        
        FPrompt = self.GetFDPromptFrac(PtCandMin, PtCandMax, PtHadMin, PtHadMax)
        
        for iBin in range(hRatioNonPrompt.GetNbinsX()):
            # Logic from C++: if (iBin > 1 && iBin < 13) -> indices 2 to 12 (1-based)
            # In Python 0-based: indices 1 to 11?
            # C++: iBin starts at 0. hRatioNonPrompt->GetNbinsX()
            # C++ loop: for (int iBin = 0; iBin < hRatioNonPrompt->GetNbinsX(); iBin++)
            # C++ check: if (iBin > 1 && iBin < 13)
            # This means iBin = 2, 3, ..., 12.
            # Bin index in ROOT is iBin + 1. So bins 3 to 13.
            
            if iBin > 1 and iBin < 13:
                bin_idx = iBin + 1
                recoKineVal = hRatioNonPrompt.GetBinContent(bin_idx) - (fitVal - 1)
                
                recP = hRecPrompt.GetBinContent(bin_idx)
                recNP = hRecNonPrompt.GetBinContent(bin_idx)
                denom = (recP * FPrompt + recNP * (1 - FPrompt))
                
                if denom != 0:
                    relAmplC = recP / denom
                    relAmplB = recNP / denom
                    
                    if recoKineVal != 0:
                        modul = relAmplC * FPrompt + relAmplB * (1 - FPrompt) / recoKineVal
                        hModul.SetBinContent(bin_idx, modul)
                        hModul.SetBinError(bin_idx, 0.)
            else:
                hModul.SetBinContent(iBin + 1, 1.)
                hModul.SetBinError(iBin + 1, 0.)
                
        return hModul

    def ExtractCorrelations(self, PtCandMin, PtCandMax, PtHadMin, PtHadMax, InvMassMin, InvMassMax, codeName):
        if self.fSubtractSoftPiME:
            print("[INFO] Fake softPi subtraction in ME via extraction code is enabled!")

        if not self.fCorrectPoolsSeparately:
            self.fNpools = 1

        hSE_Sign = [None] * self.fNpools
        hME_Sign = [None] * self.fNpools
        hME_Original = [None] * self.fNpools
        hME_Sign_SoftPi = [None] * self.fNpools
        hCorr_Sign = [None] * self.fNpools
        
        h2D_Sign = None
        h2D_SE = None
        h2D_ME = None
        
        for iPool in range(self.fNpools):
            print("Start getting correlation histograms")
            hSE_Sign[iPool] = self.GetCorrelHisto(self.kSE, iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax, InvMassMin, InvMassMax)
            print("Got SE histogram region")
            hME_Sign[iPool] = self.GetCorrelHisto(self.kME, iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax, InvMassMin, InvMassMax)
            print("Got ME histogram region")
            
            hSE_Sign[iPool].Sumw2()
            hME_Sign[iPool].Sumw2()
            
            if self.fRebinAngCorr:
                hSE_Sign[iPool].Rebin2D(self.fRebinAxisDeltaEta, self.fRebinAxisDeltaPhi)
                hME_Sign[iPool].Rebin2D(self.fRebinAxisDeltaEta, self.fRebinAxisDeltaPhi)
                print("SE and ME histograms rebinned")
                
            hME_Original[iPool] = hME_Sign[iPool].Clone(f"hME_Original_Pool{iPool}")
            
            self.NormalizeMEplot(hME_Sign[iPool], hME_Sign_SoftPi[iPool])
            
            hCorr_Sign[iPool] = hSE_Sign[iPool].Clone(f"hCorr_Sign_Pool{iPool}")
            hCorr_Sign[iPool].Sumw2()
            hCorr_Sign[iPool].Divide(hME_Sign[iPool])
            
            N_SEsign = 0
            N_sign = 0
            for i in range(1, hCorr_Sign[iPool].GetXaxis().GetNbins() + 1):
                for j in range(1, hCorr_Sign[iPool].GetYaxis().GetNbins() + 1):
                    N_SEsign += hSE_Sign[iPool].GetBinContent(i, j)
                    N_sign += hCorr_Sign[iPool].GetBinContent(i, j)
            
            hSE_Sign[iPool].SetEntries(N_SEsign)
            hCorr_Sign[iPool].SetEntries(N_sign)
            
            if iPool == 0:
                h2D_Sign = hCorr_Sign[0].Clone("h2D_Sign")
                h2D_SE = hSE_Sign[0].Clone("h2D_SE")
                h2D_ME = hME_Original[0].Clone("h2D_ME")
                h2D_Sign.Sumw2()
                h2D_SE.Sumw2()
                h2D_ME.Sumw2()
            else:
                h2D_Sign.Add(hCorr_Sign[iPool])
                h2D_SE.Add(hSE_Sign[iPool])
                h2D_ME.Add(hME_Original[iPool])

        self.fCorrHisto_2D_SE = h2D_SE.Clone(f"hCorrSE_PtCand{PtCandMin:.0f}to{PtCandMax:.0f}_PtAssoc{PtHadMin:.0f}to{PtHadMax:.0f}_InvMassBin{self.fBinInvMass}")
        self.fCorrHisto_2D_ME = h2D_ME.Clone(f"hCorrME_PtCand{PtCandMin:.0f}to{PtCandMax:.0f}_PtAssoc{PtHadMin:.0f}to{PtHadMax:.0f}_InvMassBin{self.fBinInvMass}")
        self.fCorrectedCorrHisto_2D = h2D_Sign.Clone(f"hCorrSE_MEcorrected_PtCand{PtCandMin:.0f}to{PtCandMax:.0f}_PtAssoc{PtHadMin:.0f}to{PtHadMax:.0f}_InvMassBin{self.fBinInvMass}")

        h2D_FDTemplatePrompt = None
        h2D_FDTemplateNonPrompt = None
        FDPromptFrac = 0.
        
        if self.fFDsubtraction:
            h2D_FDTemplatePrompt = self.GetFDTemplateHisto(self.kPrompt, PtCandMin, PtCandMax, PtHadMin, PtHadMax)
            h2D_FDTemplateNonPrompt = self.GetFDTemplateHisto(self.kFD, PtCandMin, PtCandMax, PtHadMin, PtHadMax)
            FDPromptFrac = self.GetFDPromptFrac(PtCandMin, PtCandMax, PtHadMin, PtHadMax)
            
            if self.fRebinFDCorr:
                h2D_FDTemplatePrompt.Rebin2D(self.fRebinAxisDeltaEta, self.fRebinAxisDeltaPhi)
                h2D_FDTemplateNonPrompt.Rebin2D(self.fRebinAxisDeltaEta, self.fRebinAxisDeltaPhi)

        h1D_Sign = h2D_Sign.ProjectionY("h1D_Sign")
        h1D_Subtr = h1D_Sign.Clone("h1D_Subtr")
        h1D_Subtr.Sumw2()
        
        h1D_SubtrNorm = h1D_Subtr.Clone("h1D_SubtrNorm")
        h1D_SubtrNorm.Sumw2()
        
        hModul = None
        if self.fCorrBiasBtoD:
            hModul = self.EvaluateMCClosModulations(PtCandMin, PtCandMax, PtHadMin, PtHadMax)
            h1D_SubtrNorm.Multiply(hModul)
            
        h1D_PrimaryPartCorr = None
        h1D_AllPartCorr = None
        h1D_SecPartFrac = None
        h1D_SubtrNorm_SecPart = None
        
        if self.fSecPartContamination:
            h1D_PrimaryPartCorr = self.GetCorrelHistoSecondaryPart(self.kPrimaryPart, PtCandMin, PtCandMax, PtHadMin, PtHadMax)
            h1D_AllPartCorr = self.GetCorrelHistoSecondaryPart(self.kAllPart, PtCandMin, PtCandMax, PtHadMin, PtHadMax)
            h1D_PrimaryPartCorr.Sumw2()
            h1D_AllPartCorr.Sumw2()
            
            if self.fRebinSecPart:
                h1D_PrimaryPartCorr.RebinX(self.fRebinAxisDeltaPhi)
                h1D_AllPartCorr.RebinX(self.fRebinAxisDeltaPhi)
                print("Secondary particle histogram rebinned")
                
            h1D_SecPartFrac = h1D_PrimaryPartCorr.Clone(f"hCorrRatio_PtD{PtCandMin:.0f}to{PtCandMax:.0f}_PtHad{PtHadMin:.0f}to{PtHadMax:.0f}")
            h1D_SecPartFrac.Sumw2()
            h1D_SecPartFrac.Divide(h1D_PrimaryPartCorr, h1D_AllPartCorr, 1., 1., "B")
            
            h1D_SubtrNorm_SecPart = h1D_SubtrNorm.Clone("h1D_SubtrNorm_SecPart")
            h1D_SubtrNorm_SecPart.Sumw2()
            if h1D_SubtrNorm_SecPart.GetNbinsX() != h1D_SecPartFrac.GetNbinsX():
                print("[ERROR]: nBinsPhi different between h1D_SubtrNorm and h1D_SecPartFrac")
                return False
            h1D_SubtrNorm_SecPart.Multiply(h1D_SecPartFrac)

        h1D_FDTemplatePrompt = None
        h1D_FDTemplateNonPrompt = None
        h1D_TemplateTotal = None
        h1D_SubtrFDNorm = None
        
        if self.fFDsubtraction:
            h1D_FDTemplatePrompt = h2D_FDTemplatePrompt.ProjectionY("h1D_FDTemplatePrompt")
            h1D_FDTemplateNonPrompt = h2D_FDTemplateNonPrompt.ProjectionY("h1D_FDTemplateNonPrompt")
            
            h1D_FDTemplatePrompt.Scale(1. / h1D_FDTemplatePrompt.GetXaxis().GetBinWidth(1))
            h1D_FDTemplateNonPrompt.Scale(1. / h1D_FDTemplateNonPrompt.GetXaxis().GetBinWidth(1))
            
            h1D_TemplateTotal = h1D_FDTemplatePrompt.Clone("h1D_TemplateTotal")
            h1D_TemplateTotal.Sumw2()
            h1D_TemplateTotal.Scale(FDPromptFrac)
            h1D_TemplateTotal.Add(h1D_FDTemplateNonPrompt, 1 - FDPromptFrac)
            
            BaselineFD = self.CalculateBaseline(h1D_TemplateTotal, True, False)
            BaselineData = 0.
            if self.fSecPartContamination:
                BaselineData = self.CalculateBaseline(h1D_SubtrNorm_SecPart, True, False)
            else:
                BaselineData = self.CalculateBaseline(h1D_SubtrNorm, True, False)
                
            print("=====================")
            print(f"Baseline FD: {BaselineFD}")
            print(f"Baseline Data: {BaselineData}")
            print("=====================")
            
            Baselinediff = BaselineData - BaselineFD
            hBaselineDiff = h1D_FDTemplateNonPrompt.Clone("hBaselineDiff")
            for iBin in range(hBaselineDiff.GetNbinsX()):
                hBaselineDiff.SetBinContent(iBin + 1, Baselinediff)
                
            h1D_FDTemplateNonPrompt.Add(hBaselineDiff)
            h1D_TemplateTotal.Add(hBaselineDiff)
            
            if self.fSecPartContamination:
                h1D_SubtrFDNorm = h1D_SubtrNorm_SecPart.Clone("h1D_SubtrFDNorm")
            else:
                h1D_SubtrFDNorm = h1D_SubtrNorm.Clone("h1D_SubtrFDNorm")
                
            h1D_FDTemplateNonPrompt.Scale(1 - FDPromptFrac)
            h1D_SubtrFDNorm.Add(h1D_FDTemplateNonPrompt, -1)
            h1D_SubtrFDNorm.Scale(1. / FDPromptFrac)

        if self.fFDsubtraction:
            self.fCorrectedCorrHisto = h1D_SubtrFDNorm.Clone(f"hCorrectedCorr_PtCand{PtCandMin:.0f}to{PtCandMax:.0f}_PtAssoc{PtHadMin:.0f}to{PtHadMax:.0f}_InvMassBin{self.fBinInvMass}")
        elif self.fSecPartContamination:
            self.fCorrectedCorrHisto = h1D_SubtrNorm_SecPart.Clone(f"hCorrectedCorr_PtCand{PtCandMin:.0f}to{PtCandMax:.0f}_PtAssoc{PtHadMin:.0f}to{PtHadMax:.0f}_InvMassBin{self.fBinInvMass}")
        else:
            self.fCorrectedCorrHisto = h1D_SubtrNorm.Clone(f"hCorrectedCorr_PtCand{PtCandMin:.0f}to{PtCandMax:.0f}_PtAssoc{PtHadMin:.0f}to{PtHadMax:.0f}_InvMassBin{self.fBinInvMass}")

        print("Analysis steps completed - baseline subtraction missing")

        # Baseline subtraction
        BaselineData = 0.
        BaselineDataErr = 0.
        hBaseline = h1D_SubtrNorm.Clone("hBaseline")
        hBaseline.Sumw2()
        
        target_histo = None
        if self.fFDsubtraction:
            target_histo = h1D_SubtrFDNorm
        elif self.fSecPartContamination:
            target_histo = h1D_SubtrNorm_SecPart
        else:
            target_histo = h1D_SubtrNorm
            
        BaselineData = self.CalculateBaseline(target_histo, True, False)
        BaselineDataErr = self.CalculateBaselineError(target_histo, True, False)
        
        for iBin in range(hBaseline.GetNbinsX()):
            hBaseline.SetBinContent(iBin + 1, BaselineData)
            hBaseline.SetBinError(iBin + 1, BaselineDataErr)
            
        self.fCorrectedCorrHisto_BaselineSubtr = target_histo.Clone("h1D_BaselineSubtr")
        self.fCorrectedCorrHisto_BaselineSubtr.Add(hBaseline, -1.)
        self.fCorrectedCorrHisto_BaselineSubtr.SetName(f"hCorrectedCorrBaselineSubtr_PtCand{PtCandMin:.0f}to{PtCandMax:.0f}_PtAssoc{PtHadMin:.0f}to{PtHadMax:.0f}_InvMassBin{self.fBinInvMass}")

        # Reflected histograms
        h1D_ReflCorr = self.ReflectCorrHistogram(target_histo)
        
        hBaseline_Refl = h1D_ReflCorr.Clone("hBaseline_Refl")
        hBaseline_Refl.Sumw2()
        BaselineData = self.CalculateBaseline(h1D_ReflCorr, False, True)
        BaselineDataErr = self.CalculateBaselineError(h1D_ReflCorr, False, True)
        
        for iBin in range(hBaseline_Refl.GetNbinsX()):
            hBaseline_Refl.SetBinContent(iBin + 1, BaselineData)
            hBaseline_Refl.SetBinError(iBin + 1, BaselineDataErr)
            
        h1D_ReflCorr_BaselineSubtr = h1D_ReflCorr.Clone("h1D_ReflCorr_BaselineSubtr")
        h1D_ReflCorr_BaselineSubtr.Sumw2()
        h1D_ReflCorr_BaselineSubtr.Add(hBaseline_Refl, -1.)
        
        self.fCorrectedCorrHisto_Reflected = h1D_ReflCorr.Clone(f"hCorrectedCorrReflected_PtCand{PtCandMin:.0f}to{PtCandMax:.0f}_PtAssoc{PtHadMin:.0f}to{PtHadMax:.0f}_InvMassBin{self.fBinInvMass}")
        self.fCorrectedCorrHisto_Reflected_BaselineSubtr = h1D_ReflCorr_BaselineSubtr.Clone(f"hCorrectedCorrReflected_BaselineSubtr_PtCand{PtCandMin:.0f}to{PtCandMax:.0f}_PtAssoc{PtHadMin:.0f}to{PtHadMax:.0f}_InvMassBin{self.fBinInvMass}")

        return True

    def GetSignalAndBackgroundForNorm(self, PtCandMin, PtCandMax):
        hMassFitSgnYield = self.fFileMass.Get(self.fMassHistoNameSgn)
        hMassFitBkgYield = self.fFileMass.Get(self.fMassHistoNameBkg)
        hMassFitSBsYield = self.fFileMass.Get(self.fMassHistoNameSBs)
        hMassFitSBLYield = self.fFileMass.Get("hBackgroundSidebandLeft")
        hMassFitSBRYield = self.fFileMass.Get("hBackgroundSidebandRight")

        PtCandBin = hMassFitSgnYield.FindBin(PtCandMin + 0.01)
        if PtCandBin != hMassFitSgnYield.FindBin(PtCandMax - 0.01):
            print("[ERROR] Pt bin in invariant mass histogram not univocally defined ")

        SgnYield = hMassFitSgnYield.GetBinContent(PtCandBin)
        BkgYield = hMassFitBkgYield.GetBinContent(PtCandBin)
        SBsYield = hMassFitSBsYield.GetBinContent(PtCandBin)
        SBLYield = hMassFitSBLYield.GetBinContent(PtCandBin)
        SBRYield = hMassFitSBRYield.GetBinContent(PtCandBin)

        print("================================= ")
        print("Getting invariant mass parameters ")
        print(f"Pt cand {PtCandMin} - {PtCandMax}")
        print(f"Signal yield    = {SgnYield}")
        print(f"Bkg yield    = {BkgYield}")
        print(f"Sideband yield    = {SBsYield}")
        print(f"Sideband left yield    = {SBLYield}")
        print(f"Sideband right yield    = {SBRYield}")
        print("================================= ")
        print(" ")

        self.SetSignalYieldforNorm(SgnYield)
        self.SetBkgYield(BkgYield)
        if self.fUseSidebLeft and self.fUseSidebRight:
            self.SetBkgScaleFactor(BkgYield / SBsYield)
            self.SetSBYield(SBsYield)
        elif self.fUseSidebLeft and not self.fUseSidebRight:
            self.SetBkgScaleFactor(BkgYield / SBLYield)
            self.SetSBYield(SBLYield)
        elif not self.fUseSidebLeft and self.fUseSidebRight:
            self.SetBkgScaleFactor(BkgYield / SBRYield)
            self.SetSBYield(SBRYield)

    def GetSignalAndBackgroundForNorm_v2(self, PtCandMin, PtCandMax, InvMassMin, InvMassMax):
        hMassSparse = self.fFileMass.Get(self.fMassHistoNameSgn)
        hMassFitSgnYield = hMassSparse.Projection(1, 0)
        hMassFitSgnYield.SetName("hMassVsPt")

        self.fMassVsPt = hMassFitSgnYield

        ptBinMin = hMassFitSgnYield.GetYaxis().FindBin(PtCandMin + 1e-6)
        ptBinMax = hMassFitSgnYield.GetYaxis().FindBin(PtCandMax - 1e-6)

        hMass = hMassFitSgnYield.ProjectionX(f"{hMassFitSgnYield.GetName()}_proj_mass", ptBinMin, ptBinMax)

        massBinMin = hMass.GetXaxis().FindBin(InvMassMin)
        massBinMax = hMass.GetXaxis().FindBin(InvMassMax)

        SgnYield = hMass.Integral(massBinMin, massBinMax)

        print("================================= ")
        print("Getting invariant mass parameters ")
        print(f"Pt cand {PtCandMin} - {PtCandMax}")
        print(f"Inv mass yield for norm   = {SgnYield}")
        print("================================= ")
        print(" ")

        self.SetSignalYieldforNorm(SgnYield)
