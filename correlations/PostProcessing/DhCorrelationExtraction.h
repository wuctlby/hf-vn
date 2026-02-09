// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DhCorrelationExtraction.h
/// \brief Class for D-h correlation extraction
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>

#ifndef PWGHF_HFC_MACROS_DHCORRELATIONEXTRACTION_H_
#define PWGHF_HFC_MACROS_DHCORRELATIONEXTRACTION_H_

#include <TAttMarker.h>
#include <TDirectoryFile.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TObject.h>
#include <TString.h>

#include <Rtypes.h>
#include <RtypesCore.h>
#include <vector>

class DhCorrelationExtraction : public TObject
{

 public:
  enum DmesonSpecie { kD0toKpi = 0,
                      kDplusKpipi,
                      kDsToKKPi,
                      kDStarD0pi };
  enum selectAnalysisType { kSE,
                            kME };
  enum selectParticleType { kPrimaryPart,
                            kAllPart };
  enum method {
    kMassBinning,
    kDeltaPhiBinning
  };
  enum axis { kPool = 0,
              kPtCand = 1,
              kPtHad = 2,
              kDeltaEta = 3,
              kDeltaPhi = 4,
              kMass = 5 };
  struct interval {
    Double_t low;
    Double_t up;
  };

  DhCorrelationExtraction(); // default constructor
  DhCorrelationExtraction(const DhCorrelationExtraction& source);
  virtual ~DhCorrelationExtraction();
  static DhCorrelationExtraction CreateDefault();
  static DhCorrelationExtraction* CreateCopy(const DhCorrelationExtraction& source);

  /// Input files, directories and thnsparse
  // D meson specie
  Bool_t SetDmesonSpecie(DmesonSpecie k);
  // mass
  void SetInputFilenameMass(TString filenameMass) { fFileNameMass = filenameMass; }
  void SetMassSparseName(TString massSparseName) { fMassSparseName = massSparseName; }
  // same event
  void SetInputFilenameSE(TString filenameSE) { fFileNameSE = filenameSE; }
  void SetDirNameSE(TString dirNameSE) { fDirNameSE = dirNameSE; }
  void SetCorrelSparseNameSE(TString correlSparseNameSE) { fCorrelSparseNameSE = correlSparseNameSE; }
  // mixed event
  void SetInputFilenameME(TString filenameME) { fFileNameME = filenameME; }
  void SetDirNameME(TString dirNameME) { fDirNameME = dirNameME; }
  void SetCorrelSparseNameME(TString correlSparseNameME) { fCorrelSparseNameME = correlSparseNameME; }
  // secondary particle
  void SetInputFilenameSecPart(TString filenameSecPart) { fFileSecPartName = filenameSecPart; }
  void SetDirNameSecPart(TString dirNameSecPart) { fDirSecPartName = dirNameSecPart; }
  void SetHistoSecPartName(TString histoPrimaryPartName, TString histoAllPartName)
  {
    fHistoPrimaryPartName = histoPrimaryPartName;
    fHistoAllPartName = histoAllPartName;
  }
  // method
  void SetMethod(method m) { fMethod = m; }
  // Debug level
  void SetDebugLevel(Int_t debug) { fDebug = debug; }

  /// Input settings
  // number of pools and pool by pool option
  void SetPoolSettings(Int_t nPools, Bool_t doPoolByPool)
  {
    fNpools = nPools;
    fDoPoolByPool = doPoolByPool;
  }
  // deltaEta bin for correlation histograms
  void SetBinDeltaEtaLeft(Double_t etaMin, Double_t etaMax)
  {
    fDeltaEtaLeftMin = etaMin;
    fDeltaEtaLeftMax = etaMax;
  }
  void SetBinDeltaEtaRight(Double_t etaMin, Double_t etaMax)
  {
    fDeltaEtaRightMin = etaMin;
    fDeltaEtaRightMax = etaMax;
  }
  // deltaPhi and deltaEta values for ME normalization
  void SetBinDeltaPhiEtaForMEnorm(Double_t valDeltaPhiMEnorm, Double_t valDeltaEtaMEnorm)
  { 
    fValDeltaPhiMEnorm = valDeltaPhiMEnorm; 
    fValDeltaEtaMEnorm = valDeltaEtaMEnorm;
  }

  /// Input conditions: PtCand, PtHad, PoolBins
  // void SetBinCandAndHad(Int_t binCand, Int_t binHad, Int_t binInvMass)
  // {
  //   fBinPtCand = binCand;
  //   fBinPtHad = binHad;
  //   fBinInvMass = binInvMass;
  // }
  void SetCandAndHadBins(std::vector<Double_t> PtCandBins, std::vector<Double_t> PtHadBins)
  {
    fPtCandBins = PtCandBins;
    fPtHadBins = PtHadBins;
  }
  void SetInvMassBins(std::vector<Double_t> InvMassBins)
  {
    fInvMassBins = InvMassBins;
  }
  void SetDeltaPhiBins(std::vector<Double_t> DeltaPhiBins)
  {
    fDeltaPhiBins = DeltaPhiBins;
  }
  void SetRebin2DcorrelHisto(Int_t rebinDeltaEta, Int_t rebinDeltaPhi)
  {
    fRebinAxisDeltaEta = rebinDeltaEta;
    fRebinAxisDeltaPhi = rebinDeltaPhi;
  }
  // optional settings
  void SetSubtractSoftPiInMEdistr(Bool_t subtractSoftPiME) { fdoSubtractSoftPiME = subtractSoftPiME; }
  void SetdoRebinSecondaryPart(Bool_t rebinSecPart) { fdoRebinSecPart = rebinSecPart; }
  void SetSecPartContamination(Bool_t secPartContamination) { fdoSecPartContamination = secPartContamination; }

  /// Analysis methods
  Bool_t ExtractCorrelations();
  Bool_t Init();
  Bool_t ReadInputSEandME();
  void ProjMassVsPt();
  TH2D* ProjCorrelHisto(Int_t SEorME, Int_t pool);
  TH1D* ProjCorrelHistoSecondaryPart(Int_t PrimaryPart, Double_t PtCandMin, Double_t PtCandMax, Double_t PtHadMin, Double_t PtHadMax);
  void NormalizeMEplot(TH2D*& histoME, TH2D*& histoMEsoftPi);
  Double_t CalculateTriggerNormalizationFactor(TH2D* hMassVsPt, Double_t PtCandMin, Double_t PtCandMax, Double_t InvMassMin, Double_t InvMassMax);
  TH1D* ReflectCorrHistogram(TH1D*& histo);
  TH1D* ReflectHistoRun2(TH1D* h, Double_t scale);
  Double_t CalculateBaseline(TH1D*& histo, Bool_t totalRange = kTRUE, Bool_t reflected = kFALSE);
  Double_t CalculateBaselineError(TH1D*& histo, Bool_t totalRange = kTRUE, Bool_t reflected = kFALSE);
  Bool_t ReadInputSecondaryPartContamination();

  /// Getters for the results
  TH2D* GetMassVsPtHist2D() { return fMassVsPt_2D; }
  // final results, debug level 0
  TH1D* GetCorrectedCorrHisto() { return fCorrectedCorrHisto; }
  TH1D* GetCorrectedMassPairs() { return fCorrectedMassPairs; }
  TH1D* GetCorrectedRatioVsDeltaEta() { return fCorrectedRatioVsDeltaEta; }
  TH1D* GetCorrectedCorrHisto_BaselineSubtr() { return fCorrectedCorrHisto_BaselineSubtr; }
  TH1D* GetCorrectedCorrHisto_Reflected() { return fCorrectedCorrHisto_Reflected; }
  TH1D* GetCorrectedCorrHisto_Reflected_BaselineSubtr() { return fCorrectedCorrHisto_Reflected_BaselineSubtr; }
  // intermediate results, debug level 1
  TH2D* GetCorrel_SE_2D() { return fCorrel_SE_2D; }
  TH2D* GetCorrel_ME_2D() { return fCorrel_ME_2D; }
  TH2D* GetCorrectedCorrel_2D() { return fCorrectedCorrel_2D; }
  TH1D* GetNonNormalizedCorrHisto() { return fNonNormalizedCorrHisto; }
  std::vector<TH1D*> GetVecCorrectedRatioVsDeltaEta() { return fVecCorrectedRatioVsDeltaEta; }
  std::vector<std::vector<TH1D*>> GetVecVecCorrectedMassPairsVsDeltaEta() { return fVecVecCorrectedMassPairsVsDeltaEta; }
  TH1D* GetBaselineHisto() { return fBaselineHisto; }
  TH1D* GetBaselineHisto_Reflected() { return fBaselineHisto_Reflected; }
  TH1D* GetCorrel_PrimaryPart() { return fCorrel_PrimaryPart; }
  TH1D* GetCorrel_AllPart() { return fCorrel_AllPart; }
  TH1D* GetFracSecondaryPart() { return fFracSecondaryPart; }
  TH1D* GetCorrectedCorrHisto_Before_SecPart() { return fCorrectedCorrHisto_Before_SecPart; }
  // original data histograms, debug level 2
  std::vector<TH2D*> GetVecOriginalCorrel_SE_2D() { return fVecOriginalCorrel_SE_2D; }
  std::vector<TH2D*> GetVecOriginalCorrel_ME_2D() { return fVecOriginalCorrel_ME_2D; }
  std::vector<TH2D*> GetVecCorrel_SE_2D() { return fVecCorrel_SE_2D; }
  std::vector<TH2D*> GetVecCorrel_ME_2D() { return fVecCorrel_ME_2D; }
  std::vector<TH2D*> GetVecCorrel_ME_norm_2D() { return fVecCorrel_ME_norm_2D; }
  std::vector<TH2D*> GetVecCorrectedCorrel_2D() { return fVecCorrectedCorrel_2D; }
  std::vector<TH2D*> GetVecMassPairsVsDeltaEta_2D() { return fVecMassPairsVsDeltaEta_2D; }
  TH1D* GetOriginalCorrel_PrimaryPart() { return fOriginalCorrel_PrimaryPart; }
  TH1D* GetOriginalCorrel_AllPart() { return fOriginalCorrel_AllPart; }

  /// Histogram style
  void SetTH1HistoStyle(TH1D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle, Style_t markerStyle = kFullCircle, Color_t markerColor = kRed + 1, Double_t markerSize = 1.4, Color_t lineColor = kRed + 1, Int_t lineWidth = 3, Float_t hTitleXaxisOffset = 1.0, Float_t hTitleYaxisOffset = 1.0, Float_t hTitleXaxisSize = 0.060, Float_t hTitleYaxisSize = 0.060, Float_t hLabelXaxisSize = 0.060, Float_t hLabelYaxisSize = 0.060, Bool_t centerXaxisTitle = false, Bool_t centerYaxisTitle = false);
  void SetTH2HistoStyle(TH2D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle, TString hZaxisTitle, Float_t hTitleXaxisOffset = 1.8, Float_t hTitleYaxisOffset = 1.8, Float_t hTitleZaxisOffset = 1.2, Float_t hTitleXaxisSize = 0.060, Float_t hTitleYaxisSize = 0.060, Float_t hTitleZaxisSize = 0.060, Float_t hLabelXaxisSize = 0.060, Float_t hLabelYaxisSize = 0.060, Float_t hLabelZaxisSize = 0.060, Bool_t centerXaxisTitle = true, Bool_t centerYaxisTitle = true);

 private:
  //Copyable properties-----------------------------------------------------------------------------------------------------//
  // Input configuration parameters，general
  DmesonSpecie fDmesonSpecies;           // D meson specie
  TString fDmesonLabel;                  // D meson label
  TString fFileNameMass;                 // File name contaning invariant mass output
  TString fFileNameSE;                   // File name contaning Same Event (SE) output
  TString fFileNameME;                   // File name contaning Mixed Event (ME) output
  TString fFileSecPartName;              // File name contaning secondary particle correction output
  TString fDirNameSE;                    // Directory in the file containing SE output
  TString fDirNameME;                    // Directory in the file containing ME output
  TString fDirSecPartName;               // Directory in the file containing secondary particle correction output
  TString fMassSparseName;                   // Inv. mass THnSparse name and directory
  TString fCorrelSparseNameSE;               // THnSparse name containing SE
  TString fCorrelSparseNameME;               // THnSparse name containing ME
  TString fHistoPrimaryPartName;         // Primary particle histogram (to be used for secondary particle contamination correction)
  TString fHistoAllPartName;             // All particle histogram (to be used for secondary particle contamination correction)
  
  Int_t fNpools;                         // Number of pools used for the ME correction
  Bool_t fDoPoolByPool;                  // Possibility to do the ME correction pool-by-pool (kTRUE) or merging all pools (kFALSE)
  Double_t fDeltaEtaLeftMin;                // DeltaEta min value
  Double_t fDeltaEtaLeftMax;                 // DeltaEta max value
  Double_t fDeltaEtaRightMin;                // DeltaEta min value
  Double_t fDeltaEtaRightMax;                 // DeltaEta max value
  Double_t fValDeltaPhiMEnorm;                // Delta phi value to ME normalisation if (0,0) bin is empty
  Double_t fValDeltaEtaMEnorm;                // Delta eta value to ME normalisation if (0,0) bin is empty
  
  Bool_t fdoSubtractSoftPiME;       // Soft pion subtraction (for D0 case)
  Bool_t fdoRebinSecPart;           // Rebin secondary particle contamination histogram
  Bool_t fdoSecPartContamination;   // Enable seconday particle contamination correction

  Int_t fDebug;                          // Debug level

  //Non-copyable properties-----------------------------------------------------------------------------------------------------//
  TFile* fFileSE;              // File containing the Same Event (SE) output
  TFile* fFileME;              // File containing the Mixed Event (ME) output
  TFile* fFileSecPart;         // File containing secondary particle contaminaion teplates
  TDirectoryFile* fDirSE;      // TDirectory for SE info
  TDirectoryFile* fDirME;      // TDirectory for ME info
  TDirectoryFile* fDirSecPart; // TDirectory for seondary particle correction

  std::vector<Double_t> fPtCandBins;   // Pt bins of the candidate
  std::vector<Double_t> fPtHadBins;    // Pt bins of the hadron
  std::vector<Double_t> fInvMassBins;  // Inv mass bins, whole range for deltaPhi binning method
  std::vector<Double_t> fDeltaPhiBins; // Delta phi bins, only for deltaPhi binning method
  Int_t fRebinAxisDeltaEta; // Rebin deltaEta axis value
  Int_t fRebinAxisDeltaPhi; // Rebin deltaPhi axis value

  Double_t fNpairs;         // Number of pairs per trigger in a given candidate and hadron pt bin
  Double_t fNpairsError;    // Error on number of pairs
  method fMethod; // Method to be used for the analysis

  // Results---------------------------------------------------------------------------------------------------------------//
  TH2D* fMassVsPt_2D;                                 // Invariant mass vs pt histo
  // final results, debug level 0
  TH1D* fCorrectedCorrHisto;                          // Corrected correlation histogram ME corrected
  TH1D* fCorrectedMassPairs;                      // Corrected pairs mass vs deltaEta histogram
  TH1D* fCorrectedRatioVsDeltaEta;                      // Ratio of corrected SE / SE vs deltaEta histogram
  TH1D* fCorrectedCorrHisto_BaselineSubtr;            // Corrected correlation histogram with baseline subtracted
  TH1D* fCorrectedCorrHisto_Reflected;                // Corrected correlation histogram reflected in azimuth
  TH1D* fCorrectedCorrHisto_Reflected_BaselineSubtr;  // Corrected correlation histogram reflected in azimuth with baseline subtracted
  // intermediate results, debug level 1
  TH2D* fCorrel_SE_2D;                                // Processed same Event correlation histogram 2D
  TH2D* fCorrel_ME_2D;                                // Processed Mixed Event correlation histogram 2D
  TH2D* fCorrel_ME_norm_2D;                           // Normalized Mixed Event correlation histogram 2D
  TH2D* fCorrectedCorrel_2D;                          // Corrected correlation histogram 2D
  TH1D* fNonNormalizedCorrHisto;                      //Corrected correlation histogram 1D, not normalized by the trigger, projected in deltaPhi
  std::vector<std::vector<TH1D*>> fVecVecCorrectedMassPairsVsDeltaEta;         // Vector of vector corrected mass distribution for each deltaEta bin and pool integrated
  std::vector<TH1D*> fVecCorrectedRatioVsDeltaEta;            // Vector of mass vs deltaEta histograms for each deltaEta bin and pool
  TH1D* fBaselineHisto;                               // Baseline histogram 1D, projected in deltaPhi
  TH1D* fBaselineHisto_Reflected;                     // Baseline histogram reflected in azimuth 1D, projected in deltaPhi
  TF1* fTFConstZero;                                  // Constant zero function
  TH1D* fCorrel_PrimaryPart;                          // Processed secondary particle contamination histogram for primary particles
  TH1D* fCorrel_AllPart;                              // Processed secondary particle contamination histogram for all particles
  TH1D* fFracSecondaryPart;                           // Fraction of secondary particles histogram
  TH1D* fCorrectedCorrHisto_Before_SecPart;           // Corrected correlation histogram before secondary particle contamination correction
  // original data histograms, debug level 2
  std::vector<TH2D*> fVecOriginalCorrel_SE_2D;        // Vector of original SE correlation histograms 2D
  std::vector<TH2D*> fVecOriginalCorrel_ME_2D;        // Vector of original ME correlation histograms 2D
  std::vector<TH2D*> fVecCorrel_SE_2D;                // Vector of processed SE correlation histograms 2D
  std::vector<TH2D*> fVecCorrel_ME_2D;                // Vector of processed ME correlation histograms 2D
  std::vector<TH2D*> fVecCorrel_ME_norm_2D;           // Vector of normalized ME correlation histograms 2D
  std::vector<TH2D*> fVecCorrectedCorrel_2D;          // Vector of corrected correlation histograms 2D
  std::vector<TH2D*> fVecMassPairsVsDeltaEta_2D;                  // Vector of mass vs deltaEta histograms
  TH1D* fOriginalCorrel_PrimaryPart;                  // Original secondary particle contamination histogram for primary particles
  TH1D* fOriginalCorrel_AllPart;                      // Original secondary particle contamination histogram for all particles
  ClassDef(DhCorrelationExtraction, 1);
};

#endif // PWGHF_HFC_MACROS_DHCORRELATIONEXTRACTION_H_
