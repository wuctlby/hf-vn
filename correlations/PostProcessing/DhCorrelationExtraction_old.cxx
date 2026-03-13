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

/// \file DhCorrelationExtraction.cxx
/// \brief class for D-h correlation extraction
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>

#include "DhCorrelationExtraction.h"
ClassImp(DhCorrelationExtraction);

#include <TCanvas.h>
#include <TDirectoryFile.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TString.h>

#include <RtypesCore.h>

#include <cstdio>
#include <iostream>

DhCorrelationExtraction::DhCorrelationExtraction() : // default constructor
  // Copyable properties-----------------------------------------------------------------------------------------------------//
  fDmesonSpecies(kDsToKKPi),
  fDmesonLabel("Ds"),
  fFileNameSE(""),
  fFileNameME(""),
  fFileSecPartName(""),
  fDirNameSE(""),
  fDirNameME(""),
  fDirSecPartName(""),
  fMassSparseName(""),
  fCorrelSparseNameSE(""),
  fCorrelSparseNameME(""),
  fHistoPrimaryPartName(""),
  fHistoAllPartName(""),
  fNpools(10),
  fDebug(0),
  fDoPoolByPool(kTRUE),
  fDeltaEtaLeftMin(-1.),
  fDeltaEtaLeftMax(1.),
  fDeltaEtaRightMin(-1.),
  fDeltaEtaRightMax(1.),
  fValDeltaPhiMEnorm(1.),
  fValDeltaEtaMEnorm(1.),
  fdoRebinSecPart(kFALSE),
  fdoSubtractSoftPiME(kFALSE),
  fdoSecPartContamination(0),
  // Non-copyable properties-----------------------------------------------------------------------------------//
  fFileSE(0x0),
  fFileME(0x0),
  fFileSecPart(0x0),
  fDirSE(0x0),
  fDirME(0x0),
  fDirSecPart(0x0),
  fPtCandBins(),
  fPtHadBins(),
  fInvMassBins(),
  fDeltaPhiBins(),
  fFactorsNormME(),
  fMethod(kMassBinning),
  fRebinAxisDeltaEta(1),
  fRebinAxisDeltaPhi(1),
  fNpairs(0.0),
  // Results---------------------------------------------------------------------------------------------------------------//
  fMassVsPt_2D(0x0), // copyable properties, for trigger normalization
  // final results, debug level 0
  fCorrectedCorrHisto(0x0),
  fCorrectedMassPairs(0x0),
  fCorrectedRatioVsDeltaEta(0x0),
  fCorrectedCorrHisto_BaselineSubtr(0x0),
  fCorrectedCorrHisto_Reflected(0x0),
  fCorrectedCorrHisto_Reflected_BaselineSubtr(0x0),
  // intermediate results, debug level 1
  fCorrel_SE_2D(0x0),
  fCorrel_ME_2D(0x0),
  fCorrel_ME_norm_2D(0x0),
  fCorrectedCorrel_2D(0x0),
  fNonNormalizedCorrHisto(0x0),
  fVecVecCorrectedMassPairsVsDeltaEta(),
  fVecCorrectedRatioVsDeltaEta(),
  fBaselineHisto(0x0),
  fBaselineHisto_Reflected(0x0),
  fTFConstZero(0x0),
  fCorrel_PrimaryPart(0x0),         // secondary particle contamination histograms
  fCorrel_AllPart(0x0),
  fFracSecondaryPart(0x0),
  fCorrectedCorrHisto_Before_SecPart(0x0),
  // original data histograms, debug level 2
  fVecOriginalCorrel_SE_2D(),
  fVecOriginalCorrel_ME_2D(),
  fVecCorrel_SE_2D(),
  fVecCorrel_ME_2D(),
  fVecCorrel_ME_norm_2D(),
  fVecCorrectedCorrel_2D(),
  fVecMassPairsVsDeltaEta_2D(),
  fOriginalCorrel_PrimaryPart(0x0), // secondary particle contamination histograms
  fOriginalCorrel_AllPart(0x0)
{
}

DhCorrelationExtraction::DhCorrelationExtraction(const DhCorrelationExtraction& source) : // copy constructor
  TObject(source),
  // Copyable properties-----------------------------------------------------------------------------------------------------//
  fDmesonSpecies(source.fDmesonSpecies),
  fDmesonLabel(source.fDmesonLabel),
  fFileNameSE(source.fFileNameSE),
  fFileNameME(source.fFileNameME),
  fFileSecPartName(source.fFileSecPartName),
  fDirNameSE(source.fDirNameSE),
  fDirNameME(source.fDirNameME),
  fDirSecPartName(source.fDirSecPartName),
  fMassSparseName(source.fMassSparseName),
  fCorrelSparseNameSE(source.fCorrelSparseNameSE),
  fCorrelSparseNameME(source.fCorrelSparseNameME),
  fHistoPrimaryPartName(source.fHistoPrimaryPartName),
  fHistoAllPartName(source.fHistoAllPartName),
  fNpools(source.fNpools),
  fDebug(source.fDebug),
  fDoPoolByPool(source.fDoPoolByPool),
  fDeltaEtaLeftMin(source.fDeltaEtaLeftMin),
  fDeltaEtaLeftMax(source.fDeltaEtaLeftMax),
  fDeltaEtaRightMin(source.fDeltaEtaRightMin),
  fDeltaEtaRightMax(source.fDeltaEtaRightMax),
  fValDeltaPhiMEnorm(source.fValDeltaPhiMEnorm),
  fValDeltaEtaMEnorm(source.fValDeltaEtaMEnorm),
  fdoRebinSecPart(source.fdoRebinSecPart),
  fdoSubtractSoftPiME(source.fdoSubtractSoftPiME),
  fdoSecPartContamination(source.fdoSecPartContamination),
  // Non-copyable properties-----------------------------------------------------------------------------------//
  fFileSE(0x0),
  fFileME(0x0),
  fFileSecPart(0x0),
  fDirSE(0x0),
  fDirME(0x0),
  fDirSecPart(0x0),
  fPtCandBins(),
  fPtHadBins(),
  fInvMassBins(),
  fDeltaPhiBins(),
  fFactorsNormME(),
  fMethod(kMassBinning),
  fRebinAxisDeltaEta(1),
  fRebinAxisDeltaPhi(1),
  fNpairs(0.0),
  // Results---------------------------------------------------------------------------------------------------------------//
  fMassVsPt_2D(source.fMassVsPt_2D), // copyable properties, for trigger normalization
  // final results, debug level 0
  fCorrectedCorrHisto(0x0),
  fCorrectedMassPairs(0x0),
  fCorrectedRatioVsDeltaEta(0x0),
  fCorrectedCorrHisto_BaselineSubtr(0x0),
  fCorrectedCorrHisto_Reflected(0x0),
  fCorrectedCorrHisto_Reflected_BaselineSubtr(0x0),
  // intermediate results, debug level 1
  fCorrel_SE_2D(0x0),
  fCorrel_ME_2D(0x0),
  fCorrel_ME_norm_2D(0x0),
  fCorrectedCorrel_2D(0x0),
  fNonNormalizedCorrHisto(0x0),
  fVecVecCorrectedMassPairsVsDeltaEta(),  
  fVecCorrectedRatioVsDeltaEta(),
  fBaselineHisto(0x0),
  fBaselineHisto_Reflected(0x0),
  fTFConstZero(0x0),
  fCorrel_PrimaryPart(0x0),         // secondary particle contamination histograms
  fCorrel_AllPart(0x0),
  fFracSecondaryPart(0x0),
  fCorrectedCorrHisto_Before_SecPart(0x0),
  // original data histograms, debug level 2
  fVecOriginalCorrel_SE_2D(),
  fVecOriginalCorrel_ME_2D(),
  fVecCorrel_SE_2D(),
  fVecCorrel_ME_2D(),
  fVecCorrel_ME_norm_2D(),
  fVecCorrectedCorrel_2D(),
  fVecMassPairsVsDeltaEta_2D(),
  fOriginalCorrel_PrimaryPart(0x0), // secondary particle contamination histograms
  fOriginalCorrel_AllPart(0x0)
{
}

DhCorrelationExtraction::~DhCorrelationExtraction()
// destructor
{
}

DhCorrelationExtraction DhCorrelationExtraction::CreateDefault() {
    static DhCorrelationExtraction instance;
    instance.SetDirNameSE("hf-correlator-flow-charm-hadrons-reduced");
    instance.SetDirNameME("hf-correlator-flow-charm-hadrons-reduced");
    instance.SetCorrelSparseNameSE("hSparseCorrelationsSECharmHad");
    instance.SetCorrelSparseNameME("hSparseCorrelationsMECharmHad");
    instance.SetMassSparseName("hf-correlator-flow-charm-hadrons-reduced/hSparseTrigCandsCharm");
    return instance;
}

DhCorrelationExtraction* DhCorrelationExtraction::CreateCopy(const DhCorrelationExtraction& source) {
    return new DhCorrelationExtraction(source);
}

Bool_t DhCorrelationExtraction::SetDmesonSpecie(DmesonSpecie k)
{

  if (k < 0 || k > 3) {
    printf("[ERROR] D meson specie not correctly set!\n");
    return kFALSE;
  } else if (k == 0) {
    fDmesonLabel = "Dzero";
  } else if (k == 1) {
    fDmesonLabel = "Dplus";
  } else if (k == 2) {
    fDmesonLabel = "Ds";
  } else {
    fDmesonLabel = "Dstar";
  }
  fDmesonSpecies = k;
  return kTRUE;
}

Bool_t DhCorrelationExtraction::Init()
{
  // Open input SE and ME files
  // Bool_t success = ReadInputSEandME();
  if (!ReadInputSEandME()) {
    return kFALSE;
  }
  return kTRUE;
}

Bool_t DhCorrelationExtraction::ExtractCorrelations()
{
  TH1::AddDirectory(kFALSE);

  if (Init() == kFALSE) {
    return kFALSE;
  }

  if (fdoSubtractSoftPiME) {
    printf("[INFO] Fake softPi subtraction in ME via extraction code is enabled!\n");
  }

  if (!fDoPoolByPool)
    fNpools = 1; // single histogram with integrated pools

  fFactorsNormME.resize(fNpools, 1.0);
  // Histograms definition
  std::vector<TH2D*> hSE_Sign(fNpools, nullptr);
  std::vector<TH2D*> hME_Sign(fNpools, nullptr);
  std::vector<TH2D*> hME_Normalized(fNpools, nullptr);
  std::vector<TH2D*> hME_Sign_SoftPi(fNpools, nullptr);
  std::vector<TH2D*> hCorrected_Sign(fNpools, nullptr);
  
  // todo: sideband region
  std::vector<TH2D*> hSE_Sideb(fNpools, nullptr);
  std::vector<TH2D*> hME_Sideb(fNpools, nullptr);
  std::vector<TH2D*> hME_Sideb_SoftPi(fNpools, nullptr);
  std::vector<TH2D*> hCorrected_Sideb(fNpools, nullptr);

  std::vector<TH1D*> hCorrectedMassPairs(fNpools, nullptr);
  TH1D* h1D_correctedMassPairs = nullptr;
  TH1D* h1D_correctedRatioVsDeltaEta = nullptr;

  TH2D* h2D_SE = nullptr;
  TH2D* h2D_ME = nullptr;
  TH2D* h2D_ME_norm = nullptr;
  TH2D* h2D_Sign = nullptr;

  TH1D* h1D_Sign = nullptr;
  TH1D* h1D_Norm = nullptr;
  TH1D* h1D_PrimaryPartCorr = nullptr;
  TH1D* h1D_AllPartCorr = nullptr;
  TH1D* h1D_SecPartFrac = nullptr;
  TH1D* h1D_Norm_SecPart = nullptr;
  TH1D* h1D_BaselineSubtr = nullptr;
  TH1D* h1D_ReflCorr = nullptr;
  TH1D* h1D_ReflCorr_BaselineSubtr = nullptr;


  for (int iPool = 0; iPool < fNpools; iPool++) {
    // Retrieve 2D plots for SE and ME, signal and bkg regions, for each pTbin and pool
    hSE_Sign[iPool] = ProjCorrelHisto(kSE, iPool);
    hME_Sign[iPool] = ProjCorrelHisto(kME, iPool);

    hSE_Sign[iPool]->Sumw2();
    hME_Sign[iPool]->Sumw2();

    // rebin axes deltaEta and deltaPhi
    if (fRebinAxisDeltaEta > 1 || fRebinAxisDeltaPhi > 1) {
      hSE_Sign[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi); // Xaxis: deltaEta, Yaxis: deltaPhi
      hME_Sign[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
    }

    hME_Normalized[iPool] = reinterpret_cast<TH2D*>(hME_Sign[iPool]->Clone(Form("hNormalizedCorrel_ME_Pool%d", iPool)));
    // Normalize ME plots for the entries in (deltaEta, deltaPhi) = (0, 0)
    NormalizeMEplot(hME_Normalized[iPool], hME_Sign_SoftPi[iPool], iPool);

    // Apply Event Mixing Correction
    hCorrected_Sign[iPool] = reinterpret_cast<TH2D*>(hSE_Sign[iPool]->Clone(Form("hCorrected_Sign_Pool%d", iPool)));
    hCorrected_Sign[iPool]->Sumw2();
    hCorrected_Sign[iPool]->Divide(hME_Normalized[iPool]);

    // Apply the ME correction on the Mass by the ratio of SE/ME integrated over deltaPhi bins for each deltaEta bin
    if (fMethod == kDeltaPhiBinning) {
      TH1D* hTempCorrectedRatioVsDeltaEta = nullptr;
      std::vector<TH1D*> VecTempCorrectedMassVsDeltaEta;

      // Apply the correction in the whole deltaPhi range
      Double_t original_content = hSE_Sign[iPool]->Integral(hSE_Sign[iPool]->GetXaxis()->FindBin(fDeltaEtaLeftMin), hSE_Sign[iPool]->GetXaxis()->FindBin(fDeltaEtaRightMax), 
                                                            hSE_Sign[iPool]->GetYaxis()->FindBin(fDeltaPhiBins.front()), hSE_Sign[iPool]->GetYaxis()->FindBin(fDeltaPhiBins.back()));
      if (original_content == 0) {
        printf("[WARNING] Original content in SE signal histogram for pool %d is 0 in the selected deltaPhi range! Skipping ME correction for mass histogram in this pool.\n", iPool);
        continue;
      }
      Double_t corrected_content = hCorrected_Sign[iPool]->Integral(hCorrected_Sign[iPool]->GetXaxis()->FindBin(fDeltaEtaLeftMin), hCorrected_Sign[iPool]->GetXaxis()->FindBin(fDeltaEtaRightMax), 
                                                                    hCorrected_Sign[iPool]->GetYaxis()->FindBin(fDeltaPhiBins.front()), hCorrected_Sign[iPool]->GetYaxis()->FindBin(fDeltaPhiBins.back()));
      Double_t ratio_correction = corrected_content / original_content;
      TH2D* tempClone = reinterpret_cast<TH2D*>(fVecMassPairsVsDeltaEta_2D[iPool]->Clone(Form("hCorrectedMassPairs_Pool%d", iPool)));
      hCorrectedMassPairs[iPool] = tempClone->ProjectionY(Form("hCorrectedMassPairs_ProjY_Pool%d", iPool), 1, fVecMassPairsVsDeltaEta_2D[iPool]->GetXaxis()->GetNbins());
      hCorrectedMassPairs[iPool]->Scale(ratio_correction);
      if (fDebug > 0) {
        VecTempCorrectedMassVsDeltaEta.push_back(reinterpret_cast<TH1D*>(hCorrectedMassPairs[iPool]->Clone(Form("hCorrectedMassPairs_Pool%d", iPool))));
        VecTempCorrectedMassVsDeltaEta.back()->SetTitle(Form("Corrected associated pairs mass, Pool %d, DeltaEta %f - %f", iPool, hCorrected_Sign[iPool]->GetXaxis()->GetBinLowEdge(1), hCorrected_Sign[iPool]->GetXaxis()->GetBinUpEdge(hCorrected_Sign[iPool]->GetXaxis()->GetNbins())));
        hTempCorrectedRatioVsDeltaEta = new TH1D(Form("hCorrectedRatioVsDeltaEta_Pool%d", iPool), Form("Correction ratio vs DeltaEta, DeltaPhi %f - %f, Pool %d", fDeltaPhiBins.front(), fDeltaPhiBins.back(), iPool),
                                              hCorrected_Sign[iPool]->GetXaxis()->GetNbins(), hCorrected_Sign[iPool]->GetXaxis()->GetXmin(), hCorrected_Sign[iPool]->GetXaxis()->GetXmax());
        hTempCorrectedRatioVsDeltaEta->SetTitle(Form("Correction ratio vs DeltaEta, DeltaPhi %f - %f, Pool %d", fDeltaPhiBins.front(), fDeltaPhiBins.back(), iPool));
        // for (int iEtaBin = 1; iEtaBin <= hCorrected_Sign[iPool]->GetXaxis()->GetNbins(); iEtaBin++) {
        //   Double_t original_content_bin = hSE_Sign[iPool]->Integral(iEtaBin, iEtaBin, hSE_Sign[iPool]->GetYaxis()->FindBin(fDeltaPhiBins.front()), hSE_Sign[iPool]->GetYaxis()->FindBin(fDeltaPhiBins.back()));
        //   if (original_content_bin == 0) {
        //     hTempCorrectedRatioVsDeltaEta->SetBinContent(iEtaBin, 0);
        //     continue;
        //   }
        //   Double_t corrected_content_bin = hCorrected_Sign[iPool]->Integral(iEtaBin, iEtaBin, hCorrected_Sign[iPool]->GetYaxis()->FindBin(fDeltaPhiBins.front()), hCorrected_Sign[iPool]->GetYaxis()->FindBin(fDeltaPhiBins.back()));
        //   Double_t ratio_correction_bin = corrected_content_bin / original_content_bin;
          hTempCorrectedRatioVsDeltaEta->SetBinContent(1, ratio_correction);
        // }
        fVecCorrectedRatioVsDeltaEta.push_back(reinterpret_cast<TH1D*>(hTempCorrectedRatioVsDeltaEta->Clone(Form("hCorrectedRatioVsDeltaEta_Pool%d", iPool))));
        fVecVecCorrectedMassPairsVsDeltaEta.push_back(VecTempCorrectedMassVsDeltaEta);
        delete hTempCorrectedRatioVsDeltaEta;
        VecTempCorrectedMassVsDeltaEta.clear();
      }
      // Apply the correction in each deltaEta bin and sum up the corrected mass histograms
      // for (int iEtaBin = 1; iEtaBin <= hCorrected_Sign[iPool]->GetXaxis()->GetNbins(); iEtaBin++) {
      //   TH1D* hTempCorrectedMass = nullptr;
      //   Double_t original_content = hSE_Sign[iPool]->Integral(iEtaBin, iEtaBin, 
      //                                                         hCorrected_Sign[iPool]->GetYaxis()->FindBin(fDeltaPhiBins.front()),
      //                                                         hCorrected_Sign[iPool]->GetYaxis()->FindBin(fDeltaPhiBins.back())
      //                                                       );
      //   if (original_content == 0) {
      //     continue;
      //   }
      //   Double_t corrected_content = hCorrected_Sign[iPool]->Integral(iEtaBin, iEtaBin, 
      //                                                                 hCorrected_Sign[iPool]->GetYaxis()->FindBin(fDeltaPhiBins.front()), 
      //                                                                 hCorrected_Sign[iPool]->GetYaxis()->FindBin(fDeltaPhiBins.back())
      //                                                               );
      //   Double_t ratio_correction = corrected_content / original_content;
      //   TH2D* tempClone = reinterpret_cast<TH2D*>(fVecMassPairsVsDeltaEta_2D[iPool]->Clone(Form("hCorrectedMassPairs_Pool%d_EtaBin%d", iPool, iEtaBin)));
      //   hTempCorrectedMass = tempClone->ProjectionY(Form("hCorrectedMassPairs_ProjY_Pool%d_EtaBin%d", iPool, iEtaBin), iEtaBin, iEtaBin);
      //   hTempCorrectedMass->Scale(ratio_correction);
      //   if (iEtaBin == 1) {
      //     hCorrectedMassPairs[iPool] = reinterpret_cast<TH1D*>(hTempCorrectedMass->Clone(Form("hCorrectedMassPairs_Pool%d", iPool)));
      //   } else {
      //     hCorrectedMassPairs[iPool]->Add(hTempCorrectedMass);
      //   }
      //   if (fDebug > 0) {
      //     VecTempCorrectedMassVsDeltaEta.push_back(reinterpret_cast<TH1D*>(hTempCorrectedMass->Clone(Form("hCorrectedMassPairs_Pool%d_EtaBin%d", iPool, iEtaBin))));
      //     VecTempCorrectedMassVsDeltaEta.back()->SetTitle(Form("Corrected associated pairs mass, Pool %d, DeltaEta %f - %f", iPool, hCorrected_Sign[iPool]->GetXaxis()->GetBinLowEdge(iEtaBin), hCorrected_Sign[iPool]->GetXaxis()->GetBinUpEdge(iEtaBin)));
      //     if (iEtaBin == 1) {
      //       hTempCorrectedRatioVsDeltaEta = new TH1D(Form("hCorrectedRatioVsDeltaEta_Pool%d", iPool), Form("Correction ratio vs DeltaEta, DeltaPhi %f - %f, Pool %d", fDeltaPhiBins.front(), fDeltaPhiBins.back(), iPool),
      //                                             hCorrected_Sign[iPool]->GetXaxis()->GetNbins(), hCorrected_Sign[iPool]->GetXaxis()->GetXmin(), hCorrected_Sign[iPool]->GetXaxis()->GetXmax());
      //       hTempCorrectedRatioVsDeltaEta->SetTitle(Form("Correction ratio vs DeltaEta, DeltaPhi %f - %f, Pool %d", fDeltaPhiBins.front(), fDeltaPhiBins.back(), iPool));
      //     }
      //     hTempCorrectedRatioVsDeltaEta->SetBinContent(iEtaBin, ratio_correction);
      //   }
      //   delete tempClone;
      //   delete hTempCorrectedMass;
      // }
      // fVecCorrectedRatioVsDeltaEta.push_back(reinterpret_cast<TH1D*>(hTempCorrectedRatioVsDeltaEta->Clone(Form("hCorrectedRatioVsDeltaEta_Pool%d", iPool))));
      // fVecVecCorrectedMassPairsVsDeltaEta.push_back(VecTempCorrectedMassVsDeltaEta);

      // delete hTempCorrectedRatioVsDeltaEta;
      // VecTempCorrectedMassVsDeltaEta.clear();
    }

    // Set proper number of entries after ME correction
    Double_t N_SEsign = 0, N_sign = 0, N_massPairs = 0;
    for (int i = 1; i <= hCorrected_Sign[iPool]->GetXaxis()->GetNbins(); i++) {
      for (int j = 1; j <= hCorrected_Sign[iPool]->GetYaxis()->GetNbins(); j++) {
        N_SEsign += hSE_Sign[iPool]->GetBinContent(i, j);
        N_sign += hCorrected_Sign[iPool]->GetBinContent(i, j);
      }
      for (int k = 1; k <= hCorrectedMassPairs[iPool]->GetXaxis()->GetNbins(); k++) {
        N_massPairs += hCorrectedMassPairs[iPool]->GetBinContent(i, k);
      }
    }
    hSE_Sign[iPool]->SetEntries(N_SEsign);
    hCorrected_Sign[iPool]->SetEntries(N_sign);
    hCorrectedMassPairs[iPool]->SetEntries(N_massPairs);
    std::cout << "[Debug] ramdomly check the number of entries after ME correction: " << hCorrectedMassPairs[iPool]->GetBinContent(5) << std::endl;

    // debug: normalized ME and corrected correlation histos pool by pool
    if (fDebug > 0) {
      fVecCorrel_ME_norm_2D.push_back(reinterpret_cast<TH2D*>(hME_Normalized[iPool]->Clone(Form("hNorm_Correl_ME_Pool%s_2D", fDoPoolByPool ? Form("%d", iPool) : "All"))));
      fVecCorrectedCorrel_2D.push_back(reinterpret_cast<TH2D*>(hCorrected_Sign[iPool]->Clone(Form("hCorrected_Correl_2D_Pool%s", fDoPoolByPool ? Form("%d", iPool) : "All"))));
    }
 
    // Pools integration
    if (iPool == 0) {
      h2D_SE = reinterpret_cast<TH2D*>(hSE_Sign[0]->Clone("h2D_SE"));
      h2D_ME = reinterpret_cast<TH2D*>(hME_Sign[0]->Clone("h2D_ME"));
      h2D_ME_norm = reinterpret_cast<TH2D*>(hME_Normalized[0]->Clone("h2D_ME_norm"));
      h2D_Sign = reinterpret_cast<TH2D*>(hCorrected_Sign[0]->Clone("h2D_Sign"));
      h1D_correctedMassPairs = reinterpret_cast<TH1D*>(hCorrectedMassPairs[0]->Clone("h1D_correctedMassPairs"));
      h1D_correctedMassPairs->SetTitle("Corrected associated pairs mass distribution");
      h1D_correctedRatioVsDeltaEta = reinterpret_cast<TH1D*>(fVecCorrectedRatioVsDeltaEta[0]->Clone("h1D_correctedRatioVsDeltaEta"));
      h1D_correctedRatioVsDeltaEta->SetTitle(Form("Correction ratio vs DeltaEta, DeltaPhi %f - %f", fDeltaPhiBins.front(), fDeltaPhiBins.back()));
      h2D_SE->Sumw2();
      h2D_ME->Sumw2();
      h2D_ME_norm->Sumw2();
      h2D_Sign->Sumw2();
    } else {
      h2D_Sign->Add(hCorrected_Sign[iPool]);
      h2D_SE->Add(hSE_Sign[iPool]);
      h2D_ME->Add(hME_Sign[iPool]);
      h2D_ME_norm->Add(hME_Normalized[iPool]);
      h1D_correctedMassPairs->Add(hCorrectedMassPairs[iPool]);
      h1D_correctedRatioVsDeltaEta->Add(fVecCorrectedRatioVsDeltaEta[iPool]);
    }
  } // end pool loop

  // clean up pool histos
  for (int iPool = 0; iPool < fNpools; iPool++) {
    delete hSE_Sign[iPool];
    hSE_Sign[iPool] = nullptr;
    delete hME_Sign[iPool];
    hME_Sign[iPool] = nullptr;
    delete hME_Normalized[iPool];
    hME_Normalized[iPool] = nullptr;
    delete hCorrected_Sign[iPool];
    hCorrected_Sign[iPool] = nullptr;
    delete hME_Sign_SoftPi[iPool];
    hME_Sign_SoftPi[iPool] = nullptr;
    delete hCorrectedMassPairs[iPool];
    hCorrectedMassPairs[iPool] = nullptr;
  }

  std::cout << "[Debug] ramdomly check the bin content after pool integration: " << h1D_correctedMassPairs->GetBinContent(5) << std::endl;
  fCorrectedMassPairs = reinterpret_cast<TH1D*>(h1D_correctedMassPairs->Clone("hCorrectedMassPairs_AllPools"));
  fCorrectedRatioVsDeltaEta = reinterpret_cast<TH1D*>(h1D_correctedRatioVsDeltaEta->Clone("hCorrectedRatioVsDeltaEta_AllPools"));

  // debug: integrated SE, ME, normalized ME and corrected correlation histos
  if (fDebug > 0) {
    fCorrel_SE_2D = reinterpret_cast<TH2D*>(h2D_SE->Clone(Form("hCorrel_SE_allPools")));
    fCorrel_ME_2D = reinterpret_cast<TH2D*>(h2D_ME->Clone(Form("hCorrel_ME_allPools")));
    fCorrel_ME_norm_2D = reinterpret_cast<TH2D*>(h2D_ME_norm->Clone(Form("hCorrel_ME_norm_allPools")));
    fCorrectedCorrel_2D = reinterpret_cast<TH2D*>(h2D_Sign->Clone(Form("hCorrectedCorrel_allPools")));
  }

  // clean up integrated SE, ME and normalized ME histos
  delete h2D_SE;
  h2D_SE = nullptr;
  delete h2D_ME;
  h2D_ME = nullptr;
  delete h2D_ME_norm;
  h2D_ME_norm = nullptr;
  delete h1D_correctedMassPairs;
  h1D_correctedMassPairs = nullptr;

  //==========================================================================================================================
  // 1D projection
  h1D_Sign = reinterpret_cast<TH1D*>(h2D_Sign->ProjectionY("h1D_Sign")); // projection on deltaPhi axis
  h1D_Sign->Sumw2();

  // Apply normalization to number of triggers - NOT DONE
  h1D_Norm = reinterpret_cast<TH1D*>(h1D_Sign->Clone("h1D_triggerNorm_Sign"));
  h1D_Norm->Sumw2();
  if (!fMassVsPt_2D) {
    ProjMassVsPt();
  }
  Double_t N_triggers = CalculateTriggerNormalizationFactor(fMassVsPt_2D, fPtCandBins[0], fPtCandBins[1], fInvMassBins[0], fInvMassBins[1]);
  if (fDebug > 0) {
    TH1D* hDebug_CorrBeforeNorm = reinterpret_cast<TH1D*>(h1D_Sign->Clone("hCorrectedCorr_BeforeNorm"));
    fNonNormalizedCorrHisto = reinterpret_cast<TH1D*>(h1D_Sign->Clone("hNonNormalizedCorrHisto"));
    SetTH1HistoStyle(fNonNormalizedCorrHisto, Form("%s-h Correlation - Before trigger [%.0f], %.0f < p_{T} < %.0f GeV/c", fDmesonLabel.Data(), N_triggers, fPtCandBins[0], fPtCandBins[1]), "#Delta#phi [rad]", "dN/d#Delta#phi");
    hDebug_CorrBeforeNorm = nullptr;
    delete hDebug_CorrBeforeNorm;
  }
  h1D_Norm->Scale(1. / N_triggers);
  h1D_Norm->SetTitle(Form("%s-h Correlation - Normalized to number of triggers, %.0f", fDmesonLabel.Data(), N_triggers));
  SetTH1HistoStyle(h1D_Norm, Form("%.0f < p_{T} < %.0f GeV/c", fPtCandBins[0], fPtCandBins[1]), "#Delta#phi [rad]", "1/N_{trig} dN/d#Delta#phi");

  // clean up
  delete h2D_Sign;
  h2D_Sign = nullptr;
  delete h1D_Sign;
  h1D_Sign = nullptr;

  // Secondary particle contamination
  if (fdoSecPartContamination) {
    h1D_PrimaryPartCorr = ProjCorrelHistoSecondaryPart(kPrimaryPart, fPtCandBins[0], fPtCandBins[1], fPtHadBins[0], fPtHadBins[1]);
    h1D_AllPartCorr = ProjCorrelHistoSecondaryPart(kAllPart, fPtCandBins[0], fPtCandBins[1], fPtHadBins[0], fPtHadBins[1]);
    h1D_PrimaryPartCorr->Sumw2();
    h1D_AllPartCorr->Sumw2();
    if (fdoRebinSecPart && (fRebinAxisDeltaPhi > 1 || fRebinAxisDeltaEta > 1)) {
      h1D_PrimaryPartCorr->RebinX(fRebinAxisDeltaPhi); // Xaxis: deltaPhi
      h1D_AllPartCorr->RebinX(fRebinAxisDeltaPhi);
    }

    h1D_SecPartFrac = reinterpret_cast<TH1D*>(h1D_PrimaryPartCorr->Clone(Form("hCorrRatio_PtD%.0fto%.0f_PtHad%.0fto%.0f", fPtCandBins[0], fPtCandBins[1], fPtHadBins[0], fPtHadBins[1])));
    h1D_SecPartFrac->Sumw2();
    h1D_SecPartFrac->Divide(h1D_PrimaryPartCorr, h1D_AllPartCorr, 1., 1., "B");
    SetTH1HistoStyle(h1D_SecPartFrac, Form("%.0f < p_{T} < %.0f GeV/c", fPtCandBins[0], fPtCandBins[1]), "#Delta#phi [rad]", "#frac{primary part.}{part. selected}");

    if (fDebug >= 1) {
      fFracSecondaryPart = reinterpret_cast<TH1D*>(h1D_SecPartFrac->Clone(Form("hFracSecondaryPart_PtD%.0fto%.0f_PtHad%.0fto%.0f", fPtCandBins[0], fPtCandBins[1], fPtHadBins[0], fPtHadBins[1])));
      fCorrectedCorrHisto_Before_SecPart = reinterpret_cast<TH1D*>(h1D_Norm->Clone(Form("hCorrectedCorr_Before_SecPart_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f_InvMass%.0fto%.0f", fPtCandBins[0], fPtCandBins[1], fPtHadBins[0], fPtHadBins[1], fInvMassBins[0], fInvMassBins[1])));
    }
  
    h1D_Norm_SecPart = reinterpret_cast<TH1D*>(h1D_Norm->Clone("h1D_Norm_SecPart"));
    h1D_Norm_SecPart->Sumw2();
    Int_t nBinsPhi = h1D_Norm_SecPart->GetNbinsX();
    if (nBinsPhi != h1D_SecPartFrac->GetNbinsX()) {
      std::cerr << "[ERROR]: nBinsPhi different between h1D_Norm and h1D_SecPartFrac" << std::endl;
      return kFALSE;
    }
    h1D_Norm_SecPart->Multiply(h1D_SecPartFrac);

    // clean up
    delete h1D_PrimaryPartCorr;
    h1D_PrimaryPartCorr = nullptr;
    delete h1D_AllPartCorr;
    h1D_AllPartCorr = nullptr;
    delete h1D_SecPartFrac;
    h1D_SecPartFrac = nullptr;
  }

  // set 1D plots (Signal region, normalized)
  h1D_Norm->SetLineColor(kBlue + 1);
  h1D_Norm->SetMarkerColor(kBlue + 1);
  h1D_Norm->SetMarkerStyle(kFullCircle);
  h1D_Norm->SetMinimum(0);
  if (fdoSecPartContamination) {
    h1D_Norm_SecPart->SetLineColor(kRed + 1);
    h1D_Norm_SecPart->SetMarkerColor(kRed + 1);
    h1D_Norm_SecPart->SetMarkerStyle(kFullCircle);
  }

  if (fdoSecPartContamination) {
    fCorrectedCorrHisto = reinterpret_cast<TH1D*>(h1D_Norm_SecPart->Clone("hCorrectedCorr"));
  } else {
    fCorrectedCorrHisto = reinterpret_cast<TH1D*>(h1D_Norm->Clone("hCorrectedCorr"));
  }

  // Baseline subtraction
  Double_t BaselineData, BaselineDataErr;
  TH1D* hBaseline = reinterpret_cast<TH1D*>(h1D_Norm->Clone("hBaseline"));
  hBaseline->Sumw2();

  if (fdoSecPartContamination) {
    BaselineData = CalculateBaseline(h1D_Norm_SecPart, kTRUE, kFALSE);
    BaselineDataErr = CalculateBaselineError(h1D_Norm_SecPart, kTRUE, kFALSE);
    for (int iBin = 0; iBin < hBaseline->GetNbinsX(); iBin++) {
      hBaseline->SetBinContent(iBin + 1, BaselineData);
      hBaseline->SetBinError(iBin + 1, BaselineDataErr);
    }
    h1D_BaselineSubtr = reinterpret_cast<TH1D*>(h1D_Norm_SecPart->Clone("h1D_BaselineSubtr"));
    h1D_BaselineSubtr->Add(hBaseline, -1.);
  } else {
    BaselineData = CalculateBaseline(h1D_Norm, kTRUE, kFALSE);
    BaselineDataErr = CalculateBaselineError(h1D_Norm, kTRUE, kFALSE);
    for (int iBin = 0; iBin < hBaseline->GetNbinsX(); iBin++) {
      hBaseline->SetBinContent(iBin + 1, BaselineData);
      hBaseline->SetBinError(iBin + 1, BaselineDataErr);
    }
    h1D_BaselineSubtr = reinterpret_cast<TH1D*>(h1D_Norm->Clone("h1D_BaselineSubtr"));
    h1D_BaselineSubtr->Add(hBaseline, -1.);
  }

  h1D_BaselineSubtr->SetMarkerColor(kOrange + 8);
  h1D_BaselineSubtr->SetLineColor(kOrange + 8);
  h1D_BaselineSubtr->GetYaxis()->SetRangeUser(-0.2, 8.);
  fCorrectedCorrHisto_BaselineSubtr = reinterpret_cast<TH1D*>(h1D_BaselineSubtr->Clone(Form("hCorrectedCorrBaselineSubtr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f_InvMass%.0fto%.0f", fPtCandBins[0], fPtCandBins[1], fPtHadBins[0], fPtHadBins[1], fInvMassBins[0], fInvMassBins[1])));

  hBaseline->SetMarkerColor(kPink - 6);
  hBaseline->SetMarkerStyle(kFullSquare);
  hBaseline->SetLineColor(kPink - 6);
  fBaselineHisto = reinterpret_cast<TH1D*>(hBaseline->Clone(Form("hBaseline_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f_InvMass%.0fto%.0f", fPtCandBins[0], fPtCandBins[1], fPtHadBins[0], fPtHadBins[1], fInvMassBins[0], fInvMassBins[1])));

  // clean up
  delete hBaseline;
  hBaseline = nullptr;
  delete h1D_BaselineSubtr;
  h1D_BaselineSubtr = nullptr;

  if (fdoSecPartContamination) {
    h1D_ReflCorr = ReflectCorrHistogram(h1D_Norm_SecPart);
  } else {
    h1D_ReflCorr = ReflectCorrHistogram(h1D_Norm);
  }

  /* used as control using Run2 reflection function
  if (fFDsubtraction) {
    h1D_ReflCorr = ReflectHistoRun2(h1D_SubtrFDNorm, 0.5);
  } else if (fdoSecPartContamination) {
    h1D_ReflCorr = ReflectHistoRun2(h1D_SubtrNorm_SecPart, 0.5);
  } else {
    h1D_ReflCorr = ReflectHistoRun2(h1D_SubtrNorm, 0.5);
  }*/

  h1D_ReflCorr->SetMinimum(0);
  SetTH1HistoStyle(h1D_ReflCorr, Form("%.0f < p_{T} < %.0f GeV/c", fPtCandBins[0], fPtCandBins[1]), "#Delta#phi [rad]", "#frac{1}{N_{D}}#frac{dN^{assoc}}{d#Delta#phi} [rad^{-1}]", kFullCircle, kOrange + 8, 1.6, kOrange + 8, 3);
  fCorrectedCorrHisto_Reflected = reinterpret_cast<TH1D*>(h1D_ReflCorr->Clone(Form("hCorrectedCorrReflected_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f_InvMass%.0fto%.0f", fPtCandBins[0], fPtCandBins[1], fPtHadBins[0], fPtHadBins[1], fInvMassBins[0], fInvMassBins[1])));

  // Reflected histograms baseline subtracted
  TH1D* hBaseline_Refl = reinterpret_cast<TH1D*>(h1D_ReflCorr->Clone("hBaseline_Refl"));
  hBaseline_Refl->Sumw2();
  BaselineData = CalculateBaseline(h1D_ReflCorr, kFALSE, kTRUE);
  BaselineDataErr = CalculateBaselineError(h1D_ReflCorr, kFALSE, kTRUE);

  for (int iBin = 0; iBin < hBaseline_Refl->GetNbinsX(); iBin++) {
    hBaseline_Refl->SetBinContent(iBin + 1, BaselineData);
    hBaseline_Refl->SetBinError(iBin + 1, BaselineDataErr);
  }
  h1D_ReflCorr_BaselineSubtr = reinterpret_cast<TH1D*>(h1D_ReflCorr->Clone("h1D_ReflCorr_BaselineSubtr"));
  h1D_ReflCorr_BaselineSubtr->Sumw2();
  h1D_ReflCorr_BaselineSubtr->Add(hBaseline_Refl, -1.);

  TF1* fConstZero = new TF1("fConstZero", "[0]", 0., TMath::Pi());
  fConstZero->SetParameter(0, 0.);
  fConstZero->SetLineColor(kMagenta);
  fConstZero->SetLineStyle(9);
  fConstZero->SetLineWidth(4);
  fConstZero->SetTitle("");

  hBaseline_Refl->SetMarkerColor(kOrange);
  hBaseline_Refl->SetMarkerStyle(kFullSquare);
  hBaseline_Refl->SetLineColor(kOrange);
  fBaselineHisto_Reflected = reinterpret_cast<TH1D*>(hBaseline_Refl->Clone(Form("hBaselineReflected_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f_InvMass%.0fto%.0f", fPtCandBins[0], fPtCandBins[1], fPtHadBins[0], fPtHadBins[1], fInvMassBins[0], fInvMassBins[1])));

  SetTH1HistoStyle(h1D_ReflCorr_BaselineSubtr, Form("%.0f < p_{T} < %.0f GeV/c", fPtCandBins[0], fPtCandBins[1]), "#Delta#phi [rad]", "#frac{1}{N_{D}}#frac{dN^{assoc}}{d#Delta#phi} [rad^{-1}]", kFullCircle, kRed + 1, 1.6, kRed + 1, 3);
  h1D_ReflCorr_BaselineSubtr->SetStats(0);
  fCorrectedCorrHisto_Reflected_BaselineSubtr = reinterpret_cast<TH1D*>(h1D_ReflCorr_BaselineSubtr->Clone(Form("hCorrectedCorrReflected_BaselineSubtr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f_InvMass%.0fto%.0f", fPtCandBins[0], fPtCandBins[1], fPtHadBins[0], fPtHadBins[1],fInvMassBins[0], fInvMassBins[1])));
  fTFConstZero = reinterpret_cast<TF1*>(fConstZero->Clone(Form("fConstZero_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f_InvMass%.0fto%.0f", fPtCandBins[0], fPtCandBins[1], fPtHadBins[0], fPtHadBins[1], fInvMassBins[0], fInvMassBins[1])));

  // clean up
  delete h1D_Norm;
  h1D_Norm = nullptr;
  if (fdoSecPartContamination) {
    delete h1D_Norm_SecPart;
    h1D_Norm_SecPart = nullptr;
  }
  delete h1D_ReflCorr;
  h1D_ReflCorr = nullptr;
  delete hBaseline_Refl;
  hBaseline_Refl = nullptr;

  return kTRUE;
}

Bool_t DhCorrelationExtraction::ReadInputSEandME()
{

  fFileSE = TFile::Open(fFileNameSE.Data());
  if (!fFileSE) {
    std::cerr << "[ERROR] File " << fFileNameSE << " cannot be opened! check your file path!";
    return kFALSE;
  }

  fFileME = TFile::Open(fFileNameME.Data());
  if (!fFileME) {
    std::cerr << "[ERROR] File " << fFileNameME << " cannot be opened! check your file path!";
    return kFALSE;
  }

  fDirSE = reinterpret_cast<TDirectoryFile*>(fFileSE->Get(fDirNameSE.Data()));
  fDirME = reinterpret_cast<TDirectoryFile*>(fFileME->Get(fDirNameME.Data()));

  std::cout << "===================== " << std::endl;
  std::cout << "Read inputs SE and ME" << std::endl;
  std::cout << "TFile SE    = " << fFileNameSE << std::endl;
  std::cout << "TFile ME    = " << fFileNameME << std::endl;
  std::cout << "===================== " << std::endl;
  std::cout << " " << std::endl;

  return kTRUE;
}

Bool_t DhCorrelationExtraction::ReadInputSecondaryPartContamination()
{

  fFileSecPart = TFile::Open(fFileSecPartName.Data());
  if (!fFileSecPart) {
    std::cerr << "[ERROR] File " << fFileSecPartName << " cannot be opened! check your file path!" << std::endl;
    return kFALSE;
  }

  fDirSecPart = reinterpret_cast<TDirectoryFile*>(fFileSecPart->Get(fDirSecPartName.Data()));

  if (!fDirSecPart) {
    std::cerr << "[ERROR] Directory " << fDirSecPart << " cannot be opened! check your file path!" << std::endl;
    return kFALSE;
  }

  std::cout << "===================== " << std::endl;
  std::cout << "Read inputs SE and ME" << std::endl;
  std::cout << "TFile Sec. part.    = " << fFileSecPartName << std::endl;
  std::cout << "===================== " << std::endl;
  std::cout << " " << std::endl;

  return kTRUE;
}

TH2D* DhCorrelationExtraction::ProjCorrelHisto(Int_t SEorME, Int_t pool)
{
  // TODO: Subtraction of softpion
  TH2D* h2D = new TH2D(); // pointer to be returned
  TH2D* h2DOrignal = nullptr;

  THnSparseF* hSparse = 0x0;
  if (SEorME == kSE) { // Same Event
      hSparse = reinterpret_cast<THnSparseF*>(fDirSE->Get(fCorrelSparseNameSE.Data()));
  } else { // Mixed Event
      hSparse = reinterpret_cast<THnSparseF*>(fDirME->Get(fCorrelSparseNameME.Data()));
  }

  // Check pointer
  if (!hSparse) {
    std::cerr << "[ERROR] hSparse is null! Check that the object name exists in the directory and the file is open." << std::endl;
    throw std::runtime_error("hSparse is null");
  }

  Int_t binExtPoolMin;
  Int_t binExtPoolMax;
  if (fDoPoolByPool) {
    binExtPoolMin = (Int_t)hSparse->GetAxis(kPool)->FindBin(pool + 0.01); // axis1: pool bin
    binExtPoolMax = (Int_t)hSparse->GetAxis(kPool)->FindBin(pool + 0.99);
  } else { // merge all pools in one
    binExtPoolMin = 1;
    binExtPoolMax = (Int_t)hSparse->GetAxis(kPool)->GetNbins();
  }

  if (fDeltaEtaLeftMin < hSparse->GetAxis(kDeltaEta)->GetXmin()) fDeltaEtaLeftMin = hSparse->GetAxis(kDeltaEta)->GetXmin();
  if (fDeltaEtaRightMax > hSparse->GetAxis(kDeltaEta)->GetXmax()) fDeltaEtaRightMax = hSparse->GetAxis(kDeltaEta)->GetXmax();

  // debug: get original histogram before any range is set
  if (fDebug > 1) {
    h2DOrignal = reinterpret_cast<TH2D*>(hSparse->Projection(kDeltaPhi, kDeltaEta)); // axis4: deltaPhi, axis3: deltaEta
    h2DOrignal->SetDirectory(0);
    h2DOrignal->SetName(Form("hOrignal_Correl_%s_Pool%s_2D", (SEorME == kSE) ? "SE" : "ME", fDoPoolByPool ? Form("%d", pool) : "All"));
    h2DOrignal->SetTitle(Form("Original Correlation %s for Pool %s", (SEorME == kSE) ? "SE" : "ME", fDoPoolByPool ? Form("%d", pool) : "All"));
    if (SEorME == kSE) {
      fVecOriginalCorrel_SE_2D.push_back(reinterpret_cast<TH2D*>(h2DOrignal->Clone(Form("hOrignal_Correl_SE_Pool%s_2D", fDoPoolByPool ? Form("%d", pool) : "All"))));
    } else {
      fVecOriginalCorrel_ME_2D.push_back(reinterpret_cast<TH2D*>(h2DOrignal->Clone(Form("hOrignal_Correl_ME_Pool%s_2D", fDoPoolByPool ? Form("%d", pool) : "All"))));
    }
    delete h2DOrignal;
    h2DOrignal = nullptr;
  }

  // set ranges
  hSparse->GetAxis(kPool)->SetRangeUser(pool+0.01, fDoPoolByPool ? pool+0.99 : hSparse->GetAxis(kPool)->GetXmax()); // axis0: pool bin
  hSparse->GetAxis(kPtCand)->SetRangeUser(fPtCandBins[0], fPtCandBins[1]);     // axis1: ptCand
  hSparse->GetAxis(kPtHad)->SetRangeUser(fPtHadBins[0], fPtHadBins[1]);       // axis2: ptHad
  if (SEorME == kSE) { // Same Event
    TH1D* temp = reinterpret_cast<TH1D*>(hSparse->Projection(kPtCand)); // axis1: ptCand
    fNpairs = temp->GetEntries();
    // fNpairsErr = TMath::Sqrt(temp->GetEntries());
    delete temp;
    temp = nullptr;
  }

  // mass selecton for kMassBinning method, but whole range for kDeltaPhiBinning method
  hSparse->GetAxis(kMass)->SetRangeUser(fInvMassBins[0], fInvMassBins[1]); // axis5: invMass
  hSparse->GetAxis(kDeltaEta)->SetRangeUser(fDeltaEtaLeftMin+0.01, fDeltaEtaRightMax-0.01); // axis3: deltaEta
  TH2D* hFinal = (TH2D*)hSparse->Projection(kDeltaPhi, kDeltaEta);            // axis4: deltaPhi, axis3: deltaEta
  if (SEorME == kME) CalculateNormaliztionFactorME(hFinal, pool);
  if(fMethod == kDeltaPhiBinning) hSparse->GetAxis(kDeltaPhi)->SetRangeUser(fDeltaPhiBins[0], fDeltaPhiBins[1]); // axis4: deltaPhi
  TH2D* hFinalMass = (TH2D*)hSparse->Projection(kMass, kDeltaEta);            // axis5: invMass, axis3: deltaEta

  // set to 0 for the inner deltaEta gap
  Int_t leftBinMax = hFinal->GetXaxis()->FindBin(fDeltaEtaLeftMax - 0.01);
  Int_t rightBinMin = hFinal->GetXaxis()->FindBin(fDeltaEtaRightMin + 0.01);
  Int_t nBinsY = hFinal->GetNbinsY();
  bool doFinalMass = (fMethod == kDeltaPhiBinning && SEorME == kSE);
  Int_t nBinsYMass = doFinalMass ? hFinalMass->GetNbinsY() : 0;
  for (Int_t ix = leftBinMax + 1; ix < rightBinMin; ++ix) {
      for (Int_t iy = 1; iy <= nBinsY; ++iy) {
          hFinal->SetBinContent(ix, iy, 0);
          hFinal->SetBinError(ix, iy, 0);
      }
      if (doFinalMass) {
          for (Int_t iy = 1; iy <= nBinsYMass; ++iy) {
              hFinalMass->SetBinContent(ix, iy, 0);
              hFinalMass->SetBinError(ix, iy, 0);
          }
      }
  }

  h2D = (TH2D*)hFinal->Clone();
  delete hFinal;
  hFinal = nullptr;
  if (fMethod == kDeltaPhiBinning && SEorME == kSE) {
    fVecMassPairsVsDeltaEta_2D.push_back(reinterpret_cast<TH2D*>(hFinalMass->Clone(Form("hMassPairsVsEta_%s_Pool%s_2D", (SEorME == kSE) ? "SE" : "ME", fDoPoolByPool ? Form("%d", pool) : "All"))));
    delete hFinalMass;
    hFinalMass = nullptr;
  }

  // project to 2D histogram
  // h2D = reinterpret_cast<TH2D*>(hSparse->Projection(4, 3));            // axis4: deltaPhi, axis3: deltaEta
  h2D->SetDirectory(0);
  h2D->SetName(Form("hCorrel_%s_Pool%s_2D", (SEorME == kSE) ? "SE" : "ME", fDoPoolByPool ? Form("%d", pool) : "All"));
  h2D->SetTitle(Form("Correlation %s for Pool %s", (SEorME == kSE) ? "SE" : "ME", fDoPoolByPool ? Form("%d", pool) : "All"));

  // debug: projected histogram after ranges are set
  if (fDebug > 0) {
    if (SEorME == kSE) {
      fVecCorrel_SE_2D.push_back(reinterpret_cast<TH2D*>(h2D->Clone(Form("hCorrel_SE_Pool%s_2D", fDoPoolByPool ? Form("%d", pool) : "All"))));
    } else {
      fVecCorrel_ME_2D.push_back(reinterpret_cast<TH2D*>(h2D->Clone(Form("hCorrel_ME_Pool%s_2D", fDoPoolByPool ? Form("%d", pool) : "All"))));
    }
  }

  // clean up
  delete hSparse;
  hSparse = nullptr;

  return h2D;
}

void DhCorrelationExtraction::CalculateNormaliztionFactorME(TH2D* histoME, Int_t pool)
{

  Int_t bin0phi = histoME->GetYaxis()->FindBin(0.);
  Int_t bin0eta = histoME->GetXaxis()->FindBin(0.);

  // evaluate the normalization (from ALL tracks, including possible fake softpions) -> **histoME indeed includes bin1+bin2 of THnSparse, i.e. all the tracks**
  Double_t factorNorm = 0;
  for (int in = -1; in <= 0; in++) {
    factorNorm += histoME->GetBinContent(bin0eta, bin0phi + in);
  }
  for (int in = -1; in <= 0; in++) {
    factorNorm += histoME->GetBinContent(bin0eta - 1, bin0phi + in);
  }
  factorNorm /= 4.;

  if (factorNorm == 0) {
    throw std::runtime_error("Normalization factor is zero (Only bin (0,0) is used for normalization)");
  }

  fFactorsNormME[pool] = factorNorm;
}

TH1D* DhCorrelationExtraction::ProjCorrelHistoSecondaryPart(Int_t PartType, Double_t PtCandMin, Double_t PtCandMax, Double_t PtHadMin, Double_t PtHadMax)
{

  TH1D* h1D = new TH1D(); // pointer to be returned
  TH1D* h1DOrig = nullptr;

  THnSparseD* hSparse = 0x0;

  if (PartType == kPrimaryPart) { // primary particles
    hSparse = reinterpret_cast<THnSparseD*>(fDirSecPart->Get(fHistoPrimaryPartName.Data()));
  } else { // all selected particles
    hSparse = reinterpret_cast<THnSparseD*>(fDirSecPart->Get(fHistoAllPartName.Data()));
  }

  // Check pointer
  if (!hSparse) {
    std::cerr << "[ERROR] hSparse is null! Check that the object name exists in the directory and the file is open." << std::endl;
    throw std::runtime_error("hSparse is null");
  }

  // get bin ranges
  Int_t binExtPtCandMin = (Int_t)hSparse->GetAxis(2)->FindBin(PtCandMin + 0.01); // axis2: ptCand, the 0.01 to avoid bin edges!
  Int_t binExtPtCandMax = (Int_t)hSparse->GetAxis(2)->FindBin(PtCandMax - 0.01);
  Int_t binExtPtHadMin = (Int_t)hSparse->GetAxis(3)->FindBin(PtHadMin + 0.01); // axis3: ptHad
  Int_t binExtPtHadMax = (Int_t)hSparse->GetAxis(3)->FindBin(PtHadMax - 0.01);
  Int_t binExtPoolMin;
  Int_t binExtPoolMax;
  if (PartType == kAllPart) {
    binExtPoolMin = 1;
    binExtPoolMax = (Int_t)hSparse->GetAxis(4)->GetNbins();
  }
  // possibility to select a certain eta region
  Int_t binExtEtaLeftMin = (Int_t)hSparse->GetAxis(1)->FindBin(fDeltaEtaLeftMin + 0.0001);
  Int_t binExtEtaLeftMax = (Int_t)hSparse->GetAxis(1)->FindBin(fDeltaEtaLeftMax - 0.0001);
  if (binExtEtaLeftMax > hSparse->GetAxis(1)->GetNbins())
    binExtEtaLeftMax = hSparse->GetAxis(1)->GetNbins();
  if (binExtEtaLeftMin < 1)
    binExtEtaLeftMin = 1;
  Int_t binExtEtaRightMin = (Int_t)hSparse->GetAxis(1)->FindBin(fDeltaEtaRightMin + 0.0001);
  Int_t binExtEtaRightMax = (Int_t)hSparse->GetAxis(1)->FindBin(fDeltaEtaRightMax - 0.0001);
  if (binExtEtaRightMax > hSparse->GetAxis(1)->GetNbins())
    binExtEtaRightMax = hSparse->GetAxis(1)->GetNbins();
  if (binExtEtaRightMin < 1)
    binExtEtaRightMin = 1;

  if (fDebug > 1) {
    if (PartType == kPrimaryPart) {
      h1DOrig = reinterpret_cast<TH1D*>(hSparse->Projection(0)); // axis0: deltaPhi
      h1DOrig->SetName("hPrimaryPartCorr_Orig");
      fOriginalCorrel_PrimaryPart = reinterpret_cast<TH1D*>(h1DOrig->Clone("hPrimaryPartCorr_Orig"));
      h1DOrig = nullptr;
    } else {
      h1DOrig = reinterpret_cast<TH1D*>(hSparse->Projection(0)); // axis0: deltaPhi
      h1DOrig->SetName("hAllPartCorr_Orig");
      fOriginalCorrel_AllPart = reinterpret_cast<TH1D*>(h1DOrig->Clone("hAllPartCorr_Orig"));
      h1DOrig = nullptr;
    }
  }

  // set ranges
  hSparse->GetAxis(1)->SetRange(binExtEtaLeftMin, binExtEtaLeftMax);       // axis1: deltaEta
  hSparse->GetAxis(1)->SetRange(binExtEtaRightMin, binExtEtaRightMax); // axis1: deltaEta
  hSparse->GetAxis(2)->SetRange(binExtPtCandMin, binExtPtCandMax); // axis2: ptCand
  hSparse->GetAxis(3)->SetRange(binExtPtHadMin, binExtPtHadMax);   // axis3: ptHad
  if (PartType == kAllPart) {
    hSparse->GetAxis(4)->SetRange(binExtPoolMin, binExtPoolMax); // axis4: pool bin
  }

  h1D = reinterpret_cast<TH1D*>(hSparse->Projection(0)); // axis0: deltaPhi
  if (PartType == kPrimaryPart) {                        // primary particles
    h1D->SetName(Form("hPrimaryPartCorr_PtD%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  } else { // all selected particles
    h1D->SetName(Form("hAllPartCorr_PtD%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  }

  if (fDebug > 0) {
    if (PartType == kPrimaryPart) {
      fCorrel_PrimaryPart = reinterpret_cast<TH1D*>(h1D->Clone(Form("hPrimaryPartCorr_PtD%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax)));
    } else {
      fCorrel_AllPart = reinterpret_cast<TH1D*>(h1D->Clone(Form("hAllPartCorr_PtD%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax)));
    }
  }

  // clean up
  delete hSparse;
  hSparse = nullptr;

  return h1D;
}

// load and project mass THnSparse to get mass vs pt histogram
void DhCorrelationExtraction::ProjMassVsPt() {
  TFile* fileMass = TFile::Open(fFileNameMass.Data());
  if (!fileMass || fileMass->IsZombie()) {
    std::cerr << "[ERROR] Could not open file: " << fFileNameMass.Data() << std::endl;
    return;
  }

  THnSparseF* sparseMass = reinterpret_cast<THnSparseF*>(fileMass->Get(fMassSparseName.Data()));
  TH2D* hMassVsPt = reinterpret_cast<TH2D*>(sparseMass->Projection(1,0));
  hMassVsPt -> SetDirectory(0);
  hMassVsPt -> SetName("hMassVsPt");
  fMassVsPt_2D = hMassVsPt;

  fileMass -> Close();
  delete sparseMass;
}

// calculate the number of trigger D mesons in given pt range from mass vs pt histogram
Double_t DhCorrelationExtraction::CalculateTriggerNormalizationFactor(TH2D* hMassVsPt, Double_t ptMin, Double_t ptMax, Double_t massMin, Double_t massMax) {
  if (!hMassVsPt) {
    std::cerr << "[ERROR] hMassVsPt is null!" << std::endl;
    return 0.0;
  }

  Int_t ptBinMin = hMassVsPt->GetYaxis()->FindBin(ptMin + 1e-6);
  Int_t ptBinMax = hMassVsPt->GetYaxis()->FindBin(ptMax - 1e-6);

  TH1D* hMassProj = hMassVsPt->ProjectionX("hMassProj", ptBinMin, ptBinMax);

  Int_t lowEdgeBin = hMassProj->GetXaxis()->FindBin(1.72);
  Int_t upEdgeBin = hMassProj->GetXaxis()->FindBin(2.02);
  Double_t Ntrigegr = hMassProj->Integral(lowEdgeBin, upEdgeBin);

  Int_t massBinMin = hMassProj->GetXaxis()->FindBin(massMin);
  Int_t massBinMax = hMassProj->GetXaxis()->FindBin(massMax);

  // Double_t triggerCount = hMassProj->Integral(massBinMin, massBinMax) * fNpairs / Ntrigger;
  Double_t triggerCount = hMassProj->Integral(massBinMin, massBinMax);
  delete hMassProj;

  return triggerCount;
}

TH1D* DhCorrelationExtraction::ReflectCorrHistogram(TH1D*& histo)
{

  // nBinsPhi must be a multple of 4 in order to reflect correcty the histogram
  Int_t nBinsPhi = histo->GetNbinsX();
  Int_t nBinsPhiRefl = nBinsPhi / 2;
  Int_t bin0Phi = nBinsPhi / 4 + 1;
  Int_t binPiPhi = 3 * nBinsPhi / 4;

  TH1D* h1D = new TH1D("h1D_Reflected", "", nBinsPhiRefl, 0., TMath::Pi()); // pointer to be returned
  h1D->Sumw2();
  // TH1D* h1D = reinterpret_cast<TH1D*> histo -> Clone("h1D_Reflected");
  // h1D -> GetXaxis() -> SetRange(bin0Phi, binPiPhi);

  // reflection
  Double_t reflectedContent, reflectedContentError;
  for (int iBin = 0; iBin < nBinsPhiRefl / 2; iBin++) {
    reflectedContent = (histo->GetBinContent(bin0Phi - iBin - 1) + histo->GetBinContent(bin0Phi + iBin)) / 2;
    reflectedContentError = 0.5 * TMath::Sqrt(TMath::Power(histo->GetBinError(iBin + 1), 2) + TMath::Power(histo->GetBinError(bin0Phi + iBin), 2));
    h1D->SetBinContent(iBin + 1, reflectedContent);
    h1D->SetBinError(iBin + 1, reflectedContentError);
  }
  for (int iBin = nBinsPhiRefl / 2; iBin < nBinsPhiRefl; iBin++) {
    reflectedContent = (histo->GetBinContent(bin0Phi + iBin) + histo->GetBinContent(binPiPhi + 2 * bin0Phi - iBin - 2)) / 2;
    reflectedContentError = 0.5 * TMath::Sqrt(TMath::Power(histo->GetBinError(bin0Phi + iBin), 2) + TMath::Power(histo->GetBinError(binPiPhi + 2 * bin0Phi - iBin - 2), 2));
    h1D->SetBinContent(iBin + 1, reflectedContent);
    h1D->SetBinError(iBin + 1, reflectedContentError);
  }

  return h1D;
}

TH1D* DhCorrelationExtraction::ReflectHistoRun2(TH1D* h, Double_t scale)
{

  TH1D* h2 = new TH1D(Form("%sReflected", h->GetName()), Form("%sReflected", h->GetName()), h->GetNbinsX() / 2., 0., TMath::Pi());
  for (Int_t j = 1; j <= h->GetNbinsX(); j++) {
    Double_t x = h->GetBinCenter(j);
    Double_t y0 = h->GetBinContent(j);
    Double_t ey0 = h->GetBinError(j);
    Int_t j2;
    if (x > 0 && x < TMath::Pi()) {
      j2 = h2->FindBin(x);
    } else if (x < 0) {
      j2 = h2->FindBin(-1. * x);
    } else if (x > TMath::Pi()) {
      j2 = h2->FindBin(2. * TMath::Pi() - x);
    } else {
      printf("Point %d excluded \n", j);
      continue;
    }
    Double_t y = h2->GetBinContent(j2);
    Double_t ey = h2->GetBinError(j2);
    h2->SetBinContent(j2, (y + y0));
    h2->SetBinError(j2, TMath::Sqrt(ey0 * ey0 + ey * ey));
  }
  h2->Scale(scale);

  return h2;
}

void DhCorrelationExtraction::NormalizeMEplot(TH2D*& histoME, TH2D*& histoMEsoftPi, Int_t pool)
{
  // apply the normalization
  histoME->Scale(1. / fFactorsNormME[pool]);
  histoME->SetTitle(Form("ME normalized to %.2f", fFactorsNormME[pool]));
  return;
}

Double_t DhCorrelationExtraction::CalculateBaseline(TH1D*& histo, Bool_t totalRange, Bool_t reflected)
{

  // total range = 2*Pi
  // half range = Pi , for histogram reflected under symmetric assumption

  Double_t baseline, errBaseline;
  Int_t nBinsPhi = histo->GetNbinsX();
  Int_t binPhiHalf = nBinsPhi / 2;
  Int_t binPhiHalfMinus1 = nBinsPhi / 2 - 1;
  Int_t binPhiHalfPlus1 = nBinsPhi / 2 + 1;
  Int_t binPhiHalfPlus2 = nBinsPhi / 2 + 1;

  if (totalRange) {
    // baseline evaluated considering: the two first points, the last two points and four points in the middle (corresponding to the outer points)
    if (nBinsPhi >= 32) {
      baseline =
        ((histo->GetBinContent(1)) * (1. / TMath::Power(histo->GetBinError(1), 2)) +
         (histo->GetBinContent(2)) * (1. / TMath::Power(histo->GetBinError(2), 2)) +
         (histo->GetBinContent(binPhiHalfMinus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
         (histo->GetBinContent(binPhiHalf)) * (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
         (histo->GetBinContent(binPhiHalfPlus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
         (histo->GetBinContent(binPhiHalfPlus2)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)) +
         (histo->GetBinContent(nBinsPhi - 1)) * (1. / TMath::Power(histo->GetBinError(nBinsPhi - 1), 2)) +
         (histo->GetBinContent(nBinsPhi)) * (1. / TMath::Power(histo->GetBinError(nBinsPhi), 2))) /
        ((1. / TMath::Power(histo->GetBinError(1), 2)) +
         (1. / TMath::Power(histo->GetBinError(2), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)) +
         (1. / TMath::Power(histo->GetBinError(nBinsPhi - 1), 2)) +
         (1. / TMath::Power(histo->GetBinError(nBinsPhi), 2)));
    } else {
      baseline =
        ((histo->GetBinContent(1)) * (1. / TMath::Power(histo->GetBinError(1), 2)) +
         (histo->GetBinContent(binPhiHalf)) * (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
         (histo->GetBinContent(binPhiHalfPlus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
         (histo->GetBinContent(nBinsPhi)) * (1. / TMath::Power(histo->GetBinError(nBinsPhi), 2))) /
        ((1. / TMath::Power(histo->GetBinError(1), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
         (1. / TMath::Power(histo->GetBinError(nBinsPhi), 2)));
    }
  } else {
    if (reflected) {
      baseline =
        ((histo->GetBinContent(binPhiHalfMinus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
         (histo->GetBinContent(binPhiHalf)) * (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
         (histo->GetBinContent(binPhiHalfPlus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
         (histo->GetBinContent(binPhiHalfPlus2)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2))) /
        ((1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)));
    } else {
      // baseline evaluated using the 4 middle points in the transverese region
      if (nBinsPhi >= 32) {
        baseline =
          ((histo->GetBinContent(binPhiHalfMinus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
           (histo->GetBinContent(binPhiHalf)) * (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
           (histo->GetBinContent(binPhiHalfPlus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
           (histo->GetBinContent(binPhiHalfPlus2)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2))) /
          ((1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
           (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
           (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
           (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)));
      } else {
        baseline =
          ((histo->GetBinContent(binPhiHalf)) * (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
           (histo->GetBinContent(binPhiHalfPlus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2))) /
          ((1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
           (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)));
      }
    }
  }

  return baseline;
}

Double_t DhCorrelationExtraction::CalculateBaselineError(TH1D*& histo, Bool_t totalRange, Bool_t reflected)
{

  // total range = 2*Pi
  // half range = Pi , for histogram reflected under symmetric assumption

  Double_t errBaseline;
  Int_t nBinsPhi = histo->GetNbinsX();
  Int_t binPhiHalf = nBinsPhi / 2;
  Int_t binPhiHalfMinus1 = nBinsPhi / 2 - 1;
  Int_t binPhiHalfPlus1 = nBinsPhi / 2 + 1;
  Int_t binPhiHalfPlus2 = nBinsPhi / 2 + 1;

  if (totalRange) {
    // baseline evaluated considering: the two first points, the last two points and four points in the middle (corresponding to the outer points)
    if (nBinsPhi >= 32) {
      errBaseline = 1. /
                    TMath::Sqrt((1. / TMath::Power(histo->GetBinError(1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(2), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)) +
                                (1. / TMath::Power(histo->GetBinError(nBinsPhi - 1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(nBinsPhi), 2)));
    } else { // fon nBinsPhi = 16 (rebin 4)
      errBaseline = 1. /
                    TMath::Sqrt((1. / TMath::Power(histo->GetBinError(1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(nBinsPhi), 2)));
    }
  } else {
    // baseline evaluated using the 4 middle points in the transverese region
    if (reflected) {
      errBaseline = 1. /
                    TMath::Sqrt((1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)));
    } else {
      if (nBinsPhi >= 32) {
        errBaseline = 1. /
                      TMath::Sqrt((1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
                                  (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
                                  (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
                                  (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)));
      } else {
        errBaseline = 1. /
                      TMath::Sqrt((1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
                                  (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)));
      }
    }
  }

  return errBaseline;
}

void DhCorrelationExtraction::SetTH1HistoStyle(TH1D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle,
                                               Style_t markerStyle, Color_t markerColor, Double_t markerSize,
                                               Color_t lineColor, Int_t lineWidth, Float_t hTitleXaxisOffset, Float_t hTitleYaxisOffset,
                                               Float_t hTitleXaxisSize, Float_t hTitleYaxisSize, Float_t hLabelXaxisSize, Float_t hLabelYaxisSize,
                                               Bool_t centerXaxisTitle, Bool_t centerYaxisTitle)
{

  histo->SetTitle(hTitle.Data());
  histo->GetXaxis()->SetTitle(hXaxisTitle.Data());
  histo->GetYaxis()->SetTitle(hYaxisTitle.Data());
  histo->SetMarkerStyle(markerStyle);
  histo->SetMarkerColor(markerColor);
  histo->SetMarkerSize(markerSize);
  histo->SetLineColor(lineColor);
  histo->SetLineWidth(lineWidth);
  histo->GetXaxis()->SetTitleOffset(hTitleXaxisOffset);
  histo->GetYaxis()->SetTitleOffset(hTitleYaxisOffset);
  histo->GetXaxis()->SetTitleSize(hTitleXaxisSize);
  histo->GetYaxis()->SetTitleSize(hTitleYaxisSize);
  histo->GetXaxis()->SetLabelSize(hLabelXaxisSize);
  histo->GetYaxis()->SetLabelSize(hLabelYaxisSize);
  histo->GetXaxis()->CenterTitle(centerXaxisTitle);
  histo->GetYaxis()->CenterTitle(centerYaxisTitle);

  return;
}

void DhCorrelationExtraction::SetTH2HistoStyle(TH2D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle, TString hZaxisTitle,
                                               Float_t hTitleXaxisOffset, Float_t hTitleYaxisOffset, Float_t hTitleZaxisOffset,
                                               Float_t hTitleXaxisSize, Float_t hTitleYaxisSize, Float_t hTitleZaxisSize,
                                               Float_t hLabelXaxisSize, Float_t hLabelYaxisSize, Float_t hLabelZaxisSize,
                                               Bool_t centerXaxisTitle, Bool_t centerYaxisTitle)
{

  histo->SetTitle(hTitle.Data());
  histo->GetXaxis()->SetTitle(hXaxisTitle.Data());
  histo->GetYaxis()->SetTitle(hYaxisTitle.Data());
  histo->GetZaxis()->SetTitle(hZaxisTitle.Data());
  histo->GetXaxis()->SetTitleOffset(hTitleXaxisOffset);
  histo->GetYaxis()->SetTitleOffset(hTitleYaxisOffset);
  histo->GetZaxis()->SetTitleOffset(hTitleZaxisOffset);
  histo->GetXaxis()->SetTitleSize(hTitleXaxisSize);
  histo->GetYaxis()->SetTitleSize(hTitleYaxisSize);
  histo->GetZaxis()->SetTitleSize(hTitleZaxisSize);
  histo->GetXaxis()->SetLabelSize(hLabelXaxisSize);
  histo->GetYaxis()->SetLabelSize(hLabelYaxisSize);
  histo->GetZaxis()->SetLabelSize(hLabelZaxisSize);
  histo->GetXaxis()->CenterTitle(centerXaxisTitle);
  histo->GetYaxis()->CenterTitle(centerYaxisTitle);

  return;
}
