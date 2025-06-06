#include "VnVsMassFitter.h"

#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TColor.h>
#include <TLegend.h>
#include <TList.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <TDatabasePDG.h>
#include <TPaveText.h>
#include "Fit/BinData.h"
#include "HFitInterface.h"
#include <vector>

/// \cond CLASSIMP
ClassImp(VnVsMassFitter);
/// \endcond

//________________________________________________________________
VnVsMassFitter::VnVsMassFitter()
  :TObject()
  ,fMassHisto(0x0)
  ,fVnVsMassHisto(0x0)
  ,fMassSgnFuncType(kGaus)
  ,fMassBkgFuncType(kExpo)
  ,fVnBkgFuncType(kLin)
  ,fMassFuncFromPrefit(0x0)
  ,fMassBkgFunc(0x0)
  ,fMassSgnFunc(0x0)
  ,fMassTemplFunc(0x0)
  ,fMassTotFunc(0x0)
  ,fVnBkgFuncSb(0x0)
  ,fVnBkgFunc(0x0)
  ,fVnTotFunc(0x0)
  ,fMassFitter(0x0)
  ,fMassMin(1.69)
  ,fMassMax(2.05)
  ,fVn(0.)
  ,fVnUncertainty(0.)
  ,fSigma(0.)
  ,fSigmaUncertainty(0.)
  ,fMean(0.)
  ,fMeanUncertainty(0.)
  ,fRawYield(0.)
  ,fRawYieldUncertainty(0.)
  ,fChiSquare(0.)
  ,fNDF(0)
  ,fProb(0.)
  ,fSBVnPrefitChiSquare(0.)
  ,fSBVnPrefitNDF(0)
  ,fSBVnPrefitProb(0.)
  ,fMassPrefitChiSquare(0.)
  ,fMassPrefitNDF(0)
  ,fMassPrefitProb(0.)
  ,fNSigmaForSB(3.)
  ,fSigmaInit(0.012)
  ,fMeanInit(1.870)
  ,fSigma2GausInit(0.012)
  ,fFrac2GausInit(0.2)
  ,fMeanFixedFromMassFit(kFALSE)
  ,fSigmaFixedFromMassFit(kFALSE)
  ,fSigma2GausFixedFromMassFit(kFALSE)
  ,fFrac2GausFixedFromMassFit(kFALSE)
  ,fMassParticle(1.870)
  ,fNParsMassSgn(3)
  ,fNParsMassBkg(2)
  ,fNParsVnBkg(2)
  ,fNParsVnSgn(1)
  ,fNParsVnSecPeak(0)
  ,fNParsVnRfl(0)
  ,fSigmaFixed(0)
  ,fMeanFixed(0)
  ,fSigma2GausFixed(0)
  ,fFrac2GausFixed(0)
  ,fPolDegreeBkg(3)
  ,fPolDegreeVnBkg(3)
  ,fReflections(kFALSE)
  ,fNParsRfl(0)
  ,fTemplates(kFALSE)
  ,fNParsTempls(0)
  ,fRflOverSig(0.)
  ,fFixRflOverSig(kFALSE)
  ,fHistoTemplRfl(0x0)
  ,fHistoTemplRflInit(0x0)
  ,fMassRflFunc(0x0)
  ,fMassBkgRflFunc(0x0)
  ,fRflOpt("1gaus")
  ,fMinRefl(0.)
  ,fMaxRefl(0.)
  ,fSmoothRfl(kFALSE)
  ,fRawYieldHelp(0.)
  ,fVnRflOpt(0)
  ,fVnRflLimited(kFALSE)
  ,fVnRflMin(-1.)
  ,fVnRflMax(1.)
  ,fSecondPeak(kFALSE)
  ,fMassSecPeakFunc(0x0)
  ,fVnSecPeakFunc(0x0)
  ,fNParsSec(0)
  ,fSecMass(-999.)
  ,fSecWidth(9999.)
  ,fFixSecMass(kFALSE)
  ,fFixSecWidth(kFALSE)
  ,fDoSecondPeakVn(kFALSE)
  ,fFixVnSecPeakToSgn(kFALSE)
  ,fAnchorTemplsMode(Free)
  ,fHarmonic(2) {

    //default constructor
}

//________________________________________________________________
VnVsMassFitter::VnVsMassFitter(TH1F* hMass, TH1F* hvn, Double_t min, Double_t max, Int_t funcMassBkg, Int_t funcMassSgn, Int_t funcvnBkg)
  :TObject()
  ,fMassSgnFuncType(funcMassSgn)
  ,fMassBkgFuncType(funcMassBkg)
  ,fVnBkgFuncType(funcvnBkg)
  ,fMassFuncFromPrefit(0x0)
  ,fMassBkgFunc(0x0)
  ,fMassSgnFunc(0x0)
  ,fMassTemplFunc(0x0)
  ,fMassTotFunc(0x0)
  ,fVnBkgFuncSb(0x0)
  ,fVnBkgFunc(0x0)
  ,fVnTotFunc(0x0)
  ,fMassMin(min)
  ,fMassMax(max)
  ,fVn(0.)
  ,fVnUncertainty(0.)
  ,fSigma(0.)
  ,fSigmaUncertainty(0.)
  ,fMean(0.)
  ,fMeanUncertainty(0.)
  ,fRawYield(0.)
  ,fRawYieldUncertainty(0.)
  ,fChiSquare(0.)
  ,fNDF(0)
  ,fProb(0.)
  ,fSBVnPrefitChiSquare(0.)
  ,fSBVnPrefitNDF(0)
  ,fSBVnPrefitProb(0.)
  ,fMassPrefitChiSquare(0.)
  ,fMassPrefitNDF(0)
  ,fMassPrefitProb(0.)
  ,fNSigmaForSB(3.)
  ,fSigmaInit(0.012)
  ,fMeanInit(1.870)
  ,fSigma2GausInit(0.012)
  ,fFrac2GausInit(0.2)
  ,fMeanFixedFromMassFit(kFALSE)
  ,fSigmaFixedFromMassFit(kFALSE)
  ,fSigma2GausFixedFromMassFit(kFALSE)
  ,fFrac2GausFixedFromMassFit(kFALSE)
  ,fMassParticle(1.870)
  ,fNParsMassSgn(3)
  ,fNParsMassBkg(2)
  ,fNParsVnBkg(2)
  ,fNParsVnSgn(1)
  ,fNParsVnSecPeak(0)
  ,fNParsVnRfl(0)
  ,fSigmaFixed(0)
  ,fMeanFixed(0)
  ,fSigma2GausFixed(0)
  ,fFrac2GausFixed(0)
  ,fPolDegreeBkg(3)
  ,fPolDegreeVnBkg(3)
  ,fReflections(kFALSE)
  ,fNParsRfl(0)
  ,fTemplates(kFALSE)
  ,fNParsTempls(0)
  ,fRflOverSig(0.)
  ,fFixRflOverSig(kFALSE)
  ,fHistoTemplRfl(0x0)
  ,fHistoTemplRflInit(0x0)
  ,fMassRflFunc(0x0)
  ,fMassBkgRflFunc(0x0)
  ,fRflOpt("1gaus")
  ,fMinRefl(0.)
  ,fMaxRefl(0.)
  ,fSmoothRfl(kFALSE)
  ,fRawYieldHelp(0.)
  ,fVnRflOpt(0)
  ,fVnRflLimited(kFALSE)
  ,fVnRflMin(-1.)
  ,fVnRflMax(1.)
  ,fSecondPeak(kFALSE)
  ,fMassSecPeakFunc(0x0)
  ,fVnSecPeakFunc(0x0)
  ,fNParsSec(0)
  ,fSecMass(-999.)
  ,fSecWidth(9999.)
  ,fFixSecMass(kFALSE)
  ,fFixSecWidth(kFALSE)
  ,fDoSecondPeakVn(kFALSE)
  ,fFixVnSecPeakToSgn(kFALSE)
  ,fAnchorTemplsMode(Free)
  ,fHarmonic(2)
  ,fMassVar("mass", "Mass", this->fMassMin, this->fMassMax) {

    //standard constructor
    fMassHisto = (TH1F*)hMass->Clone("fHistoInvMass");
    fMassHisto->SetDirectory(0);
    fVnVsMassHisto = (TH1F*)hvn->Clone(Form("fHistoV%dVsMass",fHarmonic));
    fVnVsMassHisto->SetDirectory(0);
    fMassFitter = new InvMassFitter(fMassHisto,fMassMin,fMassMax,fMassBkgFuncType,fMassSgnFuncType);

    DefineNumberOfParameters();
}

//________________________________________________________________
VnVsMassFitter::~VnVsMassFitter() {

  //destructor
  if(fMassHisto)          delete fMassHisto;
  if(fVnVsMassHisto)      delete fVnVsMassHisto;
  if(fMassBkgFunc)        delete fMassBkgFunc;
  if(fMassBkgRflFunc)     delete fMassBkgRflFunc;
  if(fMassSgnFunc)        delete fMassSgnFunc;
  if(fMassTemplFunc)      delete fMassTemplFunc;
  if(fMassTotFunc)        delete fMassTotFunc;
  if(fVnBkgFuncSb)        delete fVnBkgFuncSb;
  if(fVnBkgFunc)          delete fVnBkgFunc;
  if(fVnTotFunc)          delete fVnTotFunc;
  if(fMassFitter)         delete fMassFitter;
  if(fHistoTemplRfl)      delete fHistoTemplRfl;
  if(fHistoTemplRflInit)  delete fHistoTemplRflInit;
  if(fMassRflFunc)        delete fMassRflFunc;
  if(fMassSecPeakFunc)    delete fMassSecPeakFunc;
  if(fVnSecPeakFunc)      delete fVnSecPeakFunc;
}

//________________________________________________________________
TH1D* VnVsMassFitter::GetPullDistribution() {
  if(!fMassTotFunc) {
      throw std::invalid_argument("Fit not performed, pulls cannot be calculated!");
  }

  std::vector<double> pulls;
  for(int iBin=0; iBin<this->fMassHisto->GetNbinsX(); iBin++) {    
      if(this->fMassHisto->GetBinCenter(iBin+1) >= this->fMassMin &&
         this->fMassHisto->GetBinCenter(iBin+1) <= this->fMassMax) {
              pulls.push_back( (this->fMassHisto->GetBinContent(iBin+1) - this->GetMassTotFitFunc()->Eval(this->fMassHisto->GetBinCenter(iBin+1))) /         
                                this->fMassHisto->GetBinError(iBin+1));
      }
  }

  TH1D *histoPulls = new TH1D("hPulls", "hPulls;M (GeV/c); Data - fit", pulls.size(), this->fMassMin, this->fMassMax);
  for(int iBin=0; iBin<this->fMassHisto->GetNbinsX(); iBin++) {    
      histoPulls->SetBinContent(iBin+1, pulls[iBin]);         
  }

  return histoPulls;
}

//________________________________________________________________
Bool_t VnVsMassFitter::SimultaneousFit() {

  if(!fMassHisto || !fVnVsMassHisto) {printf("Histograms not set! Exit."); return kFALSE;}
  DefineNumberOfParameters();

  const Int_t nparsmass = fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsTempls;
  Int_t NvnParsSgn = 1;
  if(fSecondPeak && fDoSecondPeakVn) {NvnParsSgn+=1;}
  if(fReflections && fVnRflOpt==kFreePar) {NvnParsSgn+=1;}
  Int_t NvnParsTempls = 0;
  if(!fTemplSameVnOfSignal) {NvnParsTempls+=fNParsTempls;}
  const Int_t nparsvn = nparsmass+fNParsVnBkg+NvnParsSgn+NvnParsTempls;

  Bool_t massprefit=MassPrefit();
  if(!massprefit) {printf("Impossible to perform the mass prefit"); return kFALSE;}
  Bool_t vnprefit=VnSBPrefit();
  if(!vnprefit) {printf("Impossible to perform the bkg vn prefit"); return kFALSE;}
  std::vector<Double_t> initpars;
  for(Int_t iBkgPar=0; iBkgPar<fNParsMassBkg; iBkgPar++) {
    initpars.push_back(fMassFuncFromPrefit->GetParameter(iBkgPar));
  }
  for(Int_t iSgnPar=0; iSgnPar<fNParsMassSgn; iSgnPar++) {
    initpars.push_back(fMassFuncFromPrefit->GetParameter(iSgnPar+fNParsMassBkg));
  }
  for(Int_t iSecPeakPar=0; iSecPeakPar<fNParsSec; iSecPeakPar++) {
    initpars.push_back(fMassFuncFromPrefit->GetParameter(iSecPeakPar+fNParsMassBkg+fNParsMassSgn));
  }
  for(Int_t iReflPar=0; iReflPar<fNParsRfl; iReflPar++) {
    initpars.push_back(fMassFuncFromPrefit->GetParameter(iReflPar+fNParsMassBkg+fNParsMassSgn+fNParsSec));
  }
  for(Int_t iTemplPar=0; iTemplPar<fNParsTempls; iTemplPar++) {
    initpars.push_back(fMassFuncFromPrefit->GetParameter(iTemplPar+fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl));
  }
  for(Int_t iVnBkgPar=0; iVnBkgPar<fNParsVnBkg; iVnBkgPar++) {
    if(vnprefit) {initpars.push_back(fVnBkgFuncSb->GetParameter(iVnBkgPar));}
    else {initpars.push_back(0.05);}
  }
  initpars.push_back(0.10); //initial parameter for signal vn
  if(fSecondPeak && fDoSecondPeakVn) {initpars.push_back(0.10);} //initial parameter for second peak vn

  fMassTotFunc = new TF1("fMassTotFunc",this,&VnVsMassFitter::MassFunc,fMassMin,fMassMax,nparsmass,"VnVsMassFitter","MassFunc");
  fVnTotFunc = new TF1("fVnTotFunc",this,&VnVsMassFitter::vnFunc,fMassMin,fMassMax,nparsvn,"VnVsMassFitter","vnFunc");
  SetParNames();

  ROOT::Math::WrappedMultiTF1 wfTotMass(*fMassTotFunc,1);
  ROOT::Math::WrappedMultiTF1 wfTotVn(*fVnTotFunc,1);

  // set data options and ranges
  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange rangeMass; //same range for two functions

  rangeMass.SetRange(fMassMin,fMassMax);
  ROOT::Fit::BinData dataMass(opt,rangeMass);
  ROOT::Fit::FillData(dataMass, fMassHisto);
  ROOT::Fit::BinData dataVn(opt,rangeMass);
  ROOT::Fit::FillData(dataVn, fVnVsMassHisto);

  //define the 2 chi squares
  ROOT::Fit::Chi2Function chi2Mass(dataMass, wfTotMass);
  ROOT::Fit::Chi2Function chi2Vn(dataVn, wfTotVn);

  //define the global chi square
  GlobalChi2 globalChi2(chi2Mass, chi2Vn);

  //define fitter
  ROOT::Fit::Fitter fitter;
  // create before the parameter settings in order to fix or set range on them
  fitter.Config().SetParamsSettings(nparsvn,initpars.data()); //set initial parameters from prefits
  if(fMeanFixed==2 || fMeanFixedFromMassFit) {fitter.Config().ParSettings(fNParsMassBkg+1).Fix();}
  fitter.Config().ParSettings(fNParsMassBkg+2).SetLimits(0,1);
  if(fSigmaFixed==2 || fSigmaFixedFromMassFit) {fitter.Config().ParSettings(fNParsMassBkg+2).Fix();}
  if(fMassSgnFuncType==k2Gaus) {
    if(fFrac2GausFixed==2 || fFrac2GausFixedFromMassFit) {fitter.Config().ParSettings(fNParsMassBkg+3).Fix();}
    if(fSigma2GausFixed==2 || fSigma2GausFixedFromMassFit) {fitter.Config().ParSettings(fNParsMassBkg+4).Fix();}
  }
  if(fSecondPeak) {
    if(fFixSecMass) {fitter.Config().ParSettings(fNParsMassBkg+fNParsMassSgn+1).Fix();}
    if(fFixSecWidth) {fitter.Config().ParSettings(fNParsMassBkg+fNParsMassSgn+2).Fix();}
  }
  if(fReflections) {
    if(fFixRflOverSig) fitter.Config().ParSettings(fNParsMassBkg+fNParsMassSgn+fNParsSec).Fix();
    if(fVnRflLimited) fitter.Config().ParSettings(nparsmass+fNParsVnBkg+NvnParsSgn-1).SetLimits(fVnRflMin,fVnRflMax);
  }
  // if(fTemplates) {
    //   for(int iTemplPar=0; iTemplPar<this->fMassInitWeights.size(); iTemplPar++) {
      //     fitter.Config().ParSettings(iTemplPar+fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl).SetValue(fMassInitWeights[iTemplPar]);
      //     if(this->fMassWeightsLowerLims[iTemplPar] > this->fMassWeightsUpperLims[iTemplPar]) {
        //       fitter.Config().ParSettings(iTemplPar+fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl).Fix();
        //     } else {
          //       fitter.Config().ParSettings(iTemplPar+fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl).SetLimits(fMassWeightsLowerLims[iTemplPar],fMassWeightsUpperLims[iTemplPar]);        
          //     }
          //     if(!fTemplSameVnOfSignal) {
            //       fitter.Config().ParSettings(iTemplPar+fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsTempls+fNParsVnBkg+fNParsVnSgn+fNParsVnSecPeak+fNParsRfl).SetValue(fVnInitWeights[iTemplPar]);
            //       if(this->fVnWeightsLowerLims[iTemplPar] > this->fVnWeightsUpperLims[iTemplPar]) {
              //         fitter.Config().ParSettings(iTemplPar+fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsTempls+fNParsVnBkg+fNParsVnSgn+fNParsVnSecPeak+fNParsRfl).Fix();
              //       } else {
                //         fitter.Config().ParSettings(iTemplPar+fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsTempls+fNParsVnBkg+fNParsVnSgn+fNParsVnSecPeak+fNParsRfl).SetLimits(fVnWeightsLowerLims[iTemplPar],fVnWeightsUpperLims[iTemplPar]);        
                //       }
                //     }
                //   }
                // }
                
  if(fInitFuncPars.size()>0) {
    for (int iInitPar=0; iInitPar<fInitFuncPars.size(); iInitPar++) {
      int parIdx = fVnTotFunc->GetParNumber(std::get<0>(fInitFuncPars[iInitPar]).Data());
      if (parIdx < fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl+fNParsTempls) {
        cout << "[VnVsMassFitter] Parameter " << std::get<0>(fInitFuncPars[iInitPar]) << " at index " << parIdx << " is a mass function parameter, skipping it." << endl;
        // The parameter init value was set in the mass prefit, now the 
        // value determined in the mass prefit will be used as init value 
        continue;
      }
      fitter.Config().ParSettings(parIdx).SetValue(std::get<1>(fInitFuncPars[iInitPar]));
      if (std::get<2>(fInitFuncPars[iInitPar]) >= std::get<3>(fInitFuncPars[iInitPar])) {
        cout << "[VnVsMassFitter] Fixing parameter " << std::get<0>(fInitFuncPars[iInitPar]) << " at index " << parIdx << " to " << std::get<1>(fInitFuncPars[iInitPar]);
        cout << " with limits " << std::get<2>(fInitFuncPars[iInitPar]) << " and " << std::get<3>(fInitFuncPars[iInitPar]) << endl;
        fitter.Config().ParSettings(parIdx).Fix();
      } else {
        cout << "[VnVsMassFitter] Setting parameter " << std::get<0>(fInitFuncPars[iInitPar]) << " at index " << parIdx << " to " << std::get<1>(fInitFuncPars[iInitPar]);
        cout << " with limits " << std::get<2>(fInitFuncPars[iInitPar]) << " and " << std::get<3>(fInitFuncPars[iInitPar]) << endl;
        fitter.Config().ParSettings(parIdx).SetLimits(std::get<2>(fInitFuncPars[iInitPar]), std::get<3>(fInitFuncPars[iInitPar]));
      }
    }
  }
  
  // When sgn func is Double CrystalBall, the sgnInt, alpha and N parameters are fixed 
  // to the ones obtained from the mass prefit (alpha and N fixed to the MC from the config)
  if (fMassSgnFuncType==kDoubleCBSymm) {
    // fitter.Config().ParSettings(fNParsMassBkg).Fix();
    // fitter.Config().ParSettings(fNParsMassBkg+2).Fix();   // alpha
    // fitter.Config().ParSettings(fNParsMassBkg+3).Fix();   // alpha
    // fitter.Config().ParSettings(fNParsMassBkg+4).Fix();   // N
  }
  if (fMassSgnFuncType==kDoubleCBAsymm) {
    // fitter.Config().ParSettings(fNParsMassBkg).Fix();
    fitter.Config().ParSettings(fNParsMassBkg+2).SetLimits(initpars[fNParsMassBkg+2] - (initpars[fNParsMassBkg+2] / 5), initpars[fNParsMassBkg+2] + (initpars[fNParsMassBkg+2] / 5));   // alpha
    fitter.Config().ParSettings(fNParsMassBkg+3).SetLimits(initpars[fNParsMassBkg+3] - (initpars[fNParsMassBkg+3] / 5), initpars[fNParsMassBkg+3] + (initpars[fNParsMassBkg+3] / 5));   // alpha
    fitter.Config().ParSettings(fNParsMassBkg+4).SetLimits(initpars[fNParsMassBkg+4] - (initpars[fNParsMassBkg+4] / 5), initpars[fNParsMassBkg+4] + (initpars[fNParsMassBkg+4] / 5));   // N
    fitter.Config().ParSettings(fNParsMassBkg+5).SetLimits(initpars[fNParsMassBkg+5] - (initpars[fNParsMassBkg+5] / 5), initpars[fNParsMassBkg+5] + (initpars[fNParsMassBkg+5] / 5));   // N
    fitter.Config().ParSettings(fNParsMassBkg+6).SetLimits(initpars[fNParsMassBkg+6] - (initpars[fNParsMassBkg+6] / 5), initpars[fNParsMassBkg+6] + (initpars[fNParsMassBkg+6] / 5));   // N

    // fitter.Config().ParSettings(fNParsMassBkg+2).Fix();   // alpha
    // fitter.Config().ParSettings(fNParsMassBkg+3).Fix();   // alpha
    // fitter.Config().ParSettings(fNParsMassBkg+4).Fix();   // N
    // fitter.Config().ParSettings(fNParsMassBkg+5).Fix();   // N
    // fitter.Config().ParSettings(fNParsMassBkg+6).Fix();   // N
  }
  
  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");
  for(Int_t iPar=0; iPar<nparsvn; iPar++) {fitter.Config().ParSettings(iPar).SetName(fVnTotFunc->GetParName(iPar));}
  
  // fit FCN function directly
  // (specify optionally data size and flag to indicate that is a chi2 fit
  
  // Set limit > 0 for integrals fit parameters
  fitter.Config().ParSettings(0).SetLimits(0,10000000); // Bkg integral
  fitter.Config().ParSettings(this->fNParsMassBkg).SetLimits(0,10000000); // Sgn integral
  Bool_t isFitOk = fitter.FitFCN(nparsvn,globalChi2,0,dataMass.Size()+dataVn.Size(),kFALSE);
  if(!isFitOk) return kFALSE;

  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);
  if(fTemplates && fTemplSameVnOfSignal) {
    printf("\n --->Templates share the vn parameter with the signal! \n");
  }

  //set parameters in every function
  fVnBkgFunc = new TF1("fVnBkgFunc",this,&VnVsMassFitter::vnBkgFunc,fMassMin,fMassMax,fNParsVnBkg,"VnVsMassFitter","vnBkgFunc");
  fMassBkgFunc = new TF1("fMassBkgFunc",this,&VnVsMassFitter::MassBkg,fMassMin,fMassMax,fNParsMassBkg,"VnVsMassFitter","MassBkg");
  fMassSgnFunc = new TF1("fMassSgnFunc",this,&VnVsMassFitter::MassSignal,fMassMin,fMassMax,fNParsMassSgn,"VnVsMassFitter","MassSignal");
  if(fReflections) {fMassRflFunc = new TF1("fMassRflFunc",this,&VnVsMassFitter::MassRfl,fMassMin,fMassMax,fNParsRfl,"VnVsMassFitter","MassRfl");}
  if(fReflections) {fMassBkgRflFunc = new TF1("fMassBkgRflFunc",this,&VnVsMassFitter::MassBkgRfl,fMassMin,fMassMax,fNParsMassBkg+fNParsRfl,"VnVsMassFitter","MassBkgRfl");}
  if(fSecondPeak) {fMassSecPeakFunc = new TF1("fMassSecPeakFunc",this,&VnVsMassFitter::MassSecondPeak,fMassMin,fMassMax,fNParsSec,"VnVsMassFitter","MassSecondPeak");}
  for(Int_t iPar=0; iPar<nparsvn; iPar++) {
    fVnTotFunc->SetParameter(iPar,result.Parameter(iPar));
    fVnTotFunc->SetParError(iPar,result.ParError(iPar));
    if(iPar<nparsmass) {
      fMassTotFunc->SetParameter(iPar,result.Parameter(iPar));
      fMassTotFunc->SetParError(iPar,result.ParError(iPar));
    }
    if(iPar>=nparsmass && iPar<nparsvn-NvnParsSgn) {
      fVnBkgFunc->SetParameter(iPar-nparsmass,result.Parameter(iPar));
      fVnBkgFunc->SetParError(iPar-nparsmass,result.ParError(iPar));
    }
    if(iPar>=fNParsMassBkg && iPar<fNParsMassBkg+fNParsMassSgn) {
      fMassSgnFunc->SetParameter(iPar-fNParsMassBkg,result.Parameter(iPar));
    }
    if(iPar<fNParsMassBkg) {
      fMassBkgFunc->SetParameter(iPar,result.Parameter(iPar));
      fMassBkgFunc->SetParError(iPar,result.ParError(iPar));
      if(fReflections) {
        fMassBkgRflFunc->SetParameter(iPar,result.Parameter(iPar));
        fMassBkgRflFunc->SetParError(iPar,result.ParError(iPar));
      }
    }
    if(fReflections && (iPar>=fNParsMassBkg+fNParsMassSgn+fNParsSec && iPar<fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl)) {
      fMassRflFunc->SetParameter(iPar-(fNParsMassBkg+fNParsMassSgn+fNParsSec),result.Parameter(iPar));
      fMassRflFunc->SetParError(iPar-(fNParsMassBkg+fNParsMassSgn+fNParsSec),result.ParError(iPar));
      fMassBkgRflFunc->SetParameter(iPar-(fNParsMassSgn+fNParsSec),result.Parameter(iPar));
      fMassBkgRflFunc->SetParError(iPar-(fNParsMassSgn+fNParsSec),result.ParError(iPar));
    }
    if(fSecondPeak && (iPar>=fNParsMassBkg+fNParsMassSgn && iPar<fNParsMassBkg+fNParsMassSgn+fNParsSec)) {
      fMassSecPeakFunc->SetParameter(iPar-(fNParsMassBkg+fNParsMassSgn),result.Parameter(iPar));
      fMassSecPeakFunc->SetParError(iPar-(fNParsMassBkg+fNParsMassSgn),result.ParError(iPar));
    }
  }
  
  if(fTemplates) {
    switch (fAnchorTemplsMode) {
      case TemplAnchorMode::Free:
      fMassTemplFunc = new TF1("fMassTemplFunc",this,&VnVsMassFitter::MassTemplates,fMassMin,fMassMax,fNParsTempls,"VnVsMassFitter","MassTemplates");
      for (int iTempl=0; iTempl<fNParsTempls; iTempl++) {
        fMassTemplFunc->SetParameter(iTempl, result.Parameter(iTempl + (nparsmass - fNParsTempls)));
      }
      break;
      case TemplAnchorMode::AnchorToFirst:
        fMassTemplFunc = new TF1("fMassTemplFunc",this,&VnVsMassFitter::MassTemplates,fMassMin,fMassMax,1,"VnVsMassFitter","MassTemplates");
        fMassTemplFunc->SetParameter(0, result.Parameter(nparsmass - fNParsTempls));
        break;
        case TemplAnchorMode::AnchorToSgn:
        fMassTemplFunc = new TF1("fMassTemplFunc",this,&VnVsMassFitter::MassTemplates,fMassMin,fMassMax,1,"VnVsMassFitter","MassTemplates");
        fMassTemplFunc->SetParameter(0, result.Parameter(fNParsMassBkg));
        break;
      default:
        std::cerr << "Error: Invalid fAnchorTemplsMode value!" << std::endl;
    }
    cout << "fNParsTempls: " << fNParsTempls << endl; 
    cout << "fMassTemplFunc->Eval(1.861): " << fMassTemplFunc->Eval(1.861) << endl;
  }

  fVn = fVnTotFunc->GetParameter(fVnTotFunc->GetNpar()-NvnParsSgn);
  fVnUncertainty = fVnTotFunc->GetParError(fVnTotFunc->GetNpar()-NvnParsSgn);
  fRawYield = fVnTotFunc->GetParameter(fNParsMassBkg)/fMassHisto->GetBinWidth(10);
  fRawYieldUncertainty = fVnTotFunc->GetParError(fNParsMassBkg)/fMassHisto->GetBinWidth(10);
  fMean = fVnTotFunc->GetParameter(fNParsMassBkg+1);
  fMeanUncertainty = fVnTotFunc->GetParError(fNParsMassBkg+1);
  fSigma = fVnTotFunc->GetParameter(fNParsMassBkg+2);
  fSigmaUncertainty = fVnTotFunc->GetParError(fNParsMassBkg+2);
  fChiSquare = result.MinFcnValue();
  fNDF = result.Ndf();
  fProb = result.Prob();

  // Get Vn components to be drawn
  int idxParMassTemplsScaling = fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl;
  int idxParVnSgn = idxParMassTemplsScaling+fNParsTempls+fNParsVnBkg;
  int idxParVnTempl = idxParVnSgn+fNParsVnSgn+fNParsVnSecPeak+fNParsVnRfl;
  double vnSgn = result.Parameter(idxParVnSgn);
  fVnCompsDraw.push_back(new TF1("vnSgn",
    [&, this, vnSgn] (double *x, double *par) {
      return (vnSgn * this->fMassSgnFunc->Eval(x[0])) / (this->fMassTotFunc->Eval(x[0]));
    }, fMassMin, fMassMax, 0));
    fVnCompsDraw.push_back(new TF1("vnBkg",
      [&, this] (double *x, double *par) {
        return (this->fVnBkgFunc->Eval(x[0]) * this->fMassBkgFunc->Eval(x[0])) / (this->fMassTotFunc->Eval(x[0]));
      }, fMassMin, fMassMax, 0));
      if(fDoSecondPeakVn) {
        fVnSecPeak = fVnTotFunc->GetParameter(fVnTotFunc->GetNpar()-1);
        fVnSecPeakUncertainty = fVnTotFunc->GetParError(fVnTotFunc->GetNpar()-1);
        fVnSecPeakFunc = new TF1("vnSecPeak",
          [&, this] (double *x, double *par) {
            return (this->fVnSecPeak * this->fMassSecPeakFunc->Eval(x[0])) / (this->fMassTotFunc->Eval(x[0]));
          }, fMassMin, fMassMax, 0);
          fVnCompsDraw.push_back(fVnSecPeakFunc); 
        }
        if(fTemplates) {

          double templScalingPar = 0.;
          for(int iTempl=0; iTempl<fHistoTemplates.size(); iTempl++) {
            switch (fAnchorTemplsMode) {
              case TemplAnchorMode::Free:
                templScalingPar = result.Parameter(iTempl + idxParMassTemplsScaling);
                break;
              case TemplAnchorMode::AnchorToFirst:
                templScalingPar = result.Parameter(idxParMassTemplsScaling) * this->fRelWeights[iTempl];
                break;
              case TemplAnchorMode::AnchorToSgn:
                templScalingPar = result.Parameter(this->fNParsMassBkg) * this->fRelWeights[iTempl];
                break;
              default:
                std::cerr << "Error: Invalid fAnchorTemplsMode value!" << std::endl;
            }
            cout << "Mass scaling parameter of " << iTempl << ": " << templScalingPar << endl;
            cout << "this->fRelWeights.size(): " << this->fRelWeights.size() << endl;
            cout << "Rel weight of " << iTempl << ": " << this->fRelWeights[iTempl] << endl;
            fMassTemplatesDraw.push_back(new TF1(Form("fTempl_%i", iTempl),
            [&, this, iTempl, templScalingPar] (double *x, double *par) {
              double xval = x[0];
              this->fMassVar.setVal(xval);
              double pdfVal = fHistoTemplates[iTempl]->getVal(RooArgSet(this->fMassVar));
              // std::cout << "PDF value at mass = " << xval << " is: " << pdfVal << std::endl;
              return templScalingPar * pdfVal;
            }, fMassMin, fMassMax, 0));
            
            if(fTemplSameVnOfSignal) {
              fVnCompsDraw.push_back(new TF1(Form("fVnTempl_%i", iTempl),
              [&, this, iTempl, vnSgn] (double *x, double *par) {
                return (vnSgn * this->fRelWeights[iTempl] * this->fMassTemplatesDraw[iTempl]->Eval(x[0])) / (this->fMassTotFunc->Eval(x[0]));
              }, fMassMin, fMassMax, 0));
            } else {
        fVnCompsDraw.push_back(new TF1(Form("fVnTempl_%i", iTempl),
        [&, this, iTempl, idxParVnTempl, result] (double *x, double *par) {
          double templVnScalingPar = result.Parameter(iTempl + idxParVnTempl);
          return (templVnScalingPar * this->fRelWeights[iTempl] * this->fMassTemplatesDraw[iTempl]->Eval(x[0])) / (this->fMassTotFunc->Eval(x[0]));
        }, fMassMin, fMassMax, 0));
      }
    }
  }
  return kTRUE;
}

//________________________________________________________________
Bool_t VnVsMassFitter::MassPrefit() {

  // //define proper maxs and mins from histos
  Double_t tmpmin = TMath::Max(fMassHisto->GetBinLowEdge(1),fVnVsMassHisto->GetBinLowEdge(1));
  fMassMin=TMath::Max(fMassMin,tmpmin);
  Double_t tmpmax = TMath::Min(fMassHisto->GetBinLowEdge(fMassHisto->GetNbinsX())+fMassHisto->GetBinWidth(fMassHisto->GetNbinsX()),fVnVsMassHisto->GetBinLowEdge(fVnVsMassHisto->GetNbinsX())+fVnVsMassHisto->GetBinWidth(fVnVsMassHisto->GetNbinsX()));
  fMassMax=TMath::Min(fMassMax,tmpmax);
  
  // fMassFitter = new InvMassFitter(fMassHisto,fMassMin,fMassMax,fMassBkgFuncType,fMassSgnFuncType);
  fMassFitter->SetNSigma4SideBands(fNSigmaForSB);
  if(fSigmaFixed==1) fMassFitter->SetInitialGaussianSigma(fSigmaInit);
  else if(fSigmaFixed==2) fMassFitter->SetFixGaussianSigma(fSigmaInit);
  if(fMeanFixed==1) fMassFitter->SetInitialGaussianMean(fMeanInit);
  else if(fMeanFixed==2) fMassFitter->SetFixGaussianMean(fMeanInit);
  if(fMassSgnFuncType==k2Gaus) {
    if(fSigma2GausFixed==1) fMassFitter->SetInitialSecondGaussianSigma(fSigma2GausInit);
    else if(fSigma2GausFixed==2) fMassFitter->SetFixSecondGaussianSigma(fSigma2GausInit);
    if(fFrac2GausFixed==1) fMassFitter->SetInitialFrac2Gaus(fFrac2GausInit);
    else if(fFrac2GausFixed==2) fMassFitter->SetFixFrac2Gaus(fFrac2GausInit);
  }
  fMassFitter->SetUseLikelihoodFit();
  if(fMassBkgFuncType==kPoln) {fMassFitter->SetPolDegreeForBackgroundFit(fPolDegreeBkg);}
  if(fSecondPeak) {fMassFitter->IncludeSecondGausPeak(fSecMass,fFixSecMass,fSecWidth,fFixSecWidth);}
  if(fReflections) {
    fHistoTemplRfl = (TH1F*)fMassFitter->SetTemplateReflections(fHistoTemplRflInit,fRflOpt,fMinRefl,fMaxRefl)->Clone("fHistoTemplRfl");
    if(fRflOverSig>0) {fMassFitter->SetInitialReflOverS(fRflOverSig);}
    if(fFixRflOverSig) {fMassFitter->SetFixReflOverS(fRflOverSig);}
  }
  if (fInitFuncPars.size() > 0) {
    fMassFitter->SetInitPars(fInitFuncPars);
  }
  Bool_t status = fMassFitter->MassFitter(kFALSE);

  if(status) {
    fMassFuncFromPrefit = (TF1*)fMassFitter->GetMassFunc();
    fMassFuncFromPrefit->SetName("fMassFuncFromPrefit");
    fMassPrefitChiSquare = fMassFitter->GetChiSquare();
    fMassPrefitNDF       = fMassFitter->GetMassFunc()->GetNDF();
    fMassPrefitProb      = fMassFitter->GetFitProbability();
  }
  if(fReflections) fRawYieldHelp=fMassFitter->GetRawYield();

  return status;
}

//________________________________________________________________
Bool_t VnVsMassFitter::VnSBPrefit() {

  Double_t mean = fMassFitter->GetMean();
  Double_t sigma = fMassFitter->GetSigma();
  const Int_t nMassBins = fVnVsMassHisto->GetNbinsX();
  Double_t SBbins[nMassBins];
  Int_t nSBbins=0;
  for(Int_t iBin=0; iBin<nMassBins; iBin++) {
    Double_t min = fVnVsMassHisto->GetBinLowEdge(iBin+1);
    Double_t max = fVnVsMassHisto->GetBinLowEdge(iBin+1)+fVnVsMassHisto->GetBinWidth(iBin+1);
    if(max<mean){
      if(max<(mean-fNSigmaForSB*sigma)) {SBbins[iBin]=1; nSBbins++;}
      else {SBbins[iBin]=0;}
    }
    if(min>=mean){
      if(min>(mean+fNSigmaForSB*sigma)) {SBbins[iBin]=1; nSBbins++;}
      else {SBbins[iBin]=0;}
    }
  }
  TGraphErrors* gVnVsMassSB = new TGraphErrors(nSBbins);
  for(Int_t iBin=0; iBin<nMassBins; iBin++) {
    if(SBbins[iBin]==1) {
      gVnVsMassSB->SetPoint(iBin,fVnVsMassHisto->GetBinCenter(iBin+1),fVnVsMassHisto->GetBinContent(iBin+1));
      gVnVsMassSB->SetPointError(iBin,fVnVsMassHisto->GetBinWidth(iBin+1)/2,fVnVsMassHisto->GetBinError(iBin+1));
    }
  }
  fVnBkgFuncSb = new TF1("fVnBkgFuncSb",this,&VnVsMassFitter::vnBkgFunc,fMassMin,fMassMax,fNParsVnBkg,"VnVsMassFitter","vnBkgFunc");
  switch(fVnBkgFuncType) {
    case 1:
      fVnBkgFuncSb->SetParName(0,"ConstVnBkg");
      fVnBkgFuncSb->SetParName(1,"SlopeVnBkg");
      break;
    case 2:
      fVnBkgFuncSb->SetParName(0,"ConstVnBkg");
      fVnBkgFuncSb->SetParName(1,"Coef1VnBkg");
      fVnBkgFuncSb->SetParName(2,"Coef2VnBkg");
      break;
    default:
      printf("Error in setting signal par names: check fVnBkgFuncType");
      break;
  }
  Int_t status = gVnVsMassSB->Fit(fVnBkgFuncSb,"","",fMassMin,fMassMax);

  fSBVnPrefitChiSquare = fVnBkgFuncSb->GetChisquare();
  fSBVnPrefitNDF       = fVnBkgFuncSb->GetNDF();
  fSBVnPrefitProb      = fVnBkgFuncSb->GetProb();

  delete gVnVsMassSB;

  if(status==0) return kTRUE;
  return kFALSE;
}

//________________________________________________________________
void VnVsMassFitter::DefineNumberOfParameters() {

  switch(fMassSgnFuncType) {
    case 0: //single gaus
      fNParsMassSgn=3;
      break;
    case 1: //double gaus
      fNParsMassSgn=5;
      break;
    case 3:
      cout << "Setting number of parameters for signal to 7" << endl;
      fNParsMassSgn=7;
      break;
    case 4:
      cout << "Setting number of parameters for signal to 5" << endl;
      fNParsMassSgn=5;
      break;
    default:
      printf("Error in computing fMassSgnFuncType: check fMassSgnFuncType");
      break;
  }

  switch(fMassBkgFuncType) {
    case 0: //expo
      fNParsMassBkg=2;
      break;
    case 1: //lin
      fNParsMassBkg=2;
      break;
    case 2: //pol2
      fNParsMassBkg=3;
      break;
    case 3: //no bkg
      fNParsMassBkg=1;
      break;
    case 4: //power law
      fNParsMassBkg=2;
      break;
    case 5: //power expo
      fNParsMassBkg=2;
      break;
    case 6: //high degree pol
      fNParsMassBkg=fPolDegreeBkg+1;
      break;
    default:
      printf("Error in computing fNParsMassBkg: check fMassBkgFuncType");
      break;
  }

  switch(fVnBkgFuncType) {
    case 0: //expo
      fNParsVnBkg=2;
      break;
    case 1: //lin
      fNParsVnBkg=2;
      break;
    case 2: //pol2
      fNParsVnBkg=3;
      break;
    case 6: //high degree pol
      fNParsVnBkg=fPolDegreeVnBkg+1;
      break;
    default:
      printf("Error in computing fNParsVnBkg: check fVnBkgFuncType");
      break;
  }

  fNParsVnSgn=1;
  
  if(fReflections) {
    fNParsRfl=1;
    if(fVnRflOpt==3) fNParsVnRfl=1;
    else fNParsVnRfl=0;
  }
  else {
    fNParsRfl=0;
    fNParsVnRfl=0;
  }

  if (fTemplates) {
    switch (fAnchorTemplsMode) {
      case TemplAnchorMode::Free:
        fNParsTempls = fHistoTemplates.size();
        break;
      case TemplAnchorMode::AnchorToFirst:
        fNParsTempls = 1;
        break;
       case TemplAnchorMode::AnchorToSgn:
        fNParsTempls = 0;
        break;
      default:
        std::cerr << "Error: Invalid fAnchorTemplsMode value!" << std::endl;
    }
  }

  else {
  }
  
  if(fSecondPeak) {
    fNParsSec=3;
    fNParsVnSecPeak=1;
  }
  else {
    fNParsSec=0;
    fNParsVnSecPeak=0;
  }
}

//________________________________________________________________
void VnVsMassFitter::SetParNames() {

  switch(fMassSgnFuncType) {
    case 0: //single gaus
      fVnTotFunc->SetParName(fNParsMassBkg,"SgnInt");
      fVnTotFunc->SetParName(fNParsMassBkg+1,"Mean");
      fVnTotFunc->SetParName(fNParsMassBkg+2,"Sigma");
      break;
    case 1: //double gaus
      fVnTotFunc->SetParName(fNParsMassBkg,"SgnInt");
      fVnTotFunc->SetParName(fNParsMassBkg+1,"Mean");
      fVnTotFunc->SetParName(fNParsMassBkg+2,"Sigma1");
      fVnTotFunc->SetParName(fNParsMassBkg+3,"Frac");
      fVnTotFunc->SetParName(fNParsMassBkg+4,"Sigma2");
      break;
    case 3: //asymmetric crystalball
      fVnTotFunc->SetParName(fNParsMassBkg,"SgnInt");
      fVnTotFunc->SetParName(fNParsMassBkg+1,"Mean");
      fVnTotFunc->SetParName(fNParsMassBkg+2,"Sigma");
      fVnTotFunc->SetParName(fNParsMassBkg+3,"Alpha1");
      fVnTotFunc->SetParName(fNParsMassBkg+4,"N1");
      fVnTotFunc->SetParName(fNParsMassBkg+5,"Alpha2");
      fVnTotFunc->SetParName(fNParsMassBkg+6,"N2");
      break;
    case 4: //symmetric crystalball
      fVnTotFunc->SetParName(fNParsMassBkg,"SgnInt");
      fVnTotFunc->SetParName(fNParsMassBkg+1,"Mean");
      fVnTotFunc->SetParName(fNParsMassBkg+2,"Sigma");
      fVnTotFunc->SetParName(fNParsMassBkg+3,"Alpha");
      fVnTotFunc->SetParName(fNParsMassBkg+4,"N");
      break;

    default:
      printf("Error in setting signal par names: check fMassSgnFuncType");
      break;
  }
  
  switch(fMassBkgFuncType) {
    case 0: //expo
      fVnTotFunc->SetParName(0,"BkgInt");
      fVnTotFunc->SetParName(1,"Slope");
      break;
    case 1: //lin
      fVnTotFunc->SetParName(0,"BkgInt");
      fVnTotFunc->SetParName(1,"Slope");
      break;
    case 2: //pol2
      fVnTotFunc->SetParName(0,"BkgInt");
      fVnTotFunc->SetParName(1,"Coef1");
      fVnTotFunc->SetParName(2,"Coef2");
      break;
    case 3: //no bkg
      fVnTotFunc->SetParName(0,"Const");
      break;
    case 4: //power law
      fVnTotFunc->SetParName(0,"BkgInt");
      fVnTotFunc->SetParName(1,"Coef1");
      break;
    case 5: //power expo
      fVnTotFunc->SetParName(0,"Coef1");
      fVnTotFunc->SetParName(1,"Coef2");
      break;
    case 6: //high degree pol
      fVnTotFunc->SetParName(0,"BkgInt");
      for(Int_t iPar=1; iPar<fNParsMassBkg; iPar++) {fVnTotFunc->SetParName(iPar,Form("Coef%d",iPar));}
      break;
    default:
      printf("Error in setting signal par names: check fMassBkgFuncType");
      break;
  }

  for(Int_t iPar=0; iPar<fNParsVnBkg; iPar++) {fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl+fNParsTempls+iPar,fVnBkgFuncSb->GetParName(iPar));}

  if(fReflections) {fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+fNParsSec,"ReflOverS");}
  if(fTemplates) {
    switch (fAnchorTemplsMode) {
      case TemplAnchorMode::Free:
      for(int iTempl=0; iTempl<this->fNParsTempls; iTempl++) {
        fVnTotFunc->SetParName(iTempl+fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl,Form("wm_%s", fHistoTemplates[iTempl]->GetName()));
      }
      break;
      case TemplAnchorMode::AnchorToFirst:
      fVnTotFunc->SetParName(fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl,"wtempls_anchor_first");
      break;
      case TemplAnchorMode::AnchorToSgn:
      cout << "Templates anchored to signal, no parameter for them!" << endl;
      break;
      default:
      std::cerr << "Error: Invalid fAnchorTemplsMode value!" << std::endl;
    }
    if (!fTemplSameVnOfSignal){
      for(int iTempl=0; iTempl<this->fNParsTempls; iTempl++) {
        fVnTotFunc->SetParName(fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsTempls+fNParsVnBkg+fNParsVnSgn+fNParsVnSecPeak+fNParsRfl+iTempl,Form("wvn_%s", fHistoTemplates[iTempl]->GetName()));
      }
    }
  }

  if(fSecondPeak) {
    fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn,"SecPeakInt");
    fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+1,"SecPeakMean");
    fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+2,"SecPeakSigma");
  }
  fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+fNParsRfl+fNParsSec+fNParsTempls+fNParsVnBkg,Form("v%dSgn",fHarmonic));
  
  if(fSecondPeak && fDoSecondPeakVn) {fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+fNParsRfl+fNParsSec+fNParsTempls+fNParsVnBkg+1,Form("v%dSecPeak",fHarmonic));}
  if(fReflections && fVnRflOpt==kFreePar) {fVnTotFunc->SetParName(fNParsMassBkg+fNParsMassSgn+fNParsRfl+fNParsSec+fNParsVnBkg+fNParsTempls+1+fNParsVnSecPeak,Form("v%dRefl",fHarmonic));}
}

//_________________________________________________________________________
void VnVsMassFitter::Signal(Double_t nOfSigma,Double_t &signal,Double_t &errsignal) const {
  /// Return signal integral in mean +- n sigma
  ///

  Double_t minMass=fMean-nOfSigma*fSigma;
  Double_t maxMass=fMean+nOfSigma*fSigma;
  Signal(minMass,maxMass,signal,errsignal);
  return;
}

//_________________________________________________________________________
void VnVsMassFitter::Signal(Double_t min, Double_t max, Double_t &signal,Double_t &errsignal) const {
  /// Return signal integral in a range
  ///
  if(!fMassSgnFunc) {signal=-1; errsignal=0; return;}

  signal=fMassSgnFunc->Integral(min, max)/(Double_t)fMassHisto->GetBinWidth(1);
  errsignal=(fRawYieldUncertainty/fRawYield)*signal;/*assume relative error is the same as for total integral*/
  
  return;
}

//___________________________________________________________________________
void VnVsMassFitter::Background(Double_t nOfSigma,Double_t &background,Double_t &errbackground) const {
  /// Return background integral in mean +- n sigma
  ///

  Double_t minMass=fMean-nOfSigma*fSigma;
  Double_t maxMass=fMean+nOfSigma*fSigma;
  Background(minMass,maxMass,background,errbackground);

  return;
}

//___________________________________________________________________________
void VnVsMassFitter::Background(Double_t min, Double_t max, Double_t &background,Double_t &errbackground) const {
  /// Return background integral in a range
  ///

  if(!fMassBkgFunc) {background=-1; errbackground=0; return;}

  Double_t intB=fMassBkgFunc->GetParameter(0);
  Double_t intBerr=fMassBkgFunc->GetParError(0);
  //relative error evaluation: from histo

  Int_t leftBand=fMassHisto->FindBin(fMean-4*fSigma);
  Int_t rightBand=fMassHisto->FindBin(fMean+4*fSigma);
  intB=fMassHisto->Integral(1,leftBand)+fMassHisto->Integral(rightBand,fMassHisto->GetNbinsX());
  Double_t sum2=0;
  for(Int_t iBin=1; iBin<=leftBand; iBin++){
    sum2+=fMassHisto->GetBinError(iBin)*fMassHisto->GetBinError(iBin);
  }
  for(Int_t iBin=rightBand; iBin<=fMassHisto->GetNbinsX(); iBin++){
    sum2+=fMassHisto->GetBinError(iBin)*fMassHisto->GetBinError(iBin);
  }

  intBerr=TMath::Sqrt(sum2);

  background=fMassBkgFunc->Integral(min,max)/(Double_t)fMassHisto->GetBinWidth(1);
  errbackground=intBerr/intB*background;

  return;
}

//__________________________________________________________________________
void VnVsMassFitter::Significance(Double_t nOfSigma,Double_t &significance,Double_t &errsignificance) const  {
  /// Return significance in mean +- n sigma
  ///

  Double_t minMass=fMean-nOfSigma*fSigma;
  Double_t maxMass=fMean+nOfSigma*fSigma;
  Significance(minMass, maxMass, significance, errsignificance);

  return;
}

//__________________________________________________________________________
void VnVsMassFitter::Significance(Double_t min, Double_t max, Double_t &significance,Double_t &errsignificance) const {
  /// Return significance integral in a range
  ///

  Double_t background,errbackground;
  Background(min,max,background,errbackground);

  if (fRawYield+background <= 0.){
    significance=-1;
    errsignificance=0;
    return;
  }

  InvMassFitter::ComputeSignificance(fRawYield,fRawYieldUncertainty,background,errbackground,significance,errsignificance);

  return;
}

//__________________________________________________________________________
Double_t VnVsMassFitter::DoubleSidedCBAsymmForVn(double x, double mu, double width, double a1, double n1, double a2, double n2) {
    // Define the variable for the PDF (e.g., mass)
    // RooRealVar mass("mass", "Mass", this->fMassMin, this->fMassMax);
    this->fMassVar.setVal(x);  // Set the mass value to x

    // Define the parameters for the Crystal Ball PDF
    RooRealVar mean("mean", "Mean", mu);  // Set mean to mu
    RooRealVar sigma("sigma", "Sigma", width, 0.001, 5.0);  // Set sigma to width
    RooRealVar alphaLeft("alphaLeft", "AlphaLeft", a1, 0.0001, 10);  // Left tail parameter
    RooRealVar nLeft("nLeft", "nLeft", n1, 0.0001, 10);  // Left tail exponent
    RooRealVar alphaRight("alphaRight", "AlphaRight", a2, 0.0001, 10);  // Right tail parameter
    RooRealVar nRight("nRight", "nRight", n2, 0.0001, 10);  // Right tail exponent

    // Create the RooCrystalBall PDF
    RooCrystalBall cb("cb", "Double-Sided Crystal Ball PDF", fMassVar, mean, sigma, alphaLeft, nLeft, alphaRight, nRight);

    // Evaluate the PDF at the specified mass value and return the result
    double pdf_value = cb.getVal(RooArgSet(fMassVar));
    
    // // Print the evaluated PDF value for debugging purposes
    // std::cout << "[DoubleSidedCBAsymmForVn] PDF value at mass = " << x << " with mean = " << mu << " and width = " << width 
    //           << " is: " << pdf_value << std::endl;

    return pdf_value;
}

//__________________________________________________________________________
Double_t VnVsMassFitter::DoubleSidedCBSymmForVn(double x, double mu, double width, double a, double n) {
    // Define the variable for the PDF (e.g., mass)
    // RooRealVar mass("mass", "Mass", this->fMassMin, this->fMassMax);
    // mass.setVal(x);  // Set the mass value to x
    this->fMassVar.setVal(x);  // Set the mass value to x

    // Define the parameters for the Crystal Ball PDF
    RooRealVar mean("mean", "Mean", mu);  // Set mean to mu
    RooRealVar sigma("sigma", "Sigma", width, 0.001, 5.0);  // Set sigma to width
    RooRealVar alphaLeftRight("alphaLeftRight", "AlphaLeftRight", a, 0.0001, 10);  // Left tail parameter
    RooRealVar nLeftRight("nLeftRight", "nLeftRight", n, 0.0001, 10);  // Left tail exponent

    // Create the RooCrystalBall PDF
    RooCrystalBall cb("cb", "Double-Sided Crystal Ball PDF", fMassVar, mean, sigma, alphaLeftRight, nLeftRight, alphaLeftRight, nLeftRight);

    // Evaluate the PDF at the specified mass value and return the result
    double pdf_value = cb.getVal(RooArgSet(fMassVar));
    
    // Print all parameter values
    // std::cout << "[DoubleSidedCBSymmForVn] PDF value at mass = " << x
    //           << " with mean = " << mu 
    //           << ", width (sigma) = " << width
    //           << ", alphaLeftRight = " << a
    //           << ", nLeftRight = " << n
    //           << " is: " << pdf_value << std::endl;

    return pdf_value;
}

//________________________________________________________________
Double_t VnVsMassFitter::GetGausPDF(Double_t x, Double_t mean, Double_t sigma) {

  return TMath::Gaus(x,mean,sigma,kTRUE);
}

//________________________________________________________________
Double_t VnVsMassFitter::GetExpoPDF(Double_t x, Double_t coeff, Bool_t isnorm) {

  if(isnorm) {return TMath::Exp(x/coeff)/(coeff*(TMath::Exp(fMassMax/coeff)-TMath::Exp(fMassMin/coeff)));}
  else return TMath::Exp(x/coeff);
}

//________________________________________________________________
Double_t VnVsMassFitter::GetPolPDF(Double_t x, Double_t *pars, Int_t order, Bool_t isnorm) {

  switch(order) {
    case 0:
      if(isnorm) {return 1./(fMassMax-fMassMin);}
      else {return pars[0];}
      break;
    case 1:
      if(isnorm) {return (pars[0]-pars[1]/2*(fMassMax*fMassMax-fMassMin*fMassMin))/(fMassMax-fMassMin)+pars[1]*x;}
      else {return pars[0]+pars[1]*x;}
      break;
    case 2:
      if(isnorm) {return (pars[0]-pars[1]/2*(fMassMax*fMassMax-fMassMin*fMassMin)-pars[2]/3*(fMassMax*fMassMax*fMassMax-fMassMin*fMassMin*fMassMin))/(fMassMax-fMassMin)+pars[1]*x;}
      else {return pars[0]+pars[1]*x+pars[2]*x*x;}
  }
  return 0;
}

//________________________________________________________________
Double_t VnVsMassFitter::GetPowerFuncPDF(Double_t x, Double_t *pars) {

  Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  return pars[0]*(pars[1]+1.)/(TMath::Power(fMassMax-mpi,pars[1]+1.)-TMath::Power(fMassMin-mpi,pars[1]+1.))*TMath::Power(x-mpi,pars[1]);
}

//________________________________________________________________
Double_t VnVsMassFitter::GetPowerExpoPDF(Double_t x, Double_t *pars) {

  Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  return pars[0]*TMath::Sqrt(x - mpi)*TMath::Exp(-1.*pars[1]*(x-mpi));
}

//________________________________________________________________
Double_t VnVsMassFitter::GetHigherPolFuncPDF(Double_t x, Double_t *pars, Int_t Ndeg, Bool_t isnorm) {

  Double_t total=pars[0];
  for(Int_t iT=1; iT<=Ndeg; iT++){
    if(isnorm) total+=pars[iT]*TMath::Power(x-fMassParticle,iT)/TMath::Factorial(iT);
    else total+=pars[iT]*TMath::Power(x,iT);
  }
  return total;
}

//________________________________________________________________
Double_t VnVsMassFitter::MassSignal(Double_t *m, Double_t *pars) {

  switch(fMassSgnFuncType) {
    case 0:
      return pars[0]*GetGausPDF(m[0],pars[1],pars[2]);
      break;
    case 1:
      return pars[0]*(pars[3]*GetGausPDF(m[0],pars[1],pars[2])+(1-pars[3])*GetGausPDF(m[0],pars[1],pars[4]));
      break;
    case 3:
      return pars[0]*DoubleSidedCBAsymmForVn(m[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6]);
      break;
    case 4:
      return pars[0]*DoubleSidedCBSymmForVn(m[0],pars[1],pars[2],pars[3],pars[4]);
      break;

  }

  return 0;
}

//________________________________________________________________
Double_t VnVsMassFitter::MassBkg(Double_t *m, Double_t *pars) {

  switch(fMassBkgFuncType) {
    case 0: //exponential
      return pars[0]*GetExpoPDF(m[0],pars[1],kTRUE);
      break;
    case 1: //linear
      return GetPolPDF(m[0],pars,1,kTRUE);
      break;
    case 2: //parabolic
      return GetPolPDF(m[0],pars,2,kTRUE);
      break;
    case 3: //constant
      return GetPolPDF(m[0],pars,0,kTRUE);
      break;
    case 4: //power law
      return GetPowerFuncPDF(m[0],pars);
      break;
    case 5: //power law expo
      return GetPowerExpoPDF(m[0],pars);
      break;
    case 6: //higher order (>=3) polinomial
      return GetHigherPolFuncPDF(m[0],pars,fPolDegreeBkg,kTRUE);
      break;
  }
  return 0;
}

//_________________________________________________________________________
Double_t VnVsMassFitter::MassRfl(Double_t *m,Double_t *pars){
  /// Fit function for reflections:
  /// D0->Kpi decays with swapped mass assignment to pion and kaon decay tracks
  if(!fHistoTemplRfl) return 0;

  Int_t bin =fHistoTemplRfl->FindBin(m[0]);
  Double_t value=fHistoTemplRfl->GetBinContent(bin);
  Int_t binmin=fHistoTemplRfl->FindBin(fMassMin*1.00001);
  Int_t binmax=fHistoTemplRfl->FindBin(fMassMax*0.99999);
  Double_t norm=fHistoTemplRfl->Integral(binmin,binmax)*fHistoTemplRfl->GetBinWidth(bin);
  if(TMath::Abs(value)<1.e-14 && fSmoothRfl){// very rough, assume a constant trend, much better would be a pol1 or pol2 over a broader range
    value+=fHistoTemplRfl->GetBinContent(bin-1)+fHistoTemplRfl->GetBinContent(bin+1);
    value/=3.;
  }

  return pars[0]*value/norm*fRawYieldHelp*fMassHisto->GetBinWidth(1);
}

//_________________________________________________________________________
Double_t VnVsMassFitter::MassTemplates(Double_t *m,Double_t *pars){
  // Add the contributions of the templates loaded in fHistoTemplates, each
  // scaled by a multiplicative constant, left as free fit parameter
  this->fMassVar.setVal(m[0]);  // Set the mass value to x
  Double_t totalTemplates = 0.;
  switch (fAnchorTemplsMode) {
    case TemplAnchorMode::Free:
      for(int iTempl=0; iTempl<fHistoTemplates.size(); iTempl++) {
        totalTemplates += pars[iTempl]*fHistoTemplates[iTempl]->getVal(RooArgSet(this->fMassVar));
      }
      break;
    case TemplAnchorMode::AnchorToFirst:
      for(int iTempl=0; iTempl<fHistoTemplates.size(); iTempl++) {
        // cout << "[MassTemplates, first] fRelWeights[" << iTempl << "]: " << fRelWeights[iTempl] << ", pars[0]: " << pars[0] << ", eval templ: " << fHistoTemplates[iTempl]->getVal(RooArgSet(this->fMassVar)) << ", m[0]: " << m[0] << endl;
        totalTemplates += pars[0]*fRelWeights[iTempl]*fHistoTemplates[iTempl]->getVal(RooArgSet(this->fMassVar));
      }
      break;
    case TemplAnchorMode::AnchorToSgn:
      for(int iTempl=0; iTempl<fHistoTemplates.size(); iTempl++) {
        // cout << "[MassTemplates, signal] fRelWeights[" << iTempl << "]: " << fRelWeights[iTempl] << ", pars[0]: " << pars[0] << ", eval templ: " << fHistoTemplates[iTempl]->getVal(RooArgSet(this->fMassVar)) << ", m[0]: " << m[0] << ", " << pars[0]*fRelWeights[iTempl]*fHistoTemplates[iTempl]->getVal(RooArgSet(this->fMassVar)) << endl;
        totalTemplates += pars[0]*fRelWeights[iTempl]*fHistoTemplates[iTempl]->getVal(RooArgSet(this->fMassVar));
      }
      break;
    default:
      std::cerr << "Error: Invalid fAnchorTemplsMode value!" << std::endl;
  }
  return totalTemplates;
}

//_________________________________________________________________________
Double_t VnVsMassFitter::VnTemplates(Double_t *m,Double_t *pars){
  // Add the contributions of the templates loaded in fHistoTemplates, each
  // scaled by a multiplicative constant, left as free fit parameter
  this->fMassVar.setVal(m[0]);  // Set the mass value to x
  Double_t totalTemplates = 0.;
  switch (fAnchorTemplsMode) {
    case TemplAnchorMode::Free:
      for(int iTempl=0; iTempl<fHistoTemplates.size(); iTempl++) {
        totalTemplates += pars[iTempl+fNParsTempls]*pars[iTempl]*fHistoTemplates[iTempl]->getVal(RooArgSet(this->fMassVar));
      }
      break;
    case TemplAnchorMode::AnchorToFirst:
      for(int iTempl=0; iTempl<fHistoTemplates.size(); iTempl++) {
        totalTemplates += pars[0]*fRelWeights[iTempl]*fHistoTemplates[iTempl]->getVal(RooArgSet(this->fMassVar));
      }
      break;
    case TemplAnchorMode::AnchorToSgn:
      // cout << "fNParsTempls: " << fNParsTempls << endl;
      // cout << "[VnTemplates] Anchoring to signal" << endl;
      for(int iTempl=0; iTempl<fHistoTemplates.size(); iTempl++) { 
        totalTemplates += pars[0]*fRelWeights[iTempl]*fHistoTemplates[iTempl]->getVal(RooArgSet(this->fMassVar));
      }
      break;
    default:
      std::cerr << "Error: Invalid fAnchorTemplsMode value!" << std::endl;
  }

  return totalTemplates;
}

//_________________________________________________________________________
Double_t VnVsMassFitter::MassBkgRfl(Double_t *m,Double_t *pars){

  if(!fHistoTemplRfl) {return MassBkg(m,pars);}
  else {
    //bkg mass parameters
    const Int_t nBkgPars = fNParsMassBkg;
    Double_t bkgpars[nBkgPars];
    for(Int_t iPar=0; iPar<fNParsMassBkg; iPar++) {bkgpars[iPar] = pars[iPar];}
    //reflection parameters
    Double_t rflpars[1]; //maximum number of parameters for rfl = 1 for the implemented functions
    for(Int_t iPar=0; iPar<fNParsRfl; iPar++) {rflpars[iPar] = pars[iPar+fNParsMassBkg];}
    return MassBkg(m,bkgpars)+MassRfl(m,rflpars);
  }
}

//_________________________________________________________________________
Double_t VnVsMassFitter::MassSecondPeak(Double_t *m,Double_t *pars){
  /// Fit function for a second gaussian peak
  /// To be used, e.g., for D+->KKpi in the Ds mass spectrum

  return pars[0]*GetGausPDF(m[0],pars[1],pars[2]);
}

//________________________________________________________________
Double_t VnVsMassFitter::vnBkgFunc(Double_t *m, Double_t *pars) {

  switch(fVnBkgFuncType) {
    case 0: //expo
      return pars[0]*GetExpoPDF(m[0],pars[1],kFALSE);
      break;
    case 1: //linear
      return GetPolPDF(m[0],pars,1,kFALSE);
      break;
    case 2: //parabolic
      return GetPolPDF(m[0],pars,2,kFALSE);
      break;
    case 6: //higher order (>=3) polinomial
      return GetHigherPolFuncPDF(m[0],pars,fPolDegreeVnBkg,kFALSE);
      break;
  }
  return 0;
}

//________________________________________________________________
Double_t VnVsMassFitter::MassFunc(Double_t *m, Double_t *pars) {

  //bkg mass parameters
  const Int_t nBkgPars = fNParsMassBkg;
  Double_t bkgpars[nBkgPars];
  for(Int_t iPar=0; iPar<fNParsMassBkg; iPar++) {bkgpars[iPar] = pars[iPar];}
  //signal mass parameters
  Double_t sgnpars[5]; //maximum number of parameters for sgn = 5 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsMassSgn; iPar++) {sgnpars[iPar] = pars[iPar+fNParsMassBkg];}
  //second peak parameters
  Double_t secpeakpars[3]; //maximum number of parameters for second peak = 3 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsSec; iPar++) {secpeakpars[iPar] = pars[iPar+fNParsMassBkg+fNParsMassSgn];}
  //reflection parameters
  Double_t rflpars[1]; //maximum number of parameters for rfl = 1 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsRfl; iPar++) {rflpars[iPar] = pars[iPar+fNParsMassBkg+fNParsMassSgn+fNParsSec];}

  Double_t total = MassSignal(m,sgnpars)+MassBkg(m,bkgpars);
  if(fSecondPeak) {total += MassSecondPeak(m,secpeakpars);}
  if(fReflections) {total += MassRfl(m,rflpars);}

  if (fTemplates) {
    switch (fAnchorTemplsMode) {
      case TemplAnchorMode::Free:
        total += MassTemplates(m,&pars[fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl]);
        break;
      case TemplAnchorMode::AnchorToFirst:
        total += MassTemplates(m,&pars[fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl]);
        break;
      case TemplAnchorMode::AnchorToSgn:
        total += MassTemplates(m,sgnpars);
        // cout << "[m[0] = " << m[0] << "] MassTemplates(m,sgnpars): " << MassTemplates(m,sgnpars) << endl;
        break;
      default:
        std::cerr << "Error: Invalid fAnchorTemplsMode value!" << std::endl;
     }
  }

  return total;
}

//________________________________________________________________
Double_t VnVsMassFitter::vnFunc(Double_t *m, Double_t *pars) {

  //bkg mass parameters
  const Int_t nBkgPars = fNParsMassBkg;
  Double_t massbkgpars[nBkgPars];
  for(Int_t iPar=0; iPar<fNParsMassBkg; iPar++) {massbkgpars[iPar] = pars[iPar];}
  //signal mass parameters
  Double_t masssgnpars[5]; //maximum number of parameters for mass sgn = 5 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsMassSgn; iPar++) {masssgnpars[iPar] = pars[iPar+fNParsMassBkg];}
  //second peak parameters
  Double_t secpeakpars[3]; //maximum number of parameters for second peak = 3 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsSec; iPar++) {secpeakpars[iPar] = pars[iPar+fNParsMassBkg+fNParsMassSgn];}
  //reflection parameters
  Double_t rflpars[1]; //maximum number of parameters for rfl = 1 for the implemented functions
  for(Int_t iPar=0; iPar<fNParsRfl; iPar++) {rflpars[iPar] = pars[iPar+fNParsMassBkg+fNParsMassSgn+fNParsSec];}
  
  //bkg vn parameters
  const Int_t nVnBkgPars = fNParsVnBkg;
  Double_t vnbkgpars[nVnBkgPars];
  for(Int_t iPar=0; iPar<fNParsVnBkg; iPar++) {
    vnbkgpars[iPar] = pars[iPar+fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsTempls];
  }
  //signal vn parameter
  Double_t vnSgn = pars[fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsTempls+fNParsVnBkg];
  //second peak vn parameter
  Double_t vnSecPeak = 0;
  if(fSecondPeak && fDoSecondPeakVn && !fFixVnSecPeakToSgn) {vnSecPeak = pars[fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsTempls+fNParsVnBkg+fNParsVnSgn];}
  if(fSecondPeak && fDoSecondPeakVn && fFixVnSecPeakToSgn) {vnSecPeak = pars[fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsTempls+fNParsVnBkg];}
  //refl vn parameter
  Double_t vnRefl = 0;
  if(fReflections) {
    switch(fVnRflOpt) {
      case 0:
        vnRefl = pars[fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsVnBkg];
        break;
      case 1:
        vnRefl = -pars[fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsVnBkg];
        break;
      case 2:
        vnRefl = 0; //not used
        break;
      case 3:
        vnRefl = pars[fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsVnBkg+fNParsVnSgn+fNParsVnSecPeak];
        break;
      default:
        printf("Error in setting reflection vn option: check fVnRflOpt");
        break;
    }
  }
  Double_t vnBkg = vnBkgFunc(m,vnbkgpars);
  Double_t Sgn = MassSignal(m,masssgnpars);
  Double_t Bkg = MassBkg(m,massbkgpars);
  Double_t SecPeak = 0;
  if(fSecondPeak) {
    if(fDoSecondPeakVn) SecPeak += MassSecondPeak(m,secpeakpars);
    else Bkg += MassSecondPeak(m,secpeakpars);
  }
  Double_t Refl=0;
  if(fReflections) {
    if(fVnRflOpt==kSameVnBkg) Bkg += MassRfl(m,rflpars);
    else Refl += MassRfl(m,rflpars);
  }

  Double_t TemplatesVn = 0;
  Double_t TemplatesMass = 0;
  if (fTemplates) {
    switch (fAnchorTemplsMode) {
      case TemplAnchorMode::Free:
        TemplatesMass += MassTemplates(m,&pars[fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl]);
        break;
      case TemplAnchorMode::AnchorToFirst:
        TemplatesMass += MassTemplates(m,&pars[fNParsMassBkg+fNParsMassSgn+fNParsSec+fNParsRfl]);
        break;
      case TemplAnchorMode::AnchorToSgn:
        // cout << "Anchoring to signal in vnFunc" << endl;
        TemplatesMass += MassTemplates(m,masssgnpars);
        break;
      default:
        std::cerr << "Error: Invalid fAnchorTemplsMode value!" << std::endl;
    }

    if(fTemplSameVnOfSignal){
      TemplatesVn += vnSgn*TemplatesMass;
    } else {
      TemplatesVn += VnTemplates(m,&pars[fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsTempls+fNParsVnBkg+fNParsVnSgn+fNParsVnSecPeak+fNParsRfl]);
    }

  }

  return (vnSgn*Sgn+vnBkg*Bkg+vnSecPeak*SecPeak+vnRefl*Refl+TemplatesVn)/(Sgn+Bkg+SecPeak+Refl+TemplatesMass);
}
