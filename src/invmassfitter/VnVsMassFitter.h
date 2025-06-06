#ifndef VNVSMASSFITTER_H
#define VNVSMASSFITTER_H
  /// \class VnVsMassFitter
  /// \class that performs the vn vs mass simultaneus fit for D mesons
  
#include <TObject.h>
#include <Riostream.h>
#include <TVirtualPad.h>
#include <TSpline.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TKDE.h>
#include "Fit/Fitter.h"
#include "Fit/Chi2FCN.h"
#include "Math/WrappedMultiTF1.h"
#include "InvMassFitter.h"

class VnVsMassFitter : public TObject {

public:
  VnVsMassFitter();
  VnVsMassFitter(TH1F* hMass, TH1F* hvn, Double_t min, Double_t max, Int_t funcMassBkg, Int_t funcMassSgn, Int_t funcvnBkg);
  ~VnVsMassFitter();

  enum ETypeOfBkg{kExpo=0, kLin=1, kPol2=2, kNoBk=3, kPow=4, kPowEx=5, kPoln=6};
  enum ETypeOfSgn{kGaus=0, k2Gaus=1, kDoubleCBAsymm=3, kDoubleCBSymm=4};
  enum ETypeOfVnRfl{kSameVnSignal=0, kOppVnSignal=1, kSameVnBkg=2, kFreePar=3};
  enum TemplAnchorMode{Free=0, AnchorToFirst=1, AnchorToSgn=2};

  Bool_t SimultaneousFit();

  //setters
  void SetInitialGaussianSigma(Double_t sigma, Int_t opt) {fSigmaInit=sigma; fSigmaFixed=opt;}
  void SetInitialGaussianMean(Double_t mean, Int_t opt) {fMeanInit=mean; fMeanFixed=opt;}
  void SetInitialGaussianSigma2Gaus(Double_t sigma, Int_t opt) {fSigma2GausInit=sigma; fSigma2GausFixed=opt;}
  void SetInitialFrac2Gaus(Double_t frac, Int_t opt) {fFrac2GausInit=frac; fFrac2GausFixed=opt;}
  void SetParticlePdgMass(Double_t mass){fMassParticle=mass;}
  void SetMassSgnFunc(Int_t functype) {fMassSgnFuncType=functype;}
  void SetMassBkgFunc(Int_t functype) {fMassBkgFuncType=functype;}
  void SetVnBkgFunc(Int_t functype) {fVnBkgFuncType=functype;}
  void FixSigmaFromMassFit() {fSigmaFixedFromMassFit=kTRUE;}
  void FixMeanFromMassFit() {fMeanFixedFromMassFit=kTRUE;}
  void FixSigma2GausFromMassFit() {fSigma2GausFixedFromMassFit=kTRUE;}
  void FixFrac2GausFromMassFit() {fFrac2GausFixedFromMassFit=kTRUE;}
  void SetNSigmaForVnSB(Int_t nsigma=4) {
    cout << "SetNSigmaForVnSB: " << nsigma << endl;
    fNSigmaForSB=nsigma;
  }
  void SetPolDegreeForBackgroundFit(Int_t deg){
    if(fMassBkgFuncType!=6) printf("fMassBkgFuncType should be set to 6 to use higher order polynomials\n");
    fPolDegreeBkg=deg;
  }
  void SetPolDegreeForVnBackgroundFit(Int_t deg){
    if(fVnBkgFuncType!=6) printf("fVnBkgFuncType should be set to 6 to use higher order polynomials\n");
    fPolDegreeVnBkg=deg;
  }
  void SetTemplateReflections(const TH1 *h, TString opt, Double_t minRange, Double_t maxRange) {
    fHistoTemplRflInit=(TH1F*)h->Clone();
    /// option could be:
    ///    "template"                use MC histograms
    ///    "1gaus" ot "singlegaus"   single gaussian function fit to MC templates
    ///    "2gaus" ot "doublegaus"   double gaussian function fit to MC templates
    ///    "pol3"                    3rd order polynomial fit to MC templates
    ///    "pol6"                    6th order polynomial fit to MC templates
    fRflOpt=opt;
    fMinRefl=minRange;
    fMaxRefl=maxRange;
    fReflections=kTRUE;
  }
  void SetTemplatesHisto(int anchormode, std::vector<Double_t> relcombweights, std::vector<std::string> templsnames, std::vector<const TH1*> histotempl,
                         std::vector<Double_t> initweights, std::vector<Double_t> minweights, std::vector<Double_t> maxweights, 
                         std::vector<Double_t> vninitweights, std::vector<Double_t> vnminweights, std::vector<Double_t> vnmaxweights, 
                         Bool_t samevnofsignal) {

    fTemplates=kTRUE;
    TFile* file = new TFile("templates_from_roofit.root", "RECREATE");

    for (int iTempl = 0; iTempl < histotempl.size(); ++iTempl) {
      file->mkdir(Form("Template_%i", iTempl));
      file->cd(Form("Template_%i", iTempl));

      // Set fMassVar range and binning to match histogram
      Double_t xmin = histotempl[iTempl]->GetXaxis()->GetXmin();
      Double_t xmax = histotempl[iTempl]->GetXaxis()->GetXmax();
      Int_t nbins = histotempl[iTempl]->GetNbinsX();

      this->fMassVar.setRange("fullRange", xmin, xmax);
      this->fMassVar.setBins(nbins);
      this->fMassVar.setMin(xmin);
      this->fMassVar.setMax(xmax);

      // Clone and normalize histogram to unit area (PDF style)
      TH1D* histPdf = (TH1D*)histotempl[iTempl]->Clone("histPdf");
      histPdf->Scale(1.0 / histPdf->Integral("width"));

      // Use the normalized histogram for RooDataHist
      RooDataHist* data_hist = new RooDataHist(Form("templ_%i", iTempl), Form("templ_%i", iTempl),
                                              RooArgList(this->fMassVar), histPdf);

      RooHistPdf* pdf = new RooHistPdf(Form("templ_%i_pdf", iTempl), Form("templ_%i_pdf", iTempl),
                                      RooArgSet(this->fMassVar), *data_hist);
      fHistoTemplates.push_back(pdf);

      // Check normalization of the RooHistPdf
      RooAbsReal* integral = pdf->createIntegral(RooArgSet(this->fMassVar), RooFit::NormSet(this->fMassVar));
      std::cout << "PDF integral over fMassVar (templ " << iTempl << "): " << integral->getVal() << std::endl;

      // Plot the PDF
      RooPlot* frame = this->fMassVar.frame();
      frame->SetName(Form("frame_%i", iTempl));
      pdf->plotOn(frame);

      TCanvas* c = new TCanvas(Form("canvas_%i", iTempl), Form("canvas_%i", iTempl), 800, 600);
      frame->Draw();

      // Save original and normalized histograms
      histotempl[iTempl]->Write();
      histPdf->Write();  // Proper PDF histogram

      // Optional: Save histogram with just normalized counts
      TH1D* histNormCounts = (TH1D*)histotempl[iTempl]->Clone("histNormCounts");
      histNormCounts->Scale(1.0 / histNormCounts->Integral());
      histNormCounts->Write();

      // Save canvas and RooFit objects
      c->Write();
      data_hist->Write();

      double xval = 1.751;  // example value in the domain of fMassVar
      this->fMassVar.setVal(xval);  // set the value to evaluate

      double pdfVal = pdf->getVal(RooArgSet(this->fMassVar));
      std::cout << "PDF value at mass = " << xval << " is: " << pdfVal << std::endl;

    }

    file->Close();
    std::cout << "SetTemplatesHisto VnVsMassFitter ended" << std::endl;

    fMassInitWeights=initweights;
    fMassWeightsLowerLims=minweights;
    fMassWeightsUpperLims=maxweights;
    fVnInitWeights=vninitweights;
    fVnWeightsLowerLims=vnminweights;
    fVnWeightsUpperLims=vnmaxweights;
    for(int iFunc=0; iFunc<fHistoTemplates.size(); iFunc++) {
      fHistoTemplates[iFunc]->SetName(Form("TemplFlag_%s", templsnames[iFunc].c_str()));
      fHistoTemplates[iFunc]->SetTitle(Form("TemplFlag_%s", templsnames[iFunc].c_str()));
    }
    if(samevnofsignal) {printf("WARNING: Vn parameter of templates will be the same as the one of the signal! \n");}
    fTemplSameVnOfSignal=samevnofsignal;
    fTemplates=kTRUE;
    fRelWeights=relcombweights;
    fAnchorTemplsMode=static_cast<TemplAnchorMode>(anchormode);

    // Set templates for fMassFitter
    cout << "Setting templates for mass fitter" << endl;
    fMassFitter->SetTemplates(anchormode, relcombweights, templsnames, histotempl, initweights, minweights, maxweights);
    cout << "Templates for mass fitter set" << endl;

  }
  void SetInitialReflOverS(Double_t rovers){fRflOverSig=rovers;}
  void SetFixReflOverS(Double_t rovers){
    SetInitialReflOverS(rovers);
    fFixRflOverSig=kTRUE;
  }
  void SetReflVnOption(Int_t opt) {fVnRflOpt=opt;}
  void SetReflVnParLimits(Double_t min, Double_t max) {
    fVnRflLimited=kTRUE;
    fVnRflMin=min;
    fVnRflMax=max;
  }
  void IncludeSecondGausPeak(Double_t mass, Bool_t fixm, Double_t width, Bool_t fixw, Bool_t doVn, Bool_t fixtosgn){
    fSecondPeak=kTRUE; fSecMass=mass; fSecWidth=width;
    fFixSecMass=fixm;  fFixSecWidth=fixw;
    fDoSecondPeakVn=doVn;
    fFixVnSecPeakToSgn=fixtosgn;
  }
  void SetInitPars(std::vector<std::tuple<TString, double, double, double>>  initFuncPars) {
    cout << "SetInitPars VnVsMassfitter" << endl;
    fInitFuncPars = initFuncPars;
  }
  void ApplyInitPars();
  void SetHarmonic(Int_t harmonic=2) {fHarmonic=harmonic;}

  // Double-sided crystal ball functions
  Double_t DoubleSidedCBAsymmForVn(double x, double mu, double width, double a1, double n1, double a2, double n2);
  Double_t DoubleSidedCBSymmForVn(double x, double mu, double width, double a, double n);

  TH1D *GetPullDistribution();

  //getters
  Double_t GetVn() const {return fVn;}
  Double_t GetVnUncertainty() const {return fVnUncertainty;}
  Double_t GetMean() const {return fMean;}
  Double_t GetMeanUncertainty() const {return fMeanUncertainty;}
  Double_t GetSigma() const {return fSigma;}
  Double_t GetSigmaUncertainty() const {return fSigmaUncertainty;}
  Double_t GetRawYield() const {return fRawYield;}
  Double_t GetRawYieldUncertainty() const {return fRawYieldUncertainty;}
  Double_t GetChiSquare() const {return fChiSquare;}
  Int_t GetNDF() const {return fNDF;}
  Double_t GetReducedChiSquare() const {return fChiSquare/fNDF;}
  Double_t GetFitProbability() const {return fProb;}
  Double_t GetSBVnPrefitChiSquare() const {return fSBVnPrefitChiSquare;}
  Int_t GetSBVnPrefitNDF() const {return fSBVnPrefitNDF;}
  Double_t GetSBVnPrefitReducedChiSquare() const {return fSBVnPrefitChiSquare/fSBVnPrefitNDF;}
  Double_t GetSBVnPrefitProbability() const {return fSBVnPrefitProb;}
  InvMassFitter* GetMassPrefitObject() const {return fMassFitter;}
  Double_t GetMassPrefitChiSquare() const {return fMassPrefitChiSquare;}
  Int_t GetMassPrefitNDF() const {return fMassPrefitNDF;}
  Double_t GetMassPrefitReducedChiSquare() const {return fMassPrefitChiSquare/fMassPrefitNDF;}
  Double_t GetMassPrefitProbability() const {return fMassPrefitProb;}
  Double_t GetParticlePdgMass() const {return fMassParticle;}
  TH1F* GetTemplateReflections() {
    if(fHistoTemplRfl) {return (TH1F*)fHistoTemplRfl->Clone("fHistoTemplRfl");}
    else if(fHistoTemplRflInit) {return (TH1F*)fHistoTemplRflInit->Clone("fHistoTemplRflInit");}
    else {return 0;}
  }
  void Signal(Double_t nOfSigma,Double_t &signal,Double_t &errsignal) const;
  void Signal(Double_t min,Double_t max,Double_t &signal,Double_t &errsignal) const;
  void Background(Double_t nOfSigma, Double_t &background,Double_t &errbackground) const;
  void Background(Double_t min, Double_t max, Double_t &background,Double_t &errbackground) const;
  void Significance(Double_t nOfSigma, Double_t &significance,Double_t &errsignificance) const;
  void Significance(Double_t min, Double_t max, Double_t &significance,Double_t &errsignificance) const;
  TF1* GetMassTotFitFunc() const {
    if(fMassTotFunc) return fMassTotFunc;
    else return nullptr;
  }
  Int_t GetNMassBkgPars() const {
    return fNParsMassBkg;
  }
  Int_t GetNMassSgnPars() const {
    return fNParsMassSgn;
  }
  Int_t GetNMassSecPeakPars() const {
    // only for gaussian case
    return 3;
  }
  Int_t GetNMassReflPars() const {
    // only for gaussian case
    return fNParsRfl;
  }
  Int_t GetNVnBkgPars() const {
    return fNParsVnBkg;
  }
  Int_t GetNVnSgnPars() const {
    return fNParsVnSgn;
  }
  Int_t GetNVnSecPeakPars() const {
    return fNParsVnSecPeak;
  }
  Int_t GetNVnReflPars() const {
    return fNParsVnRfl;
  }
  TF1* GetMassSignalFitFunc() const {
    if(fMassSgnFunc) return fMassSgnFunc;
    else return nullptr;
  }
  TF1* GetMassBkgFitFunc() const {
    if(fMassBkgFunc) return fMassBkgFunc;
    else return nullptr;
  }
  TF1* GetMassTemplFitFunc() const {
    if(fMassTemplFunc) return fMassTemplFunc;
    else return nullptr;
  }
  TF1* GetVnVsMassTotFitFunc() const {
    if(fVnTotFunc) return fVnTotFunc;
    else return nullptr;
  }
  TF1* GetVnVsMassBkgFitFunc() const {
    if(fVnBkgFunc) return fVnBkgFunc;
    else return nullptr;
  }
  TF1* GetMassRflFunc() const {
    if(fReflections) return fMassRflFunc;
    else return nullptr;
  }
  TF1* GetMassBkgRflFunc() const {
    if(fReflections) return fMassBkgRflFunc;
    else return nullptr;
  }
  TF1* GetMassSecPeakFunc() const {
    if(fSecondPeak) return fMassSecPeakFunc;
    else return nullptr;
  }
  TF1* GetVnSecPeakFunc() const {
    if(fSecondPeak) return fVnSecPeakFunc;
    else return nullptr;
  }
  std::vector<TF1*> GetMassTemplFuncts() const {
    if(fTemplates) return fMassTemplatesDraw;
    else return {};
  }
  double GetTemplOverSig() const {
    if(fMassTemplFunc && fMassSgnFunc) return fMassTemplFunc->Integral(this->fMassMin, this->fMassMax) / fMassSgnFunc->Integral(this->fMassMin, this->fMassMax);
    else return 0;
  }
  std::vector<TF1*> GetVnCompsFuncts() const {
    return fVnCompsDraw;
  }
  std::vector<double> GetVnTemplates() const {
    std::vector<double> vnPars;
    if(!fTemplSameVnOfSignal) {
      for(int iFunc=0; iFunc<fHistoTemplates.size(); iFunc++) {
        vnPars.push_back(fVnTotFunc->GetParameter(iFunc+fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsTempls+fNParsVnBkg+fNParsVnSgn+fNParsVnSecPeak+fNParsRfl));
      }
    } else {
      for(int iFunc=0; iFunc<fHistoTemplates.size(); iFunc++) {
        vnPars.push_back(GetVn());
      }
    }
    return vnPars;
  }
  std::vector<double> GetVnTemplatesUncertainties() const {
    std::vector<double> vnPars;
    if(!fTemplSameVnOfSignal) {
      for(int iFunc=0; iFunc<fHistoTemplates.size(); iFunc++) {
        vnPars.push_back(fVnTotFunc->GetParError(iFunc+fNParsMassSgn+fNParsMassBkg+fNParsSec+fNParsRfl+fNParsTempls+fNParsVnBkg+fNParsVnSgn+fNParsVnSecPeak+fNParsRfl));
      }
    } else {
      for(int iFunc=0; iFunc<fHistoTemplates.size(); iFunc++) {
        vnPars.push_back(GetVnUncertainty());
      }
    }
    return vnPars;
  }
  //struct for global chi2 (for simultaneus fit)
  struct GlobalChi2 {
    GlobalChi2(ROOT::Math::IMultiGenFunction & f1,ROOT::Math::IMultiGenFunction & f2) : fChi2_1(&f1), fChi2_2(&f2) {}

    double operator() (const double *par) const {
        return (*fChi2_1)(par) + (*fChi2_2)(par);
    }
    const  ROOT::Math::IMultiGenFunction * fChi2_1;
    const  ROOT::Math::IMultiGenFunction * fChi2_2;
  };

private:

    ///fit functions
  Double_t GetGausPDF(Double_t x, Double_t mean, Double_t sigma);
  Double_t GetExpoPDF(Double_t x, Double_t slope, Bool_t isnorm=kTRUE);
  Double_t GetPolPDF(Double_t x, Double_t *pars, Int_t order, Bool_t isnorm=kTRUE);
  Double_t GetPowerFuncPDF(Double_t x, Double_t *pars);
  Double_t GetPowerExpoPDF(Double_t x, Double_t *pars);
  Double_t GetHigherPolFuncPDF(Double_t x, Double_t *pars, Int_t Ndeg, Bool_t isnorm=kTRUE);
  Double_t MassSignal(Double_t *m, Double_t *pars);
  Double_t MassBkg(Double_t *m, Double_t *pars);
  Double_t MassRfl(Double_t *m,Double_t *par);
  Double_t MassBkgRfl(Double_t *m,Double_t *par);
  Double_t MassTemplates(Double_t *m,Double_t *pars);
  Double_t MassSecondPeak(Double_t *m,Double_t *par);
  Double_t MassFunc(Double_t *m, Double_t *pars);
  Double_t vnBkgFunc(Double_t *m, Double_t *pars);
  Double_t vnFunc(Double_t *m, Double_t *pars);
  Double_t VnTemplates(Double_t *m,Double_t *pars);

    ///private methods
  void DefineNumberOfParameters();
  Bool_t MassPrefit();
  Bool_t VnSBPrefit();
  void SetParNames();

    ///data members
  TH1F*                 fMassHisto;                     /// mass histogram to fit
  TH1F*                 fVnVsMassHisto;                 /// vn vs. mass histogram to fit
  Int_t                 fMassSgnFuncType;               /// type of mass signal fit function
  Int_t                 fMassBkgFuncType;               /// type of mass bkg fit function
  Int_t                 fVnBkgFuncType;                 /// type of vn bkg fit function
  TF1*                  fMassFuncFromPrefit;            /// mass fit function (1st step, from prefit)
  TF1*                  fMassBkgFunc;                   /// mass bkg fit function (final, after simultaneus fit)
  TF1*                  fMassSgnFunc;                   /// mass signal fit function (final, after simultaneus fit)
  TF1*                  fMassTemplFunc;                 /// mass signal fit function (final, after simultaneus fit)
  TF1*                  fMassTotFunc;                   /// mass fit function (final, after simultaneus fit)
  TF1*                  fVnBkgFuncSb;                   /// vn bkg fit function (1st step from SB prefit)
  TF1*                  fVnBkgFunc;                     /// vn bkg fit function (final, after simultaneus fit)
  TF1*                  fVnTotFunc;                     /// vn fit function (final, after simultaneus fit)
  InvMassFitter*        fMassFitter;                    /// mass fitter for mass prefit
  Double_t              fMassMin;                       /// upper mass limit
  Double_t              fMassMax;                       /// lower mass limit
  Double_t              fVn;                            /// vn of the signal from fit
  Double_t              fVnUncertainty;                 /// uncertainty on vn of the signal from simultaneus fit
  Double_t              fSigma;                         /// mass peak width from simultaneus fit
  Double_t              fSigmaUncertainty;              /// uncertainty on mass peak width from simultaneus fit
  Double_t              fMean;                          /// mass peak position from simultaneus fit
  Double_t              fMeanUncertainty;               /// uncertainty on mass peak position from simultaneus fit
  Double_t              fRawYield;                      /// raw yield from simultaneus fit
  Double_t              fRawYieldUncertainty;           /// uncertainty raw yield from simultaneus fit
  Double_t              fChiSquare;                     /// simultaneus fit chi square
  Int_t                 fNDF;                           /// simultaneus fit number of degree of freedom
  Double_t              fProb;                          /// simultaneus fit probability
  Double_t              fSBVnPrefitChiSquare;           /// vn SB prefit chi square
  Int_t                 fSBVnPrefitNDF;                 /// vn SB prefit number of degree of freedom
  Double_t              fSBVnPrefitProb;                /// vn SB prefit probability
  Double_t              fMassPrefitChiSquare;           /// Mass prefit chi square
  Int_t                 fMassPrefitNDF;                 /// Mass prefit number of degree of freedom
  Double_t              fMassPrefitProb;                /// Mass prefit probability
  Int_t                 fNSigmaForSB;                   /// number of sigma for sidebands region (vn bkg prefit)
  Double_t              fSigmaInit;                     /// initialization for peak width
  Double_t              fMeanInit;                      /// initialization for peak position
  Double_t              fSigma2GausInit;                /// initialization for second peak width in case of k2Gaus
  Double_t              fFrac2GausInit;                 /// initialization for fraction of second gaussian in case of k2Gaus
  Bool_t                fMeanFixedFromMassFit;          /// flag to fix peak position from mass prefit
  Bool_t                fSigmaFixedFromMassFit;         /// flag to fix peak width from mass prefit
  Bool_t                fSigma2GausFixedFromMassFit;    /// flag to fix second peak width from mass prefit in case of k2Gaus
  Bool_t                fFrac2GausFixedFromMassFit;     /// flag to fix fraction of second gaussian in case of k2Gaus
  Double_t              fMassParticle;                  /// mass of selected particle
  Int_t                 fNParsMassSgn;                  /// number of parameters in mass signal fit function
  Int_t                 fNParsMassBkg;                  /// number of parameters in mass bkg fit function
  Int_t                 fNParsVnBkg;                    /// number of parameters in vn bkg fit function
  Int_t                 fNParsVnSgn;                    /// number of parameters in vn sgn fit function (1)
  Int_t                 fNParsVnSecPeak;                /// number of parameters in vn sec peak fit function (1 if included, 0 otherwise)
  Int_t                 fNParsVnRfl;                    /// number of parameters in vn refl fit function (1 if included, 0 otherwise)
  Int_t                 fSigmaFixed;                    /// flag to fix peak width
  Int_t                 fMeanFixed;                     /// flag to fix peak position
  Int_t                 fSigma2GausFixed;               /// flag to fix second peak width in case of k2Gaus
  Int_t                 fFrac2GausFixed;                /// flag to fix fraction of second gaussian in case of k2Gaus
  Int_t                 fPolDegreeBkg;                  /// degree of polynomial expansion for back fit (option 6 for back)
  Int_t                 fPolDegreeVnBkg;                /// degree of polynomial expansion for vn back fit (option 6 for back)
  Bool_t                fReflections;                   /// flag use/not use reflections
  Int_t                 fNParsRfl;                      /// fit parameters in reflection fit function
  Double_t              fRflOverSig;                    /// reflection/signal
  Bool_t                fFixRflOverSig;                 /// switch for fix refl/signal
  TH1F*                 fHistoTemplRfl;                 /// histogram with reflection template
  TH1F*                 fHistoTemplRflInit;             /// initial histogram with reflection template
  TF1*                  fMassRflFunc;                   /// fit function for reflections
  TF1*                  fMassBkgRflFunc;                /// mass bkg fit function plus reflections (final, after simultaneus fit)
  TString               fRflOpt;                        /// refelction option
  Double_t              fMinRefl;                       /// minimum for refelction histo
  Double_t              fMaxRefl;                       /// maximum for refelction histo
  Bool_t                fSmoothRfl;                     /// switch for smoothing of reflection template
  Double_t              fRawYieldHelp;                  /// internal variable for fit with reflections
  Int_t                 fVnRflOpt;                      /// option for reflection vn type
  Bool_t                fVnRflLimited;                  /// flag to limit or not the vn of reflections
  Double_t              fVnRflMin;                      /// minimum vn of reflections
  Double_t              fVnRflMax;                      /// maximum vn of reflections
  Bool_t                fSecondPeak;                    /// switch off/on second peak (for D+->KKpi in Ds)
  TF1*                  fMassSecPeakFunc;               /// fit function for second peak
  TF1*                  fVnSecPeakFunc;                 /// fit function for second peak
  Int_t                 fNParsSec;                      /// number of parameters in second peak fit function
  Double_t              fSecMass;                       /// position of the 2nd peak
  Double_t              fSecWidth;                      /// width of the 2nd peak
  Bool_t                fFixSecMass;                    /// flag to fix the position of the 2nd peak
  Bool_t                fFixSecWidth;                   /// flag to fix the width of the 2nd peak
  Double_t              fVnSecPeak;                     /// vn of second peak from fit
  Bool_t                fDoSecondPeakVn;                /// flag to introduce second peak vn in the vn vs. mass fit
  Bool_t                fFixVnSecPeakToSgn;             /// flag to fix the vn of the second peak to the one of signal
  Double_t              fVnSecPeakUncertainty;          /// vn uncertainty of second peak from fit
  Int_t                 fHarmonic;                      /// harmonic number for drawing
  Bool_t                fTemplates;                     /// flag use/not use templates
  Int_t                 fNParsTempls;                   /// fit parameters to include templates
  std::vector<TF1 *>    fVnCompsDraw;                   /// vector to store TKDE to be added as templates to the fit function 
  std::vector<TF1 *>    fMassTemplatesDraw;          /// vector to store TKDE to be added as templates to the fit function 
  std::vector<Double_t> fRelWeights;                    /// relative weights of templates 
  std::vector<Double_t> fMassWeightsUpperLims;          /// upper limit of the templates' weights
  std::vector<Double_t> fMassWeightsLowerLims;          /// lower limit of the templates' weights
  std::vector<Double_t> fVnWeightsUpperLims;            /// upper limit of the templates' weights
  std::vector<Double_t> fVnWeightsLowerLims;            /// lower limit of the templates' weights
  std::vector<Double_t> fMassInitWeights;               /// init values of the templates' weights
  std::vector<Double_t> fVnInitWeights;                 /// init values of the templates' weights
  Bool_t                fTemplSameVnOfSignal;           /// init values of the templates' weights
  TemplAnchorMode       fAnchorTemplsMode;              /// init values of the templates' weights
  std::vector<std::tuple<TString, double, double, double>> fInitFuncPars;  /// init values of total fit function
  RooRealVar fMassVar;
  std::vector<RooHistPdf*> fHistoTemplates;  /// vector to store TKDE to be added as templates to the fit function

    /// \cond CLASSDEF
  ClassDef(VnVsMassFitter,5);
    /// \endcond
};
#endif //VNVSMASSFITTER
