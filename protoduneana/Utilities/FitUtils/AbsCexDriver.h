#ifndef ABSCEXDRIVER_hh
#define ABSCEXDRIVER_hh

#include "ThinSliceDriver.h"
#include "TH2D.h"
#include "TFile.h"
#include "TSpline.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TRandom3.h"
#include <map>
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"

namespace protoana {
class AbsCexDriver : public ThinSliceDriver {
 public:
  AbsCexDriver(const fhicl::ParameterSet & extra_options);
  virtual ~AbsCexDriver();

  void BuildDataHists(
    TTree * tree, ThinSliceDataSet & data_set, double & flux,
    int split_val = 0) override;
  void BuildFakeData(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0) override;
  void FakeDataSampleScales(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);
  void FakeDataBinnedScales(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);
  void FakeDataG4RW(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);
  void FakeDataG4RWGrid(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);
  void FakeDataEffVar(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);

  void FakeDatadEdX(
    TTree * tree,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);

  void BuildMCSamples(
      //TTree * tree,
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::map<int, double> & nominal_fluxes,
      std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
      std::vector<double> & beam_energy_bins) override;

  void RefillMCSamples(
      //TTree * tree,
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins,
      const std::map<int, std::vector<double>> & signal_pars,
      const std::map<int, double> & flux_pars,
      const std::map<std::string, ThinSliceSystematic> & syst_pars,
      bool fill_incident = false) override;

  /*void BuildSystSamples(
      TTree * tree,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins) override;*/
  
  void SetupSyst_G4RW(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);

  void SetupSyst_dEdX_Cal(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);

  std::pair<double, size_t> CalculateChi2(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set) override;
  void CompareSelections(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set,
      TFile & output_file,
      std::vector<std::pair<int, int>> plot_style,
      bool plot_rebinned,
      bool post_fit, int nPars,
      TDirectory * plot_dir) override;

  void GetCurrentHists(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set,
      std::map<int, std::vector<TH1*>> & throw_hists,
      bool plot_rebinned) override;

  virtual void GetCurrentTruthHists(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      std::map<int, std::vector<TH1*>> & hists,
      std::map<int, std::vector<TH1*>> & inc_hists,
      std::map<int, std::vector<TH1*>> & xsec_hists,
      const std::vector<int> & incident_samples,
      const std::map<int, std::vector<double>> & signal_bins) override;

  void PlotThrows(
      ThinSliceDataSet & data_set, std::map<int, std::vector<TH1*>> & throw_hists,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      size_t nThrows,
      std::map<int, std::vector<TH1*>> & truth_throw_hists,
      std::map<int, std::vector<TH1*>> & truth_inc_hists,
      std::map<int, std::vector<TH1*>> & truth_xsec_hists,
      std::map<int, TH1*> & best_fit_incs,
      std::map<int, TH1*> & best_fit_xsecs,
      std::map<int, TH1*> & nominal_incs,
      std::map<int, TH1*> & nominal_xsecs,
      TFile & output_file, bool plot_rebinned,
      std::map<int, std::vector<double>> * sample_scales = 0x0) override;

  void SetupSysts(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file) override;

  /*void SetupSyst_BeamRes(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);*/
  void SetupSyst_BeamShift(
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);
  void SetupSyst_BeamShift2D(
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);
  void SetupSyst_EffVar(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);
  void SetupSyst_EffVarWeight(
      const std::map<std::string, ThinSliceSystematic> & pars);
  void SetupSyst_EDivWeight(
      const std::map<std::string, ThinSliceSystematic> & pars);

  /*double GetSystWeight_BeamRes(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);*/
  double GetSystWeight_BeamShift(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_BeamShift2D(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_G4RW(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars,
      const ThinSliceSample & sample,
      int selection_ID, double val);
  double GetSystWeight_EffVar(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_EDiv(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  void WrapUpSysts(TFile & output_file) override;
  int RecalculateSelectionID(
      const ThinSliceEvent & event,
      double C_cal,
      TProfile * prot_template);
  double TruncatedMean(const std::vector<double> & dEdX);
 private:
   TH2D * fEndSlices;
   TFile * fIn;
   std::map<int, double> fMeans;

   double fEnergyFix;
   bool fDoEnergyFix;

   double fPitch;
   double fZ0;
   double fEndZCut;
   double fTrajZStart;
   std::string fSliceMethod;
   int fSliceCut;

   double fBetaP, fRho, fWion, fAlpha, fNominalCCal;

   //bool fStaticBeamResWidth = false;
   //bool fStaticBeamResMean = false;
   //double fBeamResMeanVal = 1.;
   //double fBeamResWidthVal = 1.;
   TTree /** fSystBeamResTree, */* fSystBeamShiftTree, * fSystBeamShift2DTree;
   //double fSystBeamResWeight, fSystBeamResMeanOutput, fSystBeamResWidthOutput;
   //double fSystBeamResWeightCap, fSystBeamResOutput;
   double fSystBeamShiftWeight, fSystBeamShiftVal, fSystBeamShiftR;
   bool /*fSetupSystBeamRes = false,*/ fSetupSystBeamShift = false,
        fSetupSystBeamShift2D = false, fSetupSystEffVar = false,
        fSystBeamShiftTreeSave = false;
   double fSystBeamShift2DWeight, fSystBeamShift2DBVal, fSystBeamShift2DVal,
          fSystBeamShift2DR;
  // double fEffVarSystVal;
   TGraph2D * fSystBeamShiftMap, * fSystBeam2DMeans, * fSystBeam2DStdDevs;
   TGraph * fSystBeamShiftMeans, * fSystBeamShiftWidths;
   std::pair<double, double> fSystBeamShiftLimits;

   std::map<std::string, std::map<int, std::vector<TH1D*>>> fFullSelectionVars;
   std::map<std::string, std::map<int, std::vector<TSpline3*>>> fFullSelectionSplines;

   std::map<std::string, std::map<int, std::vector<TH1D*>>> fG4RWSelectionVarsPlus;
   std::map<std::string, std::map<int, std::vector<TH1D*>>> fG4RWSelectionVarsMinus;
   std::vector<std::string> fActiveG4RWSysts;
   TRandom3 fRNG = TRandom3(0);

   double fEffVarF, fEffVarCut;
   double fEDivF, fEDivCut;
   ProtoDUNETrackUtils fTrackUtil;
};
}
#endif
