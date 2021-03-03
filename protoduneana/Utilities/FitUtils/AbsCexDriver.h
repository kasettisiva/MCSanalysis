#ifndef ABSCEXDRIVER_hh
#define ABSCEXDRIVER_hh

#include "ThinSliceDriver.h"
namespace protoana {
class AbsCexDriver : public ThinSliceDriver {
 public:
  AbsCexDriver(const fhicl::ParameterSet & extra_options);
  virtual ~AbsCexDriver();

  void BuildDataHists(
    TTree * tree, ThinSliceDataSet & data_set, double & flux) override;
  void BuildFakeData(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales) override;
  void FakeDataSampleScales(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales);
  void FakeDataBinnedScales(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales);
  void FakeDataG4RW(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales);

  void BuildMCSamples(
      TTree * tree,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::map<int, double> & nominal_fluxes,
      std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
      std::vector<double> & beam_energy_bins/*,
      std::map<int, std::pair<TH1D, TH1D>> & signal_eff_parts*/) override;

  void BuildSystSamples(
      TTree * tree,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks/*,
      std::map<int, double> & nominal_fluxes,
      std::map<int, std::vector<double>> & fluxes_by_sample*/) override;
  
  void SystRoutineG4RW(
      TTree * tree,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      const fhicl::ParameterSet & routine);

  std::pair<double, size_t> CalculateChi2(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set) override;
  void CompareSelections(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set,
      TFile & output_file,
      std::vector<std::pair<int, int>> plot_style,
      bool plot_rebinned,
      bool post_fit) override;

  void GetCurrentHists(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set,
      std::map<int, std::vector<TH1*>> & throw_hists,
      bool plot_rebinned) override;

  virtual void GetCurrentTruthHists(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      std::map<int, std::vector<TH1*>> & hists) override;

  void PlotThrows(
    ThinSliceDataSet & data_set, std::map<int, std::vector<TH1*>> & throw_hists,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    std::map<int, std::vector<TH1*>> & truth_throw_hists,
    TFile & output_file, bool plot_rebinned,
    std::map<int, std::vector<double>> * sample_scales = 0x0) override;

  /*void PostFitThrows(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set,
      TH1D & pars,
      TH2D & cov,
      TFile & output_file,
      std::vector<std::pair<int, int>> plot_style,
      bool plot_rebinned) override;*/
};
}
#endif
