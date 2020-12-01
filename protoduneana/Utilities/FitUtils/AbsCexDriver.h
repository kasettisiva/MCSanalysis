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
    TTree * tree, std::map<int, std::vector<ThinSliceSample>> & samples,
    ThinSliceDataSet & data_set, double & flux) override;
  void FakeDataSampleScales(
    TTree * tree, std::map<int, std::vector<ThinSliceSample>> & samples,
    ThinSliceDataSet & data_set, double & flux);
  void BuildMCSamples(
      TTree * tree,
      std::map<int, std::vector<ThinSliceSample>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::map<int, double> & nominal_fluxes,
      std::map<int, std::vector<double>> & fluxes_by_sample/*,
      std::map<int, std::pair<TH1D, TH1D>> & signal_eff_parts*/) override;

  void BuildSystSamples(
      TTree * tree,
      std::map<int, std::vector<ThinSliceSample>> & samples,
      const std::map<int, bool> & signal_sample_checks/*,
      std::map<int, double> & nominal_fluxes,
      std::map<int, std::vector<double>> & fluxes_by_sample*/) override;
  
  void SystRoutineG4RW(
      TTree * tree,
      std::map<int, std::vector<ThinSliceSample>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      const fhicl::ParameterSet & routine);

  std::pair<double, size_t> CalculateChi2(
      std::map<int, std::vector<ThinSliceSample>> & samples,
      ThinSliceDataSet & data_set) override;
  void CompareSelections(
      std::map<int, std::vector<ThinSliceSample>> & samples,
      ThinSliceDataSet & data_set,
      TFile & output_file,
      std::vector<std::pair<int, int>> plot_style,
      bool plot_rebinned,
      bool post_fit) override;
};
}
#endif
