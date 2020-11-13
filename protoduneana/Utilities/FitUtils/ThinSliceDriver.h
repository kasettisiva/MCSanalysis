#ifndef THINSLICEDRIVER_hh
#define THINSLICEDRIVER_hh

#include <map>
#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"

#include "ThinSliceSample.h"
#include "ThinSliceDataSet.h"

#include "fhiclcpp/ParameterSet.h"

namespace protoana {
class ThinSliceDriver {
 public:
  ThinSliceDriver(const fhicl::ParameterSet & extra_options);
  virtual ~ThinSliceDriver();
  virtual void BuildDataHists(
      TTree * tree, TH1D & incident_hist,
      std::map<int, TH1 *> & selected_hists) = 0;
  virtual void BuildMCSamples(
      TTree * tree,
      std::map<int, std::vector<ThinSliceSample>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::map<int, double> & nominal_fluxes,
      std::map<int, std::vector<double>> & fluxes_by_sample) = 0;
  virtual std::pair<double, size_t> CalculateChi2(
      std::map<int, std::vector<ThinSliceSample>> & samples,
      ThinSliceDataSet & data_set) = 0;

  virtual void CompareSelections(
      std::map<int, std::vector<ThinSliceSample>> & samples,
      ThinSliceDataSet & data_set,
      TFile & output_file,
      std::vector<std::pair<int, int>> plot_style,
      bool plot_rebinned,
      bool post_fit) = 0;

  void CompareDataMC(
      std::map<int, std::vector<ThinSliceSample>> & samples,
      ThinSliceDataSet & data_set,
      TFile & output_file,
      std::vector<std::pair<int, int>> plot_style,
      bool plot_rebinned = false,
      bool post_fit = false);
  /*virtual void BuildAndSaveStacks(
      std::map<int, std::vector<ThinSliceSample>> & samples,
      ThinSliceDataSet & data_set,
      TFile & output_file,
      THStack * incident_stack, 
      std::map<int, THStack *> * selection_stacks_x,
      std::map<int, THStack *> * selection_stacks_y = 0x0,
      std::map<int, THStack *> * selection_stacks_z = 0x0,
      bool post_fit = false) = 0;*/

  std::pair<int, int> GetColorAndStyle(
      size_t i, const std::vector<std::pair<int, int>> & plot_style);
 protected:
  fhicl::ParameterSet fExtraOptions;
 private:
};
}
#endif
