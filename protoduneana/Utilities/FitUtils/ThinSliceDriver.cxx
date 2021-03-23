#include "ThinSliceDriver.h"

#include "THStack.h"
#include "TCanvas.h"

protoana::ThinSliceDriver::ThinSliceDriver(
    const fhicl::ParameterSet & extra_options)
    : fExtraOptions(extra_options) {}

protoana::ThinSliceDriver::~ThinSliceDriver() {}

void protoana::ThinSliceDriver::CompareDataMC(
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    ThinSliceDataSet & data_set, TFile & output_file,
    std::vector<std::pair<int, int>> plot_style,
    int nPars,
    bool plot_rebinned,
    bool post_fit) {

  for (auto it = samples.begin(); it != samples.end(); ++it) {
    std::vector<ThinSliceSample> & vec = it->second[0];
    for (size_t i = 0; i < vec.size(); ++i) {
      vec[i].RefillRebinnedHists();

    }

    for (size_t i = 1; i < it->second.size(); ++i) {
      std::vector<ThinSliceSample> & vec = it->second[i];
      for (size_t j = 0; j < vec.size(); ++j) {
        vec[j].RefillRebinnedHists();
      }
    }
  }

  CompareSelections(
      samples, data_set, output_file, plot_style, plot_rebinned,
      post_fit, nPars);
}


std::pair<int, int> protoana::ThinSliceDriver::GetColorAndStyle(
    size_t i,
    const std::vector<std::pair<int, int>> & plot_style) {
  return {plot_style.at(i % plot_style.size()).first,
          (i < plot_style.size() ? 1001: 3244)};
}
