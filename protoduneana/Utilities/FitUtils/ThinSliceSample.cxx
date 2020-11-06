#include "ThinSliceSample.h"

protoana::ThinSliceSample::ThinSliceSample(
    std::string name, const std::map<int, std::string> & selections,
    const std::vector<double> & incident_bins,
    const std::vector<double> & selected_bins,
    bool is_signal, std::pair<double, double> range)
    : fSampleName(name),
      fIsSignal(is_signal),
      fRange(range) {

  std::string inc_name = "";
  if (is_signal) {
    inc_name = "sample_" + name + "_" + std::to_string(range.first) + "_" +
               std::to_string(range.second) + "_incident_hist";
  }
  else {
    inc_name = "sample_" + name + "_incident_hist";
  }
  fIncidentHist = TH1D(inc_name.c_str(), "", incident_bins.size() - 1,
                       &incident_bins[0]);

  for (auto it = selections.begin(); it != selections.end(); ++it) {
    std::string sel_name = "";
    if (is_signal) {
      sel_name = "sample_" + name + "_" + std::to_string(range.first) + "_" +
                 std::to_string(range.second) + "_selected_" + it->second +
                 "_hist";
    }
    else {
      sel_name = "sample_" + name + "_selected_" +
                 it->second + "_hist";
    }
    fSelectionHists[it->first] = TH1D(sel_name.c_str(), "",
                                      selected_bins.size() - 1,
                                      &selected_bins[0]);
  }
}
