#include "ThinSliceSample.h"

std::string protoana::PreciseToString(const double val, const int n) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << val;
  return out.str();
}

protoana::ThinSliceSample::ThinSliceSample(
    std::string name, int flux_type,
    const std::map<int, std::string> & selections,
    const std::vector<double> & incident_bins,
    const std::vector<double> & selected_bins,
    bool is_signal, std::pair<double, double> range)
    : fSampleName(name),
      fFluxType(flux_type),
      fIsSignal(is_signal),
      fRange(range) {

  std::string inc_name = "";
  std::string title = name + (is_signal ?
                              ("(" + protoana::PreciseToString(range.first) +
                               " " + protoana::PreciseToString(range.second) + ")") :
                              "");
  if (is_signal) {
    inc_name = "sample_" + name + "_" + protoana::PreciseToString(range.first) + "_" +
               protoana::PreciseToString(range.second) + "_incident_hist";
  }
  else {
    inc_name = "sample_" + name + "_incident_hist";
  }
  fIncidentHist = TH1D(inc_name.c_str(), title.c_str(), incident_bins.size() - 1,
                       &incident_bins[0]);

  for (auto it = selections.begin(); it != selections.end(); ++it) {
    std::string sel_name = "";
    if (is_signal) {
      sel_name = "sample_" + name + "_" + protoana::PreciseToString(range.first) + "_" +
                 protoana::PreciseToString(range.second) + "_selected_" + it->second +
                 "_hist";
    }
    else {
      sel_name = "sample_" + name + "_selected_" +
                 it->second + "_hist";
    }
    fSelectionHists[it->first] = TH1D(sel_name.c_str(), title.c_str(),
                                      selected_bins.size() - 1,
                                      &selected_bins[0]);
  }
}
