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
                              "") + ";Reconstructed KE (MeV)";
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
  MakeRebinnedHists();
}

void protoana::ThinSliceSample::MakeRebinnedHists() {
  if (!fMadeRebinned) {
    std::string inc_name = fIncidentHist.GetName();
    inc_name += "Rebinned";
    fIncidentHistRebinned = TH1D(inc_name.c_str(), fIncidentHist.GetTitle(),
                                 fIncidentHist.GetNbinsX(), 0, fIncidentHist.GetNbinsX());
    for (int i = 1; i <= fIncidentHist.GetNbinsX(); ++i) {
      fIncidentHistRebinned.SetBinContent(i, fIncidentHist.GetBinContent(i));

      double low_edge = fIncidentHist.GetXaxis()->GetBinLowEdge(i);
      double up_edge = fIncidentHist.GetXaxis()->GetBinUpEdge(i);
      std::string bin_label = (low_edge < 0. ? "< 0." :
                               (protoana::PreciseToString(low_edge, 0) + " - " +
                                protoana::PreciseToString(up_edge, 0)));
      fIncidentHistRebinned.GetXaxis()->SetBinLabel(i, bin_label.c_str());
    }

    for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
      TH1D & sel_hist = it->second;
      std::string name = sel_hist.GetName();
      name += "Rebinned";
      fSelectionHistsRebinned[it->first] = TH1D(
          name.c_str(), sel_hist.GetTitle(), sel_hist.GetNbinsX(), 0,
          sel_hist.GetNbinsX());
      for (int i = 1; i <= sel_hist.GetNbinsX(); ++i) {
        fSelectionHistsRebinned[it->first].SetBinContent(
            i, sel_hist.GetBinContent(i));
        double low_edge = sel_hist.GetXaxis()->GetBinLowEdge(i);
        double up_edge = sel_hist.GetXaxis()->GetBinUpEdge(i);
        std::string bin_label = (low_edge < 0. ? "< 0." :
                                 (protoana::PreciseToString(low_edge, 0) + " - " +
                                  protoana::PreciseToString(up_edge, 0)));
        fSelectionHistsRebinned[it->first].GetXaxis()->SetBinLabel(
            i, bin_label.c_str());
      }
    }

    fMadeRebinned = true;
  }
}

void protoana::ThinSliceSample::RefillRebinnedHists() {
  for (int i = 1; i <= fIncidentHist.GetNbinsX(); ++i) {
    fIncidentHistRebinned.SetBinContent(i, fIncidentHist.GetBinContent(i));
  }

  for (auto it = fSelectionHistsRebinned.begin();
       it != fSelectionHistsRebinned.end(); ++it) {
    for (int i = 1; i <= it->second.GetNbinsX(); ++i) {
      it->second.SetBinContent(i, fSelectionHists[it->first].GetBinContent(i));
    }
  }
}
