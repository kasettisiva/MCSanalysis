#include "ThinSliceSample.h"

std::string protoana::PreciseToString(const double val, const int n) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << val;
  return out.str();
}

protoana::ThinSliceSample::ThinSliceSample(
    std::string name, int flux_type,
    const std::vector<fhicl::ParameterSet> & selections,
    const std::vector<double> & incident_bins,
    const std::vector<double> & true_incident_bins,
    bool is_signal, std::pair<double, double> range)
    : fSampleName(name),
      fFluxType(flux_type),
      fIsSignal(is_signal),
      fRange(range) {

  std::string inc_name = "";
  std::string title = name +
                      (is_signal ?
                          ("(" + protoana::PreciseToString(range.first) + " "
                           + protoana::PreciseToString(range.second) + ")") :
                           "") +
                      ";Reconstructed KE (MeV)";
  if (is_signal) {
    inc_name = "sample_" + name + "_" +
               protoana::PreciseToString(range.first) + "_" +
               protoana::PreciseToString(range.second) + "_incident_hist";
  }
  else {
    inc_name = "sample_" + name + "_incident_hist";
  }
  fIncidentHist = TH1D(inc_name.c_str(), title.c_str(), incident_bins.size() - 1,
                       &incident_bins[0]);

  inc_name += "_true";
  fTrueIncidentHist = TH1D(inc_name.c_str(), title.c_str(),
                           true_incident_bins.size() - 1,
                           &true_incident_bins[0]);

  for (auto it = selections.begin(); it != selections.end(); ++it) {
    std::string sel_name = "";
    if (is_signal) {
      sel_name = "sample_" + name + "_" +
                 protoana::PreciseToString(range.first) + "_" +
                 protoana::PreciseToString(range.second) + "_selected_" +
                 it->get<std::string>("Name") + "_hist";
    }
    else {
      sel_name = "sample_" + name + "_selected_" +
                 it->get<std::string>("Name") + "_hist";
    }

    std::vector<std::vector<double>> selected_bins
        = it->get<std::vector<std::vector<double>>>("RecoBins");
    if (selected_bins.size() == 1) {
      fSelectionHists[it->get<int>("ID")] = new TH1D(
          sel_name.c_str(), title.c_str(), selected_bins[0].size() - 1,
          &selected_bins[0][0]);
    }
    else if (selected_bins.size() == 2) {
      fSelectionHists[it->get<int>("ID")] = new TH2D(
          sel_name.c_str(), title.c_str(),
          selected_bins[0].size() - 1, &selected_bins[0][0],
          selected_bins[1].size() - 1, &selected_bins[1][0]);
    }
    else if (selected_bins.size() == 3) {
      fSelectionHists[it->get<int>("ID")] = new TH3D(
          sel_name.c_str(), title.c_str(),
          selected_bins[0].size() - 1, &selected_bins[0][0],
          selected_bins[1].size() - 1, &selected_bins[1][0],
          selected_bins[2].size() - 1, &selected_bins[2][0]);
    }
    /*else {
      throw
    }*/
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
      TH1 * sel_hist = (TH1 *)it->second;
      std::string name = sel_hist->GetName();
      name += "Rebinned";
      
      size_t nAxes = 1;
      if (sel_hist->GetNbinsY() > 1) ++nAxes;
      if (sel_hist->GetNbinsZ() > 1) ++nAxes;

      if (nAxes == 1) {
        fSelectionHistsRebinned[it->first] = new TH1D(
            name.c_str(), sel_hist->GetTitle(),
            sel_hist->GetNbinsX(), 0, sel_hist->GetNbinsX());
        Rebin1D(sel_hist, fSelectionHistsRebinned[it->first]);
      }
      else if (nAxes == 2) {
        fSelectionHistsRebinned[it->first] = new TH2D(
            name.c_str(), sel_hist->GetTitle(),
            sel_hist->GetNbinsX(), 0, sel_hist->GetNbinsX(),
            sel_hist->GetNbinsY(), 0, sel_hist->GetNbinsY());
        Rebin2D(sel_hist, fSelectionHistsRebinned[it->first]);
      }
      else if (nAxes == 3) {
        fSelectionHistsRebinned[it->first] = new TH3D(
            name.c_str(), sel_hist->GetTitle(),
            sel_hist->GetNbinsX(), 0, sel_hist->GetNbinsX(),
            sel_hist->GetNbinsY(), 0, sel_hist->GetNbinsY(),
            sel_hist->GetNbinsZ(), 0, sel_hist->GetNbinsZ());
        Rebin3D(sel_hist, fSelectionHistsRebinned[it->first]);
      }
    }

    fMadeRebinned = true;
  }
}

void protoana::ThinSliceSample::Rebin1D(TH1 * sel_hist, TH1 * rebinned) {
  for (int i = 1; i <= sel_hist->GetNbinsX(); ++i) {
    double low_x = sel_hist->GetXaxis()->GetBinLowEdge(i);
    double up_x = sel_hist->GetXaxis()->GetBinUpEdge(i);
    std::string bin_label = (low_x < 0. ? "< 0." :
                             (protoana::PreciseToString(low_x, 0) + " - " +
                              protoana::PreciseToString(up_x, 0)));
    rebinned->GetXaxis()->SetBinLabel(i, bin_label.c_str());

    rebinned->SetBinContent(i, sel_hist->GetBinContent(i));
  }
}

void protoana::ThinSliceSample::Rebin2D(TH1 * sel_hist, TH1 * rebinned) {
  for (int i = 1; i <= sel_hist->GetNbinsX(); ++i) {
    double low_x = sel_hist->GetXaxis()->GetBinLowEdge(i);
    double up_x = sel_hist->GetXaxis()->GetBinUpEdge(i);
    std::string bin_label = (low_x < 0. ? "< 0." :
                             (protoana::PreciseToString(low_x, 0) + " - " +
                              protoana::PreciseToString(up_x, 0)));
    rebinned->GetXaxis()->SetBinLabel(
        i, bin_label.c_str());
    for (int j = 1; j <= sel_hist->GetNbinsY(); ++j) {
      double low_y = sel_hist->GetYaxis()->GetBinLowEdge(j);
      double up_y = sel_hist->GetYaxis()->GetBinUpEdge(j);
      std::string y_label = (low_y < 0. ? "< 0." :
                             (protoana::PreciseToString(low_y, 0) + " - " +
                              protoana::PreciseToString(up_y, 0)));
      rebinned->GetYaxis()->SetBinLabel(j, bin_label.c_str());
      rebinned->SetBinContent(i, j, sel_hist->GetBinContent(i, j));
    }
  }
}

void protoana::ThinSliceSample::Rebin3D(TH1 * sel_hist, TH1 * rebinned) {
  for (int i = 1; i <= sel_hist->GetNbinsX(); ++i) {
    double low_x = sel_hist->GetXaxis()->GetBinLowEdge(i);
    double up_x = sel_hist->GetXaxis()->GetBinUpEdge(i);
    std::string bin_label = (low_x < 0. ? "< 0." :
                             (protoana::PreciseToString(low_x, 0) + " - " +
                              protoana::PreciseToString(up_x, 0)));
    rebinned->GetXaxis()->SetBinLabel(i, bin_label.c_str());
    for (int j = 1; j <= sel_hist->GetNbinsY(); ++j) {
      double low_y = sel_hist->GetYaxis()->GetBinLowEdge(j);
      double up_y = sel_hist->GetYaxis()->GetBinUpEdge(j);
      std::string y_label = (low_y < 0. ? "< 0." :
                             (protoana::PreciseToString(low_y, 0) + " - " +
                              protoana::PreciseToString(up_y, 0)));
      rebinned->GetYaxis()->SetBinLabel(j, bin_label.c_str());

      for (int k = 1; k <= sel_hist->GetNbinsY(); ++k) {
        double low_z = sel_hist->GetYaxis()->GetBinLowEdge(k);
        double up_z = sel_hist->GetYaxis()->GetBinUpEdge(k);
        std::string y_label = (low_z < 0. ? "< 0." :
                               (protoana::PreciseToString(low_z, 0) + " - " +
                                protoana::PreciseToString(up_z, 0)));
        rebinned->GetZaxis()->SetBinLabel(k, bin_label.c_str());

        rebinned->SetBinContent(i, j, k, sel_hist->GetBinContent(i, j, k));
      }
    }
  }
}

void protoana::ThinSliceSample::RefillRebinnedHists() {
  for (int i = 1; i <= fIncidentHist.GetNbinsX(); ++i) {
    fIncidentHistRebinned.SetBinContent(i, fIncidentHist.GetBinContent(i));
  }

  for (auto it = fSelectionHistsRebinned.begin();
       it != fSelectionHistsRebinned.end(); ++it) {
    for (int i = 1; i <= it->second->GetNbinsX(); ++i) {
      it->second->SetBinContent(i, fSelectionHists[it->first]->GetBinContent(i));
    }
  }
}
