#include "ThinSliceDataSet.h"
#include "ThinSliceSample.h"

protoana::ThinSliceDataSet::ThinSliceDataSet(
    const std::vector<double> & incident_bins,
    const std::vector<fhicl::ParameterSet> & selections) {
  fIncidentHist = TH1D("Data_incident_hist",
                           "Data;Reconstructed KE (MeV)",
                           incident_bins.size() - 1,
                           &incident_bins[0]);
  for (auto it = selections.begin(); it != selections.end(); ++it) {
    fSelectionNames[it->get<int>("ID")] = it->get<std::string>("Name");
    std::string sel_name = "Data_selected_" + it->get<std::string>("Name") +
                           "_hist";
    std::vector<std::vector<double>> selected_bins =
        it->get<std::vector<std::vector<double>>>("RecoBins");

    std::vector<std::string> titles =
        it->get<std::vector<std::string>>("AxisTitles");
    TString title = "Data";
    for (auto & t : titles) {
      title += ";" + t; 
    }

    if (selected_bins.size() == 1) {
      fSelectionHists[it->get<int>("ID")] = new TH1D(
          sel_name.c_str(), title/*.c_str()"Data;Reconstructed KE (MeV)"*/,
          selected_bins[0].size() - 1, &selected_bins[0][0]);
    }
    else if (selected_bins.size() == 2) {
      fSelectionHists[it->get<int>("ID")] = new TH2D(
          sel_name.c_str(), title/*.c_str()"Data;Reconstructed KE (MeV)"*/,
          selected_bins[0].size() - 1, &selected_bins[0][0],
          selected_bins[1].size() - 1, &selected_bins[1][0]);
    }
    else if (selected_bins.size() == 3) {
      fSelectionHists[it->get<int>("ID")] = new TH3D(
          sel_name.c_str(), title/*.c_str()"Data;Reconstructed KE (MeV)"*/,
          selected_bins[0].size() - 1, &selected_bins[0][0],
          selected_bins[1].size() - 1, &selected_bins[1][0],
          selected_bins[2].size() - 1, &selected_bins[2][0]);
    }
    /*else {
     * throw
     * }*/
  }
}

void protoana::ThinSliceDataSet::MakeRebinnedHists() {
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
        TString title = sel_hist->GetTitle();
        title += ";";
        title += sel_hist->GetXaxis()->GetTitle();
        fSelectionHistsRebinned[it->first] = new TH1D(
            name.c_str(), title/*.c_str()sel_hist->GetTitle()*/,
            sel_hist->GetNbinsX(), 0, sel_hist->GetNbinsX());
        Rebin1D(sel_hist, fSelectionHistsRebinned[it->first]);
      }
      else if (nAxes == 2) {
        std::string title = sel_hist->GetTitle();
        title += ";";
        title += sel_hist->GetXaxis()->GetTitle();
        title += ";";
        title += sel_hist->GetYaxis()->GetTitle();

        fSelectionHistsRebinned[it->first] = new TH2D(
            name.c_str(), title.c_str()/*sel_hist->GetTitle()*/,
            sel_hist->GetNbinsX(), 0, sel_hist->GetNbinsX(),
            sel_hist->GetNbinsY(), 0, sel_hist->GetNbinsY());
        Rebin2D(sel_hist, fSelectionHistsRebinned[it->first]);
      }
      else if (nAxes == 3) {
        std::string title = sel_hist->GetTitle();
        title += ";";
        title += sel_hist->GetXaxis()->GetTitle();
        title += ";";
        title += sel_hist->GetYaxis()->GetTitle();
        title += ";";
        title += sel_hist->GetZaxis()->GetTitle();

        fSelectionHistsRebinned[it->first] = new TH3D(
            name.c_str(), title.c_str()/*sel_hist->GetTitle()*/,
            sel_hist->GetNbinsX(), 0, sel_hist->GetNbinsX(),
            sel_hist->GetNbinsY(), 0, sel_hist->GetNbinsY(),
            sel_hist->GetNbinsZ(), 0, sel_hist->GetNbinsZ());
        Rebin3D(sel_hist, fSelectionHistsRebinned[it->first]);
      }
    }

    fMadeRebinned = true;
  }
}

void protoana::ThinSliceDataSet::Rebin1D(TH1 * sel_hist, TH1 * rebinned) {
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

void protoana::ThinSliceDataSet::Rebin2D(TH1 * sel_hist, TH1 * rebinned) {
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

void protoana::ThinSliceDataSet::Rebin3D(TH1 * sel_hist, TH1 * rebinned) {
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
