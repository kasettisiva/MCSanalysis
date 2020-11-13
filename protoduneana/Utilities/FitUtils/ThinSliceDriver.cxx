#include "ThinSliceDriver.h"

#include "THStack.h"
#include "TCanvas.h"

protoana::ThinSliceDriver::ThinSliceDriver(
    const fhicl::ParameterSet & extra_options)
    : fExtraOptions(extra_options) {}

protoana::ThinSliceDriver::~ThinSliceDriver() {}

void protoana::ThinSliceDriver::CompareDataMC(
    std::map<int, std::vector<ThinSliceSample>> & samples,
    ThinSliceDataSet & data_set, TFile & output_file,
    std::vector<std::pair<int, int>> plot_style,
    bool plot_rebinned,
    bool post_fit) {
  //Build the incident stack and compare
  THStack incident_stack(
      (post_fit ? "PostFitIncidentStack" : "NominalIncidentStack"), "");
  size_t iColor = 0;
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      it->second.at(i).RefillRebinnedHists();

      TH1D & inc_hist = (plot_rebinned ?
                         it->second.at(i).GetRebinnedIncidentHist() :
                         it->second.at(i).GetIncidentHist());

      std::pair<int, int> color_fill = GetColorAndStyle(iColor, plot_style);
      inc_hist.SetFillColor(color_fill.first);
      inc_hist.SetFillStyle(color_fill.second);
      inc_hist.SetLineColor(kBlack);
      incident_stack.Add(&inc_hist);
      ++iColor;
    }
  }
  output_file.cd();
  incident_stack.Write();

  TH1D & inc_data_hist = (plot_rebinned ?
                          data_set.GetRebinnedIncidentHist() :
                          data_set.GetIncidentHist());
  inc_data_hist.SetLineColor(kBlack);
  inc_data_hist.SetMarkerColor(kBlack);
  inc_data_hist.SetMarkerStyle(20);

  TCanvas cIncident((post_fit ? "cPostFitIncident" : "cNominalIncident"), "");
  cIncident.SetTicks();
  //THStack * incident_stack = (post_fit ? fPostFitIncidentMCStack :
  //                                       fNominalIncidentMCStack);
  incident_stack.Draw("hist");
  incident_stack.SetTitle("Incident Sample;Reconstructed KE (MeV)");
  incident_stack.GetHistogram()->SetTitleSize(.04, "X");
  incident_stack.Draw("hist");
  inc_data_hist.Draw("e1 same");
  output_file.cd();
  cIncident.Write();

  //Get the full incident hist from stack
  TList * l = (TList*)incident_stack.GetHists();
  TH1D * hIncidentMC = (TH1D*)l->At(0)->Clone();
  for (int i = 1; i < l->GetSize(); ++i) {
    hIncidentMC->Add((TH1D*)l->At(i));
  }

  std::string inc_ratio_name = "IncidentRatio";
  inc_ratio_name += (post_fit ? "PostFit" : "Nominal");
  TH1D * hIncidentRatio
      = (TH1D*)inc_data_hist.Clone(inc_ratio_name.c_str());
  hIncidentRatio->Divide(hIncidentMC);
  hIncidentRatio->Write(); 

  CompareSelections(
      samples, data_set, output_file, plot_style, plot_rebinned,
      post_fit);
}

/*
void protoana::ThinSliceDriver::BuildAndSaveStacks(
    std::map<int, std::vector<ThinSliceSample>> & samples,
    ThinSliceDataSet & data_set,
    TFile & output_file,
    THStack * incident_stack, 
    std::map<int, THStack *> * selection_stacks_x,
    std::map<int, THStack *> * selection_stacks_y,
    std::map<int, THStack *> * selection_stacks_z,
    bool post_fit) {
   
}*/

std::pair<int, int> protoana::ThinSliceDriver::GetColorAndStyle(
    size_t i,
    const std::vector<std::pair<int, int>> & plot_style) {
  return {plot_style.at(i % plot_style.size()).first,
          (i < plot_style.size() ? 1001: 3244)};
}
