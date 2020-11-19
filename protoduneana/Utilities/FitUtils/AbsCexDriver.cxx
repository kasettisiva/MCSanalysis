#include "AbsCexDriver.h"
#include "TCanvas.h"
#include "THStack.h"

#include "ThinSliceDriverFactory.h"
DECLARE_THINSLICEDRIVER_FACTORY_NS(protoana::AbsCexDriver, protoana, AbsCexDriver)

protoana::AbsCexDriver::~AbsCexDriver() {}

protoana::AbsCexDriver::AbsCexDriver(
    const fhicl::ParameterSet & extra_options)
    : ThinSliceDriver(extra_options) {}

void protoana::AbsCexDriver::BuildMCSamples(
    TTree * tree,
    std::map<int, std::vector<ThinSliceSample>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::map<int, double> & nominal_fluxes,
    std::map<int, std::vector<double>> & fluxes_by_sample/*,
    std::map<int, std::pair<TH1D, TH1D>> & signal_eff_parts*/) {

  int sample_ID, selection_ID; 
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  double reco_beam_endZ;
  std::vector<double> * reco_beam_incidentEnergies = 0x0,
                      * true_beam_incidentEnergies = 0x0;
  std::vector< int > * true_beam_slices = 0x0;
  tree->SetBranchAddress("new_interaction_topology", &sample_ID);
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("true_beam_interactingEnergy",
                         &true_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                         &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                         &reco_beam_incidentEnergies);
  tree->SetBranchAddress("true_beam_incidentEnergies",
                         &true_beam_incidentEnergies);
  tree->SetBranchAddress("true_beam_slices",
                         &true_beam_slices);

  //Determine the slice cut
  //based on the endZ cut and wire pitch 
  double pitch = fExtraOptions.get<double>("WirePitch");
  double z0 = fExtraOptions.get<double>("Z0");
  double endZ_cut = fExtraOptions.get<double>("EndZCut");
  int slice_cut = std::floor((endZ_cut - (z0 - pitch/2.))/pitch);

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    if (samples.find(sample_ID) == samples.end()) {
      std::cout << "Warning: skipping sample " << sample_ID << std::endl;
      continue;
    }

    //Build the true incident energy vector based on the slice cut
    std::vector<double> good_true_incEnergies;
    for (size_t j = 0; j < true_beam_incidentEnergies->size(); ++j) {
      int slice = (*true_beam_slices)[j]; 
      if (slice > slice_cut) continue;
      good_true_incEnergies.push_back((*true_beam_incidentEnergies)[j]);
    }

    std::vector<ThinSliceSample> & samples_vec = samples.at(sample_ID);
    bool is_signal = signal_sample_checks.at(sample_ID);

    if (!is_signal) {
      ThinSliceSample & sample = samples_vec.at(0);
      int flux_type = sample.GetFluxType();
      nominal_fluxes[flux_type] += 1.;
      fluxes_by_sample[sample_ID][0] += 1.;
      sample.AddFlux();
      if (reco_beam_incidentEnergies->size()) {
        sample.FillIncidentHist(*reco_beam_incidentEnergies);
        if (selection_ID != 4) {
          double energy[1] = {reco_beam_interactingEnergy};
          sample.FillSelectionHist(selection_ID, energy);
        }
        else {
          double endZ[1] = {reco_beam_endZ};
          sample.FillSelectionHist(selection_ID, endZ);
        }
      }

      //Fill the total incident hist with truth info
      //sample.FillTrueIncidentHist(*true_beam_incidentEnergies);
      sample.FillTrueIncidentHist(good_true_incEnergies);
    }
    else {
      //Iterate through the true bins and find the correct one
      for (size_t j = 0; j < samples_vec.size(); ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(true_beam_interactingEnergy)) {
          int flux_type = sample.GetFluxType();
          nominal_fluxes[flux_type] += 1.;
          fluxes_by_sample[sample_ID][j] += 1.;
          sample.AddFlux();
          if (reco_beam_incidentEnergies->size()) {
            sample.FillIncidentHist(*reco_beam_incidentEnergies);
            if (selection_ID != 4) {
              double energy[1] = {reco_beam_interactingEnergy};
              sample.FillSelectionHist(selection_ID, energy);
            }
            else {
              double endZ[1] = {reco_beam_endZ};
              sample.FillSelectionHist(selection_ID, endZ);
            }
          }

          //Fill the total incident hist with truth info
         //sample.FillTrueIncidentHist(*true_beam_incidentEnergies);
          sample.FillTrueIncidentHist(good_true_incEnergies);

          break;
        }
      }
      
      //Fill the total hist with truth info
      //signal_eff_parts[sample_ID].first.Fill(true_beam_interactingEnergy);
    }

  }

  /*
  //Determine the efficiencies for the signals
  for (auto it = signal_eff_parts.begin(); it != signal_eff_parts.end(); ++it) {
    int sample_ID = it->first;
    //TH1D & total = it->second.first;
    TH1D & passed = it->second.second;

    //Fill the passed with the hist contents from above
    std::vector<ThinSliceSample> & samples_vec = samples.at(sample_ID);
    for (size_t i = 0; i < samples_vec.size(); ++i) {
      //Need to worry bout over/underflow?
      passed.SetBinContent(i+1, samples_vec[i].GetNominalFlux());
    }
  }*/
}

void protoana::AbsCexDriver::BuildDataHists(
    TTree * tree, ThinSliceDataSet & data_set, double & flux) {
  int selection_ID; 
  double reco_beam_interactingEnergy, reco_beam_endZ;
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                        &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                        &reco_beam_incidentEnergies);

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  flux = tree->GetEntries();

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4) {
          selected_hists[selection_ID]->Fill(reco_beam_interactingEnergy);
        }
        else {
          selected_hists[selection_ID]->Fill(reco_beam_endZ);
        }
      }
    }
  }
}

void protoana::AbsCexDriver::BuildFakeData(
    TTree * tree, std::map<int, std::vector<ThinSliceSample>> & samples,
    ThinSliceDataSet & data_set, double & flux) {
  std::string routine = fExtraOptions.get<std::string>("FakeDataRoutine");
  if (routine == "SampleScales") {
    FakeDataSampleScales(tree, samples, data_set, flux);
  }
}

void protoana::AbsCexDriver::FakeDataSampleScales(
    TTree * tree, std::map<int, std::vector<ThinSliceSample>> & samples,
    ThinSliceDataSet & data_set, double & flux) {

  //Build the map for fake data scales
  std::vector<std::pair<int, double>> temp_vec
      = fExtraOptions.get<std::vector<std::pair<int, double>>>(
          "FakeDataScales");
  std::map<int, double> fake_data_scales(temp_vec.begin(), temp_vec.end()); 

  int sample_ID, selection_ID; 
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  double reco_beam_endZ;
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  tree->SetBranchAddress("new_interaction_topology", &sample_ID);
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("true_beam_interactingEnergy",
                         &true_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                         &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                         &reco_beam_incidentEnergies);

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  flux = 0.;

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    if (samples.find(sample_ID) == samples.end())
      continue;

    double scale = (fake_data_scales.find(sample_ID) != fake_data_scales.end() ?
                    fake_data_scales.at(sample_ID) : 1.);
    flux += scale; //1 or scaled
    if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j], scale);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4) {
          selected_hists[selection_ID]->Fill(reco_beam_interactingEnergy, scale);
        }
        else {
          selected_hists[selection_ID]->Fill(reco_beam_endZ, scale);
        }
      }
    }
  }
}

std::pair<double, size_t> protoana::AbsCexDriver::CalculateChi2(
    std::map<int, std::vector<ThinSliceSample>> & samples,
    ThinSliceDataSet & data_set) {

  double chi2 = 0.;
  size_t nPoints = 0;

  std::map<int, TH1 *> & selected_data_hists = data_set.GetSelectionHists();
  for (auto it = selected_data_hists.begin();
       it != selected_data_hists.end(); ++it) {
    TH1D * data_hist = (TH1D*)it->second;
    int selection_ID = it->first;
    for (int i = 1; i <= data_hist->GetNbinsX(); ++i) {
      double data_val = data_hist->GetBinContent(i);
      double data_err = data_hist->GetBinError(i);

      double mc_val = 0.;
      //Go through all the samples and get the values from mc
      for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
        std::vector<ThinSliceSample> & samples_vec = it2->second;
        for (size_t j = 0; j < samples_vec.size(); ++j) {
          ThinSliceSample & sample = samples_vec[j];
          mc_val += sample.GetSelectionHist(selection_ID)->GetBinContent(i);
        }
      }
      chi2 += (std::pow((data_val - mc_val), 2) / std::pow(data_err, 2));
      ++nPoints;
    }
  }

  //Go through the incident hist. Calculate it as its own point.
  //Then add reduced incident chi2 to full chi2 
  if (!fExtraOptions.get<bool>("SkipIncidentChi2")) {
    double incident_chi2 = 0.;
    size_t incident_points = 0;
    TH1D & incident_data_hist = data_set.GetIncidentHist();
    for (int i = 1; i <= incident_data_hist.GetNbinsX(); ++i) {
      double data_val = incident_data_hist.GetBinContent(i);
      double data_err = incident_data_hist.GetBinError(i);

      double mc_val = 0.;
      //Go through all the samples and get the values from mc
      for (auto it = samples.begin(); it != samples.end(); ++it) {
        std::vector<ThinSliceSample> & samples_vec = it->second;
        for (size_t j = 0; j < samples_vec.size(); ++j) {
          ThinSliceSample & sample = samples_vec[j];
          mc_val += sample.GetIncidentHist().GetBinContent(i);
        }
      }

      incident_chi2 += (std::pow((data_val - mc_val), 2) / std::pow(data_err, 2));
      ++incident_points;
    }

    if (fExtraOptions.get<bool>("ReducedChi2")) {
      //Add reduced incident chi2 to full chi2, increment total nPoints
      chi2 += incident_chi2 / incident_points;
      ++nPoints;
    }
    else {
      chi2 += incident_chi2;
      nPoints += incident_points;
    }
  }

  return {chi2, nPoints};
}

void protoana::AbsCexDriver::CompareSelections(
    std::map<int, std::vector<ThinSliceSample>> & samples,
    ThinSliceDataSet & data_set, TFile & output_file,
    std::vector<std::pair<int, int>> plot_style,
    bool plot_rebinned,
    bool post_fit) {

  output_file.cd();
  std::map<int, TH1*> data_hists
      = (plot_rebinned ?
         data_set.GetRebinnedSelectionHists() :
         data_set.GetSelectionHists());
  for (auto it = data_hists.begin(); it != data_hists.end(); ++it) {
    int selection_ID = it->first;
    TH1D * data_hist = (TH1D*)it->second;
    data_hist->SetLineColor(kBlack);
    data_hist->SetMarkerColor(kBlack);
    data_hist->SetMarkerStyle(20);   

    std::string canvas_name = "c";
    canvas_name += (post_fit ? "PostFit" : "Nominal") +
                   data_set.GetSelectionName(selection_ID);
    TCanvas cSelection(canvas_name.c_str(), "");
    cSelection.SetTicks();

    //Build the mc stack
    std::string stack_name = (post_fit ? "PostFit" : "Nominal") +
                             data_set.GetSelectionName(selection_ID) +
                             "Stack";
    THStack mc_stack(stack_name.c_str(), "");
    size_t iColor = 0;
    for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
      for (size_t i = 0; i < it2->second.size(); ++i) {
        TH1D * sel_hist
            = (TH1D*)(plot_rebinned ?
                      it2->second.at(i).GetRebinnedSelectionHist(selection_ID) :
                      it2->second.at(i).GetSelectionHist(selection_ID));

        std::pair<int, int> color_fill = GetColorAndStyle(iColor, plot_style);
        sel_hist->SetFillColor(color_fill.first);
        sel_hist->SetFillStyle(color_fill.second);
        sel_hist->SetLineColor(kBlack);
        mc_stack.Add(sel_hist);
        ++iColor;
      }
    }
    mc_stack.Write();


    mc_stack.Draw("hist");
    double max_mc = mc_stack.GetHistogram()->GetMaximum();
    int max_data_bin = data_hist->GetMaximumBin();
    double max_data = data_hist->GetBinContent(max_data_bin) +
                      data_hist->GetBinError(max_data_bin);
    std::string title = data_set.GetSelectionName(selection_ID) +
                        (selection_ID != 4 ? ";Reconstructed KE (MeV)" :
                                             ";Reconstructed End Z (cm)");
    mc_stack.SetTitle(title.c_str());
    mc_stack.GetHistogram()->SetTitleSize(.04, "X");
    mc_stack.SetMaximum((max_data > max_mc ? max_data : max_mc));
    mc_stack.Draw("hist");
    data_hist->Draw("e1 same");

    cSelection.Write();

    //Get the full incident hist from stack
    TList * l = (TList*)mc_stack.GetHists();
    TH1D * hMC = (TH1D*)l->At(0)->Clone();
    for (int i = 1; i < l->GetSize(); ++i) {
      hMC->Add((TH1D*)l->At(i));
    }

    std::string ratio_name = data_set.GetSelectionName(selection_ID) + "Ratio" +
                             (post_fit ? "PostFit" : "Nominal");
    TH1D * hRatio
        = (TH1D*)data_hist->Clone(ratio_name.c_str());
    hRatio->Divide(hMC);
    hRatio->Write(); 
  }
}
