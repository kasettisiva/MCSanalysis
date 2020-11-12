#include "AbsCexDriver.h"

#include "ThinSliceDriverFactory.h"
DECLARE_THINSLICEDRIVER_FACTORY_NS(protoana::AbsCexDriver, protoana, AbsCexDriver)

protoana::AbsCexDriver::AbsCexDriver(const std::string & analysis)
    : ThinSliceDriver(analysis) {}

void protoana::AbsCexDriver::BuildMCSamples(
    TTree * tree,
    std::map<int, std::vector<ThinSliceSample>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::map<int, double> & nominal_fluxes,
    std::map<int, std::vector<double>> & fluxes_by_sample) {
  int sample_ID, selection_ID; 
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  tree->SetBranchAddress("new_interaction_topology", &sample_ID);
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("true_beam_interactingEnergy",
                            &true_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                            &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                            &reco_beam_incidentEnergies);

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    if (samples.find(sample_ID) == samples.end())
      continue;

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
        double energy[1] = {reco_beam_interactingEnergy};
        sample.FillSelectionHist(selection_ID, energy);
      }
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
            sample.FillSelectionHist(selection_ID, reco_beam_interactingEnergy);
          }
          break;
        }
      }
    }
  }
}

void protoana::AbsCexDriver::BuildDataHists(
    TTree * tree, TH1D & incident_hist,
    std::map<int, TH1 *> & selected_hists) {
  int selection_ID; 
  double reco_beam_interactingEnergy;
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                        &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                        &reco_beam_incidentEnergies);

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        selected_hists[selection_ID]->Fill(reco_beam_interactingEnergy);
      }
    }
  }
}

/*
std::pair<double, size_t> protoana::AbsCexDriver::CalculateChi2() {
  //Go through the selection hists and fit mc to data
  double chi2 = 0.;
  size_t nPoints = 0;

  for (auto it = fSelectedDataHists.begin();
       it != fSelectedDataHists.end(); ++it) {
    TH1D & data_hist = it->second;
    int selection_ID = it->first;
    for (int i = 1; i <= data_hist.GetNbinsX(); ++i) {
      double data_val = data_hist.GetBinContent(i);
      double data_err = data_hist.GetBinError(i);

      double mc_val = 0.;
      //Go through all the samples and get the values from mc
      for (auto it2 = fSamples.begin(); it2 != fSamples.end(); ++it2) {
        std::vector<ThinSliceSample> & samples = it2->second;
        for (size_t j = 0; j < samples.size(); ++j) {
          ThinSliceSample & sample = samples[j];
          mc_val += sample.GetSelectionHist(selection_ID).GetBinContent(i);
        }
      }
      chi2 += (std::pow((data_val - mc_val), 2) / std::pow(data_err, 2));
      ++nPoints;
    }
  }

  //Go through the incident hist. Calculate it as its own point.
  //Then add reduced incident chi2 to full chi2 
  double incident_chi2 = 0.;
  size_t incident_points = 0;
  for (int i = 1; i <= fIncidentDataHist.GetNbinsX(); ++i) {
    double data_val = fIncidentDataHist.GetBinContent(i);
    double data_err = fIncidentDataHist.GetBinError(i);

    double mc_val = 0.;
    //Go through all the samples and get the values from mc
    for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
      std::vector<ThinSliceSample> & samples = it->second;
      for (size_t j = 0; j < samples.size(); ++j) {
        ThinSliceSample & sample = samples[j];
        mc_val += sample.GetIncidentHist().GetBinContent(i);
      }
    }

    incident_chi2 += (std::pow((data_val - mc_val), 2) / std::pow(data_err, 2));
    ++incident_points;
  }

  if (fReducedIncidentChi2) {
    //Add reduced incident chi2 to full chi2, increment total nPoints
    chi2 += incident_chi2 / incident_points;
    ++nPoints;
  }
  else {
    chi2 += incident_chi2;
    nPoints += incident_points;
  }

  return {chi2, nPoints};

}*/

protoana::AbsCexDriver::~AbsCexDriver() {}
