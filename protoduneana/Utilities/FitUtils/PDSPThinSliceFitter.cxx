// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"

// ROOT includes
#include "TFile.h"
#include "THStack.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"


#include "PDSPThinSliceFitter.h"

protoana::PDSPThinSliceFitter::PDSPThinSliceFitter(std::string fcl_file,
                                                   std::string output_file)
    : fOutputFile(output_file.c_str(), "RECREATE") {
  Configure(fcl_file);
}

protoana::PDSPThinSliceFitter::~PDSPThinSliceFitter() {
  fOutputFile.Close();
}

void protoana::PDSPThinSliceFitter::InitializeMCSamples() {
  for (size_t i = 0; i < fSampleSets.size(); ++i) {
    fhicl::ParameterSet sample_set = fSampleSets[i];
    std::string sample_name = sample_set.get<std::string>("Name");
    int sample_ID = sample_set.get<int>("ID");

    if (fSamples.find(sample_set.get<int>("ID")) == fSamples.end()) {
      fSamples[sample_ID] = std::vector<ThinSliceSample>();
    }
    fIsSignalSample[sample_ID] = sample_set.get<bool>("IsSignal");

    if (sample_set.get<bool>("IsSignal")) {
      std::vector<double> bins = sample_set.get<std::vector<double>>("SignalBins");

      for (size_t j = 1; j < bins.size(); ++j) {
        ThinSliceSample sample(sample_name, fSelectionIDs, fIncidentRecoBins,
                               fSelectedRecoBins, true,
                               {bins[j-1], bins[j]});
        fSamples[sample_ID].push_back(sample);
      }
    }
    else {
      ThinSliceSample sample(sample_name, fSelectionIDs, fIncidentRecoBins,
                             fSelectedRecoBins);
      fSamples[sample_ID].push_back(sample);
    }
  }
}

void protoana::PDSPThinSliceFitter::BuildMCSamples() {
  //Open the MC file and set branches
  TFile fMCFile(fMCFileName.c_str(), "OPEN");
  fMCTree = (TTree*)fMCFile.Get(fTreeName.c_str());

  int true_sample, selection_ID; 
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  //std::vector<double> * true_beam_incidentEnergies = 0x0
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  fMCTree->SetBranchAddress("new_interaction_topology", &true_sample);
  fMCTree->SetBranchAddress("selection_ID", &selection_ID);
  fMCTree->SetBranchAddress("true_beam_interactingEnergy",
                            &true_beam_interactingEnergy);
  fMCTree->SetBranchAddress("reco_beam_interactingEnergy",
                            &reco_beam_interactingEnergy);
  //fMCTree->SetBranchAddress("true_beam_incidentEnergies",
  //                          &true_beam_incidentEnergies);
  fMCTree->SetBranchAddress("reco_beam_incidentEnergies",
                            &reco_beam_incidentEnergies);

  for (int i = 0; i < fMCTree->GetEntries(); ++i) {
    fMCTree->GetEntry(i);

    if (fSamples.find(true_sample) == fSamples.end())
      continue;

    std::vector<ThinSliceSample> & samples = fSamples[true_sample];
    bool is_signal = fIsSignalSample[true_sample];

    //Get the selection

    if (!is_signal) {
      ThinSliceSample & sample = samples.at(0);
      if (reco_beam_incidentEnergies->size()) {
        sample.FillIncidentHist(*reco_beam_incidentEnergies);
        sample.FillSelectionHist(selection_ID, reco_beam_interactingEnergy);
      }
    }
    else {
      //Iterate through the true bins and find the correct one
      for (size_t j = 0; j < samples.size(); ++j) {
        ThinSliceSample & sample = samples.at(j);
        if (sample.CheckInSignalRange(true_beam_interactingEnergy)) {
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

void protoana::PDSPThinSliceFitter::SaveMCSamples() {
  fOutputFile.cd();
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      it->second.at(i).GetIncidentHist().Write();
      const std::map<int, TH1D> & hists = it->second.at(i).GetSelectionHists();
      for (auto it2 = hists.begin(); it2 != hists.end(); ++it2) {
        it2->second.Write();
      }
    }
  }
}

void protoana::PDSPThinSliceFitter::BuildAndSaveStacks() {

  //Incident
  THStack sIncident("IncidentStack", "");
  int iColor = 1;
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      TH1D & inc_hist = it->second.at(i).GetIncidentHist();
      inc_hist.SetFillColor(iColor);
      inc_hist.SetLineColor(iColor);
      sIncident.Add(&inc_hist);
      ++iColor;
    }
  }

  fOutputFile.cd();
  sIncident.Write();

  for (auto it1 = fSelectionIDs.begin(); it1 != fSelectionIDs.end(); ++it1) {
    std::string s_name = it1->second + "Stack";
    THStack sSelection(s_name.c_str(), "");
    iColor = 1;
    for (auto it2 = fSamples.begin(); it2 != fSamples.end(); ++it2) {
      for (size_t i = 0; i < it2->second.size(); ++i) {
        TH1D & sel_hist = it2->second.at(i).GetSelectionHist(it1->first);
        sel_hist.SetFillColor(iColor);
        sel_hist.SetLineColor(iColor);
        sSelection.Add(&sel_hist);
        ++iColor;
      }
    }
    sSelection.Write();
  }
}

void protoana::PDSPThinSliceFitter::Configure(std::string fcl_file) {
  fhicl::ParameterSet pset;

  char const* fhicl_env = getenv("FHICL_FILE_PATH");
  std::string search_path;

  if (!fhicl_env) {
    std::cerr << "Expected environment variable FHICL_FILE_PATH is missing " <<
                 "or empty: using \".\"\n";
    search_path = ".";
  }
  else {
    search_path = std::string{fhicl_env};
  }

  cet::filepath_first_absolute_or_lookup_with_dot lookupPolicy{search_path};
  fhicl::make_ParameterSet(fcl_file, lookupPolicy, pset);

  //Getting the configurable parameters 
  fMCFileName = pset.get<std::string>("MCFileName");
  fDataFileName = pset.get<std::string>("DataFileName");
  fTreeName = pset.get<std::string>("TreeName");

  std::vector<std::pair<int, std::string>> temp_vec =
      pset.get<std::vector<std::pair<int, std::string>>>("SelectionIDs");
  fSelectionIDs = std::map<int, std::string>(temp_vec.begin(), temp_vec.end());

  fSelectedRecoBins = pset.get<std::vector<double>>("SelectedRecoBins");
  fIncidentRecoBins = pset.get<std::vector<double>>("IncidentRecoBins");
  
  fSampleSets = pset.get<std::vector<fhicl::ParameterSet>>("Samples");
}
