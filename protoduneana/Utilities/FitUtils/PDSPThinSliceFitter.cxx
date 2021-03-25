// C++ language includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include "math.h"
#include "stdio.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"

#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TList.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TDecompChol.h"
#include "TF1.h"

#include "PDSPThinSliceFitter.h"
#include "AbsCexDriver.h"

#include "ThinSliceDriverRegistry.h"

protoana::PDSPThinSliceFitter::PDSPThinSliceFitter(std::string fcl_file,
                                                   std::string output_file)
    : fOutputFile(output_file.c_str(), "RECREATE") {
  Configure(fcl_file);
  gStyle->SetOptStat(0);
  gROOT->SetBatch(1);

  try {
    fThinSliceDriver = protoana::ThinSliceDriverRegistry::Instance()->GetDriver(
        fDriverName, fAnalysisOptions);
  }
  catch (const std::runtime_error & e) {
    protoana::ThinSliceDriverRegistry::Instance()->PrintAvailableDrivers();
    throw e;
  }

}

protoana::PDSPThinSliceFitter::~PDSPThinSliceFitter() {
  fOutputFile.Close();
}

//Good
void protoana::PDSPThinSliceFitter::MakeMinimizer() {
  fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
    (ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  fMinimizer->SetMaxFunctionCalls(fMaxCalls);
  fMinimizer->SetTolerance(fTolerance);

  size_t total_parameters = fTotalSignalParameters + fTotalFluxParameters + fTotalSystParameters;
  TH1D parsHist("preFitPars", "", total_parameters, 0, total_parameters);

  size_t n_par = 0;
  for (auto it = fSignalParameters.begin();
       it != fSignalParameters.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (fRandomStart) {
        it->second[i] = fRNG.Gaus(1., .1);
      }
      fMinimizer->SetVariable(n_par, fSignalParameterNames[it->first][i].c_str(),
                              it->second[i], 0.01);
      fMinimizer->SetVariableLimits(n_par, fLowerLimit, fUpperLimit);
      parsHist.SetBinContent(n_par+1, it->second[i]);
      parsHist.GetXaxis()->SetBinLabel(
          n_par+1, fSignalParameterNames[it->first][i].c_str());
      ++n_par;

    }
  }

  for (auto it = fFluxParameters.begin();
       it != fFluxParameters.end(); ++it) {
    if (fRandomStart) {
      it->second = fRNG.Gaus(1., .1);
    }
    fMinimizer->SetVariable(n_par, fFluxParameterNames[it->first].c_str(),
                            it->second, 0.01);
    fMinimizer->SetVariableLimits(n_par, fLowerLimit, fUpperLimit);
    parsHist.SetBinContent(n_par+1, it->second);
    parsHist.GetXaxis()->SetBinLabel(
        n_par+1, fFluxParameterNames[it->first].c_str());
    ++n_par;
  }

  for (auto it = fSystParameters.begin(); it != fSystParameters.end(); ++it) {
    std::cout << "Adding parameter " << it->second.GetName().c_str() << std::endl;
    if (fRandomStart) {
      it->second.SetValue(fRNG.Gaus(1., .1)*it->second.GetCentral());
    }
    fMinimizer->SetVariable(n_par, it->second.GetName().c_str(),
                            it->second.GetValue(), 0.01);
    fMinimizer->SetVariableLimits(n_par, it->second.GetLowerLimit(),
                                  it->second.GetUpperLimit());
    parsHist.SetBinContent(n_par+1, it->second.GetValue()/it->second.GetCentral());
    parsHist.GetXaxis()->SetBinLabel(
        n_par+1, it->second.GetName().c_str());
    ++n_par;
  }

  parsHist.SetMarkerColor(kBlue);
  parsHist.SetMarkerStyle(20);
  fOutputFile.cd();
  parsHist.Write();

}

//Good
void protoana::PDSPThinSliceFitter::InitializeMCSamples() {
  for (size_t i = 0; i < fSampleSets.size(); ++i) {
    fhicl::ParameterSet sample_set = fSampleSets[i];
    std::string sample_name = sample_set.get<std::string>("Name");
    int sample_ID = sample_set.get<int>("ID");

    if (fSamples.find(sample_ID) == fSamples.end()) {
      fSamples[sample_ID] = std::vector<std::vector<ThinSliceSample>>();
      fFluxesBySample[sample_ID] = std::vector<std::vector<double>>();
    }
    fIsSignalSample[sample_ID] = sample_set.get<bool>("IsSignal");

    int flux_type = sample_set.get<int>("FluxType");
    std::vector<double> bins
        = sample_set.get<std::vector<double>>("SignalBins");
    fSignalBins[sample_ID] = bins;

    for (size_t j = 1; j < fBeamEnergyBins.size(); ++j) {
      fSamples[sample_ID].push_back(std::vector<ThinSliceSample>());
      fFluxesBySample[sample_ID].push_back(std::vector<double>());
      if (sample_set.get<bool>("IsSignal")) {

        if (j == 1) {
          fSignalParameters[sample_ID] = std::vector<double>();
          fSignalParameterNames[sample_ID] = std::vector<std::string>();
        }

        ThinSliceSample underflow_sample(
            sample_name + "Underflow",
            flux_type, fSelectionSets,
            fIncidentRecoBins, fTrueIncidentBins, j);
        fSamples[sample_ID].back().push_back(underflow_sample);
        fFluxesBySample[sample_ID].back().push_back(0.);

        for (size_t k = 1; k < bins.size(); ++k) {
          ThinSliceSample sample(sample_name, flux_type, fSelectionSets,
                                 fIncidentRecoBins, fTrueIncidentBins,
                                 j, true,
                                 {bins[k-1], bins[k]});
          fSamples[sample_ID].back().push_back(sample);

          if (j == 1) {
            fSignalParameters[sample_ID].push_back(1.);

            std::string par_name = "par_" + sample_name + "_" +
                                   PreciseToString(bins[k-1]) + "_" +
                                   PreciseToString(bins[k]);
            fSignalParameterNames[sample_ID].push_back(par_name);
            ++fTotalSignalParameters;
          }
          fFluxesBySample[sample_ID].back().push_back(0.);
        }

        ThinSliceSample overflow_sample(sample_name + "Overflow",
                               flux_type, fSelectionSets,
                               fIncidentRecoBins, fTrueIncidentBins, j);
        fSamples[sample_ID].back().push_back(overflow_sample);
        fFluxesBySample[sample_ID].back().push_back(0.);
      }
      else {
        ThinSliceSample sample(sample_name, flux_type, fSelectionSets,
                                       fIncidentRecoBins, fTrueIncidentBins, j);
        fSamples[sample_ID].back().push_back(sample);
        fFluxesBySample[sample_ID].back().push_back(0.);
        if (fFluxParameters.find(flux_type) != fFluxParameters.end() &&
            j == 1) {
          fFluxParsToSamples[flux_type].push_back(sample_ID);
        }
      }
    }
  }

  MakeMinimizer();
}

//Good
void protoana::PDSPThinSliceFitter::BuildMCSamples() {
  //Open the MC file and set branches
  TFile fMCFile(fMCFileName.c_str(), "OPEN");
  fMCTree = (TTree*)fMCFile.Get(fTreeName.c_str());

  FillMCEvents();
  fThinSliceDriver->BuildMCSamples(/*fMCTree*/fEvents, fSamples, fIsSignalSample,
                                   fNominalFluxes, fFluxesBySample,
                                   fBeamEnergyBins);
  fMCFile.Close();

  if (fDoSysts) {
    fThinSliceDriver->SetupSysts(fEvents, fSamples, fIsSignalSample,
                                 fBeamEnergyBins, fSystParameters,
                                 fOutputFile);
  }

}

void protoana::PDSPThinSliceFitter::FillMCEvents() {

  std::cout << "Filling MC Events" << std::endl;

  int sample_ID, selection_ID, event, run, subrun;
  int true_beam_PDG;
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  double true_beam_endP, true_beam_mass;
  double reco_beam_endZ, true_beam_startP;
  double beam_inst_P;
  std::vector<double> * reco_beam_incidentEnergies = 0x0,
                      * true_beam_incidentEnergies = 0x0,
                      * true_beam_traj_Z = 0x0,
                      * true_beam_traj_KE = 0x0;
  std::vector<int> * true_beam_slices = 0x0;
  fMCTree->SetBranchAddress("event", &event);
  fMCTree->SetBranchAddress("subrun", &subrun);
  fMCTree->SetBranchAddress("run", &run);

  fMCTree->SetBranchAddress("true_beam_PDG", &true_beam_PDG);

  fMCTree->SetBranchAddress("new_interaction_topology", &sample_ID);
  fMCTree->SetBranchAddress("selection_ID", &selection_ID);
  fMCTree->SetBranchAddress("true_beam_interactingEnergy",
                         &true_beam_interactingEnergy);
  fMCTree->SetBranchAddress("true_beam_endP", &true_beam_endP);
  fMCTree->SetBranchAddress("true_beam_mass", &true_beam_mass);
  fMCTree->SetBranchAddress("reco_beam_interactingEnergy",
                         &reco_beam_interactingEnergy);
  fMCTree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  fMCTree->SetBranchAddress("reco_beam_incidentEnergies",
                         &reco_beam_incidentEnergies);
  fMCTree->SetBranchAddress("true_beam_incidentEnergies",
                         &true_beam_incidentEnergies);
  fMCTree->SetBranchAddress("true_beam_slices",
                         &true_beam_slices);
  fMCTree->SetBranchAddress("true_beam_startP", &true_beam_startP);
  fMCTree->SetBranchAddress("true_beam_traj_Z", &true_beam_traj_Z);
  fMCTree->SetBranchAddress("true_beam_traj_KE", &true_beam_traj_KE);
  std::vector<double> * calibrated_dQdX = 0x0, * beam_EField = 0x0,
                      * track_pitch = 0x0;
  fMCTree->SetBranchAddress("reco_beam_calibrated_dQdX_SCE", &calibrated_dQdX);
  fMCTree->SetBranchAddress("reco_beam_EField_SCE", &beam_EField);
  fMCTree->SetBranchAddress("reco_beam_TrkPitch_SCE", &track_pitch);
  fMCTree->SetBranchAddress("beam_inst_P", &beam_inst_P);

  for (int i = 0; i < fMCTree->GetEntries(); ++i) {
    fMCTree->GetEntry(i);

    fEvents.push_back(ThinSliceEvent(event, subrun, run));
    fEvents.back().SetSampleID(sample_ID);
    fEvents.back().SetSelectionID(selection_ID);
    fEvents.back().SetTrueInteractingEnergy(true_beam_interactingEnergy);
    fEvents.back().SetRecoInteractingEnergy(reco_beam_interactingEnergy);
    fEvents.back().SetTrueEndP(true_beam_endP);
    fEvents.back().SetTrueStartP(true_beam_startP);
    fEvents.back().SetTrueMass(true_beam_mass);
    fEvents.back().SetRecoEndZ(reco_beam_endZ);

    fEvents.back().SetRecoIncidentEnergies(*reco_beam_incidentEnergies);
    fEvents.back().SetTrueIncidentEnergies(*true_beam_incidentEnergies);
    fEvents.back().SetTrueTrajZ(*true_beam_traj_Z);
    fEvents.back().SetTrueTrajKE(*true_beam_traj_KE);
    fEvents.back().SetTrueSlices(*true_beam_slices);
    fEvents.back().SetdQdXCalibrated(*calibrated_dQdX);
    fEvents.back().SetEField(*beam_EField);
    fEvents.back().SetTrackPitch(*track_pitch);
    fEvents.back().SetBeamInstP(beam_inst_P);
    fEvents.back().SetPDG(true_beam_PDG);
  }

  std::cout << "Filled MC Events" << std::endl;
}

//Good
void protoana::PDSPThinSliceFitter::ScaleMCToData() {
  double total_nominal = 0.;
  for (auto it = fNominalFluxes.begin(); it != fNominalFluxes.end(); ++it) {
    total_nominal += it->second;
  }

  fMCDataScale = fDataFlux/total_nominal;
  for (auto it = fNominalFluxes.begin(); it != fNominalFluxes.end(); ++it) {
    it->second *= fMCDataScale;
  }

  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        it->second[i][j].SetDataMCScale(fMCDataScale);
      }
    }
  }
}

void protoana::PDSPThinSliceFitter::BuildDataHists() {
  //Create the data hists
  fDataSet = ThinSliceDataSet(fIncidentRecoBins, fSelectionSets);

  TFile inputFile((!fDoFakeData ? fDataFileName.c_str() : fMCFileName.c_str()),
                  "OPEN");
  TTree * tree = (TTree*)inputFile.Get(fTreeName.c_str());
  if (!fDoFakeData) {
    fThinSliceDriver->BuildDataHists(tree, fDataSet, fDataFlux);
  }
  else {
    for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
      std::vector<ThinSliceSample> & samples_vec = it->second[0];
      fFakeDataScales[it->first] = std::vector<double>(samples_vec.size(), 0.);
    }
    fThinSliceDriver->BuildFakeData(tree, fSamples, fIsSignalSample, fDataSet,
                                    fDataFlux, fFakeDataScales);
    std::cout << "Sample scales: "; 
    for (auto it = fFakeDataScales.begin(); it != fFakeDataScales.end(); ++it) {
      std::cout << fSamples[it->first][0][0].GetName() << " ";
      for (size_t i = 0; i < it->second.size(); ++i) {
        std::cout << it->second[i] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    BuildFakeDataXSecs();
  }

  if (fDoFakeData) {
    std::cout << "Testing" << std::endl;
    for (auto it = fFakeDataXSecs.begin(); it != fFakeDataXSecs.end(); ++it) {
      std::cout << it->first << " " << it->second->GetNbinsX() << std::endl;
    }
  }

  inputFile.Close();
  if (fDoFakeData) {
    std::cout << "Testing" << std::endl;
    for (auto it = fFakeDataXSecs.begin(); it != fFakeDataXSecs.end(); ++it) {
      std::cout << it->first << " " << it->second->GetNbinsX() << std::endl;
    }
  }

  fOutputFile.cd();
  TDirectory * out = (TDirectory *)fOutputFile.mkdir("Data");
  out->cd();
  fDataSet.GetIncidentHist().Write();
  std::map<int, TH1 *> & selected_hists = fDataSet.GetSelectionHists();
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Write();
  }

  fDataSet.MakeRebinnedHists();
}


void protoana::PDSPThinSliceFitter::SaveMCSamples() {
  fOutputFile.cd();
  fOutputFile.mkdir("MC_Samples");
  fOutputFile.cd("MC_Samples");
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      auto vec = it->second.at(i);
      for (size_t j = 0; j < vec.size(); ++j) {
        //vec[j].GetIncidentHist().Write();
        const std::map<int, TH1 *> & hists = vec[j].GetSelectionHists();
        for (auto it2 = hists.begin(); it2 != hists.end(); ++it2) {
          it2->second->Write();
        }
      }
    }
  }

  TDirectory * top_dir = fOutputFile.GetDirectory("");
  TDirectory * xsec_dir = fOutputFile.mkdir("PreFitXSec");
  std::string extra_name = "PreFit";
  CompareDataMC(extra_name, xsec_dir, top_dir);
}

void protoana::PDSPThinSliceFitter::CompareDataMC(
    std::string extra_name, TDirectory * xsec_dir, TDirectory * plot_dir,
    bool post_fit) {
  fThinSliceDriver->CompareDataMC(
      fSamples, fDataSet, fOutputFile,
      fPlotStyle, (fTotalFluxParameters + fTotalSignalParameters),
      plot_dir, fPlotRebinned, post_fit);

  xsec_dir->cd();

  //Get the incident histogram from all of the relevant samples
  for (size_t i = 0; i < fMeasurementSamples.size(); ++i) {
    int sample_ID = fMeasurementSamples[i];
    auto & samples_vec_2D
        = fSamples[sample_ID]; 

    std::vector<double> signal_bins = fSignalBins[sample_ID];
    if (fDrawXSecUnderflow) signal_bins.insert(signal_bins.begin(), 0.);

    std::string signal_name = extra_name;
    signal_name += samples_vec_2D[0][1].GetName() + "Hist";
    TH1D signal_hist(signal_name.c_str(), "", signal_bins.size() - 1,
                     &signal_bins[0]);

    for (size_t j = 0; j < samples_vec_2D.size(); ++j) {
      auto & samples_vec = samples_vec_2D[j];
      for (size_t k = (fDrawXSecUnderflow? 0 : 1);
           k < samples_vec.size()-1; ++k) {
        ThinSliceSample & sample = samples_vec[k];
        if (fDrawXSecUnderflow) {
          signal_hist.AddBinContent(k+1, sample.GetVariedFlux());
        }
        else {
          signal_hist.AddBinContent(k, sample.GetVariedFlux());
        }
      }
    }
    signal_hist.Write();

    //Get the incident histogram from all of the relevant samples
    std::string inc_name = extra_name;
    inc_name += "TotalIncident" + samples_vec_2D[0][0].GetName();

    TH1D * total_incident_hist = new TH1D(inc_name.c_str(), "",
                                          signal_bins.size() - 1,
                                          &signal_bins[0]);
    std::map<int, std::vector<TH1D*>> temp_hists;
    for (size_t i = 0; i < fIncidentSamples.size(); ++i) {
      auto & vec_2D = fSamples[fIncidentSamples[i]];
      temp_hists[fIncidentSamples[i]] = std::vector<TH1D*>();
      for (size_t j = 0; j < vec_2D.size(); ++j) {
        auto & samples_vec = vec_2D[j];
        for (size_t k = 0; k < samples_vec.size(); ++k) {
          if (j == 0) {
            std::string name = extra_name;
            name += "Incident" + samples_vec_2D[0][1].GetName() +
                     samples_vec[k].GetName() + std::to_string(k);
            std::string title = samples_vec[k].GetSelectionHists().begin()->second->GetTitle();
            temp_hists[fIncidentSamples[i]].push_back(
                new TH1D(name.c_str(),
                         title.c_str(),
                         signal_bins.size() - 1,
                         &signal_bins[0]));
          }
          if (fSliceMethod == "E") {
            samples_vec[k].FillESliceHist(*total_incident_hist);
            samples_vec[k].FillESliceHist(
                *(temp_hists[fIncidentSamples[i]][k]));
          }
          else {
            samples_vec[k].FillHistFromIncidentEnergies(*total_incident_hist);
            samples_vec[k].FillHistFromIncidentEnergies(
                *(temp_hists[fIncidentSamples[i]][k]));
          }
        }
      }
    }
    total_incident_hist->Write();
    if (post_fit) {
      fBestFitIncs[sample_ID] = total_incident_hist;
    }
    else {
      fNominalIncs[sample_ID] = total_incident_hist;
    }

    std::string xsec_name = extra_name + samples_vec_2D[0][1].GetName() +
                            "XSec";
    TH1D * xsec_hist = (TH1D*)signal_hist.Clone(xsec_name.c_str());
    xsec_hist->Sumw2();
    xsec_hist->Divide(total_incident_hist);
    if (fSliceMethod == "E") {
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
        xsec_hist->SetBinContent(i, -1.*log(1. - xsec_hist->GetBinContent(i)));
      }
      xsec_hist->Scale(1.E27*39.948/(1.4 * 6.022E23));
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {

        double bethe_val = BetheBloch(xsec_hist->GetBinCenter(i), 139.57); 

        xsec_hist->SetBinContent(i, (bethe_val*
                                     xsec_hist->GetBinContent(i)/
                                     xsec_hist->GetBinWidth(i)));
      }
    }
    else if (fSliceMethod == "Alt") {
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
        xsec_hist->SetBinContent(i, -1.*log(1. - xsec_hist->GetBinContent(i)));
      }
      xsec_hist->Scale(1.E27*39.948/(fPitch * 1.4 * 6.022E23));
    }
    else {
      xsec_hist->Scale(1.E27/ (fPitch * 1.4 * 6.022E23 / 39.948 ));
    }
    xsec_hist->Write();
    if (post_fit) {
      fBestFitXSecs[sample_ID] = xsec_hist;
    }
    else {
      fNominalXSecs[sample_ID] = xsec_hist;
    }

    std::string stack_name = extra_name;
    stack_name += "IncidentStack" + samples_vec_2D[0][1].GetName();
    THStack inc_stack(stack_name.c_str(), ";Reconstructed KE (MeV)");
    int iColor = 0;
    for (auto it = temp_hists.begin(); it != temp_hists.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        std::pair<int, int> color_fill =
            fThinSliceDriver->GetColorAndStyle(iColor, fPlotStyle);
        it->second[i]->SetFillColor(color_fill.first);
        it->second[i]->SetFillStyle(color_fill.second);
        it->second[i]->SetLineColor(1);
        inc_stack.Add(it->second[i]);
        ++iColor;
      }
    }
    inc_stack.Write();
  }

  //Get the Best fit selection hists
  std::map<int, std::vector<TH1*>> temp_map;
  fThinSliceDriver->GetCurrentHists(fSamples, fDataSet, temp_map,
                                    fPlotRebinned);
  for (auto it = temp_map.begin(); it != temp_map.end(); ++it) {
    fBestFitSelectionHists[it->first] = (TH1D*)it->second[0];
  }

  //Get the best fit truth vals
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    fBestFitTruthVals[it->first]
        = std::vector<double>(it->second[0].size(), 0.);
   
    for (size_t i = 0; i < it->second[0].size(); ++i) {
      double best_fit_val_i = 0.;
      for (size_t j = 0; j < it->second.size(); ++j) {
        best_fit_val_i += it->second[j][i].GetVariedFlux();
      }

      fBestFitTruthVals[it->first][i] = best_fit_val_i;
    }
  }
}

void protoana::PDSPThinSliceFitter::RunFitAndSave() {
  TFile fMCFile(fMCFileName.c_str(), "OPEN");
  fMCTree = (TTree*)fMCFile.Get(fTreeName.c_str());


  switch (fFitFunctionType) {
    case 1:
      DefineFitFunction();
      break;
    case 2:
      DefineFitFunction2();
      break;
    default:
      DefineFitFunction();
      break;
  }

  fMinimizer->SetFunction(fFitFunction);
  int fit_status = fMinimizer->Minimize();

  if (!fit_status) {
    std::cout << "Failed to find minimum" << std::endl;
    ParameterScans();
  }
  else {
    std::vector<double> vals;
    std::cout << "Found minimimum: " << std::endl;
    size_t total_parameters = fTotalSignalParameters + fTotalFluxParameters + fTotalSystParameters;
    for (size_t i = 0; i < total_parameters; ++i) {
      std::cout << fMinimizer->VariableName(i) << " " << fMinimizer->X()[i] <<
                   std::endl;
      vals.push_back(fMinimizer->X()[i]);
    }

    TH2D covHist("covHist", "", total_parameters, 0,
                 total_parameters, total_parameters, 0,
                 total_parameters);
    TH2D corrHist("corrHist", "", total_parameters, 0,
                  total_parameters, total_parameters, 0,
                  total_parameters);

    TH1D parsHist("postFitPars", "", total_parameters, 0,
                  total_parameters);
    TH1D parsHist_normal_val("postFitParsNormal", "", 
                             total_parameters, 0, total_parameters);
    int n_par = 0;
    for (auto it = fSignalParameters.begin();
         it != fSignalParameters.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        parsHist.SetBinContent(n_par+1, fMinimizer->X()[n_par]);
        parsHist.SetBinError(n_par+1,
                             sqrt(fMinimizer->CovMatrix(n_par, n_par)));
        parsHist.GetXaxis()->SetBinLabel(
            n_par+1, fSignalParameterNames[it->first][i].c_str());

        parsHist_normal_val.SetBinContent(n_par+1, fMinimizer->X()[n_par]);
        parsHist_normal_val.SetBinError(
            n_par+1, sqrt(fMinimizer->CovMatrix(n_par, n_par)));
        parsHist_normal_val.GetXaxis()->SetBinLabel(
            n_par+1, fSignalParameterNames[it->first][i].c_str());


        covHist.GetXaxis()->SetBinLabel(
            n_par+1, fSignalParameterNames[it->first][i].c_str());
        corrHist.GetXaxis()->SetBinLabel(
            n_par+1, fSignalParameterNames[it->first][i].c_str());
        covHist.GetYaxis()->SetBinLabel(
            n_par+1, fSignalParameterNames[it->first][i].c_str());
        corrHist.GetYaxis()->SetBinLabel(
            n_par+1, fSignalParameterNames[it->first][i].c_str());

        ++n_par;
      }
    }
    for (auto it = fFluxParameters.begin();
         it != fFluxParameters.end(); ++it) {
      parsHist.SetBinContent(n_par+1, fMinimizer->X()[n_par]);
      parsHist.SetBinError(n_par+1,
                           sqrt(fMinimizer->CovMatrix(n_par, n_par)));
      parsHist.GetXaxis()->SetBinLabel(
          n_par+1, fFluxParameterNames[it->first].c_str());

      parsHist_normal_val.SetBinContent(n_par+1, fMinimizer->X()[n_par]);
      parsHist_normal_val.SetBinError(
          n_par+1, sqrt(fMinimizer->CovMatrix(n_par, n_par)));
      parsHist_normal_val.GetXaxis()->SetBinLabel(
          n_par+1, fFluxParameterNames[it->first].c_str());



      covHist.GetXaxis()->SetBinLabel(
          n_par+1, fFluxParameterNames[it->first].c_str());
      corrHist.GetXaxis()->SetBinLabel(
          n_par+1, fFluxParameterNames[it->first].c_str());
      covHist.GetYaxis()->SetBinLabel(
          n_par+1, fFluxParameterNames[it->first].c_str());
      corrHist.GetYaxis()->SetBinLabel(
          n_par+1, fFluxParameterNames[it->first].c_str());

      ++n_par;
    }

    for (auto it = fSystParameters.begin();
         it != fSystParameters.end(); ++it) {

      parsHist.SetBinContent(n_par+1, fMinimizer->X()[n_par]/it->second.GetCentral());
      parsHist.SetBinError(n_par+1,
                           sqrt(fMinimizer->CovMatrix(n_par, n_par))/it->second.GetCentral());
      parsHist.GetXaxis()->SetBinLabel(
          n_par+1, it->second.GetName().c_str());

      parsHist_normal_val.SetBinContent(n_par+1, fMinimizer->X()[n_par]);
      parsHist_normal_val.SetBinError(
          n_par+1, sqrt(fMinimizer->CovMatrix(n_par, n_par)));
      parsHist_normal_val.GetXaxis()->SetBinLabel(
          n_par+1, it->second.GetName().c_str());



      covHist.GetXaxis()->SetBinLabel(
          n_par+1, it->second.GetName().c_str());
      corrHist.GetXaxis()->SetBinLabel(
          n_par+1, it->second.GetName().c_str());
      covHist.GetYaxis()->SetBinLabel(
          n_par+1, it->second.GetName().c_str());
      corrHist.GetYaxis()->SetBinLabel(
          n_par+1, it->second.GetName().c_str());

      ++n_par;
    }

    TMatrixD * cov = new TMatrixD(total_parameters, total_parameters);

    for (size_t i = 0; i < total_parameters; ++i) {
      for (size_t j = 0; j < total_parameters; ++j) {
        double cov_val = fMinimizer->CovMatrix(i, j);
        covHist.SetBinContent(i+1, j+1, cov_val);
        corrHist.SetBinContent(i+1, j+1, fMinimizer->Correlation(i, j));
        (*cov)(i, j) = fMinimizer->CovMatrix(i, j);
      }
    }

    fOutputFile.cd();
    parsHist.SetFillColor(kRed);
    parsHist.SetFillStyle(3144);
    parsHist.SetMarkerColor(kBlack);
    parsHist.SetMarkerStyle(20);
    parsHist.Write();

    fBestFitSignalPars = fSignalParameters; 
    fBestFitFluxPars = fFluxParameters; 
    fBestFitSystPars = fSystParameters;

    //Drawing pre + post fit pars
    TCanvas cPars("cParameters", "");
    cPars.SetTicks();
    parsHist.GetYaxis()->SetTitle("Parameter Values");
    parsHist.SetTitleSize(.05, "Y");
    parsHist.Draw("e2");
    TH1D * preFitParsHist = (TH1D*)fOutputFile.Get("preFitPars");
    preFitParsHist->Draw("p same");
    TLine l(0., 1., parsHist.GetXaxis()->GetBinUpEdge(parsHist.GetNbinsX()), 1.);
    l.SetLineColor(kBlack);
    l.Draw("same");
    TLegend leg(.15, .65, .45, .85);
    leg.AddEntry(preFitParsHist, "Pre-Fit", "p");
    leg.AddEntry(&parsHist, "Post-Fit #pm 1 #sigma", "lpf");
    leg.Draw();
    cPars.Write();

    covHist.Write();
    cov->Write("CovMatrix");
    corrHist.Write();

    ParameterScans();



    fFillIncidentInFunction = true;
    fFitFunction(&vals[0]);

    //SetBestFit();  //Deprecated

    //save post fit stacks
    TDirectory * top_dir = fOutputFile.GetDirectory("");
    TDirectory * xsec_dir = fOutputFile.mkdir("PostFitXSec");
    std::string extra_name = "PostFit";
    CompareDataMC(extra_name, xsec_dir, top_dir, true);
    if (fDoThrows) 
      DoThrows(parsHist_normal_val, cov);

    if (fDo1DShifts)
      Do1DShifts(parsHist_normal_val);
  }

  if (fDoSysts)
    fThinSliceDriver->WrapUpSysts(fOutputFile);

  fMCFile.Close();
}

std::vector<double> protoana::PDSPThinSliceFitter::GetBestFitParsVec() {
  std::vector<double> results;

  for (auto it = fSignalParameters.begin();
       it != fSignalParameters.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      results.push_back(it->second.at(i));
    }
  }
  for (auto it = fFluxParameters.begin();
       it != fFluxParameters.end(); ++it) {
    results.push_back(it->second);
  }

  for (auto it = fSystParameters.begin();
       it != fSystParameters.end(); ++it) {
    results.push_back(it->second.GetValue());
  }

  return results;
}

void protoana::PDSPThinSliceFitter::ParameterScans() {
  fOutputFile.cd();
  TDirectory * out = (TDirectory *)fOutputFile.mkdir("Scans");
  out->cd();

  size_t total_parameters = fTotalSignalParameters + fTotalFluxParameters + fTotalSystParameters;
  std::cout << "Total parameters " << total_parameters << std::endl;

  double * x = new double[fNScanSteps] {};
  double * y = new double[fNScanSteps] {};
  for (size_t i = 0; i < total_parameters; ++i) {
    std::cout << "\tParameter " << fMinimizer->VariableName(i) << std::endl;
    bool scanned = fMinimizer->Scan(i, fNScanSteps, x, y);
    if (scanned) {
      TGraph gr(fNScanSteps - 1, x, y);
      gr.Write(fMinimizer->VariableName(i).c_str());
    }
  }

  delete[] x;
  delete[] y;
}

void protoana::PDSPThinSliceFitter::DoThrows(const TH1D & pars, const TMatrixD * cov) {
  std::vector<std::vector<double>> vals;
  TVectorD rand(pars.GetNbinsX());

  TDecompChol chol(*cov);
  bool success = chol.Decompose();
  if (!success) {
    std::cout << "Error" << std::endl;
    return;
  }

  fOutputFile.mkdir("Throws");
  fOutputFile.cd("Throws");

  TCanvas cPars("cParametersThrown", "");
  cPars.SetTicks();

  //Build up the map of selection hists
  std::map<int, std::vector<TH1*>> throw_hists;
  for (auto it = fDataSet.GetSelectionHists().begin();
       it != fDataSet.GetSelectionHists().end(); ++it) {
    throw_hists[it->first] = std::vector<TH1*>();
  }

  std::map<int, std::vector<TH1*>> truth_throw_hists;
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    truth_throw_hists[it->first] = std::vector<TH1*>();
  }

  std::map<int, std::vector<TH1*>> truth_inc_throw_hists,
                                   truth_xsec_throw_hists;
  for (auto it = fMeasurementSamples.begin();
       it != fMeasurementSamples.end(); ++it) {
    truth_inc_throw_hists[*it] = std::vector<TH1*>();
    truth_xsec_throw_hists[*it] = std::vector<TH1*>();
  }

  int n_signal_flux_pars = fTotalSignalParameters + fTotalFluxParameters;

  for (size_t i = 0; i < fNThrows; ++i) {
    if (! (i%100) ) {
      std::cout << "Throw " << i << std::endl;
    //auto start = std::chrono::high_resolution_clock::now();
    }
    vals.push_back(std::vector<double>(pars.GetNbinsX(), 1.));
    
    bool rethrow = true;
    size_t nRethrows = 0;
    while (rethrow && nRethrows < fMaxRethrows) {
      bool all_pos = true;
      for (int j = 1; j <= pars.GetNbinsX(); ++j) {
        rand[j-1] = fRNG.Gaus();
      }

      TVectorD rand_times_chol = chol.GetU()*rand;

      for (int j = 1; j <= pars.GetNbinsX(); ++j) {
        //std::cout << rand_times_chol[j-1] << " ";
        vals.back()[j-1] = pars.GetBinContent(j) + rand_times_chol[j-1];

        //Configure this to truncate at 0 or rethrow
        if ((j-1 < n_signal_flux_pars) && (vals.back()[j-1] < 0.)) {
          all_pos = false;
        }
        rethrow = !all_pos;

      }
      //std::cout << std::endl;
      ++nRethrows;
    }
    
    //for (double v : vals.back()) {
    //  std::cout << v << " ";
    //}
    //std::cout << std::endl;

    //Applies the variations according to the thrown parameters
    fFillIncidentInFunction = true;
    fFitFunction(&(vals.back()[0]));
    fThinSliceDriver->GetCurrentHists(fSamples, fDataSet, throw_hists, fPlotRebinned);
    GetCurrentTruthHists(truth_throw_hists,
                         truth_inc_throw_hists,
                         truth_xsec_throw_hists);
    for(auto it = truth_xsec_throw_hists.begin();
        it != truth_xsec_throw_hists.end(); ++it) { 
      TH1 * xsec_hist = it->second.back();
      if (fSliceMethod == "E") {
        for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
          xsec_hist->SetBinContent(i, -1.*log(1. - xsec_hist->GetBinContent(i)));
        }
        xsec_hist->Scale(1.E27*39.948/(1.4 * 6.022E23));
        for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {

          double bethe_val = BetheBloch(xsec_hist->GetBinCenter(i), 139.57); 

          xsec_hist->SetBinContent(i, (bethe_val*
                                       xsec_hist->GetBinContent(i)/
                                       xsec_hist->GetBinWidth(i)));
        }
      }
      else if (fSliceMethod == "Alt") {
        for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
          xsec_hist->SetBinContent(i, -1.*log(1. - xsec_hist->GetBinContent(i)));
        }
        xsec_hist->Scale(1.E27*39.948/(fPitch * 1.4 * 6.022E23));
      }
      else {
        xsec_hist->Scale(1.E27/ (fPitch * 1.4 * 6.022E23 / 39.948 ));
      }
    }

    //auto end = std::chrono::high_resolution_clock::now();
    //auto difference = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    //std::cout << "Time taken: " << difference << std::endl;
  }

  PlotThrows(throw_hists, truth_throw_hists,
             truth_inc_throw_hists, truth_xsec_throw_hists);
}

void protoana::PDSPThinSliceFitter::Do1DShifts(const TH1D & pars) {
  TDirectory * shift_dir = fOutputFile.mkdir("1DShifts");
  fFillIncidentInFunction = false;

  //Iterate over the bins in the parameter hist
  for (int i = 1; i <= pars.GetNbinsX(); ++i) {
    
    std::string par_name = pars.GetXaxis()->GetBinLabel(i);

    std::vector<double> vals;

    //Set to +1 sigma for parameter i
    for (int j = 1; j <= pars.GetNbinsX(); ++j) {
      if (j == i) {
        vals.push_back(pars.GetBinContent(j) + pars.GetBinError(i));
      }
      else {
        vals.push_back(pars.GetBinContent(j));
      }
    }
    fFitFunction(&vals[0]);

    shift_dir->cd();
    std::string dir_name = par_name + "_PlusSigma";
    TDirectory * plus_dir = shift_dir->mkdir(dir_name.c_str()); 
    TDirectory * xsec_dir = plus_dir->mkdir("XSec");
    CompareDataMC(dir_name, xsec_dir, plus_dir);

    vals.clear();

    //Set to -1 sigma for parameter i
    for (int j = 1; j <= pars.GetNbinsX(); ++j) {
      if (j == i) {
        vals.push_back(pars.GetBinContent(j) - pars.GetBinError(i));
      }
      else {
        vals.push_back(pars.GetBinContent(j));
      }
    }
    fFitFunction(&vals[0]);

    shift_dir->cd();
    dir_name = par_name + "_MinusSigma";
    TDirectory * minus_dir = shift_dir->mkdir(dir_name.c_str()); 
    xsec_dir = minus_dir->mkdir("XSec");
    CompareDataMC(dir_name, xsec_dir, minus_dir);
  }
}

void protoana::PDSPThinSliceFitter::DefineFitFunction() {

  fFitFunction = ROOT::Math::Functor(
      [&](double const * coeffs) {

        //Set all the parameters
        //coeffs are ordered according to
        //keys of signal par map
        //----vector for given signal
        //then keys of flux par map
        size_t par_position = 0;
        for (auto it = fSignalParameters.begin();
             it != fSignalParameters.end(); ++it) {
          for (size_t i = 0; i < it->second.size(); ++i) {
            it->second.at(i) = coeffs[par_position];
            ++par_position;
          }
        }
        for (auto it = fFluxParameters.begin();
             it != fFluxParameters.end(); ++it) {
          it->second = coeffs[par_position];
          ++par_position;
        }

        double total_varied_flux = 0.;
        double total_nominal_flux = 0.;
        std::vector<double> varied_flux(fBeamEnergyBins.size() - 1, 0.);
        std::vector<double> nominal_flux(fBeamEnergyBins.size() - 1, 0.);

        //Go through and determine the scaled flux
        for (auto it = fFluxesBySample.begin();
             it != fFluxesBySample.end(); ++it) {
          int sample_ID = it->first;
          std::vector<std::vector<ThinSliceSample>> & samples_2D
              = fSamples[sample_ID];
          for (size_t i = 0; i < samples_2D.size(); ++i) {
            std::vector<ThinSliceSample> & samples = samples_2D[i];
            for (size_t j = 0; j < samples.size(); ++j) {
              ThinSliceSample & sample = samples.at(j);
              int flux_type = sample.GetFluxType();
              if (flux_type == 1) {
                nominal_flux[i] += fFluxesBySample[sample_ID][i][j];
              }
              total_nominal_flux += fFluxesBySample[sample_ID][i][j];
              
              if (fSignalParameters.find(sample_ID) != fSignalParameters.end()) {
                //Scale the flux for this signal channel
                //check for under/overflow bin here. Don't scale by signal pars
                if (j == 0 || j == samples.size() - 1) {
                  if (flux_type == 1) {
                    varied_flux[i]
                        += (fFluxesBySample[sample_ID][i][j]);
                  }
                  total_varied_flux += fFluxesBySample[sample_ID][i][j];
                }
                else {
                  if (flux_type == 1) {
                    varied_flux[i]
                        += (fFluxesBySample[sample_ID][i][j] *
                            fSignalParameters[sample_ID][j-1]);
                  }
                  total_varied_flux += (fFluxesBySample[sample_ID][i][j]*
                                        fSignalParameters[sample_ID][j-1]);
                }
              }
              else if (fFluxParameters.find(flux_type) != fFluxParameters.end()) {
                if (flux_type == 1) {
                  varied_flux[i]
                      += (fFluxesBySample[sample_ID][i][j] *
                          fFluxParameters[flux_type]);
                }
                total_varied_flux += (fFluxesBySample[sample_ID][i][j]*
                                      fFluxParameters[flux_type]);
              }
              else {
                if (flux_type == 1) {
                  varied_flux[i]
                      += fFluxesBySample[sample_ID][i][j];
                }
                total_varied_flux += fFluxesBySample[sample_ID][i][j];
              }
            }
          }
        }

        double total_flux_factor = total_nominal_flux / total_varied_flux;

        double nominal_primary = 0., varied_primary = 0.;
        for (size_t i = 0; i < nominal_flux.size(); ++i) {
          nominal_primary += nominal_flux[i];
          varied_primary += varied_flux[i];
        }

        std::vector<double> flux_factor = nominal_flux;
        for (size_t i = 0; i < flux_factor.size(); ++i) {
          flux_factor[i] = (nominal_flux[i] > 1.e-7 ?
                            flux_factor[i] / varied_flux[i] :
                            0.)*(varied_primary/nominal_primary);
        }

        for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
          //Check if it's a signal
          int sample_ID = it->first;
          std::vector<std::vector<ThinSliceSample>> & samples_2D = it->second;

          if (fSignalParameters.find(sample_ID) != fSignalParameters.end()) {
            for (size_t i = 0; i < samples_2D.size(); ++i) {
              std::vector<ThinSliceSample> & samples = samples_2D[i]; 
              for (size_t j = 0; j < samples.size(); ++j) {
                //check for under/overflow bin here. Don't scale by signal pars
                double val = ((j == 0 || j == samples.size() - 1) ?
                              1. : fSignalParameters[sample_ID][j-1]);
                ThinSliceSample & sample = samples.at(j);
                int flux_type = sample.GetFluxType();
                if (flux_type == 1) {
                  sample.SetFactorAndScale(flux_factor[i]*val*
                                           total_flux_factor);
                }
                else {
                  sample.SetFactorAndScale(val*total_flux_factor);
                }
              }
            }
          }
          else {
            bool flux_sample = false; 
            //Check if this is a sample that should be scaled by a flux parameter
            for (auto it2 = fFluxParsToSamples.begin();
                 it2 != fFluxParsToSamples.end(); ++it2) {
              int flux_type = it2->first;
              std::vector<int> flux_samples = it2->second;
              double val = fFluxParameters[flux_type];
              if (std::find(flux_samples.begin(), flux_samples.end(),
                            sample_ID) != flux_samples.end()) {

                for (size_t i = 0; i < samples_2D.size(); ++i) {
                  std::vector<ThinSliceSample> & samples = samples_2D[i];
                  for (size_t j = 0; j < samples.size(); ++j) {
                    ThinSliceSample & sample = samples.at(j);
                    if (flux_type == 1) {
                      sample.SetFactorAndScale(flux_factor[i]*val*
                                               total_flux_factor);
                    }
                    else {
                      sample.SetFactorAndScale(val*total_flux_factor);
                    }
                  }
                }

                flux_sample = true;
                break;
              }
            }

            if (!flux_sample) {
              for (size_t i = 0; i < samples_2D.size(); ++i) {
                std::vector<ThinSliceSample> & samples = samples_2D[i];
                for (size_t j = 0; j < samples.size(); ++j) {
                  ThinSliceSample & sample = samples.at(j);
                  int flux_type = sample.GetFluxType();
                  if (flux_type == 1) {
                    sample.SetFactorAndScale(flux_factor[i]*
                                             total_flux_factor);
                  }
                  else {
                    sample.SetFactorAndScale(total_flux_factor);
                  }
                }
              }
            }
          }
        }

        /*
        double final_nominal_total = 0., final_varied_total = 0.;
        for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
          for (size_t i = 0; i < it->second.size(); ++i) {
            for (size_t j = 0; j < it->second[i].size(); ++j) {
              final_nominal_total += it->second[i][j].GetNominalFlux();
              final_varied_total += it->second[i][j].GetVariedFlux();
            }
          }
        }*/

        std::pair<double, size_t> chi2_points
            = fThinSliceDriver->CalculateChi2(fSamples, fDataSet);
        std::cout << "Fit step: " << fNFitSteps << " " << chi2_points.first << std::endl;
        ++fNFitSteps;

        return (chi2_points.first);
      },
      fTotalSignalParameters + fTotalFluxParameters);

}

void protoana::PDSPThinSliceFitter::DefineFitFunction2() {
  std::cout << "FF2" << std::endl;

  fFitFunction = ROOT::Math::Functor(
      [&](double const * coeffs) {

        //Set all the parameters
        //coeffs are ordered according to
        //keys of signal par map
        //----vector for given signal
        //then keys of flux par map
        //then values of syst parameters
        size_t par_position = 0;
        for (auto it = fSignalParameters.begin();
             it != fSignalParameters.end(); ++it) {
          for (size_t i = 0; i < it->second.size(); ++i) {
            it->second.at(i) = coeffs[par_position];
            ++par_position;
          }
        }
        for (auto it = fFluxParameters.begin();
             it != fFluxParameters.end(); ++it) {
          it->second = coeffs[par_position];
          ++par_position;
        }

        for (auto it = fSystParameters.begin();
             it != fSystParameters.end(); ++it) {
          it->second.SetValue(coeffs[par_position]);
          ++par_position;
        }

        //Refill the hists 
        fThinSliceDriver->RefillMCSamples(
            fEvents, fSamples, fIsSignalSample,
            fBeamEnergyBins, fSignalParameters,
            fFluxParameters, fSystParameters,
            fFillIncidentInFunction);

        //Scale to Data
        for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
          for (size_t i = 0; i < it->second.size(); ++i) {
            for (size_t j = 0; j < it->second[i].size(); ++j) {
              it->second[i][j].ScaleToDataMC();
            }
          }
        }


        double total_varied_flux = 0.;
        double total_nominal_flux = 0.;
        std::vector<double> varied_flux(fBeamEnergyBins.size() - 1, 0.);
        std::vector<double> nominal_flux(fBeamEnergyBins.size() - 1, 0.);

        //Go through and determine the scaled flux
        for (auto it = fFluxesBySample.begin();
             it != fFluxesBySample.end(); ++it) {
          int sample_ID = it->first;
          std::vector<std::vector<ThinSliceSample>> & samples_2D
              = fSamples[sample_ID];
          for (size_t i = 0; i < samples_2D.size(); ++i) {
            std::vector<ThinSliceSample> & samples = samples_2D[i];
            for (size_t j = 0; j < samples.size(); ++j) {
              ThinSliceSample & sample = samples.at(j);
              int flux_type = sample.GetFluxType();
              if (flux_type == 1) {
                nominal_flux[i] += fFluxesBySample[sample_ID][i][j];
                varied_flux[i] += sample.GetVariedFlux();
              }
              total_nominal_flux += fFluxesBySample[sample_ID][i][j];
              total_varied_flux += sample.GetVariedFlux();
              
            }
          }
        }

        double total_flux_factor = total_nominal_flux / total_varied_flux;
        //std::cout << "Flux: " << total_nominal_flux << " " <<
        //             total_varied_flux << " " << total_flux_factor << std::endl;

        double nominal_primary = 0., varied_primary = 0.;
        for (size_t i = 0; i < nominal_flux.size(); ++i) {
          nominal_primary += nominal_flux[i];
          varied_primary += varied_flux[i];
        }

        std::vector<double> flux_factor = nominal_flux;
        for (size_t i = 0; i < flux_factor.size(); ++i) {
          flux_factor[i] = (nominal_flux[i] > 1.e-7 ?
                            flux_factor[i] / varied_flux[i] :
                            0.)*(varied_primary/nominal_primary);
        }

        for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
          for (size_t i = 0; i < it->second.size(); ++i) {
            for (size_t j = 0; j < it->second[i].size(); ++j) {
              int flux_type = it->second[i][j].GetFluxType();
              std::cout << std::fixed;
              //double total = 0.;
              //std::cout << "Varied flux1 : " << std::setprecision(8) <<it->second[i][j].GetVariedFlux() << std::endl;
              //std::cout << "scaling by " << total_flux_factor << " " << (flux_type == 1 ? flux_factor[i] : 1.) <<
              //             " " << total_flux_factor*(flux_type == 1 ? flux_factor[i] : 1.) << std::endl;
              it->second[i][j].SetFactorAndScale(
                  total_flux_factor*(flux_type == 1 ? flux_factor[i] : 1.));
              //std::cout << "Varied flux2 : " << std::setprecision(8) <<it->second[i][j].GetVariedFlux() << std::endl;
            }
          }
        }

        std::pair<double, size_t> chi2_points
            = fThinSliceDriver->CalculateChi2(fSamples, fDataSet);
        //std::cout << "Fit step: " << fNFitSteps << " " << chi2_points.first << std::endl;
        //for (size_t i = 0;
        //     i < fTotalSignalParameters + fTotalFluxParameters + fTotalSystParameters; ++i) {
        //  std::cout << coeffs[i] << " ";
        //}
        //std::cout << std::endl;
        ++fNFitSteps;

        return (chi2_points.first);
      },
      fTotalSignalParameters + fTotalFluxParameters + fTotalSystParameters);

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

  fSelectionSets = pset.get<std::vector<fhicl::ParameterSet>>("Selections");

  fIncidentRecoBins = pset.get<std::vector<double>>("IncidentRecoBins");
  fTrueIncidentBins = pset.get<std::vector<double>>("TrueIncidentBins");
  fBeamEnergyBins = pset.get<std::vector<double>>("BeamEnergyBins");
  fIncidentSamples = pset.get<std::vector<int>>("IncidentSamples");
  fMeasurementSamples = pset.get<std::vector<int>>("MeasurementSamples");

  fDrawXSecUnderflow = pset.get<bool>("DrawXSecUnderflow");
  
  fSampleSets = pset.get<std::vector<fhicl::ParameterSet>>("Samples");
  std::vector<std::pair<int, std::string>> temp_vec
      = pset.get<std::vector<std::pair<int, std::string>>>("FluxTypes");
  fFluxTypes = std::map<int, std::string>(temp_vec.begin(), temp_vec.end());

  fFitFlux = pset.get<bool>("FitFlux");
  if (fFitFlux) {
    for (size_t i = 0; i < temp_vec.size() - 1; ++i) {
      fFluxParameterNames[temp_vec[i].first] = "par_" + temp_vec[i].second +
                                               "_flux";
      fFluxParameters[temp_vec[i].first] = 1.;
    }
    fTotalFluxParameters = temp_vec.size() - 1;
  }

  fMaxCalls = pset.get<int>("MaxCalls");
  fNScanSteps = pset.get<unsigned int>("NScanSteps") + 1;
  fTolerance = pset.get<double>("Tolerance");
  fLowerLimit = pset.get<double>("LowerLimit");
  fUpperLimit = pset.get<double>("UpperLimit");
  fPlotStyle = pset.get<std::vector<std::pair<int,int>>>("PlotStyle");
  fPlotRebinned = pset.get<bool>("PlotRebinned");
  fRandomStart = pset.get<bool>("RandomStart");

  fDriverName = pset.get<std::string>("DriverName");
  fAnalysisOptions = pset.get<fhicl::ParameterSet>("AnalysisOptions");
  fPitch = fAnalysisOptions.get<double>("WirePitch");
  fSliceMethod = fAnalysisOptions.get<std::string>("SliceMethod");
  fDoFakeData = pset.get<bool>("DoFakeData");
  fDoThrows = pset.get<bool>("DoThrows");
  fDo1DShifts = pset.get<bool>("Do1DShifts");
  fDoSysts = pset.get<bool>("DoSysts");
  fFitFunctionType = pset.get<int>("FitFunctionType");

  fNThrows = pset.get<size_t>("NThrows");
  fMaxRethrows = pset.get<size_t>("MaxRethrows");

  //Initialize systematics
  if (fDoSysts) {
    std::vector<fhicl::ParameterSet> par_vec
        = pset.get<std::vector<fhicl::ParameterSet>>("Systematics");
    for (size_t i = 0; i < par_vec.size(); ++i) {
      ThinSliceSystematic syst(par_vec[i]);
      fSystParameters[par_vec[i].get<std::string>("Name")] = syst;
      ++fTotalSystParameters;
    }
  }
}

void protoana::PDSPThinSliceFitter::SetBestFit() {
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        it->second[i][j].SetBestFit();
      }
    }
  }
}


void protoana::PDSPThinSliceFitter::PlotThrows(
    std::map<int, std::vector<TH1*>> & throw_hists,
    std::map<int, std::vector<TH1*>> & truth_throw_hists,
    std::map<int, std::vector<TH1*>> & truth_inc_hists,
    std::map<int, std::vector<TH1*>> & truth_xsec_hists) {
  std::map<int, TH1*> data_hists
      = (fPlotRebinned ?
         fDataSet.GetRebinnedSelectionHists() :
         fDataSet.GetSelectionHists());

  //Build best fit hists and get bins for covariance 
  //std::map<int, TH1D*> best_fit_selection_hists;
  int nBins = 0;
  for (auto it = data_hists.begin(); it != data_hists.end(); ++it ) {
    //TH1D * best_fit_hist = (TH1D*)it->second->Clone();
    //best_fit_hist->Reset();
    //for (auto it2 = fSamples.begin(); it2 != fSamples.end(); ++it2) {
    //  for (size_t i = 0; i < it2->second.size(); ++i) {
    //    for (size_t j = 0; j < it2->second[i].size(); ++j) {
    //      it2->second[i][j].SetFactorToBestFit();
    //      best_fit_hist->Add(
    //          (TH1D*)(fPlotRebinned ?
    //                  it2->second[i][j].GetRebinnedSelectionHist(it->first) :
    //                  it2->second[i][j].GetSelectionHist(it->first)));
    //    }
    //  }
    //}
    //best_fit_selection_hists[it->first] = best_fit_hist;
    //nBins += best_fit_hist->GetNbinsX();
    nBins += it->second->GetNbinsX();
  }

  TH2D selection_cov("SelectionCov", "", nBins, 0, nBins, nBins, 0, nBins);

  nBins = 0;
  std::map<int, size_t> sample_bins;
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    nBins += it->second[0].size();
    sample_bins[it->first] = it->second[0].size();
  }

  //std::map<int, std::vector<double>> best_fit_truth;
  std::map<int, std::vector<double>> best_fit_errs;

  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    //best_fit_truth[it->first]
    //    = std::vector<double>(sample_bins[it->first], 0.);
    best_fit_errs[it->first]
        = std::vector<double>(sample_bins[it->first], 0.);
   
    //for (size_t i = 0; i < sample_bins[it->first]; ++i) {
    //  double best_fit_val_i = 0.;
    //  for (size_t j = 0; j < it->second.size(); ++j) {
    //    best_fit_val_i += it->second[j][i].GetVariedFlux();
    //  }

    //  best_fit_truth[it->first][i] = best_fit_val_i;
    //}
  }

  TH2D interaction_cov("interaction_cov", "", nBins, 0, nBins, nBins, 0, nBins);
  std::map<int, std::vector<double>> best_fit_inc_truth;
  std::map<int, std::vector<double>> best_fit_xsec_truth;
  std::map<int, std::vector<double>> best_fit_inc_errs;
  std::map<int, std::vector<double>> best_fit_xsec_errs;

  nBins = 0;
  std::map<int, size_t> xsec_bins;
  for (auto it = fBestFitIncs.begin(); it != fBestFitIncs.end(); ++it) {
    int s = it->first;
    nBins += it->second->GetNbinsX();
    xsec_bins[s] = it->second->GetNbinsX();

    best_fit_inc_truth[s] = std::vector<double>(xsec_bins[s], 0.);
    best_fit_xsec_truth[s] = std::vector<double>(xsec_bins[s], 0.);
    best_fit_inc_errs[s] = std::vector<double>(xsec_bins[s], 0.);
    best_fit_xsec_errs[s] = std::vector<double>(xsec_bins[s], 0.);
    
    for (size_t i = 0; i < xsec_bins[s]; ++i) {
      best_fit_inc_truth[s][i] = it->second->GetBinContent(i+1);
      best_fit_xsec_truth[s][i] = fBestFitXSecs[s]->GetBinContent(i+1);
    }
  }

  //TH2D incident_cov("incident_cov", "", nBins, 0, nBins, nBins, 0, nBins);
  TH2D xsec_cov("xsec_cov", "", nBins, 0, nBins, nBins, 0, nBins);

  for (size_t z = 0; z < fNThrows; ++z) {
    int bin_i = 1;
    for (auto it = fBestFitSelectionHists.begin();
         it != fBestFitSelectionHists.end(); ++it) {
      TH1D * best_fit = it->second;
      int selection_ID = it->first;
      std::vector<TH1*> & temp_throws = throw_hists[selection_ID];
      for (int i = 1; i <= best_fit->GetNbinsX(); ++i) {
        double best_fit_val_i = best_fit->GetBinContent(i);
        int bin_j = 1;
        for (auto it2 = fBestFitSelectionHists.begin();
             it2 != fBestFitSelectionHists.end(); ++it2) {

          TH1D * best_fit_2 = it2->second;
          int selection_ID_2 = it2->first;
          std::vector<TH1*> & temp_throws_2 = throw_hists[selection_ID_2];
          for (int j = 1; j <= best_fit_2->GetNbinsX(); ++j) {
            double best_fit_val_j = best_fit_2->GetBinContent(j);
            double val = (best_fit_val_i - temp_throws[z]->GetBinContent(i))*
                         (best_fit_val_j - temp_throws_2[z]->GetBinContent(j));
            selection_cov.SetBinContent(
                bin_i, bin_j, (val/temp_throws.size() +
                               selection_cov.GetBinContent(bin_i, bin_j)));
            ++bin_j;
          }
        }
        ++bin_i;
      }
    }

    bin_i = 1;
    for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
      std::vector<TH1 *> throw_hists_i = truth_throw_hists[it->first];
     
      for (size_t i = 0; i < sample_bins[it->first]; ++i) {
        double best_fit_val_i = fBestFitTruthVals[it->first][i];

        int bin_j = 1;
        for (auto it2 = fSamples.begin(); it2 != fSamples.end(); ++it2) {
          std::vector<TH1 *> throw_hists_j = truth_throw_hists[it2->first];
          for (size_t j = 0; j < sample_bins[it2->first]; ++j) {
            double best_fit_val_j = fBestFitTruthVals[it2->first][j];

            double val
                = (throw_hists_i[z]->GetBinContent(i+1) - best_fit_val_i)*
                  (throw_hists_j[z]->GetBinContent(j+1) - best_fit_val_j);
            interaction_cov.SetBinContent(
                bin_i, bin_j,
                (interaction_cov.GetBinContent(bin_i, bin_j) +
                 val/throw_hists_i.size()));
            if (bin_i == bin_j && (z == fNThrows - 1)) {
              best_fit_errs[it->first][i]
                  = sqrt(interaction_cov.GetBinContent(bin_i, bin_j));
            }
            ++bin_j;
          }
        }

        ++bin_i;
      }
    }

    bin_i = 1;
    for (auto it = truth_inc_hists.begin(); it != truth_inc_hists.end(); ++it) {
      //std::vector<TH1 *> inc_hists_i = it->second;
      std::vector<TH1 *> xsec_hists_i = truth_xsec_hists[it->first];

      for (size_t i = 0; i < xsec_bins[it->first]; ++i) {
        //double best_fit_inc_i = best_fit_inc_truth[it->first][i];
        double best_fit_xsec_i = best_fit_xsec_truth[it->first][i];

        int bin_j = 1;
        for (auto it2 = truth_inc_hists.begin(); it2 != truth_inc_hists.end();
             ++it2) {
          std::vector<TH1 *> xsec_hists_j = truth_xsec_hists[it2->first];
          for (size_t j = 0; j < xsec_bins[it2->first]; ++j) {
            double best_fit_xsec_j = best_fit_xsec_truth[it2->first][j];

            double val
                = (xsec_hists_i[z]->GetBinContent(i+1) - best_fit_xsec_i)*
                  (xsec_hists_j[z]->GetBinContent(j+1) - best_fit_xsec_j);
            xsec_cov.SetBinContent(
                bin_i, bin_j,
                (xsec_cov.GetBinContent(bin_i, bin_j) +
                 val/fNThrows));
            if (bin_i == bin_j && (z == fNThrows - 1)) {
              best_fit_xsec_errs[it->first][i]
                  = sqrt(xsec_cov.GetBinContent(bin_i, bin_j));
            }
            ++bin_j;
          }
        }
        ++bin_i;
      }
    }
  }


  fOutputFile.cd("Throws");
  selection_cov.Write();
  interaction_cov.Write();
  xsec_cov.Write();

  int bin_count = 0;
  for (auto it = data_hists.begin(); it != data_hists.end(); ++it) {
    int selection_ID = it->first;
    std::vector<TH1*> hists = throw_hists.at(selection_ID);

    std::string canvas_name = "cThrow" +
                              fDataSet.GetSelectionName(selection_ID);
    TCanvas cThrow(canvas_name.c_str(), "");
    cThrow.SetTicks();

    std::string name = "Throw" + fDataSet.GetSelectionName(selection_ID);
    auto data_hist = it->second;
    std::vector<double> xs, xs_width;
    std::vector<double> ys, errs;
    for (int i = 1;
         i <= fBestFitSelectionHists[it->first]->GetNbinsX(); ++i) {
      ys.push_back(
          fBestFitSelectionHists[it->first]->GetBinContent(i));
      errs.push_back(
          sqrt(selection_cov.GetBinContent(bin_count+i, bin_count+i)));
      xs.push_back(data_hist->GetBinCenter(i));
      xs_width.push_back(data_hist->GetBinWidth(i)/2.);
    } 

    TGraphAsymmErrors throw_gr(data_hist->GetNbinsX(),
                               &xs[0], &ys[0], 
                               &xs_width[0], &xs_width[0], &errs[0], &errs[0]);

    throw_gr.SetTitle(fDataSet.GetSelectionName(selection_ID).c_str());
    throw_gr.GetXaxis()->SetTitle(data_hist->GetXaxis()->GetTitle());
    throw_gr.SetFillStyle(3144);
    throw_gr.SetFillColor(kRed);
    data_hist->Draw();
    throw_gr.Draw("same a2");
    data_hist->Draw("same e1");
    fOutputFile.cd("Throws");
    cThrow.Write();

    bin_count += data_hist->GetNbinsX();
  }

  bin_count = 0;
  for (auto it = truth_throw_hists.begin(); it != truth_throw_hists.end(); ++it) {
    int sample_ID = it->first;

    std::vector<double> xs, xs_width;
    for (size_t i = 0; i < sample_bins[it->first]; ++i) {
      xs.push_back(i + 0.5);
      xs_width.push_back(.5);
    }

    std::string name = "hNominal" + fSamples[sample_ID][0][0].GetName();
    TH1D temp_nominal(name.c_str(), "", xs.size(), 0, xs.size());
    std::string dummy_name = "hDummy" + fSamples[sample_ID][0][0].GetName();
    TH1D dummy(dummy_name.c_str(), "", xs.size(), 0, xs.size());
    if (fIsSignalSample[sample_ID]) {
      dummy.GetXaxis()->SetBinLabel(1, "Underflow");
      dummy.GetXaxis()->SetBinLabel(dummy.GetNbinsX(), "Overflow");
      for (int i = 2; i < dummy.GetNbinsX(); ++i) {
        std::ostringstream low_stream, high_stream;
        low_stream.precision(2);
        high_stream.precision(2);
        low_stream << std::fixed <<
                      fSamples[sample_ID][0][i-1].RangeLowEnd();
        high_stream
            << std::fixed <<
               fSamples[sample_ID][0][i-1].RangeHighEnd();
        std::string label = low_stream.str();
        label += " - ";
        label += high_stream.str();
        //std::string label
        //    = std::to_string(fSamples[sample_ID][0][i-1].RangeLowEnd());
        //label += " - ";
        //label += std::to_string(fSamples[sample_ID][0][i-1].RangeHighEnd());
        dummy.GetXaxis()->SetBinLabel(i, label.c_str());
      }
    }
    else {
      dummy.GetXaxis()->SetBinLabel(1, "");
    }


    std::vector<std::vector<ThinSliceSample>> & samples_vec_2D
        = fSamples[sample_ID];
    for (size_t i = 0; i < samples_vec_2D.size(); ++i) {
      for (size_t j = 0; j < samples_vec_2D[i].size(); ++j) {
        temp_nominal.AddBinContent(j+1, samples_vec_2D[i][j].GetNominalFlux());
      }
    }

    double max = -999.;
    for (size_t i = 0; i < sample_bins[it->first]; ++i) {
      if ((fBestFitTruthVals[sample_ID][i] + best_fit_errs[sample_ID][i]) > max)
        max = (fBestFitTruthVals[sample_ID][i] + best_fit_errs[sample_ID][i]);

      if (temp_nominal.GetBinContent(i+1) > max)
        max = temp_nominal.GetBinContent(i+1);
    }

    fOutputFile.cd("Throws");
    std::string canvas_name = "cTruthThrow" + fSamples[sample_ID][0][0].GetName();
    TCanvas cThrow(canvas_name.c_str(), "");
    cThrow.SetTicks();
    TGraphAsymmErrors throw_gr(xs.size(),
                                &xs[0], &fBestFitTruthVals[it->first][0], 
                                &xs_width[0], &xs_width[0],
                                &best_fit_errs[it->first][0],
                                &best_fit_errs[it->first][0]);
    throw_gr.SetFillStyle(3144);
    throw_gr.SetFillColor(kRed);
    throw_gr.SetMinimum(0.);
    throw_gr.SetMaximum(1.5*max);
    dummy.SetMaximum(1.5*max);
    dummy.Draw();
    //throw_gr.Draw("a2");
    throw_gr.Draw("same 2");
    throw_gr.Draw("p");

    temp_nominal.SetMarkerColor(kBlue);
    temp_nominal.SetMarkerStyle(20);
    temp_nominal.Draw("same p");

    TLegend leg;
    leg.AddEntry(&throw_gr, "Throws", "lpf");
    leg.AddEntry(&temp_nominal, "Nominal", "p");

    if (fDoFakeData) {
      name = "hVaried" + fSamples[sample_ID][0][0].GetName();
      TH1D * temp_varied = (TH1D*)temp_nominal.Clone(name.c_str());
      for (size_t i = 0; i < xs.size(); ++i) {
        temp_varied->SetBinContent(
            i+1, temp_varied->GetBinContent(i+1)*fFakeDataScales[sample_ID][i]);
      }
      temp_varied->SetMarkerColor(kBlack);
      temp_varied->SetMarkerStyle(20);
      temp_varied->Draw("same p");
      leg.AddEntry(temp_varied, "Fake Data", "p");
    }

    leg.Draw();
    cThrow.Write();

    bin_count += xs.size();
  }

  for (auto it = best_fit_xsec_truth.begin(); it != best_fit_xsec_truth.end();
       ++it) {
    int sample_ID = it->first;

    std::vector<double> xs, xs_width;
    for (size_t i = 0; i < xsec_bins[sample_ID]; ++i) {
      xs.push_back(i + 0.5);
      xs_width.push_back(.5);
    }

    fOutputFile.cd("Throws");
    std::string canvas_name = "cXSecThrow" + fSamples[sample_ID][0][0].GetName();
    TCanvas cThrow(canvas_name.c_str(), "");
    cThrow.SetTicks();

    std::string dummy_name = "hDummyXSec" + fSamples[sample_ID][0][0].GetName();
    TH1D dummy(dummy_name.c_str(), "", xs.size(), 0, xs.size());
    for (int i = 1; i <= dummy.GetNbinsX(); ++i) {
      std::ostringstream low_stream, high_stream;
      low_stream.precision(2);
      high_stream.precision(2);
      low_stream << std::fixed <<
                    fSamples[sample_ID][0][i].RangeLowEnd();
      high_stream
          << std::fixed <<
             fSamples[sample_ID][0][i].RangeHighEnd();
      std::string label = low_stream.str();
      label += " - ";
      label += high_stream.str();
      //std::string label
      //    = std::to_string(fSamples[sample_ID][0][i-1].RangeLowEnd());
      //label += " - ";
      //label += std::to_string(fSamples[sample_ID][0][i-1].RangeHighEnd());
      dummy.GetXaxis()->SetBinLabel(i, label.c_str());
    }

    TGraphAsymmErrors throw_gr(xs.size(),
                               &xs[0], &best_fit_xsec_truth[it->first][0], 
                               &xs_width[0], &xs_width[0],
                               &best_fit_xsec_errs[it->first][0],
                               &best_fit_xsec_errs[it->first][0]);
    double max = -999.;
    for (size_t i = 0; i < best_fit_xsec_truth[it->first].size(); ++i) {
      if ((best_fit_xsec_truth[it->first][i] +
           best_fit_xsec_errs[it->first][0]) > max) {
        max = (best_fit_xsec_truth[it->first][i] +
               best_fit_xsec_errs[it->first][0]);
      }
    }
    throw_gr.SetFillStyle(3144);
    throw_gr.SetFillColor(kRed);
    throw_gr.SetMinimum(0.);
    dummy.SetMaximum(1.5*max);
    dummy.SetMinimum(0.);
    dummy.Draw();
    dummy.GetXaxis()->SetTitle("Kinetic Energy (MeV)");
    dummy.GetYaxis()->SetTitle("#sigma (mb)");
    throw_gr.Draw("same 2");
    throw_gr.Draw("p");

    std::vector<double> nominal_xsec_vals;
    for (size_t i = 0; i < xs.size(); ++i) {
      nominal_xsec_vals.push_back(
          fNominalXSecs[it->first]->GetBinContent(i+1));
    }
    TGraph nominal_gr(xs.size(), &xs[0], &nominal_xsec_vals[0]);
    nominal_gr.SetMarkerColor(kBlue);
    nominal_gr.SetMarkerStyle(20);
    nominal_gr.Draw("same p");

    TLegend leg;
    leg.AddEntry(&throw_gr, "Throws", "lpf");
    leg.AddEntry(&nominal_gr, "Nominal", "p");

    if (fDoFakeData) {
      //std::cout << "Plotting fake data" << std::endl;
      std::vector<double> fake_xsec_vals;
      //std::cout << xs.size() << std::endl;
      for (size_t i = 0; i < xs.size(); ++i) {

        //std::cout << fFakeDataXSecs[it->first] << std::endl;
        //std::cout << "Testing" << std::endl;
        //std::cout << fFakeDataXSecs[it->first]->GetNbinsX() << std::endl;
        fake_xsec_vals.push_back(
            fFakeDataXSecs[it->first]->GetBinContent(i+1));
        //std::cout << "Adding " << fake_xsec_vals.back() << std::endl;
      }
      TGraph * fake_gr = new TGraph(xs.size(), &xs[0], &fake_xsec_vals[0]);
      fake_gr->Write();
      fake_gr->SetMarkerColor(kBlack);
      fake_gr->SetMarkerStyle(20);
      fake_gr->Draw("same p");
      leg.AddEntry(fake_gr, "Fake Data", "p");
    }

    leg.Draw("same");
    cThrow.Write();
  }
}

void protoana::PDSPThinSliceFitter::GetCurrentTruthHists(
    std::map<int, std::vector<TH1*>> & throw_hists,
    std::map<int, std::vector<TH1*>> & throw_inc_hists,
    std::map<int, std::vector<TH1*>> & throw_xsec_hists) {
  //Loop over the samples
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    //Get the number of bins from the first entry of the beam energy bins
    std::vector<std::vector<ThinSliceSample>> & samples_vec_2D = it->second;
    size_t nBins = samples_vec_2D[0].size();
    std::string name = it->second[0][0].GetName() + "Throw" +
                       std::to_string(throw_hists[it->first].size());
    TH1D * temp_hist = new TH1D(name.c_str(), "", nBins, 0, nBins); 
    for (size_t i = 0; i < samples_vec_2D.size(); ++i) {
      std::vector<ThinSliceSample> & samples_vec = samples_vec_2D[i];
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        temp_hist->AddBinContent(j+1, samples_vec[j].GetVariedFlux());
      }
    }
    throw_hists[it->first].push_back(temp_hist);
  }

  for (auto it = throw_inc_hists.begin(); it != throw_inc_hists.end(); ++it) {
    int s = it->first;
    auto & samples_vec_2D = fSamples[s];
    const std::vector<double> & bins = fSignalBins[s];
    std::string name = samples_vec_2D[0][0].GetName();
    name += "IncidentThrow" +
             std::to_string(throw_inc_hists[it->first].size());
    TH1D * temp_inc_hist = new TH1D(name.c_str(), "", bins.size() - 1, &bins[0]); 
    

    name = samples_vec_2D[0][0].GetName();
    name += "XSecThrow" +
             std::to_string(throw_inc_hists[it->first].size());
    TH1D * temp_xsec_hist = new TH1D(name.c_str(), "", bins.size() - 1,
                                     &bins[0]);
    for (auto i_s : fIncidentSamples) {
      auto & incident_vec_2D = fSamples[i_s];
      for (size_t i = 0; i < incident_vec_2D.size(); ++i) {
        for (size_t j = 0; j < incident_vec_2D[i].size(); ++j) {
          if (fSliceMethod == "E") {
            incident_vec_2D[i][j].FillESliceHist(*temp_inc_hist);
          }
          else {
            incident_vec_2D[i][j].FillHistFromIncidentEnergies(*temp_inc_hist);
          }
        }
      }
    }
    throw_inc_hists[s].push_back(temp_inc_hist);

    for (int i = 1; i <= temp_xsec_hist->GetNbinsX(); ++i) {
      temp_xsec_hist->SetBinContent(
          i, throw_hists[s].back()->GetBinContent(i+1));
    }
    temp_xsec_hist->Divide(temp_inc_hist);
    throw_xsec_hists[s].push_back(temp_xsec_hist);
  }  
}

void protoana::PDSPThinSliceFitter::BuildFakeDataXSecs() {
  //First, set all samples to fake data scales
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        it->second[i][j].SetFactorAndScale(fFakeDataScales[it->first][j]);
      }
    }
  }

  for (auto s : fMeasurementSamples) {
    auto & samples_vec_2D = fSamples[s];
    std::vector<double> & bins = fSignalBins[s];
    std::string xsec_name = "FakeData" +
                             samples_vec_2D[0][1].GetName() + "XSec";
    TH1D * temp_xsec = new TH1D(xsec_name.c_str(), "",
                                fSignalBins[s].size() - 1, &bins[0]);
    for (size_t i = 0; i < samples_vec_2D.size(); ++i) {
      for (size_t j = 1; j < samples_vec_2D[i].size() - 1; ++j) {
        temp_xsec->AddBinContent(j, samples_vec_2D[i][j].GetVariedFlux());
      }
    }

    std::string inc_name = "FakeData" +
                             samples_vec_2D[0][1].GetName() + "Inc";
    TH1D * temp_inc = new TH1D(inc_name.c_str(), "", fSignalBins[s].size() - 1,
                               &fSignalBins[s][0]);
    for (size_t i = 0; i < fIncidentSamples.size(); ++i) {
      auto & vec_2D = fSamples[fIncidentSamples[i]];
      for (size_t j = 0; j < vec_2D.size(); ++j) {
        auto & samples_vec = vec_2D[j];
        for (size_t k = 0; k < samples_vec.size(); ++k) {
          if (fSliceMethod == "E") {
            samples_vec[k].FillESliceHist(*temp_inc);
          }
          else {
            samples_vec[k].FillHistFromIncidentEnergies(*temp_inc);
          }
        }
      }
    }

    temp_xsec->Divide(temp_inc);

    if (fSliceMethod == "E") {
      for (int i = 1; i <= temp_xsec->GetNbinsX(); ++i) {
        temp_xsec->SetBinContent(i, -1.*log(1. - temp_xsec->GetBinContent(i)));
      }
      temp_xsec->Scale(1.E27*39.948/(1.4 * 6.022E23));
      for (int i = 1; i <= temp_xsec->GetNbinsX(); ++i) {

        double bethe_val = BetheBloch(temp_xsec->GetBinCenter(i), 139.57); 

        temp_xsec->SetBinContent(i, (bethe_val*
                                     temp_xsec->GetBinContent(i)/
                                     temp_xsec->GetBinWidth(i)));
      }
    }
    else if (fSliceMethod == "Alt") {
      for (int i = 1; i <= temp_xsec->GetNbinsX(); ++i) {
        temp_xsec->SetBinContent(i, -1.*log(1. - temp_xsec->GetBinContent(i)));
      }
      temp_xsec->Scale(1.E27*39.948/(fPitch * 1.4 * 6.022E23));
    }
    else {
      temp_xsec->Scale(1.E27/ (fPitch * 1.4 * 6.022E23 / 39.948 ));
    }
    fFakeDataXSecs[s] = temp_xsec;
    fFakeDataIncs[s] = temp_inc;

    fFakeDataXSecs[s]->SetDirectory(0);
    fFakeDataIncs[s]->SetDirectory(0);

    //std::cout << "Xsec vals " << fFakeDataXSecs[s] << std::endl;
    //for (int i = 1; i <= temp_xsec->GetNbinsX(); ++i) {
    //  std::cout << temp_xsec->GetBinContent(i) << " ";
    //}
    //std::cout << std::endl;
  }

  //Now, reset all samples
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        it->second[i][j].ResetFactor();
      }
    }
  }
}

