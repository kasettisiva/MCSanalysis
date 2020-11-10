// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
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
#include "TList.h"
#include "TDirectory.h"
#include "TH2D.h"
#include "TGraph.h"

#include "PDSPThinSliceFitter.h"

protoana::PDSPThinSliceFitter::PDSPThinSliceFitter(std::string fcl_file,
                                                   std::string output_file)
    : fOutputFile(output_file.c_str(), "RECREATE") {
  Configure(fcl_file);
  gStyle->SetOptStat(0);
  gROOT->SetBatch(1);
}

protoana::PDSPThinSliceFitter::~PDSPThinSliceFitter() {
  fOutputFile.Close();
}

void protoana::PDSPThinSliceFitter::MakeMinimizer() {
  fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
    (ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  fMinimizer->SetMaxFunctionCalls(fMaxCalls);
  fMinimizer->SetTolerance(fTolerance);

  size_t total_parameters = fTotalSignalParameters + fTotalFluxParameters;
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

  parsHist.SetMarkerColor(kBlue);
  parsHist.SetMarkerStyle(20);
  fOutputFile.cd();
  parsHist.Write();

}

void protoana::PDSPThinSliceFitter::InitializeMCSamples() {
  for (size_t i = 0; i < fSampleSets.size(); ++i) {
    fhicl::ParameterSet sample_set = fSampleSets[i];
    std::string sample_name = sample_set.get<std::string>("Name");
    int sample_ID = sample_set.get<int>("ID");

    if (fSamples.find(sample_ID) == fSamples.end()) {
      fSamples[sample_ID] = std::vector<ThinSliceSample>();
      fFluxesBySample[sample_ID] = std::vector<double>();
    }
    fIsSignalSample[sample_ID] = sample_set.get<bool>("IsSignal");

    int flux_type = sample_set.get<int>("FluxType");

    if (sample_set.get<bool>("IsSignal")) {
      std::vector<double> bins
          = sample_set.get<std::vector<double>>("SignalBins");
      fSignalParameters[sample_ID] = std::vector<double>();
      fSignalParameterNames[sample_ID] = std::vector<std::string>();
      for (size_t j = 1; j < bins.size(); ++j) {
        ThinSliceSample sample(sample_name, flux_type, fSelectionIDs,
                               fIncidentRecoBins,
                               fSelectedRecoBins, true,
                               {bins[j-1], bins[j]});
        fSamples[sample_ID].push_back(sample);
        fSignalParameters[sample_ID].push_back(1.);

        std::string par_name = "par_" + sample_name + "_" +
                               protoana::PreciseToString(bins[j-1]) + "_" +
                               protoana::PreciseToString(bins[j]);
        fSignalParameterNames[sample_ID].push_back(par_name);
        ++fTotalSignalParameters;
        fFluxesBySample[sample_ID].push_back(0.);
      }
    }
    else {
      ThinSliceSample sample(sample_name, flux_type, fSelectionIDs,
                             fIncidentRecoBins,
                             fSelectedRecoBins);
      fSamples[sample_ID].push_back(sample);
      fFluxesBySample[sample_ID].push_back(0.);
    }
  }

  MakeMinimizer();
}

void protoana::PDSPThinSliceFitter::BuildMCSamples() {
  //Open the MC file and set branches
  TFile fMCFile(fMCFileName.c_str(), "OPEN");
  fMCTree = (TTree*)fMCFile.Get(fTreeName.c_str());

  int true_sample, selection_ID; 
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  fMCTree->SetBranchAddress("new_interaction_topology", &true_sample);
  fMCTree->SetBranchAddress("selection_ID", &selection_ID);
  fMCTree->SetBranchAddress("true_beam_interactingEnergy",
                            &true_beam_interactingEnergy);
  fMCTree->SetBranchAddress("reco_beam_interactingEnergy",
                            &reco_beam_interactingEnergy);
  fMCTree->SetBranchAddress("reco_beam_incidentEnergies",
                            &reco_beam_incidentEnergies);

  for (int i = 0; i < fMCTree->GetEntries(); ++i) {
    fMCTree->GetEntry(i);

    if (fSamples.find(true_sample) == fSamples.end())
      continue;

    std::vector<ThinSliceSample> & samples = fSamples[true_sample];
    bool is_signal = fIsSignalSample[true_sample];

    if (!is_signal) {
      ThinSliceSample & sample = samples.at(0);
      int flux_type = sample.GetFluxType();
      fNominalFluxes[flux_type] += 1.;
      fFluxesBySample[true_sample][0] += 1.;
      sample.AddFlux();
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
          int flux_type = sample.GetFluxType();
          fNominalFluxes[flux_type] += 1.;
          fFluxesBySample[true_sample][j] += 1.;
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
  fMCFile.Close();
}

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
      it->second[i].SetDataMCScale(fMCDataScale);
    }
  }
}

void protoana::PDSPThinSliceFitter::BuildDataHists() {

  //Create the data hists
  fIncidentDataHist = TH1D("Data_incident_hist",
                           "Data;Reconstructed KE (MeV)",
                           fIncidentRecoBins.size() - 1,
                           &fIncidentRecoBins[0]);
  for (auto it = fSelectionIDs.begin(); it != fSelectionIDs.end(); ++it) {
    std::string sel_name = "Data_selected_" + it->second + "_hist";
    fSelectedDataHists[it->first] = TH1D(sel_name.c_str(),
                                         "Data;Reconstructed KE (MeV)",
                                         fSelectedRecoBins.size() - 1,
                                         &fSelectedRecoBins[0]);
  }

  //Open the Data file and set branches
  TFile fDataFile(fDataFileName.c_str(), "OPEN");
  fDataTree = (TTree*)fDataFile.Get(fTreeName.c_str());

  int selection_ID; 
  double reco_beam_interactingEnergy;
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  fDataTree->SetBranchAddress("selection_ID", &selection_ID);
  fDataTree->SetBranchAddress("reco_beam_interactingEnergy",
                            &reco_beam_interactingEnergy);
  fDataTree->SetBranchAddress("reco_beam_incidentEnergies",
                            &reco_beam_incidentEnergies);

  fDataFlux = fDataTree->GetEntries();

  for (int i = 0; i < fDataTree->GetEntries(); ++i) {
    fDataTree->GetEntry(i);
    if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        fIncidentDataHist.Fill((*reco_beam_incidentEnergies)[j]);
      }
      if (fSelectionIDs.find(selection_ID) != fSelectionIDs.end()) {
        fSelectedDataHists[selection_ID].Fill(reco_beam_interactingEnergy);
      }
    }
  }
  fOutputFile.cd();
  TDirectory * out = (TDirectory *)fOutputFile.mkdir("Data");
  out->cd();
  fIncidentDataHist.Write();
  for (auto it = fSelectedDataHists.begin(); it != fSelectedDataHists.end();
       ++it) {
    it->second.Write();
  }
  MakeRebinnedDataHists();
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

void protoana::PDSPThinSliceFitter::BuildAndSaveStacks(bool post_fit) {

  THStack * incident_stack = new THStack(
      (post_fit ? "PostFitIncidentStack" : "NominalIncidentStack"), "");
  //Incident
  size_t iColor = 0;
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      it->second.at(i).RefillRebinnedHists();

      TH1D & inc_hist = (fPlotRebinned ?
                         it->second.at(i).GetRebinnedIncidentHist() :
                         it->second.at(i).GetIncidentHist());

      std::pair<int, int> color_fill = GetColorAndStyle(iColor);
      inc_hist.SetFillColor(color_fill.first);
      inc_hist.SetFillStyle(color_fill.second);
      inc_hist.SetLineColor(kBlack);
      incident_stack->Add(&inc_hist);
      ++iColor;
    }
  }

  fOutputFile.cd();
  incident_stack->Write();
  if (post_fit) {
    fPostFitIncidentMCStack = incident_stack;
  }
  else {
    fNominalIncidentMCStack = incident_stack;
  }

  for (auto it1 = fSelectionIDs.begin(); it1 != fSelectionIDs.end(); ++it1) {
    std::string s_name = (post_fit ? "PostFit" : "Nominal") + it1->second +
                         "Stack";
    THStack * selected_stack = new THStack(s_name.c_str(), "");
    iColor = 0;
    for (auto it2 = fSamples.begin(); it2 != fSamples.end(); ++it2) {
      for (size_t i = 0; i < it2->second.size(); ++i) {
        TH1D & sel_hist = (fPlotRebinned ?
                           it2->second.at(i).GetRebinnedSelectionHist(
                               it1->first) :
                           it2->second.at(i).GetSelectionHist(it1->first));

        std::pair<int, int> color_fill = GetColorAndStyle(iColor);
        sel_hist.SetFillColor(color_fill.first);
        sel_hist.SetFillStyle(color_fill.second);
        sel_hist.SetLineColor(kBlack);
        selected_stack->Add(&sel_hist);
        ++iColor;
      }
    }
    selected_stack->Write();
    if (post_fit) {
      fPostFitSelectedMCStacks[it1->first] = selected_stack;
    }
    else {
      fNominalSelectedMCStacks[it1->first] = selected_stack;
    }
  }
}

void protoana::PDSPThinSliceFitter::CompareDataMC(bool post_fit) {
  //Incident
  TH1D & inc_data_hist = (fPlotRebinned ? fRebinnedIncidentDataHist :
                                          fIncidentDataHist);
  inc_data_hist.SetLineColor(kBlack);
  inc_data_hist.SetMarkerColor(kBlack);
  inc_data_hist.SetMarkerStyle(20);

  TCanvas cIncident((post_fit ? "cPostFitIncident" : "cNominalIncident"), "");
  cIncident.SetTicks();
  THStack * incident_stack = (post_fit ? fPostFitIncidentMCStack :
                                         fNominalIncidentMCStack);
  incident_stack->Draw("hist");
  incident_stack->SetTitle("Incident Sample;Reconstructed KE (MeV)");
  incident_stack->GetHistogram()->SetTitleSize(.04, "X");
  incident_stack->Draw("hist");
  inc_data_hist.Draw("e1 same");
  fOutputFile.cd();
  cIncident.Write();

  //Get the full incident hist from stack
  TList * l = (TList*)incident_stack->GetHists();
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

  //Selected hists
  for (auto it = fSelectionIDs.begin(); it != fSelectionIDs.end(); ++it) {
    int selection_ID = it->first;
    TH1D & sel_data_hist = (fPlotRebinned ?
                            fRebinnedSelectedDataHists[selection_ID] :
                            fSelectedDataHists[selection_ID]);
    sel_data_hist.SetLineColor(kBlack);
    sel_data_hist.SetMarkerColor(kBlack);
    sel_data_hist.SetMarkerStyle(20);   

    std::string canvas_name = "c";
    canvas_name += (post_fit ? "PostFit" : "Nominal") + it->second;
    TCanvas cSelection(canvas_name.c_str(), "");
    cSelection.SetTicks();
    THStack * selected_stack
        = (post_fit ? fPostFitSelectedMCStacks[selection_ID] :
                      fNominalSelectedMCStacks[selection_ID]);
    selected_stack->Draw("hist");
    std::string title = it->second + ";Reconstructed KE (MeV)";
    selected_stack->SetTitle(title.c_str());
    selected_stack->GetHistogram()->SetTitleSize(.04, "X");
    selected_stack->Draw("hist");
    sel_data_hist.Draw("e1 same");

    fOutputFile.cd();
    cSelection.Write();

    //Get the full incident hist from stack
    TList * l = (TList*)selected_stack->GetHists();
    TH1D * hMC = (TH1D*)l->At(0)->Clone();
    for (int i = 1; i < l->GetSize(); ++i) {
      hMC->Add((TH1D*)l->At(i));
    }

    std::string ratio_name = it->second + "Ratio" + (post_fit ? "PostFit" :
                                                                "Nominal");
    TH1D * hRatio
        = (TH1D*)sel_data_hist.Clone(ratio_name.c_str());
    hRatio->Divide(hMC);
    hRatio->Write(); 
  }
}

void protoana::PDSPThinSliceFitter::RunFitAndSave() {
  DefineFitFunction();
  fMinimizer->SetFunction(fFitFunction);
  int fit_status = fMinimizer->Minimize();

  if (!fit_status) {
    std::cout << "Failed to find minimum: " << std::endl;
  }
  else {
    std::cout << "Found minimimum: " << std::endl;
    size_t total_parameters = fTotalSignalParameters + fTotalFluxParameters;
    for (size_t i = 0; i < total_parameters; ++i) {
      std::cout << fMinimizer->VariableName(i) << " " << fMinimizer->X()[i] <<
                   std::endl;
    }

    TH2D covHist("covHist", "", total_parameters, 0,
                 total_parameters, total_parameters, 0,
                 total_parameters);
    TH2D corrHist("corrHist", "", total_parameters, 0,
                  total_parameters, total_parameters, 0,
                  total_parameters);

    TH1D parsHist("postFitPars", "", total_parameters, 0,
                  total_parameters);

    int n_par = 0;
    for (auto it = fSignalParameters.begin();
         it != fSignalParameters.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        parsHist.SetBinContent(n_par+1, fMinimizer->X()[n_par]);
        parsHist.SetBinError(n_par+1,
                             sqrt(fMinimizer->CovMatrix(n_par, n_par)));
        parsHist.GetXaxis()->SetBinLabel(
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

    for (size_t i = 0; i < (fTotalSignalParameters + fTotalFluxParameters);
         ++i) {
      for (size_t j = 0; j < (fTotalSignalParameters + fTotalFluxParameters);
           ++j) {
        covHist.SetBinContent(i+1, j+1, fMinimizer->CovMatrix(i, j));
        corrHist.SetBinContent(i+1, j+1, fMinimizer->Correlation(i, j));
      }
    }

    fOutputFile.cd();
    parsHist.SetFillColor(kRed);
    parsHist.SetFillStyle(3144);
    parsHist.SetMarkerColor(kBlack);
    parsHist.SetMarkerStyle(20);
    parsHist.Write();

    covHist.Write();
    corrHist.Write();

    //save post fit stacks
    BuildAndSaveStacks(true);
    CompareDataMC(true);
    ParameterScans();
  }
}

void protoana::PDSPThinSliceFitter::ParameterScans() {
  fOutputFile.cd();
  TDirectory * out = (TDirectory *)fOutputFile.mkdir("Scans");
  out->cd();

  std::cout << "Scanning parameters" << std::endl;
  size_t total_parameters = fTotalSignalParameters + fTotalFluxParameters;

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

        //Then go through and vary the fluxes by the flux parameters
        //Remember: 1 less parameter than number of flux types
        std::map<int, double> fluxes, scaled_fluxes;
        double total_nominal_flux = 0.;
        double total_varied_flux = 0.;
        for (auto it = fNominalFluxes.begin();
             it != fNominalFluxes.end(); ++it) {
          int flux_id = it->first;
          fluxes[flux_id]
              = it->second * ((fFluxParameters.find(flux_id) != fFluxParameters.end()) ?
                              fFluxParameters[flux_id] : 1.);

          //Determine total of fluxes, both nominal and varied, in order to 
          //rescale later on
          total_nominal_flux += it->second;
          total_varied_flux += fluxes[flux_id];

          scaled_fluxes[flux_id] = 0.;
        }

        for (auto it = fluxes.begin(); it != fluxes.end(); ++it) {
          it->second *= (total_nominal_flux/total_varied_flux);
        }

        //Go through and determine the scaled flux
        for (auto it = fFluxesBySample.begin();
             it != fFluxesBySample.end(); ++it) {
          int sample_ID = it->first;
          std::vector<ThinSliceSample> & samples
              = fSamples[sample_ID];

          for (size_t i = 0; i < samples.size(); ++i) {
            ThinSliceSample & sample = samples.at(i);
            int flux_type = sample.GetFluxType();
            
            if (fSignalParameters.find(sample_ID) == fSignalParameters.end()) {
              scaled_fluxes[flux_type]
                  += fFluxesBySample[sample_ID][i];
            }
            else {
              //Scale the flux for this signal channel
              scaled_fluxes[flux_type]
                  += (fFluxesBySample[sample_ID][i] *
                      fSignalParameters[sample_ID][i]);
            }
          }
        }

        //Use the fluxes from above to figure out the correct
        //scaling factor for the signal bins, and use this to 
        //scale the hists in the signal samples accordingly
        //along with the parameter
        for (auto it = fSignalParameters.begin();
             it != fSignalParameters.end(); ++it) {
          int sample_ID = it->first;
          std::vector<ThinSliceSample> & samples
              = fSamples[sample_ID];
          for (size_t i = 0; i < it->second.size(); ++i) {
            ThinSliceSample & sample = samples.at(i);
            int flux_type = sample.GetFluxType();
            double flux_factor = fluxes[flux_type]/scaled_fluxes[flux_type];
            sample.SetFactorAndScale(flux_factor*it->second.at(i));
          }
        }

        std::pair<double, size_t> chi2_points = CalculateChi2();
        return (chi2_points.first/chi2_points.second);
      },
      fTotalSignalParameters + fTotalFluxParameters);

}

std::pair<double, size_t> protoana::PDSPThinSliceFitter::CalculateChi2() {
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

  //fFluxTypes = pset.get<std::vector<int>>("FluxTypes");
  temp_vec = pset.get<std::vector<std::pair<int, std::string>>>(
      "FluxTypes");
  fFluxTypes = std::map<int, std::string>(temp_vec.begin(), temp_vec.end());
  for (size_t i = 0; i < temp_vec.size() - 1; ++i) {
    fFluxParameterNames[temp_vec[i].first] = "par_" + temp_vec[i].second +
                                             "_flux";
    fFluxParameters[temp_vec[i].first] = 1.;
  }
  fTotalFluxParameters = temp_vec.size() - 1;

  fMaxCalls = pset.get<int>("MaxCalls");
  fNScanSteps = pset.get<unsigned int>("NScanSteps") + 1;
  fTolerance = pset.get<double>("Tolerance");
  fLowerLimit = pset.get<double>("LowerLimit");
  fUpperLimit = pset.get<double>("UpperLimit");
  fReducedIncidentChi2 = pset.get<bool>("ReducedIncidentChi2");
  fPlotStyle = pset.get<std::vector<std::pair<int,int>>>("PlotStyle");
  fPlotRebinned = pset.get<bool>("PlotRebinned");
  fRandomStart = pset.get<bool>("RandomStart");
}

std::pair<int, int> protoana::PDSPThinSliceFitter::GetColorAndStyle(size_t i) {
  return {fPlotStyle[i % fPlotStyle.size()].first,
          (i < fPlotStyle.size() ? 1001 : 3244)};
  //return (i < fPlotStyle.size() ? fPlotStyle[i].first : i);
}

void protoana::PDSPThinSliceFitter::MakeRebinnedDataHists() {
  fRebinnedIncidentDataHist = TH1D("Data_incident_hist_rebinned",
                                   "Data;Reconstructed KE (MeV)",
                                   fIncidentRecoBins.size() - 1, 0,
                                   fIncidentRecoBins.size() - 1);
  for (int i = 1; i <= fIncidentDataHist.GetNbinsX(); ++i) {
    fRebinnedIncidentDataHist.SetBinContent(
        i, fIncidentDataHist.GetBinContent(i));
    double low_edge = fIncidentDataHist.GetXaxis()->GetBinLowEdge(i);
    double up_edge = fIncidentDataHist.GetXaxis()->GetBinUpEdge(i);
    std::string bin_label = (low_edge < 0. ?
                             "< 0." :
                             (PreciseToString(low_edge, 0) + " - " + 
                              PreciseToString(up_edge, 0)));
    fRebinnedIncidentDataHist.GetXaxis()->SetBinLabel(i, bin_label.c_str());
  }
  for (auto it = fSelectionIDs.begin(); it != fSelectionIDs.end(); ++it) {
    TH1D & sel_hist = fSelectedDataHists[it->first];
    std::string sel_name = "Data_selected_" + it->second + "_hist_rebinned";
    fRebinnedSelectedDataHists[it->first] = TH1D(sel_name.c_str(),
                                                 "Data;Reconstructed KE (MeV)",
                                                 fSelectedRecoBins.size() - 1, 0,
                                                 fSelectedRecoBins.size() - 1);
    TH1D & sel_hist_rebinned = fRebinnedSelectedDataHists[it->first];
    for (int i = 1; i <= sel_hist.GetNbinsX(); ++i) {
      sel_hist_rebinned.SetBinContent(
          i, sel_hist.GetBinContent(i));
      double low_edge = sel_hist.GetXaxis()->GetBinLowEdge(i);
      double up_edge = sel_hist.GetXaxis()->GetBinUpEdge(i);
      std::string bin_label = (low_edge < 0. ?
                               "< 0." :
                               (PreciseToString(low_edge, 0) + " - " + 
                                PreciseToString(up_edge, 0)));
      sel_hist_rebinned.GetXaxis()->SetBinLabel(i, bin_label.c_str());
    }
  }
}
