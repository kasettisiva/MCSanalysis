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
#include "TLine.h"
#include "TList.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraph.h"

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

//Good
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
        ThinSliceSample sample(sample_name, flux_type, fSelectionSets,
                               fIncidentRecoBins,
                               true,
                               {bins[j-1], bins[j]});
        fSamples[sample_ID].push_back(sample);
        fSignalParameters[sample_ID].push_back(1.);

        std::string par_name = "par_" + sample_name + "_" +
                               PreciseToString(bins[j-1]) + "_" +
                               PreciseToString(bins[j]);
        fSignalParameterNames[sample_ID].push_back(par_name);
        ++fTotalSignalParameters;
        fFluxesBySample[sample_ID].push_back(0.);
      }
    }
    else {
      ThinSliceSample sample(sample_name, flux_type, fSelectionSets,
                             fIncidentRecoBins);
      fSamples[sample_ID].push_back(sample);
      fFluxesBySample[sample_ID].push_back(0.);
      if (fFluxParameters.find(flux_type) != fFluxParameters.end()) {
        fFluxParsToSamples[flux_type].push_back(sample_ID);
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
  fThinSliceDriver->BuildMCSamples(fMCTree, fSamples, fIsSignalSample,
                                   fNominalFluxes, fFluxesBySample);
  fMCFile.Close();
}

//Good
void protoana::PDSPThinSliceFitter::ScaleMCToData() {
  double total_nominal = 0.;
  for (auto it = fNominalFluxes.begin(); it != fNominalFluxes.end(); ++it) {
    total_nominal += it->second;
  }

  /*
  std::cout << "Fluxes: " << std::endl;
  for (auto it = fFluxesBySample.begin(); it != fFluxesBySample.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      std::cout << "\t" << it->first << " " << i << " " << it->second[i] <<
                   std::endl;
    }
  }

  std::cout << "Total MC: " << total_nominal << std::endl;
  std::cout << "Total Data: " << fDataFlux << std::endl;
  */

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
  fDataSet = ThinSliceDataSet(fIncidentRecoBins, fSelectionSets);

  //Open the Data file and set branches
  //fDataFlux = fDataTree->GetEntries();

  if (!fDoFakeData) {
    TFile fDataFile(fDataFileName.c_str(), "OPEN");
    fDataTree = (TTree*)fDataFile.Get(fTreeName.c_str());
    fThinSliceDriver->BuildDataHists(fDataTree, fDataSet, fDataFlux);
  }
  else {
    TFile fMCFile(fMCFileName.c_str(), "OPEN");
    fMCTree = (TTree*)fMCFile.Get(fTreeName.c_str());
    fThinSliceDriver->BuildFakeData(fMCTree, fSamples, fDataSet, fDataFlux);
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
      it->second.at(i).GetIncidentHist().Write();
      const std::map<int, TH1 *> & hists = it->second.at(i).GetSelectionHists();
      for (auto it2 = hists.begin(); it2 != hists.end(); ++it2) {
        it2->second->Write();
      }
    }
  }
}

void protoana::PDSPThinSliceFitter::CompareDataMC(bool post_fit) {
  fThinSliceDriver->CompareDataMC(fSamples, fDataSet, fOutputFile,
                                  fPlotStyle, fPlotRebinned, post_fit);
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

    //Drawing pre + post fit pars
    TCanvas cPars("cParameters", "");
    cPars.SetTicks();
    parsHist.Draw("e2");
    TH1D * preFitParsHist = (TH1D*)fOutputFile.Get("preFitPars");
    preFitParsHist->Draw("p same");
    TLine l(0., 1., parsHist.GetXaxis()->GetBinUpEdge(parsHist.GetNbinsX()), 1.);
    l.SetLineColor(kBlack);
    l.Draw("same");
    TLegend leg(.15, .65, .45, .85);
    leg.AddEntry(preFitParsHist, "Pre-Fit", "p");
    leg.AddEntry(&parsHist, "Post-Fit", "lp");
    leg.Draw();
    cPars.Write();

    covHist.Write();
    corrHist.Write();

    //save post fit stacks
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

        double varied_flux = 0.;
        double nominal_flux = 0.;

        //Go through and determine the scaled flux
        for (auto it = fFluxesBySample.begin();
             it != fFluxesBySample.end(); ++it) {
          int sample_ID = it->first;
          std::vector<ThinSliceSample> & samples
              = fSamples[sample_ID];

          for (size_t i = 0; i < samples.size(); ++i) {
            ThinSliceSample & sample = samples.at(i);
            int flux_type = sample.GetFluxType();
            nominal_flux += fFluxesBySample[sample_ID][i];
            
            if (fSignalParameters.find(sample_ID) != fSignalParameters.end()) {
              //Scale the flux for this signal channel
              varied_flux
                  += (fFluxesBySample[sample_ID][i] *
                      fSignalParameters[sample_ID][i]);
            }
            else if (fFluxParameters.find(flux_type) != fFluxParameters.end()) {
              varied_flux
                  += (fFluxesBySample[sample_ID][i] *
                      fFluxParameters[flux_type]);
            }
            else {
              varied_flux
                  += fFluxesBySample[sample_ID][i];
            }
          }
        }

        double flux_factor = nominal_flux / varied_flux;

        for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
          //Check if it's a signal
          int sample_ID = it->first;
          std::vector<ThinSliceSample> & samples = it->second;
          if (fSignalParameters.find(sample_ID) != fSignalParameters.end()) {
            for (size_t i = 0; i < it->second.size(); ++i) {
              double val = fSignalParameters[sample_ID][i];
              ThinSliceSample & sample = samples.at(i);
              sample.SetFactorAndScale(flux_factor*val);
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
                for (size_t i = 0; i < it->second.size(); ++i) {
                  ThinSliceSample & sample = samples.at(i);
                  sample.SetFactorAndScale(flux_factor*val);
                }
                flux_sample = true;
                break;
              }
            }
            if (!flux_sample) {
              for (size_t i = 0; i < it->second.size(); ++i) {
                ThinSliceSample & sample = samples.at(i);
                sample.SetFactorAndScale(flux_factor);
              }
            }
          }
        }

        std::pair<double, size_t> chi2_points
            = fThinSliceDriver->CalculateChi2(fSamples, fDataSet);
        return (chi2_points.first/chi2_points.second);
      },
      fTotalSignalParameters + fTotalFluxParameters);

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
  fReducedIncidentChi2 = pset.get<bool>("ReducedIncidentChi2");
  fPlotStyle = pset.get<std::vector<std::pair<int,int>>>("PlotStyle");
  fPlotRebinned = pset.get<bool>("PlotRebinned");
  fRandomStart = pset.get<bool>("RandomStart");

  fDriverName = pset.get<std::string>("DriverName");
  fAnalysisOptions = pset.get<fhicl::ParameterSet>("AnalysisOptions");
  fDoFakeData = pset.get<bool>("DoFakeData");
}

