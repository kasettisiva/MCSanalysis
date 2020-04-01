// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"

// ROOT includes
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TAxis.h>
#include <TDirectory.h>
#include <Math/MinimizerOptions.h>

#include <RooFitResult.h>
#include <RooWorkspace.h>

#include <RooStats/ModelConfig.h>
#include <RooStats/HistFactory/Sample.h>
#include <RooStats/HistFactory/Systematics.h>
#include <RooStats/HistFactory/HistoToWorkspaceFactory.h>
#include <RooStats/HistFactory/HistoToWorkspaceFactoryFast.h>
#include <RooStats/HistFactory/MakeModelAndMeasurementsFast.h>
#include <RooStats/HistFactory/PiecewiseInterpolation.h>

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ProtoDUNEFit.h"
#include "ProtoDUNEFitUtils.h"
#include "ProtoDUNESelectionUtils.h"
#include "MCToyGenerationAndFit.h"


//********************************************************************
protoana::ProtoDUNEFit::ProtoDUNEFit(){
  //********************************************************************

}

//********************************************************************
protoana::ProtoDUNEFit::ProtoDUNEFit(std::string configPath){
  //********************************************************************

  Configure(configPath);

}

//********************************************************************
protoana::ProtoDUNEFit::~ProtoDUNEFit(){
  //********************************************************************

}

//********************************************************************
void protoana::ProtoDUNEFit::BuildWorkspace(TString Outputfile, int analysis){
  //********************************************************************

  bool hfilled = false;

  // Pion analysis
  if(analysis == 1)
    hfilled = FillHistogramVectors_Pions();

  if(!hfilled) return;

  // Get pion flux
  //TH1* pionflux_ereco_histo  =  protoana::ProtoDUNESelectionUtils::FillMCFlux_Pions(_MCFileNames[0], _TruthTreeName, _TruthBinning, 1);
  //TH1* pionflux_etruth_histo =  protoana::ProtoDUNESelectionUtils::FillMCFlux_Pions(_MCFileNames[0], _TruthTreeName, _TruthBinning, 2);
  //TH1* pionflux_eff_histo = (TH1*)pionflux_ereco_histo->Clone();
  //pionflux_eff_histo->Divide(pionflux_etruth_histo);
  //pionflux_eff_histo->SetNameTitle("pionflux_eff_histo","Efficiency for incident particles");
  //TH1* pionrecoflux_ereco_histo  =  protoana::ProtoDUNESelectionUtils::FillMCFlux_Pions(_MCFileNames[0], _TruthTreeName, _RecoBinning, 3);
  //TH1* pionrecoflux_eff_histo = (TH1*)pionrecoflux_ereco_histo->Clone();
  //pionrecoflux_eff_histo->Divide(pionflux_etruth_histo);
  //pionrecoflux_eff_histo->SetNameTitle("pionrecoflux_eff_histo","Efficiency for incident particles");

  // Create measurement object
  RooStats::HistFactory::Measurement meas("ProtoDUNEFitExample","ProtoDUNE fit example");

  // Since this is not LHC, don't care about luminosity. Set to constant
  meas.SetLumi(1.0);
  meas.SetLumiRelErr(0.001);
  meas.AddConstantParam("Lumi");

  // SetExportOnly to false meaning that we will fit the model
  meas.SetExportOnly(false);

  // Add samples and channels to measurement
  AddSamplesAndChannelsToMeasurement(meas);
  if( _AddIncidentToMeasurement )
    AddIncidentSamplesAndChannelsToMeasurement(meas);
  // AddSidebandSamplesAndChannelsToMeasurement(meas);

  // Print info about meas
  meas.PrintTree();
  
  // Export meas to workspace
  RooStats::HistFactory::HistoToWorkspaceFactoryFast h2w(meas);
  RooWorkspace* ws = h2w.MakeCombinedModel(meas);

  // Define interpolation code: 6th order polynomial interpolation and linear extrapolation
  protoana::ProtoDUNEFitUtils::SetInterpolationCode(ws, 4);

  // Save initial snapshot
  protoana::ProtoDUNEFitUtils::SaveSnapshot(ws, Form("%s_initial_snapshot",ws->GetName()));

  // Fit and/or generate data using the workspace
  protoana::MCToyGenerationAndFit* fitandgen = new protoana::MCToyGenerationAndFit();

  // Fit options - must be defined first
  fitandgen->SetFitStrategy(_FitStrategy);
  fitandgen->SetMinimiser(_Minimizer);
  if(_EnableMinosError)
    fitandgen->EnableMinosError();

  RooFitResult *fitresult = NULL;
  TTree *toys_tree = NULL;
  TString treename("protodUNE_dataresults");

  // Plots before fit
  std::vector<TString> truebinsnameVec;
  
  for(unsigned int l = 0; l < _BackgroundTopologyName.size(); l++){
    TString str = Form("%s", _BackgroundTopologyName[l].c_str());
    truebinsnameVec.push_back(str);
  }

  if(!_FitInReco){
    for(unsigned int l = 1; l < _TruthBinning.size(); l++){
      TString str = Form("Signal %.1f-%.1f", _TruthBinning[l-1], _TruthBinning[l]);
      truebinsnameVec.push_back(str);
    }
  }

  std::vector< TString > incidentNameVec;
  for( size_t l = 0; l < _IncidentTopologyName.size(); ++l ){
    TString str = Form("%s", _IncidentTopologyName[l].c_str());
    incidentNameVec.push_back(str); 
  }

  std::vector<TCanvas*> bfplots = protoana::ProtoDUNEFitUtils::PlotDatasetsAndPdfs(ws, "beforefit", "Poisson", "ratio", truebinsnameVec, _RecoBinning, incidentNameVec, "Before Fit", _DoNegativeReco);

  RooAbsData *asimovdata = ws->data("asimovData");
  std::vector<TCanvas*> bfAsimovplots = protoana::ProtoDUNEFitUtils::PlotDatasetsAndPdfs(ws, "asimov", "Poisson", "ratio", truebinsnameVec, _RecoBinning, incidentNameVec, "Asimov Dataset", _DoNegativeReco, asimovdata);

  // ----------------------------------------------------------------------------------------------------
  // Check if this is MC toys case
  // ----------------------------------------------------------------------------------------------------

  if(_NToys > 0){
    toys_tree = fitandgen->GenerateAndFit(ws, _NToys);
    toys_tree->SetNameTitle("protodUNE_mctoysresults", "protodUNE_mctoysresults");

    // For each set of plots need to clone the output tree to avoid ROOT crashes
    TTree *mctoys_results1 = (TTree*)toys_tree->Clone();
    TCanvas* nuisancecanvas =  protoana::ProtoDUNEFitUtils::PlotNuisanceParameters(mctoys_results1, ws);

    TTree *mctoys_results2 = (TTree*)toys_tree->Clone();
    TCanvas* pullscanvas =  protoana::ProtoDUNEFitUtils::PlotParametersPull(mctoys_results2, ws);

    TTree *mctoys_results3 = (TTree*)toys_tree->Clone();
    TCanvas* nuispullscanvas =  protoana::ProtoDUNEFitUtils::PlotNuisanceParametersPull(mctoys_results3, ws);

    TTree *mctoys_results4 = (TTree*)toys_tree->Clone();
    TCanvas* avresultcanvas = protoana::ProtoDUNEFitUtils::PlotAverageResultsFromToys(mctoys_results4, ws, "Channel0", "SigTopo");

    TFile *f = new TFile(Outputfile.Data(), "UPDATE");

    // Save all histograms used in the measurement
    meas.writeToFile(f);

    toys_tree->Write();
    nuisancecanvas->Write();
    pullscanvas->Write();
    nuispullscanvas->Write();
    avresultcanvas->Write();

    for(unsigned int i=0; i < bfAsimovplots.size(); i++){
      bfAsimovplots[i]->Write();
    }
    for(unsigned int i=0; i < bfplots.size(); i++){
      bfplots[i]->Write();
    }

    // Histograms below are also saved from the measurement
    /*
    // Save all input histograms
    TDirectory *HistoDir = f->mkdir("OriginalHistograms");
    HistoDir->cd();
    
    for(unsigned int k = 0; k < _sighistos.size(); k++){
      _sighistos[k]->Write();
    }
    for(unsigned int k = 0; k < _bkghistos.size(); k++){
      _bkghistos[k]->Write();
    }
    for(unsigned int k = 0; k < _datahistos.size(); k++){
      _datahistos[k]->Write();
    }
    for(unsigned int k = 0; k < _truthsighistos.size(); k++){
      _truthsighistos[k]->Write();
    }
    
    // Decorate and save efficiency graphs
    for(unsigned int i=0; i < _efficiencyGraphs.size(); i++){
      TGraphAsymmErrors* effgraph = (TGraphAsymmErrors*)_efficiencyGraphs.at(i);

      TCanvas* ceffgraph = new TCanvas(Form("ceffgraph%i",i),Form("Efficiency for channel %i",i));

      TAxis *ax = effgraph->GetHistogram()->GetXaxis();
      Double_t x1 = ax->GetBinLowEdge(1); 
      //Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
      effgraph->GetHistogram()->GetXaxis()->Set(_TruthBinning.size(),x1,_TruthBinning.size());
      for(unsigned int i = 1; i < _TruthBinning.size(); i++){
	TString ibinstr = Form("%.1f-%.1f",_TruthBinning[i-1],_TruthBinning[i]);
	effgraph->GetHistogram()->GetXaxis()->SetBinLabel(i, ibinstr.Data());
      }

      effgraph->Draw("*a");
      effgraph->SetMarkerStyle(20);
      effgraph->SetMarkerColor(1);
      effgraph->SetTitle("Efficiency");
      effgraph->GetXaxis()->SetTitle("E_{true} at vertex [MeV]");
      effgraph->GetYaxis()->SetTitle("Efficiency");
      effgraph->GetYaxis()->SetTitleOffset(1.25);
      
      ceffgraph->Write();
    }
    
    HistoDir->cd("..");
    */

    f->Close();

    // Nothing else to do here
    return;
  }

  // ----------------------------------------------------------------------------------------------------
  // Continue with either asimov or data fit
  // ----------------------------------------------------------------------------------------------------

  if(_DoAsimovFit){
    fitresult = fitandgen->FitAsimovData(ws);
    treename = "protodUNE_asimovresults";
  }
  else{
    fitresult = fitandgen->FitData(ws);
  }

  // Try different fit configurations if the fit failed
  // Reset all values and errors
  //protoana::ProtoDUNEFitUtils::ResetAllValuesAndErrors(ws);

  if(fitresult->status() != 0){
    std::cout << "WARNING::Fit failed - trying with Scan" << std::endl;
    fitandgen->SetAlgorithm("Scan");
    fitresult = fitandgen->FitData(ws);
  }

  if(fitresult->status() != 0){
    if(_FitStrategy == 0){
      std::cout << "WARNING::Fit failed with fit strategy 0 - trying with fit strategy 1." << std::endl;
      fitandgen->SetAlgorithm(ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
      fitandgen->SetFitStrategy(1);
      fitresult = fitandgen->FitData(ws);
    }
  }

  if(fitresult->status() != 0){
    std::cout << "WARNING::Fit failed - trying with improve" << std::endl;
    fitandgen->SetAlgorithm("migradimproved");
    fitandgen->SetMinimiser("Minuit");
    fitresult = fitandgen->FitData(ws);
  }

  // Save results from data or asimov fit
  if(!fitresult){
    std::cerr << "No fit resuls found. Will exit!" << std::endl;
    return;
  }

  // Import fit result in workspace
  ws->import(*fitresult,kTRUE);

  // Print fit results on the screen
  //fitresult->Print();

  std::vector<TCanvas*> afplots = protoana::ProtoDUNEFitUtils::PlotDatasetsAndPdfs(ws, "afterfit", "Poisson", "ratio", truebinsnameVec, _RecoBinning, incidentNameVec, "After fit", _DoNegativeReco, NULL, fitresult);

  // Save post-fit workspace snapshot
  protoana::ProtoDUNEFitUtils::SaveSnapshot(ws, Form("%s_postfit_snapshot",ws->GetName()));
  protoana::ProtoDUNEFitUtils::SaveSnapshot(ws, Form("%s_postfitForPlots_snapshot",ws->GetName()));

  TH2* fitcov = protoana::ProtoDUNEFitUtils::GetFitCovariance(fitresult);
  TH2* fitcor = protoana::ProtoDUNEFitUtils::GetFitCorrelationMatrix(fitresult);

  toys_tree = fitandgen->RooFitResultToTTree(ws, fitresult);
  toys_tree->SetNameTitle(treename.Data(),treename.Data());

  // Plot NLL
  //std::vector<TCanvas*> nllplots = protoana::ProtoDUNEFitUtils::PlotNLL(ws,"PDFit",fitresult);

  TTree *mctoys_results1 = (TTree*)toys_tree->Clone();
  TCanvas* nuisancecanvas = protoana::ProtoDUNEFitUtils::PlotNuisanceParameters(mctoys_results1, ws);

  TTree *mctoys_results2 = (TTree*)toys_tree->Clone();
  TCanvas* avresultcanvas = protoana::ProtoDUNEFitUtils::PlotAverageResultsFromToys(mctoys_results2, ws, "POI", "POI");

  //TTree *mctoys_results3 = (TTree*)toys_tree->Clone();
  //TCanvas* avresultcanvas_inc = protoana::ProtoDUNEFitUtils::PlotAverageResultsFromToys(mctoys_results3, ws, "POI", "POI");
  
  // Fit fixing nuisance parameters
  //protoana::ProtoDUNEFitUtils::LoadSnapshot(ws, Form("%s_postfitForPlots_snapshot",ws->GetName()));
  //protoana::ProtoDUNEFitUtils::MakeNuisanceParamsConstant(ws,"");
  //RooFitResult *fitresultnosysts = fitandgen->FitData(ws);
  //if(fitresultnosysts){
  //ws->import(*fitresultnosysts,kTRUE);
  //}
  //protoana::ProtoDUNEFitUtils::SaveSnapshot(ws, Form("%s_postfitNosysts_snapshot",ws->GetName()));

  // Save workspace
  protoana::ProtoDUNEFitUtils::SaveWorkspace(ws, Outputfile);

  TFile *f = new TFile(Outputfile.Data(), "UPDATE");

  // Save all histograms used in the measurement
  meas.writeToFile(f);

  fitcov->Write();
  fitcor->Write();
  toys_tree->Write();
  if(nuisancecanvas)
    nuisancecanvas->Write();
  if(avresultcanvas)
    avresultcanvas->Write();
  //if(avresultcanvas_inc)
  //avresultcanvas_inc->Write();
  
  for(unsigned int i=0; i < bfAsimovplots.size(); i++){
    bfAsimovplots[i]->Write();
  }
  for(unsigned int i=0; i < bfplots.size(); i++){
    bfplots[i]->Write();
  }
  for(unsigned int i=0; i < afplots.size(); i++){
    afplots[i]->Write();
  }

  //pionflux_etruth_histo->Write();
  //pionflux_ereco_histo->Write();
  //pionflux_eff_histo->Write();
  //pionrecoflux_ereco_histo->Write();
  //pionrecoflux_eff_histo->Write();
  
  TDirectory *HistoDir = f->mkdir("OriginalHistograms");
  HistoDir->cd();

  for(unsigned int k = 0; k < _sighistos.size(); k++){
    _sighistos[k]->Write();
  }
  for(unsigned int k = 0; k < _bkghistos.size(); k++){
    _bkghistos[k]->Write();
  }
  for(unsigned int k = 0; k < _datahistos.size(); k++){
    _datahistos[k]->Write();
  }
  for(unsigned int k = 0; k < _truthsighistos.size(); k++){
    _truthsighistos[k]->Write();
  }

  // Decorate and save efficiency graphs
  for(unsigned int i=0; i < _efficiencyGraphs.size(); i++){
    TGraphAsymmErrors* effgraph = (TGraphAsymmErrors*)_efficiencyGraphs.at(i);
    
    TCanvas* ceffgraph = new TCanvas(Form("ceffgraph%i",i),Form("Efficiency for channel %i",i));
    
    DecorateEfficiency( effgraph );
    effgraph->Draw("*a");

    ceffgraph->Write();
  }

  for(unsigned int k = 0; k < _incsighistos.size(); k++){
    _incsighistos[k]->Write();
  }
  for(unsigned int k = 0; k < _incbkghistos.size(); k++){
    _incbkghistos[k]->Write();
  }
  for(unsigned int k = 0; k < _incdatahistos.size(); k++){
    _incdatahistos[k]->Write();
  }

  //Incident efficiency
  TCanvas* ceffgraph = new TCanvas( "ceffgraph", "Efficiency" );

  _incidentEfficiency->Draw("a");
  DecorateEfficiency(_incidentEfficiency);

  ceffgraph->Write();

  _incidentEfficiency->Write();
  _incidentEfficiencyNum->Write();
  _incidentEfficiencyDenom->Write();
  ////////////////////////////
  
  //Interacting efficiency
  for( size_t i = 0; i < _interactingEfficiencyDenoms.size(); ++i ){
    _interactingEfficiencyDenoms[i]->Write();
    _interactingEfficiencyNums[i]->Write();

    DecorateEfficiency(_interactingEfficiencies[i]);
    _interactingEfficiencies[i]->Write();

  }

  HistoDir->cd("..");

  f->Close();
  
  return;

}

//********************************************************************
void protoana::ProtoDUNEFit::AddSamplesAndChannelsToMeasurement(RooStats::HistFactory::Measurement& meas){
  //********************************************************************

  const int nmcchannels = _MCFileNames.size();
  for(int i=0; i < nmcchannels; i++){
    TString channelname = Form("Channel%s", _ChannelNames[i].c_str());
    RooStats::HistFactory::Channel channel(channelname.Data());
    // Add data to channel
    RooStats::HistFactory::Data data;
    for(unsigned int j=0; j < _datahistos.size(); j++){
      // Clone histogram
      TH1D* htemp = (TH1D*)(_datahistos.at(j)->Clone());

      TString hname(htemp->GetName());
      if(hname.Contains(channelname.Data()) && hname.Contains("Data")){
	mf::LogInfo("AddSamplesAndChannelsToMeasurement") << "Adding dataset " << hname.Data() << " to channel " << channelname.Data() << " with " << htemp->Integral() << " events";
	data.SetHisto(htemp);
	channel.SetData(data);
      }
    }
    
    // Get histograms with the systematics 
    std::vector<TH1*> systvec;
    if (_EnableSystematicError) {
      std::cout << "Attempting to get systs from " <<  _SystFileNames[i] << std::endl; 
      systvec = protoana::ProtoDUNEFitUtils::GetSystHistograms(_SystFileNames[i]); 
      std::cout << "Got " << systvec.size() << " systematic histograms" << std::endl;
    }
    
    // Add bkg samples to channel
    for(unsigned int j=0; j < _bkghistos.size(); j++){
      // Clone histogram
      TH1D* htemp = (TH1D*)(_bkghistos.at(j)->Clone());

      TString hname(htemp->GetName());
      //if(htemp->GetEntries() == 0 || htemp->Integral() == 0) continue;
      if(hname.Contains(channelname.Data()) && hname.Contains("MC")){
	mf::LogInfo("AddSamplesAndChannelsToMeasurement") << "Adding MC sample " << hname.Data() << " to channel " << channelname.Data() << " with " << htemp->Integral() << " events";
	TString samplename = hname + TString("_sample");
	RooStats::HistFactory::Sample sample(samplename.Data());
	sample.SetNormalizeByTheory(true);
	
	// Check to enable statistical uncertainty
	if(_EnableStatisticalError){
	  //mf::LogInfo("AddSamplesAndChannelsToMeasurement") << "Enable statistical uncertainty";
	  sample.ActivateStatError();
	  TH1* staterrorhisto = protoana::ProtoDUNEFitUtils::GetStatsSystHistogram(htemp);
	  RooStats::HistFactory::StatError& staterror = sample.GetStatError();
	  staterror.SetUseHisto();
	  staterror.SetErrorHist(staterrorhisto);
	}

	if(_EnableSystematicError){
	  ApplySystematicToSample(sample, htemp, systvec, false, false);
	}

	// Set histogram for sample
	sample.SetHisto(htemp);
	
	if(htemp->Integral() > 0){
	  //TString poiname = hname;
	  //poiname.ReplaceAll("MC","POI");
	  //poiname.ReplaceAll("_Histo","");
	  //poiname.ReplaceAll(".","");
	  //poiname.ReplaceAll("-","_");
	  //poiname.ReplaceAll(channelname.Data(),"");
	  //poiname.ReplaceAll("__","_");
	  //meas.SetPOI(poiname.Data());
	  //sample.AddNormFactor(poiname.Data(), 1.0, 0.0, 2.0);
	  //mf::LogInfo("BuildMeasurement") << "Sample " << sample.GetName() << " has normalisation parameter " << poiname.Data();
	}

	// Add sample to channel
	channel.AddSample(sample);
      }
    }

    // Add signal samples to channel
    for(unsigned int j=0; j < _sighistos.size(); j++){
      // Clone histogram
      TH1D* htemp = (TH1D*)(_sighistos.at(j)->Clone());

      TString hname(htemp->GetName());
      //if(htemp->GetEntries() == 0 || htemp->Integral() == 0) continue;
      if(hname.Contains(channelname.Data()) && hname.Contains("MC")){
	mf::LogInfo("AddSamplesAndChannelsToMeasurement") << "Adding MC sample" << hname.Data() << " to channel " << channelname.Data() << " with " << htemp->Integral() << " events";
	TString samplename = hname + TString("_sample");
	RooStats::HistFactory::Sample sample(samplename.Data());
	sample.SetNormalizeByTheory(true);
	
	// Check to enable statistical uncertainty
	if(_EnableStatisticalError){
	  //mf::LogInfo("AddSamplesAndChannelsToMeasurement") << "Enable statistical uncertainty";
	  sample.ActivateStatError();
	  TH1* staterrorhisto = protoana::ProtoDUNEFitUtils::GetStatsSystHistogram(htemp);
	  RooStats::HistFactory::StatError& staterror = sample.GetStatError();
	  staterror.SetUseHisto();
	  staterror.SetErrorHist(staterrorhisto);
	}

	if(_EnableSystematicError){
	  ApplySystematicToSample(sample, htemp, systvec, true, _NormalisedSystematic);
	}

	// Set histogram for sample
	sample.SetHisto(htemp);

	// Add POI to this sample
	TString poiname = hname;
	poiname.ReplaceAll("MC","POI");
	poiname.ReplaceAll("_Histo","");
	poiname.ReplaceAll(".0","");
	poiname.ReplaceAll("-","_");
	poiname.ReplaceAll(channelname.Data(),"");
	poiname.ReplaceAll("__","_");
	poiname.ReplaceAll("_0_1000000","");
	meas.SetPOI(poiname.Data()); // AddPOI would also work
	sample.AddNormFactor(poiname.Data(), 1.0, 0.0, 100.0);
	mf::LogInfo("AddSamplesAndChannelsToMeasurement") << "Sample " << sample.GetName() << " has normalisation parameter " << poiname.Data();

	// Add sample to channel
	channel.AddSample(sample);
      }
    }

    // Statistical uncertainty less than 1% is ignored
    channel.SetStatErrorConfig(_IgnoreStatisticalErrorBelow,"Poisson"); // Poisson or Gaussian

    // Add channel to measurement
    mf::LogInfo("AddSamplesAndChannelsToMeasurement") <<  "Adding channel " << channel.GetName() << " to measurement " << meas.GetName();
    meas.AddChannel(channel);
  }

}

//********************************************************************
void protoana::ProtoDUNEFit::AddIncidentSamplesAndChannelsToMeasurement(RooStats::HistFactory::Measurement& meas){
  //********************************************************************

  TString channelname("ChannelIncident");
  RooStats::HistFactory::Channel channel(channelname.Data());
  // Add data to channel
  RooStats::HistFactory::Data data;
  for(unsigned int j=0; j < _incdatahistos.size(); j++){
    // Clone histogram
    TH1D* htemp = (TH1D*)(_incdatahistos.at(j)->Clone());

    mf::LogInfo("AddIncidentSamplesAndChannelsToMeasurement") << "Adding Incident dataset " << htemp->GetName() << " to channel " << channelname.Data() << " with " << htemp->Integral() << " events";
    data.SetHisto(htemp);
    channel.SetData(data);
  }
    
  // Get histograms with the systematics 
  //std::vector<TH1*> systvec;
  //if(_EnableSystematicError)
  //systvec = protoana::ProtoDUNEFitUtils::GetSystHistograms(_SystFileNames[i]); 
    
  // Add bkg samples to channel
  for(unsigned int j=0; j < _incbkghistos.size(); j++){
    // Clone histogram
    TH1D* htemp = (TH1D*)(_incbkghistos.at(j)->Clone());
    
    TString hname(htemp->GetName());
    //if(htemp->GetEntries() == 0 || htemp->Integral() == 0) continue;
   
    mf::LogInfo("AddIncidentSamplesAndChannelsToMeasurement") << "Adding Incident MC sample " << hname.Data() << " to channel " << channelname.Data() << " with " << htemp->Integral() << " events";
    TString samplename = hname + TString("_sample");
    RooStats::HistFactory::Sample sample(samplename.Data());
    sample.SetNormalizeByTheory(true);
    
    // Check to enable statistical uncertainty
    if(_EnableStatisticalError){
      //mf::LogInfo("AddIncidentSamplesAndChannelsToMeasurement") << "Enable statistical uncertainty";
      sample.ActivateStatError();
      TH1* staterrorhisto = protoana::ProtoDUNEFitUtils::GetStatsSystHistogram(htemp);
      RooStats::HistFactory::StatError& staterror = sample.GetStatError();
      staterror.SetUseHisto();
      staterror.SetErrorHist(staterrorhisto);
    }
    
    //if(_EnableSystematicError){
    //ApplySystematicToSample(sample, htemp, systvec, false, false);
    //}

    // Set histogram for sample
    sample.SetHisto(htemp);
    
    if(htemp->Integral() > 0){
      //TString poiname = hname;
      //poiname.ReplaceAll("MC","POI");
      //poiname.ReplaceAll("_Histo","");
      //poiname.ReplaceAll(".","");
      //poiname.ReplaceAll("-","_");
      //poiname.ReplaceAll(channelname.Data(),"");
      //poiname.ReplaceAll("__","_");
      //meas.SetPOI(poiname.Data());
      //sample.AddNormFactor(poiname.Data(), 1.0, 0.0, 2.0);
      //mf::LogInfo("BuildMeasurement") << "Sample " << sample.GetName() << " has normalisation parameter " << poiname.Data();
    }
    
    // Add sample to channel
    channel.AddSample(sample);
  }

  // Add signal samples to channel
  for(unsigned int j=0; j < _incsighistos.size(); j++){
    // Clone histogram
    TH1D* htemp = (TH1D*)(_incsighistos.at(j)->Clone());
    
    TString hname(htemp->GetName());
    //if(htemp->GetEntries() == 0 || htemp->Integral() == 0) continue;

    mf::LogInfo("AddSamplesAndChannelsToMeasurement") << "Adding Incident MC sample" << hname.Data() << " to channel " << channelname.Data() << " with " << htemp->Integral() << " events";
    TString samplename = hname + TString("_sample");
    RooStats::HistFactory::Sample sample(samplename.Data());
    sample.SetNormalizeByTheory(true);
    
    // Check to enable statistical uncertainty
    if(_EnableStatisticalError){
      //mf::LogInfo("AddIncidentSamplesAndChannelsToMeasurement") << "Enable statistical uncertainty";
      sample.ActivateStatError();
      TH1* staterrorhisto = protoana::ProtoDUNEFitUtils::GetStatsSystHistogram(htemp);
      RooStats::HistFactory::StatError& staterror = sample.GetStatError();
      staterror.SetUseHisto();
      staterror.SetErrorHist(staterrorhisto);
    }
    
    //if(_EnableSystematicError){
    //ApplySystematicToSample(sample, htemp, systvec, true, _NormalisedSystematic);
    //}
    
    // Set histogram for sample
    sample.SetHisto(htemp);
    
    // Add POI to this sample
    TString poiname = hname;
    poiname.ReplaceAll("MC","POI");
    poiname.ReplaceAll("_Histo","");
    poiname.ReplaceAll(".0","");
    poiname.ReplaceAll("-","_");
    meas.SetPOI(poiname.Data()); // AddPOI would also work
    sample.AddNormFactor(poiname.Data(), 1.0, 0.0, 100.0);
    mf::LogInfo("AddIncidentSamplesAndChannelsToMeasurement") << "Sample " << sample.GetName() << " has normalisation parameter " << poiname.Data();
    
    // Add sample to channel
    channel.AddSample(sample);
  }

  // Statistical uncertainty less than 1% is ignored
  channel.SetStatErrorConfig(_IgnoreStatisticalErrorBelow,"Poisson"); // Poisson or Gaussian
  
  // Add channel to measurement
  mf::LogInfo("AddIncidentSamplesAndChannelsToMeasurement") <<  "Adding channel " << channel.GetName() << " to measurement " << meas.GetName();
  meas.AddChannel(channel);

}

//********************************************************************
void protoana::ProtoDUNEFit::AddSidebandSamplesAndChannelsToMeasurement(RooStats::HistFactory::Measurement& meas){
  //********************************************************************

  const int nsidmcchannels = _MCControlSampleFiles.size();

  for(int i=0; i < nsidmcchannels; i++){
    TString toponame = Form("Topo%i",_BackgroundTopology[i]);
    TString channelname = Form("SidebandChannel%i", i);
    RooStats::HistFactory::Channel channel(channelname.Data());
    // Add data to channel
    RooStats::HistFactory::Data data;
    for(unsigned int j=0; j < _siddatahistos.size(); j++){
      TH1D* htemp = (TH1D*)(_siddatahistos.at(j)->Clone());
      TString hname(htemp->GetName());
      if(hname.Contains(channelname.Data()) && hname.Contains("Data")){
	mf::LogInfo("AddSidebandSamplesAndChannelsToMeasurement") << "Adding sideband dataset " << hname.Data() << " to channel " << channelname.Data() << " with " << htemp->Integral() << " events";
	data.SetHisto(htemp);
	channel.SetData(data);
      }
    }

    // Add bkg samples to channel
    for(unsigned int j=0; j < _sidhistos.size(); j++){
      TH1D* htemp = (TH1D*)(_sidhistos.at(j)->Clone());

      TString hname(htemp->GetName());
      if(htemp->GetEntries() == 0 || htemp->Integral() == 0) continue;
      if(hname.Contains(channelname.Data()) && hname.Contains("MC")){
	mf::LogInfo("AddSidebandSamplesAndChannelsToMeasurement") << "Adding MC sideband sample " << hname.Data() << " to channel " << channelname.Data() << " with " <<  htemp->Integral() << " events";
	TString samplename = hname + TString("_sample");
	RooStats::HistFactory::Sample sample(samplename.Data());
	sample.SetNormalizeByTheory(true);
	
	// Check to enable statistical uncertainty
	if(_EnableStatisticalError){
	  //mf::LogInfo("AddSidebandSamplesAndChannelsToMeasurement") << "Enable statistical uncertainty";
	  sample.ActivateStatError();
	  TH1* staterrorhisto = protoana::ProtoDUNEFitUtils::GetStatsSystHistogram(htemp);
	  RooStats::HistFactory::StatError& staterror = sample.GetStatError();
	  staterror.SetUseHisto();
	  staterror.SetErrorHist(staterrorhisto);
	}

	//if(_EnableSystematicError){
	//ApplySystematicToSample(sample, htemp, systvec, false, false);
	//}

	// Set histogram for sample
	sample.SetHisto(htemp);

	// Add bkg normalisation parameter to this sample
	if(hname.Contains(toponame.Data())){
	  TString poiname = hname;
	  poiname.ReplaceAll("MC","POI");
	  poiname.ReplaceAll("_Histo","");
	  poiname.ReplaceAll(".","");
	  poiname.ReplaceAll("-","_");
	  poiname.ReplaceAll(channelname.Data(),"");
	  poiname.ReplaceAll("__","_");
	  meas.SetPOI(poiname.Data());
	  sample.AddNormFactor(poiname.Data(), 1.0, 0.0, 100.0);
	  mf::LogInfo("AddSidebandSamplesAndChannelsToMeasurement") << "Sample " << sample.GetName() << " has normalisation parameter " << poiname.Data();
	}

	// Add sample to channel
	channel.AddSample(sample);
      }
    }

    // Statistical uncertainty less than 1% is ignored
    channel.SetStatErrorConfig(_IgnoreStatisticalErrorBelow,"Poisson"); // Poisson or Gaussian

    // Add channel to measurement
    mf::LogInfo("AddSidebandSamplesAndChannelsToMeasurement") <<  "Adding sideband channel " << channel.GetName() << " to measurement " << meas.GetName();
    meas.AddChannel(channel);
  }

}

//********************************************************************
bool protoana::ProtoDUNEFit::FillHistogramVectors_Pions(){
  //********************************************************************

  const int nmcchannels   = _MCFileNames.size();
  const int ndatachannels = _DataFileNames.size();
  const int nchannelnames = _ChannelNames.size();
  const int ninctoponames = _IncidentTopologyName.size();

  const int nbkgtopo      = _BackgroundTopology.size();
  const int nbkgtoponames = _BackgroundTopologyName.size();
  const int nsigtopo      = _SignalTopology.size();
  const int nsigtoponames = _SignalTopologyName.size();
  const int ninctopo      = _IncidentTopology.size();

  if(nmcchannels != ndatachannels){
    mf::LogError("FillHistogramVectors_Pions") << "The number of data and MC channels is not the same. Check MCFileNames and DataFileNames in the fcl file. Time to die!";
    return false;
  }

  if(nmcchannels != nchannelnames || ndatachannels != nchannelnames){
    mf::LogError("FillHistogramVectors_Pions") << "The channel names do not correspond to the input files. Check MCFileNames, DataFileNames and ChannelNames in the fcl file. Time to die!";
    return false;
  }

  if(nbkgtopo != nbkgtoponames){
    mf::LogError("FillHistogramVectors_Pions") << "Background topologies and background names do not have the same length. Check BackgroundTopology and BackgroundTopologyName in the fcl file. Time to die!";
    return false;
  }

  if(nsigtopo != nsigtoponames){
    mf::LogError("FillHistogramVectors_Pions") << "Signal topologies and background name vectors do not have the same size. Check SignalTopology and SignalTopologyName in the fcl file. Time to die!";
    return false;
  }

  if(ninctopo != ninctoponames){
    mf::LogError("FillHistogramVectors_Pions") << "Incident topologies and name vectors do not have the same size. Check  IncidentTopology and IncidentTopologyName in the fcl file. Time to die!";
    return false;
  }

  // Get total number of data and MC triggers
  //int nmctriggers = protoana::ProtoDUNESelectionUtils::GetNTriggers_Pions(_MCFileNames[0], _RecoTreeName);
  //int ndatatriggers = protoana::ProtoDUNESelectionUtils::GetNTriggers_Pions(_DataFileNames[0], _RecoTreeName, false);
  //double mcnorm = (double)ndatatriggers/nmctriggers;
  //mf::LogInfo("FillHistogramVectors_Pions") << "Total number of MC triggers = " << nmctriggers << ", total number of data triggers = " << ndatatriggers << " , data/MC = " << mcnorm;

  Double_t tmin = _TruthBinning[0];
  Double_t tmax = _TruthBinning[_TruthBinning.size()-1];
  if(_FitInReco){
    tmin = 0.0;
    tmax = 1000000.0;
  }

  for(int i=0; i < nmcchannels; i++){
    for(int j=0; j < nbkgtopo; j++){
      int topo = _BackgroundTopology[j];
      _bkghistos.push_back(
          protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(
              _MCFileNames[i], _RecoTreeName, _RecoBinning, _ChannelNames[i],
              _BackgroundTopologyName[j], topo, _EndZCut, tmin, tmax,
              _DoNegativeReco));
      //_incbkghistos.push_back( protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(_MCFileNames[0], _RecoTreeName, _RecoBinning, i, topo, true) );
    }

    //TH1D* sigevenshisto = new TH1D(Form("sigevenshisto_channel%i",i),Form("sigevenshisto_channel%i",i),_TruthBinning.size()-1,0,_TruthBinning.size()-1);
    //TH1D* truevenshisto = new TH1D(Form("truevenshisto_channel%i",i),Form("truevenshisto_channel%i",i),_TruthBinning.size()-1,0,_TruthBinning.size()-1);

    for(int j=0; j < nsigtopo; j++){
      int topo = _SignalTopology[j];
      for(unsigned int k=1; k < _TruthBinning.size(); k++){
	if(_FitInReco){
          TH1* hsignal =
              protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(
                  _MCFileNames[i], _RecoTreeName, _RecoBinning,
                  _ChannelNames[i], _SignalTopologyName[j], topo, tmin, tmax,
                  _EndZCut);

	  // Make one histogram per bin
	  for(int k=1; k <= hsignal->GetNbinsX(); k++){
	    TString hname = Form("%s_RecoBin%i",hsignal->GetName(), k);
	    TH1 *hsignal_clone = (TH1*)hsignal->Clone();
	    hsignal_clone->SetName(hname.Data());
	    hsignal_clone->SetBinContent(k, hsignal->GetBinContent(k));
	    for(int kk=1; kk <= hsignal_clone->GetNbinsX(); kk++){
	      if(kk == k) continue;
	      hsignal_clone->SetBinContent(kk, 0.0);
	    }
	    _sighistos.push_back(hsignal_clone);
	  }
	}
	else{
          _sighistos.push_back(
              protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(
                  _MCFileNames[i], _RecoTreeName, _RecoBinning,
                  _ChannelNames[i], _SignalTopologyName[j], topo,
                  _TruthBinning[k-1], _TruthBinning[k], _EndZCut,
                  _DoNegativeReco));
	//sigevenshisto->SetBinContent(k, (_sighistos.back())->Integral() + sigevenshisto->GetBinContent(k));
	}
      }
       //_incsighistos.push_back( protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(_MCFileNames[0], _RecoTreeName, _RecoBinning, i, topo, 0, 0, true) );
      //TH1* inchistoall = protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(_MCFileNames[0], _RecoTreeName, _RecoBinning, i, topo, 0, 0, true);
      //for(Int_t k=1; k <= inchistoall->GetNbinsX(); k++){
      //TH1 *newhist = (TH1*)inchistoall->Clone();
      //newhist->SetNameTitle(Form("%s_RecoBin%i",inchistoall->GetName(),k), Form("%s in Reco Bin %i",inchistoall->GetName(),k));
      //for(Int_t kk=1; kk <= newhist->GetNbinsX(); kk++){
      //  if(kk != k){
      //    newhist->SetBinContent(kk,0);
      //    newhist->SetBinError(kk,0);
      //  }
      //}
      //_incsighistos.push_back(newhist);
      //}
    }

    //_truthsighistos.push_back( protoana::ProtoDUNESelectionUtils::FillMCTruthSignalHistogram_Pions( _MCFileNames[0], _TruthTreeName, _TruthBinning, i) );
    //truevenshisto->Add(_truthsighistos.back());

    //TGraphAsymmErrors* effgraph = new TGraphAsymmErrors(sigevenshisto,truevenshisto);
    //effgraph->SetNameTitle(Form("Efficiency_channel%i",i), Form("Efficiency for channel %i",i));
    //_efficiencyGraphs.push_back(effgraph); 
			       
  }

  for(int i=0; i < ndatachannels; i++){
    _datahistos.push_back( protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(_DataFileNames[i], _RecoTreeName, _RecoBinning, _ChannelNames[i]) );
    //_incdatahistos.push_back( protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(_DataFileNames[0], _RecoTreeName, _RecoBinning, i, true) );
  }

  for(int i=0; i < ninctopo; i++){
    TH1* inchisto = protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(_MCFileNames[0], _RecoTreeName, _RecoBinning, _ChannelNames[0], _IncidentTopologyName[i], _IncidentTopology[i], _EndZCut);
    for(int j=1; j < nmcchannels; j++){
      inchisto->Add(protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(_MCFileNames[j], _RecoTreeName, _RecoBinning, _ChannelNames[j], _IncidentTopologyName[i], _IncidentTopology[i], _EndZCut));
    }

    inchisto->SetNameTitle(Form("MC_ChannelIncident_%s_Histo",_IncidentTopologyName[i].c_str()), Form("Incident MC for topology %s", _IncidentTopologyName[i].c_str()));
    /*
    if(i == 0){
      // Split into multiple histograms
      for(int j=1; j <= inchisto->GetNbinsX(); j++){
	TString hname = Form("%s_IncBin%i",inchisto->GetName(),j);
	TH1* inchisto_h = (TH1*)inchisto->Clone();
	inchisto_h->Reset();
	inchisto_h->SetName(hname.Data());
	for(int k=1; k <= inchisto->GetNbinsX(); k++){
	  if(k == j) inchisto_h->SetBinContent(k, inchisto->GetBinContent(k));
	  else inchisto_h->SetBinContent(k, 0.0);
	}
	_incsighistos.push_back(inchisto_h);
      }   
    }
    else{*/
      _incbkghistos.push_back(inchisto);
    //}
  }

  TH1* incdatahisto = protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(_DataFileNames[0], _RecoTreeName, _RecoBinning, _ChannelNames[0], true);
  for(int i=1; i < ndatachannels; i++){
    incdatahisto->Add(protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(_DataFileNames[i], _RecoTreeName, _RecoBinning, _ChannelNames[i], true));
  }
  _incdatahistos.push_back(incdatahisto);

  std::pair<TH1 *, TH1 *> inc_eff_num_denom = 
        protoana::ProtoDUNESelectionUtils::GetMCIncidentEfficiency(
        _IncidentMCFileNames[0], _TruthTreeName, _TruthBinning, _EndZCut);

  _incidentEfficiencyNum = inc_eff_num_denom.first;
  _incidentEfficiencyDenom = inc_eff_num_denom.second;

  _incidentEfficiency = new TGraphAsymmErrors(_incidentEfficiencyNum,
        _incidentEfficiencyDenom);

  _incidentEfficiency->SetNameTitle("MC_Incident_Efficiency", "Efficiency");

  //_incidentEfficiency = protoana::ProtoDUNESelectionUtils::GetMCIncidentEfficiency( _IncidentMCFileNames[0], _TruthTreeName, _TruthBinning );

  //Interacting Efficiencies
  for( int i = 0; i < nsigtopo; ++i ){
    int topo = _SignalTopology[i];
    //_interactingEfficiencyDenoms.push_back( protoana::ProtoDUNESelectionUtils::GetMCInteractingEfficiencyDenominator( _IncidentMCFileNames[0], _TruthTreeName, _TruthBinning, _ChannelNames[i], _SignalTopologyName[i], topo ) );
    
    std::pair< TH1 *, TH1 * > eff_num_denom = 
        protoana::ProtoDUNESelectionUtils::GetMCInteractingEfficiency( 
            _IncidentMCFileNames[0], _TruthTreeName, _TruthBinning,
            _ChannelNames[i], _SignalTopologyName[i], topo, _EndZCut);

    _interactingEfficiencyNums.push_back(eff_num_denom.first); 
    _interactingEfficiencyDenoms.push_back(eff_num_denom.second); 

    TGraphAsymmErrors * eff = new TGraphAsymmErrors(eff_num_denom.first,
        eff_num_denom.second);

    std::string name = "MC_Channel_" + _ChannelNames[i] + "_" + 
        _SignalTopologyName[i] + "_Interacting_Efficiency";

    eff->SetNameTitle(name.c_str(), "Efficiency");

    _interactingEfficiencies.push_back(eff);
    
  }



  return true;

}

//********************************************************************
bool protoana::ProtoDUNEFit::ApplySystematicToSample(RooStats::HistFactory::Sample& sample, TH1* histo, std::vector<TH1*> systvec, bool hasnormfactor, bool isnorm){
  //********************************************************************

  if(!histo){
    std::cerr << "ERROR::No input histogram found! Will not apply systematics!" << std::endl;
    return false;
  }

  if(systvec.empty()){
    std::cout << "INFO::No detector systematic vector found! Will not apply systematics!" << std::endl;
    return false;
  }

  TString hname(histo->GetName());
  hname.ReplaceAll("_Histo","");

  if(_SystToConsider.size() != _SystType.size()){
    std::cerr << "Vectors SystToConsider and SystType do not have the same size. Will not apply systematics. Check your fcl file." << std::endl;
    return false;
  }
  
  for (size_t i = 0; i < _SystToConsider.size(); ++i) {
    TString systname(_SystToConsider[i].c_str());
    TString systtype(_SystType[i].c_str());
    
    TH1* highhist = NULL; TH1* lowhist = NULL;
    for (size_t j = 0; j < systvec.size(); ++j) {
      TH1* systhisto = (TH1*)systvec[j];
      TString systhistoname(systhisto->GetName());
      
      if (!systhistoname.Contains(hname.Data())) {
        continue;
      }
      if (!systhistoname.Contains(systname.Data())) {
        continue;
      }
      
      if(systhistoname.Contains("LOW") || systhistoname.Contains("Low") || systhistoname.Contains("low")){ 
	lowhist = systhisto;
      }
      if(systhistoname.Contains("HIGH") || systhistoname.Contains("High") || systhistoname.Contains("high")){ 
	highhist = systhisto;
      }
    }
    
    if(!highhist || !lowhist){
      std::cerr << "ERROR::Stage1: Systematic histograms not found! " <<
                   "Will not apply systematic " << systname.Data() <<
                   " to histogram " << histo->GetName() << ". Will skip!" <<
                   std::endl;
      continue;
    } 

    TH1* highsyst = protoana::ProtoDUNEFitUtils::GetSystematicHistoFromNominal(
        histo, highhist, "UP");
    TH1* lowsyst  = protoana::ProtoDUNEFitUtils::GetSystematicHistoFromNominal(
        histo, lowhist, "DOWN");

    if(!highsyst || !lowsyst){
      std::cerr << "ERROR::Stage2: Systematic histograms not found! " <<
      "Will not apply systematic " << systname.Data() << 
      " to histogram " << histo->GetName() << ". Will skip!" << std::endl;
      continue;
    } 

    Double_t highval = highsyst->Integral()/histo->Integral();
    Double_t lowval  = lowsyst->Integral()/histo->Integral();

    if(1.-lowsyst->Integral()/histo->Integral() < _IgnoreSystematicErrorBelow || highsyst->Integral()/histo->Integral()-1. < _IgnoreSystematicErrorBelow ){
      std::cout << "WARNING::Ignoring systematic uncertainty " << systname << " for sample " << hname << " as it is below the threshold " << _IgnoreSystematicErrorBelow <<  " , syst fraction = " << 1.-lowsyst->Integral()/histo->Integral() << " , " << highsyst->Integral()/histo->Integral()-1. << std::endl;
      continue;
    }

    // Case for normalisation systematic only
    if(systtype == "NormOnly"){
      if(isnorm && hname.Contains("Signal")){
	if(highval > 0. && lowval > 0.){
	  std::cout << "INFO::Normalised systematic " << systname.Data() << " is considered." << std::endl;
	  highsyst->Scale(1./highval);
	  lowsyst->Scale(1./lowval);
	}
      }
      sample.AddOverallSys(systname.Data(), lowval, highval);
      std::cout << "INFO::Adding norm systematic by default with name " << systname << " , to sample " << hname << " , with high value = " << highval << " , and low value = " << lowval << " , with nominal integral = " << histo->Integral() << " , high histogram integral = " << highsyst->Integral() << " , low histogram integral = " << lowsyst->Integral() << " , and syst fraction = " << 1.-lowsyst->Integral()/histo->Integral() << " , " << highsyst->Integral()/histo->Integral()-1. << std::endl;
      continue;
    }

    // Mixed (norm and/or shape) systematic case
    if(!hasnormfactor){ // Sample without normalisation parameter
      if(protoana::ProtoDUNEFitUtils::IsSingleBinHisto(histo)){ // Single bin samples and without norm factor are considered as norm syst
	sample.AddOverallSys(systname.Data(), lowval, highval);
	std::cout << "INFO::Adding single bin and without norm parameter norm systematic " << systname << " , to sample " << hname << " , with high value = " << highval << " , and low value = " << lowval << " , with nominal integral = " << histo->Integral() << " , high histogram integral = " << highsyst->Integral() << " , low histogram integral = " << lowsyst->Integral() << " , and syst fraction = " << 1.-lowsyst->Integral()/histo->Integral() << " , " << highsyst->Integral()/histo->Integral()-1. << std::endl;
      }
      else{ // multi-bin region and without norm parameter
	RooStats::HistFactory::HistoSys syst;
	syst.SetName(systname.Data());
	syst.SetHistoHigh(highsyst);
	syst.SetHistoLow(lowsyst);
	sample.AddHistoSys(syst);
	std::cout << "INFO::Adding multi-bin and without norm parameter shape+norm systematic " << systname << " , to sample " << hname << " , with nominal integral = " << histo->Integral() << " , high histogram integral = " << highsyst->Integral() << " , low histogram integral = " << lowsyst->Integral() <<  ", and syst fraction = " << 1.-lowsyst->Integral()/histo->Integral() << " , " << highsyst->Integral()/histo->Integral()-1. << std::endl;
      }
    }
    else{ // Sample with normalisation parameter
      // Apply normalised systematic uncertainty
      if(isnorm && hname.Contains("Signal")){
	if(highval > 0. && lowval > 0.){
	  std::cout << "INFO::Normalised systematic " << systname.Data() << " is considered." << std::endl;
	  highsyst->Scale(1./highval);
	  lowsyst->Scale(1./lowval);
	}
      }
      
      if(protoana::ProtoDUNEFitUtils::IsSingleBinHisto(histo)){ // Single bin samples and with norm factor are considered as norm syst
	sample.AddOverallSys(systname.Data(), lowval, highval);
	std::cout << "INFO::Adding norm systematic " << systname << " , to sample " << hname << " , with high value = " << highval << " , and low value = " << lowval << " , with nominal integral = " << histo->Integral() << " , high histogram integral = " << highsyst->Integral() << " , low histogram integral = " << lowsyst->Integral() << " , and syst fraction = " << 1.-lowsyst->Integral()/histo->Integral() << " , " << highsyst->Integral()/histo->Integral()-1. << std::endl;
      }
      else{ // multi-bin region and with norm parameter
	RooStats::HistFactory::HistoSys syst;
	syst.SetName(systname.Data());
	syst.SetHistoHigh(highsyst);
	syst.SetHistoLow(lowsyst);
	sample.AddHistoSys(syst);
	std::cout << "INFO::Adding shape+norm systematic " << systname << " , to sample " << hname << " , with nominal integral = " << histo->Integral() << " , high histogram integral = " << highsyst->Integral() << " , low histogram integral = " << lowsyst->Integral() <<  ", and syst fraction = " << 1.-lowsyst->Integral()/histo->Integral() << " , " << highsyst->Integral()/histo->Integral()-1. << std::endl;
      }
    }
  }

  return true;

}

//********************************************************************
bool protoana::ProtoDUNEFit::Configure(std::string configPath){
  //********************************************************************

  cet::filepath_first_absolute_or_lookup_with_dot cetfilelook("FHICL_FILE_PATH");

  // Parse configuration file
  fhicl::intermediate_table tbl;
  fhicl::parse_document(configPath, cetfilelook, tbl);

  // Convert to ParameterSet
  fhicl::ParameterSet pset;
  fhicl::make_ParameterSet(tbl, pset);

  _RecoTreeName                = pset.get<std::string>("RecoTreeName");
  _Minimizer                   = pset.get<std::string>("Minimizer");
  _TruthTreeName               = pset.get<std::string>("TruthTreeName");

  _DataFileNames               = pset.get< std::vector<std::string> >("DataFileNames");
  _MCFileNames                 = pset.get< std::vector<std::string> >("MCFileNames");
  _MCControlSampleFiles        = pset.get< std::vector<std::string> >("MCControlSampleFiles");
  _DataControlSampleFiles      = pset.get< std::vector<std::string> >("DataControlSampleFiles");
  _IncidentMCFileNames         = pset.get< std::vector<std::string> >("IncidentMCFileNames");
  _SystFileNames               = pset.get< std::vector<std::string> >("SystFileNames");
  _SystToConsider              = pset.get< std::vector<std::string> >("SystToConsider");
  _SystType                    = pset.get< std::vector<std::string> >("SystType");
  _BackgroundTopologyName      = pset.get< std::vector<std::string> >("BackgroundTopologyName");
  _SignalTopologyName          = pset.get< std::vector<std::string> >("SignalTopologyName");
  _IncidentTopologyName        = pset.get< std::vector<std::string> >("IncidentTopologyName");
  _ChannelNames                = pset.get< std::vector<std::string> >("ChannelNames");
  
  _FitStrategy                 = pset.get<int>("FitStrategy");
  _NToys                       = pset.get<int>("NToys");

  _IgnoreStatisticalErrorBelow = pset.get<double>("IgnoreStatisticalErrorBelow");
  _IgnoreSystematicErrorBelow  = pset.get<double>("IgnoreSystematicErrorBelow");
  
  _EnableMinosError            = pset.get<bool>("EnableMinosError");
  _DoAsimovFit                 = pset.get<bool>("DoAsimovFit");
  _EnableStatisticalError      = pset.get<bool>("EnableStatisticalError");
  _EnableSystematicError       = pset.get<bool>("EnableSystematicError");
  _NormalisedSystematic        = pset.get<bool>("NormalisedSystematic");

  _RecoBinning                 = pset.get< std::vector<double> >("RecoBinning");
  _TruthBinning                = pset.get< std::vector<double> >("TruthBinning");
  _FitInReco                   = pset.get<bool>("FitInReco");

  _SignalTopology              = pset.get< std::vector<int> >("SignalTopology");
  _BackgroundTopology          = pset.get< std::vector<int> >("BackgroundTopology");
  _IncidentTopology            = pset.get< std::vector<int> >("IncidentTopology");

  _AddIncidentToMeasurement    = pset.get<bool>("AddIncidentToMeasurement");
  _DoNegativeReco              = pset.get<bool>("DoNegativeReco"); 
  _EndZCut                     = pset.get<double>("EndZCut");

  return true;

}

void protoana::ProtoDUNEFit::DecorateEfficiency( TGraphAsymmErrors * eff ){
  TAxis * ax = eff->GetHistogram()->GetXaxis();
  double x1 = ax->GetBinLowEdge(1); 
  /*eff->GetHistogram()->GetXaxis()*/ax->Set((_TruthBinning.size() - 1), 
                                             x1, (_TruthBinning.size() - 1));

  for(size_t i = 1; i < _TruthBinning.size(); ++i) {
    TString label = Form("%.1f-%.1f", _TruthBinning[i-1], _TruthBinning[i]);
    /*eff->GetHistogram()->GetXaxis()*/ax->SetBinLabel(i, label.Data());
  }
  
  //eff->Draw("*a");
  eff->SetMarkerStyle(20);
  eff->SetMarkerColor(1);
  eff->SetTitle("Efficiency");
  eff->GetXaxis()->SetTitle("E_{true} at vertex [MeV]");
  eff->GetYaxis()->SetTitle("Efficiency");
  eff->GetYaxis()->SetTitleOffset(1.25);

}
