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
#include <TDirectory.h>
#include <RooFitResult.h>
#include "RooWorkspace.h"
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

  // Create measurement object
  RooStats::HistFactory::Measurement meas("ProtoDUNEFitExample","ProtoDUNE fit example");

  // Since this is not LHC, don't care about luminosity. Set to constant
  meas.SetLumi(1.0);
  meas.SetLumiRelErr(0.001);
  meas.AddConstantParam("Lumi");

  // SetExportOnly to false meaning that we will fit the model
  meas.SetExportOnly(false);

  // Add samples and channels to measurement
  this->BuildMeasurement(meas);
  
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

  for(unsigned int l = 1; l <= _TruthBinning.size(); l++){
    TString str = Form("Signal %.1f-%.1f", _TruthBinning[l-1], _TruthBinning[l]);
    truebinsnameVec.push_back(str);
  }

  Double_t totalchi2bf = 0.0;
  std::vector<TCanvas*> bfplots = protoana::ProtoDUNEFitUtils::PlotDatasetsAndPdfs(ws, "beforefit", "Poisson", "pull", totalchi2bf, truebinsnameVec, _RecoBinning, "Before Fit");

  RooAbsData *asimovdata = ws->data("asimovData");
  totalchi2bf = 0.0;
  std::vector<TCanvas*> bfAsimovplots = protoana::ProtoDUNEFitUtils::PlotDatasetsAndPdfs(ws, "asimovbeforefit", "Poisson", "pull", totalchi2bf, truebinsnameVec, _RecoBinning, "Asimov Dataset", asimovdata);

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
    
    // Decorate and save efficiency graphs
    for(unsigned int i=0; i < _efficiencyGraphs.size(); i++){
      TGraphAsymmErrors* effgraph = (TGraphAsymmErrors*)_efficiencyGraphs.at(i);
      TCanvas* ceffgraph = new TCanvas(Form("ceffgraph%i",i),Form("Efficiency for channel %i",i));
      effgraph->Draw("*a");
      effgraph->SetMarkerStyle(20);
      effgraph->SetMarkerColor(1);
      effgraph->SetTitle("Efficiency");
      effgraph->GetXaxis()->SetTitle("True bin");
      effgraph->GetYaxis()->SetTitle("Efficiency");
      effgraph->GetYaxis()->SetTitleOffset(1.25);
      
      ceffgraph->Write();
    }
    
    HistoDir->cd("..");
    
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

  // Save results from data or asimov fit
  if(!fitresult){
    std::cerr << "No fit resuls found. Will exit!" << std::endl;
    return;
  }

  // Import fit result in workspace
  ws->import(*fitresult,kTRUE);

  // Print fit results on the screen
  //fitresult->Print();

  totalchi2bf = 0.0;
  std::vector<TCanvas*> afplots = protoana::ProtoDUNEFitUtils::PlotDatasetsAndPdfs(ws, "afterfit", "Poisson", "pull", totalchi2bf, truebinsnameVec, _RecoBinning, "After fit", NULL, fitresult);

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
  TCanvas* avresultcanvas = protoana::ProtoDUNEFitUtils::PlotAverageResultsFromToys(mctoys_results2, ws, "Channel0", "SigTopo");
  
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
  fitcov->Write();
  fitcor->Write();
  toys_tree->Write();
  nuisancecanvas->Write();
  avresultcanvas->Write();
  
  for(unsigned int i=0; i < bfAsimovplots.size(); i++){
    bfAsimovplots[i]->Write();
  }
  for(unsigned int i=0; i < bfplots.size(); i++){
    bfplots[i]->Write();
  }
  for(unsigned int i=0; i < afplots.size(); i++){
    afplots[i]->Write();
  }
  
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

  // Decorate and save efficiency graphs
  for(unsigned int i=0; i < _efficiencyGraphs.size(); i++){
    TGraphAsymmErrors* effgraph = (TGraphAsymmErrors*)_efficiencyGraphs.at(i);
    TCanvas* ceffgraph = new TCanvas(Form("ceffgraph%i",i),Form("Efficiency for channel %i",i));
     effgraph->Draw("*a");
     effgraph->SetMarkerStyle(20);
     effgraph->SetMarkerColor(1);
     effgraph->SetTitle("Efficiency");
     effgraph->GetXaxis()->SetTitle("True bin");
     effgraph->GetYaxis()->SetTitle("Efficiency");
     effgraph->GetYaxis()->SetTitleOffset(1.25);

     ceffgraph->Write();
  }

  HistoDir->cd("..");

  f->Close();
  
  return;

}

//********************************************************************
void protoana::ProtoDUNEFit::BuildMeasurement(RooStats::HistFactory::Measurement& meas){
  //********************************************************************

  const int nmcchannels = _MCFileNames.size();
  for(int i=0; i < nmcchannels; i++){
    TString channelname = Form("Channel%i", i);
    RooStats::HistFactory::Channel channel(channelname.Data());
    // Add data to channel
    RooStats::HistFactory::Data data;
    for(unsigned int j=0; j < _datahistos.size(); j++){
      // Clone histogram
      TH1D* htemp = (TH1D*)(_datahistos.at(j)->Clone());

      TString hname(htemp->GetName());
      if(hname.Contains(channelname.Data()) && hname.Contains("Data")){
	mf::LogInfo("BuildMeasurement") << "Adding dataset " << hname.Data() << " to channel " << channelname.Data() << " with " << htemp->Integral() << " events";
	data.SetHisto(htemp);
	channel.SetData(data);
      }
    }
    
    // Get histograms with the systematics 
    std::vector<TH1*> systvec;
    if(_EnableSystematicError)
      systvec = protoana::ProtoDUNEFitUtils::GetSystHistograms(_SystFileNames[i]); 
    
    // Add bkg samples to channel
    for(unsigned int j=0; j < _bkghistos.size(); j++){
      // Clone histogram
      TH1D* htemp = (TH1D*)(_bkghistos.at(j)->Clone());

      TString hname(htemp->GetName());
      //if(htemp->GetEntries() == 0 || htemp->Integral() == 0) continue;
      if(hname.Contains(channelname.Data()) && hname.Contains("MC")){
	mf::LogInfo("BuildMeasurement") << "Adding MC sample " << hname.Data() << " to channel " << channelname.Data() << " with " << htemp->Integral() << " events";
	TString samplename = hname + TString("_sample");
	RooStats::HistFactory::Sample sample(samplename.Data());
	sample.SetNormalizeByTheory(true);
	
	// Check to enable statistical uncertainty
	if(_EnableStatisticalError){
	  mf::LogInfo("BuildMeasurement") << "Enable statistical uncertainty";
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
	
	//TString poiname = hname;
	//poiname.ReplaceAll("MC","POI");
	//poiname.ReplaceAll("_Histo","");
	//poiname.ReplaceAll(".","");
	//poiname.ReplaceAll("-","_");
	//poiname.ReplaceAll(channelname.Data(),"");
	//poiname.ReplaceAll("__","_");
	//meas.SetPOI(poiname.Data());
	//sample.AddNormFactor(poiname.Data(), 1.0, 0.0, 100.0);
	//mf::LogInfo("BuildMeasurement") << "Sample " << sample.GetName() << " has normalisation parameter " << poiname.Data();

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
	mf::LogInfo("BuildMeasurement") << "Adding MC sample" << hname.Data() << " to channel " << channelname.Data() << " with " << htemp->Integral() << " events";
	TString samplename = hname + TString("_sample");
	RooStats::HistFactory::Sample sample(samplename.Data());
	sample.SetNormalizeByTheory(true);
	
	// Check to enable statistical uncertainty
	if(_EnableStatisticalError){
	  mf::LogInfo("BuildMeasurement") << "Enable statistical uncertainty";
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
	poiname.ReplaceAll(".","");
	poiname.ReplaceAll("-","_");
	poiname.ReplaceAll("0000","00");
	poiname.ReplaceAll("000_","0_");
	meas.SetPOI(poiname.Data());
	sample.AddNormFactor(poiname.Data(), 1.0, 0.0, 100.0);
	mf::LogInfo("BuildMeasurement") << "Sample " << sample.GetName() << " has normalisation parameter " << poiname.Data();

	// Add sample to channel
	channel.AddSample(sample);
      }
    }

    // Statistical uncertainty less than 1% is ignored
    channel.SetStatErrorConfig(_IgnoreStatisticalErrorBelow,"Poisson"); // Poisson or Gaussian

    // Add channel to measurement
    mf::LogInfo("BuildMeasurement") <<  "Adding channel " << channel.GetName() << " to measurement " << meas.GetName();
    meas.AddChannel(channel);
  }

  // Now the sidebands
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
	mf::LogInfo("BuildMeasurement") << "Adding sideband dataset " << hname.Data() << " to channel " << channelname.Data() << " with " << htemp->Integral() << " events";
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
	mf::LogInfo("BuildMeasurement") << "Adding MC sideband sample " << hname.Data() << " to channel " << channelname.Data() << " with " <<  htemp->Integral() << " events";
	TString samplename = hname + TString("_sample");
	RooStats::HistFactory::Sample sample(samplename.Data());
	sample.SetNormalizeByTheory(true);
	
	// Check to enable statistical uncertainty
	if(_EnableStatisticalError){
	  mf::LogInfo("BuildMeasurement") << "Enable statistical uncertainty";
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
	  mf::LogInfo("BuildMeasurement") << "Sample " << sample.GetName() << " has normalisation parameter " << poiname.Data();
	}

	// Add sample to channel
	channel.AddSample(sample);
      }
    }

    // Statistical uncertainty less than 1% is ignored
    channel.SetStatErrorConfig(_IgnoreStatisticalErrorBelow,"Poisson"); // Poisson or Gaussian

    // Add channel to measurement
    mf::LogInfo("BuildMeasurement") <<  "Adding channel " << channel.GetName() << " to measurement " << meas.GetName();
    meas.AddChannel(channel);
  }

}

//********************************************************************
bool protoana::ProtoDUNEFit::FillHistogramVectors_Pions(){
  //********************************************************************

  const int nmcchannels   = _MCFileNames.size();
  const int ndatachannels = _DataFileNames.size();

  if(nmcchannels != ndatachannels){
    mf::LogError("FillHistogramVectors_Pions") << "The number of data and MC channels is not the same. Check MCFileNames and DataFileNames in the fcl file. Time to die!";
    return false;
  }

  for(int i=0; i < nmcchannels; i++){
    for(int j=2; j < 7; j++){
      _bkghistos.push_back( protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(_MCFileNames[0], _RecoTreeName, _RecoBinning, i, j) );
    }

    for(int j=0; j < 1; j++){
      for(unsigned int k=1; k < _TruthBinning.size(); k++){
	_sighistos.push_back( protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(_MCFileNames[0], _RecoTreeName, _RecoBinning, i, 1, _TruthBinning[k-1], _TruthBinning[k]) );

      }
    }
  }

  for(int i=0; i < ndatachannels; i++){
    _datahistos.push_back( protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(_DataFileNames[0], _RecoTreeName, _RecoBinning, i) );
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
  
  for(unsigned int j=0; j<_SystToConsider.size(); j++){
    TString systname(_SystToConsider[j].c_str());
    TString systtype(_SystType[j].c_str());
    
    TH1* highhist = NULL; TH1* lowhist = NULL;
    for(unsigned int i = 0; i < systvec.size(); i++){
      TH1* systhisto = (TH1*)systvec[i];
      TString systhistoname(systhisto->GetName());
      
      if(!systhistoname.Contains(hname.Data())) continue;
      if(!systhistoname.Contains(systname.Data())) continue;
      
      if(systhistoname.Contains("LOW") || systhistoname.Contains("Low") || systhistoname.Contains("low")){ 
	lowhist = systhisto;
      }
      if(systhistoname.Contains("HIGH") || systhistoname.Contains("High") || systhistoname.Contains("high")){ 
	highhist = systhisto;
      }
    }
    
    if(!highhist || !lowhist){
      std::cerr << "ERROR::Stage1: Systematic histograms not found! Will not apply systematic " << systname.Data() << " to histogram " << histo->GetName() << ". Will skip!" << std::endl;
      continue;
    } 

    TH1* highsyst = protoana::ProtoDUNEFitUtils::GetSystematicHistoFromNominal(histo, highhist, "UP");
    TH1* lowsyst  = protoana::ProtoDUNEFitUtils::GetSystematicHistoFromNominal(histo, lowhist, "DOWN");

    if(!highsyst || !lowsyst){
      std::cerr << "ERROR::Stage2: Systematic histograms not found! Will not apply systematic " << systname.Data() << " to histogram " << histo->GetName() << ". Will skip!" << std::endl;
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
  _SystFileNames               = pset.get< std::vector<std::string> >("SystFileNames");
  _SystToConsider              = pset.get< std::vector<std::string> >("SystToConsider");
  _SystType                    = pset.get< std::vector<std::string> >("SystType");
  _BackgroundTopologyName      = pset.get< std::vector<std::string> >("BackgroundTopologyName");
  
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

  _SignalTopology              = pset.get< std::vector<int> >("SignalTopology");
  _BackgroundTopology          = pset.get< std::vector<int> >("BackgroundTopology");

  return true;

}
