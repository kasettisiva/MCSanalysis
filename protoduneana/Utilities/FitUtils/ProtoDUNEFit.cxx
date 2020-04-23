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


/*namespace protoana {
enum HistType { 
 kSignal,
 kBackground,
 kIncident
};
}*/

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

  if (_DoScaleMCToData) {
    ScaleMCToData();
  }

  if (_OnlyDrawXSecs) {
    mf::LogInfo("BuildWorkspace") << "Only drawing XSecs. Exiting";

    TFile *f = new TFile(Outputfile.Data(), "RECREATE");
    TDirectory *HistoDir = f->mkdir("OriginalHistograms");
    HistoDir->cd();

    std::vector<TH1 *> pre_xsecs = DrawXSecs(0x0);
    for (size_t k = 0; k < pre_xsecs.size(); ++k) {
       pre_xsecs[k]->Write();
    }

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


    f->Close();

    return;
  }

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
  if (_AddIncidentToMeasurement)
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

  if (_DistinguishIncidentSignal) {
    for (size_t l = 1; l < _IncidentTopologyName.size(); ++l) {
      TString str = Form("%s", _IncidentTopologyName[l].c_str());
      std::cout << "NameVec " << l << " " << str << std::endl;
      incidentNameVec.push_back(str); 
    }

    for (size_t l = 1; l < _TruthBinning.size(); ++l) {
      incidentNameVec.push_back(Form("%s_%.1f-%.1f",
                                     _IncidentTopologyName[0].c_str(), _TruthBinning[l-1], _TruthBinning[l]));
      std::cout << incidentNameVec.back() << std::endl;
    }
  }
  else {
    for (size_t l = 0; l < _IncidentTopologyName.size(); ++l) {
      TString str = Form("%s", _IncidentTopologyName[l].c_str());
      std::cout << "NameVec " << l << " " << str << std::endl;
      incidentNameVec.push_back(str); 
    }
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

  if (_syst_hists.size()) {
    TDirectory * SystDir = f->mkdir("SystematicHists");
    SystDir->cd();
    
    for (size_t i = 0; i < _syst_hists.size(); ++i) {
      _syst_hists[i]->Write();
    }
    HistoDir->cd();
  }
  // Decorate and save efficiency graphs
  for(unsigned int i=0; i < _efficiencyGraphs.size(); i++){
    TGraphAsymmErrors* effgraph = (TGraphAsymmErrors*)_efficiencyGraphs.at(i);
    
    TCanvas* ceffgraph = new TCanvas(Form("ceffgraph%i",i),Form("Efficiency for channel %i",i));
    
    DecorateEfficiency( effgraph );
    effgraph->Draw("*a");

    ceffgraph->Write();
  }

  // Try drawing the xsecs
  std::vector<TH1 *> post_xsecs = DrawXSecs(fitresult, true);
  for (size_t k = 0; k < post_xsecs.size(); ++k) {
     post_xsecs[k]->Write();
  }
  std::cout << "Finished XSec" << std::endl;

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
	mf::LogInfo("AddSamplesAndChannelsToMeasurement") << "Adding dataset " <<
            hname.Data() << " to channel " << channelname.Data() << " with " <<
            htemp->Integral() << " events";

	data.SetHisto(htemp);
	channel.SetData(data);
      }
    }
    
    // Get histograms with the systematics 
    
    std::vector<TH1*> systvec;
    if (_EnableSystematicError && _UseComputedSysts) {
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

	if (_EnableSystematicError){
          if (_UseComputedSysts) {
	    ApplySystematicToSample(sample, htemp, systvec, false, false);
          }
          else {
            BuildBackgroundSystThenApplyToSample(sample, htemp, false, false,
                _bkg_chan_index[j], _bkg_topo_index[j]);
            //BuildSystThenApplyToSample(sample, htemp, false, false,
            //    kBackground, _bkg_chan_index[j], _bkg_topo_index[j]);
          }
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

	if (_EnableSystematicError){
          if (_UseComputedSysts) { //Make this a vector?
	    ApplySystematicToSample(sample, htemp, systvec, true, _NormalisedSystematic);
          }
          else { 
            BuildSignalSystThenApplyToSample(sample, htemp, true,
                _NormalisedSystematic, _sig_chan_index[j], _sig_topo_index[j],
                _sig_truth_index[j]);
          }
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
  
  std::vector<TH1*> systvec;
  if (_EnableSystematicError && _UseComputedSysts) {
    //Just use the first for now
    systvec = protoana::ProtoDUNEFitUtils::GetSystHistograms(_SystFileNames[0]);
    std::cout << "Got " << systvec.size() << " syst hists" << std::endl;
  }
    
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
    
    
    
    if (_EnableSystematicError) {
      if (_UseComputedSysts) {
        ApplySystematicToSample(sample, htemp, systvec, false, false);
      }
      else{
        BuildIncidentSystThenApplyToSample(sample, htemp, false, false,
            _inc_bkg_topo_index[j]);
      }
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

  //Needed to make the systs for the bin-by-bin incident sample
  //std::pair<std::vector<TH1 *>, std::vector<TH1 *>> inc_sig_systs;
  std::vector<std::vector<std::pair<TH1*, TH1*>>> inc_sig_systs;
  if (_DistinguishIncidentSignal && _EnableSystematicError && !_UseComputedSysts) {
     std::cout << "Incident systs for topo: " << _IncidentTopology[0] << " " <<
         _IncidentTopologyName[0] << std::endl;
     inc_sig_systs = BuildIncidentSignalSyst(0);
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
    
    
    if (_EnableSystematicError) {
      if (_UseComputedSysts) {
        ApplySystematicToSample(sample, htemp, systvec, false, false);
      }
      else{ //Move this up above and create a vector -- for each inc sig bin
        //BuildIncidentSystThenApplyToSample(sample, htemp, false, false,
        //    _inc_sig_topo_index[j]); 
        for (size_t k = 0; k < _SystToConsider.size(); ++k) {
          std::cout << "Applying syst " << _SystToConsider[k] << " to " <<
                       htemp->GetName() << " " << inc_sig_systs[k][j].first <<
                       inc_sig_systs[k][j].second << std::endl;
          ApplyBuiltSystToSample(htemp, inc_sig_systs[k][j].first,
                                 inc_sig_systs[k][j].second, sample,
                                 _SystToConsider[k], _SystType[k], false);
        }
      }
    }
    
    
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

      //For systs later
      _bkg_chan_index.push_back(i);
      _bkg_topo_index.push_back(j);

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

          //For systs later
          _sig_chan_index.push_back(i);
          _sig_topo_index.push_back(j);
          _sig_truth_index.push_back(k);

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
    _datahistos.push_back(
        protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(
            _DataFileNames[i], _RecoTreeName, _RecoBinning, _ChannelNames[i],
            _DoNegativeReco));
    //_incdatahistos.push_back( protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(_DataFileNames[0], _RecoTreeName, _RecoBinning, i, true) );
  }

  for(int i=0; i < ninctopo; i++){
    /*
    TH1* inchisto =
        protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
            _IncidentMCFileNames[0], _RecoTreeName, _RecoBinning,
            _IncidentTopologyName[i], _IncidentTopology[i], _EndZCut,
            _DoNegativeReco);

    for (size_t j = 1; j < _IncidentMCFileNames.size(); ++j) {
      inchisto->Add(
          protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
              _IncidentMCFileNames[j], _RecoTreeName, _RecoBinning,
              _IncidentTopologyName[i], _IncidentTopology[i], _EndZCut,
              _DoNegativeReco));
    }

    inchisto->SetNameTitle(Form("MC_ChannelIncident_%s_Histo",
                                _IncidentTopologyName[i].c_str()),
                           Form("Incident MC for topology %s",
                                _IncidentTopologyName[i].c_str()));


    
    if (i == 0 && _DistinguishIncidentSignal) {
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
        _inc_sig_topo_index.push_back(i);
      }   
    }
    else {
      _incbkghistos.push_back(inchisto);
      _inc_bkg_topo_index.push_back(i);
    }
  }
*/
    if (i == 0) {
      for (size_t k = 1; k < _TruthBinning.size(); ++k) {
        TH1 * inchisto =
            protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
                _IncidentMCFileNames[0], _RecoTreeName, _RecoBinning,
                _IncidentTopologyName[i], _IncidentTopology[i],
                _EndZCut, _TruthBinning[k-1], _TruthBinning[k],
                _DoNegativeReco);
        for (size_t j = 1; j < _IncidentMCFileNames.size(); ++j) {
          inchisto->Add(
              protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
                  _IncidentMCFileNames[j], _RecoTreeName, _RecoBinning,
                  _IncidentTopologyName[i], _IncidentTopology[i],
                  _EndZCut, _TruthBinning[k-1], _TruthBinning[k],
                  _DoNegativeReco));
        }
        TString hname = Form("%s_TrueBin_%.1f-%.1f", inchisto->GetName(),
            _TruthBinning[k-1], _TruthBinning[k]);
        inchisto->SetName(hname.Data());
        _incsighistos.push_back(inchisto);
        _inc_sig_topo_index.push_back(i);
      }
    }
    else {
      TH1* inchisto =
          protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
              _IncidentMCFileNames[0], _RecoTreeName, _RecoBinning,
              _IncidentTopologyName[i], _IncidentTopology[i], _EndZCut,
              0., 10000., _DoNegativeReco);

      for (size_t j = 1; j < _IncidentMCFileNames.size(); ++j) {
        inchisto->Add(
            protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
                _IncidentMCFileNames[j], _RecoTreeName, _RecoBinning,
                _IncidentTopologyName[i], _IncidentTopology[i], _EndZCut,
                0., 10000., _DoNegativeReco));
      }

      inchisto->SetNameTitle(Form("MC_ChannelIncident_%s_Histo",
                                  _IncidentTopologyName[i].c_str()),
                             Form("Incident MC for topology %s",
                                  _IncidentTopologyName[i].c_str()));

      _incbkghistos.push_back(inchisto);
      _inc_bkg_topo_index.push_back(i);
    }
  }

  TH1* incdatahisto =
      protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(
          //_DataFileNames[0], _RecoTreeName, _RecoBinning, _ChannelNames[0],
          _IncidentDataFileNames[0], _RecoTreeName, _RecoBinning, ""/*_ChannelNames[0]*/,
          _DoNegativeReco, true);

  for (size_t i = 1; i < _IncidentDataFileNames.size(); ++i) {
    incdatahisto->Add(
      protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(
          //_DataFileNames[i], _RecoTreeName, _RecoBinning, _ChannelNames[i],
          _IncidentDataFileNames[i], _RecoTreeName, _RecoBinning, ""/*_ChannelNames[i]*/,
          _DoNegativeReco, true));
  }
  _incdatahistos.push_back(incdatahisto);

  std::pair<TH1 *, TH1 *> inc_eff_num_denom = 
        protoana::ProtoDUNESelectionUtils::GetMCIncidentEfficiency(
        _IncidentMCFileNames[0], _TruthTreeName, _TruthBinning, _EndZCut);

  _incidentEfficiencyNum = inc_eff_num_denom.first;
  _incidentEfficiencyDenom = inc_eff_num_denom.second;
  for (int i = 1; i < _incidentEfficiencyDenom->GetNbinsX(); ++i) {
    if (_incidentEfficiencyDenom->GetBinContent(i) < 1.e-5)
      _incidentEfficiencyDenom->SetBinContent(i, 1.);
  }

  //_incidentEfficiency = new TGraphAsymmErrors(_incidentEfficiencyNum,
  //      _incidentEfficiencyDenom);

  //Need to build up the hists for events the pass the selection
  TH1D * incident_numerator = new TH1D("incident_numerator", "",
                                       _TruthBinning.size()-1, 0,
                                       _TruthBinning.size()-1);

  for (size_t i = 0; i < _TruthBinning.size() - 1; ++i) {
    incident_numerator->SetBinContent(i+1, _incsighistos[i]->Integral());
  }
  std::cout << "n bins: " << inc_eff_num_denom.first->GetNbinsX() << " " <<
               inc_eff_num_denom.second->GetNbinsX() << " " <<
               incident_numerator->GetNbinsX() << std::endl;
  _incidentEfficiency = new TGraphAsymmErrors(incident_numerator,
                                              _incidentEfficiencyDenom);

  _incidentEfficiency->SetNameTitle("MC_Incident_Efficiency", "Efficiency");

  std::map<std::string, std::vector<TH1 *>> signal_hists_by_topo;
  for (std::string topo : _SignalTopologyName) {
    for (size_t i = 0; i < _sighistos.size(); ++i) {
      std::string name = _sighistos[i]->GetName();
      //std::cout << name << std::endl;
      if (name.find("_" + topo + "_") == std::string::npos) continue;
      if (name.find("MC_Channel" + topo) != std::string::npos) {
        signal_hists_by_topo[topo].push_back(_sighistos[i]); 
      }
    }
  }

  //Interacting Efficiencies
  for( int i = 0; i < nsigtopo; ++i ){
    int topo = _SignalTopology[i];
    
    std::pair< TH1 *, TH1 * > eff_num_denom = 
        protoana::ProtoDUNESelectionUtils::GetMCInteractingEfficiency( 
            _IncidentMCFileNames[0], _TruthTreeName, _TruthBinning,
            _ChannelNames[i], _SignalTopologyName[i], topo, _EndZCut);

    _interactingEfficiencyNums.push_back(eff_num_denom.first); 
    for (int i = 1; i < eff_num_denom.second->GetNbinsX(); ++i) {
      if (eff_num_denom.second->GetBinContent(i) < 1.e-5)
        eff_num_denom.second->SetBinContent(i, 1.);
    }
    _interactingEfficiencyDenoms.push_back(eff_num_denom.second); 

    //TGraphAsymmErrors * eff = new TGraphAsymmErrors(eff_num_denom.first,
    //    eff_num_denom.second);

    TH1D * interacting_numerator = new TH1D("incident_numerator", "",
                                            _TruthBinning.size()-1, 0,
                                            _TruthBinning.size()-1);
    
    std::cout << "n bins: " << eff_num_denom.first->GetNbinsX() << " " <<
                 eff_num_denom.second->GetNbinsX() << " " <<
                 interacting_numerator->GetNbinsX() << std::endl;

    for (size_t j = 0; j < signal_hists_by_topo[_SignalTopologyName[i]].size(); ++j) {
      interacting_numerator->SetBinContent(
          j+1, signal_hists_by_topo[_SignalTopologyName[i]][j]->Integral());
    }
    TGraphAsymmErrors * eff = new TGraphAsymmErrors(interacting_numerator, eff_num_denom.second);

    std::string name = "MC_Channel_" + _ChannelNames[i] + "_" + 
        _SignalTopologyName[i] + "_Interacting_Efficiency";

    eff->SetNameTitle(name.c_str(), "Efficiency");

    _interactingEfficiencies.push_back(eff);
    
  }



  return true;

}

void protoana::ProtoDUNEFit::ScaleMCToData() {
  //Get the number of incident pions from MC
  double nIncidentPionsMC = 0.;
  for (size_t i = 0; i < _IncidentMCFileNames.size(); ++i) {
    //Get the tree
    TFile incidentFile(_IncidentMCFileNames[i].c_str(), "OPEN");
    TTree * incident_tree = (TTree*)incidentFile.Get(_RecoTreeName.c_str());
    nIncidentPionsMC += incident_tree->GetEntries();
    incidentFile.Close();
  }
  std::cout << "Got " << nIncidentPionsMC << " incident MC Pions" << std::endl;
  
  //Get the number of incident pions from Data
  double nIncidentPionsData = 0.;
  for (size_t i = 0; i < _IncidentDataFileNames.size(); ++i) {
    //Get the tree
    TFile incidentFile(_IncidentDataFileNames[i].c_str(), "OPEN");
    TTree * incident_tree = (TTree*)incidentFile.Get(_RecoTreeName.c_str());
    nIncidentPionsData += incident_tree->GetEntries();
    incidentFile.Close();
  }
  std::cout << "Got " << nIncidentPionsData << " incident Data Pions" << std::endl;

  double scale_factor = nIncidentPionsData/nIncidentPionsMC;
  std::cout << "Scale factor: " << scale_factor << std::endl;

  for(size_t i = 0; i < _sighistos.size(); i++){
    _sighistos[i]->Scale(scale_factor);
  }
  for(size_t i = 0; i < _bkghistos.size(); i++){
    _bkghistos[i]->Scale(scale_factor);
  }
  for(size_t i = 0; i < _incsighistos.size(); i++){
    _incsighistos[i]->Scale(scale_factor);
  }
  for(size_t i = 0; i < _incbkghistos.size(); i++){
    _incbkghistos[i]->Scale(scale_factor);
  }

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

    if (fabs(1. - lowval) < _IgnoreSystematicErrorBelow ||
        fabs(1. - highval) < _IgnoreSystematicErrorBelow) {

      std::cout << "WARNING::Ignoring systematic uncertainty " << systname
                << " for sample " << hname << " as it is below the threshold "
                << _IgnoreSystematicErrorBelow <<  " , syst fraction = "
                << 1. - lowval << " , " << highval - 1. << std::endl;

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
  _IncidentDataFileNames       = pset.get< std::vector<std::string> >("IncidentDataFileNames");
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
  _UseComputedSysts            = pset.get<bool>("UseComputedSysts");
  _NormalisedSystematic        = pset.get<bool>("NormalisedSystematic");

  _RecoBinning                 = pset.get< std::vector<double> >("RecoBinning");
  _TruthBinning                = pset.get< std::vector<double> >("TruthBinning");
  _FitInReco                   = pset.get<bool>("FitInReco");

  _SignalTopology              = pset.get< std::vector<int> >("SignalTopology");
  _BackgroundTopology          = pset.get< std::vector<int> >("BackgroundTopology");
  _IncidentTopology            = pset.get< std::vector<int> >("IncidentTopology");

  _AddIncidentToMeasurement    = pset.get<bool>("AddIncidentToMeasurement");
  _DoNegativeReco              = pset.get<bool>("DoNegativeReco"); 
  _DistinguishIncidentSignal   = pset.get<bool>("DistinguishIncidentSignal"); 
  _EndZCut                     = pset.get<double>("EndZCut");

  _DoScaleMCToData             = pset.get<bool>("DoScaleMCToData");
  _OnlyDrawXSecs               = pset.get<bool>("OnlyDrawXSecs");
  _WirePitch                   = pset.get<double>("WirePitch");

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


bool protoana::ProtoDUNEFit::BuildSignalSystThenApplyToSample(
    RooStats::HistFactory::Sample& sample, TH1* histo,
    bool hasnormfactor, bool isnorm,
    size_t iChan, size_t iTopo, size_t iTruthBin) {

  TH1 * low_hist = 0x0;
  TH1 * high_hist = 0x0;

  for (size_t i = 0; i < _SystToConsider.size(); ++i) {
    std::string syst_name = _SystToConsider[i];
    std::string syst_type = _SystType[i];

    if (iTruthBin > _TruthBinning.size() - 1) {
      std::cout << "Error! Requesting truth bin " << iTruthBin <<
                   " from binning vector of size " << 
                   _TruthBinning.size() << std::endl;
      return false;
    }

    low_hist =
        protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(
            _MCFileNames[iChan], _RecoTreeName, _RecoBinning,
            _ChannelNames[iChan], _SignalTopologyName[iTopo],
            _SignalTopology[iTopo], _TruthBinning[iTruthBin-1],
            _TruthBinning[iTruthBin], _EndZCut, false, -1, syst_name);

    high_hist =
        protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(
            _MCFileNames[iChan], _RecoTreeName, _RecoBinning,
            _ChannelNames[iChan], _SignalTopologyName[iTopo],
            _SignalTopology[iTopo], _TruthBinning[iTruthBin-1],
            _TruthBinning[iTruthBin], _EndZCut, false, 1, syst_name);

    if (!(low_hist && high_hist)) {
      continue;
    }

    std::cout << "Adding systs: " << low_hist->GetName() << " " << 
                 high_hist->GetName() << std::endl;
    _syst_hists.push_back((TH1*)low_hist->Clone()); 
    _syst_hists.push_back((TH1*)high_hist->Clone()); 

    ApplyBuiltSystToSample(histo, high_hist, low_hist, sample, syst_name,
                           syst_type, hasnormfactor);
  }

  return true;
}

bool protoana::ProtoDUNEFit::BuildBackgroundSystThenApplyToSample(
    RooStats::HistFactory::Sample& sample, TH1* histo,
    bool hasnormfactor, bool isnorm,
    size_t iChan, size_t iTopo) {

  TH1 * low_hist = 0x0;
  TH1 * high_hist = 0x0;

  for (size_t i = 0; i < _SystToConsider.size(); ++i) {
    std::string syst_name = _SystToConsider[i];
    std::string syst_type = _SystType[i];

    double tmin = _TruthBinning[0];
    double tmax = _TruthBinning.back();
    if(_FitInReco){
      tmin = 0.0;
      tmax = 1000000.0;
    }

    low_hist =
        protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(
            _MCFileNames[iChan], _RecoTreeName, _RecoBinning,
            _ChannelNames[iChan], _BackgroundTopologyName[iTopo],
            _BackgroundTopology[iTopo], _EndZCut, tmin, tmax, false, -1,
            syst_name);

    high_hist =
        protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(
            _MCFileNames[iChan], _RecoTreeName, _RecoBinning,
            _ChannelNames[iChan], _BackgroundTopologyName[iTopo],
            _BackgroundTopology[iTopo], _EndZCut, tmin, tmax, false, 1,
            syst_name);

    if (!(low_hist && high_hist)) {
      continue;
    }

    std::cout << "Adding systs: " << low_hist->GetName() << " " << 
                 high_hist->GetName() << std::endl;
    _syst_hists.push_back((TH1*)low_hist->Clone()); 
    _syst_hists.push_back((TH1*)high_hist->Clone()); 

    ApplyBuiltSystToSample(histo, high_hist, low_hist, sample, syst_name,
                           syst_type, hasnormfactor);
  }

  return true;
}

bool protoana::ProtoDUNEFit::BuildIncidentSystThenApplyToSample(
    RooStats::HistFactory::Sample& sample, TH1* histo,
    bool hasnormfactor, bool isnorm, size_t iTopo) {

  TH1 * low_hist = 0x0;
  TH1 * high_hist = 0x0;

  for (size_t i = 0; i < _SystToConsider.size(); ++i) {
    std::string syst_name = _SystToConsider[i];
    std::string syst_type = _SystType[i];

    low_hist =
        protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
            _IncidentMCFileNames[0], _RecoTreeName, _RecoBinning,
            _IncidentTopologyName[iTopo], _IncidentTopology[iTopo], 0., 10000.,
            _EndZCut, _DoNegativeReco, -1, _SystToConsider[i]);

    for (size_t j = 1; j < _IncidentMCFileNames.size(); ++j) {
      low_hist->Add(
          protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
              _IncidentMCFileNames[j], _RecoTreeName, _RecoBinning,
              _IncidentTopologyName[iTopo], _IncidentTopology[iTopo], 0., 10000.,
              _EndZCut, _DoNegativeReco, -1, _SystToConsider[i]));
    }


    high_hist =
        protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
            _IncidentMCFileNames[0], _RecoTreeName, _RecoBinning,
            _IncidentTopologyName[iTopo], _IncidentTopology[iTopo], 0., 10000.,
            _EndZCut, _DoNegativeReco, +1, _SystToConsider[i]);

    for (size_t j = 1; j < _IncidentMCFileNames.size(); ++j) {
      high_hist->Add(
          protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
              _IncidentMCFileNames[j], _RecoTreeName, _RecoBinning,
              _IncidentTopologyName[iTopo], _IncidentTopology[iTopo], 0., 10000.,
              _EndZCut, _DoNegativeReco, +1, _SystToConsider[i]));
    }

    if (!(low_hist && high_hist)) {
      continue;
    }

    std::cout << "Adding systs: " << low_hist->GetName() << " " << 
                 high_hist->GetName() << std::endl;
    _syst_hists.push_back((TH1*)low_hist->Clone()); 
    _syst_hists.push_back((TH1*)high_hist->Clone()); 

    ApplyBuiltSystToSample(histo, high_hist, low_hist, sample, syst_name,
                           syst_type, hasnormfactor);
  }

  return true;
}

std::vector<std::vector<std::pair<TH1*, TH1*>>>
//std::pair<std::vector<TH1*>, std::vector<TH1*>>
    protoana::ProtoDUNEFit::BuildIncidentSignalSyst(size_t iTopo) {
  
  //std::vector<TH1 *> low_hists, high_hists;
  std::vector<std::vector<std::pair<TH1*, TH1*>>> results;

  TH1 * low_hist = 0x0;
  TH1 * high_hist = 0x0;

  for (size_t i = 0; i < _SystToConsider.size(); ++i) {
    std::string syst_name = _SystToConsider[i];
    std::string syst_type = _SystType[i];

    low_hist =
        protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
            _IncidentMCFileNames[0], _RecoTreeName, _RecoBinning,
            _IncidentTopologyName[iTopo], _IncidentTopology[iTopo], 0., 10000.,
            _EndZCut, _DoNegativeReco, -1, _SystToConsider[i]);

    for (size_t j = 1; j < _IncidentMCFileNames.size(); ++j) {
      low_hist->Add(
          protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
              _IncidentMCFileNames[j], _RecoTreeName, _RecoBinning,
              _IncidentTopologyName[iTopo], _IncidentTopology[iTopo], 0., 10000.,
              _EndZCut, _DoNegativeReco, -1, _SystToConsider[i]));
    }


    high_hist =
        protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
            _IncidentMCFileNames[0], _RecoTreeName, _RecoBinning,
            _IncidentTopologyName[iTopo], _IncidentTopology[iTopo], 0., 10000.,
            _EndZCut, _DoNegativeReco, +1, _SystToConsider[i]);

    for (size_t j = 1; j < _IncidentMCFileNames.size(); ++j) {
      high_hist->Add(
          protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
              _IncidentMCFileNames[j], _RecoTreeName, _RecoBinning,
              _IncidentTopologyName[iTopo], _IncidentTopology[iTopo], 0., 10000.,
              _EndZCut, _DoNegativeReco, +1, _SystToConsider[i]));
    }

    if (!(low_hist && high_hist)) {
      return results;
    }
    
    //Create a new vector for this syst
    results.push_back(std::vector<std::pair<TH1*, TH1*>>());

    //Split up the histos 
    for (int j = 1; j <= low_hist->GetNbinsX(); j++) {
      TString low_name = Form("%s_IncBin%i", low_hist->GetName(),j);
      TString high_name = Form("%s_IncBin%i", high_hist->GetName(),j);

      TH1* temp_low_hist = (TH1*)low_hist->Clone();
      temp_low_hist->Reset();
      temp_low_hist->SetName(low_name.Data());

      TH1* temp_high_hist = (TH1*)high_hist->Clone();
      temp_high_hist->Reset();
      temp_high_hist->SetName(high_name.Data());

      for (int k = 1; k <= low_hist->GetNbinsX(); k++) {
        if (k == j) {
          temp_low_hist->SetBinContent(k, low_hist->GetBinContent(k));
          temp_high_hist->SetBinContent(k, high_hist->GetBinContent(k));
        }
        else {
          temp_low_hist->SetBinContent(k, 0.0);
          temp_high_hist->SetBinContent(k, 0.0);
        }
      }
      results.back().push_back({temp_high_hist, temp_low_hist});
      //low_hists.push_back(temp_low_hist);
      //high_hists.push_back(temp_high_hist);
      std::cout << "Adding systs: " << temp_low_hist->GetName() << " " << 
                   temp_high_hist->GetName() << std::endl;
      _syst_hists.push_back((TH1*)temp_low_hist->Clone()); 
      _syst_hists.push_back((TH1*)temp_high_hist->Clone()); 
    }
  }

  std::cout << "Systs: " << _SystToConsider.size() << " " << results.size() << std::endl;

  for (size_t i = 0; i < _SystToConsider.size(); ++i) {
    std::cout << i << " " << _SystToConsider[i] << " " << results[i].size() << std::endl;
    for (size_t j = 0; j < results[i].size(); ++j) {
      std::cout << results[i][j].first << " " << results[i][j].second << std::endl;
    }
  }

  return results;

}

bool protoana::ProtoDUNEFit::BuildSystThenApplyToSample(
    RooStats::HistFactory::Sample& sample, TH1* histo,
    bool hasnormfactor, bool isnorm, protoana::HistType this_histType,
    size_t iChan, size_t iTopo, size_t iTruthBin) {

  TH1 * low_hist = 0x0;
  TH1 * high_hist = 0x0;

  for (size_t i = 0; i < _SystToConsider.size(); ++i) {
    std::string syst_name = _SystToConsider[i];
    std::string syst_type = _SystType[i];

    switch (this_histType) {
      case kSignal: {
        
        if (iTruthBin > _TruthBinning.size() - 1) {
          std::cout << "Error! Requesting truth bin " << iTruthBin <<
                       " from binning vector of size " << 
                       _TruthBinning.size() << std::endl;
          return false;
        }

        low_hist =
            protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(
                _MCFileNames[iChan], _RecoTreeName, _RecoBinning,
                _ChannelNames[iChan], _SignalTopologyName[iTopo],
                _SignalTopology[iTopo], _TruthBinning[iTruthBin-1],
                _TruthBinning[iTruthBin], _EndZCut, false, -1, syst_name);

        high_hist =
            protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(
                _MCFileNames[iChan], _RecoTreeName, _RecoBinning,
                _ChannelNames[iChan], _SignalTopologyName[iTopo],
                _SignalTopology[iTopo], _TruthBinning[iTruthBin-1],
                _TruthBinning[iTruthBin], _EndZCut, false, 1, syst_name);

        break;
      }
      case kBackground: {

        double tmin = _TruthBinning[0];
        double tmax = _TruthBinning.back();
        if(_FitInReco){
          tmin = 0.0;
          tmax = 1000000.0;
        }

        low_hist =
            protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(
                _MCFileNames[iChan], _RecoTreeName, _RecoBinning,
                _ChannelNames[iChan], _BackgroundTopologyName[iTopo],
                _BackgroundTopology[iTopo], _EndZCut, tmin, tmax, false, -1,
                syst_name);

        high_hist =
            protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(
                _MCFileNames[iChan], _RecoTreeName, _RecoBinning,
                _ChannelNames[iChan], _BackgroundTopologyName[iTopo],
                _BackgroundTopology[iTopo], _EndZCut, tmin, tmax, false, 1,
                syst_name);
        break;
      }
      case kIncident: {
        break;
        /*low_hist = 
            protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
                _IncidentMCFileNames[0], _RecoTreeName, _RecoBinning,
                _ChannelNames[0], _IncidentTopologyName[iTopo], _IncidentTopology[iTopo],
                _EndZCut);*/
      }
      default: {
        return false;
      }
    }

    if (!(low_hist && high_hist)) {
      continue;
    }

    //ApplyBuiltSystToSample(histo, high_hist, low_hist, sample, syst_name,
                           //syst_type, hasnormfactor);

    double nominal_integral = histo->Integral();
    double high_integral = high_hist->Integral();
    double high_val = high_integral/nominal_integral;
    double low_integral = low_hist->Integral();
    double low_val = low_integral/nominal_integral;
    std::string hname = histo->GetName();

    if ((fabs(1. - low_val) < _IgnoreSystematicErrorBelow) ||
        (fabs(high_val - 1.) < _IgnoreSystematicErrorBelow)) {
      std::cout << "WARNING::Ignoring systematic uncertainty " << syst_name <<
                   " for sample " << hname << " as it is below the threshold " <<
                   _IgnoreSystematicErrorBelow <<  " , syst fraction = " <<
                   fabs(1. - low_val) << " , " << fabs(high_val - 1.) << std::endl;
      continue;
    }
    
    if (nominal_integral < 1.e-4 ||
        high_integral < 1.e-4 ||
        low_integral < 1.e-4 ) {
      continue;
    }

    _syst_hists.push_back((TH1*)low_hist->Clone()); 
    _syst_hists.push_back((TH1*)high_hist->Clone()); 
    

    //need this above
    if (syst_type == "NormOnly") {
      continue;
    }

    // Mixed (norm and/or shape) systematic case
    if(!hasnormfactor){ // Sample without normalisation parameter

      // Single bin samples and without norm factor are considered as norm syst
      if (protoana::ProtoDUNEFitUtils::IsSingleBinHisto(histo)) {
        sample.AddOverallSys(syst_name.c_str(), low_val, high_val);//systname from above

        std::cout << "INFO::Adding single bin and without norm parameter " <<
                     "norm systematic " << syst_name << " , to sample " << hname <<
                     " , with high value = " << high_val <<
                     " , and low value = " << low_val <<
                     " , with nominal integral = " << nominal_integral <<
                     " , high histogram integral = " << high_integral <<
                     " , low histogram integral = " << low_integral <<
                     " , and syst fraction = " << 1. - low_val <<
                     " , " << high_val - 1. << std::endl;
        continue; 
      }
      else{ // multi-bin region and without norm parameter
        RooStats::HistFactory::HistoSys syst;
        syst.SetName(syst_name.c_str());
        syst.SetHistoHigh(high_hist);
        syst.SetHistoLow(low_hist);
        sample.AddHistoSys(syst);

        std::cout << "INFO::Adding multi-bin and without norm parameter " <<
                     "shape+norm systematic " << syst_name <<
                     " , to sample " << hname <<
                     " , with nominal integral = " << nominal_integral <<
                     " , high histogram integral = " << high_integral <<
                     " , low histogram integral = " << low_integral <<
                     ", and syst fraction = " << 1. - low_val <<
                     " , " << high_val - 1. << std::endl;
        continue; 
      }
    }
    
  }

  return true;
}

bool protoana::ProtoDUNEFit::ApplyBuiltSystToSample(
    TH1 * histo, TH1 * high_hist, TH1 * low_hist,
    RooStats::HistFactory::Sample& sample,
    std::string syst_name, std::string syst_type,
    bool hasnormfactor) {

  double nominal_integral = histo->Integral();
  double high_integral = high_hist->Integral();
  double high_val = high_integral/nominal_integral;
  double low_integral = low_hist->Integral();
  double low_val = low_integral/nominal_integral;
  std::string hname = histo->GetName();

  if ((fabs(1. - low_val) < _IgnoreSystematicErrorBelow) ||
      (fabs(high_val - 1.) < _IgnoreSystematicErrorBelow)) {
    std::cout << "WARNING::Ignoring systematic uncertainty " << syst_name <<
                 " for sample " << hname << " as it is below the threshold " <<
                 _IgnoreSystematicErrorBelow <<  " , syst fraction = " <<
                 1. - low_val << " , " << high_val - 1. << std::endl;
    return false; 
  }

  if (nominal_integral < 1.e-4 ||
      high_integral < 1.e-4 ||
      low_integral < 1.e-4 ) {
    return false;
  }

  //need this above
  if (syst_type == "NormOnly") {
    return false;
  }

  // Mixed (norm and/or shape) systematic case
  if(!hasnormfactor){ // Sample without normalisation parameter

    // Single bin samples and without norm factor are considered as norm syst
    if (protoana::ProtoDUNEFitUtils::IsSingleBinHisto(histo)) {
      sample.AddOverallSys(syst_name.c_str(), low_val, high_val);//systname from above

      std::cout << "INFO::Adding single bin and without norm parameter " <<
                   "norm systematic " << syst_name << " , to sample " << hname <<
                   " , with high value = " << high_val <<
                   " , and low value = " << low_val <<
                   " , with nominal integral = " << nominal_integral <<
                   " , high histogram integral = " << high_integral <<
                   " , low histogram integral = " << low_integral <<
                   " , and syst fraction = " << 1. - low_val <<
                   " , " << high_val - 1. << std::endl;
      return true; 
    }
    else{ // multi-bin region and without norm parameter
      RooStats::HistFactory::HistoSys syst;
      syst.SetName(syst_name.c_str());
      syst.SetHistoHigh(high_hist);
      syst.SetHistoLow(low_hist);
      sample.AddHistoSys(syst);

      std::cout << "INFO::Adding multi-bin and without norm parameter " <<
                   "shape+norm systematic " << syst_name <<
                   " , to sample " << hname <<
                   " , with nominal integral = " << nominal_integral <<
                   " , high histogram integral = " << high_integral <<
                   " , low histogram integral = " << low_integral <<
                   ", and syst fraction = " << 1. - low_val <<
                   " , " << high_val - 1. << std::endl;
      return true; 
    }
  }
  
  return true;
}

std::vector<TH1 *> protoana::ProtoDUNEFit::DrawXSecs(RooFitResult *fitresult, bool doPostFit) {
  
  std::cout << "Drawing XSec" << std::endl;
  
  std::vector<TH1 *> xsecs;
  
  std::map<std::string, std::vector<double>> POI_vals;

  if (doPostFit) {
    RooArgList floatParsList = fitresult->floatParsFinal();
    TIterator* itr = floatParsList.createIterator();
    RooRealVar * var = 0x0;
    while ( (var = (RooRealVar*)itr->Next()) ) {
      std::string name = var->GetName();
      if (name.find("POI") == std::string::npos) continue;
      std::cout << var->GetName() << " " << var->getVal() << std::endl;

      bool found = false;
      for (size_t i = 0; i < _SignalTopologyName.size(); ++i) {
        if (name.find(_SignalTopologyName[i]) != std::string::npos) {
          POI_vals[_SignalTopologyName[i]].push_back(var->getVal()); 
          found = true;
          break;
        }
      }
      if (!found)
        POI_vals["INC"].push_back(var->getVal()); 
    }
  }
  
  //Make the incident hists
  
  std::string incident_name = (doPostFit ? "IncidentSignalPostFit" : "IncidentSignalPreFit");
  TH1 * inc_signal_hist = new TH1D(incident_name.c_str(), "",
                                   _TruthBinning.size()-1, 0,
                                   _TruthBinning.size()-1);
  //TH1 * inc_signal_hist_prefit = new TH1D("IncidentSignalPreFit", "",
  //                                        _TruthBinning.size()-1, 0,
  //                                        _TruthBinning.size()-1);
  for (size_t i = 0; i < _TruthBinning.size()-1; ++i) {
    std::cout << "Making bin for Truth Bin: " << _TruthBinning[i] << " " <<
                 _TruthBinning[i+1] <<  std::endl;
    std::cout << "Inc hist: " << _incsighistos[i]->GetName() << std::endl;

    //Scale by the post fit parameter value
    inc_signal_hist->SetBinContent(i+1, 
        _incsighistos[i]->Integral()*(doPostFit ? POI_vals["INC"][i] : 1.));
    //inc_signal_hist_prefit->SetBinContent(i+1, _incsighistos[i]->Integral());
  }

  //xsecs.push_back((TH1*)inc_signal_hist_prefit->Clone("IncidentSignalPreFitNoEff"));

  //Efficiency correction
  for (int i = 0; i < _incidentEfficiency->GetN(); ++i) {
    if (_incidentEfficiency->GetY()[i] < 1.e-5) continue;
    inc_signal_hist->SetBinContent(i+1,
        inc_signal_hist->GetBinContent(i+1) / _incidentEfficiency->GetY()[i]);

    //inc_signal_hist_prefit->SetBinContent(i+1,
    //    inc_signal_hist_prefit->GetBinContent(i+1) / _incidentEfficiency->GetY()[i]);
  }

  xsecs.push_back(inc_signal_hist);

  std::map<std::string, std::vector<TH1 *>> signal_hists_by_topo;
  std::map<std::string, std::vector<TH1 *>> mixed_hists_by_topo;
  for (std::string topo : _SignalTopologyName) {
    for (size_t i = 0; i < _sighistos.size(); ++i) {
      std::string name = _sighistos[i]->GetName();
      if (name.find("_" + topo + "_") == std::string::npos) continue;
      if (name.find("MC_Channel" + topo) == std::string::npos) {
        mixed_hists_by_topo[topo].push_back(_sighistos[i]); 
      }
      else {
        signal_hists_by_topo[topo].push_back(_sighistos[i]); 
      }
    }
  }

  
  for (std::string topo : _SignalTopologyName) {
    std::string signal_name = "Signal" + topo +
                              (doPostFit ? "PostFit" : "PreFit");
    TH1 * signal_hist = new TH1D(signal_name.c_str(), "",
                                 inc_signal_hist->GetNbinsX(), 0,
                                 inc_signal_hist->GetNbinsX());
    //TH1 * signal_hist_prefit = new TH1D(("PreFitSignal" + topo).c_str(), "",
    //                                    inc_signal_hist->GetNbinsX(), 0,
    //                                    inc_signal_hist->GetNbinsX());

    TGraphAsymmErrors * eff = 0x0;
    for (size_t i = 0; i < _interactingEfficiencies.size(); ++i) {
      std::string name = _interactingEfficiencies[i]->GetName();
      if (name.find(topo) != std::string::npos) {
        eff = _interactingEfficiencies[i];
        break;
      }
    }

    for (size_t i = 0; i < signal_hists_by_topo[topo].size(); ++i) {
      //correct by eff
      double from_signal = signal_hists_by_topo[topo][i]->Integral();
      from_signal = from_signal / (eff->GetY()[i] > 1.e-5 ? eff->GetY()[i] : 1.);

      //double from_mixed = mixed_hists_by_topo[topo][i]->Integral();

      double combined = (from_signal /*+ from_mixed*/) /*/ eff->GetY()[i]*/;

      signal_hist->SetBinContent(i+1, (doPostFit ? POI_vals[topo][i] : 1.)*combined);
      //signal_hist_prefit->SetBinContent(i+1, combined);
    }

    xsecs.push_back(signal_hist);
    //xsecs.push_back(signal_hist_prefit);


    //Now divide the interacting hist by the incident hist
    std::string xsec_name = "XSEC_" + topo +
                            (doPostFit ? "_PostFit" : "_PreFit");
    TH1 * xsec = (TH1*)signal_hist->Clone(xsec_name.c_str());
    xsec->Sumw2();
    inc_signal_hist->Sumw2();
    xsec->Divide(inc_signal_hist);
    xsec->Scale(1.E27/ (_WirePitch * 1.390 * 6.022E23 / 39.95 ));

    for (size_t i = 1; i < _TruthBinning.size(); ++i) {
      TString ibinstr = Form("%.1f-%.1f",_TruthBinning[i-1],_TruthBinning[i]);
      xsec->GetXaxis()->SetBinLabel(i, ibinstr);
    }
    
    //TCanvas * xsec_canvas = new TCanvas(("Canvas_XSec_"+chan).c_str(), "xsec", 500, 400);
    xsec->SetMinimum(0.);
    xsec->SetMarkerColor(1);
    xsec->SetLineColor(1);
    xsec->SetMarkerStyle(20);
    xsec->SetMarkerSize(.75);
    xsec->SetLineWidth(2);
    xsecs.push_back(xsec);
 
    //PREFIT
    //Now divide the interacting hist by the incident hist
    /*
    TH1 * xsec_prefit = (TH1*)signal_hist_prefit->Clone(
        ("PreFit_XSEC_"+topo).c_str());
    xsec_prefit->Sumw2();
    inc_signal_hist_prefit->Sumw2();
    xsec_prefit->Divide(inc_signal_hist_prefit);
    xsec_prefit->Scale(1.E27/ (.4792 * 1.390 * 6.022E23 / 39.95 ));

    for (size_t i = 1; i < _RecoBinning.size(); ++i) {
      TString ibinstr = Form("%.1f-%.1f",_RecoBinning[i-1],_RecoBinning[i]);
      xsec_prefit->GetXaxis()->SetBinLabel(i, ibinstr);
    }
    
    //TCanvas * xsec_canvas = new TCanvas(("Canvas_XSec_"+chan).c_str(), "xsec", 500, 400);
    xsec_prefit->SetMinimum(0.);
    xsec_prefit->SetMarkerColor(1);
    xsec_prefit->SetLineColor(1);
    xsec_prefit->SetMarkerStyle(20);
    xsec_prefit->SetMarkerSize(.75);
    xsec_prefit->SetLineWidth(2);
    xsecs.push_back(xsec_prefit);
    */
  }
  
  return xsecs;
  //Get the signal hists to determine XSec
  /*
  for (std::string chan : {"ABS", "CEX"}) {

    TH1 * signal_hist = new TH1D(("Signal" + chan).c_str(), "", inc_signal_hist->GetNbinsX(), 0, inc_signal_hist->GetNbinsX());
    TH1 * signal_hist_prefit = new TH1D(("PreFitSignal" + chan).c_str(), "", inc_signal_hist->GetNbinsX(), 0, inc_signal_hist->GetNbinsX());

    for (size_t i = 0; i < _sighistos.size(); ++i ) {
      std::string name = _sighistos[i]->GetName();

      if (name.find("MC_Channel" + chan) == std::string::npos) continue;
      
      TH1 * temp = (TH1*)_sighistos[i]->Clone();
      TH1 * temp_prefit = (TH1*)temp->Clone();

      //Scale the bin by the post fit parameter value
      temp->SetBinContent(_sig_truth_index[i],
          temp->GetBinContent(_sig_truth_index[i])*
          POI_vals[chan][_sig_truth_index[i]-1]);

      total_hist->Add(temp);
      total_hist_prefit->Add(temp_prefit);

      auto find_signal = name.find("MC_Channel" + chan + "_" + chan);
      std::cout << _sig_truth_index[i] << " " << _sighistos[i]->GetName() 
                << " " << find_signal << std::endl;

      if (find_signal != std::string::npos) {
        signal_hist->SetBinContent(_sig_truth_index[i],
                                   temp->GetBinContent(_sig_truth_index[i]));
        signal_hist_prefit->SetBinContent(_sig_truth_index[i],
                                          temp_prefit->GetBinContent(
                                              _sig_truth_index[i]));
      }

    }

    for (size_t i = 0; i < _bkghistos.size(); ++i) {
      std::string name = _bkghistos[i]->GetName();
      auto find_channel = name.find("MC_Channel" + chan);
      if (find_channel != std::string::npos) {
        total_hist->Add(_bkghistos[i]);
        total_hist_prefit->Add(_bkghistos[i]);
      }
    }

    //Make correction
    TH1 * correction_hist = (TH1*)signal_hist->Clone(("Correction" + chan).c_str());
    TH1 * correction_hist_prefit = (TH1*)signal_hist_prefit->Clone(("PreFitCorrection" + chan).c_str());
    correction_hist->Divide(total_hist);
    correction_hist_prefit->Divide(total_hist_prefit);

    TH1 * data_hist_chan = 0x0;
    for (size_t i = 0; i < _datahistos.size(); ++i) {
      std::string name = _datahistos[i]->GetName();
      if (name.find("Data_Channel" + chan) != std::string::npos) {
        data_hist_chan = (TH1*)_datahistos[i]->Clone();
        break;
      }
    }
    data_hist_chan->Multiply(correction_hist);

    TH1 * data_hist_chan_prefit = 0x0;
    for (size_t i = 0; i < _datahistos.size(); ++i) {
      std::string name = _datahistos[i]->GetName();
      if (name.find("Data_Channel" + chan) != std::string::npos) {
        data_hist_chan_prefit = (TH1*)_datahistos[i]->Clone();
        break;
      }
    }
    data_hist_chan_prefit->Multiply(correction_hist_prefit);

    for (size_t i = 0; i < _interactingEfficiencies.size(); ++i) {
      std::string name = _interactingEfficiencies[i]->GetName();
      if (name.find(chan) != std::string::npos) {
        TGraphAsymmErrors * eff = _interactingEfficiencies[i];

        for (int j = 0; j < eff->GetN(); ++j) {
          data_hist_chan->SetBinContent(j+1,
              data_hist_chan->GetBinContent(j+1)/eff->GetY()[j]);

          data_hist_chan_prefit->SetBinContent(j+1,
              data_hist_chan_prefit->GetBinContent(j+1)/eff->GetY()[j]);
        }
        break;
      }
    }
    
    //Now divide the interacting hist by the incident hist
    TH1 * xsec = (TH1*)data_hist_chan->Clone(("XSEC_"+chan).c_str());
    xsec->Sumw2();
    data_hist_inc->Sumw2();
    xsec->Divide(data_hist_inc);
    xsec->Scale(1.E27/ (.4792 * 1.390 * 6.022E23 / 39.95 ));

    for (size_t i = 1; i < _RecoBinning.size(); ++i) {
      TString ibinstr = Form("%.1f-%.1f",_RecoBinning[i-1],_RecoBinning[i]);
      xsec->GetXaxis()->SetBinLabel(i, ibinstr);
    }
    
    //TCanvas * xsec_canvas = new TCanvas(("Canvas_XSec_"+chan).c_str(), "xsec", 500, 400);
    xsec->SetMinimum(0.);
    xsec->SetMarkerColor(1);
    xsec->SetLineColor(1);
    xsec->SetMarkerStyle(20);
    xsec->SetMarkerSize(.75);
    xsec->SetLineWidth(2);
    xsecs.push_back(xsec);
 
    //PREFIT
    //Now divide the interacting hist by the incident hist
    TH1 * xsec_prefit = (TH1*)data_hist_chan_prefit->Clone(
        ("PreFit_XSEC_"+chan).c_str());
    xsec_prefit->Sumw2();
    data_hist_inc_prefit->Sumw2();
    xsec_prefit->Divide(data_hist_inc_prefit);
    xsec_prefit->Scale(1.E27/ (.4792 * 1.390 * 6.022E23 / 39.95 ));

    for (size_t i = 1; i < _RecoBinning.size(); ++i) {
      TString ibinstr = Form("%.1f-%.1f",_RecoBinning[i-1],_RecoBinning[i]);
      xsec_prefit->GetXaxis()->SetBinLabel(i, ibinstr);
    }
    
    //TCanvas * xsec_canvas = new TCanvas(("Canvas_XSec_"+chan).c_str(), "xsec", 500, 400);
    xsec_prefit->SetMinimum(0.);
    xsec_prefit->SetMarkerColor(1);
    xsec_prefit->SetLineColor(1);
    xsec_prefit->SetMarkerStyle(20);
    xsec_prefit->SetMarkerSize(.75);
    xsec_prefit->SetLineWidth(2);
    xsecs.push_back(xsec_prefit);

  }

  return xsecs;
  */
  /*
  if (_DistinguishIncidentSignal) {
    inc_signal_hist = (TH1*)_incsighistos[0]->Clone();
    for (size_t i = 1; i < _incsighistos.size(); ++i) {
      inc_signal_hist->Add(_incsighistos[i]);
    }

    inc_signal_hist_prefit = (TH1*)_incsighistos[0]->Clone();
    for (size_t i = 1; i < _incsighistos.size(); ++i) {
      inc_signal_hist_prefit->Add(_incsighistos[i]);
    }
  }
  else {
    inc_signal_hist = (TH1*)_incbkghistos[0]->Clone();
    inc_signal_hist_prefit = (TH1*)_incbkghistos[0]->Clone();
  }
  

  TH1 * total_incident_hist = new TH1D("Total_Incident", "",
                                       inc_signal_hist->GetNbinsX(), 0,
                                       inc_signal_hist->GetNbinsX());

  TH1 * total_incident_hist_prefit = new TH1D("Total_Incident_PreFit", "",
                                              inc_signal_hist->GetNbinsX(), 0,
                                              inc_signal_hist->GetNbinsX());

  for (size_t i = 0; i < _incbkghistos.size(); ++i) {
    total_incident_hist->Add(_incbkghistos[i]);
    total_incident_hist_prefit->Add(_incbkghistos[i]);
  }
  if (_DistinguishIncidentSignal) {
    total_incident_hist->Add(inc_signal_hist);
    total_incident_hist_prefit->Add(inc_signal_hist_prefit);
  }

  //Get Data incident hist
  TH1 * data_hist_inc = (TH1*)_incdatahistos[0]->Clone();
  TH1 * data_hist_inc_prefit = (TH1*)_incdatahistos[0]->Clone();
  for (size_t i = 1; i < _incdatahistos.size(); ++i) {
    data_hist_inc->Add(_incdatahistos[i]);
    data_hist_inc_prefit->Add(_incdatahistos[i]);
  }


  //Make correction
  TH1 * correction_incident_hist = (TH1*)inc_signal_hist->Clone("CorrectionIncident");
  correction_incident_hist->Divide(total_incident_hist);

  TH1 * correction_incident_hist_prefit = (TH1*)inc_signal_hist_prefit->Clone("CorrectionIncident");
  correction_incident_hist_prefit->Divide(total_incident_hist_prefit);

  data_hist_inc->Multiply(correction_incident_hist);
  data_hist_inc_prefit->Multiply(correction_incident_hist_prefit);

  for (int i = 0; i < _incidentEfficiency->GetN(); ++i) {
    data_hist_inc->SetBinContent(i+1,
        data_hist_inc->GetBinContent(i+1) / _incidentEfficiency->GetY()[i]);

    data_hist_inc_prefit->SetBinContent(i+1,
        data_hist_inc_prefit->GetBinContent(i+1) / _incidentEfficiency->GetY()[i]);

    if (_DistinguishIncidentSignal)
      data_hist_inc->SetBinContent(i+1,
          data_hist_inc->GetBinContent(i+1) * POI_vals["INC"][i]);
  }

  //Get the signal hists to determine
  for (std::string chan : {"ABS", "CEX"}) {


    TH1 * signal_hist = new TH1D(("Signal" + chan).c_str(), "", inc_signal_hist->GetNbinsX(), 0, inc_signal_hist->GetNbinsX());
    TH1 * total_hist = new TH1D(("Total" + chan).c_str(), "", inc_signal_hist->GetNbinsX(), 0, inc_signal_hist->GetNbinsX());

    TH1 * signal_hist_prefit = new TH1D(("PreFitSignal" + chan).c_str(), "", inc_signal_hist->GetNbinsX(), 0, inc_signal_hist->GetNbinsX());
    TH1 * total_hist_prefit = new TH1D(("PreFitTotal" + chan).c_str(), "", inc_signal_hist->GetNbinsX(), 0, inc_signal_hist->GetNbinsX());

    for (size_t i = 0; i < _sighistos.size(); ++i ) {
      std::string name = _sighistos[i]->GetName();

      if (name.find("MC_Channel" + chan) == std::string::npos) continue;
      
      TH1 * temp = (TH1*)_sighistos[i]->Clone();
      TH1 * temp_prefit = (TH1*)temp->Clone();

      //Scale the bin by the post fit parameter value
      temp->SetBinContent(_sig_truth_index[i],
          temp->GetBinContent(_sig_truth_index[i])*
          POI_vals[chan][_sig_truth_index[i]-1]);

      total_hist->Add(temp);
      total_hist_prefit->Add(temp_prefit);

      auto find_signal = name.find("MC_Channel" + chan + "_" + chan);
      std::cout << _sig_truth_index[i] << " " << _sighistos[i]->GetName() 
                << " " << find_signal << std::endl;

      if (find_signal != std::string::npos) {
        signal_hist->SetBinContent(_sig_truth_index[i],
                                   temp->GetBinContent(_sig_truth_index[i]));
        signal_hist_prefit->SetBinContent(_sig_truth_index[i],
                                          temp_prefit->GetBinContent(
                                              _sig_truth_index[i]));
      }

    }

    for (size_t i = 0; i < _bkghistos.size(); ++i) {
      std::string name = _bkghistos[i]->GetName();
      auto find_channel = name.find("MC_Channel" + chan);
      if (find_channel != std::string::npos) {
        total_hist->Add(_bkghistos[i]);
        total_hist_prefit->Add(_bkghistos[i]);
      }
    }

    //Make correction
    TH1 * correction_hist = (TH1*)signal_hist->Clone(("Correction" + chan).c_str());
    TH1 * correction_hist_prefit = (TH1*)signal_hist_prefit->Clone(("PreFitCorrection" + chan).c_str());
    correction_hist->Divide(total_hist);
    correction_hist_prefit->Divide(total_hist_prefit);

    TH1 * data_hist_chan = 0x0;
    for (size_t i = 0; i < _datahistos.size(); ++i) {
      std::string name = _datahistos[i]->GetName();
      if (name.find("Data_Channel" + chan) != std::string::npos) {
        data_hist_chan = (TH1*)_datahistos[i]->Clone();
        break;
      }
    }
    data_hist_chan->Multiply(correction_hist);

    TH1 * data_hist_chan_prefit = 0x0;
    for (size_t i = 0; i < _datahistos.size(); ++i) {
      std::string name = _datahistos[i]->GetName();
      if (name.find("Data_Channel" + chan) != std::string::npos) {
        data_hist_chan_prefit = (TH1*)_datahistos[i]->Clone();
        break;
      }
    }
    data_hist_chan_prefit->Multiply(correction_hist_prefit);

    for (size_t i = 0; i < _interactingEfficiencies.size(); ++i) {
      std::string name = _interactingEfficiencies[i]->GetName();
      if (name.find(chan) != std::string::npos) {
        TGraphAsymmErrors * eff = _interactingEfficiencies[i];

        for (int j = 0; j < eff->GetN(); ++j) {
          data_hist_chan->SetBinContent(j+1,
              data_hist_chan->GetBinContent(j+1)/eff->GetY()[j]);

          data_hist_chan_prefit->SetBinContent(j+1,
              data_hist_chan_prefit->GetBinContent(j+1)/eff->GetY()[j]);
        }
        break;
      }
    }
    
    //Now divide the interacting hist by the incident hist
    TH1 * xsec = (TH1*)data_hist_chan->Clone(("XSEC_"+chan).c_str());
    xsec->Sumw2();
    data_hist_inc->Sumw2();
    xsec->Divide(data_hist_inc);
    xsec->Scale(1.E27/ (.4792 * 1.390 * 6.022E23 / 39.95 ));

    for (size_t i = 1; i < _RecoBinning.size(); ++i) {
      TString ibinstr = Form("%.1f-%.1f",_RecoBinning[i-1],_RecoBinning[i]);
      xsec->GetXaxis()->SetBinLabel(i, ibinstr);
    }
    
    //TCanvas * xsec_canvas = new TCanvas(("Canvas_XSec_"+chan).c_str(), "xsec", 500, 400);
    xsec->SetMinimum(0.);
    xsec->SetMarkerColor(1);
    xsec->SetLineColor(1);
    xsec->SetMarkerStyle(20);
    xsec->SetMarkerSize(.75);
    xsec->SetLineWidth(2);
    xsecs.push_back(xsec);
 
    //PREFIT
    //Now divide the interacting hist by the incident hist
    TH1 * xsec_prefit = (TH1*)data_hist_chan_prefit->Clone(
        ("PreFit_XSEC_"+chan).c_str());
    xsec_prefit->Sumw2();
    data_hist_inc_prefit->Sumw2();
    xsec_prefit->Divide(data_hist_inc_prefit);
    xsec_prefit->Scale(1.E27/ (.4792 * 1.390 * 6.022E23 / 39.95 ));

    for (size_t i = 1; i < _RecoBinning.size(); ++i) {
      TString ibinstr = Form("%.1f-%.1f",_RecoBinning[i-1],_RecoBinning[i]);
      xsec_prefit->GetXaxis()->SetBinLabel(i, ibinstr);
    }
    
    //TCanvas * xsec_canvas = new TCanvas(("Canvas_XSec_"+chan).c_str(), "xsec", 500, 400);
    xsec_prefit->SetMinimum(0.);
    xsec_prefit->SetMarkerColor(1);
    xsec_prefit->SetLineColor(1);
    xsec_prefit->SetMarkerStyle(20);
    xsec_prefit->SetMarkerSize(.75);
    xsec_prefit->SetLineWidth(2);
    xsecs.push_back(xsec_prefit);

  }

  return xsecs;
*/
}
