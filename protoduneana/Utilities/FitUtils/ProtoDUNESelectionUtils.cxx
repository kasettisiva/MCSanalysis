#include "ProtoDUNESelectionUtils.h"

#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "messagefacility/MessageLogger/MessageLogger.h"

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, int channel, int topo, bool IsIncidentHisto, double weight){
  //********************************************************************

  TFile *file = new TFile(filename.c_str(), "READ");
  TTree *defaultTree  = (TTree*)file->Get(treename.c_str());

  Int_t reco_beam_type; // 13 -> track-like, 11 -> shower-like
  Int_t reco_beam_nTrackDaughters, reco_beam_nShowerDaughters;
  Double_t reco_beam_len, reco_beam_vtxX, reco_beam_vtxY, reco_beam_vtxZ, reco_beam_startX, reco_beam_startY, reco_beam_startZ, reco_beam_trackDirZ, reco_beam_interactingEnergy, reco_beam_Chi2_proton;
  bool reco_beam_true_byHits_matched; // Does the true particle contributing most to the reconstructed beam track coincide with the actual beam particle that generated the event
  Int_t reco_beam_true_byHits_origin, reco_beam_true_byHits_PDG;  // Origin and PDG of the reconstructed beam track
  Int_t true_beam_PDG;
  Double_t true_beam_interactingEnergy;
  std::string *true_beam_endProcess = 0;
  std::string *reco_beam_true_byHits_endProcess = 0;
  std::vector<double> *reco_beam_incidentEnergies = 0;

  defaultTree->SetBranchAddress("reco_beam_type",                   &reco_beam_type);
  defaultTree->SetBranchAddress("reco_beam_len",                    &reco_beam_len);
  defaultTree->SetBranchAddress("reco_beam_vtxX",                   &reco_beam_vtxX);
  defaultTree->SetBranchAddress("reco_beam_vtxY",                   &reco_beam_vtxY);
  defaultTree->SetBranchAddress("reco_beam_vtxZ",                   &reco_beam_vtxZ);
  defaultTree->SetBranchAddress("reco_beam_startX",                 &reco_beam_startX);
  defaultTree->SetBranchAddress("reco_beam_startY",                 &reco_beam_startY);
  defaultTree->SetBranchAddress("reco_beam_startZ",                 &reco_beam_startZ);
  defaultTree->SetBranchAddress("reco_beam_trackDirZ",              &reco_beam_trackDirZ);
  defaultTree->SetBranchAddress("reco_beam_nTrackDaughters",        &reco_beam_nTrackDaughters);
  defaultTree->SetBranchAddress("reco_beam_nShowerDaughters",       &reco_beam_nShowerDaughters);
  defaultTree->SetBranchAddress("reco_beam_interactingEnergy",      &reco_beam_interactingEnergy);
  defaultTree->SetBranchAddress("reco_beam_Chi2_proton",            &reco_beam_Chi2_proton);
  defaultTree->SetBranchAddress("reco_beam_incidentEnergies",       &reco_beam_incidentEnergies);

  defaultTree->SetBranchAddress("reco_beam_true_byHits_matched",    &reco_beam_true_byHits_matched);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_origin",     &reco_beam_true_byHits_origin);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_PDG",        &reco_beam_true_byHits_PDG);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_endProcess", &reco_beam_true_byHits_endProcess);

  defaultTree->SetBranchAddress("true_beam_interactingEnergy",      &true_beam_interactingEnergy);
  defaultTree->SetBranchAddress("true_beam_PDG",                    &true_beam_PDG);
  defaultTree->SetBranchAddress("true_beam_endProcess",             &true_beam_endProcess);

  //defaultTree->SetBranchAddress("",        &);

  const int nrecobins = recoBins.size();
  TH1D* mchisto = new TH1D(Form("MC_Channel%i_BkgTopo%i_Histo",channel,topo), Form("MC Background for channel %i and topology %i", channel,topo), nrecobins-1, 0, nrecobins-1);
  if(IsIncidentHisto)
    mchisto->SetNameTitle( Form("MC_ChannelIncident%i_BkgTopo%i_Histo",channel,topo), Form("Incident MC Background for channel %i and topology %i", channel,topo) );
  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCBackgroundHistogram_Pions") << "Filling MC background histogram " << mchisto->GetName() << " from file " << filename.c_str() << " for channel " << channel << " with topology " << topo;

  for(Int_t k=0; k < defaultTree->GetEntries(); k++){
    defaultTree->GetEntry(k);

    // Temporary selection of events. To speed up the fit the selection should be done in advance.

    // Different background topologies - temporary: cosmics: 2, 
    Int_t topology = -1;
    std::string reco_beam_true_byHits_endProcess_str = *reco_beam_true_byHits_endProcess;

    if(reco_beam_true_byHits_origin == 2){ // cosmic origin
      topology = 2; // cosmics background
    }
    else if(reco_beam_true_byHits_origin == 4){ // beam origin
      if(true_beam_PDG != 211) continue; // look only for pion beam events

      if(abs(reco_beam_true_byHits_PDG) == 2212) topology = 3; // protons
      else if(abs(reco_beam_true_byHits_PDG) == 13) topology = 4; // muons
      else if(reco_beam_true_byHits_PDG == 211){ // background pions
	if(reco_beam_true_byHits_endProcess_str == "pi+Inelastic") continue; // this is signal - skip
	topology = 5; // other pion
      }
      else{	
	topology = 6; // other backgrounds
      }
    }
    else{
      topology = 6; // other backgrounds
    }

    // Select only the correct topology
    if(topology != topo) continue;
    
    // Track-like beam particles
    if(reco_beam_type != 13) continue;
    
    // Position cuts
    if(reco_beam_startZ < 28.0  || reco_beam_startZ > 32.0) continue;
    if(reco_beam_startY < 410.0 || reco_beam_startY > 445.0) continue;
    if(reco_beam_startX < -45.0 || reco_beam_startX > 0.0) continue;

    // Direction cuts
    if(reco_beam_trackDirZ < 0.9) continue;
    
    // Remove secondary protons from pion interacting outside the TPC
    if(reco_beam_Chi2_proton < 450.0) continue;

    // Minimum and maximum track length
    if(reco_beam_len < 5.0 || reco_beam_len > 270.0) continue;

    // Sometimes the reco energy at vertex is mis-reconstructed
    if(reco_beam_interactingEnergy < 0.0) continue;

    // For inelastic scattering at least one daughter is expected
    Int_t allDaughters = reco_beam_nTrackDaughters + reco_beam_nShowerDaughters;
    if(allDaughters <= 0) continue;

    if(IsIncidentHisto){
      for(UInt_t l = 0; l < reco_beam_incidentEnergies->size(); l++){
	for(Int_t m = 1; m <= nrecobins; m++){
	  if(reco_beam_incidentEnergies->at(l) > recoBins[m-1] && reco_beam_incidentEnergies->at(l) <= recoBins[m]){
	    mchisto->SetBinContent(m, mchisto->GetBinContent(m) + weight);
	    break;
	  }
	}
      }
      continue;
    }

    // Fill histogram in reco energy
    for(Int_t l = 1; l <= nrecobins; l++){
      if(reco_beam_interactingEnergy > recoBins[l-1] && reco_beam_interactingEnergy <= recoBins[l]){
	mchisto->SetBinContent(l, mchisto->GetBinContent(l) + weight);
	break;
      }
    }
  }

  file->Close();

  return mchisto;

}

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, int channel, int topo, double minval, double maxval, bool IsIncidentHisto, double weight){
  //********************************************************************

  TFile *file = new TFile(filename.c_str(), "READ");
  TTree *defaultTree  = (TTree*)file->Get(treename.c_str());

  Int_t reco_beam_type; // 13 -> track-like, 11 -> shower-like
  Int_t reco_beam_nTrackDaughters, reco_beam_nShowerDaughters;
  Double_t reco_beam_len, reco_beam_vtxX, reco_beam_vtxY, reco_beam_vtxZ, reco_beam_startX, reco_beam_startY, reco_beam_startZ, reco_beam_trackDirZ, reco_beam_interactingEnergy, reco_beam_Chi2_proton;
  bool reco_beam_true_byHits_matched; // Does the true particle contributing most to the reconstructed beam track coincide with the actual beam particle that generated the event
  Int_t reco_beam_true_byHits_origin, reco_beam_true_byHits_PDG;  // Origin and PDG of the reconstructed beam track
  Int_t true_beam_PDG;
  Double_t true_beam_interactingEnergy;
  std::string *true_beam_endProcess = 0;
  std::string *reco_beam_true_byHits_endProcess = 0;
  std::vector<double> *reco_beam_incidentEnergies = 0;

  defaultTree->SetBranchAddress("reco_beam_type",                   &reco_beam_type);
  defaultTree->SetBranchAddress("reco_beam_len",                    &reco_beam_len);
  defaultTree->SetBranchAddress("reco_beam_vtxX",                   &reco_beam_vtxX);
  defaultTree->SetBranchAddress("reco_beam_vtxY",                   &reco_beam_vtxY);
  defaultTree->SetBranchAddress("reco_beam_vtxZ",                   &reco_beam_vtxZ);
  defaultTree->SetBranchAddress("reco_beam_startX",                 &reco_beam_startX);
  defaultTree->SetBranchAddress("reco_beam_startY",                 &reco_beam_startY);
  defaultTree->SetBranchAddress("reco_beam_startZ",                 &reco_beam_startZ);
  defaultTree->SetBranchAddress("reco_beam_trackDirZ",              &reco_beam_trackDirZ);
  defaultTree->SetBranchAddress("reco_beam_nTrackDaughters",        &reco_beam_nTrackDaughters);
  defaultTree->SetBranchAddress("reco_beam_nShowerDaughters",       &reco_beam_nShowerDaughters);
  defaultTree->SetBranchAddress("reco_beam_interactingEnergy",      &reco_beam_interactingEnergy);
  defaultTree->SetBranchAddress("reco_beam_Chi2_proton",            &reco_beam_Chi2_proton);
  defaultTree->SetBranchAddress("reco_beam_incidentEnergies",       &reco_beam_incidentEnergies);

  defaultTree->SetBranchAddress("reco_beam_true_byHits_matched",    &reco_beam_true_byHits_matched);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_origin",     &reco_beam_true_byHits_origin);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_PDG",        &reco_beam_true_byHits_PDG);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_endProcess", &reco_beam_true_byHits_endProcess);

  defaultTree->SetBranchAddress("true_beam_interactingEnergy",      &true_beam_interactingEnergy);
  defaultTree->SetBranchAddress("true_beam_PDG",                    &true_beam_PDG);
  defaultTree->SetBranchAddress("true_beam_endProcess",             &true_beam_endProcess);

  const int nrecobins = recoBins.size();
  TH1D* mchisto = new TH1D(Form("MC_Channel%i_SigTopo%i_%.1f-%.1f_Histo",channel,topo,minval,maxval), Form("MC Signal for channel %i and topology %i and true region %.1f-%.1f", channel,topo,minval,maxval), nrecobins-1, 0, nrecobins-1);
  if(IsIncidentHisto)
    mchisto->SetNameTitle( Form("MC_ChannelIncident%i_SigTopo%i_Histo",channel,topo), Form("Incident MC Signal for channel %i and topology %i", channel,topo) );
  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCSignalHistogram_Pions") << "Filling MC signal histogram " << mchisto->GetName() << " from file " << filename.c_str() << " for channel " << channel << " with topology " << topo << " in the true region " << minval << "-" << maxval;

  for(Int_t k=0; k < defaultTree->GetEntries(); k++){
    defaultTree->GetEntry(k);

    // Temporary selection of events. To speed up the fit the selection should be done in advance.
    std::string reco_beam_true_byHits_endProcess_str = *reco_beam_true_byHits_endProcess;
    
    // Selection of signal beam pions
    if(reco_beam_true_byHits_origin != 4) continue;
    if(true_beam_PDG != 211) continue;
    // Secondary tracks non-pion coming from beam pions interacting outsided the TPC are considered as background - should be revised
    if(reco_beam_true_byHits_PDG != 211) continue;
    // Interaction should be pi+ inelastic
    if(reco_beam_true_byHits_endProcess_str != "pi+Inelastic") continue;
    
    // Track-like beam particles
    if(reco_beam_type != 13) continue;
    
    // Position cuts
    if(reco_beam_startZ < 28.0  || reco_beam_startZ > 32.0) continue;
    if(reco_beam_startY < 410.0 || reco_beam_startY > 445.0) continue;
    if(reco_beam_startX < -45.0 || reco_beam_startX > 0.0) continue;

    // Direction cuts
    if(reco_beam_trackDirZ < 0.9) continue;
    
    // Remove secondary protons from pion interacting outside the TPC
    if(reco_beam_Chi2_proton < 450.0) continue;

    // Minimum and maximum track length
    if(reco_beam_len < 5.0 || reco_beam_len > 270.0) continue;

    // Sometimes the reco energy at vertex is mis-reconstructed
    if(reco_beam_interactingEnergy < 0.0) continue;

    // For inelastic scattering at least one daughter is expected
    Int_t allDaughters = reco_beam_nTrackDaughters + reco_beam_nShowerDaughters;
    if(allDaughters <= 0) continue;

    if(IsIncidentHisto){
      for(UInt_t l = 0; l < reco_beam_incidentEnergies->size(); l++){
	for(Int_t m = 1; m <= nrecobins; m++){
	  if(reco_beam_incidentEnergies->at(l) > recoBins[m-1] && reco_beam_incidentEnergies->at(l) <= recoBins[m]){
	    mchisto->SetBinContent(m, mchisto->GetBinContent(m) + weight);
	    break;
	  }
	}
      }
      continue;
    }

    // True energy bin
    if(true_beam_interactingEnergy < minval) continue;
    if(true_beam_interactingEnergy >= maxval) continue;

    // Fill histogram in reco energy
    for(Int_t l = 1; l <= nrecobins; l++){
      if(reco_beam_interactingEnergy > recoBins[l-1] && reco_beam_interactingEnergy <= recoBins[l]){
	mchisto->SetBinContent(l, mchisto->GetBinContent(l) + weight);
	break;
      }
    }
  }

  file->Close();

  return mchisto;

}

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, int channel, bool IsIncidentHisto){
  //********************************************************************

  TFile *file = new TFile(filename.c_str(), "READ");
  TTree *defaultTree  = (TTree*)file->Get(treename.c_str());

  Int_t reco_beam_type; // 13 -> track-like, 11 -> shower-like
  Int_t reco_beam_nTrackDaughters, reco_beam_nShowerDaughters;
  Double_t reco_beam_len, reco_beam_vtxX, reco_beam_vtxY, reco_beam_vtxZ, reco_beam_startX, reco_beam_startY, reco_beam_startZ, reco_beam_trackDirZ, reco_beam_interactingEnergy, reco_beam_Chi2_proton;
  bool reco_beam_true_byHits_matched; // Does the true particle contributing most to the reconstructed beam track coincide with the actual beam particle that generated the event
  Int_t reco_beam_true_byHits_origin, reco_beam_true_byHits_PDG;  // Origin and PDG of the reconstructed beam track
  Int_t true_beam_PDG;
  Double_t true_beam_interactingEnergy;
  std::string *true_beam_endProcess = 0;
  std::string *reco_beam_true_byHits_endProcess = 0;
  std::vector<double> *reco_beam_incidentEnergies = 0;

  defaultTree->SetBranchAddress("reco_beam_type",                   &reco_beam_type);
  defaultTree->SetBranchAddress("reco_beam_len",                    &reco_beam_len);
  defaultTree->SetBranchAddress("reco_beam_vtxX",                   &reco_beam_vtxX);
  defaultTree->SetBranchAddress("reco_beam_vtxY",                   &reco_beam_vtxY);
  defaultTree->SetBranchAddress("reco_beam_vtxZ",                   &reco_beam_vtxZ);
  defaultTree->SetBranchAddress("reco_beam_startX",                 &reco_beam_startX);
  defaultTree->SetBranchAddress("reco_beam_startY",                 &reco_beam_startY);
  defaultTree->SetBranchAddress("reco_beam_startZ",                 &reco_beam_startZ);
  defaultTree->SetBranchAddress("reco_beam_trackDirZ",              &reco_beam_trackDirZ);
  defaultTree->SetBranchAddress("reco_beam_nTrackDaughters",        &reco_beam_nTrackDaughters);
  defaultTree->SetBranchAddress("reco_beam_nShowerDaughters",       &reco_beam_nShowerDaughters);
  defaultTree->SetBranchAddress("reco_beam_interactingEnergy",      &reco_beam_interactingEnergy);
  defaultTree->SetBranchAddress("reco_beam_Chi2_proton",            &reco_beam_Chi2_proton);
  defaultTree->SetBranchAddress("reco_beam_incidentEnergies",       &reco_beam_incidentEnergies);

  defaultTree->SetBranchAddress("reco_beam_true_byHits_matched",    &reco_beam_true_byHits_matched);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_origin",     &reco_beam_true_byHits_origin);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_PDG",        &reco_beam_true_byHits_PDG);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_endProcess", &reco_beam_true_byHits_endProcess);

  defaultTree->SetBranchAddress("true_beam_interactingEnergy",      &true_beam_interactingEnergy);
  defaultTree->SetBranchAddress("true_beam_PDG",                    &true_beam_PDG);
  defaultTree->SetBranchAddress("true_beam_endProcess",             &true_beam_endProcess);

  const int nrecobins = recoBins.size();
  TH1D* datahisto = new TH1D(Form("Data_Channel%i_Histo",channel), Form("Data for channel %i",channel), nrecobins-1, 0, nrecobins-1);
  if(IsIncidentHisto)
    datahisto->SetNameTitle( Form("Data_ChannelIncident%i_Histo",channel), Form("Incident Data for channel %i", channel) );
  datahisto->SetDirectory(0);

  mf::LogInfo("FillDataHistogram_Pions") << "Filling data histogram " << datahisto->GetName() << " from file " << filename.c_str() << " for channel " << channel;

  for(Int_t k=0; k < defaultTree->GetEntries(); k++){
    defaultTree->GetEntry(k);

    // Temporary selection of events. To speed up the fit the selection should be done in advance.
    if(true_beam_PDG != 211) continue; // should be removed when proper data input is used

    // Track-like beam particles
    if(reco_beam_type != 13) continue;
    
    // Position cuts
    if(reco_beam_startZ < 28.0  || reco_beam_startZ > 32.0) continue;
    if(reco_beam_startY < 410.0 || reco_beam_startY > 445.0) continue;
    if(reco_beam_startX < -45.0 || reco_beam_startX > 0.0) continue;

    // Direction cuts
    if(reco_beam_trackDirZ < 0.9) continue;
    
    // Remove secondary protons from pion interacting outside the TPC
    if(reco_beam_Chi2_proton < 450.0) continue;

    // Minimum and maximum track length
    if(reco_beam_len < 5.0 || reco_beam_len > 270.0) continue;

    // Sometimes the reco energy at vertex is mis-reconstructed
    if(reco_beam_interactingEnergy < 0.0) continue;

    // For inelastic scattering at least one daughter is expected
    Int_t allDaughters = reco_beam_nTrackDaughters + reco_beam_nShowerDaughters;
    if(allDaughters <= 0) continue;

    if(IsIncidentHisto){
      for(UInt_t l = 0; l < reco_beam_incidentEnergies->size(); l++){
	for(Int_t m = 1; m <= nrecobins; m++){
	  if(reco_beam_incidentEnergies->at(l) > recoBins[m-1] && reco_beam_incidentEnergies->at(l) <= recoBins[m]){
	    datahisto->SetBinContent(m, datahisto->GetBinContent(m) + 1);
	    break;
	  }
	}
      }
      continue;
    }

    // Fill histogram in reco energy
    for(Int_t l = 1; l <= nrecobins; l++){
      if(reco_beam_interactingEnergy > recoBins[l-1] && reco_beam_interactingEnergy <= recoBins[l]){
	datahisto->SetBinContent(l, datahisto->GetBinContent(l) + 1);
	break;
      }
    }
    
  }

  file->Close();

  return datahisto;

}

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillMCTruthSignalHistogram_Pions(std::string filename, std::string treename, std::vector<double> truthBins, int channel, double weight){
  //********************************************************************

  TFile *file = new TFile(filename.c_str(), "READ");
  TTree *truthTree  = (TTree*)file->Get(treename.c_str());

  Int_t true_beam_PDG;
  Double_t true_beam_interactingEnergy, true_beam_endZ;
  std::string *true_beam_endProcess = 0;

  truthTree->SetBranchAddress("true_beam_interactingEnergy",      &true_beam_interactingEnergy);
  truthTree->SetBranchAddress("true_beam_PDG",                    &true_beam_PDG);
  truthTree->SetBranchAddress("true_beam_endZ",                   &true_beam_endZ);
  truthTree->SetBranchAddress("true_beam_endProcess",             &true_beam_endProcess);

  const int ntruthbins = truthBins.size();
  TH1D* mchisto = new TH1D(Form("MC_Channel%i_TruthSig_Histo",channel), Form("MC Truth Signal for channel %i",channel), ntruthbins-1, 0, ntruthbins-1);
  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCTruthSignalHistogram_Pions") << "Filling MC truth signal histogram " << mchisto->GetName() << " from file " << filename.c_str() << " for channel " << channel;

  for(Int_t k=0; k < truthTree->GetEntries(); k++){
    truthTree->GetEntry(k);

    // Pion beam
    if(true_beam_PDG != 211) continue;

    // If the pion does not make at the front face of the TPC, then don't count. Probably need a better estimate.
    if(true_beam_endZ < 0.0) continue;
    
    // True energy bin
    //if(true_beam_interactingEnergy < minval) continue;
    //if(true_beam_interactingEnergy >= maxval) continue;

    for(Int_t l = 1; l <= ntruthbins; l++){
      if(true_beam_interactingEnergy > truthBins[l-1] && true_beam_interactingEnergy <= truthBins[l]){
	mchisto->SetBinContent(l, mchisto->GetBinContent(l) + weight);
	break;
      }
    }
  }

  file->Close();

  return mchisto;
}

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillMCFlux_Pions(std::string filename, std::string treename, std::vector<double> Bins, int mode, double weight){
  //********************************************************************

  TFile *file = new TFile(filename.c_str(), "READ");
  TTree *truthTree  = (TTree*)file->Get(treename.c_str());

  Int_t reco_beam_type; // 13 -> track-like, 11 -> shower-like
  Double_t reco_beam_len, reco_beam_startX, reco_beam_startY, reco_beam_startZ, reco_beam_trackDirZ, reco_beam_interactingEnergy, reco_beam_Chi2_proton;
  Int_t true_beam_PDG, reco_beam_true_byHits_origin, reco_beam_true_byHits_PDG, reco_beam_nTrackDaughters, reco_beam_nShowerDaughters;
  Double_t true_beam_interactingEnergy, true_beam_endZ;
  std::vector<double> *true_beam_incidentEnergies = 0;  std::vector<double> *reco_beam_incidentEnergies = 0;

  truthTree->SetBranchAddress("reco_beam_type",                   &reco_beam_type);
  truthTree->SetBranchAddress("reco_beam_len",                    &reco_beam_len);
  truthTree->SetBranchAddress("reco_beam_startX",                 &reco_beam_startX);
  truthTree->SetBranchAddress("reco_beam_startY",                 &reco_beam_startY);
  truthTree->SetBranchAddress("reco_beam_startZ",                 &reco_beam_startZ);
  truthTree->SetBranchAddress("reco_beam_trackDirZ",              &reco_beam_trackDirZ);
  truthTree->SetBranchAddress("reco_beam_interactingEnergy",      &reco_beam_interactingEnergy);
  truthTree->SetBranchAddress("reco_beam_Chi2_proton",            &reco_beam_Chi2_proton);
  truthTree->SetBranchAddress("reco_beam_true_byHits_origin",     &reco_beam_true_byHits_origin);
  truthTree->SetBranchAddress("reco_beam_true_byHits_PDG",        &reco_beam_true_byHits_PDG);
  truthTree->SetBranchAddress("reco_beam_incidentEnergies",       &reco_beam_incidentEnergies);
  truthTree->SetBranchAddress("reco_beam_nTrackDaughters",        &reco_beam_nTrackDaughters);
  truthTree->SetBranchAddress("reco_beam_nShowerDaughters",       &reco_beam_nShowerDaughters);

  truthTree->SetBranchAddress("true_beam_interactingEnergy",      &true_beam_interactingEnergy);
  truthTree->SetBranchAddress("true_beam_PDG",                    &true_beam_PDG);
  truthTree->SetBranchAddress("true_beam_endZ",                   &true_beam_endZ);
  truthTree->SetBranchAddress("true_beam_incidentEnergies",       &true_beam_incidentEnergies);
  
  const int nbins = Bins.size();
  TH1D* mchisto = new TH1D(Form("MC_PionFlux%i_Histo",mode), Form("MC Pion Flux - %i",mode), nbins-1, 0, nbins-1);
  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCFlux_Pions") << "Filling MC pion flux histogram " << mchisto->GetName() << " from file " << filename.c_str() << " for mode " << mode;

  for(Int_t k=0; k < truthTree->GetEntries(); k++){
    truthTree->GetEntry(k);

    // Pion beam
    if(true_beam_PDG != 211) continue;

    // If the pion does not make at the front face of the TPC, then don't count. Probably need a better estimate.
    if(true_beam_endZ < 0.0) continue;

    if(true_beam_incidentEnergies->size() == 0 && reco_beam_incidentEnergies->size() == 0) continue;

    if(mode == 1){
      for(UInt_t l = 0; l < true_beam_incidentEnergies->size(); l++){
	for(Int_t m = 1; m <= nbins; m++){
	  if(true_beam_incidentEnergies->at(l) > Bins[m-1] && true_beam_incidentEnergies->at(l) <= Bins[m]){
	    mchisto->SetBinContent(m, mchisto->GetBinContent(m) + weight);
	    break;
	  }
	}
      }
    }
    else{
      if(reco_beam_true_byHits_origin != 4) continue;
      if(reco_beam_true_byHits_PDG != 211) continue;
      if(reco_beam_incidentEnergies->size() == 0) continue;

      // Track-like beam particles
      if(reco_beam_type != 13) continue;
      
      // Position cuts
      if(reco_beam_startZ < 28.0  || reco_beam_startZ > 32.0) continue;
      if(reco_beam_startY < 410.0 || reco_beam_startY > 445.0) continue;
      if(reco_beam_startX < -45.0 || reco_beam_startX > 0.0) continue;
      
      // Direction cuts
      if(reco_beam_trackDirZ < 0.9) continue;
      
      // Remove secondary protons from pion interacting outside the TPC
      if(reco_beam_Chi2_proton < 450.0) continue;
      
      // Minimum and maximum track length
      if(reco_beam_len < 5.0 || reco_beam_len > 270.0) continue;
      
      // Sometimes the reco energy at vertex is mis-reconstructed
      if(reco_beam_interactingEnergy < 0.0) continue;

      // For inelastic scattering at least one daughter is expected
      Int_t allDaughters = reco_beam_nTrackDaughters + reco_beam_nShowerDaughters;
      if(allDaughters <= 0) continue;

      if(mode == 2){
	for(UInt_t l = 0; l < true_beam_incidentEnergies->size(); l++){
	  for(Int_t m = 1; m <= nbins; m++){
	    if(true_beam_incidentEnergies->at(l) > Bins[m-1] && true_beam_incidentEnergies->at(l) <= Bins[m]){
	      mchisto->SetBinContent(m, mchisto->GetBinContent(m) + weight);
	      break;
	    }
	  }
	}
      }
      else if(mode == 3){
	for(UInt_t l = 0; l < reco_beam_incidentEnergies->size(); l++){
	  for(Int_t m = 1; m <= nbins; m++){
	    if(reco_beam_incidentEnergies->at(l) > Bins[m-1] && reco_beam_incidentEnergies->at(l) <= Bins[m]){
	      mchisto->SetBinContent(m, mchisto->GetBinContent(m) + weight);
	      break;
	    }
	  }
	}
      }
    }
  }
  
  file->Close();

  return mchisto;

}

//********************************************************************
int protoana::ProtoDUNESelectionUtils::GetNTriggers_Pions(std::string filename, std::string treename, bool IsMC){
  //********************************************************************

  TFile *file = new TFile(filename.c_str(), "READ");
  TTree *defaultTree  = (TTree*)file->Get(treename.c_str());

  int ntriggers = 0;

  Int_t true_beam_PDG, reco_beam_true_byHits_origin;
  
  defaultTree->SetBranchAddress("true_beam_PDG",                    &true_beam_PDG);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_origin",     &reco_beam_true_byHits_origin);

  for(Int_t k=0; k < defaultTree->GetEntries(); k++){
    defaultTree->GetEntry(k);

    if(true_beam_PDG != 211) continue; // should be removed when proper data input is used
    if(reco_beam_true_byHits_origin != 4) continue;
    ntriggers++;

  }

  return ntriggers;

}
