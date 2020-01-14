#include "ProtoDUNESelectionUtils.h"

#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "messagefacility/MessageLogger/MessageLogger.h"

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, int channel, int topo){
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

  defaultTree->SetBranchAddress("reco_beam_true_byHits_matched",    &reco_beam_true_byHits_matched);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_origin",     &reco_beam_true_byHits_origin);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_PDG",        &reco_beam_true_byHits_PDG);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_endProcess", &reco_beam_true_byHits_endProcess);

  defaultTree->SetBranchAddress("true_beam_interactingEnergy",      &true_beam_interactingEnergy);
  defaultTree->SetBranchAddress("true_beam_PDG",                    &true_beam_PDG);
  defaultTree->SetBranchAddress("true_beam_endProcess",             &true_beam_endProcess);

  //defaultTree->SetBranchAddress("",        &);

  const int nrecobins = recoBins.size();
  TH1D* mchisto = new TH1D(Form("MC_Channel%i_BkgTopo%i_Histo",channel,topo), Form("MC Background channel %i and topology %i", channel,topo), nrecobins-1, 0, nrecobins-1);
  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCBackgroundHistogram") << "Filling MC background histogram from file " << filename.c_str() << " for channel " << channel << " with topology " << topo;

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

    // Fill histogram in reco energy
    for(Int_t l = 1; l <= nrecobins; l++){
      if(reco_beam_interactingEnergy > recoBins[l-1] && reco_beam_interactingEnergy <= recoBins[l]){
	mchisto->SetBinContent(l, mchisto->GetBinContent(l) + 1);
	break;
      }
    }
  }

  file->Close();

  return mchisto;

}

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, int channel, int topo, double minval, double maxval){
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

  defaultTree->SetBranchAddress("reco_beam_true_byHits_matched",    &reco_beam_true_byHits_matched);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_origin",     &reco_beam_true_byHits_origin);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_PDG",        &reco_beam_true_byHits_PDG);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_endProcess", &reco_beam_true_byHits_endProcess);

  defaultTree->SetBranchAddress("true_beam_interactingEnergy",      &true_beam_interactingEnergy);
  defaultTree->SetBranchAddress("true_beam_PDG",                    &true_beam_PDG);
  defaultTree->SetBranchAddress("true_beam_endProcess",             &true_beam_endProcess);

  const int nrecobins = recoBins.size();
  TH1D* mchisto = new TH1D(Form("MC_Channel%i_SigTopo%i_%.2f-%.2f_Histo",channel,topo,minval,maxval), Form("MC Signal for channel %i and topology %i and true region %.1f-%.1f", channel,topo,minval,maxval), nrecobins-1, 0, nrecobins-1);
  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCSignalHistogram") << "Filling MC signal histogram from file " << filename.c_str() << " for channel " << channel << " with topology " << topo << " in the true region " << minval << "-" << maxval;

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

    // True energy bin
    if(true_beam_interactingEnergy < minval) continue;
    if(true_beam_interactingEnergy >= maxval) continue;

    // Fill histogram in reco energy
    for(Int_t l = 1; l <= nrecobins; l++){
      if(reco_beam_interactingEnergy > recoBins[l-1] && reco_beam_interactingEnergy <= recoBins[l]){
	mchisto->SetBinContent(l, mchisto->GetBinContent(l) + 1);
	break;
      }
    }
  }

  file->Close();

  return mchisto;

}

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, int channel){
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

  defaultTree->SetBranchAddress("reco_beam_true_byHits_matched",    &reco_beam_true_byHits_matched);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_origin",     &reco_beam_true_byHits_origin);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_PDG",        &reco_beam_true_byHits_PDG);
  defaultTree->SetBranchAddress("reco_beam_true_byHits_endProcess", &reco_beam_true_byHits_endProcess);

  defaultTree->SetBranchAddress("true_beam_interactingEnergy",      &true_beam_interactingEnergy);
  defaultTree->SetBranchAddress("true_beam_PDG",                    &true_beam_PDG);
  defaultTree->SetBranchAddress("true_beam_endProcess",             &true_beam_endProcess);

  const int nrecobins = recoBins.size();
  TH1D* datahisto = new TH1D(Form("Data_Channel%i_Histo",channel), Form("Data histogram for channel %i",channel), nrecobins-1, 0, nrecobins-1);
  datahisto->SetDirectory(0);

  mf::LogInfo("FillDataHistogram") << "Filling data histogram from file " << filename.c_str() << " for channel " << channel;

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
