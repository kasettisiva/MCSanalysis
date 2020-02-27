#include "ProtoDUNESelectionUtils.h"

#include <algorithm>

// ROOT
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "messagefacility/MessageLogger/MessageLogger.h"

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, std::string channel, std::string topo, int toponum, double weight){
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

  Int_t true_chexSignal, true_absSignal, true_backGround, true_nPi0Signal;

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

  defaultTree->SetBranchAddress("true_chexSignal",                  &true_chexSignal);
  defaultTree->SetBranchAddress("true_nPi0Signal",                  &true_nPi0Signal);
  defaultTree->SetBranchAddress("true_absSignal",                   &true_absSignal);
  defaultTree->SetBranchAddress("true_backGround",                  &true_backGround);

  //New: backgrounds within the tree
  bool primaryMuon, isCosmic, isExtraBeam, upstreamInt, isDecay;
  defaultTree->SetBranchAddress("primaryMuon",                  &primaryMuon);
  defaultTree->SetBranchAddress("isCosmic",                     &isCosmic);
  defaultTree->SetBranchAddress("isExtraBeam",                  &isExtraBeam);
  defaultTree->SetBranchAddress("upstreamInt",                  &upstreamInt);
  defaultTree->SetBranchAddress("isDecay",                      &isDecay);


  channel.erase(std::remove(channel.begin(), channel.end(), '.'), channel.end());
  channel.erase(std::remove(channel.begin(), channel.end(), ' '), channel.end());
  topo.erase(std::remove(topo.begin(), topo.end(), '.'), topo.end());
  topo.erase(std::remove(topo.begin(), topo.end(), ' '), topo.end());

  const int nrecobins = recoBins.size();
  TH1D* mchisto = new TH1D(Form("MC_Channel%s_%s_Histo",channel.c_str(),topo.c_str()), Form("MC Background for channel %s and topology %s", channel.c_str(),topo.c_str()), nrecobins-1, 0, nrecobins-1);
  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCBackgroundHistogram_Pions") << "Filling MC background histogram " << mchisto->GetName() << " from file " << filename.c_str() << " for channel " << channel.c_str() << " with topology " << topo.c_str();

  for(Int_t k=0; k < defaultTree->GetEntries(); k++){
    defaultTree->GetEntry(k);

    // Different background topologies
    Int_t topology = -1;
    std::string reco_beam_true_byHits_endProcess_str = *reco_beam_true_byHits_endProcess;


    // New From Jake
    if( upstreamInt )
      topology = 4;
    else if( primaryMuon )
      topology = 5;
    else if( isDecay )
      topology = 6;
    else if( isExtraBeam )
      //topology = 7;
      topology = 7;
    else if( isCosmic )
      //topology = 7;
      topology = 7;
    else{
      if( reco_beam_true_byHits_matched ){
        if(  true_backGround  )
          topology = 3;
        //Else: this is signal
        //should be caught by topology != toponum check
      }
      else
        topology = 7;
    }
       

    /////////////////////


    // Select only the correct topology
    if(topology != toponum) continue;

    // Sometimes the reco energy at vertex is mis-reconstructed
    if(reco_beam_interactingEnergy < 0.0) continue;

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
TH1* protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, std::string channel, std::string topo, int toponum, double minval, double maxval, double weight){
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
  Int_t true_chexSignal, true_absSignal, true_backGround, true_nPi0Signal;

  int event, run;

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

  defaultTree->SetBranchAddress("true_chexSignal",                  &true_chexSignal);
  defaultTree->SetBranchAddress("true_nPi0Signal",                  &true_nPi0Signal);
  defaultTree->SetBranchAddress("true_absSignal",                   &true_absSignal);
  defaultTree->SetBranchAddress("true_backGround",                  &true_backGround);

  defaultTree->SetBranchAddress("event", &event);
  defaultTree->SetBranchAddress("run", &run);

  bool primaryMuon, isCosmic, isExtraBeam, upstreamInt, isDecay;
  defaultTree->SetBranchAddress("primaryMuon",                  &primaryMuon);
  defaultTree->SetBranchAddress("isCosmic",                     &isCosmic);
  defaultTree->SetBranchAddress("isExtraBeam",                  &isExtraBeam);
  defaultTree->SetBranchAddress("upstreamInt",                  &upstreamInt);
  defaultTree->SetBranchAddress("isDecay",                      &isDecay);


  channel.erase(std::remove(channel.begin(), channel.end(), '.'), channel.end());
  channel.erase(std::remove(channel.begin(), channel.end(), ' '), channel.end());
  topo.erase(std::remove(topo.begin(), topo.end(), '.'), topo.end());
  topo.erase(std::remove(topo.begin(), topo.end(), ' '), topo.end());

  const int nrecobins = recoBins.size();
  TH1D* mchisto = new TH1D(Form("MC_Channel%s_%s_%.1f-%.1f_Histo",channel.c_str(),topo.c_str(),minval,maxval), Form("MC Signal for channel %s and topology %s and true region %.1f-%.1f", channel.c_str(),topo.c_str(),minval,maxval), nrecobins-1, 0, nrecobins-1);
  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCSignalHistogram_Pions") << "Filling MC signal histogram " << mchisto->GetName() << " from file " << filename.c_str() << " for channel " << channel.c_str() << " with topology " << topo.c_str() << " in the true region " << minval << "-" << maxval;

  for(Int_t k=0; k < defaultTree->GetEntries(); k++){
    defaultTree->GetEntry(k);

    /*
     *  Jake: Should we define the true signal for Charge Exchange
     *        as any number of pi0? If so, then we'll have to slightly
     *        change the definition here, and we'll have to make Cex
     *        (topology 2) include true_nPi0Signal
     *
     *
     *
     *  NEED TO EDIT IN THE NEW DEFINITIONS FOR SIGNAL/BACKGROUND
     */


    Int_t topology = -1;

    if( !( upstreamInt || primaryMuon || isDecay || isExtraBeam || isCosmic ) && reco_beam_true_byHits_matched && !true_backGround ){
      if( true_absSignal )
        topology = 1;
      else if( true_chexSignal || true_nPi0Signal )
        topology = 2;
      else{
        std::cout << "WARNING!!! OTHER CASE CAUGHT. NEED TO INVESTIGATE" << std::endl;
        std::cout << "Event: " << event << " Run: " << run << std::endl;
      }
    }

    /*
    if(true_absSignal == 1) topology = 1;
    //else if(true_chexSignal == 1) topology = 2;
    //New: including nPi0
    else if( (true_chexSignal == 1) || (true_nPi0Signal == 1) ) topology = 2;
    */

    // Remove events with the wrong topology
    if(topology != toponum) continue;

    // Classified as ABS and CEX signal => Bkg
    // This never happens
    //if(true_chexSignal == 1 && true_absSignal == 1) continue;

    // Sometimes the reco energy at vertex is mis-reconstructed
    if(reco_beam_interactingEnergy < 0.0) continue;

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
TH1* protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, std::string channel, bool IsIncidentHisto){
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

  channel.erase(std::remove(channel.begin(), channel.end(), '.'), channel.end());
  channel.erase(std::remove(channel.begin(), channel.end(), ' '), channel.end());

  const int nrecobins = recoBins.size();
  TH1D* datahisto = new TH1D(Form("Data_Channel%s_Histo",channel.c_str()), Form("Data for channel %s",channel.c_str()), nrecobins-1, 0, nrecobins-1);
  if(IsIncidentHisto)
    datahisto->SetNameTitle( Form("Data_ChannelIncident%s_Histo",channel.c_str()), Form("Incident Data for channel %s", channel.c_str()) );
  datahisto->SetDirectory(0);

  mf::LogInfo("FillDataHistogram_Pions") << "Filling data histogram " << datahisto->GetName() << " from file " << filename.c_str() << " for channel " << channel.c_str();

  for(Int_t k=0; k < defaultTree->GetEntries(); k++){
    defaultTree->GetEntry(k);

    // Sometimes the reco energy at vertex is mis-reconstructed
    if(reco_beam_interactingEnergy < 0.0) continue;

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
TH1* protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, std::string channel, std::string topo, int toponum, double weight){
  //********************************************************************

  TFile *file = new TFile(filename.c_str(), "READ");
  TTree *defaultTree  = (TTree*)file->Get(treename.c_str());

  Int_t reco_beam_type; // 13 -> track-like, 11 -> shower-like
  Int_t reco_beam_nTrackDaughters, reco_beam_nShowerDaughters;
  Double_t reco_beam_len, reco_beam_vtxX, reco_beam_vtxY, reco_beam_vtxZ, reco_beam_startX, reco_beam_startY, reco_beam_startZ, reco_beam_trackDirZ, reco_beam_interactingEnergy, reco_beam_Chi2_proton;
  bool reco_beam_true_byHits_matched; // Does the true particle contributing most to the reconstructed beam track coincide with the actual beam particle that generated the event
  Int_t reco_beam_true_byHits_origin, reco_beam_true_byHits_PDG;  // Origin and PDG of the reconstructed beam track
  Int_t true_beam_PDG, true_beam_ID;
  Double_t true_beam_interactingEnergy;
  std::string *true_beam_endProcess = 0;
  std::string *reco_beam_true_byHits_endProcess = 0;
  std::vector<double> *reco_beam_incidentEnergies = 0;

  Int_t true_chexSignal, true_absSignal, true_backGround;

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
  defaultTree->SetBranchAddress("true_beam_ID",                    &true_beam_ID);
  defaultTree->SetBranchAddress("true_beam_endProcess",             &true_beam_endProcess);

  defaultTree->SetBranchAddress("true_chexSignal",                  &true_chexSignal);
  defaultTree->SetBranchAddress("true_absSignal",                   &true_absSignal);
  defaultTree->SetBranchAddress("true_backGround",                  &true_backGround);

  std::vector< int > * true_beam_slices = 0x0;
  std::vector< int > * reco_beam_hit_true_ID = 0x0;
  std::vector< int > * reco_beam_hit_true_origin = 0x0;
  std::vector< int > * reco_beam_hit_true_slice = 0x0;
  std::vector< double > * true_beam_slices_deltaE = 0x0;

  defaultTree->SetBranchAddress("true_beam_slices", &true_beam_slices);
  defaultTree->SetBranchAddress("true_beam_slices_deltaE", &true_beam_slices_deltaE);
  defaultTree->SetBranchAddress("reco_beam_hit_true_ID", &reco_beam_hit_true_ID);
  defaultTree->SetBranchAddress("reco_beam_hit_true_slice", &reco_beam_hit_true_slice);
  defaultTree->SetBranchAddress("reco_beam_hit_true_origin", &reco_beam_hit_true_origin);

  std::vector< std::string > * true_beam_processes = 0x0;
  std::vector< int > * true_beam_process_matched = 0x0;
  std::vector< int > * true_beam_process_dSlice = 0x0;
  std::vector< int > * true_beam_process_slice = 0x0;
  defaultTree->SetBranchAddress("true_beam_processes", &true_beam_processes);
  defaultTree->SetBranchAddress("true_beam_process_matched", &true_beam_process_matched);
  defaultTree->SetBranchAddress("true_beam_process_dSlice", &true_beam_process_dSlice);
  defaultTree->SetBranchAddress("true_beam_process_slice", &true_beam_process_slice);

  //std::replace(channel.begin(), channel.end(), ' ', '-');
  channel.erase(std::remove(channel.begin(), channel.end(), '.'), channel.end());
  channel.erase(std::remove(channel.begin(), channel.end(), ' '), channel.end());
  //std::replace(topo.begin(), topo.end(), ' ', '-');
  topo.erase(std::remove(topo.begin(), topo.end(), '.'), topo.end());
  topo.erase(std::remove(topo.begin(), topo.end(), ' '), topo.end());

  size_t nrecobins = recoBins.size();
  TH1D* mchisto = new TH1D(Form("MC_ChannelIncident%s_%s_Histo",channel.c_str(),topo.c_str()), Form("Incident MC for channel %s and topology %s", channel.c_str(),topo.c_str()), nrecobins-1, 0, nrecobins-1);
  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCBackgroundHistogram_Pions") << "Filling MC background histogram " << mchisto->GetName() << " from file " << filename.c_str() << " for channel " << channel.c_str() << " with topology " << topo.c_str();

  for(Int_t k=0; k < defaultTree->GetEntries(); k++){
    defaultTree->GetEntry(k);

    // Different background topologies
    Int_t topology = -1;
 
 /*
    if(true_chexSignal == 1 || true_absSignal == 1 || reco_beam_true_byHits_PDG == 211) topology = 3; // pions
    else if(abs(reco_beam_true_byHits_PDG) == 2212) topology = 4; // protons
    else if(abs(reco_beam_true_byHits_PDG) == 13) topology = 5; // muons
    else topology = 6; // other

    // Select only the correct topology
    if(topology != toponum) continue;
*/
    // Sometimes the reco energy at vertex is mis-reconstructed
    if(reco_beam_interactingEnergy < 0.0) continue;

    int max_slice_found = -999;
    //Go through the reco beam hits and find the max slice found
    for( size_t l = 0; l < reco_beam_hit_true_ID->size(); ++l ){
      int true_id = (*reco_beam_hit_true_ID)[l]; 
      int true_slice = (*reco_beam_hit_true_slice)[l]; 

      //Found the slice
      if( true_id == true_beam_ID && true_slice != -999 ){
        if( true_slice > max_slice_found ) max_slice_found = true_slice;
      }
    }

    //Determine if any process was matched 
    size_t max_matched = 999;
    bool any_matched = false;
    for( size_t l = 0; l < true_beam_process_matched->size(); ++l ){
      if( (*true_beam_process_matched)[l] ){
        max_matched = l; 
        any_matched = true;
      }
    }


    bool found_last_proc = false;
    if( any_matched ){
      //Check if the max_matched is the last
      if( max_matched == (true_beam_process_matched->size()-1) ){
        //This means we should treat all of the slices up to the last process
        found_last_proc = true;
      }
    }



    bool passed_last_found = false;
    //Go through the incident energies in the beamline
    for(size_t l = 0; l < reco_beam_incidentEnergies->size(); ++l){

      //Check the ID, origin, true_slice
      int true_id = (*reco_beam_hit_true_ID)[l]; 
      int true_origin = (*reco_beam_hit_true_origin)[l]; 
      int true_slice = (*reco_beam_hit_true_slice)[l]; 

      if( max_slice_found > -999 && !passed_last_found && found_last_proc ){
        if( true_slice >= max_slice_found ) passed_last_found = true;
      }

      if( true_origin == 2 ){
        topology = 5;
      }
      else{
        if( true_id != true_beam_ID ){
          if( true_id == -999 ) topology = 8; 
          else{
            if( passed_last_found ) topology = 7;
            else topology = 8;
          }
        }
        else{
          //True id matched to hit, but no slice matched
          if( true_slice == -999 ){
            topology = 6;
          }
          else{
            if( true_beam_PDG == 211 ) topology = 3;
            else if( true_beam_PDG == -13 ) topology = 4;
          }
        }
      }


      if(topology != toponum) continue;

      for(size_t m = 1; m <= nrecobins; m++){
	if(reco_beam_incidentEnergies->at(l) > recoBins[m-1] && reco_beam_incidentEnergies->at(l) <= recoBins[m]){
	  mchisto->SetBinContent(m, mchisto->GetBinContent(m) + weight);
	  break;
	}
      }
    }
  }

  file->Close();

  return mchisto;

}

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillMCTruthSignalHistogram_Pions(std::string filename, std::string treename, std::vector<double> truthBins, std::string channel, double weight){
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

  channel.erase(std::remove(channel.begin(), channel.end(), '.'), channel.end());
  channel.erase(std::remove(channel.begin(), channel.end(), ' '), channel.end());

  const int ntruthbins = truthBins.size();
  TH1D* mchisto = new TH1D(Form("MC_Channel%s_TruthSig_Histo",channel.c_str()), Form("MC Truth Signal for channel %s",channel.c_str()), ntruthbins-1, 0, ntruthbins-1);
  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCTruthSignalHistogram_Pions") << "Filling MC truth signal histogram " << mchisto->GetName() << " from file " << filename.c_str() << " for channel " << channel.c_str();

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
      // Jake: This is deprecated in later versions
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
