#include "ProtoDUNESelectionUtils.h"

#include <algorithm>

// ROOT
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "messagefacility/MessageLogger/MessageLogger.h"

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(
    std::string filename, std::string treename, std::vector<double> recoBins,
    std::string channel, std::string topo, int toponum, double endZ_cut,
    double minval, double maxval, bool doNegativeReco, int doSyst, std::string systName, double weight) {
//********************************************************************

  TFile *file = new TFile(filename.c_str(), "READ");
  TTree *defaultTree  = (TTree*)file->Get(treename.c_str());

  int reco_beam_type; // 13 -> track-like, 11 -> shower-like
  int reco_beam_nTrackDaughters, reco_beam_nShowerDaughters;
  double reco_beam_len, reco_beam_vtxX, reco_beam_vtxY, reco_beam_vtxZ,
      reco_beam_startX, reco_beam_startY, reco_beam_startZ, reco_beam_trackDirZ,
      reco_beam_interactingEnergy, reco_beam_Chi2_proton;

  // Does the true particle contributing most to the reconstructed beam track 
  // coincide with the actual beam particle that generated the event
  bool reco_beam_true_byHits_matched; 

  // Origin and PDG of the reconstructed beam track
  int reco_beam_true_byHits_origin, reco_beam_true_byHits_PDG;
  int true_beam_PDG;
  double true_beam_endZ;
  double true_beam_interactingEnergy;
  std::string *true_beam_endProcess = 0;
  std::string *reco_beam_true_byHits_endProcess = 0;
  std::vector<double> *reco_beam_incidentEnergies = 0;

  int true_chexSignal, true_absSignal, true_backGround, true_nPi0Signal;

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
  defaultTree->SetBranchAddress("true_beam_endZ",                    &true_beam_endZ);
  defaultTree->SetBranchAddress("true_beam_endProcess",             &true_beam_endProcess);

  defaultTree->SetBranchAddress("true_chexSignal",                  &true_chexSignal);
  defaultTree->SetBranchAddress("true_nPi0Signal",                  &true_nPi0Signal);
  defaultTree->SetBranchAddress("true_absSignal",                   &true_absSignal);
  defaultTree->SetBranchAddress("true_backGround",                  &true_backGround);

  //Testing
  std::vector<double> * reco_beam_calibrated_dEdX = 0x0;
  defaultTree->SetBranchAddress("reco_beam_calibrated_dEdX", &reco_beam_calibrated_dEdX);
  std::vector<double> * reco_beam_TrkPitch = 0x0;
  defaultTree->SetBranchAddress("reco_beam_TrkPitch", &reco_beam_TrkPitch);

  double reco_beam_endX, reco_beam_endY, reco_beam_endZ;
  defaultTree->SetBranchAddress("reco_beam_endX",                 &reco_beam_endX);
  defaultTree->SetBranchAddress("reco_beam_endY",                 &reco_beam_endY);
  defaultTree->SetBranchAddress("reco_beam_endZ",                 &reco_beam_endZ);

  int run, event;
  defaultTree->SetBranchAddress("event", &event);
  defaultTree->SetBranchAddress("run", &run);
  /////////////////////////////////

  std::vector< int > * reco_beam_hit_true_origin = 0x0;
  std::vector< int > * reco_beam_hit_true_ID = 0x0;
  std::vector< int > * true_beam_daughter_ID = 0x0;
  std::vector< int > * true_beam_grand_daughter_ID = 0x0;
  int true_beam_ID;
  int true_daughter_nPiPlus, true_daughter_nPiMinus, true_daughter_nPi0;
  defaultTree->SetBranchAddress( "reco_beam_hit_true_origin", &reco_beam_hit_true_origin );
  defaultTree->SetBranchAddress( "reco_beam_hit_true_ID", &reco_beam_hit_true_ID );
    defaultTree->SetBranchAddress( "true_beam_ID", &true_beam_ID );
  defaultTree->SetBranchAddress( "true_daughter_nPiPlus", &true_daughter_nPiPlus );
  defaultTree->SetBranchAddress( "true_daughter_nPiMinus", &true_daughter_nPiMinus );
  defaultTree->SetBranchAddress( "true_daughter_nPi0", &true_daughter_nPi0 );
  defaultTree->SetBranchAddress( "true_beam_daughter_ID", &true_beam_daughter_ID );
  defaultTree->SetBranchAddress( "true_beam_grand_daughter_ID", &true_beam_grand_daughter_ID );

  //For vertex type 
  std::vector< std::string > * true_beam_processes = 0x0;
  std::vector< std::vector < double > > * reco_beam_vertex_dRs = 0x0;
  std::vector< int > * reco_beam_vertex_hits_slices = 0x0;
  defaultTree->SetBranchAddress( "true_beam_processes", &true_beam_processes );
  defaultTree->SetBranchAddress( "reco_beam_vertex_dRs", &reco_beam_vertex_dRs );
  defaultTree->SetBranchAddress( "reco_beam_vertex_hits_slices", &reco_beam_vertex_hits_slices );

  //For systs
  std::vector<double> * g4rw_primary_plus_sigma_weight = 0x0;
  std::vector<double> * g4rw_primary_minus_sigma_weight = 0x0;
  std::vector<std::string> * g4rw_primary_var = 0x0;
  if (doSyst != 0) {
    defaultTree->SetBranchAddress("g4rw_primary_plus_sigma_weight", &g4rw_primary_plus_sigma_weight);
    defaultTree->SetBranchAddress("g4rw_primary_minus_sigma_weight", &g4rw_primary_minus_sigma_weight);
    defaultTree->SetBranchAddress("g4rw_primary_var", &g4rw_primary_var);
  }

  channel.erase(std::remove(channel.begin(), channel.end(), '.'), channel.end());
  channel.erase(std::remove(channel.begin(), channel.end(), ' '), channel.end());
  topo.erase(std::remove(topo.begin(), topo.end(), '.'), topo.end());
  topo.erase(std::remove(topo.begin(), topo.end(), ' '), topo.end());

  const int nrecobins = recoBins.size();
  std::string hist_name = "MC_Channel" + channel + "_" + topo + "_Histo";
  std::string hist_title = "MC Background for channel " + channel + 
                           " and topology " + topo;
  
  if (doSyst == 1) {
    hist_name += "_high_" + systName;
    hist_title += " +1 sigma";
  }
  else if (doSyst == -1) {
    hist_name += "_low_" + systName;
    hist_title += " -1 sigma";
  }

  TH1D* mchisto = new TH1D(hist_name.c_str(), hist_title.c_str(), 
                           (doNegativeReco ? nrecobins : (nrecobins-1)),
                           (doNegativeReco ? -1 : 0), nrecobins-1);

  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCBackgroundHistogram_Pions") << 
      "Filling MC background histogram " << mchisto->GetName() << 
      " from file " << filename.c_str() << " for channel " << 
      channel.c_str() << " with topology " << topo.c_str();

  bool done_check = false;
  for(Int_t k=0; k < defaultTree->GetEntries(); k++){

    defaultTree->GetEntry(k);

    if (!done_check && doSyst != 0) {
      if (g4rw_primary_var->size() > 0) {
        std::cout << "Checking" << std::endl;
        auto syst_check = std::find(g4rw_primary_var->begin(), g4rw_primary_var->end(), systName);
        if (syst_check == g4rw_primary_var->end()) {
          std::cout << "Error! Could not find syst named " << systName << std::endl;
          std::exception e;
          throw e;
        }
        done_check = true;
      }
    }

    double syst_weight = 1.;
    std::map<std::string, double> weights;
    if (doSyst == 1) {
      for (size_t i = 0; i < g4rw_primary_var->size(); ++i) {
        weights[(*g4rw_primary_var)[i]] = (*g4rw_primary_plus_sigma_weight)[i];
      }
    }
    else if (doSyst == -1) {
      for (size_t i = 0; i < g4rw_primary_var->size(); ++i) {
        weights[(*g4rw_primary_var)[i]] = (*g4rw_primary_minus_sigma_weight)[i];
      }
    }

    // Different background topologies
    Int_t topology = -1;

    //Bin edges cases. If there are signal events that fall out the truth bins,
    //but are present in the reco bins then count them as backgrounds
    if (true_backGround == 0 && (true_chexSignal == 1 || true_absSignal == 1 ||
                                 true_nPi0Signal == 1)){
      if (true_beam_interactingEnergy < minval &&
          reco_beam_interactingEnergy > recoBins[0]){
        true_backGround = 1;
      }

      if (true_beam_interactingEnergy > maxval &&
          reco_beam_interactingEnergy < recoBins[nrecobins-1]){
        true_backGround = 1;
      }
    }

    /////////////////////
    
    if( reco_beam_true_byHits_origin == 2 ){
      topology = 7;  
    }
    else if( true_beam_PDG == 211 ){
      
      //make sure there is any hits at all
      if( !reco_beam_hit_true_ID->size() ){ 
        //std::cout << "no true_ID for hits " << reco_beam_interactingEnergy <<
        //             " " << reco_beam_calibrated_dEdX->size() << " " << 
        //             reco_beam_TrkPitch->size() << std::endl;
        //std::cout << "Checking: " << event << " " << run << std::endl;
        continue;
      }

      if (true_beam_endZ < 0.) {
        topology = 4;
      }
      else if (true_beam_endZ > endZ_cut) {
        topology = 7;
      }
      else if ( *true_beam_endProcess == "pi+Inelastic" ){
        if ( true_daughter_nPiPlus == 0 && true_daughter_nPiMinus == 0 ) {
          topology = 1;
        }
        else {
          topology = 3;
        }
      }
      else {
        topology = 7;
      }

      if (doSyst == 1 || doSyst == -1) { //Do +1 sigma
        syst_weight = weights[systName];
      }
    }
    else if ( true_beam_PDG == -13 ){
      //First check that if the last hit is from a cosmic
      if ( reco_beam_hit_true_origin->back() == 2 )
        topology = 7; //Other for now
      else if( reco_beam_hit_true_ID->back() == true_beam_ID )
        topology = 5;
      else
        topology = 7;
    }
    else{ //Shouldn't be any but w/e
      topology = 7;
    }

    // Select only the correct topology
    if (topology != toponum) continue;

    if (doNegativeReco) {
      if (reco_beam_interactingEnergy < 0.) {
        mchisto->AddBinContent(1, weight*syst_weight);
      }
      else{
        for (int l = 1; l < nrecobins; ++l) {
          if (reco_beam_interactingEnergy > recoBins[l-1] &&
              reco_beam_interactingEnergy <= recoBins[l]) {
            //Fill +1 because the first bin is negative
            mchisto->AddBinContent(l+1, weight*syst_weight);
            break;
          }
        }       
      }
    }
    else {
      // Sometimes the reco energy at vertex is mis-reconstructed
      if (reco_beam_interactingEnergy < 0.0) continue;

      for (int l = 1; l < nrecobins; ++l) {
        if (reco_beam_interactingEnergy > recoBins[l-1] &&
            reco_beam_interactingEnergy <= recoBins[l]) {
          mchisto->AddBinContent(l/*+1*/, weight*syst_weight);
          break;
        }
      }
    }
  }

  file->Close();

  return mchisto;

}

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(
    std::string filename, std::string treename, std::vector<double> recoBins,
    std::string channel, std::string topo, int toponum, double minval,
    double maxval, double endZ_cut, bool doNegativeReco, int doSyst,
    std::string systName, double weight) {
//********************************************************************

  TFile *file = new TFile(filename.c_str(), "READ");
  TTree *defaultTree  = (TTree*)file->Get(treename.c_str());

  int reco_beam_type; // 13 -> track-like, 11 -> shower-like
  int reco_beam_nTrackDaughters, reco_beam_nShowerDaughters;
  double reco_beam_len, reco_beam_vtxX, reco_beam_vtxY, reco_beam_vtxZ, reco_beam_startX, reco_beam_startY, reco_beam_startZ, reco_beam_trackDirZ, reco_beam_interactingEnergy, reco_beam_Chi2_proton;
  bool reco_beam_true_byHits_matched; // Does the true particle contributing most to the reconstructed beam track coincide with the actual beam particle that generated the event
  int reco_beam_true_byHits_origin, reco_beam_true_byHits_PDG;  // Origin and PDG of the reconstructed beam track
  int true_beam_PDG;
  double true_beam_endZ;
  double true_beam_interactingEnergy;
  std::string *true_beam_endProcess = 0;
  std::string *reco_beam_true_byHits_endProcess = 0;
  std::vector<double> *reco_beam_incidentEnergies = 0;
  int true_chexSignal, true_absSignal, true_backGround, true_nPi0Signal;

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
  defaultTree->SetBranchAddress("true_beam_endZ",                   &true_beam_endZ);
  defaultTree->SetBranchAddress("true_beam_endProcess",             &true_beam_endProcess);

  defaultTree->SetBranchAddress("true_chexSignal",                  &true_chexSignal);
  defaultTree->SetBranchAddress("true_nPi0Signal",                  &true_nPi0Signal);
  defaultTree->SetBranchAddress("true_absSignal",                   &true_absSignal);
  defaultTree->SetBranchAddress("true_backGround",                  &true_backGround);

  defaultTree->SetBranchAddress("event",                            &event);
  defaultTree->SetBranchAddress("run",                              &run);

  //Testing
  std::vector<double> * reco_beam_calibrated_dEdX = 0x0;
  defaultTree->SetBranchAddress("reco_beam_calibrated_dEdX", &reco_beam_calibrated_dEdX);
  std::vector<double> * reco_beam_TrkPitch = 0x0;
  defaultTree->SetBranchAddress("reco_beam_TrkPitch", &reco_beam_TrkPitch);

  /////////////////////////////////


  int true_daughter_nPiPlus, true_daughter_nPiMinus, true_daughter_nPi0, true_beam_ID;
  std::vector< std::string > * true_beam_processes = 0x0;
  std::vector< int > * reco_beam_hit_true_ID = 0x0;
  double new_true_beam_interactingEnergy;
  defaultTree->SetBranchAddress( "true_daughter_nPiPlus", &true_daughter_nPiPlus );
  defaultTree->SetBranchAddress( "true_daughter_nPiMinus", &true_daughter_nPiMinus );
  defaultTree->SetBranchAddress( "true_daughter_nPi0", &true_daughter_nPi0 );
  defaultTree->SetBranchAddress( "true_beam_ID", &true_beam_ID );
  defaultTree->SetBranchAddress( "true_beam_processes", &true_beam_processes );
  defaultTree->SetBranchAddress( "reco_beam_hit_true_ID", &reco_beam_hit_true_ID );
  defaultTree->SetBranchAddress("new_true_beam_interactingEnergy",      &new_true_beam_interactingEnergy);


  std::vector< int > * reco_beam_vertex_hits_slices = 0x0; 
  std::vector< std::vector< double > > * reco_beam_vertex_dRs = 0x0;
  defaultTree->SetBranchAddress( "reco_beam_vertex_hits_slices", &reco_beam_vertex_hits_slices );
  defaultTree->SetBranchAddress( "reco_beam_vertex_dRs", &reco_beam_vertex_dRs );

  //For systs
  std::vector<double> * g4rw_primary_plus_sigma_weight = 0x0;
  std::vector<double> * g4rw_primary_minus_sigma_weight = 0x0;
  std::vector<std::string> * g4rw_primary_var = 0x0;
  if (doSyst != 0) {
    defaultTree->SetBranchAddress("g4rw_primary_plus_sigma_weight", &g4rw_primary_plus_sigma_weight);
    defaultTree->SetBranchAddress("g4rw_primary_minus_sigma_weight", &g4rw_primary_minus_sigma_weight);
    defaultTree->SetBranchAddress("g4rw_primary_var", &g4rw_primary_var);
  }

  channel.erase(std::remove(channel.begin(), channel.end(), '.'), channel.end());
  channel.erase(std::remove(channel.begin(), channel.end(), ' '), channel.end());
  topo.erase(std::remove(topo.begin(), topo.end(), '.'), topo.end());
  topo.erase(std::remove(topo.begin(), topo.end(), ' '), topo.end());

  const int nrecobins = recoBins.size();
  TString hist_name = Form("MC_Channel%s_%s_%.1f-%.1f_Histo",
                           channel.c_str(), topo.c_str(), minval, maxval);
  TString hist_title = Form("MC Signal for channel %s and topology %s and true region %.1f-%.1f",
                            channel.c_str(), topo.c_str(), minval, maxval);

  if (doSyst == 1) {
    hist_name += "_high_" + systName;
    hist_title += " +1 sigma";
  }
  else if (doSyst == -1) {
    hist_name += "_low_" + systName;
    hist_title += " -1 sigma";
  }

  TH1D* mchisto = new TH1D(hist_name, hist_title,
                           (doNegativeReco ? nrecobins : (nrecobins-1)),
                           (doNegativeReco ? -1 : 0), nrecobins-1);
  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCSignalHistogram_Pions") << "Filling MC signal histogram " <<
              mchisto->GetName() << " from file " << filename.c_str() <<
              " for channel " << channel.c_str() << " with topology " <<
              topo.c_str() << " in the true region " << minval << "-" << maxval;

  bool done_check = false;
  for (int k=0; k < defaultTree->GetEntries(); k++) {

    defaultTree->GetEntry(k);

    if (!done_check && doSyst != 0) {
      if (g4rw_primary_var->size() > 0) {
        auto syst_check = std::find(g4rw_primary_var->begin(), g4rw_primary_var->end(), systName);
        if (syst_check == g4rw_primary_var->end()) {
          std::cout << "Error! Could not find syst named " << systName << std::endl;
          std::exception e;
          throw e;
        }
        done_check = true;
      }
    }

    double syst_weight = 1.;
    std::map<std::string, double> weights;
    if (doSyst == 1) {
      for (size_t i = 0; i < g4rw_primary_var->size(); ++i) {
        weights[(*g4rw_primary_var)[i]] = (*g4rw_primary_plus_sigma_weight)[i];
      }
    }
    else if (doSyst == -1) {
      for (size_t i = 0; i < g4rw_primary_var->size(); ++i) {
        weights[(*g4rw_primary_var)[i]] = (*g4rw_primary_minus_sigma_weight)[i];
      }
    }

    int topology = -1;

    if ( reco_beam_true_byHits_origin == 2 )
      continue; 

    //First determine if this is a pion ending in abs/cex/nPi0
    if ( true_beam_PDG == 211 && *true_beam_endProcess == "pi+Inelastic" &&
        true_daughter_nPiPlus == 0 && true_daughter_nPiMinus == 0 ){
      
      //make sure there is any hits at all
      if ( !reco_beam_hit_true_ID->size() ) {
        //std::cout << "no true_ID for hits " << reco_beam_interactingEnergy <<
        //             " " << reco_beam_calibrated_dEdX->size() << " " << 
        //             reco_beam_TrkPitch->size() << std::endl;
        continue;
      }

      if (true_beam_endZ < 0. || true_beam_endZ > endZ_cut) continue;

      if (doSyst == 1 || doSyst == -1 ) {
        syst_weight = weights[systName];
      }

      // Absorption
      if ( true_daughter_nPi0 == 0 ) 
        topology = 1;
      // Charge Exchange + nPi0
      else 
        topology = 2;
    }
    else continue;

   
    // Remove events with the wrong topology
    if(topology != toponum) continue;

    // True energy bin
    //if(true_beam_interactingEnergy < minval) continue;
    //if(true_beam_interactingEnergy >= maxval) continue;
    if (new_true_beam_interactingEnergy < minval ||
        new_true_beam_interactingEnergy >= maxval) continue;

    if (doNegativeReco) {
      if (reco_beam_interactingEnergy < 0.0){
        mchisto->AddBinContent(1, weight*syst_weight);        
      }
      else{
        for (int l = 1; l < nrecobins; l++) {
          if (reco_beam_interactingEnergy > recoBins[l-1] && 
              reco_beam_interactingEnergy <= recoBins[l]) {
            mchisto->AddBinContent(l+1, weight*syst_weight);
            break;
          }
        }
      }
    }
    else {
      // Sometimes the reco energy at vertex is mis-reconstructed
      if (reco_beam_interactingEnergy < 0.0) continue;

      // Fill histogram in reco energy
      for (int l = 1; l < nrecobins; l++) {
        if (reco_beam_interactingEnergy > recoBins[l-1] && 
            reco_beam_interactingEnergy <= recoBins[l]) {
          mchisto->AddBinContent(l, weight*syst_weight);
          break;
        }
      }
    }
  }

  file->Close();

  return mchisto;

}

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillDataHistogram_Pions(
    std::string filename, std::string treename, std::vector<double> recoBins,
    std::string channel, bool doNegativeReco, bool IsIncidentHisto) {
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
  TString hist_name = Form("Data_Channel%s_Histo",channel.c_str());
  TString hist_title = Form("Data for channel %s",channel.c_str());
  TH1D* datahisto = new TH1D(hist_name, hist_title,
                             (doNegativeReco ? nrecobins : nrecobins-1),
                             (doNegativeReco ? -1 : 0),
                             nrecobins-1);
  if(IsIncidentHisto)
    datahisto->SetNameTitle("Data_ChannelIncident_Histo", "Incident Data");
  datahisto->SetDirectory(0);

  mf::LogInfo("FillDataHistogram_Pions") << "Filling data histogram " << datahisto->GetName() << " from file " << filename.c_str() << " for channel " << channel.c_str();

  for(Int_t k=0; k < defaultTree->GetEntries(); k++){
    defaultTree->GetEntry(k);

    // Sometimes the reco energy at vertex is mis-reconstructed
    if (!doNegativeReco && reco_beam_interactingEnergy < 0.0) continue;

    if (reco_beam_interactingEnergy == -999.) continue;

    if (IsIncidentHisto) {
      for (size_t l = 0; l < reco_beam_incidentEnergies->size(); l++) {
        double energy = (*reco_beam_incidentEnergies)[l];
        if (doNegativeReco) {
          if (energy < 0.) {
            datahisto->AddBinContent(1, 1);
          }
          else {
	    for (int m = 1; m < nrecobins; m++) {
	      if (energy > recoBins[m-1] && energy <= recoBins[m]) {
	        datahisto->AddBinContent(m+1, 1);
	        break;
	      }
	    }
          }
        }
        else {
	  for (int m = 1; m < nrecobins; m++) {
	    if (energy > recoBins[m-1] && energy <= recoBins[m]) {
	      datahisto->AddBinContent(m, 1);
	      break;
	    }
	  }
        }
      }
      continue;
    }

    // Fill histogram in reco energy
    if (doNegativeReco) {
      if (reco_beam_interactingEnergy < 0.) {
        datahisto->AddBinContent(1, 1);
      }
      else {
        for (int l = 1; l < nrecobins; l++) {
          if (reco_beam_interactingEnergy > recoBins[l-1] &&
              reco_beam_interactingEnergy <= recoBins[l]) {
            datahisto->AddBinContent(l+1, 1);
            break;
          }
        }
      }
    }
    else {
      for (int l = 1; l < nrecobins; l++) {
        if (reco_beam_interactingEnergy > recoBins[l-1] &&
            reco_beam_interactingEnergy <= recoBins[l]) {
          datahisto->AddBinContent(l, 1);
          break;
        }
      }
    }
    
  }

  file->Close();

  return datahisto;

}

//********************************************************************
TH1* protoana::ProtoDUNESelectionUtils::FillMCIncidentHistogram_Pions(
    std::string filename, std::string treename, std::vector<double> recoBins,
    std::string topo, int toponum,
    double reco_beam_endZ_cut, double minval, double maxval,
    bool doNegativeReco, int doSyst, std::string systName, double weight) {
  //********************************************************************

  TFile *file = new TFile(filename.c_str(), "READ");
  TTree *defaultTree  = (TTree*)file->Get(treename.c_str());

  Int_t reco_beam_type; // 13 -> track-like, 11 -> shower-like
  Int_t reco_beam_nTrackDaughters, reco_beam_nShowerDaughters;
  Double_t reco_beam_len, reco_beam_vtxX, reco_beam_vtxY, reco_beam_vtxZ, reco_beam_startX, reco_beam_startY, reco_beam_startZ, reco_beam_trackDirZ, reco_beam_interactingEnergy, reco_beam_Chi2_proton;
  bool reco_beam_true_byHits_matched; // Does the true particle contributing most to the reconstructed beam track coincide with the actual beam particle that generated the event
  Int_t reco_beam_true_byHits_origin, reco_beam_true_byHits_PDG;  // Origin and PDG of the reconstructed beam track
  Int_t true_beam_PDG, true_beam_ID;
  double true_beam_endZ;
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
  defaultTree->SetBranchAddress("true_beam_endZ",                    &true_beam_endZ);
  defaultTree->SetBranchAddress("true_beam_ID",                    &true_beam_ID);
  defaultTree->SetBranchAddress("true_beam_endProcess",             &true_beam_endProcess);

  defaultTree->SetBranchAddress("true_chexSignal",                  &true_chexSignal);
  defaultTree->SetBranchAddress("true_absSignal",                   &true_absSignal);
  defaultTree->SetBranchAddress("true_backGround",                  &true_backGround);

  std::vector< int > * reco_beam_hit_true_ID = 0x0;
  std::vector< int > * reco_beam_hit_true_origin = 0x0;
  std::vector< int > * reco_beam_hit_true_slice = 0x0;

  defaultTree->SetBranchAddress("reco_beam_hit_true_ID", &reco_beam_hit_true_ID);
  defaultTree->SetBranchAddress("reco_beam_hit_true_slice", &reco_beam_hit_true_slice);
  defaultTree->SetBranchAddress("reco_beam_hit_true_origin", &reco_beam_hit_true_origin);

  std::vector< std::string > * true_beam_processes = 0x0;
  defaultTree->SetBranchAddress("true_beam_processes", &true_beam_processes);

  std::vector< int > * reco_beam_vertex_hits_slices = 0x0;
  std::vector< std::vector<double> > * reco_beam_vertex_dRs = 0x0;
  std::vector< int > * true_beam_daughter_ID = 0x0;
  std::vector< int > * true_beam_grand_daughter_ID = 0x0;
  defaultTree->SetBranchAddress( "reco_beam_vertex_hits_slices", &reco_beam_vertex_hits_slices );
  defaultTree->SetBranchAddress( "reco_beam_vertex_dRs", &reco_beam_vertex_dRs );
  defaultTree->SetBranchAddress( "true_beam_daughter_ID", &true_beam_daughter_ID );
  defaultTree->SetBranchAddress( "true_beam_grand_daughter_ID", &true_beam_grand_daughter_ID );

  //For systs
  std::vector<double> * g4rw_primary_plus_sigma_weight = 0x0;
  std::vector<double> * g4rw_primary_minus_sigma_weight = 0x0;
  std::vector<std::string> * g4rw_primary_var = 0x0;
  if (doSyst != 0) {
    defaultTree->SetBranchAddress("g4rw_primary_plus_sigma_weight", &g4rw_primary_plus_sigma_weight);
    defaultTree->SetBranchAddress("g4rw_primary_minus_sigma_weight", &g4rw_primary_minus_sigma_weight);
    defaultTree->SetBranchAddress("g4rw_primary_var", &g4rw_primary_var);
  }

  //Splitting by true energy
  std::vector<int> * true_beam_slices = 0x0;
  std::vector<double> * new_true_beam_incidentEnergies = 0x0;
  defaultTree->SetBranchAddress( "true_beam_slices", &true_beam_slices );
  defaultTree->SetBranchAddress("new_true_beam_incidentEnergies",
                                &new_true_beam_incidentEnergies);

  //std::replace(channel.begin(), channel.end(), ' ', '-');
  //channel.erase(std::remove(channel.begin(), channel.end(), '.'), channel.end());
  //channel.erase(std::remove(channel.begin(), channel.end(), ' '), channel.end());
  //std::replace(topo.begin(), topo.end(), ' ', '-');
  topo.erase(std::remove(topo.begin(), topo.end(), '.'), topo.end());
  topo.erase(std::remove(topo.begin(), topo.end(), ' '), topo.end());

  size_t nrecobins = recoBins.size();

  TString hist_name = Form("MC_ChannelIncident_%s_Histo", topo.c_str());
                      //Form("MC_ChannelIncident%s_%s_Histo",
                      //     channel.c_str(),topo.c_str());
  TString hist_title = Form("Incident MC for topology %s", topo.c_str());
                       //Form("Incident MC for channel %s and topology %s",
                       //     channel.c_str(),topo.c_str());
  
  if (doSyst == 1) {
    hist_name += "_high_" + systName;
    hist_title += " +1 sigma";
  }
  else if (doSyst == -1) {
    hist_name += "_low_" + systName;
    hist_title += " -1 sigma";
  }

  //TH1D* mchisto = new TH1D(hist_name, hist_title, nrecobins-1, 0, nrecobins-1);
  TH1D* mchisto = new TH1D(hist_name, hist_title, 
                           (doNegativeReco ? nrecobins : (nrecobins-1)),
                           (doNegativeReco ? -1 : 0), nrecobins-1);
  mchisto->SetDirectory(0);

  mf::LogInfo("FillMCIncidentHistogram_Pions") <<
      "Filling MC incident histogram " << mchisto->GetName() <<
      " from file " << filename.c_str() << /*" for channel " <<
      channel.c_str() <<*/ " with topology " << topo.c_str();

  double pitch = 0.4792;
  double z0 = 0.56035;
  int slice_cut = std::floor((reco_beam_endZ_cut - (z0 - pitch/2.)) / pitch);

  bool done_check = false;
  for(Int_t k=0; k < defaultTree->GetEntries(); k++){

    defaultTree->GetEntry(k);

    if (!done_check && doSyst != 0) {
      if (g4rw_primary_var->size() > 0) {
        std::cout << "Checking" << std::endl;
        auto syst_check = std::find(g4rw_primary_var->begin(), g4rw_primary_var->end(), systName);
        if (syst_check == g4rw_primary_var->end()) {
          std::cout << "Error! Could not find syst named " << systName << std::endl;
          std::exception e;
          throw e;
        }
        done_check = true;
      }
    }

    double syst_weight = 1.;
    std::map<std::string, double> weights;
    if (doSyst == 1) {
      for (size_t i = 0; i < g4rw_primary_var->size(); ++i) {
        weights[(*g4rw_primary_var)[i]] = (*g4rw_primary_plus_sigma_weight)[i];
      }
    }
    else if (doSyst == -1) {
      for (size_t i = 0; i < g4rw_primary_var->size(); ++i) {
        weights[(*g4rw_primary_var)[i]] = (*g4rw_primary_minus_sigma_weight)[i];
      }
    }

    // Different background topologies
    Int_t topology = -1;
 
    // Sometimes the reco energy at vertex is mis-reconstructed
    if (!doNegativeReco && reco_beam_interactingEnergy < 0.0) continue;

    if (reco_beam_interactingEnergy  == -999.) continue;

    bool check_cosmics = ( reco_beam_true_byHits_origin == 2 );

    // Go through the incident energies from the beam particle 
    for(size_t l = 0; l < reco_beam_incidentEnergies->size(); ++l){

      // Check the ID, origin, true_slice
      int true_id = (*reco_beam_hit_true_ID)[l]; 
      int true_origin = (*reco_beam_hit_true_origin)[l]; 
      int true_slice = (*reco_beam_hit_true_slice)[l]; 

      if ( check_cosmics && true_origin != 2 ) {
        std::cout << "Notice! Beam matched to cosmic, with non-cosmic hit" 
                  << std::endl;
      }                  

      double true_energy = 0.;
        

      // Cosmic
      if( true_origin == 2 ){
        topology = 5;
      }
      else{
        if( true_id != true_beam_ID ){
          if ( true_id == -999 ) 
            topology = 8; 
          else {
            if ( true_beam_endZ < 0. ) //Beam particle ends before TPC 
              topology = 7; //Just consider it downstream
            else {
              // Consider it downstream. This is slightly different than just above
              // i.e. if the beam particle interacts/decays before the TPC
              auto daughter_ID_check = std::find( 
                  true_beam_daughter_ID->begin(), true_beam_daughter_ID->end(), 
                  true_id);

              auto g_daughter_ID_check = std::find( 
                  true_beam_grand_daughter_ID->begin(), 
                  true_beam_grand_daughter_ID->end(), 
                  true_id);

              if (daughter_ID_check != true_beam_daughter_ID->end()) {
                topology = 7;
              }
              else if (g_daughter_ID_check != true_beam_grand_daughter_ID->end()){
                topology = 7; 
              }
              else {
                topology = 8;
              }

            }
          }
        }
        else{
          // True id matched to hit, but no slice matched
          // i.e. This is a 'messy' event
          if (true_slice == -999) {
            topology = 6;
          }
          else{
            if (true_beam_PDG == 211) {
              if (true_slice > slice_cut) {
                topology = 8; // Past endZ cut 
              }
              else {
                topology = 3; // Is Pion
                
                //Go through the true slices to find the 
                //true incident energy
                for (size_t i = 0; i < true_beam_slices->size(); ++i) {
                  int check_slice = (*true_beam_slices)[i];  
                  true_energy = (*new_true_beam_incidentEnergies)[i];
                  if (true_slice == check_slice) {
                    std::cout << "Reco inc energy: " <<
                                (*reco_beam_incidentEnergies)[l] << " " <<
                                "True inc energy: " << true_energy <<
                                std::endl;
                    break;
                  }
                }
              }
            }
            else if (true_beam_PDG == -13) {
              topology = 4; // Is muon
            }

            /*
            if (doSyst == 1) {
              syst_weight = g4rw_primary_plus_sigma_weight;
            }
            else if (doSyst == -1) {
              syst_weight = g4rw_primary_minus_sigma_weight;
            }
            */
            if (doSyst == 1 || doSyst == -1) { //Do +1 sigma
              syst_weight = weights[systName];
            }
          }
        }
      }


      if (topology != toponum) continue;
      if (true_energy < minval || true_energy >= maxval) {
        std::cout << "Wrong true energy bin " << minval << " " << maxval << std::endl;
        continue;
      }

/*
      for (size_t m = 1; m < nrecobins; m++) {
	if ((*reco_beam_incidentEnergies)[l] > recoBins[m-1] &&
            (*reco_beam_incidentEnergies)[l] <= recoBins[m]) {
	  mchisto->AddBinContent(m, weight*syst_weight);
	  break;
	}
      }
      */

      double energy = (*reco_beam_incidentEnergies)[l];
      if (doNegativeReco) {
        if (energy < 0.) {
          mchisto->AddBinContent(1, weight*syst_weight);
        }
        else{
          for (size_t m = 1; m < nrecobins; ++m) {
            if (energy > recoBins[m-1] &&
                energy <= recoBins[m]) {
              //Fill +1 because the first bin is negative
              mchisto->AddBinContent(m+1, weight*syst_weight);
              break;
            }
          }       
        }
      }
      else {
        // Sometimes the reco energy at vertex is mis-reconstructed
        if (energy < 0.0) continue;

        for (size_t m = 1; m < nrecobins; ++m) {
          if (energy > recoBins[m-1] &&
              energy <= recoBins[m]) {
            mchisto->AddBinContent(m/*+1*/, weight*syst_weight);
            break;
          }
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

    for(Int_t l = 1; l < ntruthbins; l++){
      if(true_beam_interactingEnergy > truthBins[l-1] && true_beam_interactingEnergy <= truthBins[l]){
	mchisto->AddBinContent(l, weight);
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
	  for(Int_t m = 1; m < nbins; m++){
	    if(true_beam_incidentEnergies->at(l) > Bins[m-1] && true_beam_incidentEnergies->at(l) <= Bins[m]){
	      mchisto->SetBinContent(m, mchisto->GetBinContent(m) + weight);
	      break;
	    }
	  }
	}
      }
      else if(mode == 3){
	for(UInt_t l = 0; l < reco_beam_incidentEnergies->size(); l++){
	  for(Int_t m = 1; m < nbins; m++){
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

int protoana::ProtoDUNESelectionUtils::GetVertexType( const std::vector< std::string > & processes, const std::vector< int > & vertex_hits_slices, const std::vector< std::vector< double > > & vertex_dRs, double cut, int max_slices ){
  
  bool matched = false;
  bool inel = false;
  bool el = false;
  bool other = false;
  
  //int matched_proc = -1;

  for( size_t i = 0; i < processes.size(); ++i ){
    std::vector< double > the_dRs;     
    for( size_t j = 0; j < vertex_dRs[i].size(); ++j ){
      int slice = vertex_hits_slices[j];
      if( slice < max_slices ) 
        the_dRs.push_back( vertex_dRs[i][j] );
      else break;
    }

    if( !the_dRs.size() ) break;

    double min_dR = 99999.;
    for( size_t j = 0; j < the_dRs.size(); ++j ){
      if( the_dRs[j] < min_dR ) 
        min_dR = the_dRs[j];
    }

    if( min_dR < cut ){
      matched = true;
      //matched_proc = i; // Will return the last proc matched
      if( processes[i] == "hadElastic" )     
        el = true;
      else if( processes[i] == "pi+Inelastic" ) 
        inel = true;
      else
        other = true;
    }

  }

  if( matched && inel && !other && !el ) 
    return 1;
  else if( matched && el && !other && !inel ) 
    return 2;
  else if( matched && el && inel && !other )  
    return 3;
  else if( matched && other ) 
    return 4;
  else 
    return 0;
  
}

std::pair< TH1 *, TH1 * > 
    protoana::ProtoDUNESelectionUtils::GetMCIncidentEfficiency(
        std::string fileName, std::string treeName, 
        std::vector< double > bins, double reco_beam_endZ_cut,
        bool doNegativeReco, int doSyst, double weight) {

  TFile * file = new TFile(fileName.c_str(), "READ");
  TTree * defaultTree  = (TTree*)file->Get(treeName.c_str());

   
  const size_t nBins = bins.size();
  TString hist_name = "MC_Incident_Efficiency_Denominator";
  TString hist_title = "Incident MC Efficiency Denominator";

  if (doSyst == 1) {
    hist_name += "_high";
    hist_title += " +1 sigma";
  }
  else if (doSyst == -1) {
    hist_name += "_low";
    hist_title += " -1 sigma";
  }

  TH1D * denominator = new TH1D(hist_name, hist_title, nBins-1, 0, nBins-1 );
  denominator->SetDirectory(0);

  hist_name = "MC_Incident_Efficiency_Numerator";
  hist_title = "Incident MC Efficiency Numerator";

  if (doSyst == 1) {
    hist_name += "_high";
    hist_title += " +1 sigma";
  }
  else if (doSyst == -1) {
    hist_name += "_low";
    hist_title += " -1 sigma";
  }

  TH1D * numerator = new TH1D(hist_name, hist_title, nBins-1, 0, nBins-1 );
  numerator->SetDirectory(0);


  double reco_beam_interactingEnergy, true_beam_endZ; 
  int true_beam_ID, true_beam_PDG;
  std::vector< int > * reco_beam_vertex_hits_slices = 0x0;
  std::vector< std::vector< double > > * reco_beam_vertex_dRs = 0x0;
  std::vector< int > * reco_beam_hit_true_ID = 0x0;
  std::vector< int > * reco_beam_hit_true_slice = 0x0;

  std::vector< std::string > * true_beam_processes = 0x0;
  std::vector< int > * true_beam_slices = 0x0;
  std::vector< double > * new_true_beam_incidentEnergies = 0x0;

  defaultTree->SetBranchAddress( "reco_beam_interactingEnergy", &reco_beam_interactingEnergy );
  defaultTree->SetBranchAddress( "reco_beam_vertex_hits_slices", &reco_beam_vertex_hits_slices );
  defaultTree->SetBranchAddress( "reco_beam_vertex_dRs", &reco_beam_vertex_dRs );
  defaultTree->SetBranchAddress( "reco_beam_hit_true_ID", &reco_beam_hit_true_ID );
  defaultTree->SetBranchAddress( "reco_beam_hit_true_slice", &reco_beam_hit_true_slice );
  defaultTree->SetBranchAddress( "true_beam_processes", &true_beam_processes );
  defaultTree->SetBranchAddress( "true_beam_slices", &true_beam_slices );
  defaultTree->SetBranchAddress( "true_beam_ID", &true_beam_ID );
  defaultTree->SetBranchAddress( "true_beam_PDG", &true_beam_PDG );
  defaultTree->SetBranchAddress( "true_beam_endZ", &true_beam_endZ );
  defaultTree->SetBranchAddress( "new_true_beam_incidentEnergies", &new_true_beam_incidentEnergies );

  double g4rw_primary_plus_sigma_weight, g4rw_primary_minus_sigma_weight;
  if (doSyst != 0) {
    defaultTree->SetBranchAddress("g4rw_primary_plus_sigma_weight", &g4rw_primary_plus_sigma_weight);
    defaultTree->SetBranchAddress("g4rw_primary_minus_sigma_weight", &g4rw_primary_minus_sigma_weight);
  }

  for( int k=0; k < defaultTree->GetEntries(); ++k ){
    double syst_weight = 1.;    
    defaultTree->GetEntry(k);

    if ( true_beam_PDG != 211 ) 
      continue;

    // Sometimes the reco energy at vertex is mis-reconstructed
    if (!doNegativeReco && reco_beam_interactingEnergy < 0.0) 
      continue;

    if (doSyst == 1) {
      syst_weight = g4rw_primary_plus_sigma_weight;
    }
    else if (doSyst == -1) {
      syst_weight = g4rw_primary_minus_sigma_weight;
    }

    if ( true_beam_endZ < 0. && true_beam_slices->size() )
      std::cout << "NOTICE: endZ < 0. but has true slices " 
          << true_beam_slices->size() << std::endl;

    //if ( true_beam_endZ > 225. ) 
    //  continue;

    double pitch = 0.4792;
    double z0 = 0.56035;
    int slice_cut = std::floor((reco_beam_endZ_cut - (z0 - pitch/2.)) / pitch);

    std::vector< int > denom_slices;

    //First, fill the denominator histogram
    //
    for (size_t i = 0; i < new_true_beam_incidentEnergies->size(); ++i) {

      int the_slice = (*true_beam_slices)[i];
      if (the_slice > slice_cut) continue;

      for ( size_t j = 1; j < nBins; ++j ) {
        if ( new_true_beam_incidentEnergies->at(i) > bins[j-1] && 
            new_true_beam_incidentEnergies->at(i) <= bins[j] ) {
          denominator->AddBinContent(j, weight*syst_weight);
          break;
        }
      }
      denom_slices.push_back(the_slice);
    }


    // Go through the incident energies from the beam particle 
    for ( size_t i = 0; i < reco_beam_hit_true_ID->size(); ++i ) {

      // Check the ID, origin, true_slice
      int true_id = (*reco_beam_hit_true_ID)[i]; 
      int true_slice = (*reco_beam_hit_true_slice)[i]; 
     
      if ( true_id == true_beam_ID ) {
        if (true_slice > slice_cut) continue;
      
        //search through the true beam slices
        for ( size_t j = 0; j < true_beam_slices->size(); ++j ) {
          if ( true_slice == (*true_beam_slices)[j] ) {

            auto slice_check = std::find(denom_slices.begin(), 
                                         denom_slices.end(), 
                                         true_slice);

            if ( slice_check == denom_slices.end() ) {
              std::cout << "Warning found true slice from hit"
                        << "that was not in denominator" << std::endl;
            }
            for ( size_t m = 1; m < nBins; ++m ) {
              if ( (*new_true_beam_incidentEnergies)[j] > bins[m-1] 
                  && (*new_true_beam_incidentEnergies)[j] <= bins[m] ) {
                numerator->AddBinContent(m, weight*syst_weight);
                break;
              }
            }
          }
        }
      }
    }
  }

  file->Close();

  return {numerator, denominator};
}

std::pair< TH1 *, TH1 *>
    protoana::ProtoDUNESelectionUtils::GetMCInteractingEfficiency(
        std::string fileName, std::string treeName,
        std::vector< double > bins, std::string channel,
        std::string topo, int toponum, double endZ_cut, int doSyst,
        double weight) {
     
  TFile * file = new TFile(fileName.c_str(), "READ");
  TTree * defaultTree  = (TTree*)file->Get(treeName.c_str());

   
  const size_t nBins = bins.size();
  TString hist_name = Form("MC_Channel%s_%s_Interacting_Denominator",channel.c_str(),topo.c_str());
  TString hist_title = Form( "Interacting MC Efficiency Denominator for channel %s and topology %s", channel.c_str(), topo.c_str());

  if (doSyst == 1) {
    hist_name += "_high";
    hist_title += " +1 sigma";
  }
  else if (doSyst == -1) {
    hist_name += "_low";
    hist_title += " -1 sigma";
  }

  TH1D * denominator = new TH1D(hist_name, hist_title, nBins-1, 0, nBins-1);
  denominator->SetDirectory(0);

  hist_name = Form("MC_Channel%s_%s_Interacting_Numerator",channel.c_str(),topo.c_str());
  hist_title = Form("Interacting MC Efficiency Numerator for channel %s and topology %s", channel.c_str(), topo.c_str());

  if (doSyst == 1) {
    hist_name += "_high";
    hist_title += " +1 sigma";
  }
  else if (doSyst == -1) {
    hist_name += "_low";
    hist_title += " -1 sigma";
  }


  TH1D * numerator = new TH1D(hist_name, hist_title, nBins-1, 0, nBins-1 );
  numerator->SetDirectory(0);

  //Init stuff here
  int true_beam_PDG, true_daughter_nPiPlus, true_daughter_nPiMinus, true_daughter_nPi0;
  std::string * true_beam_endProcess = 0x0;
  std::vector< std::string > * true_beam_processes = 0x0;
  std::vector< int > * reco_beam_vertex_hits_slices = 0x0;
  std::vector< std::vector< double > > * reco_beam_vertex_dRs = 0x0;
  double new_true_beam_interactingEnergy, true_beam_endZ;
  bool has_noPion_daughter, has_shower_nHits_distance;

  defaultTree->SetBranchAddress( "true_beam_PDG", &true_beam_PDG );
  defaultTree->SetBranchAddress( "true_beam_endZ", &true_beam_endZ );
  defaultTree->SetBranchAddress( "true_beam_endProcess", &true_beam_endProcess );
  defaultTree->SetBranchAddress( "true_daughter_nPiPlus", &true_daughter_nPiPlus );
  defaultTree->SetBranchAddress( "true_daughter_nPiMinus", &true_daughter_nPiMinus );
  defaultTree->SetBranchAddress( "true_daughter_nPi0", &true_daughter_nPi0 );
  defaultTree->SetBranchAddress( "true_beam_processes", &true_beam_processes );
  defaultTree->SetBranchAddress( "reco_beam_vertex_hits_slices", &reco_beam_vertex_hits_slices );
  defaultTree->SetBranchAddress( "reco_beam_vertex_dRs", &reco_beam_vertex_dRs );
  defaultTree->SetBranchAddress( "new_true_beam_interactingEnergy", &new_true_beam_interactingEnergy );
  defaultTree->SetBranchAddress( "has_noPion_daughter", &has_noPion_daughter );
  defaultTree->SetBranchAddress( "has_shower_nHits_distance", &has_shower_nHits_distance );

  double g4rw_primary_plus_sigma_weight, g4rw_primary_minus_sigma_weight;
  if (doSyst != 0) {
    defaultTree->SetBranchAddress("g4rw_primary_plus_sigma_weight", &g4rw_primary_plus_sigma_weight);
    defaultTree->SetBranchAddress("g4rw_primary_minus_sigma_weight", &g4rw_primary_minus_sigma_weight);
  }

  int topology = -1;

  for( int k=0; k < defaultTree->GetEntries(); ++k ){
    double syst_weight = 1.;
    defaultTree->GetEntry(k);

    //First, check if it's a pion
    if ( true_beam_PDG != 211 ) 
      continue;
    
    if ( true_beam_endZ < 0. ) 
      continue;

    //Then if it ends in abs/cex/nPi0 in the FV
    if (!(*true_beam_endProcess == "pi+Inelastic" &&
          true_daughter_nPiPlus == 0 && true_daughter_nPiMinus == 0 &&
          true_beam_endZ < endZ_cut)){
      continue;
    }

    //Abs
    if (true_daughter_nPi0 == 0) {
      topology = 1;
    }
    //Cex/nPi0
    else {
      topology = 2;
    }

    if (topology != toponum){
      continue;
    }

    if (doSyst == 1) {
      syst_weight = g4rw_primary_plus_sigma_weight;
    }
    else if (doSyst == -1) {
      syst_weight = g4rw_primary_minus_sigma_weight;
    }

    //Fill the denominator
    for ( size_t i = 1; i < nBins; ++i ) {
      if (new_true_beam_interactingEnergy > bins[i-1] 
          && new_true_beam_interactingEnergy <= bins[i]){
        denominator->AddBinContent(i, weight*syst_weight);

        if ( has_noPion_daughter ) { //Abs/Cex selection
          if (topology == 1 && !has_shower_nHits_distance) {
            numerator->AddBinContent(i, weight*syst_weight); //Abs
          }
          else if (topology == 2 && has_shower_nHits_distance) {
            numerator->AddBinContent(i, weight*syst_weight); //Cex
          }
        } 
        break;

      }
    }
  }


  file->Close();
  return {numerator, denominator};
}
