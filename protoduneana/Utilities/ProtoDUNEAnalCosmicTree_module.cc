////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEAnalCosmicTree
// File:        ProtoDUNEAnalCosmicTree_module.cc
//
// Extract protoDUNE useful information, do a firs tpre-selection and save output to a flat tree
// 
// Georgios Christodoulou - georgios.christodoulou at cern.ch
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TTimeStamp.h"

// C++ Includes
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

// Maximum number of beam particles to save
const int NMAXCOSMICPARTICLES = 100;

namespace protoana {
  class ProtoDUNEAnalCosmicTree;
}


class protoana::ProtoDUNEAnalCosmicTree : public art::EDAnalyzer {
public:

  explicit ProtoDUNEAnalCosmicTree(fhicl::ParameterSet const & p);

  ProtoDUNEAnalCosmicTree(ProtoDUNEAnalCosmicTree const &) = delete;
  ProtoDUNEAnalCosmicTree(ProtoDUNEAnalCosmicTree &&) = delete;
  ProtoDUNEAnalCosmicTree & operator = (ProtoDUNEAnalCosmicTree const &) = delete;
  ProtoDUNEAnalCosmicTree & operator = (ProtoDUNEAnalCosmicTree &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & evt) override;

private:

  // Helper utility functions
  protoana::ProtoDUNEDataUtils dataUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETrackUtils trackUtil;
  protoana::ProtoDUNETruthUtils truthUtil;

  // Track momentum algorithm calculates momentum based on track range
  trkf::TrackMomentumCalculator trmom;

  // Initialise tree variables
  void Initialise();

  // Fill cosmics tree
  void FillCosmicsTree(art::Event const & evt);

  // fcl parameters
  const art::InputTag fBeamModuleLabel;
  std::string fCalorimetryTag;
  std::string fParticleIDTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  double fMinimumTrackLength;

  TTree *fPandoraCosmics;

  // Tree variables
  int fRun;
  int fSubRun;
  int fevent;
  double fTimeStamp;
  int fNactivefembs[6];

  // Cosmic tree variables 
  int fNCOSMICS;
  int fcosmicIsClear[NMAXCOSMICPARTICLES];
  int fcosmicNHits[NMAXCOSMICPARTICLES];
  double fcosmicBDTScore[NMAXCOSMICPARTICLES];
  int fiscosmictrack[NMAXCOSMICPARTICLES];
  int fiscosmicshower[NMAXCOSMICPARTICLES];
  double fcosmicID[NMAXCOSMICPARTICLES];
  double fcosmicTheta[NMAXCOSMICPARTICLES];
  double fcosmicPhi[NMAXCOSMICPARTICLES];
  double fcosmicLength[NMAXCOSMICPARTICLES];
  //double fcosmicMomentum[NMAXCOSMICPARTICLES];
  //double fcosmicEndMomentum[NMAXCOSMICPARTICLES];
  double fcosmicMomentumByRangeMuon[NMAXCOSMICPARTICLES];
  double fcosmicMomentumByRangeProton[NMAXCOSMICPARTICLES];
  int fcosmicShowerBestPlane[NMAXCOSMICPARTICLES];
  double  fcosmicOpeningAngle[NMAXCOSMICPARTICLES];
    
  double fcosmicEndPosition[NMAXCOSMICPARTICLES][3];
  double fcosmicStartPosition[NMAXCOSMICPARTICLES][3];
  double fcosmicEndDirection[NMAXCOSMICPARTICLES][3];
  double fcosmicStartDirection[NMAXCOSMICPARTICLES][3];
  double fcosmicKineticEnergy[NMAXCOSMICPARTICLES][3];
  double fcosmicRange[NMAXCOSMICPARTICLES][3];
  double fcosmicTrkPitchC[NMAXCOSMICPARTICLES][3];

  double fcosmicT0[NMAXCOSMICPARTICLES];

  int fcosmicPID_Pdg[NMAXCOSMICPARTICLES][3];       
  int fcosmicPID_Ndf[NMAXCOSMICPARTICLES][3];       
  double fcosmicPID_MinChi2[NMAXCOSMICPARTICLES][3];   
  double fcosmicPID_DeltaChi2[NMAXCOSMICPARTICLES][3]; 
  double fcosmicPID_Chi2Proton[NMAXCOSMICPARTICLES][3];
  double fcosmicPID_Chi2Kaon[NMAXCOSMICPARTICLES][3];  
  double fcosmicPID_Chi2Pion[NMAXCOSMICPARTICLES][3];  
  double fcosmicPID_Chi2Muon[NMAXCOSMICPARTICLES][3];  
  double fcosmicPID_MissingE[NMAXCOSMICPARTICLES][3];  
  double fcosmicPID_MissingEavg[NMAXCOSMICPARTICLES][3];
  double fcosmicPID_PIDA[NMAXCOSMICPARTICLES][3];  

};


protoana::ProtoDUNEAnalCosmicTree::ProtoDUNEAnalCosmicTree(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fParticleIDTag(p.get<std::string>("ParticleIDTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fMinimumTrackLength(p.get<double>("MinimumTrackLength"))
{

}

void protoana::ProtoDUNEAnalCosmicTree::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;

  fPandoraCosmics = tfs->make<TTree>("PandoraCosmics", "Cosmic tracks reconstructed with Pandora");
  fPandoraCosmics->Branch("run",                          &fRun,                           "run/I");
  fPandoraCosmics->Branch("subrun",                       &fSubRun,                        "subrun/I");
  fPandoraCosmics->Branch("event",                        &fevent,                         "event/I");
  fPandoraCosmics->Branch("timestamp",                    &fTimeStamp,                     "timestamp/D");
  fPandoraCosmics->Branch("Nactivefembs",                 &fNactivefembs,                  "Nactivefembs[5]/I");
  fPandoraCosmics->Branch("NCOSMICS",                     &fNCOSMICS,                      "NCOSMICS/I");
  fPandoraCosmics->Branch("cosmicIsClear",                &fcosmicIsClear,                 "cosmicIsClear[NCOSMICS]/I");
  fPandoraCosmics->Branch("cosmicBDTScore",               &fcosmicBDTScore,                "cosmicBDTScore[NCOSMICS]/D");
  fPandoraCosmics->Branch("cosmicNHits",                  &fcosmicNHits,                   "cosmicNHits[NCOSMICS]/I");
  fPandoraCosmics->Branch("iscosmictrack",                &fiscosmictrack,                 "iscosmictrack[NCOSMICS]/I");
  fPandoraCosmics->Branch("iscosmicshower",               &fiscosmicshower,                "iscosmicshower[NCOSMICS]/I");
  fPandoraCosmics->Branch("cosmicID",                     &fcosmicID,                      "cosmicID[NCOSMICS]/I");
  fPandoraCosmics->Branch("cosmicTheta",                  &fcosmicTheta,                   "cosmicTheta[NCOSMICS]/D");
  fPandoraCosmics->Branch("cosmicPhi",                    &fcosmicPhi,                     "cosmicPhi[NCOSMICS]/D");
  fPandoraCosmics->Branch("cosmicLength",                 &fcosmicLength,                  "cosmicLength[NCOSMICS]/D");
  //fPandoraCosmics->Branch("cosmicMomentum",               &fcosmicMomentum,              "cosmicMomentum[NCOSMICS]/D");
  //fPandoraCosmics->Branch("cosmicEndMomentum",            &fcosmicEndMomentum,           "cosmicEndMomentum[NCOSMICS]/D");
  fPandoraCosmics->Branch("cosmicOpeningAngle",           &fcosmicOpeningAngle,             "cosmicOpeningAngle[NCOSMICS]/D");
  fPandoraCosmics->Branch("cosmicShowerBestPlane",        &fcosmicShowerBestPlane,          "cosmicShowerBestPlane[NCOSMICS]/I");
  fPandoraCosmics->Branch("cosmicEndPosition",            &fcosmicEndPosition,              "cosmicEndPosition[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicStartPosition",          &fcosmicStartPosition,            "cosmicStartPosition[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicEndDirection",           &fcosmicEndDirection,             "cosmicEndDirection[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicStartDirection",         &fcosmicStartDirection,           "cosmicStartDirection[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicMomentumByRangeMuon",    &fcosmicMomentumByRangeMuon,      "cosmicMomentumByRangeMuon[NCOSMICS]/D");
  fPandoraCosmics->Branch("cosmicMomentumByRangeProton",  &fcosmicMomentumByRangeProton,    "cosmicMomentumByRangeProton[NCOSMICS]/D");
  fPandoraCosmics->Branch("cosmicKineticEnergy",          &fcosmicKineticEnergy,            "cosmicKineticEnergy[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicRange",                  &fcosmicRange,                    "cosmicRange[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicTrkPitchC",              &fcosmicTrkPitchC,                "cosmicTrkPitchC[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicT0",                     &fcosmicT0,                       "cosmicT0[NCOSMICS]/D");

  fPandoraCosmics->Branch("cosmicPID_Pdg",                &fcosmicPID_Pdg,                  "cosmicPID_Pdg[NCOSMICS][3]/I");
  fPandoraCosmics->Branch("cosmicPID_Ndf",                &fcosmicPID_Ndf,                  "cosmicPID_Ndf[NCOSMICS][3]/I");
  fPandoraCosmics->Branch("cosmicPID_MinChi2",            &fcosmicPID_MinChi2,              "cosmicPID_MinChi2[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicPID_DeltaChi2",          &fcosmicPID_DeltaChi2,            "cosmicPID_DeltaChi2[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicPID_Chi2Proton",         &fcosmicPID_Chi2Proton,           "cosmicPID_Chi2Proton[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicPID_Chi2Kaon",           &fcosmicPID_Chi2Kaon,             "cosmicPID_Chi2Kaon[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicPID_Chi2Pion",           &fcosmicPID_Chi2Pion,             "cosmicPID_Chi2Pion[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicPID_Chi2Muon",           &fcosmicPID_Chi2Muon,             "cosmicPID_Chi2Muon[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicPID_MissingE",           &fcosmicPID_MissingE,             "cosmicPID_MissingE[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicPID_MissingEavg",        &fcosmicPID_MissingEavg,          "cosmicPID_MissingEavg[NCOSMICS][3]/D");
  fPandoraCosmics->Branch("cosmicPID_PIDA",               &fcosmicPID_PIDA,                 "cosmicPID_PIDA[NCOSMICS][3]/D");

}

void protoana::ProtoDUNEAnalCosmicTree::analyze(art::Event const & evt){

  // Initialise tree parameters
  Initialise();

  fRun = evt.run();
  fSubRun = evt.subRun();
  fevent  = evt.id().event(); 
  art::Timestamp ts = evt.time();
  if (ts.timeHigh() == 0){
    TTimeStamp ts2(ts.timeLow());
    fTimeStamp = ts2.AsDouble();
  }
  else{
    TTimeStamp ts2(ts.timeHigh(), ts.timeLow());
    fTimeStamp = ts2.AsDouble();
  }
  
  // Get number of active fembs
  int allactivefembs = 0;
  if(!evt.isRealData()){
    allactivefembs = 120;
    for(int k=0; k < 6; k++)
      fNactivefembs[0] = 20;
  }
  else{
    for(int k=0; k < 6; k++){
     fNactivefembs[k] = dataUtil.GetNActiveFembsForAPA(evt, k);
     allactivefembs += fNactivefembs[k];
    }
  }

  //trmom.SetMinLength(100);
  
  // For cosmics only save events if the all fembs are active
  if(allactivefembs != 120) return;
  
  FillCosmicsTree(evt);

}

void protoana::ProtoDUNEAnalCosmicTree::endJob(){

}

void protoana::ProtoDUNEAnalCosmicTree::FillCosmicsTree(art::Event const & evt){

  unsigned int npfParticles = pfpUtil.GetNumberPrimaryPFParticle(evt,fPFParticleTag);

  if(npfParticles == 0){
    std::cout << "WARNING::No PfParticles found! Skipping event!" << std::endl;
    return;
  }

  // Do not process more than NMAXCOSMICPARTICLES
  if(fNCOSMICS > NMAXCOSMICPARTICLES) return;

  // Get all PFParticles
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  for(unsigned int i = 0; i < recoParticles->size(); i++){
    recob::PFParticle* pfparticle = const_cast<recob::PFParticle*>(&(recoParticles->at(i)));
    //const recob::PFParticle* pfparticle = &(recoParticles->at(i));

    //  Only consider primary particles
    if(!pfparticle->IsPrimary()) continue;
    
    // Do not consider beam particles
    if(pfpUtil.IsBeamParticle(*pfparticle, evt, fPFParticleTag)) continue;

    // Get the T0 for this pfParticle
    std::vector<anab::T0> pfT0vec = pfpUtil.GetPFParticleT0(*pfparticle,evt,fPFParticleTag);
    //if(pfT0vec.empty()) continue;

    const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*pfparticle,evt,fPFParticleTag,fTrackerTag);
    const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*pfparticle,evt,fPFParticleTag,fShowerTag);

    // Do not save very short cosmic tracks
    if(thisTrack != 0x0  && thisTrack->Length()  < fMinimumTrackLength) continue;
    if(thisShower != 0x0 && thisShower->Length() < fMinimumTrackLength) continue;

    // Pandora's BDT beam-cosmic score
    fcosmicBDTScore[fNCOSMICS] = (double)pfpUtil.GetBeamCosmicScore(*pfparticle,evt,fPFParticleTag);
    
    // NHits associated with this pfParticle
    fcosmicNHits[fNCOSMICS] = (pfpUtil.GetPFParticleHits(*pfparticle,evt,fPFParticleTag)).size();

    if(!pfT0vec.empty())
      fcosmicT0[fNCOSMICS] = pfT0vec[0].Time();

    // Check if it is clear cosmic
    fcosmicIsClear[fNCOSMICS] = pfpUtil.IsClearCosmic(*pfparticle,evt,fPFParticleTag);

    if(thisTrack != 0x0){
      fiscosmictrack[fNCOSMICS]               = 1;
      fiscosmicshower[fNCOSMICS]              = 0;
      fcosmicID[fNCOSMICS]                    = thisTrack->ParticleId();
      fcosmicTheta[fNCOSMICS]                 = thisTrack->Theta();
      fcosmicPhi[fNCOSMICS]                   = thisTrack->Phi();
      fcosmicLength[fNCOSMICS]                = thisTrack->Length();
      //fcosmicMomentum[fNCOSMICS]              = thisTrack->StartMomentum();
      //fcosmicEndMomentum[fNCOSMICS]           = thisTrack->EndMomentum();
      fcosmicEndPosition[fNCOSMICS][0]        = thisTrack->Trajectory().End().X();
      fcosmicEndPosition[fNCOSMICS][1]        = thisTrack->Trajectory().End().Y();
      fcosmicEndPosition[fNCOSMICS][2]        = thisTrack->Trajectory().End().Z();
      fcosmicStartPosition[fNCOSMICS][0]      = thisTrack->Trajectory().Start().X();
      fcosmicStartPosition[fNCOSMICS][1]      = thisTrack->Trajectory().Start().Y();
      fcosmicStartPosition[fNCOSMICS][2]      = thisTrack->Trajectory().Start().Z();
      fcosmicEndDirection[fNCOSMICS][0]       = thisTrack->Trajectory().EndDirection().X();
      fcosmicEndDirection[fNCOSMICS][1]       = thisTrack->Trajectory().EndDirection().Y();
      fcosmicEndDirection[fNCOSMICS][2]       = thisTrack->Trajectory().EndDirection().Z();
      fcosmicStartDirection[fNCOSMICS][0]     = thisTrack->Trajectory().StartDirection().X();
      fcosmicStartDirection[fNCOSMICS][1]     = thisTrack->Trajectory().StartDirection().Y();
      fcosmicStartDirection[fNCOSMICS][2]     = thisTrack->Trajectory().StartDirection().Z();

      fcosmicMomentumByRangeMuon[fNCOSMICS]   = trmom.GetTrackMomentum(thisTrack->Length(),13);
      fcosmicMomentumByRangeProton[fNCOSMICS] = trmom.GetTrackMomentum(thisTrack->Length(),2212);

      // Calorimetry
      std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
      if(calovector.size() != 3)
	std::cerr << "WARNING::Calorimetry vector size for cosmic is = " << calovector.size() << std::endl;

      for(size_t k = 0; k < calovector.size() && k<3; k++){
	int plane = calovector[k].PlaneID().Plane;
	if(plane < 0) continue;
	if(plane > 2) continue;
	fcosmicKineticEnergy[fNCOSMICS][plane] = calovector[k].KineticEnergy();
	fcosmicRange[fNCOSMICS][plane]         = calovector[k].Range();
	fcosmicTrkPitchC[fNCOSMICS][plane]     = calovector[k].TrkPitchC();
      }

      // PID
      std::vector<anab::ParticleID> pids = trackUtil.GetRecoTrackPID(*thisTrack, evt, fTrackerTag, fParticleIDTag);
      if(pids.size() != 3)
	std::cerr << "WARNING::PID vector size for primary is = " << pids.size() << std::endl;
    
      for(size_t k = 0; k < pids.size() && k<3; k++){
	int plane = pids[k].PlaneID().Plane;
	if(plane < 0) continue;
	if(plane > 2) continue;
	fcosmicPID_Pdg[fNCOSMICS][plane]            = pids[plane].Pdg();
	fcosmicPID_Ndf[fNCOSMICS][plane]            = pids[plane].Ndf();
	fcosmicPID_MinChi2[fNCOSMICS][plane]        = pids[plane].MinChi2();
	fcosmicPID_DeltaChi2[fNCOSMICS][plane]      = pids[plane].DeltaChi2();
	fcosmicPID_Chi2Proton[fNCOSMICS][plane]     = pids[plane].Chi2Proton();
	fcosmicPID_Chi2Kaon[fNCOSMICS][plane]       = pids[plane].Chi2Kaon();
	fcosmicPID_Chi2Pion[fNCOSMICS][plane]       = pids[plane].Chi2Pion();
	fcosmicPID_Chi2Muon[fNCOSMICS][plane]       = pids[plane].Chi2Muon();
	fcosmicPID_MissingE[fNCOSMICS][plane]       = pids[plane].MissingE();
	fcosmicPID_MissingEavg[fNCOSMICS][plane]    = pids[plane].MissingEavg();
	fcosmicPID_PIDA[fNCOSMICS][plane]           = pids[plane].PIDA();
      }

    }
    else if(thisShower != 0x0){
      fiscosmictrack[fNCOSMICS]               = 0;
      fiscosmicshower[fNCOSMICS]              = 1;

      fcosmicID[fNCOSMICS]                    = thisShower->ID();
      fcosmicLength[fNCOSMICS]                = thisShower->Length();
      fcosmicShowerBestPlane[fNCOSMICS]       = thisShower->best_plane();
      fcosmicOpeningAngle[fNCOSMICS]          = thisShower->OpenAngle();
      fcosmicStartPosition[fNCOSMICS][0]      = thisShower->ShowerStart().X();
      fcosmicStartPosition[fNCOSMICS][1]      = thisShower->ShowerStart().Y();
      fcosmicStartPosition[fNCOSMICS][2]      = thisShower->ShowerStart().Z();
      fcosmicStartDirection[fNCOSMICS][0]     = thisShower->Direction().X();
      fcosmicStartDirection[fNCOSMICS][1]     = thisShower->Direction().Y();
      fcosmicStartDirection[fNCOSMICS][2]     = thisShower->Direction().Z();
    }
    else{
      std::cout << "INFO::Cosmic pfParticle is not track or shower. Skip!" << std::endl;
      continue;
    }

    // Increment number of cosmic tracks
    fNCOSMICS++;

  }

  fPandoraCosmics->Fill();

}

void protoana::ProtoDUNEAnalCosmicTree::Initialise(){
  
  fRun = -999;
  fSubRun = -999;
  fevent = -999;
  fTimeStamp = -999.0;
  for(int k=0; k < 5; k++)
    fNactivefembs[k] = -999;

  // Cosmics tree
  fNCOSMICS = 0;
  for(int k=0; k < NMAXCOSMICPARTICLES; k++){
    fcosmicIsClear[k]               = -999;
    fcosmicNHits[k]                 = -999;
    fcosmicBDTScore[k]              = -999.0;
    fiscosmictrack[k]               = -999;
    fiscosmicshower[k]              = -999;
    fcosmicID[k]                    = -999;
    fcosmicTheta[k]                 = -999.0;
    fcosmicPhi[k]                   = -999.0;
    fcosmicLength[k]                = -999.0;
    //fcosmicMomentum[k]              = -999.0;
    //fcosmicEndMomentum[k]           = -999.0;
    fcosmicMomentumByRangeMuon[k]   = -999.0;
    fcosmicMomentumByRangeProton[k] = -999.0;
    fcosmicT0[k]                    = -999.0;
    fcosmicShowerBestPlane[k]       = -999.0;
    fcosmicOpeningAngle[k]          = -999.0;
    
    for(int l=0; l < 3; l++){
      fcosmicEndPosition[k][l]        = -999.0;
      fcosmicStartPosition[k][l]      = -999.0;
      fcosmicEndDirection[k][l]       = -999.0;
      fcosmicStartDirection[k][l]     = -999.0;
      fcosmicKineticEnergy[k][l]      = -999.0;
      fcosmicRange[k][l]              = -999.0;
      fcosmicTrkPitchC[k][l]          = -999.0;

      fcosmicPID_Pdg[k][l]            = -999.0;
      fcosmicPID_Ndf[k][l]            = -999.0;
      fcosmicPID_MinChi2[k][l]        = -999.0;
      fcosmicPID_DeltaChi2[k][l]      = -999.0;
      fcosmicPID_Chi2Proton[k][l]     = -999.0;
      fcosmicPID_Chi2Kaon[k][l]       = -999.0;
      fcosmicPID_Chi2Pion[k][l]       = -999.0;
      fcosmicPID_Chi2Muon[k][l]       = -999.0;
      fcosmicPID_MissingE[k][l]       = -999.0;
      fcosmicPID_MissingEavg[k][l]    = -999.0;
      fcosmicPID_PIDA[k][l]           = -999.0;
    }
  }

}

DEFINE_ART_MODULE(protoana::ProtoDUNEAnalCosmicTree)

