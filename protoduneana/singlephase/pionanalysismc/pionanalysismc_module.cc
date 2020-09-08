////////////////////////////////////////////////////////////////////////
// Class:       pionanalysismc
// File:        pionanalysismc_module.cc
//
// Extract protoDUNE useful information, do a firs tpre-selection and save output to a flat tree
// 
// Some parts are copied from the beam module example
//
// Georgios Christodoulou - georgios.christodoulou at cern.ch
// Heng-Ye Liao modified for his protonanalysis - liao@phys.ksu.edu
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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
//#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEbeamsim.h"
//#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEBeamInstrument.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"

#include "larsim/MCCheater/ParticleInventoryService.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TTimeStamp.h"

// C++ Includes
#include <stdio.h>
#include <stdlib.h> 
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

// Maximum number of beam particles to save
const int NMAXDAUGTHERS = 15;

//const int kMaxTracks  = 1000;
//const int kMaxHits = 10000;

namespace protoana {
  class pionanalysismc;
}


class protoana::pionanalysismc : public art::EDAnalyzer {
public:

  explicit pionanalysismc(fhicl::ParameterSet const & p);

  pionanalysismc(pionanalysismc const &) = delete;
  pionanalysismc(pionanalysismc &&) = delete;
  pionanalysismc & operator = (pionanalysismc const &) = delete;
  pionanalysismc & operator = (pionanalysismc &&) = delete;

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

  // Initialise tree variables
  void Initialise();

  // Fill cosmics tree
  void FillCosmicsTree(art::Event const & evt, std::string pfParticleTag);

  // fcl parameters
  const art::InputTag fBeamModuleLabel;
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  bool fVerbose;

  TTree *fPandoraBeam;
  //TTree *fPandoraCosmics;

  // Tree variables
  int fRun;
  int fSubRun;
  int fevent;
  double fTimeStamp;
  int fNactivefembs[6];//why is the size=6???????

  ////// Beam track information
  int fbeamtrigger;
  double ftof;
  int fcerenkov1;
  int fcerenkov2;
  double fbeamtrackMomentum;
  double fbeamtrackP[3]; //Px/Py/Pz
  double fbeamtrackEnergy;
  double fbeamtrackPos[3];
  double fbeamtrackDir[3];
  double fbeamtrackTime;
  int fbeamtrackPdg;
  int fbeamtrackID;

 // Reconstructed tracks/showers information
  double fvertex[3];
  double fsecvertex[3];

  int fisprimarytrack;
  int fisprimaryshower;
  double fprimaryBDTScore;
  int fprimaryNHits;
  double fprimaryTheta;
  double fprimaryPhi;
  double fprimaryLength;
  double fprimaryMomentum;
  double fprimaryEndMomentum;
  double fprimaryEndPosition[3];
  double fprimaryStartPosition[3];
  double fprimaryEndDirection[3];
  double fprimaryStartDirection[3];
  double fprimaryOpeningAngle;
  int fprimaryShowerBestPlane;
  double fprimaryShowerEnergy;
  double fprimaryShowerMIPEnergy;
  double fprimaryShowerdEdx;
  double fprimaryMomentumByRangeProton;
  double fprimaryMomentumByRangeMuon;
  double fprimaryKineticEnergy[3];
  double fprimaryRange[3];
  int fprimaryID;
  double fprimaryT0;



  //*********************truth information*******************************************//
  int fprimary_truth_TrackId;
  int fprimary_truth_Pdg;
  double fprimary_truth_StartPosition[4];
  double fprimary_truth_EndPosition[4];
  std::string fprimary_truth_EndProcess;
  double fprimary_truth_P;
  double fprimary_truth_Momentum[4];
  double fprimary_truth_EndMomentum[4];
  double fprimary_truth_Pt;
  double fprimary_truth_Mass;
  double fprimary_truth_Theta;
  double fprimary_truth_Phi;
  int fprimary_truth_Process;
  int fprimary_truth_Isbeammatched;
  int fprimary_truth_NDaughters;

  //calo info
  std::vector< std::vector<double> > primtrk_dqdx;
  std::vector< std::vector<double> > primtrk_resrange;
  std::vector< std::vector<double> > primtrk_dedx;
  std::vector<double> primtrk_range;
  std::vector< std::vector<double> > primtrk_hitx;
  std::vector< std::vector<double> > primtrk_hity;
  std::vector< std::vector<double> > primtrk_hitz;
  std::vector< std::vector<double> > primtrk_pitch;
  //
  //MC Truth information
  /*  std::vector<double> truth_startx;
  std::vector<double> truth_starty;
  std::vector<double> truth_startz;
  std::vector<double> truth_endx;
  std::vector<double> truth_endy;
  std::vector<double> truth_endz;
  std::vector<double> truth_Ndaughter;
  std::vector<double> truth_pdg;
  std::vector<string> truth_endprocess;
  std::vector<string> truth_process;*/

  //////////////////////////////



  // Daughters from primary
  int fNDAUGHTERS;
  int fisdaughtertrack[NMAXDAUGTHERS];
  int fisdaughtershower[NMAXDAUGTHERS];
  int fdaughterNHits[NMAXDAUGTHERS];
  double fdaughterTheta[NMAXDAUGTHERS];
  double fdaughterPhi[NMAXDAUGTHERS];
  double fdaughterLength[NMAXDAUGTHERS];
  double fdaughterMomentum[NMAXDAUGTHERS];
  double fdaughterEndMomentum[NMAXDAUGTHERS];
  double fdaughterEndPosition[NMAXDAUGTHERS][3];
  double fdaughterStartPosition[NMAXDAUGTHERS][3];
  double fdaughterEndDirection[NMAXDAUGTHERS][3];
  double fdaughterStartDirection[NMAXDAUGTHERS][3];
  double fdaughterOpeningAngle[NMAXDAUGTHERS];
  double fdaughterShowerEnergy[NMAXDAUGTHERS];
  double fdaughterShowerMIPEnergy[NMAXDAUGTHERS];
  double fdaughterShowerdEdx[NMAXDAUGTHERS];
  int fdaughterShowerBestPlane[NMAXDAUGTHERS];
  double fdaughterMomentumByRangeProton[NMAXDAUGTHERS];
  double fdaughterMomentumByRangeMuon[NMAXDAUGTHERS];
  double fdaughterKineticEnergy[NMAXDAUGTHERS][3];
  double fdaughterRange[NMAXDAUGTHERS][3];
  int fdaughterID[NMAXDAUGTHERS];
  double fdaughterT0[NMAXDAUGTHERS];

  int fdaughter_truth_TrackId[NMAXDAUGTHERS];
  int fdaughter_truth_Pdg[NMAXDAUGTHERS];
  double fdaughter_truth_StartPosition[NMAXDAUGTHERS][4];
  double fdaughter_truth_EndPosition[NMAXDAUGTHERS][4];
  double fdaughter_truth_P[NMAXDAUGTHERS];
  double fdaughter_truth_Momentum[NMAXDAUGTHERS][4];
  double fdaughter_truth_EndMomentum[NMAXDAUGTHERS][4];
  double fdaughter_truth_Pt[NMAXDAUGTHERS];
  double fdaughter_truth_Mass[NMAXDAUGTHERS];
  double fdaughter_truth_Theta[NMAXDAUGTHERS];
  double fdaughter_truth_Phi[NMAXDAUGTHERS];
  int fdaughter_truth_Process[NMAXDAUGTHERS];

};


protoana::pionanalysismc::pionanalysismc(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fVerbose(p.get<bool>("Verbose"))
{

}

void protoana::pionanalysismc::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;

  fPandoraBeam = tfs->make<TTree>("PandoraBeam", "Beam events reconstructed with Pandora");
  fPandoraBeam->Branch("run",                           &fRun,                          "run/I");
  fPandoraBeam->Branch("subrun",                        &fSubRun,                       "subrun/I");
  fPandoraBeam->Branch("event",                         &fevent,                        "event/I");
  fPandoraBeam->Branch("timestamp",                     &fTimeStamp,                    "timestamp/D");
  fPandoraBeam->Branch("Nactivefembs",                  &fNactivefembs,                 "Nactivefembs[5]/I");

  fPandoraBeam->Branch("beamtrigger",                   &fbeamtrigger,                  "beamtrigger/I");
  fPandoraBeam->Branch("tof",                           &ftof,                          "tof/D");
  fPandoraBeam->Branch("cerenkov1",                     &fcerenkov1,                    "cerenkov1/I");
  fPandoraBeam->Branch("cerenkov2",                     &fcerenkov2,                    "cerenkov2/I");
  fPandoraBeam->Branch("beamtrackMomentum",             &fbeamtrackMomentum,            "beamtrackMomentum/D");
  fPandoraBeam->Branch("beamtrackP",                    &fbeamtrackP,                   "beamtrackP[3]/D");
  fPandoraBeam->Branch("beamtrackEnergy",               &fbeamtrackEnergy,              "beamtrackEnergy/D");
  fPandoraBeam->Branch("beamtrackPos",                  &fbeamtrackPos,                 "beamtrackPos[3]/D");
  fPandoraBeam->Branch("beamtrackDir",                  &fbeamtrackDir,                 "beamtrackDir[3]/D");
  fPandoraBeam->Branch("beamtrackTime",                 &fbeamtrackTime,                "beamtrackTime/D");
  fPandoraBeam->Branch("beamtrackPdg",                  &fbeamtrackPdg,                 "beamtrackPdg/I");
  fPandoraBeam->Branch("beamtrackID",                   &fbeamtrackID,                  "beamtrackID/I");

  fPandoraBeam->Branch("vertex",                        &fvertex,                       "vertex[3]/D");
  fPandoraBeam->Branch("secvertex",                     &fsecvertex,                    "secvertex[3]/D");
  fPandoraBeam->Branch("isprimarytrack",                &fisprimarytrack,               "isprimarytrack/I");
  fPandoraBeam->Branch("isprimaryshower",               &fisprimaryshower,              "isprimaryshower/I");
  fPandoraBeam->Branch("primaryBDTScore",               &fprimaryBDTScore,              "primaryBDTScore/D");
  fPandoraBeam->Branch("primaryNHits",                  &fprimaryNHits,                 "fprimaryNHits/I");
  fPandoraBeam->Branch("primaryTheta",                  &fprimaryTheta,                 "primaryTheta/D");
  fPandoraBeam->Branch("primaryPhi",                    &fprimaryPhi,                   "primaryPhi/D");
  fPandoraBeam->Branch("primaryLength",                 &fprimaryLength,                "primaryLength/D");
  fPandoraBeam->Branch("primaryMomentum",               &fprimaryMomentum,              "primaryMomentum/D");
  fPandoraBeam->Branch("primaryEndMomentum",            &fprimaryEndMomentum,           "primaryEndMomentum/D");
  fPandoraBeam->Branch("primaryEndPosition",            &fprimaryEndPosition,           "primaryEndPosition[3]/D");
  fPandoraBeam->Branch("primaryStartPosition",          &fprimaryStartPosition,         "primaryStartPosition[3]/D");
  fPandoraBeam->Branch("primaryEndDirection",           &fprimaryEndDirection,          "primaryEndDirection[3]/D");
  fPandoraBeam->Branch("primaryStartDirection",         &fprimaryStartDirection,        "primaryStartDirection[3]/D");
  fPandoraBeam->Branch("primaryOpeningAngle",           &fprimaryOpeningAngle,          "primaryOpeningAngle/D");
  fPandoraBeam->Branch("primaryID",                     &fprimaryID,                    "primaryID/I");
  fPandoraBeam->Branch("primaryShowerBestPlane",        &fprimaryShowerBestPlane,       "primaryShowerBestPlane/I");
  fPandoraBeam->Branch("primaryShowerEnergy",           &fprimaryShowerEnergy,          "primaryShowerEnergy/I");
  fPandoraBeam->Branch("primaryShowerMIPEnergy",        &fprimaryShowerMIPEnergy,       "primaryShowerMIPEnergy/I");
  fPandoraBeam->Branch("primaryShowerdEdx",             &fprimaryShowerdEdx,            "primaryShowerdEdx/I");
  fPandoraBeam->Branch("primaryMomentumByRangeProton",  &fprimaryMomentumByRangeProton, "primaryMomentumByRangeProton/D");
  fPandoraBeam->Branch("primaryMomentumByRangeMuon",    &fprimaryMomentumByRangeMuon,   "primaryMomentumByRangeMuon/D");
  fPandoraBeam->Branch("primaryKineticEnergy",          &fprimaryKineticEnergy,         "primaryKineticEnergy[3]/D");
  fPandoraBeam->Branch("primaryRange",                  &fprimaryRange,                 "primaryRange[3]/D");
  fPandoraBeam->Branch("primaryT0",                     &fprimaryT0,                    "primaryT0/D");

  fPandoraBeam->Branch("primary_truth_TrackId",         &fprimary_truth_TrackId,         "primary_truth_TrackId/I");
  fPandoraBeam->Branch("primary_truth_Pdg",             &fprimary_truth_Pdg,             "primary_truth_Pdg/I");
  fPandoraBeam->Branch("primary_truth_StartPosition",   &fprimary_truth_StartPosition,   "primary_truth_StartPosition[4]/D");
  fPandoraBeam->Branch("primary_truth_EndPosition",     &fprimary_truth_EndPosition,     "primary_truth_EndPosition[4]/D");
  fPandoraBeam->Branch("primary_truth_Momentum",        &fprimary_truth_Momentum,        "primary_truth_Momentum[4]/D");
  fPandoraBeam->Branch("primary_truth_EndMomentum",     &fprimary_truth_EndMomentum,     "primary_truth_EndMomentum[4]/D");
  fPandoraBeam->Branch("primary_truth_P",               &fprimary_truth_P,               "primary_truth_P/D");
  fPandoraBeam->Branch("primary_truth_Pt",              &fprimary_truth_Pt,              "primary_truth_Pt/D");
  fPandoraBeam->Branch("primary_truth_Mass",            &fprimary_truth_Mass,            "primary_truth_Mass/D");
  fPandoraBeam->Branch("primary_truth_Theta",           &fprimary_truth_Theta,           "primary_truth_Theta/D");
  fPandoraBeam->Branch("primary_truth_Phi",             &fprimary_truth_Phi,             "primary_truth_Phi/D");
  fPandoraBeam->Branch("primary_truth_Process",         &fprimary_truth_Process,         "primary_truth_Process/I");
  fPandoraBeam->Branch("primary_truth_Isbeammatched",   &fprimary_truth_Isbeammatched,   "primary_truth_Isbeammatched/I");
  fPandoraBeam->Branch("primary_truth_NDaughters",      &fprimary_truth_NDaughters,      "primary_truth_NDaughters/I");
  fPandoraBeam->Branch("fprimary_truth_EndProcess",     &fprimary_truth_EndProcess);

  fPandoraBeam->Branch("NDAUGHTERS",                    &fNDAUGHTERS,                   "NDAUGHTERS/I");
  fPandoraBeam->Branch("isdaughtertrack",               &fisdaughtertrack,              "isdaughtertrack[NDAUGHTERS]/I");
  fPandoraBeam->Branch("isdaughtershower",              &fisdaughtershower,             "isdaughtershower[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterNHits",                 &fdaughterNHits,                "daughterNHits[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterTheta",                 &fdaughterTheta,                "daughterTheta[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterPhi",                   &fdaughterPhi,                  "daughterPhi[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterLength",                &fdaughterLength,               "daughterLength[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterMomentum",              &fdaughterMomentum,             "daughterMomentum[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterEndMomentum",           &fdaughterEndMomentum,          "daughterEndMomentum[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterEndPosition",           &fdaughterEndPosition,          "daughterEndPosition[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterStartPosition",         &fdaughterStartPosition,        "daughterStartPosition[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterStartDirection",        &fdaughterStartDirection,       "daughterStartDirection[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterEndDirection",          &fdaughterEndDirection,         "daughterEndDirection[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterOpeningAngle",          &fdaughterOpeningAngle,         "daughterOpeningAngle[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerBestPlane",       &fdaughterShowerBestPlane,      "daughterShowerBestPlane[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerEnergy",          &fdaughterShowerEnergy,         "daughterShowerEnergy[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerMIPEnergy",       &fdaughterShowerMIPEnergy,      "daughterShowerMIPEnergy[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerdEdx",            &fdaughterShowerdEdx,           "daughterShowerdEdx[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterMomentumByRangeProton", &fdaughterMomentumByRangeProton,"daughterMomentumByRangeProton[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterMomentumByRangeMuon",   &fdaughterMomentumByRangeMuon,  "daughterMomentumByRangeMuon[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterKineticEnergy",         &fdaughterKineticEnergy,        "daughterKineticEnergy[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterRange",                 &fdaughterRange,                "daughterRange[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterID",                    &fdaughterID,                   "daughterID[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterT0",                    &fdaughterT0,                   "daughterT0[NDAUGHTERS]/D");

  fPandoraBeam->Branch("daughter_truth_TrackId",        &fdaughter_truth_TrackId,       "daughter_truth_TrackId[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughter_truth_Pdg",            &fdaughter_truth_Pdg,           "daughter_truth_Pdg[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughter_truth_StartPosition",  &fdaughter_truth_StartPosition, "daughter_truth_StartPosition[NDAUGHTERS][4]/D");
  fPandoraBeam->Branch("daughter_truth_EndPosition",    &fdaughter_truth_EndPosition,   "daughter_truth_EndPosition[NDAUGHTERS][4]/D");
  fPandoraBeam->Branch("daughter_truth_Momentum",       &fdaughter_truth_Momentum,      "daughter_truth_Momentum[NDAUGHTERS][4]/D");
  fPandoraBeam->Branch("daughter_truth_EndMomentum",    &fdaughter_truth_EndMomentum,   "daughter_truth_EndMomentum[NDAUGHTERS][4]/D");
  fPandoraBeam->Branch("daughter_truth_P",              &fdaughter_truth_P,             "daughter_truth_P[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Pt",             &fdaughter_truth_Pt,            "daughter_truth_Pt[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Mass",           &fdaughter_truth_Mass,          "daughter_truth_Mass[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Theta",          &fdaughter_truth_Theta,         "daughter_truth_Theta[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Phi",            &fdaughter_truth_Phi,           "daughter_truth_Phi[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Process",        &fdaughter_truth_Process,       "daughter_truth_Process[NDAUGHTERS]/I");
 

  fPandoraBeam->Branch("primtrk_dqdx",&primtrk_dqdx);
  fPandoraBeam->Branch("primtrk_dedx",&primtrk_dedx);
  fPandoraBeam->Branch("primtrk_resrange",&primtrk_resrange);
  fPandoraBeam->Branch("primtrk_range",&primtrk_range);
  fPandoraBeam->Branch("primtrk_hitx",&primtrk_hitx);
  fPandoraBeam->Branch("primtrk_hity",&primtrk_hity);
  fPandoraBeam->Branch("primtrk_hitz",&primtrk_hitz);
  fPandoraBeam->Branch("primtrk_pitch",&primtrk_pitch);


  ///////////////////////////MC truth  from BackTracker///////
  /* fPandoraBeam->Branch("truth_startx",&truth_startx);
  fPandoraBeam->Branch("truth_starty",&truth_starty);
  fPandoraBeam->Branch("truth_startz",&truth_startz);
  fPandoraBeam->Branch("truth_endx",&truth_endx);
  fPandoraBeam->Branch("truth_endy",&truth_endy);
  fPandoraBeam->Branch("truth_endz",&truth_endz);
  fPandoraBeam->Branch("truth_Ndaughter",&truth_Ndaughter);
  fPandoraBeam->Branch("truth_pdg",&truth_pdg);
  fPandoraBeam->Branch("truth_process",&truth_process);*/
  ////////////////////////////////////////////////////////////

}

void protoana::pionanalysismc::analyze(art::Event const & evt){

  // Initialise tree parameters
  Initialise();

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;



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
  if(!evt.isRealData()){
    for(int k=0; k < 6; k++)
      fNactivefembs[k] = 20;
  }
  else{
    for(int k=0; k < 6; k++)
     fNactivefembs[k] = dataUtil.GetNActiveFembsForAPA(evt, k);
  }

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

  bool beamTriggerEvent = false;
  // If this event is MC then we can check what the true beam particle is
  if(!evt.isRealData()){
    // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
    // simulation has fGeneratorTag = "generator"
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);

    // Also get the reconstructed beam information in the MC - TO DO
    //auto beamsim = evt.getValidHandle<std::vector<sim::ProtoDUNEbeamsim> >(fGeneratorTag);
    //const sim::ProtoDUNEbeamsim beamsimobj = (*beamsim)[0];
    //std::cout << beamsimobj.NInstruments() << std::endl;
    //sim::ProtoDUNEBeamInstrument beamsim_tof1 = beamsimobj.GetInstrument("TOF1");
    //sim::ProtoDUNEBeamInstrument beamsim_trig2 = beamsimobj.GetInstrument("TRIG2");
    //std::cout << beamsim_trig2.GetT() - beamsim_tof1.GetT() << " , " << beamsim_trig2.GetSmearedVar1() - beamsim_tof1.GetSmearedVar1() << std::endl;
    //std::cout << beamsimobj.GetInstrument("TRIG2").GetT() - beamsimobj.GetInstrument("TOF1").GetT() << " , " << beamsimobj.GetInstrument("TRIG2").GetSmearedVar1() - beamsimobj.GetInstrument("TOF1").GetSmearedVar1() << std::endl;
  
    // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
    // of these, so we pass the first element into the function to get the good particle
    const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);

    if(geantGoodParticle != 0x0){
      std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() 
		<< " , track id = " << geantGoodParticle->TrackId()
		<< " , Vx/Vy/Vz = " << geantGoodParticle->Vx() << "/"<< geantGoodParticle->Vy() << "/" << geantGoodParticle->Vz() 
		<< std::endl;

      beamTriggerEvent = true;
      fbeamtrigger       = 12;
      fbeamtrackPos[0]   = geantGoodParticle->Vx();
      fbeamtrackPos[1]   = geantGoodParticle->Vy();
      fbeamtrackPos[2]   = geantGoodParticle->Vz();
      fbeamtrackMomentum = geantGoodParticle->P();
      fbeamtrackP[0]     = geantGoodParticle->Px();
      fbeamtrackP[1]     = geantGoodParticle->Py();
      fbeamtrackP[2]     = geantGoodParticle->Pz();
      fbeamtrackEnergy   = geantGoodParticle->E();
      fbeamtrackPdg      = geantGoodParticle->PdgCode();
      fbeamtrackTime     = geantGoodParticle->T();
      fbeamtrackID       = geantGoodParticle->TrackId();
    }
  }
  else{
    // For data we can see if this event comes from a beam trigger
    beamTriggerEvent = dataUtil.IsBeamTrigger(evt);

    art::Handle< std::vector<beam::ProtoDUNEBeamEvent> > pdbeamHandle;
    std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;
    if(evt.getByLabel(fBeamModuleLabel, pdbeamHandle))
      art::fill_ptr_vector(beaminfo, pdbeamHandle);
  
    for(unsigned int i = 0; i < beaminfo.size(); ++i){
      //if(!beaminfo[i]->CheckIsMatched()) continue;
      fbeamtrigger = beaminfo[i]->GetTimingTrigger();
      fbeamtrackTime = (double)beaminfo[i]->GetRDTimestamp();

      // If ToF is 0-3 there was a good match corresponding to the different pair-wise combinations of the upstream and downstream channels
      if(beaminfo[i]->GetTOFChan() >= 0)
	ftof =  beaminfo[i]->GetTOF();

      // Get Cerenkov
      if(beaminfo[i]->GetBITrigger() == 1){
	fcerenkov1 = beaminfo[i]->GetCKov0Status();
	fcerenkov2 = beaminfo[i]->GetCKov1Status();
      }

      // Beam particle could have more than one tracks - for now take the first one, need to do this properly
      auto & tracks = beaminfo[i]->GetBeamTracks();
      if(!tracks.empty()){
	fbeamtrackPos[0] = tracks[0].End().X();
	fbeamtrackPos[1] = tracks[0].End().Y();
	fbeamtrackPos[2] = tracks[0].End().Z();
	fbeamtrackDir[0] = tracks[0].StartDirection().X();
	fbeamtrackDir[1] = tracks[0].StartDirection().Y();
	fbeamtrackDir[2] = tracks[0].StartDirection().Z();
      }
  
      // Beam momentum
      auto & beammom = beaminfo[i]->GetRecoBeamMomenta();
      if(!beammom.empty())
	fbeamtrackMomentum = beammom[0];

      // For now only take the first beam particle - need to add some criteria if more than one are found
      break;
 
    }
  }

  /*
  // Now we want to access the output from Pandora. This comes in the form of particle flow objects (recob::PFParticle).
  // The primary PFParticles are those we want to consider and these PFParticles then have a hierarchy of daughters that
  // describe the whole interaction of a given primary particle
  //
  //                     / daughter track
  //                    /
  //  primary track    /   
  //  ---------------- ---- daughter track
  //                   \
  //                   /\-
  //                   /\\-- daughter shower
  //
  // The above primary PFParticle will have links to three daughter particles, two track-like and one shower-like
  */

  // Track momentum algorithm calculates momentum based on track range
  trkf::TrackMomentumCalculator trmom;
  //trmom.SetMinLength(100);

  // Get all of the PFParticles, by default from the "pandora" product
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  //std::cout << "All primary pfParticles = " <<  pfpUtil.GetNumberPrimaryPFParticle(evt,fPFParticleTag) << std::endl;

  // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility 
  // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
  // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
  // PFParticles considered as primary
  // std::vector<recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  auto pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);

  // We can now look at these particles
  for(const recob::PFParticle* particle : pfParticles){

    // Pandora's BDT beam-cosmic score
    fprimaryBDTScore = (double)pfpUtil.GetBeamCosmicScore(*particle,evt,fPFParticleTag);
    
    // NHits associated with this pfParticle
    fprimaryNHits = (pfpUtil.GetPFParticleHits(*particle,evt,fPFParticleTag)).size();

    // Get the T0 for this pfParticle
    std::vector<anab::T0> pfT0vec = pfpUtil.GetPFParticleT0(*particle,evt,fPFParticleTag);
    if(!pfT0vec.empty())
      fprimaryT0 = pfT0vec[0].Time();
 
    //std::cout << "Pdg Code = " << particle->PdgCode() << std::endl;
    // "particle" is the pointer to our beam particle. The recob::Track or recob::Shower object
    // of this particle might be more helpful. These return null pointers if not track-like / shower-like
    const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
    const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
    if(thisTrack != 0x0){
      fisprimarytrack               = 1;
      fisprimaryshower              = 0;

      fprimaryID                    = thisTrack->ParticleId();
      fprimaryTheta                 = thisTrack->Theta();
      fprimaryPhi                   = thisTrack->Phi();
      fprimaryLength                = thisTrack->Length();
      fprimaryMomentum              = thisTrack->StartMomentum();
      fprimaryEndMomentum           = thisTrack->EndMomentum();

      fprimaryEndPosition[0]        = thisTrack->End().X();
      fprimaryEndPosition[1]        = thisTrack->End().Y();
      fprimaryEndPosition[2]        = thisTrack->End().Z();
      fprimaryStartPosition[0]      = thisTrack->Start().X();
      fprimaryStartPosition[1]      = thisTrack->Start().Y();
      fprimaryStartPosition[2]      = thisTrack->Start().Z();
      fprimaryEndDirection[0]       = thisTrack->Trajectory().EndDirection().X();
      fprimaryEndDirection[1]       = thisTrack->Trajectory().EndDirection().Y();
      fprimaryEndDirection[2]       = thisTrack->Trajectory().EndDirection().Z();
      fprimaryStartDirection[0]     = thisTrack->Trajectory().StartDirection().X();
      fprimaryStartDirection[1]     = thisTrack->Trajectory().StartDirection().Y();
      fprimaryStartDirection[2]     = thisTrack->Trajectory().StartDirection().Z();

      fprimaryMomentumByRangeMuon   = trmom.GetTrackMomentum(thisTrack->Length(),13);
      fprimaryMomentumByRangeProton = trmom.GetTrackMomentum(thisTrack->Length(),2212);      
      
      std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
      if(calovector.size() != 3)
	std::cerr << "WARNING::Calorimetry vector size for primary is = " << calovector.size() << std::endl;

      //HY::Get the Calorimetry(s) from thisTrack
      //std::vector<double> tmp_primtrk_dqdx;	
      //std::vector<double> tmp_primtrk_resrange;	
      //std::vector<double> tmp_primtrk_dedx;	
      std::vector<double> tmp_primtrk_dqdx;	
      std::vector<double> tmp_primtrk_resrange;	
      std::vector<double> tmp_primtrk_dedx;	
      std::vector<double> tmp_primtrk_hitx;	
      std::vector<double> tmp_primtrk_hity;	
      std::vector<double> tmp_primtrk_hitz;
      std::vector<double> tmp_primtrk_pitch;
      for (auto & calo : calovector) {
      	if (calo.PlaneID().Plane == 2){ //only collection plane
		primtrk_range.push_back(calo.Range());
        	for (size_t ihit = 0; ihit < calo.dQdx().size(); ++ihit){ //loop over hits
			tmp_primtrk_dqdx.push_back(calo.dQdx()[ihit]);
		      	tmp_primtrk_resrange.push_back(calo.ResidualRange()[ihit]);
			tmp_primtrk_dedx.push_back(calo.dEdx()[ihit]);
			tmp_primtrk_pitch.push_back(calo.TrkPitchVec()[ihit]);

			const auto &primtrk_pos=(calo.XYZ())[ihit];
		      	tmp_primtrk_hitx.push_back(primtrk_pos.X());
		      	tmp_primtrk_hity.push_back(primtrk_pos.Y());
		      	tmp_primtrk_hitz.push_back(primtrk_pos.Z());
		      	std::cout<<"dqdx="<<calo.dQdx()[ihit]<<"; resrange="<<calo.ResidualRange()[ihit]<<std::endl;
		      	std::cout<<"(X,Y,Z)="<<"("<<calo.XYZ()[ihit].X()<<","<<calo.XYZ()[ihit].Y()<<","<<calo.XYZ()[ihit].Z()<<")"<<std::endl;

                 } //loop over hits
        } //only collection plane
      }

      primtrk_dqdx.push_back(tmp_primtrk_dqdx);
      primtrk_resrange.push_back(tmp_primtrk_resrange);
      primtrk_dedx.push_back(tmp_primtrk_dedx);
      primtrk_hitx.push_back(tmp_primtrk_hitx);
      primtrk_hity.push_back(tmp_primtrk_hity);
      primtrk_hitz.push_back(tmp_primtrk_hitz);
      primtrk_pitch.push_back(tmp_primtrk_pitch);

      tmp_primtrk_dqdx.clear();
      tmp_primtrk_resrange.clear();
      tmp_primtrk_dedx.clear();
      tmp_primtrk_hitx.clear();
      tmp_primtrk_hity.clear();
      tmp_primtrk_hitz.clear();
      tmp_primtrk_pitch.clear();

      for(size_t k = 0; k < calovector.size() && k<3; k++){
	fprimaryKineticEnergy[k] = calovector[k].KineticEnergy();
	fprimaryRange[k] = calovector[k].Range();
	//const std::vector< double > & dedxvec = calovector[k].dEdx();
      }

      //Get MC truth info using BackTracker





      // Get the true mc particle
      const simb::MCParticle* mcparticle1 = truthUtil.GetMCParticleFromRecoTrack(clockData, *thisTrack, evt, fTrackerTag);
      if(mcparticle1 != 0x0){
        const simb::MCParticle *mcparticle=pi_serv->TrackIdToMotherParticle_P(mcparticle1->TrackId());
        if (mcparticle) {





	fprimary_truth_TrackId          = mcparticle->TrackId();
	fprimary_truth_Pdg              = mcparticle->PdgCode();
	fprimary_truth_StartPosition[0] = mcparticle->Vx();
	fprimary_truth_StartPosition[1] = mcparticle->Vy();
	fprimary_truth_StartPosition[2] = mcparticle->Vz();
	fprimary_truth_StartPosition[3] = mcparticle->T();
	fprimary_truth_EndPosition[0]   = mcparticle->EndX();
	fprimary_truth_EndPosition[1]   = mcparticle->EndY();
	fprimary_truth_EndPosition[2]   = mcparticle->EndZ();
	fprimary_truth_EndPosition[3]   = mcparticle->EndT();
	fprimary_truth_P                = mcparticle->P();
	fprimary_truth_Momentum[0]      = mcparticle->Px();
	fprimary_truth_Momentum[1]      = mcparticle->Py();
	fprimary_truth_Momentum[2]      = mcparticle->Pz();
	fprimary_truth_Momentum[3]      = mcparticle->E();
	fprimary_truth_Pt               = mcparticle->Pt();
	fprimary_truth_Mass             = mcparticle->Mass();
	fprimary_truth_EndMomentum[0]   = mcparticle->EndPx();
	fprimary_truth_EndMomentum[1]   = mcparticle->EndPy();
	fprimary_truth_EndMomentum[2]   = mcparticle->EndPz();
	fprimary_truth_EndMomentum[3]   = mcparticle->EndE();
	fprimary_truth_Theta            = mcparticle->Momentum().Theta();
	fprimary_truth_Phi              = mcparticle->Momentum().Phi();
	fprimary_truth_NDaughters       = mcparticle->NumberDaughters();
	fprimary_truth_EndProcess       = mcparticle->EndProcess();
	fprimary_truth_Process          = int(mcparticle->Trajectory().ProcessToKey(mcparticle->Process()));
	if(fbeamtrackID != -999 && fprimary_truth_TrackId == fbeamtrackID)
	  fprimary_truth_Isbeammatched    = 1;
	else
	  fprimary_truth_Isbeammatched    = 0;

	//std::cout << "Process = " << (mcparticle->Process()).c_str() << " , End process = " << (mcparticle->EndProcess()).c_str()
	// << " , track ID = " << mcparticle->TrackId()
	//	  << std::endl;
	std::cout<<"End Process name "<<mcparticle->EndProcess()<<std::endl;
        }
      }
    }
    else if(thisShower != 0x0){
      fisprimarytrack               = 0;
      fisprimaryshower              = 1;

      fprimaryID                    = thisShower->ID();
      fprimaryLength                = thisShower->Length();
      fprimaryShowerBestPlane       = thisShower->best_plane();
      fprimaryOpeningAngle          = thisShower->OpenAngle();
      fprimaryStartPosition[0]      = thisShower->ShowerStart().X();
      fprimaryStartPosition[1]      = thisShower->ShowerStart().Y();
      fprimaryStartPosition[2]      = thisShower->ShowerStart().Z();
      fprimaryStartDirection[0]     = thisShower->Direction().X();
      fprimaryStartDirection[1]     = thisShower->Direction().Y();
      fprimaryStartDirection[2]     = thisShower->Direction().Z();
      if( (thisShower->Energy()).size() > 0 )
	fprimaryShowerEnergy = thisShower->Energy()[0]; // thisShower->best_plane()
      if( (thisShower->MIPEnergy()).size() > 0 )
	fprimaryShowerMIPEnergy = thisShower->MIPEnergy()[0];
      if( (thisShower->dEdx()).size() > 0 )
	fprimaryShowerdEdx = thisShower->dEdx()[0];
    }
    else{
      std::cout << "INFO::Primary pfParticle is not track or shower. Skip!" << std::endl;
      continue;
    }
    
    // Find the particle vertex. We need the tracker tag here because we need to do a bit of
    // additional work if the PFParticle is track-like to find the vertex. 
    const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    fvertex[0] = vtx.X(); fvertex[1] = vtx.Y(); fvertex[2] = vtx.Z();

    // Now we can look for the interaction point of the particle if one exists, i.e where the particle
    // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
    // check this by making sure we get a sensible position
    const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    fsecvertex[0] = interactionVtx.X(); fsecvertex[1] = interactionVtx.Y(); fsecvertex[2] = interactionVtx.Z();

    // Maximum number of daugthers to be processed
    if(particle->NumDaughters() > NMAXDAUGTHERS)
      std::cout << "INFO::Number of daughters is " << particle->NumDaughters() << ". Only the first NMAXDAUGTHERS are processed." << std::endl;

    // Let's get the daughter PFParticles... we can do this simply without the utility
    for(const int daughterID : particle->Daughters()){
      // Daughter ID is the element of the original recoParticle vector
      const recob::PFParticle *daughterParticle      = &(recoParticles->at(daughterID));
      std::cout << "Daughter " << daughterID << " has " << daughterParticle->NumDaughters() << " daughters" << std::endl;
      
      const recob::Track* daughterTrack              = pfpUtil.GetPFParticleTrack(*daughterParticle,evt,fPFParticleTag,fTrackerTag);
      const recob::Shower* daughterShower            = pfpUtil.GetPFParticleShower(*daughterParticle,evt,fPFParticleTag,fShowerTag);
  
      if(daughterTrack != 0x0){
	fisdaughtertrack[fNDAUGHTERS]                = 1;
	fisdaughtershower[fNDAUGHTERS]               = 0;
	fdaughterTheta[fNDAUGHTERS]                  = daughterTrack->Theta();
	fdaughterPhi[fNDAUGHTERS]                    = daughterTrack->Phi();
	fdaughterLength[fNDAUGHTERS]                 = daughterTrack->Length();
	fdaughterMomentum[fNDAUGHTERS]               = daughterTrack->StartMomentum();
	fdaughterEndMomentum[fNDAUGHTERS]            = daughterTrack->EndMomentum();
	fdaughterStartPosition[fNDAUGHTERS][0]       = daughterTrack->Trajectory().Start().X();
	fdaughterStartPosition[fNDAUGHTERS][1]       = daughterTrack->Trajectory().Start().Y();
	fdaughterStartPosition[fNDAUGHTERS][2]       = daughterTrack->Trajectory().Start().Z();
	fdaughterEndPosition[fNDAUGHTERS][0]         = daughterTrack->Trajectory().End().X();
	fdaughterEndPosition[fNDAUGHTERS][1]         = daughterTrack->Trajectory().End().Y();
	fdaughterEndPosition[fNDAUGHTERS][2]         = daughterTrack->Trajectory().End().Z();
	fdaughterStartDirection[fNDAUGHTERS][0]      = daughterTrack->Trajectory().StartDirection().X();
	fdaughterStartDirection[fNDAUGHTERS][1]      = daughterTrack->Trajectory().StartDirection().Y();
	fdaughterStartDirection[fNDAUGHTERS][2]      = daughterTrack->Trajectory().StartDirection().Z();
	fdaughterEndDirection[fNDAUGHTERS][0]        = daughterTrack->Trajectory().EndDirection().X();
	fdaughterEndDirection[fNDAUGHTERS][1]        = daughterTrack->Trajectory().EndDirection().Y();
	fdaughterEndDirection[fNDAUGHTERS][2]        = daughterTrack->Trajectory().EndDirection().Z();

	fdaughterMomentumByRangeMuon[fNDAUGHTERS]    = trmom.GetTrackMomentum(daughterTrack->Length(),13);
	fdaughterMomentumByRangeProton[fNDAUGHTERS]  = trmom.GetTrackMomentum(daughterTrack->Length(),2212);

	std::vector<anab::Calorimetry> daughtercalovector = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackerTag, fCalorimetryTag);
	if(daughtercalovector.size() != 3)
	  std::cerr << "WARNING::Calorimetry vector size for daughter is = " << daughtercalovector.size() << std::endl;

	for(size_t k = 0; k < daughtercalovector.size() && k<3; k++){
	  fdaughterKineticEnergy[fNDAUGHTERS][k] = daughtercalovector[k].KineticEnergy();
	  fdaughterRange[fNDAUGHTERS][k] = daughtercalovector[k].Range();
	}

	// Get the true mc particle
        const simb::MCParticle* mcdaughterparticle = truthUtil.GetMCParticleFromRecoTrack(clockData, *daughterTrack, evt, fTrackerTag);
	if(mcdaughterparticle != 0x0){
	  fdaughter_truth_TrackId[fNDAUGHTERS]          = mcdaughterparticle->TrackId();
	  fdaughter_truth_Pdg[fNDAUGHTERS]              = mcdaughterparticle->PdgCode();
	  fdaughter_truth_StartPosition[fNDAUGHTERS][0] = mcdaughterparticle->Vx();
	  fdaughter_truth_StartPosition[fNDAUGHTERS][1] = mcdaughterparticle->Vy();
	  fdaughter_truth_StartPosition[fNDAUGHTERS][2] = mcdaughterparticle->Vz();
	  fdaughter_truth_StartPosition[fNDAUGHTERS][3] = mcdaughterparticle->T();
	  fdaughter_truth_EndPosition[fNDAUGHTERS][0]   = mcdaughterparticle->EndX();
	  fdaughter_truth_EndPosition[fNDAUGHTERS][1]   = mcdaughterparticle->EndY();
	  fdaughter_truth_EndPosition[fNDAUGHTERS][2]   = mcdaughterparticle->EndZ();
	  fdaughter_truth_EndPosition[fNDAUGHTERS][3]   = mcdaughterparticle->EndT();
	  fdaughter_truth_P[fNDAUGHTERS]                = mcdaughterparticle->P();
	  fdaughter_truth_Momentum[fNDAUGHTERS][0]      = mcdaughterparticle->Px();
	  fdaughter_truth_Momentum[fNDAUGHTERS][1]      = mcdaughterparticle->Py();
	  fdaughter_truth_Momentum[fNDAUGHTERS][2]      = mcdaughterparticle->Pz();
	  fdaughter_truth_Momentum[fNDAUGHTERS][3]      = mcdaughterparticle->E();
	  fdaughter_truth_Pt[fNDAUGHTERS]               = mcdaughterparticle->Pt();
	  fdaughter_truth_Mass[fNDAUGHTERS]             = mcdaughterparticle->Mass();
	  fdaughter_truth_EndMomentum[fNDAUGHTERS][0]   = mcdaughterparticle->EndPx();
	  fdaughter_truth_EndMomentum[fNDAUGHTERS][1]   = mcdaughterparticle->EndPy();
	  fdaughter_truth_EndMomentum[fNDAUGHTERS][2]   = mcdaughterparticle->EndPz();
	  fdaughter_truth_EndMomentum[fNDAUGHTERS][3]   = mcdaughterparticle->EndE();
	  fdaughter_truth_Theta[fNDAUGHTERS]            = mcdaughterparticle->Momentum().Theta();
	  fdaughter_truth_Phi[fNDAUGHTERS]              = mcdaughterparticle->Momentum().Phi();
	  fdaughter_truth_Process[fNDAUGHTERS]          = int(mcdaughterparticle->Trajectory().ProcessToKey(mcdaughterparticle->Process()));
	  std::cout << "Daughter Process = " << (mcdaughterparticle->Process()).c_str() 
		    << " , mother = " << mcdaughterparticle->Mother() 
		    << std::endl;
	}
      }
      else if(daughterShower != 0x0){
	fisdaughtertrack[fNDAUGHTERS]                = 0;
	fisdaughtershower[fNDAUGHTERS]               = 1;
	fdaughterLength[fNDAUGHTERS]                 = daughterShower->Length();
	fdaughterShowerBestPlane[fNDAUGHTERS]        = daughterShower->best_plane();
	fdaughterOpeningAngle[fNDAUGHTERS]           = daughterShower->OpenAngle();
	fdaughterStartPosition[fNDAUGHTERS][0]       = daughterShower->ShowerStart().X();
	fdaughterStartPosition[fNDAUGHTERS][1]       = daughterShower->ShowerStart().Y();
	fdaughterStartPosition[fNDAUGHTERS][2]       = daughterShower->ShowerStart().Z();
	fdaughterStartDirection[fNDAUGHTERS][0]      = daughterShower->Direction().X();
	fdaughterStartDirection[fNDAUGHTERS][1]      = daughterShower->Direction().Y();
	fdaughterStartDirection[fNDAUGHTERS][2]      = daughterShower->Direction().Z();
	if( (daughterShower->Energy()).size() > 0 )
	  fdaughterShowerEnergy[fNDAUGHTERS] = daughterShower->Energy()[0]; // thisShower->best_plane()
	if( (daughterShower->MIPEnergy()).size() > 0 )
	  fdaughterShowerMIPEnergy[fNDAUGHTERS] = daughterShower->MIPEnergy()[0];
	if( (daughterShower->dEdx()).size() > 0 )
	  fdaughterShowerdEdx[fNDAUGHTERS] = daughterShower->dEdx()[0];
      }
      else{
	std::cout << "INFO::Daughter pfParticle is not track or shower. Skip!" << std::endl;
	continue;
      }

      fdaughterID[fNDAUGHTERS]                       = daughterID;
      // NHits associated with this pfParticle
      fdaughterNHits[fNDAUGHTERS]                    = (pfpUtil.GetPFParticleHits(*daughterParticle,evt,fPFParticleTag)).size();
      // T0
      std::vector<anab::T0> pfdaughterT0vec = pfpUtil.GetPFParticleT0(*daughterParticle,evt,fPFParticleTag);
      if(!pfT0vec.empty())
	fdaughterT0[fNDAUGHTERS] = pfdaughterT0vec[0].Time();

      fNDAUGHTERS++;

      // Only process NMAXDAUGTHERS
      if(fNDAUGHTERS > NMAXDAUGTHERS) break;

    }
 
    // For actually studying the objects, it is easier to have the daughters in their track and shower forms.
    // We can use the utility to get a vector of track-like and a vector of shower-like daughters
    //const std::vector<const recob::Track*> trackDaughters = pfpUtil.GetPFParticleDaughterTracks(*particle,evt,fPFParticleTag,fTrackerTag);  
    //const std::vector<const recob::Shower*> showerDaughters = pfpUtil.GetPFParticleDaughterShowers(*particle,evt,fPFParticleTag,fShowerTag);
 
    // For now only consider the first primary track. Need a proper treatment if more than one primary particles are found.
    break;
  } 

  // Fill trees
  if(beamTriggerEvent)
    fPandoraBeam->Fill();

  //fPandoraCosmics->Fill();

}

void protoana::pionanalysismc::endJob(){

}

void protoana::pionanalysismc::FillCosmicsTree(art::Event const & evt, std::string pfParticleTag){

  // To fill

}

void protoana::pionanalysismc::Initialise(){
  
  fRun = -999;
  fSubRun = -999;
  fevent = -999;
  fTimeStamp = -999.0;
  for(int k=0; k < 5; k++)
    fNactivefembs[k] = -999;

  for(int k=0; k < 3; k++){
    fvertex[k] = -999.0;
    fsecvertex[k] = -999.0;
    fprimaryEndPosition[k] = -999.0;
    fprimaryStartPosition[k] = -999.0;
    fprimaryEndDirection[k] = -999.0;
    fprimaryStartDirection[k] = -999.0;
    fprimaryKineticEnergy[k] = -999.0;
    fprimaryRange[k] = -999.0;
  }

  fbeamtrigger = -999;
  ftof = -999.0;
  fcerenkov1 = -999;
  fcerenkov2 = -999;
  fbeamtrackMomentum = -999.0;
  fbeamtrackEnergy = 999.0;
  fbeamtrackPdg = -999;
  fbeamtrackTime = -999.0;
  fbeamtrackID = -999;
  for(int l=0; l < 3; l++){
    fbeamtrackP[l] = -999.0;
    fbeamtrackPos[l] = -999.0;
    fbeamtrackDir[l] = -999.0;
  }
 
  fisprimarytrack = 0;
  fisprimaryshower = 0;
  fNDAUGHTERS = 0;

  fprimaryBDTScore = -999.0;
  fprimaryNHits = -999;
  fprimaryTheta = -999.0;
  fprimaryPhi = -999.0;
  fprimaryLength = -999.0;
  fprimaryMomentum = -999.0;
  fprimaryEndMomentum = -999.0;
  fprimaryOpeningAngle = -999.0;
  fprimaryShowerBestPlane = -999;
  fprimaryShowerEnergy = -999;
  fprimaryShowerMIPEnergy = -999;
  fprimaryShowerdEdx = -999;
  fprimaryID = -999;
  fprimaryMomentumByRangeProton = -999.0;
  fprimaryMomentumByRangeMuon = -999.0;
  fprimaryT0 = -999.0;

  fprimary_truth_TrackId = -999;
  fprimary_truth_Pdg = -999;
  fprimary_truth_P = -999.0;
  fprimary_truth_Pt = -999.0;
  fprimary_truth_Mass = -999.0;
  fprimary_truth_Theta = -999.0;
  fprimary_truth_Phi = -999.0;
  fprimary_truth_Process = -999;
  fprimary_truth_Isbeammatched = -999;
  fprimary_truth_NDaughters = -999;
  for(int l=0; l < 4; l++){
    fprimary_truth_StartPosition[l] = -999.0;
    fprimary_truth_EndPosition[l] = -999.0;
    fprimary_truth_Momentum[l] = -999.0;
    fprimary_truth_EndMomentum[l] = -999.0;
  }

  for(int k=0; k < NMAXDAUGTHERS; k++){
    fisdaughtertrack[k] = -999;
    fisdaughtershower[k] = -999;
    fdaughterNHits[k] = -999;
    fdaughterTheta[k] = -999.0;
    fdaughterPhi[k] = -999.0;
    fdaughterLength[k] = -999.0;
    fdaughterMomentum[k] = -999.0;
    fdaughterEndMomentum[k] = -999.0;
    for(int l=0; l < 3; l++){
      fdaughterEndPosition[k][l] = -999.0;
      fdaughterStartPosition[k][l] = -999.0;
      fdaughterEndDirection[k][l] = -999.0;
      fdaughterStartDirection[k][l] = -999.0;
      fdaughterKineticEnergy[k][l] = -999.0;
      fdaughterRange[k][l] = -999.0;
    }
    fdaughterOpeningAngle[k] = -999.0;
    fdaughterShowerBestPlane[k] = -999;
    fdaughterShowerEnergy[k] = -999;
    fdaughterShowerMIPEnergy[k] = -999;
    fdaughterShowerdEdx[k] = -999;
    fdaughterMomentumByRangeProton[k] = -999.0;
    fdaughterMomentumByRangeMuon[k] = -999.0;
    fdaughterID[k] = -999;
    fdaughterT0[k] = -999;

    fdaughter_truth_TrackId[k] = -999;
    fdaughter_truth_Pdg[k] = -999;
    fdaughter_truth_P[k] = -999.0;
    fdaughter_truth_Pt[k] = -999.0;
    fdaughter_truth_Mass[k] = -999.0;
    fdaughter_truth_Theta[k] = -999.0;
    fdaughter_truth_Phi[k] = -999.0;
    fdaughter_truth_Process[k] = -999;
    for(int l=0; l < 4; l++){
      fdaughter_truth_StartPosition[k][l] = -999.0;
      fdaughter_truth_EndPosition[k][l] = -999.0;
      fdaughter_truth_Momentum[k][l] = -999.0;
      fdaughter_truth_EndMomentum[k][l] = -999.0;
    }
  }

  primtrk_dqdx.clear();
  primtrk_resrange.clear();
  primtrk_dedx.clear();
  primtrk_range.clear();
  primtrk_hitx.clear();
  primtrk_hity.clear();
  primtrk_hitz.clear();
  primtrk_pitch.clear();


}

DEFINE_ART_MODULE(protoana::pionanalysismc)
