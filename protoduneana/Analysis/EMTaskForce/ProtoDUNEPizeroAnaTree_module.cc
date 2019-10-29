////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEPizeroAnaTree
// File:        ProtoDUNEPizeroAnaTree_module.cc
//
// Extract protoDUNE useful information, do a first pre-selection and
// save output to a flat tree
//
//
// Georgios Christodoulou - georgios.christodoulou at cern.ch
// modified by Aaron Higuera ahiguera@central.uh.edu
// and Milo Vermeulen milov@nikhef.nl
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"
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
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcore/Geometry/Geometry.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "dune/Calib/XYZCalib.h"
#include "dune/CalibServices/XYZCalibService.h"
#include "dune/CalibServices/XYZCalibServiceProtoDUNE.h"
#include "lardata/ArtDataHelper/MVAReader.h"


#include "protoduneana/Analysis/ProtoDUNETrackUtils.h"
#include "protoduneana/Analysis/ProtoDUNETruthUtils.h"
#include "protoduneana/Analysis/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Analysis/ProtoDUNEDataUtils.h"
#include "protoduneana/Analysis/ProtoDUNEShowerUtils.h"
#include "protoduneana/Analysis/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Analysis/EMTaskForce/PiZeroProcess.h"

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
//const int NMAXDAUGTHERS = 25;
// Maximum number of hits to save
const int NMAXHITS = 5000;

namespace protoana {
  class ProtoDUNEPizeroAnaTree;
}


class protoana::ProtoDUNEPizeroAnaTree : public art::EDAnalyzer {
public:

  explicit ProtoDUNEPizeroAnaTree(fhicl::ParameterSet const & p);

  ProtoDUNEPizeroAnaTree(ProtoDUNEPizeroAnaTree const &) = delete;
  ProtoDUNEPizeroAnaTree(ProtoDUNEPizeroAnaTree &&) = delete;
  ProtoDUNEPizeroAnaTree & operator = (ProtoDUNEPizeroAnaTree const &) = delete;
  ProtoDUNEPizeroAnaTree & operator = (ProtoDUNEPizeroAnaTree &&) = delete;

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
  protoana::ProtoDUNEShowerUtils showerUtil;
  protoana::ProtoDUNEBeamlineUtils beamlineUtil;


  // Track momentum algorithm calculates momentum based on track range
  const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  // Initialise tree variables
  void Initialise();
  void ResetPi0Vars();
  void branchPi0TTree(TTree* tree);
  void FillPrimaryPFParticle(art::Event const & evt, const recob::PFParticle* particle);
  void setPiZeroInfo(art::Event const & evt, const pizero::PiZeroProcess& pzproc);

  // Fill cosmics tree

  // fcl parameters
  const art::InputTag fBeamModuleLabel;
  std::string fSimulationTag;
  std::string fCalorimetryTag;
  calo::CalorimetryAlg fCalorimetryAlg;
  std::string fParticleIDTag;
  std::string fTrackTag;
  std::string fShowerTag;
  std::string fShowerCaloTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  std::string fHitTag;
  int fVerbose;

  TTree *fPandoraBeam;
  TTree *fMCPi0;
  TTree *fRecoPi0;

  // Tree variables
  int fRun;
  int fSubRun;
  int fevent;
  double fTimeStamp;
  int fNactivefembs[6];

  // Beam track
  int fbeamtrigger;
  int fbeamCheckIsMatched;
  double ftof;
  int fcerenkovStatus[2];
  double fcerenkovTime[2];
  double fcerenkovPressure[2];
  double fbeamtrackMomentum;
  double fbeamtrackEnergy;
  double fbeamtrackPos[3];
  double fbeamtrackEndPos[3];
  double fbeamtrackDir[3];
  double fbeamtrackTime;
  int fbeam_ntrjPoints;
  int fbeamtrackPdg;
  int fbeamtrackID;
  double fbeamtrackPos_at[1000][3];
  double fbeamtrackMom_at[1000][4];
  std::vector<std::string> fbeamtrackEndProcess;
  int fbeamtrackNDaughters;

  // Reconstructed tracks/showers
  double fprimaryVertex[3];
  int fprimaryIstrack;
  int fprimaryIsshower;
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
  double fprimaryShowerCharge;
  double fprimaryShowerMIPEnergy;
  double fprimaryMomentumByRangeProton;
  int fprimaryIsBeamparticle;
  int    fprimaryTruth_trkID;
  int    fprimaryTruth_pdg;
  double fprimaryTruth_E;
  double fprimaryTruth_vtx[3];
  double fprimaryKineticEnergy[3];
  double fprimaryRange[3];
  int    fprimarynCal;
  double fprimarydEdx[NMAXHITS];
  double fprimarydQdx[NMAXHITS];
  double fprimary_calX[NMAXHITS];
  double fprimary_calY[NMAXHITS];
  double fprimary_calZ[NMAXHITS];
  double fprimary_cal_pitch[NMAXHITS];
  double fprimaryResidualRange[NMAXHITS];
  int    fprimaryShower_nHits; //collection only
  int    fprimaryShower_hit_w[NMAXHITS];
  double fprimaryShower_hit_q[NMAXHITS];
  double fprimaryShower_hit_t[NMAXHITS];
  double fprimaryShower_hit_X[NMAXHITS];
  double fprimaryShower_hit_Y[NMAXHITS];
  double fprimaryShower_hit_Z[NMAXHITS];
//  double fprimaryShower_hit_pitch[NMAXHITS];
  int fprimaryID;
  double fprimaryT0;

  // MC pi0 variables
  int fMCPi0ID;
  double fMCPi0Energy;
  double fMCPi0StartPosition[3];
  double fMCPi0EndPosition[3];
  int fMCPhoton1ID;
  int fMCPhoton2ID;
  double fMCPhoton1Energy;
  double fMCPhoton2Energy;
  double fMCPhoton1StartPosition[3];
  double fMCPhoton2StartPosition[3];
  double fMCPhoton1EndPosition[3];
  double fMCPhoton2EndPosition[3];
  // Reco hits related to photons
  int fMCPhoton1_numhits;
  int fMCPhoton2_numhits;
  int fMCPhoton1_hit_w[NMAXHITS];
  int fMCPhoton2_hit_w[NMAXHITS];
  double fMCPhoton1_hit_t[NMAXHITS];
  double fMCPhoton2_hit_t[NMAXHITS];
  double fMCPhoton1_hit_q[NMAXHITS];
  double fMCPhoton2_hit_q[NMAXHITS];
  double fMCPhoton1_hit_X[NMAXHITS];
  double fMCPhoton2_hit_X[NMAXHITS];
  double fMCPhoton1_hit_Y[NMAXHITS];
  double fMCPhoton2_hit_Y[NMAXHITS];
  double fMCPhoton1_hit_Z[NMAXHITS];
  double fMCPhoton2_hit_Z[NMAXHITS];

  // Reco pi0 variables
  // Showers
  int fShower1ID;
  int fShower2ID;
  double fShower1StartPosition[3];
  double fShower2StartPosition[3];
  double fShower1Direction[3];
  double fShower2Direction[3];
  double fShower1Length;
  double fShower2Length;
  // Reco hits related to showers
  int fShower1_numhits;
  double fShower1_cnn_sc[NMAXHITS];
  double fShower2_cnn_sc[NMAXHITS];
  int fShower2_numhits;
  double fShower1_cal_X[NMAXHITS];
  double fShower2_cal_X[NMAXHITS];
  double fShower1_cal_Y[NMAXHITS];
  double fShower2_cal_Y[NMAXHITS];
  double fShower1_cal_Z[NMAXHITS];
  double fShower2_cal_Z[NMAXHITS];
  double fShower1_cal_E;
  double fShower2_cal_E;
  double fShower1_E;
  double fShower2_E;
  double fShower1_cal_pitch[NMAXHITS];
  double fShower2_cal_pitch[NMAXHITS];
  double fShower1_cal_dEdx[NMAXHITS];
  double fShower2_cal_dEdx[NMAXHITS];
  double fShower1_cal_dQdx[NMAXHITS];
  double fShower2_cal_dQdx[NMAXHITS];
  // Tracks
  int fTrack1ID;
  int fTrack2ID;
  double fTrack1StartPosition[3];
  double fTrack2StartPosition[3];
  double fTrack1EndPosition[3];
  double fTrack2EndPosition[3];
  // Reco hits related to tracks
  int fTrack1_numhits;
  int fTrack2_numhits;
  double fTrack1_cnn_sc[NMAXHITS];
  double fTrack2_cnn_sc[NMAXHITS];
  double fTrack1_cal_X[NMAXHITS];
  double fTrack2_cal_X[NMAXHITS];
  double fTrack1_cal_Y[NMAXHITS];
  double fTrack2_cal_Y[NMAXHITS];
  double fTrack1_cal_Z[NMAXHITS];
  double fTrack2_cal_Z[NMAXHITS];
  double fTrack1_cal_E;
  double fTrack2_cal_E;
  double fTrack1_E;
  double fTrack2_E;
  double fTrack1_cal_pitch[NMAXHITS];
  double fTrack2_cal_pitch[NMAXHITS];
  double fTrack1_cal_dEdx[NMAXHITS];
  double fTrack2_cal_dEdx[NMAXHITS];
  double fTrack1_cal_dQdx[NMAXHITS];
  double fTrack2_cal_dQdx[NMAXHITS];
};


protoana::ProtoDUNEPizeroAnaTree::ProtoDUNEPizeroAnaTree(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  beamlineUtil(p.get<fhicl::ParameterSet>("BeamLineUtils")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fSimulationTag(p.get<std::string>("SimulationTag")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fCalorimetryAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg")),
  fParticleIDTag(p.get<std::string>("ParticleIDTag")),
  fTrackTag(p.get<std::string>("TrackTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fShowerCaloTag(p.get<std::string>("ShowerCalorimetryTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fHitTag(p.get<std::string>("HitTag")),
  fVerbose(p.get<int>("Verbose"))
{

}

void protoana::ProtoDUNEPizeroAnaTree::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;

  fPandoraBeam = tfs->make<TTree>("PandoraBeam", "Beam events reconstructed with Pandora");
  fPandoraBeam->Branch("run",                           &fRun,                          "run/I");
  fPandoraBeam->Branch("subrun",                        &fSubRun,                       "subrun/I");
  fPandoraBeam->Branch("event",                         &fevent,                        "event/I");
  fPandoraBeam->Branch("timestamp",                     &fTimeStamp,                    "timestamp/D");
  fPandoraBeam->Branch("Nactivefembs",                  &fNactivefembs,                 "Nactivefembs[5]/I");

  fPandoraBeam->Branch("beamtrigger",                   &fbeamtrigger,                  "beamtrigger/I");
  fPandoraBeam->Branch("beamCheckIsMatched",            &fbeamCheckIsMatched,           "beamCheckIsMatched/I");
  fPandoraBeam->Branch("tof",                           &ftof,                          "tof/D");
  fPandoraBeam->Branch("cerenkovStatus",                &fcerenkovStatus,               "cerenkovStatus[2]/I");
  fPandoraBeam->Branch("cerenkovTime",                  &fcerenkovTime,                 "cerenkovTime[2]/D");
  fPandoraBeam->Branch("cerenkovPressure",              &fcerenkovPressure,             "cerenkovPressure[2]/D");
  fPandoraBeam->Branch("beamtrackMomentum",             &fbeamtrackMomentum,            "beamtrackMomentum/D");
  fPandoraBeam->Branch("beamtrackEnergy",               &fbeamtrackEnergy,              "beamtrackEnergy/D");
  fPandoraBeam->Branch("beamtrackPos",                  &fbeamtrackPos,                 "beamtrackPos[3]/D");
  fPandoraBeam->Branch("beamtrackEndPos",               &fbeamtrackEndPos,              "beamtrackEndPos[3]/D");
  fPandoraBeam->Branch("beamtrackDir",                  &fbeamtrackDir,                 "beamtrackDir[3]/D");
  fPandoraBeam->Branch("beamtrackTime",                 &fbeamtrackTime,                "beamtrackTime/D");
  fPandoraBeam->Branch("beamtrackPdg",                  &fbeamtrackPdg,                 "beamtrackPdg/I");
  fPandoraBeam->Branch("beamtrackID",                   &fbeamtrackID,                  "beamtrackID/I");
  fPandoraBeam->Branch("beam_ntrjPoints",               &fbeam_ntrjPoints,              "beam_ntrjPoints/I");
  fPandoraBeam->Branch("beamtrackNDaughters",           &fbeamtrackNDaughters,          "beamtrackNDaughters/I");
  fPandoraBeam->Branch("beamtrackPos_at",               &fbeamtrackPos_at,              "beamtrackPos_at[beam_ntrjPoints][3]/D");
  fPandoraBeam->Branch("beamtrackMom_at",               &fbeamtrackMom_at,              "beamtrackMom_at[beam_ntrjPoints][4]/D");
  fPandoraBeam->Branch("beamtrackEndProcess",           &fbeamtrackEndProcess);

  fPandoraBeam->Branch("primaryVertex",                 &fprimaryVertex,                "primaryVertex[3]/D");
  fPandoraBeam->Branch("primaryIsBeamparticle",         &fprimaryIsBeamparticle,        "primaryIsBeamparticle/I");
  fPandoraBeam->Branch("primaryIstrack",                &fprimaryIstrack,               "primaryIstrack/I");
  fPandoraBeam->Branch("primaryIsshower",               &fprimaryIsshower,              "primaryIsshower/I");
  fPandoraBeam->Branch("primaryBDTScore",               &fprimaryBDTScore,              "primaryBDTScore/D");
  fPandoraBeam->Branch("primaryNHits",                  &fprimaryNHits,                 "primaryNHits/I");
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
  fPandoraBeam->Branch("primaryTruth_E",                &fprimaryTruth_E,               "primaryTruth_E/D");
  fPandoraBeam->Branch("primaryTruth_vtx",              &fprimaryTruth_vtx,             "primaryTruth_vtx[3]/D");
  fPandoraBeam->Branch("primaryTruth_pdg",              &fprimaryTruth_pdg,             "primaryTruth_pdg/I");
  fPandoraBeam->Branch("primaryTruth_trkID",            &fprimaryTruth_trkID,           "primaryTruth_trkID/I");
  fPandoraBeam->Branch("primaryShowerBestPlane",        &fprimaryShowerBestPlane,       "primaryShowerBestPlane/I");
  fPandoraBeam->Branch("primaryShowerCharge",           &fprimaryShowerCharge,          "primaryShowerCharge/D");
  fPandoraBeam->Branch("primaryShowerEnergy",           &fprimaryShowerEnergy,          "primaryShowerEnergy/D");
  fPandoraBeam->Branch("primaryShowerMIPEnergy",        &fprimaryShowerMIPEnergy,       "primaryShowerMIPEnergy/D");

  fPandoraBeam->Branch("primaryKineticEnergy",          &fprimaryKineticEnergy,         "primaryKineticEnergy[3]/D");
  fPandoraBeam->Branch("primaryRange",                  &fprimaryRange,                 "primaryRange[3]/D");
  fPandoraBeam->Branch("primarynCal",                   &fprimarynCal,                  "primarynCal/I");
  fPandoraBeam->Branch("primarydQdx",                   &fprimarydQdx,                  "primarydQdx[primarynCal]/D");
  fPandoraBeam->Branch("primary_calX",                  &fprimary_calX,                 "primary_calX[primarynCal]/D");
  fPandoraBeam->Branch("primary_calY",                  &fprimary_calY,                 "primary_calY[primarynCal]/D");
  fPandoraBeam->Branch("primary_calZ",                  &fprimary_calZ,                 "primary_calZ[primarynCal]/D");
  fPandoraBeam->Branch("primary_cal_pitch",             &fprimary_cal_pitch,            "primary_cal_pitch[primarynCal]/D");
  fPandoraBeam->Branch("primarydEdx",                   &fprimarydEdx,                  "primarydEdx[primarynCal]/D");
  fPandoraBeam->Branch("primaryResidualRange",          &fprimaryResidualRange,         "primaryResidualRange[primarynCal]/D");
  fPandoraBeam->Branch("primaryT0",                     &fprimaryT0,                    "primaryT0/D");

  fMCPi0 = tfs->make<TTree>("MCPi0", "MC pi0 information");
  fRecoPi0 = tfs->make<TTree>("RecoPi0", "Reconstructed pi0 information");
  branchPi0TTree(fMCPi0);
  branchPi0TTree(fRecoPi0);
}

void protoana::ProtoDUNEPizeroAnaTree::branchPi0TTree(TTree* tree) {
  // MC pi0 variables
  tree->Branch("MCPi0ID",                 &fMCPi0ID,                "MCPi0ID/I");
  tree->Branch("MCPi0Energy",             &fMCPi0Energy,            "MCPi0Energy/D");
  tree->Branch("MCPi0StartPosition",      &fMCPi0StartPosition,     "MCPi0StartPosition[3]/D");
  tree->Branch("MCPi0EndPosition",        &fMCPi0EndPosition,       "MCPi0EndPosition[3]/D");
  tree->Branch("MCPhoton1ID",             &fMCPhoton1ID,            "MCPhoton1ID/I");
  tree->Branch("MCPhoton2ID",             &fMCPhoton2ID,            "MCPhoton2ID/I");
  tree->Branch("MCPhoton1Energy",         &fMCPhoton1Energy,        "MCPhoton1Energy/D");
  tree->Branch("MCPhoton2Energy",         &fMCPhoton2Energy,        "MCPhoton2Energy/D");
  tree->Branch("MCPhoton1StartPosition",  &fMCPhoton1StartPosition, "MCPhoton1StartPosition[3]/D");
  tree->Branch("MCPhoton2StartPosition",  &fMCPhoton2StartPosition, "MCPhoton2StartPosition[3]/D");
  tree->Branch("MCPhoton1EndPosition",    &fMCPhoton1EndPosition,   "MCPhoton1EndPosition[3]/D");
  tree->Branch("MCPhoton2EndPosition",    &fMCPhoton2EndPosition,   "MCPhoton2EndPosition[3]/D");
  // Reco hits related to photons
  tree->Branch("MCPhoton1_numhits",       &fMCPhoton1_numhits,      "MCPhoton1_numhits/I");
  tree->Branch("MCPhoton2_numhits",       &fMCPhoton2_numhits,      "MCPhoton2_numhits/I");
  tree->Branch("MCPhoton1_hit_w",         &fMCPhoton1_hit_w,        ("MCPhoton1_hit_w[" + std::to_string(NMAXHITS) + "]/I").c_str());
  tree->Branch("MCPhoton2_hit_w",         &fMCPhoton2_hit_w,        ("MCPhoton2_hit_w[" + std::to_string(NMAXHITS) + "]/I").c_str());
  tree->Branch("MCPhoton1_hit_q",         &fMCPhoton1_hit_q,        ("MCPhoton1_hit_q[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("MCPhoton2_hit_q",         &fMCPhoton2_hit_q,        ("MCPhoton2_hit_q[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("MCPhoton1_hit_t",         &fMCPhoton1_hit_t,        ("MCPhoton1_hit_t[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("MCPhoton2_hit_t",         &fMCPhoton2_hit_t,        ("MCPhoton2_hit_t[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("MCPhoton1_hit_X",         &fMCPhoton1_hit_X,        ("MCPhoton1_hit_X[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("MCPhoton2_hit_X",         &fMCPhoton2_hit_X,        ("MCPhoton2_hit_X[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("MCPhoton1_hit_Y",         &fMCPhoton1_hit_Y,        ("MCPhoton1_hit_Y[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("MCPhoton2_hit_Y",         &fMCPhoton2_hit_Y,        ("MCPhoton2_hit_Y[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("MCPhoton1_hit_Z",         &fMCPhoton1_hit_Z,        ("MCPhoton1_hit_Z[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("MCPhoton2_hit_Z",         &fMCPhoton2_hit_Z,        ("MCPhoton2_hit_Z[" + std::to_string(NMAXHITS) + "]/D").c_str());

  // Shower variables
  tree->Branch("Shower1ID",               &fShower1ID,              "Shower1ID/I");
  tree->Branch("Shower2ID",               &fShower2ID,              "Shower2ID/I");
  tree->Branch("Shower1StartPosition",    &fShower1StartPosition,   "Shower1StartPosition[3]/D");
  tree->Branch("Shower2StartPosition",    &fShower2StartPosition,   "Shower2StartPosition[3]/D");
  tree->Branch("Shower1Direction",        &fShower1Direction,       "Shower1Direction[3]/D");
  tree->Branch("Shower2Direction",        &fShower2Direction,       "Shower2Direction[3]/D");
  tree->Branch("Shower1Length",           &fShower1Length,          "Shower1Length/D");
  tree->Branch("Shower2Length",           &fShower2Length,          "Shower2Length/D");
  // Reco hits related to showers
  tree->Branch("Shower1_numhits",         &fShower1_numhits,        "Shower1_numhits/I");
  tree->Branch("Shower2_numhits",         &fShower2_numhits,        "Shower2_numhits/I");
  tree->Branch("Shower1_cnn_sc",          &fShower1_cnn_sc,         ("Shower1_cnn_sc[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower2_cnn_sc",          &fShower2_cnn_sc,         ("Shower2_cnn_sc[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower1_cal_X",           &fShower1_cal_X,          ("Shower1_cal_X[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower2_cal_X",           &fShower2_cal_X,          ("Shower2_cal_X[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower1_cal_Y",           &fShower1_cal_Y,          ("Shower1_cal_Y[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower2_cal_Y",           &fShower2_cal_Y,          ("Shower2_cal_Y[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower1_cal_Z",           &fShower1_cal_Z,          ("Shower1_cal_Z[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower2_cal_Z",           &fShower2_cal_Z,          ("Shower2_cal_Z[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower1_cal_E",           &fShower1_cal_E,          "Shower1_cal_E/D");
  tree->Branch("Shower2_cal_E",           &fShower2_cal_E,          "Shower2_cal_E/D");
  tree->Branch("Shower1_E",               &fShower1_E,              "Shower1_E/D");
  tree->Branch("Shower2_E",               &fShower2_E,              "Shower2_E/D");
  tree->Branch("Shower1_cal_pitch",       &fShower1_cal_pitch,      ("Shower1_cal_pitch[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower2_cal_pitch",       &fShower2_cal_pitch,      ("Shower2_cal_pitch[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower1_cal_dEdx",        &fShower1_cal_dEdx,       ("Shower1_cal_dEdx[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower2_cal_dEdx",        &fShower2_cal_dEdx,       ("Shower2_cal_dEdx[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower1_cal_dQdx",        &fShower1_cal_dQdx,       ("Shower1_cal_dQdx[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Shower2_cal_dQdx",        &fShower2_cal_dQdx,       ("Shower2_cal_dQdx[" + std::to_string(NMAXHITS) + "]/D").c_str());
  // Track variables
  tree->Branch("Track1ID",               &fTrack1ID,              "Track1ID/I");
  tree->Branch("Track2ID",               &fTrack2ID,              "Track2ID/I");
  tree->Branch("Track1StartPosition",    &fTrack1StartPosition,   "Track1StartPosition[3]/D");
  tree->Branch("Track2StartPosition",    &fTrack2StartPosition,   "Track2StartPosition[3]/D");
  tree->Branch("Track1EndPosition",      &fTrack1EndPosition,     "Track1EndPosition[3]/D");
  tree->Branch("Track2EndPosition",      &fTrack2EndPosition,     "Track2EndPosition[3]/D");
  // Reco hits related to showers
  tree->Branch("Track1_numhits",         &fTrack1_numhits,        "Track1_numhits/I");
  tree->Branch("Track2_numhits",         &fTrack2_numhits,        "Track2_numhits/I");
  tree->Branch("Track1_cnn_sc",          &fTrack1_cnn_sc,         ("fTrack1_cnn_sc[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track2_cnn_sc",          &fTrack2_cnn_sc,         ("fTrack2_cnn_sc[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track1_cal_X",           &fTrack1_cal_X,          ("Track1_cal_X[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track2_cal_X",           &fTrack2_cal_X,          ("Track2_cal_X[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track1_cal_Y",           &fTrack1_cal_Y,          ("Track1_cal_Y[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track2_cal_Y",           &fTrack2_cal_Y,          ("Track2_cal_Y[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track1_cal_Z",           &fTrack1_cal_Z,          ("Track1_cal_Z[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track2_cal_Z",           &fTrack2_cal_Z,          ("Track2_cal_Z[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track1_cal_E",           &fTrack1_cal_E,          "Track1_cal_E/D");
  tree->Branch("Track2_cal_E",           &fTrack2_cal_E,          "Track2_cal_E/D");
  tree->Branch("Track1_E",               &fTrack1_E,              "Track1_E/D");
  tree->Branch("Track2_E",               &fTrack2_E,              "Track2_E/D");
  tree->Branch("Track1_cal_pitch",       &fTrack1_cal_pitch,      ("Track1_cal_pitch[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track2_cal_pitch",       &fTrack2_cal_pitch,      ("Track2_cal_pitch[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track1_cal_dEdx",        &fTrack1_cal_dEdx,       ("Track1_cal_dEdx[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track2_cal_dEdx",        &fTrack2_cal_dEdx,       ("Track2_cal_dEdx[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track1_cal_dQdx",        &fTrack1_cal_dQdx,       ("Track1_cal_dQdx[" + std::to_string(NMAXHITS) + "]/D").c_str());
  tree->Branch("Track2_cal_dQdx",        &fTrack2_cal_dQdx,       ("Track2_cal_dQdx[" + std::to_string(NMAXHITS) + "]/D").c_str());
}

void protoana::ProtoDUNEPizeroAnaTree::analyze(art::Event const & evt){
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
  if(!evt.isRealData()){
    for(int k=0; k < 6; k++)
      fNactivefembs[k] = 20;
  }
  else{
    for(int k=0; k < 6; k++)
     fNactivefembs[k] = dataUtil.GetNActiveFembsForAPA(evt, k);
  }

  bool beamTriggerEvent = false;
  if(!evt.isRealData()){
    // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);

    // Also get the reconstructed beam information in the MC - TO DO
    const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    if(geantGoodParticle != 0x0){
      beamTriggerEvent = true;
      fbeamtrigger       = 12;
      fbeamtrackPos[0]   = geantGoodParticle->Vx();
      fbeamtrackPos[1]   = geantGoodParticle->Vy();
      fbeamtrackPos[2]   = geantGoodParticle->Vz();
      fbeamtrackEndPos[0]= geantGoodParticle->EndX();
      fbeamtrackEndPos[1]= geantGoodParticle->EndY();
      fbeamtrackEndPos[2]= geantGoodParticle->EndZ();
      fbeamtrackMomentum = geantGoodParticle->P();
      fbeamtrackEnergy   = geantGoodParticle->E();
      fbeamtrackPdg      = geantGoodParticle->PdgCode();
      fbeamtrackTime     = geantGoodParticle->T();
      fbeamtrackID       = geantGoodParticle->TrackId();
      fbeamtrackEndProcess.push_back(geantGoodParticle->EndProcess());
      fbeam_ntrjPoints = geantGoodParticle->NumberTrajectoryPoints();
      fbeamtrackNDaughters = geantGoodParticle->NumberDaughters();
      for( size_t i=0; i<geantGoodParticle->NumberTrajectoryPoints() && i<1000; ++i){
         fbeamtrackPos_at[i][0] = geantGoodParticle->Position(i).X();
         fbeamtrackPos_at[i][1] = geantGoodParticle->Position(i).Y();
         fbeamtrackPos_at[i][2] = geantGoodParticle->Position(i).Z();
         fbeamtrackMom_at[i][0] = geantGoodParticle->Momentum(i).Px();
         fbeamtrackMom_at[i][1] = geantGoodParticle->Momentum(i).Py();
         fbeamtrackMom_at[i][2] = geantGoodParticle->Momentum(i).Pz();
         fbeamtrackMom_at[i][3] = geantGoodParticle->Momentum(i).E();
      }

      // Record MC and matched reco pi0 information
      // Define and fill a handle to point to a vector of the MCParticles
      art::Handle<std::vector<simb::MCParticle>> MCParticleHandle;
      if (!evt.getByLabel(fSimulationTag, MCParticleHandle)) {
        // Handle no simb::MCParticles.
        throw cet::exception("PiZeroExtract")
            << " No simb::MCParticle objects in this event - "
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }
      // Loop through pi0s
      for(const simb::MCParticle& pi0 : *MCParticleHandle) {
        if(pi0.PdgCode() != 111) continue;
        pizero::PiZeroProcess mcpzproc(pi0, evt, fShowerTag);
        setPiZeroInfo(evt, mcpzproc);
        fMCPi0->Fill();
      }

      // Record pi0s reconstructed from showers
      // Make recob::Showers accessible
      art::Handle<std::vector<recob::Shower>> rShowerHandle;
      if (!evt.getByLabel(fShowerTag, rShowerHandle)) {
        // Handle no recob::Showers.
        throw cet::exception("PiZeroExtract")
            << " No recob::Shower objects in this event - "
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }
      // Loop through showers
      for(const recob::Shower& shower : *rShowerHandle) {
        pizero::PiZeroProcess recopzproc(shower, evt, fShowerTag);
        // Only record if the showers come somewhat close
        if(recopzproc.allRecoSet() &&
           pizero::ClosestDistance(recopzproc.shower1(), recopzproc.shower2()) < 20) {
          setPiZeroInfo(evt, recopzproc);
          fRecoPi0->Fill();
        }
      }

    } //geantGoodParticle
  } // MC
  else{ //data
    // For data we can see if this event comes from a beam trigger
    beamTriggerEvent = dataUtil.IsBeamTrigger(evt);
    if( !beamTriggerEvent ) return;

    art::Handle< std::vector<beam::ProtoDUNEBeamEvent> > pdbeamHandle;
    std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;
    if(evt.getByLabel(fBeamModuleLabel, pdbeamHandle))
      art::fill_ptr_vector(beaminfo, pdbeamHandle);
    else{
       std::cout<<"No beam information from "<<fBeamModuleLabel<<std::endl;
    }

    if(beaminfo.size()){
      if( beamlineUtil.IsGoodBeamlineTrigger( evt ) ){
        fbeamCheckIsMatched = beaminfo[0]->CheckIsMatched();
        fbeamtrigger = beaminfo[0]->GetTimingTrigger();
        if(beaminfo[0]->GetTOFChan() != -1){//if TOFChan == -1, then there was not a successful match, if it's 0, 1, 2, or 3, then there was a good match corresponding to the different pair-wise combinations of the upstream and downstream channels
          ftof = beaminfo[0]->GetTOF();
        }
        // Get Cerenkov
        fcerenkovStatus[0]   = beaminfo[0]->GetCKov0Status();
        fcerenkovStatus[1]   = beaminfo[0]->GetCKov1Status();
        fcerenkovTime[0]     = beaminfo[0]->GetCKov0Time();
        fcerenkovTime[1]     = beaminfo[0]->GetCKov1Time();
        fcerenkovPressure[0] = beaminfo[0]->GetCKov0Pressure();
        fcerenkovPressure[1] = beaminfo[0]->GetCKov1Pressure();
        auto & tracks = beaminfo[0]->GetBeamTracks();
        if(!tracks.empty()){
          fbeamtrackPos[0] = tracks[0].Start().X();
          fbeamtrackPos[1] = tracks[0].Start().Y();
          fbeamtrackPos[2] = tracks[0].Start().Z();
          fbeamtrackEndPos[0] = tracks[0].End().X();
          fbeamtrackEndPos[1] = tracks[0].End().Y();
          fbeamtrackEndPos[2] = tracks[0].End().Z();
          fbeamtrackDir[0] = tracks[0].StartDirection().X();
          fbeamtrackDir[1] = tracks[0].StartDirection().Y();
          fbeamtrackDir[2] = tracks[0].StartDirection().Z();
        }
        auto & beammom = beaminfo[0]->GetRecoBeamMomenta();
        if(!beammom.empty()) fbeamtrackMomentum = beammom[0];

        // Record pi0s reconstructed from showers
        // Make recob::Showers accessible
        art::Handle<std::vector<recob::Shower>> rShowerHandle;
        if (!evt.getByLabel(fShowerTag, rShowerHandle)) {
          // Handle no recob::Showers.
          throw cet::exception("PiZeroExtract")
              << " No recob::Shower objects in this event - "
              << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        // Loop through showers
        for(const recob::Shower& shower : *rShowerHandle) {
          pizero::PiZeroProcess recopzproc(shower, evt, fShowerTag);
          // Only record if the showers come somewhat close
          if(!recopzproc.allRecoSet()) continue;
          const double shdist = pizero::ClosestDistance(recopzproc.shower1(), recopzproc.shower2());
          if(shdist > 0 && shdist < 100) {
            setPiZeroInfo(evt, recopzproc);
            fRecoPi0->Fill();
          }
        }
      } //good beam trigger
    } //good beaminfo
  }//for data

  //check for reco pandora stuff
  art::Handle<std::vector<recob::PFParticle>> recoParticleHandle;
  if( !evt.getByLabel(fPFParticleTag,recoParticleHandle) ) return;

  // Get all of the PFParticles, by default from the "pandora" product
//  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility
  // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
  // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
  // PFParticles considered as primary
  std::vector<const recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  for(const recob::PFParticle* particle : pfParticles){

    FillPrimaryPFParticle(evt, particle);
    const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackTag);

    fprimaryVertex[0] = vtx.X(); fprimaryVertex[1] = vtx.Y(); fprimaryVertex[2] = vtx.Z();
    const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackTag);


    // For now only consider the first primary track. Need a proper treatment if more than one primary particles are found.
    break;
  }
  if(beamTriggerEvent)   fPandoraBeam->Fill();

  // Put some space between events.
  std::cout << "\n\n";
}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEPizeroAnaTree::endJob(){

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEPizeroAnaTree::FillPrimaryPFParticle(art::Event const & evt, const recob::PFParticle* particle){

  // Pandora's BDT beam-cosmic score
  fprimaryBDTScore = (double)pfpUtil.GetBeamCosmicScore(*particle,evt,fPFParticleTag);
  // NHits associated with this pfParticle
  fprimaryNHits = (pfpUtil.GetPFParticleHits(*particle,evt,fPFParticleTag)).size();

  // Get the T0 for this pfParticle
  std::vector<anab::T0> pfT0vec = pfpUtil.GetPFParticleT0(*particle,evt,fPFParticleTag);
  if(!pfT0vec.empty())
    fprimaryT0 = pfT0vec[0].Time();

  // "particle" is the pointer to our beam particle. The recob::Track or recob::Shower object
  // of this particle might be more helpful. These return null pointers if not track-like / shower-like
  const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*particle, evt,fPFParticleTag,fTrackTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);

  const simb::MCParticle* mcparticle = NULL;
  if(thisTrack != 0x0){

    mcparticle = truthUtil.GetMCParticleFromRecoTrack(*thisTrack, evt, fTrackTag);
    if(mcparticle){
      fprimaryTruth_pdg = mcparticle->PdgCode();
      fprimaryTruth_trkID = mcparticle->TrackId();
      fprimaryTruth_vtx[0] =mcparticle->Vx();
      fprimaryTruth_vtx[1] =mcparticle->Vy();
      fprimaryTruth_vtx[2] =mcparticle->Vx();
      fprimaryTruth_E   = mcparticle->E();
      if( fbeamtrackID != -999 && fbeamtrackID == fprimaryTruth_trkID )
        fprimaryIsBeamparticle = 1;
    }
    fprimaryIstrack                    = 1;
    fprimaryIsshower                   = 0;
    fprimaryID                         = thisTrack->ParticleId();
    fprimaryTheta                      = thisTrack->Theta();
    fprimaryPhi                        = thisTrack->Phi();
    fprimaryLength                     = thisTrack->Length();
    fprimaryMomentum                   = thisTrack->StartMomentum();
    fprimaryEndMomentum                = thisTrack->EndMomentum();
    fprimaryEndPosition[0]             = thisTrack->Trajectory().End().X();
    fprimaryEndPosition[1]             = thisTrack->Trajectory().End().Y();
    fprimaryEndPosition[2]             = thisTrack->Trajectory().End().Z();
    fprimaryStartPosition[0]           = thisTrack->Trajectory().Start().X();
    fprimaryStartPosition[1]           = thisTrack->Trajectory().Start().Y();
    fprimaryStartPosition[2]           = thisTrack->Trajectory().Start().Z();
    fprimaryEndDirection[0]            = thisTrack->Trajectory().EndDirection().X();
    fprimaryEndDirection[1]            = thisTrack->Trajectory().EndDirection().Y();
    fprimaryEndDirection[2]            = thisTrack->Trajectory().EndDirection().Z();
    fprimaryStartDirection[0]          = thisTrack->Trajectory().StartDirection().X();
    fprimaryStartDirection[1]          = thisTrack->Trajectory().StartDirection().Y();
    fprimaryStartDirection[2]          = thisTrack->Trajectory().StartDirection().Z();


    // Calorimetry only collection plane
    std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackTag, fCalorimetryTag);
    for(size_t k = 0; k < calovector.size() && k<3; k++){
      int plane = calovector[k].PlaneID().Plane;
      if(plane !=2 ) continue;
      fprimaryKineticEnergy[plane]      = calovector[k].KineticEnergy();
      fprimaryRange[plane]              = calovector[k].Range();
      fprimarynCal = calovector[k].dEdx().size();
      for(size_t l=0; l<calovector[k].dEdx().size() && l<NMAXHITS; ++l){
         fprimarydEdx[l]= calovector[k].dEdx()[l];
         fprimarydQdx[l]= calovector[k].dQdx()[l];
         fprimaryResidualRange[l]= calovector[k].ResidualRange()[l];
         const auto &pos=(calovector[k].XYZ())[l];
         fprimary_calX[l] = pos.X();
         fprimary_calY[l] = pos.Y();
         fprimary_calZ[l] = pos.Z();
         fprimary_cal_pitch[l] =calovector[k].TrkPitchVec()[l];
      }
    }
  } // end is track
  else if(thisShower != 0x0){

    mcparticle = truthUtil.GetMCParticleFromRecoShower(*thisShower, evt, fShowerTag);
    if(mcparticle){
      fprimaryTruth_pdg = mcparticle->PdgCode();
      fprimaryTruth_trkID = mcparticle->TrackId();
      fprimaryTruth_vtx[0] =mcparticle->Vx();
      fprimaryTruth_vtx[1] =mcparticle->Vy();
      fprimaryTruth_vtx[2] =mcparticle->Vx();
      fprimaryTruth_E   = mcparticle->E();
      if( fbeamtrackID != -999 && fbeamtrackID == fprimaryTruth_trkID )
        fprimaryIsBeamparticle = 1;
    }
    fprimaryIstrack                     = 0;
    fprimaryIsshower                    = 1;
    fprimaryID                          = thisShower->ID();
    fprimaryLength                      = thisShower->Length();
    fprimaryShowerBestPlane             = thisShower->best_plane();
    fprimaryOpeningAngle                = thisShower->OpenAngle();
    fprimaryStartPosition[0]            = thisShower->ShowerStart().X();
    fprimaryStartPosition[1]            = thisShower->ShowerStart().Y();
    fprimaryStartPosition[2]            = thisShower->ShowerStart().Z();
    fprimaryStartDirection[0]           = thisShower->Direction().X();
    fprimaryStartDirection[1]           = thisShower->Direction().Y();
    fprimaryStartDirection[2]           = thisShower->Direction().Z();

    // Calorimetry only colleciton plane
    std::vector<anab::Calorimetry> calovector = showerUtil.GetRecoShowerCalorimetry(*thisShower, evt, fShowerTag, fShowerCaloTag);
    if(calovector.size() != 3 && fVerbose > 0)
      std::cerr << "WARNING::Calorimetry vector size for primary is = " << calovector.size() << std::endl;

    for(size_t k = 0; k < calovector.size() && k<3; k++){
      int plane = calovector[k].PlaneID().Plane;
      if(plane !=2 ) continue;
      fprimaryKineticEnergy[plane]      = calovector[k].KineticEnergy();
      fprimaryRange[plane]              = calovector[k].Range();
      fprimarynCal = calovector[k].dEdx().size();
      for(size_t l=0; l<calovector[k].dEdx().size() && l<NMAXHITS; ++l){
         fprimarydEdx[l]= calovector[k].dEdx()[l];
         fprimarydQdx[l]= calovector[k].dQdx()[l];
         fprimaryResidualRange[l]= calovector[k].ResidualRange()[l];
         const auto &pos=(calovector[k].XYZ())[l];
         fprimary_calX[l] = pos.X();
         fprimary_calY[l] = pos.Y();
         fprimary_calZ[l] = pos.Z();
         fprimary_cal_pitch[l] =calovector[k].TrkPitchVec()[l];
      }
    }
  } // end is shower


}

void protoana::ProtoDUNEPizeroAnaTree::setPiZeroInfo(art::Event const & evt, const pizero::PiZeroProcess& pzproc) {
  // Clean slate
  ResetPi0Vars();

  // Set MC object variables
  if(pzproc.allMCSet()) {

    const simb::MCParticle* pi0 = pzproc.pi0();
    const simb::MCParticle* photon1 = pzproc.photon1();
    const simb::MCParticle* photon2 = pzproc.photon2();

    // MC pi0 variables
    fMCPi0ID = pi0->TrackId();
    fMCPi0Energy = pi0->E();
    pi0->Position().Vect().GetXYZ(fMCPi0StartPosition);
    pi0->EndPosition().Vect().GetXYZ(fMCPi0EndPosition);
    fMCPhoton1ID = photon1->TrackId();
    fMCPhoton2ID = photon2->TrackId();
    fMCPhoton1Energy = photon1->E();
    fMCPhoton2Energy = photon2->E();
    photon1->Position().Vect().GetXYZ(fMCPhoton1StartPosition);
    photon2->Position().Vect().GetXYZ(fMCPhoton2StartPosition);
    photon1->EndPosition().Vect().GetXYZ(fMCPhoton1EndPosition);
    photon2->EndPosition().Vect().GetXYZ(fMCPhoton2EndPosition);
    // Reco hits related to photons
    std::vector<const recob::Hit*> ph1hits = truthUtil.GetMCParticleHits(*photon1, evt, fHitTag);
    std::vector<const recob::Hit*> ph2hits = truthUtil.GetMCParticleHits(*photon2, evt, fHitTag);
    // Photon 1 hits
    art::FindManyP<recob::SpacePoint> ph1sps(ph1hits, evt, fPFParticleTag);
    unsigned phi1 = 0;
    for(unsigned i=0; i<ph1hits.size(); ++i) {
      if(ph1hits[i]->WireID().Plane != 2) continue;
      fMCPhoton1_hit_w[phi1] = ph1hits[phi1]->WireID().Wire;
      fMCPhoton1_hit_t[phi1] = ph1hits[phi1]->PeakTime();
      fMCPhoton1_hit_q[phi1] = ph1hits[phi1]->Integral();
      std::vector<art::Ptr<recob::SpacePoint>> sps = ph1sps.at(i);
      if(!sps.empty()) {
        fMCPhoton1_hit_X[phi1] = sps[0]->XYZ()[0];
        fMCPhoton1_hit_Y[phi1] = sps[0]->XYZ()[1];
        fMCPhoton1_hit_Z[phi1] = sps[0]->XYZ()[2];
      }
      ++phi1;
    }
    // Photon 2 hits
    art::FindManyP<recob::SpacePoint> ph2sps(ph2hits, evt, fPFParticleTag);
    unsigned phi2 = 0;
    for(unsigned i=0; i<ph2hits.size(); ++i) {
      if(ph2hits[i]->WireID().Plane != 2) continue;
      fMCPhoton2_hit_w[phi2] = ph2hits[phi2]->WireID().Wire;
      fMCPhoton2_hit_t[phi2] = ph2hits[phi2]->PeakTime();
      fMCPhoton2_hit_q[phi2] = ph2hits[phi2]->Integral();
      std::vector<art::Ptr<recob::SpacePoint>> sps = ph2sps.at(i);
      if(!sps.empty()) {
        fMCPhoton2_hit_X[phi2] = sps[0]->XYZ()[0];
        fMCPhoton2_hit_Y[phi2] = sps[0]->XYZ()[1];
        fMCPhoton2_hit_Z[phi2] = sps[0]->XYZ()[2];
      }
      ++phi2;
    }
    fMCPhoton1_numhits = phi1;
    fMCPhoton2_numhits = phi2;

  } // if MC set

  //This is how we get cnn score for now
  auto recoShowers = evt.getValidHandle< std::vector< recob::Shower > >(fShowerTag);
  auto recoTracks = evt.getValidHandle< std::vector< recob::Track > >(fTrackTag);
  art::FindManyP<recob::Hit> findHitsFromShowers(recoShowers,evt,fShowerTag);
  art::FindManyP<recob::Hit> findHitsFromTracks(recoTracks,evt,fTrackTag);
  anab::MVAReader<recob::Hit,4> hitResults(evt, "emtrkmichelid:emtrkmichel" );
 
  // First shower
  const recob::Shower* shower1 = pzproc.shower1();
  if(shower1 != 0x0) {
    fShower1ID = showerUtil.GetShowerIndex(*shower1, evt, fShowerTag);
    std::vector< art::Ptr< recob::Hit > > sh1_hits = findHitsFromShowers.at(fShower1ID);
    for( size_t i=0; i<sh1_hits.size(); ++i){
       std::array<float,4> cnn_out = hitResults.getOutput( sh1_hits[i] ); 
       double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ]; 
       double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh; 
       fShower1_cnn_sc[i] = cnn_score;
    }
    shower1->ShowerStart().GetXYZ(fShower1StartPosition);
    shower1->Direction().GetXYZ(fShower1Direction);
    fShower1Length = shower1->Length();
    // Calorimetry related to showers
    std::vector<anab::Calorimetry> sh1calo = showerUtil.GetRecoShowerCalorimetry(*shower1, evt, fShowerTag, fShowerCaloTag);
    for(unsigned k=0; k<sh1calo.size(); ++k) {
      if(sh1calo[k].PlaneID().Plane != 2) continue;

      fShower1_numhits = std::min((int)sh1calo[k].dEdx().size(), NMAXHITS);
      fShower1_cal_E = sh1calo[k].KineticEnergy();
      std::vector<const recob::Hit*> shHitPtrs = showerUtil.GetRecoShowerHits(
                  *shower1, evt, fShowerTag);
      fShower1_E = showerUtil.EstimateEnergyFromHitCharge(shHitPtrs, fCalorimetryAlg)[2];
      for(int i=0; i<fShower1_numhits; ++i) {
        const geo::Point_t& pos = sh1calo[k].XYZ()[i];
        fShower1_cal_X[i] = pos.X();
        fShower1_cal_Y[i] = pos.Y();
        fShower1_cal_Z[i] = pos.Z();
        fShower1_cal_pitch[i] = sh1calo[k].TrkPitchVec()[i];
        fShower1_cal_dEdx[i] = sh1calo[k].dEdx()[i];
        fShower1_cal_dQdx[i] = sh1calo[k].dQdx()[i];
      }
    }
  } // if shower1 != 0x0
  // Second shower
  const recob::Shower* shower2 = pzproc.shower2();
  if(shower2 != 0x0) {
    // Reco pi0 variables
    fShower2ID = showerUtil.GetShowerIndex(*shower2, evt, fShowerTag);
    std::vector< art::Ptr< recob::Hit > > sh2_hits = findHitsFromShowers.at(fShower2ID);
    for( size_t i=0; i<sh2_hits.size(); ++i){
       std::array<float,4> cnn_out = hitResults.getOutput( sh2_hits[i] ); 
       double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ]; 
       double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh; 
       fShower2_cnn_sc[i] = cnn_score;
    }
    shower2->ShowerStart().GetXYZ(fShower2StartPosition);
    shower2->Direction().GetXYZ(fShower2Direction);
    fShower2Length = shower2->Length();
    // Calorimetry related to showers
    std::vector<anab::Calorimetry> sh2calo = showerUtil.GetRecoShowerCalorimetry(*shower2, evt, fShowerTag, fShowerCaloTag);
    for(unsigned k=0; k<sh2calo.size(); ++k) {
      if(sh2calo[k].PlaneID().Plane != 2) continue;

      fShower2_numhits = std::min((int)sh2calo[k].dEdx().size(), NMAXHITS);
      fShower2_cal_E = sh2calo[k].KineticEnergy();
      std::vector<const recob::Hit*> shHitPtrs = showerUtil.GetRecoShowerHits(
                  *shower2, evt, fShowerTag);
      fShower2_E = showerUtil.EstimateEnergyFromHitCharge(shHitPtrs, fCalorimetryAlg)[2];
      for(int i=0; i<fShower2_numhits; ++i) {
        const geo::Point_t& pos = sh2calo[k].XYZ()[i];
        fShower2_cal_X[i] = pos.X();
        fShower2_cal_Y[i] = pos.Y();
        fShower2_cal_Z[i] = pos.Z();
        fShower2_cal_pitch[i] = sh2calo[k].TrkPitchVec()[i];
        fShower2_cal_dEdx[i] = sh2calo[k].dEdx()[i];
        fShower2_cal_dQdx[i] = sh2calo[k].dQdx()[i];
      }
    }
  } // if shower2 != 0x0

  // First track
  const recob::Track* track1 = pzproc.track1();
  if(track1 != 0x0) {
    fTrack1ID = track1->ID();
    std::vector< art::Ptr< recob::Hit > > track1_hits = findHitsFromTracks.at(fTrack1ID);
    for( size_t i=0; i<track1_hits.size(); ++i){
       std::array<float,4> cnn_out = hitResults.getOutput( track1_hits[i] ); 
       double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ]; 
       double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh; 
       fTrack1_cnn_sc[i] = cnn_score;
    }
 
    track1->Start<TVector3>().GetXYZ(fTrack1StartPosition);
    track1->End<TVector3>().GetXYZ(fTrack1EndPosition);
    // Calorimetry related to tracks
    std::vector<anab::Calorimetry> tr1calo = trackUtil.GetRecoTrackCalorimetry(*track1, evt, fTrackTag, fCalorimetryTag);
    for(unsigned k=0; k<tr1calo.size(); ++k) {
      if(tr1calo[k].PlaneID().Plane != 2) continue;

      fTrack1_numhits = std::min((int)tr1calo[k].dEdx().size(), NMAXHITS);
      fTrack1_cal_E = tr1calo[k].KineticEnergy();
      std::vector<const recob::Hit*> trHitPtrs = trackUtil.GetRecoTrackHits(
                  *track1, evt, fTrackTag);
      fTrack1_E = showerUtil.EstimateEnergyFromHitCharge(trHitPtrs, fCalorimetryAlg)[2];
      for(int i=0; i<fTrack1_numhits; ++i) {
        const geo::Point_t& pos = tr1calo[k].XYZ()[i];
        fTrack1_cal_X[i] = pos.X();
        fTrack1_cal_Y[i] = pos.Y();
        fTrack1_cal_Z[i] = pos.Z();
        fTrack1_cal_pitch[i] = tr1calo[k].TrkPitchVec()[i];
        fTrack1_cal_dEdx[i] = tr1calo[k].dEdx()[i];
        fTrack1_cal_dQdx[i] = tr1calo[k].dQdx()[i];
      }
    }
  } // if track1 != 0x0
  // Second track
  const recob::Track* track2 = pzproc.track2();
  if(track2 != 0x0) {
    // Reco pi0 variables
    fTrack2ID = track2->ID();
    std::vector< art::Ptr< recob::Hit > > track2_hits = findHitsFromTracks.at(fTrack2ID);
    for( size_t i=0; i<track2_hits.size(); ++i){
       std::array<float,4> cnn_out = hitResults.getOutput( track2_hits[i] ); 
       double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ]; 
       double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh; 
       fTrack2_cnn_sc[i] = cnn_score;
    }
 
    track2->Start<TVector3>().GetXYZ(fTrack2StartPosition);
    track2->End<TVector3>().GetXYZ(fTrack2EndPosition);
    // Calorimetry related to tracks
    std::vector<anab::Calorimetry> tr2calo = trackUtil.GetRecoTrackCalorimetry(*track2, evt, fTrackTag, fCalorimetryTag);
    for(unsigned k=0; k<tr2calo.size(); ++k) {
      if(tr2calo[k].PlaneID().Plane != 2) continue;

      fTrack2_numhits = std::min((int)tr2calo[k].dEdx().size(), NMAXHITS);
      fTrack2_cal_E = tr2calo[k].KineticEnergy();
      std::vector<const recob::Hit*> trHitPtrs = trackUtil.GetRecoTrackHits(
                  *track2, evt, fTrackTag);
      fTrack2_E = showerUtil.EstimateEnergyFromHitCharge(trHitPtrs, fCalorimetryAlg)[2];
      for(int i=0; i<fTrack2_numhits; ++i) {
        const geo::Point_t& pos = tr2calo[k].XYZ()[i];
        fTrack2_cal_X[i] = pos.X();
        fTrack2_cal_Y[i] = pos.Y();
        fTrack2_cal_Z[i] = pos.Z();
        fTrack2_cal_pitch[i] = tr2calo[k].TrkPitchVec()[i];
        fTrack2_cal_dEdx[i] = tr2calo[k].dEdx()[i];
        fTrack2_cal_dQdx[i] = tr2calo[k].dQdx()[i];
      }
    }
  } // if track2 != 0x0

} // setPiZeroInfo
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void protoana::ProtoDUNEPizeroAnaTree::Initialise(){
  fRun = -999;
  fSubRun = -999;
  fevent = -999;
  fTimeStamp = -999.0;
  for(int k=0; k < 5; k++)
    fNactivefembs[k] = -999;

  for(int k=0; k < 3; k++){
    fbeamtrackPos[k] = -999.0;
    fbeamtrackEndPos[k] = -999.0;
    fbeamtrackDir[k] = -999.0;
    fprimaryTruth_vtx[k]= -999.0;
    fprimaryVertex[k] = -999.0;
    fprimaryEndPosition[k] = -999.0;
    fprimaryStartPosition[k] = -999.0;
    fprimaryEndDirection[k] = -999.0;
    fprimaryStartDirection[k] = -999.0;
    fprimaryKineticEnergy[k] = -999.0;
    fprimaryRange[k] = -999.0;
  }
  for( int l=0; l<1000; l ++) {
     fbeamtrackPos_at[l][0] = -999.;
     fbeamtrackPos_at[l][1] = -999.;
     fbeamtrackPos_at[l][2] = -999.;
     fbeamtrackMom_at[l][0] = -999.;
     fbeamtrackMom_at[l][1] = -999.;
     fbeamtrackMom_at[l][2] = -999.;
     fbeamtrackMom_at[l][3] = -999.;
  }
  fprimaryShower_nHits =0;
  for( int m=0; m<NMAXHITS; m ++){
     fprimarydEdx[m]= -999.0;
     fprimarydQdx[m]= -999.0;
     fprimary_calX[m] =-999.0;
     fprimary_calY[m] =-999.0;
     fprimary_calZ[m] =-999.0;
     fprimary_cal_pitch[m]=-999.0;
     fprimaryResidualRange[m] = -999.0;

     fprimaryShower_hit_w[m] =-999.0;
     fprimaryShower_hit_q[m] =-999.0;
     fprimaryShower_hit_t[m] =-999.0;
     fprimaryShower_hit_X[m] =-999.0;
     fprimaryShower_hit_Y[m] =-999.0;
     fprimaryShower_hit_Z[m] =-999.0;
  }
  fprimaryTruth_trkID =-999;
  fprimaryTruth_pdg = -999;
  fprimaryTruth_E = -999;
  fprimarynCal = 0;
  fbeamCheckIsMatched = -999;
  fbeamtrigger = -999;
  ftof = -999.0;
  for(int k=0; k < 2; k++){
    fcerenkovPressure[k] = -999.0;
    fcerenkovTime[k] = -999.0;
    fcerenkovStatus[k] = -999;
  }
  fbeamtrackMomentum = -999.0;
  fbeamtrackEnergy = 999.0;
  fbeamtrackPdg = -999;
  fbeamtrackTime = -999.0;
  fbeamtrackID = -999;
  fbeam_ntrjPoints =0;
  fbeamtrackNDaughters= 0;
  fprimaryIstrack = -999;
  fprimaryIsshower = -999;

  fprimaryBDTScore = -999.0;
  fprimaryNHits = 0;
  fprimaryTheta = -999.0;
  fprimaryPhi = -999.0;
  fprimaryLength = -999.0;
  fprimaryMomentum = -999.0;
  fprimaryEndMomentum = -999.0;
  fprimaryOpeningAngle = -999.0;
  fprimaryShowerBestPlane = -999;
  fprimaryShowerEnergy = -999.0;
  fprimaryShowerCharge = -999.0;
  fprimaryShowerMIPEnergy = -999.0;
  fprimaryID = -999;
  fprimaryMomentumByRangeProton = -999.0;
  fprimaryT0 = -999.0;

  ResetPi0Vars();
} // Initialise

void protoana::ProtoDUNEPizeroAnaTree::ResetPi0Vars() {
  // MC pi0 variables
  fMCPi0ID = -999;
  fMCPi0Energy = -999.0;
  for(int i=0; i<3; ++i) {
    fMCPi0StartPosition[i] = -999.0;
    fMCPi0EndPosition[i] = -999.0;
  }
  fMCPhoton1ID = -999;
  fMCPhoton2ID = -999;
  fMCPhoton1Energy = -999.0;
  fMCPhoton2Energy = -999.0;
  for(int i=0; i<3; ++i) {
    fMCPhoton1StartPosition[i] = -999.0;
    fMCPhoton2StartPosition[i] = -999.0;
    fMCPhoton1EndPosition[i] = -999.0;
    fMCPhoton2EndPosition[i] = -999.0;
  }
  // Reco hits related to photons
  fMCPhoton1_numhits = -999;
  fMCPhoton2_numhits = -999;
  for(int i=0; i < NMAXHITS; ++i) {
    fMCPhoton1_hit_w[i] = -999;
    fMCPhoton2_hit_w[i] = -999;
    fMCPhoton1_hit_t[i] = -999.0;
    fMCPhoton2_hit_t[i] = -999.0;
    fMCPhoton1_hit_q[i] = -999.0;
    fMCPhoton2_hit_q[i] = -999.0;
    fMCPhoton1_hit_X[i] = -999.0;
    fMCPhoton2_hit_X[i] = -999.0;
    fMCPhoton1_hit_Y[i] = -999.0;
    fMCPhoton2_hit_Y[i] = -999.0;
    fMCPhoton1_hit_Z[i] = -999.0;
    fMCPhoton2_hit_Z[i] = -999.0;
    // fMCPhoton1_hit_E[i] = -999.0;
    // fMCPhoton2_hit_E[i] = -999.0;
    // fMCPhoton1_hit_pitch[i] = -999.0;
    // fMCPhoton2_hit_pitch[i] = -999.0;
    // fMCPhoton1_hit_dEdx[i] = -999.0;
    // fMCPhoton2_hit_dEdx[i] = -999.0;
    // fMCPhoton1_hit_dQdx[i] = -999.0;
    // fMCPhoton2_hit_dQdx[i] = -999.0;
  }

  // Reco pi0 variables
  fShower1ID = -999;
  fShower2ID = -999;
  for(int i=0; i<3; ++i) {
    fShower1StartPosition[i] = -999.0;
    fShower2StartPosition[i] = -999.0;
    fShower1Direction[i] = -999.0;
    fShower2Direction[i] = -999.0;
  }
  fShower1Length = -999.0;
  fShower2Length = -999.0;
  // Reco hits related to showers
  fShower1_numhits = -999;
  fShower2_numhits = -999;
  fShower1_cal_E = -999.0;
  fShower2_cal_E = -999.0;
  fShower1_E = -999.0;
  fShower2_E = -999.0;
  for(int i=0; i<NMAXHITS; ++i) {
    fShower1_cnn_sc[i] = -999.0;
    fShower2_cnn_sc[i] = -999.0;
    fShower1_cal_X[i] = -999.0;
    fShower2_cal_X[i] = -999.0;
    fShower1_cal_Y[i] = -999.0;
    fShower2_cal_Y[i] = -999.0;
    fShower1_cal_Z[i] = -999.0;
    fShower2_cal_Z[i] = -999.0;
    fShower1_cal_pitch[i] = -999.0;
    fShower2_cal_pitch[i] = -999.0;
    fShower1_cal_dEdx[i] = -999.0;
    fShower2_cal_dEdx[i] = -999.0;
    fShower1_cal_dQdx[i] = -999.0;
    fShower2_cal_dQdx[i] = -999.0;
  }
  // Track variables
  fTrack1ID = -999;
  fTrack2ID = -999;
  for(int i=0; i<3; ++i) {
    fTrack1StartPosition[i] = -999.0;
    fTrack2StartPosition[i] = -999.0;
    fTrack1EndPosition[i] = -999.0;
    fTrack2EndPosition[i] = -999.0;
  }
  // Reco hits related to showers
  fTrack1_numhits = -999;
  fTrack2_numhits = -999;
  fTrack1_cal_E = -999.0;
  fTrack2_cal_E = -999.0;
  fTrack1_E = -999.0;
  fTrack2_E = -999.0;
  for(int i=0; i<NMAXHITS; ++i) {
    fTrack1_cnn_sc[i] = -999.0;
    fTrack2_cnn_sc[i] = -999.0;
    fTrack1_cal_X[i] = -999.0;
    fTrack2_cal_X[i] = -999.0;
    fTrack1_cal_Y[i] = -999.0;
    fTrack2_cal_Y[i] = -999.0;
    fTrack1_cal_Z[i] = -999.0;
    fTrack2_cal_Z[i] = -999.0;
    fTrack1_cal_pitch[i] = -999.0;
    fTrack2_cal_pitch[i] = -999.0;
    fTrack1_cal_dEdx[i] = -999.0;
    fTrack2_cal_dEdx[i] = -999.0;
    fTrack1_cal_dQdx[i] = -999.0;
    fTrack2_cal_dQdx[i] = -999.0;
  }
} // ResetPi0Vars

DEFINE_ART_MODULE(protoana::ProtoDUNEPizeroAnaTree)
