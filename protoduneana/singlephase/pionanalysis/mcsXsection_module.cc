////////////////////////////////////////////////////////////////////////
// Class:       mcsXsection                                              //
// File:        mcsXsection_module.cc                                    //
// Ajib Paudel MCS cross section module                               //   
// ajib.paudel@phys.ksu.edu                                           //
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "larcorealg/Geometry/Exceptions.h"
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

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
//#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEbeamsim.h"
//#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEBeamInstrument.h"
//#include "dune/Protodune/singlephase/DataUtils/ProtoDUNETrackUtils.h"
//#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEShowerUtils.h"
//#include "dune/Protodune/singlephase/DataUtils/ProtoDUNETruthUtils.h"
//#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEPFParticleUtils.h"
//#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamCuts.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TH2.h"
// C++ Includes
#include <stdio.h>
#include <stdlib.h> 
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include "TComplex.h"
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include <iterator>
#include <map>


// Maximum number of beam particles to save
const int NMAXDAUGTHERS = 30;
double  m_pi=0.1395;//GeV/c^2
double prim_energy=0.0;
//double geant_energy=0;


float soln(float x1, float x2, float x3, float x4, float y1, float y2, float y3, float y4, float result[2]){
  if(x3==x4){
    result[0]=x3;
    result[1]=y3;
  }
  else if((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)==0){
    result[0]=x3;
    result[1]=y3;
  }
  else{
    result[0]=((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
    result[1]=((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
  }
  return *result;
}
double angle2d(double x0, double y0, double x1, double y1, double x2, double y2){
  double theta=((x1-x0)*(x2-x1)+(y1-y0)*(y2-y1))/sqrt(((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))*((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
  return acos(theta);

} 


namespace protoana {
  class mcsXsection;
}


class protoana::mcsXsection : public art::EDAnalyzer {
public:

  explicit mcsXsection(fhicl::ParameterSet const & p);

  mcsXsection(mcsXsection const &) = delete;
  mcsXsection(mcsXsection &&) = delete;
  mcsXsection & operator = (mcsXsection const &) = delete;
  mcsXsection & operator = (mcsXsection &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & evt) override;

private:
  const art::InputTag fTrackModuleLabel;
  const art::InputTag fBeamModuleLabel;
  // Helper utility functions
  protoana::ProtoDUNEBeamCuts beam_cuts;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
  protoana::ProtoDUNEDataUtils dataUtil;

  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETrackUtils trackUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
 

  // Initialise tree variables
  void Initialise();

  // Fill cosmics tree
  void FillCosmicsTree(art::Event const & evt, std::string pfParticleTag);

  // fcl parameters
  //const art::InputTag fNNetModuleLabel; //label of the module used for CNN tagging
 
 
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  bool fVerbose;

  geo::GeometryCore const * fGeometry;

  TTree *fPandoraBeam;
  //TTree *fPandoraCosmics;
  // Tree variables
  int fRun;
  int fSubRun;
  int fevent;
  double fTimeStamp;
 

  double fbeamtrackMomentum;
  double fbeamtrackP[3]; //Px/Py/Pz
  double fbeamtrackEnergy;
  double fbeamtrackPos[3];
  double fbeamtrackDir[3];
  double fbeamtrackTime;
  int fbeamtrackPdg;
  int fbeamtrackID;

  //int NumberBeamTrajectoryPoints;
  std::vector<float> beamtrk_x;
  std::vector<float> beamtrk_y;
  std::vector<float> beamtrk_z;
  std::vector<int> beamtrk_z_wire;
  std::vector<int> beamtrk_z_tpc;

  std::vector<float> beamtrk_Px;
  std::vector<float> beamtrk_Py;
  std::vector<float> beamtrk_Pz;
  std::vector<float> beamtrk_Eng;
  std::vector<float> peakT_2;
  std::vector<float> hitz_2;
  // std::vector<float> inelscore;
  //std::vector<float> elscore;
  //std::vector<float> nonescore;


  std::vector<float> int_2;
  std::vector<float> hitz_1;
  std::vector<float> int_1;
  std::vector<float> hitz_0;
  std::vector<float> int_0;

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




  //*********************truth information*******************************************//
  int fprimary_truth_TrackId;
  int fprimary_truth_Pdg;
  int ftruthpdg;
  double fprimary_truth_StartPosition[4];
  // double fprimary_truth_Startwiretime[4];
  double fprimary_truth_EndPosition[4];
  std::string fprimary_truth_EndProcess;
  std::string truth_last_process;
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
  double fprimary_truth_tracklength;
  

  //interaction point details
  std::vector<double> interactionX;
  std::vector<double> interactionY;
  std::vector<double> interactionZ;
  std::vector<double> interactionU;
  std::vector<double> interactionV;
  std::vector<double> interactionW;
  std::vector<double> interactionT;



  std::vector<std::string> interactionProcesslist;
  std::vector<double> interactionAngles;
  std::vector<double> interactionAnglesUT;
  std::vector<double> interactionAnglesVT;
  std::vector<double> interactionAnglesZT;
  std::vector<double> Zintersection;
  std::vector<double> Zintersection1;
  std::vector<double> timeintersection;
  std::vector<double> timeintersection1;
  std::vector<double> deltaZint;
  std::vector<double> deltatimeint;
  //calo info
  std::vector< std::vector<double> > primtrk_dqdx;
  std::vector< std::vector<double> > primtrk_resrange;
  std::vector< std::vector<double> > primtrk_dedx;
  std::vector<double> primtrk_range;
  std::vector< std::vector<double> > primtrk_hitx;
  std::vector< std::vector<double> > primtrk_hity;
  std::vector< std::vector<double> > primtrk_hitz;
  std::vector< std::vector<int> > primtrk_hitz_wire;
  std::vector< std::vector<int> > primtrk_hitz_tpc;

  std::vector< std::vector<double> > primtrk_pitch;
  std::vector<std::vector<double> > primtrk_truth_Z;
  std::vector<std::vector<int> > primtrk_truth_Z_wire;
  std::vector<std::vector<int> > primtrk_truth_Z_tpc;


  std::vector<std::vector<double> > primtrk_truth_Eng;
  std::vector<std::vector<double> > primtrk_truth_trkide;
  std::vector<std::vector<float> > wireno_0;
  std::vector<std::vector<float> > wireno_1;
  std::vector<std::vector<float> > wireno_2;
  std::vector<std::vector<float> > peakTime_0;
  std::vector<std::vector<float> > peakTime_1;
  std::vector<std::vector<float> > peakTime_2;
  std::vector<std::vector<float> > dq_0;
  std::vector<std::vector<float> > dq_1;
  std::vector<std::vector<float> > dq_2;

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
  //double fdaughterT0[NMAXDAUGTHERS];

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

  double minX =  -360.0;
  double maxX = 360.0;
  double minY =0.0;
  double maxY = 600.0;
  double minZ =  0.0;
  double maxZ = 695.0;

};


protoana::mcsXsection::mcsXsection(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
 
  // fNNetModuleLabel(p.get<art::InputTag>("NNetModuleLabel")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  beam_cuts(p.get<fhicl::ParameterSet>("BeamCuts")),
  fBeamlineUtils(p.get<fhicl::ParameterSet>("BeamlineUtils")),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fVerbose(p.get<bool>("Verbose"))
{

}

void protoana::mcsXsection::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;

  fPandoraBeam = tfs->make<TTree>("PandoraBeam", "Beam events reconstructed with Pandora");
  fPandoraBeam->Branch("run",                           &fRun,                          "run/I");
  fPandoraBeam->Branch("subrun",                        &fSubRun,                       "subrun/I");
  fPandoraBeam->Branch("event",                         &fevent,                        "event/I");
  fPandoraBeam->Branch("timestamp",                     &fTimeStamp,                    "timestamp/D");
  fPandoraBeam->Branch("beamtrackMomentum",             &fbeamtrackMomentum,            "beamtrackMomentum/D");
  fPandoraBeam->Branch("beamtrackP",                    &fbeamtrackP,                   "beamtrackP[3]/D");
  fPandoraBeam->Branch("beamtrackEnergy",               &fbeamtrackEnergy,              "beamtrackEnergy/D");
  fPandoraBeam->Branch("beamtrackPos",                  &fbeamtrackPos,                 "beamtrackPos[3]/D");
  fPandoraBeam->Branch("beamtrackDir",                  &fbeamtrackDir,                 "beamtrackDir[3]/D");
  fPandoraBeam->Branch("beamtrackTime",                 &fbeamtrackTime,                "beamtrackTime/D");
  fPandoraBeam->Branch("beamtrackPdg",                  &fbeamtrackPdg,                 "beamtrackPdg/I");
  fPandoraBeam->Branch("beamtrackID",                   &fbeamtrackID,                  "beamtrackID/I");

  //fPandoraBeam->Branch("NumberBeamTrajectoryPoints",    &NumberBeamTrajectoryPoints,    "NumberBeamTrajectoryPoints/I");
  fPandoraBeam->Branch("beamtrk_x",&beamtrk_x);
  fPandoraBeam->Branch("beamtrk_y",&beamtrk_y);
  fPandoraBeam->Branch("beamtrk_z",&beamtrk_z);
  fPandoraBeam->Branch("beamtrk_z_wire",&beamtrk_z_wire);
  fPandoraBeam->Branch("beamtrk_z_tpc",&beamtrk_z_tpc);

  fPandoraBeam->Branch("beamtrk_Px",&beamtrk_Px);
  fPandoraBeam->Branch("beamtrk_Py",&beamtrk_Py);
  fPandoraBeam->Branch("beamtrk_Pz",&beamtrk_Pz);
  fPandoraBeam->Branch("beamtrk_Eng",&beamtrk_Eng);


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
 

  fPandoraBeam->Branch("primary_truth_TrackId",         &fprimary_truth_TrackId,         "primary_truth_TrackId/I");
  fPandoraBeam->Branch("primary_truth_Pdg",             &fprimary_truth_Pdg,             "primary_truth_Pdg/I");
  fPandoraBeam->Branch("truthpdg",                      &ftruthpdg,                      "truthpdg/I");


  fPandoraBeam->Branch("primary_truth_StartPosition",   &fprimary_truth_StartPosition,   "primary_truth_StartPosition[4]/D");
  // fPandoraBeam->Branch("primary_truth_Startwiretime",   &fprimary_truth_Startwiretime,   "primary_truth_Startwiretime[4]/D");
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
  fPandoraBeam->Branch("primary_truth_EndProcess",      &fprimary_truth_EndProcess);
  fPandoraBeam->Branch("truth_last_process",            &truth_last_process);
  fPandoraBeam->Branch("primary_truth_tracklength",     &fprimary_truth_tracklength,      "primary_truth_tracklength/D");

  fPandoraBeam->Branch("interactionX",&interactionX);
  fPandoraBeam->Branch("interactionY",&interactionY);
  fPandoraBeam->Branch("interactionZ",&interactionZ);
  fPandoraBeam->Branch("interactionU",&interactionU);
  fPandoraBeam->Branch("interactionV",&interactionV);
  fPandoraBeam->Branch("interactionW",&interactionW);
  fPandoraBeam->Branch("interactionT",&interactionT);
  fPandoraBeam->Branch("interactionProcesslist",&interactionProcesslist);
  fPandoraBeam->Branch("interactionAngles",&interactionAngles);
  fPandoraBeam->Branch("interactionAnglesUT",&interactionAnglesUT);
  fPandoraBeam->Branch("interactionAnglesVT",&interactionAnglesVT);
  fPandoraBeam->Branch("interactionAnglesZT",&interactionAnglesZT);

  fPandoraBeam->Branch("Zintersection",&Zintersection);
  fPandoraBeam->Branch("timeintersection",&timeintersection);
  fPandoraBeam->Branch("Zintersection1",&Zintersection1);
  fPandoraBeam->Branch("timeintersection1",&timeintersection1);
  fPandoraBeam->Branch("deltaZint",&deltaZint);
  fPandoraBeam->Branch("deltatimeint",&deltatimeint);
  fPandoraBeam->Branch("peakT_2",&peakT_2);
  fPandoraBeam->Branch("hitz_2",&hitz_2);
  // fPandoraBeam->Branch("inelscore",&inelscore);
  //fPandoraBeam->Branch("elscore",&elscore);
  //fPandoraBeam->Branch("nonescore",&nonescore);

  fPandoraBeam->Branch("int_2",&int_2);
  fPandoraBeam->Branch("hitz_1",&hitz_1);
  fPandoraBeam->Branch("int_1",&int_1);
  fPandoraBeam->Branch("hitz_0",&hitz_0);
  fPandoraBeam->Branch("int_0",&int_0);
  fPandoraBeam->Branch("wireno_0",&wireno_0);
  fPandoraBeam->Branch("wireno_1",&wireno_1);
  fPandoraBeam->Branch("wireno_2",&wireno_2);
  fPandoraBeam->Branch("peakTime_0",&peakTime_0);
  fPandoraBeam->Branch("peakTime_1",&peakTime_1);
  fPandoraBeam->Branch("peakTime_2",&peakTime_2);
  fPandoraBeam->Branch("dq_0",&dq_0);
  fPandoraBeam->Branch("dq_1",&dq_1);
  fPandoraBeam->Branch("dq_2",&dq_2);

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
  // fPandoraBeam->Branch("daughterT0",                    &fdaughterT0,                   "daughterT0[NDAUGHTERS]/D");

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
 

  // hKE_truth_reco = tfs->make<TH2D>("hKE_truth_reco","KE truth vs reco;true KE in MeV;reconstructed KE in MeV"  , 3000, -100, 1400  ,3000, -100, 1400); 
  // hEndZ_truth_reco = tfs->make<TH2D>("hEndZ_truth_reco","EndZ truth vs reco;true EndZ (cm);reconstructed EndZ(cm)"  , 695, 0, 695  ,695,0,695); 

  fPandoraBeam->Branch("primtrk_dqdx",&primtrk_dqdx);
  fPandoraBeam->Branch("primtrk_dedx",&primtrk_dedx);
  fPandoraBeam->Branch("primtrk_resrange",&primtrk_resrange);
  fPandoraBeam->Branch("primtrk_range",&primtrk_range);
  fPandoraBeam->Branch("primtrk_hitx",&primtrk_hitx);
  fPandoraBeam->Branch("primtrk_hity",&primtrk_hity);
  fPandoraBeam->Branch("primtrk_hitz",&primtrk_hitz);
  fPandoraBeam->Branch("primtrk_hitz_wire",&primtrk_hitz_wire);
  fPandoraBeam->Branch("primtrk_hitz_tpc",&primtrk_hitz_tpc);
  fPandoraBeam->Branch("primtrk_pitch",&primtrk_pitch);
  ///////////////////////////////////////////////////
  fPandoraBeam->Branch("primtrk_truth_Z",&primtrk_truth_Z);//primary track true Z positions 
  fPandoraBeam->Branch("primtrk_truth_Z_wire",&primtrk_truth_Z_wire);//primary track true Z positions 
  fPandoraBeam->Branch("primtrk_truth_Z_tpc",&primtrk_truth_Z_tpc);//primary track true Z positions 


  fPandoraBeam->Branch("primtrk_truth_Eng",&primtrk_truth_Eng);//primary track true Energy deposited for each Z position
  fPandoraBeam->Branch("primtrk_truth_trkide",&primtrk_truth_trkide);//primary track true Energy deposited for each Z position
}//begin job 

void protoana::mcsXsection::analyze(art::Event const & evt){

  // Initialise tree parameters
  Initialise();

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  fGeometry = &*(art::ServiceHandle<geo::Geometry>());
  // const sim::ParticleList& plist=pi_serv->ParticleList();

  const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // anab::MVAReader<recob::Hit,3> hitResults(evt, fNNetModuleLabel);
  art::ServiceHandle<geo::Geometry> geom;


  int beamid=-9999;
  int truthid=-999;

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

  bool beamTriggerEvent = false;
  // If this event is MC then we can check what the true beam particle is
  auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
  if(!evt.isRealData()){
    // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
    // of these, so we pass the first element into the function to get the good particle
    const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
    //std::cout<<"geantGoodParticle "<<geantGoodParticle.size()<<std::endl;
    if(geantGoodParticle != 0x0){
      std::cout<<"geant good particle loop "<<std::endl;
      std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() 
		<< " , track id = " << geantGoodParticle->TrackId()
		<< " , Vx/Vy/Vz = " << geantGoodParticle->Vx() << "/"<< geantGoodParticle->Vy() << "/" << geantGoodParticle->Vz() 

		<< std::endl;

      std::vector<double> tmp_primtrk_truth_Z;
      std::vector<double> tmp_primtrk_truth_Eng;
      std::vector<double> tmp_primtrk_truth_trkide;     
      std::vector<int> tmp_wire;
      std::vector<int> tmp_tpc;      

      beamTriggerEvent = true;
      fbeamtrackPos[0]   = geantGoodParticle->Vx();
      fbeamtrackPos[1]   = geantGoodParticle->Vy();
      fbeamtrackPos[2]   = geantGoodParticle->Vz();
      fbeamtrackMomentum = geantGoodParticle->P();
      fbeamtrackP[0]     = geantGoodParticle->Px();
      fbeamtrackP[1]     = geantGoodParticle->Py();
      fbeamtrackP[2]     = geantGoodParticle->Pz();
      // fbeamtrackEnergy   = geantGoodParticle->E();
      fbeamtrackPdg      = geantGoodParticle->PdgCode();
      fbeamtrackTime     = geantGoodParticle->T();
      fbeamtrackID       = geantGoodParticle->TrackId();
      beamid             = geantGoodParticle->TrackId();

      /* fprimary_truth_TrackId          = geantGoodParticle->TrackId();
	 fprimary_truth_StartPosition[3] = geantGoodParticle->T();
	 fprimary_truth_EndPosition[3]   = geantGoodParticle->EndT();
	 fprimary_truth_P                = geantGoodParticle->P();
	 fprimary_truth_Momentum[3]      = geantGoodParticle->E();
	 fprimary_truth_Pt               = geantGoodParticle->Pt();
	 fprimary_truth_Mass             = geantGoodParticle->Mass();
	 fprimary_truth_EndMomentum[0]   = geantGoodParticle->EndPx();
	 fprimary_truth_EndMomentum[1]   = geantGoodParticle->EndPy();
	 fprimary_truth_EndMomentum[2]   = geantGoodParticle->EndPz();
	 fprimary_truth_EndMomentum[3]   = geantGoodParticle->EndE();
	 fprimary_truth_Theta            = geantGoodParticle->Momentum().Theta();
	 fprimary_truth_Phi              = geantGoodParticle->Momentum().Phi();
	 fprimary_truth_NDaughters       = geantGoodParticle->NumberDaughters();*/

      ////////////////////////////////////////
      truth_last_process=geantGoodParticle->EndProcess();
      ////here is another new parameter added

      fprimary_truth_Process = int(geantGoodParticle->Trajectory().ProcessToKey(geantGoodParticle->Process()));
      prim_energy=0;
      for(size_t i_s=0; i_s < geantGoodParticle->NumberTrajectoryPoints(); i_s++){ //loop over beam tracks
	//	if(geantGoodParticle->Position(i_s).Z()>0) break;
	beamtrk_x.push_back(geantGoodParticle->Position(i_s).X());
	beamtrk_y.push_back(geantGoodParticle->Position(i_s).Y());
	beamtrk_z.push_back(geantGoodParticle->Position(i_s).Z());
	double pos_true[3]={geantGoodParticle->Position(i_s).X(),geantGoodParticle->Position(i_s).Y(),geantGoodParticle->Position(i_s).Z()};


	/*	if(geantGoodParticle->Position(i_s).Z()>=0 && geantGoodParticle->Position(i_s).Z()<=230){
		geo::WireID wireID = geom->NearestWireID(pos_true, 2);
		if (!wireID) wireID = geom->Plane(2).ClosestWireID(wireID);
		beamtrk_z_wire.push_back(wireID.Wire);
		beamtrk_z_tpc.push_back(wireID.TPC);
		}
		if(geantGoodParticle->Position(i_s).Z()<0 || geantGoodParticle->Position(i_s).Z()>230){
		beamtrk_z_wire.push_back(-99999);
		beamtrk_z_tpc.push_back(-99999);

		}*/
	geo::TPCID tpc = geom->FindTPCAtPosition(pos_true);
	if(tpc.isValid){
	  int tpc_no=tpc.TPC;
	  geo::PlaneID planeID = geo::PlaneID(0, tpc_no, 2);
	  geo::WireID wireID;
	  //	geo::WireID wireID = geom->NearestWireID(pos_true, planeID);
	  try{
	    wireID = geom->NearestWireID(pos_true, planeID);
	  }
	  catch(geo::InvalidWireError const& e) {
	    wireID = e.suggestedWireID(); // pick the closest valid wire
	  }
	  beamtrk_z_wire.push_back(wireID.Wire);
	  beamtrk_z_tpc.push_back(wireID.TPC);
	}
	if(!tpc.isValid){
	  beamtrk_z_wire.push_back(-9999);
	  beamtrk_z_tpc.push_back(-9999);

	}


	


	/*	geo::PlaneGeo const& plane = geom->Plane(2);
		geo::WireID wireID;
		try {
		wireID = plane.NearestWireID(pos_true);
		}
		catch (geo::InvalidWireIDError const& e) {
		if (!e.hasSuggestedWire()) throw;
		wireID = plane.ClosestWireID(e.suggestedWireID());
		}
		beamtrk_z_wire.push_back(wireID.Wire);
		beamtrk_z_tpc.push_back(wireID.TPC);*/
	///////trying new stuff
        


	beamtrk_Px.push_back(geantGoodParticle->Momentum(i_s).X());
        beamtrk_Py.push_back(geantGoodParticle->Momentum(i_s).Y());
        beamtrk_Pz.push_back(geantGoodParticle->Momentum(i_s).Z());
	
        beamtrk_Eng.push_back(geantGoodParticle->Momentum(i_s).E()-geantGoodParticle->Mass());
	if(geantGoodParticle->Position(i_s).Z()<0) prim_energy=1000*(geantGoodParticle->Momentum(i_s).E()-geantGoodParticle->Mass());
	if(geantGoodParticle->Position(i_s).Z()<0) fbeamtrackEnergy=prim_energy; //correct energy at the beginning of track
      } //loop over beam trks

      //new section 
      art::ServiceHandle<cheat::BackTrackerService> bt_serv;


      // art::ServiceHandle<geo::Geometry> geom;
      simb::MCTrajectory truetraj=geantGoodParticle->Trajectory();
      auto thisTrajectoryProcessMap1 =  truetraj.TrajectoryProcesses();
      if (thisTrajectoryProcessMap1.size()){
	for(auto const& couple: thisTrajectoryProcessMap1){
	  // int_label=truetraj.KeyToProcess(couple.second);
	  fprimary_truth_EndPosition[0]=((truetraj.at(couple.first)).first).X();
	  fprimary_truth_EndPosition[1]=((truetraj.at(couple.first)).first).Y();
	  fprimary_truth_EndPosition[2]=((truetraj.at(couple.first)).first).Z();
	  fprimary_truth_EndProcess=truetraj.KeyToProcess(couple.second);
	  fprimary_truth_Momentum[0]=((truetraj.at(couple.first)).second).X();
	  fprimary_truth_Momentum[1]= ((truetraj.at(couple.first)).second).Y();
	  fprimary_truth_Momentum[2]=((truetraj.at(couple.first)).second).Z();
	  break;
	}
      }

      ////saving the complete information of all the interactions
      std::cout<<"interaction map size "<<thisTrajectoryProcessMap1.size()<<std::endl;
      if (thisTrajectoryProcessMap1.size()){
	for(auto const& couple1: thisTrajectoryProcessMap1){
	 
	  if ((truetraj.KeyToProcess(couple1.second)).find("CoulombScat")!= std::string::npos) continue;
	  // Let's check if the interaction is in the the TPC
	  auto     interactionPos4D =  (truetraj.at(couple1.first)).first ;        
	  if      (interactionPos4D.Z() <  minZ || interactionPos4D.Z() > maxZ ) continue;
	  else if (interactionPos4D.X() <  minX || interactionPos4D.X() > maxX ) continue;
	  else if (interactionPos4D.Y() <  minY || interactionPos4D.Y() > maxY ) continue;
	  interactionX.push_back(((truetraj.at(couple1.first)).first).X());
	  interactionY.push_back(((truetraj.at(couple1.first)).first).Y());
	  interactionZ.push_back(((truetraj.at(couple1.first)).first).Z());
	  double xval=((truetraj.at(couple1.first)).first).X();
	  double zval=((truetraj.at(couple1.first)).first).Z();
	  unsigned int tpcno=1;
	  if(xval<=0 && zval<232) tpcno=1;
	  if(xval<=0 && zval>232 && zval<464) tpcno=5; 
	  if(xval<=0 && zval>=464) tpcno=9;
	  if(xval>0 && zval<232) tpcno=2; 
	  if(xval>0 && zval>232 && zval<464) tpcno=6; 
	  if(xval>0 && zval>=464) tpcno=10;

	  interactionT.push_back(detprop->ConvertXToTicks(((truetraj.at(couple1.first)).first).X(), 2, tpcno, 0));
	  interactionU.push_back(fGeometry->WireCoordinate(((truetraj.at(couple1.first)).first).Y(), ((truetraj.at(couple1.first)).first).Z(),0, tpcno, 0));
	  interactionV.push_back(fGeometry->WireCoordinate(((truetraj.at(couple1.first)).first).Y(), ((truetraj.at(couple1.first)).first).Z(),1, tpcno, 0));
	  interactionW.push_back(fGeometry->WireCoordinate(((truetraj.at(couple1.first)).first).Y(), ((truetraj.at(couple1.first)).first).Z(),2, tpcno, 0));
	  interactionProcesslist.push_back(truetraj.KeyToProcess(couple1.second));
	  std::cout<<"number of interactions "<<thisTrajectoryProcessMap1.size()<<std::endl;
	  std::cout<<"int X, Y, Z and process "<<((truetraj.at(couple1.first)).first).X()<<" "<<((truetraj.at(couple1.first)).first).Y()<<" "<<((truetraj.at(couple1.first)).first).Z()<<" "<<truetraj.KeyToProcess(couple1.second)<<std::endl;
	  ///get the interaction angle here
	  double interactionAngle = 999999.; // This needs to be changed
	  double ut=999999;
	  double vt=999999;
	  double zt=999999;
	  //--------------------- Int Angle ---------------------------
	  // Try to retreive the interaction angle
	  auto  prevInteractionPos4D = (truetraj.at(couple1.first-1)).first ;
	  auto  prevInteractionPos3D = prevInteractionPos4D.Vect() ;
	  auto  interactionPos3D     = interactionPos4D.Vect() ;
	  auto  distanceBtwPoint     = interactionPos3D - prevInteractionPos3D;
	  //Let's try to see if the next point exists
	  if (truetraj.size() > couple1.first + 1) {
	    // The particle doesn't die. No need to check for anything else.
	    auto nextInteractionPos4D =  (truetraj.at(couple1.first+1)).first ;
	    auto nextInteractionPos3D =  nextInteractionPos4D.Vect() ;
	    auto distanceBtwPointNext =  nextInteractionPos3D - interactionPos3D;
	    std::cout<<"distance between points values "<<distanceBtwPoint.X()<<" "<<distanceBtwPoint.Z()<<" "<<distanceBtwPointNext.X()<<" "<<distanceBtwPointNext.Z()<<std::endl;

	    interactionAngles.push_back(TMath::ACos(distanceBtwPointNext.Dot(distanceBtwPoint)/(distanceBtwPointNext.Mag()*distanceBtwPoint.Mag() )  ));

	   
	    double u0=0.4669*fGeometry->WireCoordinate(prevInteractionPos3D.Y(),prevInteractionPos3D.Z(),0,1,0); 
	    double v0=0.4669*fGeometry->WireCoordinate(prevInteractionPos3D.Y(),prevInteractionPos3D.Z(),1,1,0);
	    double z0=0.4792*prevInteractionPos3D.Z();
	    double x0=prevInteractionPos3D.X();

	    double u1=0.4669*fGeometry->WireCoordinate(interactionPos3D.Y(),interactionPos3D.Z(),0,1,0); 
	    double v1=0.4669*fGeometry->WireCoordinate(interactionPos3D.Y(),interactionPos3D.Z(),1,1,0);
	    double z1=0.4792*interactionPos3D.Z();
	    double x1=interactionPos3D.Z();

	    double u2=0.4669*fGeometry->WireCoordinate(nextInteractionPos3D.Y(),nextInteractionPos3D.Z(),0,1,0); 
	    double v2=0.4669*fGeometry->WireCoordinate(nextInteractionPos3D.Y(),nextInteractionPos3D.Z(),1,1,0);
	    double z2=0.4792*nextInteractionPos3D.Z();
	    double x2=nextInteractionPos3D.Z();

	    interactionAnglesUT.push_back(angle2d(u0,x0,u1,x1,u2,x2));
	    interactionAnglesVT.push_back(angle2d(v0,x0,v1,x1,v2,x2));
	    interactionAnglesZT.push_back(angle2d(z0,x0,z1,x1,z2,x2));
 

	    continue;
	  }
	  interactionAnglesUT.push_back(ut);
	  interactionAnglesVT.push_back(vt);
	  interactionAnglesZT.push_back(zt);
	  interactionAngles.push_back(interactionAngle);
	}
      }

      geo::View_t view = geom->View(2);
      auto simIDE_prim=bt_serv->TrackIdToSimIDEs_Ps(geantGoodParticle->TrackId(),view);
      std::map<double, sim::IDE> orderedSimIDE;
      for (auto& ide : simIDE_prim) orderedSimIDE[ide->z]= *ide;
      auto inTPCPoint  = truetraj.begin(); 
      auto Momentum0   = inTPCPoint->second;
      auto old_iter = orderedSimIDE.begin();
      double tlen=0.0;
      double xi=0.0;double yi=0.0;double zi=0.0;
      int count=0; 
     
    
      for ( auto iter= orderedSimIDE.begin(); iter!= orderedSimIDE.end(); iter++,old_iter++){
	auto currentIde = iter->second;
	if(currentIde.z<minZ) continue;
	else if (currentIde.x < minX || currentIde.x > maxX ) continue;
	else if (currentIde.y < minY || currentIde.y > maxY ) continue;
	tmp_primtrk_truth_Z.push_back(currentIde.z);
	tmp_primtrk_truth_Eng.push_back(currentIde.energy);
	tmp_primtrk_truth_trkide.push_back(currentIde.trackID);
	double pos_true[3]={currentIde.x,currentIde.y,currentIde.z};
	/*if(currentIde.z>=0 && currentIde.z<=230){
	  geo::WireID wireID = geom->NearestWireID(pos_true, 2);
	  if (!wireID) wireID = geom->Plane(2).ClosestWireID(wireID);
	  tmp_wire.push_back(wireID.Wire);
	  tmp_tpc.push_back(wireID.TPC);
	  }
	  if(currentIde.z<0 && currentIde.z>230){
	  tmp_wire.push_back(-99999);
	  tmp_tpc.push_back(-99999);
	  }*/
	geo::TPCID tpc1 = geom->FindTPCAtPosition(pos_true);
	if(tpc1.isValid){
	  int tpc_no=tpc1.TPC;
	  geo::PlaneID planeID = geo::PlaneID(0, tpc_no, 2);
	  //	geo::WireID wireID = geom->NearestWireID(pos_true, planeID);
	  geo::WireID wireID;
	  try{
	    wireID = geom->NearestWireID(pos_true, planeID);
	  }
	  catch(geo::InvalidWireError const& e) {
	    wireID = e.suggestedWireID(); // pick the closest valid wire
	  }
	  tmp_wire.push_back(wireID.Wire);
	  tmp_tpc.push_back(wireID.TPC);
	}
	if(!tpc1.isValid){
	  tmp_wire.push_back(-9999);
	  tmp_tpc.push_back(-9999);

	}


	if(count==0){
	  fprimary_truth_StartPosition[0] = currentIde.x;
	  fprimary_truth_StartPosition[1] = currentIde.y;
	  fprimary_truth_StartPosition[2] = currentIde.z;
	}
	if(currentIde.trackID>=0){
	  if(count>0){
	    tlen=tlen+TMath::Sqrt(std::pow(currentIde.x-xi,2)+std::pow(currentIde.y-yi,2)+std::pow(currentIde.z-zi,2));
	  }//if count
	  xi=currentIde.x;yi=currentIde.y;zi=currentIde.z;
	  count++;
	}//trackid>0 loop
      }// iter loop
      fprimary_truth_tracklength=tlen;
      primtrk_truth_Z.push_back(tmp_primtrk_truth_Z);
      primtrk_truth_Eng.push_back(tmp_primtrk_truth_Eng);
      primtrk_truth_trkide.push_back(tmp_primtrk_truth_trkide);
      primtrk_truth_Z_wire.push_back(tmp_wire);
      primtrk_truth_Z_tpc.push_back(tmp_tpc);


      tmp_primtrk_truth_Z.clear();
      tmp_primtrk_truth_Eng.clear();
      tmp_primtrk_truth_trkide.clear();
      tmp_wire.clear();
      tmp_tpc.clear();
      //new section 

    }// geantGoodParticle
  }//is not real data loop

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
  std::cout << "All primary pfParticles = " <<  pfpUtil.GetNumberPrimaryPFParticle(evt,fPFParticleTag) << std::endl;

  // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility 
  // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
  // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
  // PFParticles considered as primary
  // std::vector<recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
 


  auto pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  //cluster information
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  art::Handle< std::vector<recob::PFParticle> > PFPListHandle; 
  std::vector<art::Ptr<recob::PFParticle> > pfplist;

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if(evt.getByLabel("hitpdune",hitListHandle)) art::fill_ptr_vector(hitlist, hitListHandle);

  if(evt.getByLabel("pandoraTrack",trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
  if(evt.getByLabel("pandora",PFPListHandle)) art::fill_ptr_vector(pfplist, PFPListHandle);
  
  art::Handle< std::vector<recob::Cluster> > clusterListHandle; // to get information about the hits
  std::vector<art::Ptr<recob::Cluster>> clusterlist;
  if(evt.getByLabel("pandora", clusterListHandle))
    art::fill_ptr_vector(clusterlist, clusterListHandle);

  art::FindManyP<recob::Cluster> fmcp(PFPListHandle,evt,"pandora");
  art::FindManyP<recob::Track> pftrack(PFPListHandle,evt,"pandoraTrack");
  art::FindManyP<recob::Hit> clhit(clusterListHandle,evt,"pandora");

  std::cout<<"number of pfp_particles "<<pfplist.size()<<std::endl;
  std::cout<<" size of pfParticles "<<pfParticles.size()<<std::endl;
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt,"pandoraTrack"); // to associate tracks and hits


  std::cout<<"outside pfparticle loop"<<std::endl;
  // We can now look at these particles
  for(const recob::PFParticle* particle : pfParticles){
    // Pandora's BDT beam-cosmic score
    fprimaryBDTScore = (double)pfpUtil.GetBeamCosmicScore(*particle,evt,fPFParticleTag);
    // NHits associated with this pfParticle
    fprimaryNHits = (pfpUtil.GetPFParticleHits(*particle,evt,fPFParticleTag)).size();
    // of this particle might be more helpful. These return null pointers if not track-like / shower-like
    const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
    const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
    /////new line added here
    // std::vector<art::Ptr<recob::Cluster>> allClusters=fmcp.at(particle);
    // std::cout<<allClusters.size();
   
    std::cout<<"outside this track loop"<<std::endl;
    if(thisTrack != 0x0){
      if(!beam_cuts.IsBeamlike(*thisTrack, evt, "1")) return;
      // Get the true mc particle
      const simb::MCParticle* mcparticle = truthUtil.GetMCParticleFromRecoTrack(*thisTrack, evt, fTrackerTag);
      if(mcparticle!=0x0){
	std::cout<<"ftruth pdg "<<mcparticle->PdgCode()<<std::endl;
	ftruthpdg=mcparticle->PdgCode();
     
	truthid=mcparticle->TrackId();
	fprimary_truth_Pdg= mcparticle->PdgCode();
	std::cout<<"primary_truth_Pdg "<<fprimary_truth_Pdg<<std::endl;
	fprimary_truth_Isbeammatched=0;

	fprimary_truth_TrackId          = mcparticle->TrackId();
	fprimary_truth_StartPosition[3] = mcparticle->T();
	fprimary_truth_EndPosition[3]   = mcparticle->EndT();
	fprimary_truth_P                = mcparticle->P();
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
	if(beamid==truthid) fprimary_truth_Isbeammatched=1;
      }//mcparticle loop





      fisprimarytrack               = 1;
      fisprimaryshower              = 0;
      fprimaryID                    = thisTrack->ID();
      std::cout<<"this Track track ID "<<thisTrack->ID()<<std::endl;
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
      std::vector<double> tmp_primtrk_dqdx;	
      std::vector<double> tmp_primtrk_resrange;	
      std::vector<double> tmp_primtrk_dedx;	
      std::vector<double> tmp_primtrk_hitx;	
      std::vector<double> tmp_primtrk_hity;	
      std::vector<double> tmp_primtrk_hitz;
      std::vector<double> tmp_primtrk_pitch;
      std::vector<int> tmp_zwire;
      std::vector<int> tmp_ztpc;


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

	    double pos_true[3]={primtrk_pos.X(),primtrk_pos.Y(),primtrk_pos.Z()};
	    /* if(primtrk_pos.Z()>=0 && primtrk_pos.Z()<=230){
	       geo::WireID wireID = geom->NearestWireID(pos_true, 2);
	       if (!wireID) wireID = geom->Plane(2).ClosestWireID(wireID);
	       tmp_zwire.push_back(wireID.Wire);
	       tmp_ztpc.push_back(wireID.TPC);
	       }
	       if(primtrk_pos.Z()<0 || primtrk_pos.Z()>230){
	       tmp_zwire.push_back(-99999);
	       tmp_ztpc.push_back(-99999);
	       }
	    */
	    geo::TPCID tpc2 = geom->FindTPCAtPosition(pos_true);
	    if(tpc2.isValid){
	      int tpc_no=tpc2.TPC;
	      geo::PlaneID planeID = geo::PlaneID(0, tpc_no, 2);
	      geo::WireID wireID;
	      //geo::WireID wireID = geom->NearestWireID(pos_true, planeID);
	      try{
		wireID = geom->NearestWireID(pos_true, planeID);
	      }
	      catch(geo::InvalidWireError const& e) {
		wireID = e.suggestedWireID(); // pick the closest valid wire
	      }
	      tmp_zwire.push_back(wireID.Wire);
	      tmp_ztpc.push_back(wireID.TPC);
	    }
	    if(!tpc2.isValid){
	      tmp_zwire.push_back(-9999);
	      tmp_ztpc.push_back(-9999);
	    }
	  } //loop over hits
	} //only collection plane
      }//calovector

      if(tmp_primtrk_dqdx.size()!=0){
	if(tmp_primtrk_hitz[0]>tmp_primtrk_hitz[tmp_primtrk_hitz.size()-1]){
	  std::reverse(tmp_primtrk_hitz.begin(),tmp_primtrk_hitz.end());
	  std::reverse(tmp_primtrk_hity.begin(),tmp_primtrk_hity.end());
	  std::reverse(tmp_primtrk_hitx.begin(),tmp_primtrk_hitx.end());
	  std::reverse(tmp_primtrk_pitch.begin(),tmp_primtrk_pitch.end());
	  std::reverse(tmp_primtrk_dedx.begin(),tmp_primtrk_dedx.end());
	  std::reverse(tmp_primtrk_dqdx.begin(),tmp_primtrk_dqdx.end());
	  std::reverse(tmp_primtrk_resrange.begin(),tmp_primtrk_resrange.end());

	  std::reverse(tmp_ztpc.begin(),tmp_ztpc.end());
	  std::reverse(tmp_zwire.begin(),tmp_zwire.end());
	}
	primtrk_dqdx.push_back(tmp_primtrk_dqdx);
	primtrk_resrange.push_back(tmp_primtrk_resrange);
	primtrk_dedx.push_back(tmp_primtrk_dedx);
	primtrk_hitx.push_back(tmp_primtrk_hitx);
	primtrk_hity.push_back(tmp_primtrk_hity);
	primtrk_hitz.push_back(tmp_primtrk_hitz);
	primtrk_pitch.push_back(tmp_primtrk_pitch);

	primtrk_hitz_wire.push_back(tmp_zwire);
	primtrk_hitz_tpc.push_back(tmp_ztpc);

      }
      tmp_primtrk_dqdx.clear();
      tmp_primtrk_resrange.clear();
      tmp_primtrk_dedx.clear();
      tmp_primtrk_hitx.clear();
      tmp_primtrk_hity.clear();
      tmp_primtrk_hitz.clear();
      tmp_primtrk_pitch.clear();
      tmp_ztpc.clear();
      tmp_zwire.clear();
     
      for(size_t k = 0; k < calovector.size() && k<3; k++){
	fprimaryKineticEnergy[k] = calovector[k].KineticEnergy();
	fprimaryRange[k] = calovector[k].Range();
	//const std::vector< double > & dedxvec = calovector[k].dEdx();
      }

     

      int planenum=999;
      float zpos=-999;

      //////////Section for looking into integral dQ vs hitz position//////////// 
      //hits and calorimetry loop
      if(fmthm.isValid()){
	auto vhit=fmthm.at(fprimaryID);
	auto vmeta=fmthm.data(fprimaryID);
	for (size_t ii = 0; ii<vhit.size(); ++ii){ //loop over all meta data hit
	  bool fBadhit = false;
	  if (vmeta[ii]->Index() == std::numeric_limits<int>::max()){
	    fBadhit = true;
	    //cout<<"fBadHit"<<fBadhit<<endl;
	    continue;
	  }
	  if (vmeta[ii]->Index()>=tracklist[fprimaryID]->NumberTrajectoryPoints()){
	    throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<tracklist[fprimaryID]->NumberTrajectoryPoints()<<" for track index "<<fprimaryID<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
	  }
	  if (!tracklist[fprimaryID]->HasValidPoint(vmeta[ii]->Index())){
	    fBadhit = true;
	    // cout<<"had valid point "<<fBadhit<<endl;
	    continue;
	  }
	
	  auto loc = tracklist[fprimaryID]->LocationAtPoint(vmeta[ii]->Index());
	  // xpos=loc.X();
	  // ypos=loc.Y();
	  zpos=loc.Z();
	  //	cout<<"x, y, z "<<xpos<<"  "<<ypos<<"  "<<zpos<<endl;
	  //	cout<<"BadHit"<<fBadhit<<endl;
	  if (fBadhit) continue; //HY::If BAD hit, skip this hit and go next
	  if (zpos<-100) continue; //hit not on track
	  planenum=vhit[ii]->WireID().Plane;
	  if(planenum==2){
	    // std::cout<<"inside loop"<<std::endl;
	    // std::array<float,3> cnn_out=hitResults.getOutput(vhit[ii]);
	    peakT_2.push_back(vhit[ii]->PeakTime());
	    int_2.push_back(vhit[ii]->Integral());
	    hitz_2.push_back(zpos);
	    // inelscore.push_back(cnn_out[hitResults.getIndex("inel")]);
	    // elscore.push_back(cnn_out[hitResults.getIndex("el")]);
	    // nonescore.push_back(cnn_out[hitResults.getIndex("none")]);

	    //  std::cout<<"inel, el, none  score"<<cnn_out[hitResults.getIndex("inel")]<<" "<<cnn_out[hitResults.getIndex("el")]<<" "<<cnn_out[hitResults.getIndex("none")]<<std::endl;
	    // std::cout<<"peaktime "<<vhit[ii]->PeakTime()<<std::endl;	
	  }//planenum 2
	  if(planenum==1){
	    int_1.push_back(vhit[ii]->Integral());
	    hitz_1.push_back(zpos);	
	  }//planenum 1
	  if(planenum==0){
	    int_0.push_back(vhit[ii]->Integral());
	    hitz_0.push_back(zpos);	
	  }//planenum 0
	}//loop over vhit
      }//fmthm valid
      //hits and calorimetry loop   


      //**************Section for dE/dx correction*****************************************************************************************************************************//
      //**********************************************************************************************************************************************************************//
      std::vector<float> Stw, Endw, Stt, Endt, Stwires, Endwires, Stticks, Endticks, TPCb, TPCcl, primwireb, primtickb, wirebuf, tickbuf, primdqb;
      std::vector<std::vector<float> > primwire;
      std::vector<std::vector<float> > clwire;
      std::vector<std::vector<float> > primtick;
      std::vector<std::vector<float> >  cltick;
      Stw.clear(); Endw.clear(); Stt.clear();  Endt.clear();  Stwires.clear();  Endwires.clear(); Stticks.clear(); Endticks.clear(); TPCb.clear(); TPCcl.clear();
      float den;
      float numw, numt,wire_no,ticks_no;
      for(size_t p1=0;p1<pfplist.size();p1++){
	std::vector<art::Ptr<recob::Track>> trk=pftrack.at(p1);
	std::vector<art::Ptr<recob::Cluster>> allClusters=fmcp.at(p1);
	for(size_t c1=0;c1<allClusters.size();c1++){
	  if(trk.size() && int(trk[0].key())==fprimaryID) std::cout<<"cluster tpc and plane "<<allClusters[c1]->Plane().TPC<<"  "<<allClusters[c1]->Plane().Plane<<" size "<<clhit.at(allClusters[c1].key()).size()<<std::endl;
	  if(allClusters[c1]->Plane().TPC!=1) continue;
	  if(allClusters[c1]->Plane().Plane==0 && trk.size() && int(trk[0].key())==fprimaryID){ 
	    std::vector<art::Ptr<recob::Hit>> allHits=clhit.at(allClusters[c1].key());
	    for(size_t h1=0;h1<allHits.size();h1++){
	      primwireb.push_back(allHits[h1]->WireID().Wire);
	      primtickb.push_back(allHits[h1]->PeakTime());
	      primdqb.push_back(allHits[h1]->Integral());
	    }
	    wireno_0.push_back(primwireb);
	    peakTime_0.push_back(primtickb);
	    dq_0.push_back(primdqb);
	    primwireb.clear();
	    primtickb.clear();
	    primdqb.clear();
	  }
	  if(allClusters[c1]->Plane().Plane==1 && trk.size() && int(trk[0].key())==fprimaryID){ 
	    std::vector<art::Ptr<recob::Hit>> allHits=clhit.at(allClusters[c1].key());
	    for(size_t h1=0;h1<allHits.size();h1++){
	      primwireb.push_back(allHits[h1]->WireID().Wire);
	      primtickb.push_back(allHits[h1]->PeakTime());
	      primdqb.push_back(allHits[h1]->Integral());
	    }
	    wireno_1.push_back(primwireb);
	    peakTime_1.push_back(primtickb);
	    dq_1.push_back(primdqb);
	    primwireb.clear();
	    primtickb.clear();
	    primdqb.clear();
	  }
	  
	  if(allClusters[c1]->Plane().Plane!=2) continue;
	  if(trk.size() && int(trk[0].key())==fprimaryID){
	    std::vector<art::Ptr<recob::Hit>> allHits=clhit.at(allClusters[c1].key());
	    Stw.push_back(allClusters[c1]->StartWire());
	    Endw.push_back(allClusters[c1]->EndWire());
	    Stt.push_back(allClusters[c1]->StartTick());
	    Endt.push_back(allClusters[c1]->EndTick());
	    TPCb.push_back(allClusters[c1]->Plane().TPC);
	    for(size_t h1=0;h1<allHits.size();h1++){
	      primwireb.push_back(allHits[h1]->WireID().Wire);
	      primtickb.push_back(allHits[h1]->PeakTime());
	      primdqb.push_back(allHits[h1]->Integral());
	    }
	    wireno_2.push_back(primwireb);
	    peakTime_2.push_back(primtickb);
	    dq_2.push_back(primdqb);
	    primwire.push_back(primwireb);
	    primtick.push_back(primtickb);
	    primwireb.clear();
	    primtickb.clear();
	    primdqb.clear();
	  }
	  else if(trk.size()){
	    std::vector<art::Ptr<recob::Hit>> allHits=clhit.at(allClusters[c1].key());
	    // std::cout<<"size of hits "<<allHits.size()<<std::endl;
	    Stwires.push_back(allClusters[c1]->StartWire());
	    Endwires.push_back(allClusters[c1]->EndWire());
	    Stticks.push_back(allClusters[c1]->StartTick());
	    Endticks.push_back(allClusters[c1]->EndTick());
	    TPCcl.push_back(allClusters[c1]->Plane().TPC);
	    for(size_t h2=0;h2<allHits.size();h2++){
	      wirebuf.push_back(allHits[h2]->WireID().Wire);
	      tickbuf.push_back(allHits[h2]->PeakTime());
	    }
	    clwire.push_back(wirebuf);
	    cltick.push_back(tickbuf);
	    wirebuf.clear();
	    tickbuf.clear();
	  }
	}
      }
      //solving the equations here
      for(size_t clt=0;clt<Stw.size();clt++){
	for(size_t cl1=0;cl1<Stwires.size();cl1++){
	  if(TPCcl[cl1]!=TPCb[clt]) continue;
	  den=(Stw[clt]-Endw[clt])*(Stticks[cl1]-Endticks[cl1])-(Stt[clt]-Endt[clt])*(Stwires[cl1]-Endwires[cl1]);
	  bool testv=false;
	  if(Stwires[cl1]==Endwires[cl1]){
	    wire_no=Stwires[cl1];
	    ticks_no=Stticks[cl1];
	    testv=true;
	  }
	  if(den==0  && !testv) continue;
	  if(!testv){
	    numw=(Stw[clt]*Endt[clt]-Stt[clt]*Endw[clt])*(Stwires[cl1]-Endwires[cl1])-(Stw[clt]-Endw[clt])*(Stwires[cl1]*Endticks[cl1]-Stticks[cl1]*Endwires[cl1]);
	    numt=(Stw[clt]*Endt[clt]-Stt[clt]*Endw[clt])*(Stticks[cl1]-Endticks[cl1])-(Stt[clt]-Endt[clt])*(Stwires[cl1]*Endticks[cl1]-Stticks[cl1]*Endwires[cl1]);
	    wire_no=numw/den;
	    ticks_no=numt/den;
	  }
	  if(((Stw[clt]-20<wire_no && Endw[clt]+20>wire_no)||(Stw[clt]+20>wire_no && Endw[clt]<wire_no-20))&&((Stt[clt]-100<ticks_no && Endt[clt]+100>ticks_no)||(Stt[clt]+100>ticks_no && Endt[clt]-100<ticks_no)) && ((Stwires[cl1]-20<wire_no && Endwires[cl1]+20>wire_no)||(Stwires[cl1]+20>wire_no && Endwires[cl1]-20<wire_no)) && ((Stticks[cl1]-100<ticks_no && Endticks[cl1]+100>ticks_no)||(Stticks[cl1]+100>ticks_no && Endticks[cl1]-100<ticks_no))){
	    double xyzStart[3];
	    double xyzEnd[3];
	    unsigned int wireno=std::round(wire_no);
	    // geo::WireID wireid(0,TPCb[clt],2,wireno);
	    if(wireno>=0 && wireno<=479){
	      fGeometry->WireEndPoints(0,TPCb[clt],2,wireno, xyzStart, xyzEnd);
	      Zintersection.push_back(xyzStart[2]);
	      timeintersection.push_back(ticks_no);
	      // std::cout<<"intersecting "<<xyzStart[2]<<std::endl;
	    }
	    else{
	      Zintersection.push_back(wireno*0.47);
	      timeintersection.push_back(ticks_no);
	    }
	    //go for the second solution now
	    float pws=wireno; //p=primary w=wire s=small b=big t=tick c=other clusters
	    float pwb=wireno;
	    float pts=ticks_no;
	    float ptb=ticks_no;
	    float smallest=wireno;
	    float biggest=wireno;
	    float small, big;
	    //finding new variables for beam track
	    for(size_t hit1=0;hit1<primwire[clt].size();hit1++){
	      if(primwire[clt][hit1]>=wireno-30 && primwire[clt][hit1]<=wireno+30){
		if(primwire[clt][hit1]>=wireno-30 && primwire[clt][hit1]<=wireno){
		  small=primwire[clt][hit1];
		  if(small<smallest){
		    pws=primwire[clt][hit1];
		    pts=primtick[clt][hit1];
		  }
		  smallest=primwire[clt][hit1];
		}
		if(primwire[clt][hit1]<=wireno+30 && primwire[clt][hit1]>=wireno){
		  big=primwire[clt][hit1];
		  if(big>biggest){
		    pwb=primwire[clt][hit1];
		    ptb=primtick[clt][hit1];
		  }
		  biggest=primwire[clt][hit1];
		}
		  
	      }//if +_30 loop
		
	    }//hit1 loop
	    //finding new variables for remaining clusters
	    float cws=wireno;
	    float cwb=wireno;
	    float cts=ticks_no;
	    float ctb=ticks_no;
	    smallest=wireno;
	    biggest=wireno;
	    float csmall, cbig;
	    for(size_t hit1=0;hit1<clwire[cl1].size();hit1++){
	      if(clwire[cl1][hit1]>=wireno-30 && clwire[cl1][hit1]<=wireno+30){
		if(clwire[cl1][hit1]>=wireno-30 && clwire[cl1][hit1]<=wireno){
		  csmall=clwire[cl1][hit1];
		  if(csmall<smallest){
		    cws=clwire[cl1][hit1];
		    cts=cltick[cl1][hit1];
		  }
		  smallest=clwire[cl1][hit1];
		}
		if(clwire[cl1][hit1]<=wireno+30 && clwire[cl1][hit1]>=wireno){
		  cbig=clwire[cl1][hit1];
		  if(cbig>biggest){
		    cwb=clwire[cl1][hit1];
		    ctb=cltick[cl1][hit1];
		  }
		  biggest=clwire[cl1][hit1];
		}
			
	      }//if +_30 loop
	    }//hit1 loop
	    float ans[2];
	    soln(pws, pwb, cws, cwb, pts, ptb, cts, ctb,ans);
	    // std::cout<<"answer values "<<ans[0]<<"  "<<ans[1]<<std::endl;
	    // std::cout<<"wire number and new start end wire "<<wireno<<"  st w "<<pws<<"  end w "<<pwb<<" time and new times "<<ticks_no<<"  "<<pts<<"  "<<ptb<<" same for intersecting cluster "<<cws<<" "<<cwb<<" "<<cts<<"  "<<ctb<<" hits in a cluster "<<clwire[cl1].size()<<"start wire and ticks "<<Stwires[cl1]<<" "<<Endwires[cl1]<<" "<<Stticks[cl1]<<"  "<<Endticks[cl1]<<std::endl;    
	    int wir=round(ans[0]);
	    double xyzStart1[3];
	    double xyzEnd1[3];

	    /* fGeometry->WireEndPoints(0,2,2,10, xyzStart1, xyzEnd1);//will remove this line
	       stdd::cout<<"Wire start X in the postive direction "<<xyzStart1[0]<<std::endl;
	       fGeometry->WireEndPoints(0,1,2,10, xyzStart1, xyzEnd1);//will remove this line
	       std::cout<<"Wire start X in the negative direction "<<xyzStart1[0]<<std::endl;*/

	    if(((Stw[clt]-5<ans[0] && Endw[clt]+5>ans[0])||(Stw[clt]+5>ans[0] && Endw[clt]<ans[0]-5))&&((Stt[clt]-50<ans[1] && Endt[clt]+50>ans[1])||(Stt[clt]+50>ans[1] && Endt[clt]-50<ans[1])) && ((Stwires[cl1]-5<ans[0] && Endwires[cl1]+5>ans[0])||(Stwires[cl1]+5>ans[0] && Endwires[cl1]-5<ans[0])) && ((Stticks[cl1]-50<ans[1] && Endticks[cl1]+50>ans[1])||(Stticks[cl1]+50>ans[1] && Endticks[cl1]-50<ans[1]))){
	      if(wir>=0 && wir<=479){
		fGeometry->WireEndPoints(0,TPCb[clt],2,wir, xyzStart1, xyzEnd1);
		Zintersection1.push_back(xyzStart1[2]);
		timeintersection1.push_back(ans[1]);
		std::cout<<"wire X position is here  "<<xyzStart1[0]<<std::endl;
		std::cout<<"wire X position End  is here  "<<xyzEnd1[0]<<std::endl;
	      }
	      else{
		Zintersection1.push_back(wire_no*0.47);
		timeintersection1.push_back(ticks_no);
		ans[0]=wire_no;
		ans[1]=ticks_no;
	      }
	    }
	    else{
	      Zintersection1.push_back(99999);
	      timeintersection1.push_back(99999);
	    }
					     
	    //  if((pws<=ans[0] && pwb>=ans[0])&&(cws<=ans[0] && cwb>=ans[0])&&(pts<=ans[1] && ptb>=ans[1])&&(cts<=ans[1] && ctb>=ans[1])){
	    if(((Stw[clt]<=ans[0] && Endw[clt]>=ans[0])||(Stw[clt]>=ans[0] && Endw[clt]<=ans[0]))&&((Stt[clt]<=ans[1] && Endt[clt]>ans[1])||(Stt[clt]>ans[1] && Endt[clt]<ans[1])) && ((Stwires[cl1]<ans[0] && Endwires[cl1]>ans[0])||(Stwires[cl1]>ans[0] && Endwires[cl1]<ans[0])) && ((Stticks[cl1]<ans[1] && Endticks[cl1]>ans[1])||(Stticks[cl1]>ans[1] && Endticks[cl1]<ans[1]))){
	      deltaZint.push_back(0);
	      deltatimeint.push_back(0);
	    }
	    else if((Stw[clt]<=ans[0] && Endw[clt]>=ans[0])||(Stw[clt]>=ans[0] && Endw[clt]<=ans[0])){
	      deltaZint.push_back(abs(std::min(abs(ans[0]-Stwires[cl1]),abs(ans[0]-Endwires[cl1]))));
	      deltatimeint.push_back(abs(std::min(abs(ans[1]-Stticks[cl1]),abs(ans[1]-Endticks[cl1]))));
	    }
	    else if((Stwires[cl1]<=ans[0] && Endwires[cl1]>=ans[0])||(Stwires[cl1]>=ans[0] && Endwires[cl1]<=ans[0])){
	      deltaZint.push_back(abs(std::min(abs(ans[0]-Stw[clt]),abs(ans[0]-Endw[clt]))));
	      deltatimeint.push_back(abs(std::min(abs(ans[1]-Stt[clt]),abs(ans[1]-Endt[clt]))));
	    }
	    else{
	      deltaZint.push_back(99999);
	      deltatimeint.push_back(99999);
	    }

	  }//if solution lies within wire coordinates
	  //	}//if wire coordinates
	}//cl1 loop
      }//clt
      primwire.clear();
      clwire.clear();
      primtick.clear();
      cltick.clear();
      
      ////*****************section for dE/dx correction*******************************////

	
    }//this track
     
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
	const simb::MCParticle* mcdaughterparticle = truthUtil.GetMCParticleFromRecoTrack(*daughterTrack, evt, fTrackerTag);
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
      // std::vector<anab::T0> pfdaughterT0vec = pfpUtil.GetPFParticleT0(*daughterParticle,evt,fPFParticleTag);
      // if(!pfT0vec.empty())
      //	fdaughterT0[fNDAUGHTERS] = pfdaughterT0vec[0].Time();

      fNDAUGHTERS++;

      // Only process NMAXDAUGTHERS
      if(fNDAUGHTERS > NMAXDAUGTHERS) break;

    }

    break;
  }//particle loop pfparticle

  // Fill trees
  if(beamTriggerEvent)
    fPandoraBeam->Fill();

  //fPandoraCosmics->Fill();
}//analyzer

void protoana::mcsXsection::endJob(){

}

void protoana::mcsXsection::FillCosmicsTree(art::Event const & evt, std::string pfParticleTag){

  // To fill

}



void protoana::mcsXsection::Initialise(){
  
  fRun = -999;
  fSubRun = -999;
  fevent = -999;
  fTimeStamp = -999.0;

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
  fbeamtrackMomentum = -999.0;
  fbeamtrackEnergy = 999.0;
  fbeamtrackPdg = -9999;
  fbeamtrackTime = -999.0;
  fbeamtrackID = -999;
  for(int l=0; l < 3; l++){
    fbeamtrackP[l] = -999.0;
    fbeamtrackPos[l] = -999.0;
    fbeamtrackDir[l] = -999.0;
  }

  //NumberBeamTrajectoryPoints=0; 
  beamtrk_x.clear();
  beamtrk_y.clear();
  beamtrk_z.clear();
  beamtrk_Px.clear();
  beamtrk_Py.clear();
  beamtrk_Pz.clear();
  beamtrk_Eng.clear();
  beamtrk_z_wire.clear();
  beamtrk_z_tpc.clear();

 
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

  fprimary_truth_TrackId = -999;
  fprimary_truth_Pdg = -999;
  ftruthpdg=-999;
  fprimary_truth_P = -999.0;
  fprimary_truth_Pt = -999.0;
  fprimary_truth_Mass = -999.0;
  fprimary_truth_Theta = -999.0;
  fprimary_truth_Phi = -999.0;
  fprimary_truth_Process = -999;
  fprimary_truth_Isbeammatched = -999;
  fprimary_truth_NDaughters = -999;
  fprimary_truth_tracklength=-9999;
  fprimary_truth_EndProcess="";
  truth_last_process="";
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
    // fdaughterT0[k] = -999;

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
  int_2.clear();
  hitz_2.clear();
  // inelscore.clear();
  // elscore.clear();
  // nonescore.clear();
  int_1.clear();
  hitz_1.clear();
  int_0.clear();
  hitz_0.clear();
  peakT_2.clear();
  primtrk_dqdx.clear();
  primtrk_resrange.clear();
  primtrk_dedx.clear();
  primtrk_range.clear();
  primtrk_hitx.clear();
  primtrk_hity.clear();
  primtrk_hitz.clear();
  primtrk_hitz_wire.clear();
  primtrk_hitz_tpc.clear();

  primtrk_pitch.clear();
  primtrk_truth_Z.clear();
  primtrk_truth_Z_wire.clear();
  primtrk_truth_Z_tpc.clear();


  primtrk_truth_Eng.clear();
  primtrk_truth_trkide.clear();
  interactionX.clear();
  interactionY.clear();
  interactionZ.clear();
  interactionU.clear();
  interactionV.clear();
  interactionW.clear();
  interactionT.clear();


  interactionProcesslist.clear();
  interactionAngles.clear();
  interactionAnglesUT.clear();
  interactionAnglesVT.clear();
  interactionAnglesZT.clear();
  Zintersection.clear();
  timeintersection.clear();
  Zintersection1.clear();
  timeintersection1.clear();
  deltaZint.clear();
  deltatimeint.clear();
  wireno_0.clear();
  wireno_1.clear();
  wireno_2.clear();
  dq_0.clear();
  dq_1.clear();
  dq_2.clear();
  peakTime_0.clear();
  peakTime_1.clear();
  peakTime_2.clear();
}

DEFINE_ART_MODULE(protoana::mcsXsection)

