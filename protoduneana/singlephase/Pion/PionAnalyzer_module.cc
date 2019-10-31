////////////////////////////////////////////////////////////////////////
// Class:       PionAnalyzer
// Plugin Type: analyzer (art v3_00_00)
// File:        PionAnalyzer_module.cc
//
// Generated at Tue Jan  8 09:12:19 2019 by Jacob Calcutt using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
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

#include "dune/Protodune/singlephase/DataUtils/ProtoDUNETrackUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEShowerUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNETruthUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEBeamlineUtils.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

#include "art_root_io/TFileService.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"

namespace pionana {
  class PionAnalyzer;
}


class pionana::PionAnalyzer : public art::EDAnalyzer {
public:
  explicit PionAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PionAnalyzer(PionAnalyzer const&) = delete;
  PionAnalyzer(PionAnalyzer&&) = delete;
  PionAnalyzer& operator=(PionAnalyzer const&) = delete;
  PionAnalyzer& operator=(PionAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void reset();

private:

  // Declare member data here.
  const art::InputTag fBeamModuleLabel;
  const art::InputTag fTrackModuleLabel;

  TTree *fTree;
  // Run information
  int run;
  int subrun;
  int event;



  double startX, startY, startZ;
  double endX, endY, endZ;
  double len;
  double chi2;
  int chi2_ndof;
  double beam_costheta;
  double new_beam_costheta;
  double beamDirX, beamDirY, beamDirZ;
  double trackDirX, trackDirY, trackDirZ;
  double newDirX, newDirY, newDirZ;
  int beamTrackID;
  std::vector< int > daughter_trackID;
  std::vector< int > daughter_showerID;

  std::vector< double > dEdX, dQdX, resRange;
  std::vector< float > calibrated_dEdX;
  std::vector< std::vector< double > > daughter_dEdX, daughter_dQdX, daughter_resRange;
  double vtxX, vtxY, vtxZ; 
  std::vector< double > daughter_startX, daughter_endX;
  std::vector< double > daughter_startY, daughter_endY;
  std::vector< double > daughter_startZ, daughter_endZ;
  std::vector< double > daughter_shower_startX;
  std::vector< double > daughter_shower_startY;
  std::vector< double > daughter_shower_startZ;
  std::vector< double > daughter_len;

  int nTrackDaughters, nShowerDaughters;

  int type;
  double check_beam_endZ, check_beam_startZ;
  int nBeamParticles;

  std::string fCalorimetryTag;
  
  std::string fTrackerTag;    
  std::string fShowerTag;     
  std::string fPFParticleTag; 
  std::string dEdX_template_name;
  TFile dEdX_template_file;
  bool fVerbose;             
  fhicl::ParameterSet dataUtil;
  fhicl::ParameterSet beamlineUtil;
  double fBrokenTrackZ_low, fBrokenTrackZ_high;
  double fStitchTrackZ_low, fStitchTrackZ_high;
  double fStitchXTol, fStitchYTol;


  fhicl::ParameterSet fCalorimetryParameters;
  fhicl::ParameterSet fBrokenTrackParameters;

  TProfile * proton_template;
};


pionana::PionAnalyzer::PionAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  ,
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),

  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  dEdX_template_name(p.get<std::string>("dEdX_template_name")),
  dEdX_template_file( dEdX_template_name.c_str(), "OPEN" ),
  fVerbose(p.get<bool>("Verbose")),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  beamlineUtil( p.get<fhicl::ParameterSet>("BeamlineUtils")),

  fBrokenTrackZ_low( p.get<double>("BrokenTrackZ_low") ),
  fBrokenTrackZ_high( p.get<double>("BrokenTrackZ_high") ),
  fStitchTrackZ_low( p.get<double>("StitchTrackZ_low") ),
  fStitchTrackZ_high( p.get<double>("StitchTrackZ_high") ),
  fStitchXTol( p.get<double>("StitchXTol") ),
  fStitchYTol( p.get<double>("StitchYTol") ),

  fCalorimetryParameters( p.get< fhicl::ParameterSet > ("CalorimetryParameters") ),
  fBrokenTrackParameters( p.get< fhicl::ParameterSet > ("BrokenTrackParameters") )

{
  proton_template = (TProfile*)dEdX_template_file.Get( "dedx_range_pro" );
}

void pionana::PionAnalyzer::analyze(art::Event const& evt)
{
  //reset containers
  reset();  


  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();

  //Instantiate this here. Fill later if real data
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  
  //Beamline utils 
  protoana::ProtoDUNEBeamlineUtils BLUtil(beamlineUtil);
   
  // Get the PFParticle utility
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  nBeamParticles = beamParticles.size();
  std::cout << "Found " << nBeamParticles << " beamParticles" << std::endl;
  if( nBeamParticles < 1 ) return;

  const recob::PFParticle* particle = beamParticles.at(0);
  const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
  if(thisTrack != 0x0){
    std::cout << "Beam particle is track-like " << thisTrack->ID() << std::endl;
    type = 13;
  }
  if(thisShower != 0x0){
    std::cout << "Beam particle is shower-like" << std::endl;
    type = 11;
    return;
  }
  if(thisTrack == 0x0 && thisShower == 0x0) return;



  // Find the particle vertex. We need the tracker tag here because we need to do a bit of
  // additional work if the PFParticle is track-like to find the vertex. 
  const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
  std::cout << "Vertex: " << vtx[0] << " " << vtx[1] << " " << vtx[2] << std::endl;

  startX = thisTrack->Trajectory().Start().X();
  startY = thisTrack->Trajectory().Start().Y();
  startZ = thisTrack->Trajectory().Start().Z();
  endX = thisTrack->Trajectory().End().X();
  endY = thisTrack->Trajectory().End().Y();
  endZ = thisTrack->Trajectory().End().Z();
  len  = thisTrack->Length();    
  beamTrackID = thisTrack->ID();
  
  std::cout << "Start: " << startX << " " << startY << " " << startZ << std::endl;
  std::cout << "End: " << endX << " " << endY << " " << endZ << std::endl;
  std::cout << "len: " << len << std::endl;

  // Now we can look for the interaction point of the particle if one exists, i.e where the particle
  // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
  // check this by making sure we get a sensible position
  const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
  vtxX = interactionVtx.X();
  vtxY = interactionVtx.Y();
  vtxZ = interactionVtx.Z();
  std::cout << "Interaction Vertex: " << vtxX << " " <<  vtxY << " " <<  vtxZ << std::endl;
  

  // Let's get the daughter PFParticles... we can do this simply without the utility
  for(const int daughterID : particle->Daughters()){
    // Daughter ID is the element of the original recoParticle vector
    const recob::PFParticle *daughterParticle = &(recoParticles->at(daughterID));
    std::cout << "Daughter " << daughterID << " PDG: " << daughterParticle->PdgCode() << std::endl; 
  }
 
  // For actually studying the objects, it is easier to have the daughters in their track and shower forms.
  // We can use the utility to get a vector of track-like and a vector of shower-like daughters
  const std::vector<const recob::Track*> trackDaughters = pfpUtil.GetPFParticleDaughterTracks(*particle,evt,fPFParticleTag,fTrackerTag);  
  const std::vector<const recob::Shower*> showerDaughters = pfpUtil.GetPFParticleDaughterShowers(*particle,evt,fPFParticleTag,fShowerTag);  
  std::cout << "Beam particle has " << trackDaughters.size() << " track-like daughters and " << showerDaughters.size() << " shower-like daughters." << std::endl;
  std::cout << std::endl;

  for( size_t i = 0; i < trackDaughters.size(); ++i ){
    std::cout << "Track daughter " << i << " has len " << trackDaughters[i]->Length() << std::endl; 
    daughter_len.push_back( trackDaughters[i]->Length() );
  }

  for( size_t i = 0; i < showerDaughters.size(); ++i ){
    std::cout << "Shower daughter " << i << " Starts at " << showerDaughters[i]->ShowerStart().X() << " " << showerDaughters[i]->ShowerStart().Y() << " " << showerDaughters[i]->ShowerStart().Z() << std::endl;
    daughter_showerID.push_back( showerDaughters[i]->ID() );
    daughter_shower_startX.push_back( showerDaughters[i]->ShowerStart().X() );
    daughter_shower_startY.push_back( showerDaughters[i]->ShowerStart().Y() );
    daughter_shower_startZ.push_back( showerDaughters[i]->ShowerStart().Z() );
  }

  nTrackDaughters = trackDaughters.size();
  nShowerDaughters = showerDaughters.size();



  //Get Beamline info
  auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamModuleLabel);
  if( beamHandle.isValid()){
    art::fill_ptr_vector(beamVec, beamHandle);
  }
  //Should just have one
  const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0));
  
  std::vector< recob::Track > newTracks = BLUtil.MakeTracks( evt );
  const std::vector< recob::Track > & beamTracks = beamEvent.GetBeamTracks();
  std::cout << "Beamline tracks: " << beamTracks.size() << std::endl;
 
  if( beamTracks.size() == 1 ){
    auto trackDir = thisTrack->StartDirection();
    auto beamDir = beamTracks.at(0).StartDirection();
    
    double flip = 1.;
    if( beamDir.Z() < 0. ) flip = -1.;
    beamDirX = flip * beamDir.X(); 
    beamDirY = flip * beamDir.Y(); 
    beamDirZ = flip * beamDir.Z(); 
 
    if( trackDir.Z() < 0. ) flip = -1.;
    else flip = 1.;
    trackDirX = flip * trackDir.X(); 
    trackDirY = flip * trackDir.Y(); 
    trackDirZ = flip * trackDir.Z(); 
 
    std::cout << "beamDirX: " << beamDir.X() << std::endl;
    std::cout << "beamDirY: " << beamDir.Y() << std::endl;
    std::cout << "beamDirZ: " << beamDir.Z() << std::endl;
    std::cout << "trackDirX: " << trackDir.X() << std::endl;
    std::cout << "trackDirY: " << trackDir.Y() << std::endl;
    std::cout << "trackDirZ: " << trackDir.Z() << std::endl;
 
    beam_costheta = beamDirX*trackDirX + beamDirY*trackDirY + beamDirZ*trackDirZ;
  }
 
  if( newTracks.size() == 1 ){
    auto newDir = newTracks.at(0).StartDirection();
 
    double flip = 1.;
    if( newDir.Z() < 0. ) flip = -1.;
 
    newDirX = flip * newDir.X(); 
    newDirY = flip * newDir.Y(); 
    newDirZ = flip * newDir.Z(); 
 
    std::cout << "newDirX: " << newDir.X() << std::endl;
    std::cout << "newDirY: " << newDir.Y() << std::endl;
    std::cout << "newDirZ: " << newDir.Z() << std::endl;
    new_beam_costheta = newDirX*trackDirX + newDirY*trackDirY + newDirZ*trackDirZ;
  }
 
  //Calorimetry 
  //
  protoana::ProtoDUNETrackUtils trackUtil;
  std::vector< anab::Calorimetry> calo = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
  std::cout << "Planes: " << calo[0].PlaneID().toString() << " " << calo[1].PlaneID().toString()  << " " << calo[2].PlaneID().toString() << std::endl;
  auto calo_dQdX = calo[0].dQdx();
  auto calo_dEdX = calo[0].dEdx();
  auto calo_range = calo[0].ResidualRange();
  for( size_t i = 0; i < calo_dQdX.size(); ++i ){
 
    std::cout << calo_dQdX[i] << " " << calo_dEdX[i] << " " << calo_range[i] << std::endl;
    dQdX.push_back( calo_dQdX[i] );
    dEdX.push_back( calo_dEdX[i] );
    resRange.push_back( calo_range[i] );
  }
 
  calibrated_dEdX = trackUtil.CalibrateCalorimetry( *thisTrack, evt, fTrackerTag, fCalorimetryTag, fCalorimetryParameters);

  //Now get the chi2 from this track
  std::pair< double, int > this_chi2_ndof = trackUtil.Chi2PID( std::vector<double>(calibrated_dEdX.begin(), calibrated_dEdX.end()), resRange, proton_template );
  chi2 = this_chi2_ndof.first;
  chi2_ndof = this_chi2_ndof.second;


  //Go through the track-like daughters and save their calorimetry
  for( size_t i = 0; i < trackDaughters.size(); ++i ){
    auto daughterTrack = trackDaughters.at(i);
    
    daughter_startX.push_back( daughterTrack->Trajectory().Start().X() );
    daughter_startY.push_back( daughterTrack->Trajectory().Start().Y() );
    daughter_startZ.push_back( daughterTrack->Trajectory().Start().Z() );
    daughter_endX.push_back( daughterTrack->Trajectory().End().X() );
    daughter_endY.push_back( daughterTrack->Trajectory().End().Y() );
    daughter_endZ.push_back( daughterTrack->Trajectory().End().Z() );

    daughter_trackID.push_back( daughterTrack->ID() );


    std::vector< anab::Calorimetry > dummy_calo = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackerTag, fCalorimetryTag);
    auto dummy_dQdx = dummy_calo[0].dQdx();
    auto dummy_dEdx = dummy_calo[0].dEdx();
    auto dummy_Range = dummy_calo[0].ResidualRange();
 
    daughter_dQdX.push_back( std::vector<double>() );   
    daughter_resRange.push_back( std::vector<double>() );

    for( size_t i = 0; i < dummy_dQdx.size(); ++i ){
      daughter_dQdX.back().push_back( dummy_dQdx[i] );
      daughter_resRange.back().push_back( dummy_Range[i] );
    }
    std::vector<float> cal_daughter_dEdX = trackUtil.CalibrateCalorimetry( *daughterTrack, evt, fTrackerTag, fCalorimetryTag, fCalorimetryParameters);
    daughter_dEdX.push_back( std::vector<double>(cal_daughter_dEdX.begin(), cal_daughter_dEdX.end() ) );
  }


  fTree->Fill();
}

void pionana::PionAnalyzer::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("beamana","beam analysis tree");
  fTree->Branch("startX", &startX);
  fTree->Branch("startY", &startY);
  fTree->Branch("startZ", &startZ);
  fTree->Branch("endX", &endX);
  fTree->Branch("endY", &endY);
  fTree->Branch("endZ", &endZ);
  fTree->Branch("len", &len);
  fTree->Branch("run", &run);
  fTree->Branch("event", &event);
  fTree->Branch("type", &type);
  fTree->Branch("check_beam_startZ", &check_beam_startZ);
  fTree->Branch("check_beam_endZ", &check_beam_endZ);
  fTree->Branch("nBeamParticles", &nBeamParticles);

  fTree->Branch("beam_costheta", &beam_costheta);
  fTree->Branch("chi2", &chi2);
  fTree->Branch("chi2_ndof", &chi2_ndof);
  fTree->Branch("new_beam_costheta", &new_beam_costheta);
  fTree->Branch("beamDirX", &beamDirX);
  fTree->Branch("beamDirY", &beamDirY);
  fTree->Branch("beamDirZ", &beamDirZ);
  fTree->Branch("trackDirX", &trackDirX);
  fTree->Branch("trackDirY", &trackDirY);
  fTree->Branch("trackDirZ", &trackDirZ);
  fTree->Branch("newDirX", &newDirX);
  fTree->Branch("newDirY", &newDirY);
  fTree->Branch("newDirZ", &newDirZ);
  fTree->Branch("beamTrackID", &beamTrackID);
  fTree->Branch("daughter_trackID", &daughter_trackID);
  fTree->Branch("daughter_showerID", &daughter_showerID);

  fTree->Branch("dQdX", &dQdX);
  fTree->Branch("dEdX", &dEdX);
  fTree->Branch("calibrated_dEdX", &calibrated_dEdX);
  fTree->Branch("resRange", &resRange);
  fTree->Branch("daughter_dQdX", &daughter_dQdX);
  fTree->Branch("daughter_dEdX", &daughter_dEdX);
  fTree->Branch("daughter_resRange", &daughter_resRange);
  fTree->Branch("daughter_len", &daughter_len);
  fTree->Branch("daughter_startX", &daughter_startX);
  fTree->Branch("daughter_startY", &daughter_startY);
  fTree->Branch("daughter_startZ", &daughter_startZ);
  fTree->Branch("vtxX", &vtxX);
  fTree->Branch("vtxY", &vtxY);
  fTree->Branch("vtxZ", &vtxZ);
  fTree->Branch("daughter_endX", &daughter_endX);
  fTree->Branch("daughter_endY", &daughter_endY);
  fTree->Branch("daughter_endZ", &daughter_endZ);
  fTree->Branch("daughter_shower_startX", &daughter_shower_startX);
  fTree->Branch("daughter_shower_startY", &daughter_shower_startY);
  fTree->Branch("daughter_shower_startZ", &daughter_shower_startZ);
  fTree->Branch("nTrackDaughters", &nTrackDaughters);
  fTree->Branch("nShowerDaughters", &nShowerDaughters);
}

void pionana::PionAnalyzer::endJob()
{

}

void pionana::PionAnalyzer::reset()
{
  startX = -1;
  startY = -1;
  startZ = -1;
  endX = -1;
  endY = -1;
  endZ = -1;

  chi2 = 0.;
  chi2_ndof = 0;

  len = -1;
  type = -1;
  check_beam_endZ = 0.;
  check_beam_startZ = 0.;
  nBeamParticles = 0;
  beam_costheta = -100;
  new_beam_costheta = -100;

  nTrackDaughters = -1;
  nShowerDaughters = -1;

  dQdX.clear();
  dEdX.clear();
  calibrated_dEdX.clear();
  vtxX = -1.;
  vtxY = -1.;
  vtxZ = -1.;
  daughter_startX.clear();
  daughter_startY.clear();
  daughter_startZ.clear();
  daughter_endX.clear();
  daughter_endY.clear();
  daughter_endZ.clear();

  daughter_shower_startX.clear();
  daughter_shower_startY.clear();
  daughter_shower_startZ.clear();

  resRange.clear();

  daughter_dQdX.clear();
  daughter_dEdX.clear();
  daughter_resRange.clear();
  daughter_len.clear();

  beamTrackID = -1;
  daughter_trackID.clear();
  daughter_showerID.clear();

}

DEFINE_ART_MODULE(pionana::PionAnalyzer)
