///////////////////////////////////////////////////////////////////////
// Class:      PrimaryPionSelection
// Plugin Type: filter (art v3_02_06)
// File:        PrimaryPionSelection_module.cc
//
// Generated at Tue Jul  9 08:37:33 2019 by Jacob Calcutt using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamCuts.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"

#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "TFile.h"
#include "TProfile.h"

class PrimaryPionSelection;


class PrimaryPionSelection : public art::EDFilter {
public:
  explicit PrimaryPionSelection(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PrimaryPionSelection(PrimaryPionSelection const&) = delete;
  PrimaryPionSelection(PrimaryPionSelection&&) = delete;
  PrimaryPionSelection& operator=(PrimaryPionSelection const&) = delete;
  PrimaryPionSelection& operator=(PrimaryPionSelection&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:
  
  fhicl::ParameterSet beamlineUtil;
  //fhicl::ParameterSet BeamCuts;
  protoana::ProtoDUNEBeamCuts beam_cuts;
  std::string fNominalMomentum;
  std::string fPFParticleTag; 
  std::string fTrackerTag; 
  std::string fShowerTag; 
  std::string fCalorimetryTag;

  fhicl::ParameterSet fCalorimetryParameters;
  std::string fGeneratorTag;
  
  //std::pair< double, double > fTrackStartXCut;
  //std::pair< double, double > fTrackStartYCut;
  //std::pair< double, double > fTrackStartZCut;
  double fTrackEndZCut;
  //double fTrackDirCut;
  bool fStrictNTracks;
  bool fUseMVA;
  double fDaughterCNNCut;
  double fChi2PIDCut;

  std::string dEdX_template_name;
  TFile dEdX_template_file;
  TProfile * profile;

  bool InRange(double, double, double);
};


PrimaryPionSelection::PrimaryPionSelection(fhicl::ParameterSet const& p)
  : EDFilter{p} ,

  beamlineUtil( p.get< fhicl::ParameterSet >("BeamlineUtils")),
  beam_cuts( p.get< fhicl::ParameterSet >("BeamCuts") ),
  fNominalMomentum( p.get< std::string >("NominalMomentum") ),
  fPFParticleTag( p.get< std::string >("PFParticleTag")),
  fTrackerTag( p.get< std::string >("TrackerTag")),
  fShowerTag( p.get< std::string >("ShowerTag")),
  fCalorimetryTag( p.get< std::string >("CalorimetryTag")),
  fCalorimetryParameters( p.get< fhicl::ParameterSet >("CalorimetryParameters")),
  fGeneratorTag( p.get< std::string >("GeneratorTag") ),

  //fTrackStartXCut( p.get< std::pair< double, double> >("TrackStartXCut") ),
  //fTrackStartYCut( p.get< std::pair< double, double> >("TrackStartYCut") ),
  //fTrackStartZCut( p.get< std::pair< double, double> >("TrackStartZCut") ),
  fTrackEndZCut( p.get< double >("TrackEndZCut") ),
 // fTrackDirCut( p.get< double>("TrackDirCut") ),
  fStrictNTracks( p.get< bool >("StrictNTracks") ),
  fUseMVA( p.get< bool >("UseMVA") ),
  fDaughterCNNCut( p.get< double >("DaughterCNNCut") ),
  fChi2PIDCut( p.get< double >("Chi2PIDCut") ),

  dEdX_template_name(p.get<std::string>("dEdX_template_name")),
  dEdX_template_file( dEdX_template_name.c_str(), "OPEN" )

{
  profile = (TProfile*)dEdX_template_file.Get( "dedx_range_pro" );
}

bool PrimaryPionSelection::filter(art::Event& e)
{
  protoana::ProtoDUNEBeamlineUtils    fBeamlineUtils(beamlineUtil);  
  protoana::ProtoDUNEPFParticleUtils  pfpUtil;
  protoana::ProtoDUNETrackUtils       trackUtil;
  
  if( e.isRealData() ){
    if( !fBeamlineUtils.IsGoodBeamlineTrigger( e ) ){
      MF_LOG_INFO("PrimaryPionSelection") << "Failed Beamline Trigger Check" << "\n";
      return false;
    }
  }


  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(e,fPFParticleTag);
  if(beamParticles.size() == 0){
    MF_LOG_INFO("PrimaryPionSelection") << "We found no beam particles for this event... moving on" << "\n";
    return false;
  }

  // Get the reconstructed PFParticle tagged as beam by Pandora
  const recob::PFParticle* particle = beamParticles.at(0);

  // Determine if the beam particle is track-like or shower-like
  const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,e,fPFParticleTag,fTrackerTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,e,fPFParticleTag,fShowerTag);

  if( !thisTrack && thisShower ){
    MF_LOG_INFO("PrimaryPionSelection") << "Beam Particle Reconstructed as shower" << "\n";
    return false;
  }
  else if( !thisShower && !thisTrack ){
    MF_LOG_INFO("PrimaryPionSelection") << "Beam Particle Not Reconstructed" << "\n";
    return false;
  }
  else{
    MF_LOG_INFO("PrimaryPionSelection") << "Beam Particle Reconstructed as track" << "\n";
  }
  ////////////////////////////////////////////////////////////////////
        //The beam Particle has beam Type 13 like in PionAnalyzer module ////

  
  if( !beam_cuts.IsBeamlike( *thisTrack, e, fNominalMomentum ) ){
    MF_LOG_INFO("PrimaryPionSelection") << "Beam Particle failed Beam Cuts" << "\n";
    return false;
  }

  
  //Here add in the cuts for the position of the beam and the incident angle
  //First: need to switch reversed tracks    
  double endZ = thisTrack->Trajectory().End().Z();
  double startZ = thisTrack->Trajectory().Start().Z();
  if( startZ > endZ ){
    double tempZ = endZ;
    endZ = startZ;
    startZ = tempZ;
  }
  //Cut for track length to cut out muons/keep the track within the APA3? 
  if( endZ > fTrackEndZCut ){
    MF_LOG_INFO("PrimaryPionSelection") << "Failed End Z cut" << "\n";
    return false;
  }
  

  //**********IGNORE DAUGHTERS FOR NOW*******
  //Get some objects to use for CNN output checking later
  //anab::MVAReader< recob::Hit, 4 > * hitResults = 0x0;
  ////if( fUseMVA )  hitResults = new anab::MVAReader<recob::Hit,4>(e, "emtrkmichelid:emtrkmichel" );
  //auto recoTracks = e.getValidHandle<std::vector<recob::Track> >(fTrackerTag);
  //art::FindManyP<recob::Hit> findHits(recoTracks,e,fTrackerTag);
  ////////////////////////////
  

  //Look at the daughters and check for track-like daughters that look like showers
  //to try to pick out misreco'd pi0 gammas
  /*
  
  const std::vector< const recob::Track* > trackDaughters = pfpUtil.GetPFParticleDaughterTracks( *particle, e, fPFParticleTag, fTrackerTag );

  for( size_t i = 0; i < trackDaughters.size(); ++i ){
    auto daughterTrack = trackDaughters.at(i);

    
    auto daughterHits = findHits.at( daughterTrack->ID() ); 

    if( fUseMVA ){
      double track_total = 0.;  
      for( size_t h = 0; h < daughterHits.size(); ++h ){
        std::array<float,4> cnn_out = hitResults->getOutput( daughterHits[h] );
        track_total  += cnn_out[ hitResults->getIndex("track") ]; 
      }

      if( track_total < fDaughterCNNCut ){ 
        MF_LOG_INFO("PrimaryPionSelection") << "Found daughter track that looks like shower" << "\n";
        continue;
      }
    }

    //Now: If it's not a potential gamma, pass the calorimetry through the 
    //     Chi2 PID and see if any MIP-like daughters are associated

    auto daughter_calo = trackUtil.GetRecoTrackCalorimetry( *daughterTrack, e, fTrackerTag, fCalorimetryTag );
    std::vector<float> calo_range = daughter_calo[0].ResidualRange();
    std::vector<float> calo_dEdX;
    if( e.isRealData() ) calo_dEdX = trackUtil.CalibrateCalorimetry( *daughterTrack, e, fTrackerTag, fCalorimetryTag, fCalorimetryParameters );
    else calo_dEdX = daughter_calo[0].dEdx();

    std::vector<double> daughter_range, daughter_dEdX;
    for( size_t j = 0; j < calo_range.size(); ++j ){
      daughter_range.push_back( calo_range[i] );
      daughter_dEdX.push_back( calo_dEdX[i] );
    }
    
    std::pair< double,int > chi2_pid_results = trackUtil.Chi2PID( daughter_dEdX, daughter_range, profile );

    if( chi2_pid_results.first > fChi2PIDCut ){
      MF_LOG_INFO("PrimaryPionSelection") << "Found daughter with MIP-like Chi2 PID" << "\n"; 
      return false;
    }
  }

  */
  
    

  return true;
}

bool PrimaryPionSelection::InRange(double input, double low, double high){
  return ( (input >= low) && (input <= high) );
}

DEFINE_ART_MODULE(PrimaryPionSelection)
