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
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcore/Geometry/Geometry.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamCuts.h"
//#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"

//#include "dune/Protodune/singlephase/DataUtils/ProtoDUNECalibration.h"
#include "protoduneana/Utilities/ProtoDUNECalibration.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

#include "lardata/ArtDataHelper/MVAReader.h"

#include "geant4reweight/src/ReweightBase/G4ReweighterFactory.hh"
#include "geant4reweight/src/ReweightBase/G4Reweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightStep.hh"
#include "geant4reweight/src/PropBase/G4ReweightParameterMaker.hh"


#include "art_root_io/TFileService.h"
#include "TProfile.h"
#include "TFile.h"

// ROOT includes
#include "TTree.h"

namespace pionana {
  class PionAnalyzer;
  bool sort_IDEs( const sim::IDE * i1, const sim::IDE * i2){
    return( i1->z < i2->z ); 
  }

  enum RecoVertexType{
    kUnmatched,
    kInelastic,
    kElastic,
    kBoth,
    kOther
  };

  struct cnnOutput2D{

    cnnOutput2D();

    double track;
    double em;
    double michel;
    double none;
    size_t nHits;
  };

  struct calo_point{

    calo_point();
    calo_point( size_t w, double p, double dedx ) :
      wire(w), pitch(p), dEdX(dedx) {};

    size_t wire;
    double pitch;
    double dEdX;
  };

  cnnOutput2D GetCNNOutputFromPFParticle( const recob::PFParticle & part, const art::Event & evt, const anab::MVAReader<recob::Hit,4> & CNN_results,  protoana::ProtoDUNEPFParticleUtils & pfpUtil, std::string fPFParticleTag ){

    cnnOutput2D output;
    const std::vector< art::Ptr< recob::Hit > > daughterPFP_hits = pfpUtil.GetPFParticleHits_Ptrs( part, evt, fPFParticleTag );    

    for( size_t h = 0; h < daughterPFP_hits.size(); ++h ){
      std::array<float,4> cnn_out = CNN_results.getOutput( daughterPFP_hits[h] );
      output.track  += cnn_out[ CNN_results.getIndex("track") ];
      output.em     += cnn_out[ CNN_results.getIndex("em") ];
      output.michel += cnn_out[ CNN_results.getIndex("michel") ];
      output.none   += cnn_out[ CNN_results.getIndex("none") ];
    }

    output.nHits = daughterPFP_hits.size();

    return output;
  }


  cnnOutput2D GetCNNOutputFromPFParticleFromPlane( const recob::PFParticle & part, const art::Event & evt, const anab::MVAReader<recob::Hit,4> & CNN_results,  protoana::ProtoDUNEPFParticleUtils & pfpUtil, std::string fPFParticleTag, size_t planeID ){

    cnnOutput2D output;
    const std::vector< art::Ptr< recob::Hit > > daughterPFP_hits = pfpUtil.GetPFParticleHitsFromPlane_Ptrs( part, evt, fPFParticleTag, planeID );    

    for( size_t h = 0; h < daughterPFP_hits.size(); ++h ){
      std::array<float,4> cnn_out = CNN_results.getOutput( daughterPFP_hits[h] );
      output.track  += cnn_out[ CNN_results.getIndex("track") ];
      output.em     += cnn_out[ CNN_results.getIndex("em") ];
      output.michel += cnn_out[ CNN_results.getIndex("michel") ];
      output.none   += cnn_out[ CNN_results.getIndex("none") ];
    }

    output.nHits = daughterPFP_hits.size();

    return output;
  }
}

pionana::cnnOutput2D::cnnOutput2D() : track(0), em(0), michel(0), none(0), nHits(0) { }

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
  double lateralDist( TVector3 & n, TVector3 & x0, TVector3 & p );

private:

  
  // Declare member data here.
  const art::InputTag fTrackModuleLabel;

  TTree *fTree;
  // Run information
  int run;
  int subrun;
  int event;
  int MC;


  /////////////////////////////////////////////
  //Truth level info of the primary beam particle
  //that generated the event
  int true_beam_PDG;
  int true_beam_ID;
  std::string true_beam_endProcess;
  double true_beam_endX;
  double true_beam_endY;
  double true_beam_endZ;
  double true_beam_startX;
  double true_beam_startY;
  double true_beam_startZ;

  double true_beam_startDirX;
  double true_beam_startDirY;
  double true_beam_startDirZ;

  double true_beam_startPx;
  double true_beam_startPy;
  double true_beam_startPz;
  double true_beam_startP;

  double true_beam_endPx;
  double true_beam_endPy;
  double true_beam_endPz;
  double true_beam_endP;

  int  true_beam_nElasticScatters;
  int  true_beam_nHits;
  std::vector< double > true_beam_elastic_costheta, true_beam_elastic_X, true_beam_elastic_Y, true_beam_elastic_Z;
  double true_beam_IDE_totalDep;
  bool true_beam_IDE_found_in_recoVtx;
  std::vector< std::string > true_beam_processes;

  std::vector< std::vector< int > > true_beam_reco_byHits_PFP_ID, true_beam_reco_byHits_PFP_nHits,
                                    true_beam_reco_byHits_allTrack_ID;
  //////////////////////////////////////////////////////
  
  //////////////////////////////////////////////////////
  //Truth level info of the daughter MCParticles coming out of the 
  //true primary particle
  std::vector< int > true_beam_daughter_PDG;
  std::vector< int > true_beam_daughter_ID;
  std::vector< double > true_beam_daughter_len;
  std::vector< std::string > true_beam_daughter_Process, true_beam_daughter_endProcess;

  std::vector< double > true_beam_daughter_startX, true_beam_daughter_startY, true_beam_daughter_startZ;
  std::vector< double > true_beam_daughter_startP, true_beam_daughter_startPx, true_beam_daughter_startPy, true_beam_daughter_startPz;
  std::vector< double > true_beam_daughter_endX, true_beam_daughter_endY, true_beam_daughter_endZ;
  std::vector< int >    true_beam_daughter_nHits;


  //going from true to reco byHits
  std::vector< std::vector< int > > true_beam_daughter_reco_byHits_PFP_ID, true_beam_daughter_reco_byHits_PFP_nHits, 
                                    true_beam_daughter_reco_byHits_allTrack_ID, true_beam_daughter_reco_byHits_allShower_ID;
  std::vector< std::vector< double > > true_beam_daughter_reco_byHits_PFP_trackScore;                                    
  std::vector< std::vector< double > > true_beam_daughter_reco_byHits_allTrack_startX, true_beam_daughter_reco_byHits_allTrack_startY, true_beam_daughter_reco_byHits_allTrack_startZ;
  std::vector< std::vector< double > > true_beam_daughter_reco_byHits_allTrack_endX, true_beam_daughter_reco_byHits_allTrack_endY, true_beam_daughter_reco_byHits_allTrack_endZ;
  std::vector< std::vector< double > > true_beam_daughter_reco_byHits_allTrack_len;
  std::vector< std::vector< double > > true_beam_daughter_reco_byHits_allShower_startX, true_beam_daughter_reco_byHits_allShower_startY, true_beam_daughter_reco_byHits_allShower_startZ;
  std::vector< std::vector< double > > true_beam_daughter_reco_byHits_allShower_len;
  //////////////////////////////////////////////////////


  //Decay products from pi0s
  std::vector< int > true_beam_Pi0_decay_PDG, true_beam_Pi0_decay_ID, true_beam_Pi0_decay_parID;
  std::vector< double > true_beam_Pi0_decay_startP;
  std::vector< int > true_beam_Pi0_decay_nHits;
  std::vector< std::vector< int > > true_beam_Pi0_decay_reco_byHits_PFP_ID, true_beam_Pi0_decay_reco_byHits_PFP_nHits,
                                    true_beam_Pi0_decay_reco_byHits_allTrack_ID, true_beam_Pi0_decay_reco_byHits_allShower_ID;
  std::vector< std::vector< double > > true_beam_Pi0_decay_reco_byHits_PFP_trackScore;                                    
  std::vector< std::vector< double > > true_beam_Pi0_decay_reco_byHits_allTrack_startX, true_beam_Pi0_decay_reco_byHits_allTrack_startY, true_beam_Pi0_decay_reco_byHits_allTrack_startZ;
  std::vector< std::vector< double > > true_beam_Pi0_decay_reco_byHits_allTrack_endX, true_beam_Pi0_decay_reco_byHits_allTrack_endY, true_beam_Pi0_decay_reco_byHits_allTrack_endZ;
  std::vector< std::vector< double > > true_beam_Pi0_decay_reco_byHits_allTrack_len;
  std::vector< std::vector< double > > true_beam_Pi0_decay_reco_byHits_allShower_startX, true_beam_Pi0_decay_reco_byHits_allShower_startY, true_beam_Pi0_decay_reco_byHits_allShower_startZ;
  std::vector< std::vector< double > > true_beam_Pi0_decay_reco_byHits_allShower_len;
  //also reco nhits
  std::vector< double > true_beam_Pi0_decay_len;

  std::vector< int > true_beam_grand_daughter_PDG, true_beam_grand_daughter_ID, true_beam_grand_daughter_parID;
  std::vector< int > true_beam_grand_daughter_nHits;
  std::vector< std::string > true_beam_grand_daughter_Process, true_beam_grand_daughter_endProcess;

  //How many of each true particle came out of the true primary beam particle?
  int true_daughter_nPiPlus, true_daughter_nPiMinus, true_daughter_nPi0;
  int true_daughter_nProton, true_daughter_nNeutron, true_daughter_nNucleus;

  //Matched to vertex/slice?
  //
  int reco_beam_vertex_slice;

  std::vector< std::vector< double > > reco_beam_vertex_dRs;
  std::vector< int > reco_beam_vertex_hits_slices;
  ////////////////////////


  //Reconstructed track info
  //EDIT: STANDARDIZE
  double reco_beam_startX, reco_beam_startY, reco_beam_startZ;
  double reco_beam_endX, reco_beam_endY, reco_beam_endZ;
  double reco_beam_vtxX, reco_beam_vtxY, reco_beam_vtxZ; 
  double reco_beam_len;
  double reco_beam_trackDirX, reco_beam_trackDirY, reco_beam_trackDirZ;
  double reco_beam_trackEndDirX, reco_beam_trackEndDirY, reco_beam_trackEndDirZ;
  std::vector< double > reco_beam_dEdX, reco_beam_dQdX, reco_beam_resRange, reco_beam_TrkPitch;
  std::vector< double > reco_beam_calo_wire, reco_beam_calo_tick;
  std::vector< double > reco_beam_calibrated_dEdX;
  int reco_beam_trackID;
  bool reco_beam_flipped;
  
  //fix
  bool reco_beam_passes_beam_cuts;                       

  int reco_beam_nTrackDaughters, reco_beam_nShowerDaughters;
  int reco_beam_type;
  double reco_beam_Chi2_proton;
  int    reco_beam_Chi2_ndof;

  std::vector< double > reco_beam_cosmic_candidate_upper_hits;
  std::vector< double > reco_beam_cosmic_candidate_lower_hits;
  std::vector< int > reco_beam_cosmic_candidate_ID;
  bool beam_has_cosmic_IDE;
  std::vector< int > cosmic_has_beam_IDE;
  int n_cosmics_with_beam_IDE;
  ////////////////////////

  

  //EDIT: STANDARDIZE
  //EndProcess --> endProcess ? 
  std::string reco_beam_true_byE_endProcess, reco_beam_true_byHits_endProcess; //What process ended the reco beam particle
  std::string reco_beam_true_byE_process, reco_beam_true_byHits_process;    //What process created the reco beam particle
  int reco_beam_true_byE_PDG, reco_beam_true_byHits_PDG; 
  int reco_beam_true_byE_ID, reco_beam_true_byHits_ID;
  bool reco_beam_true_byE_matched, reco_beam_true_byHits_matched; //Does the true particle contributing most to the 
                                           //reconstructed beam track coincide with the actual
                                           //beam particle that generated the event
  int reco_beam_true_byE_origin, reco_beam_true_byHits_origin; //What is the origin of the reconstructed beam track?
  //EDIT: STANDARDIZE
  //End_P --> endP, etc.
  double reco_beam_true_byE_endPx,   reco_beam_true_byHits_endPx;
  double reco_beam_true_byE_endPy,   reco_beam_true_byHits_endPy;
  double reco_beam_true_byE_endPz,   reco_beam_true_byHits_endPz;
  double reco_beam_true_byE_endE,    reco_beam_true_byHits_endE;
  double reco_beam_true_byE_endP,    reco_beam_true_byHits_endP;
                                   
  double reco_beam_true_byE_startPx, reco_beam_true_byHits_startPx;
  double reco_beam_true_byE_startPy, reco_beam_true_byHits_startPy;
  double reco_beam_true_byE_startPz, reco_beam_true_byHits_startPz;
  double reco_beam_true_byE_startE,  reco_beam_true_byHits_startE;
  double reco_beam_true_byE_startP,  reco_beam_true_byHits_startP;
  //also throw in byE
  double reco_beam_true_byHits_purity;                      
  //////////////////////////

  std::vector< double > reco_beam_incidentEnergies;
  double reco_beam_interactingEnergy;
  std::vector< double > true_beam_incidentEnergies;
  double true_beam_interactingEnergy;

  int    reco_beam_PFP_ID;
  int    reco_beam_PFP_nHits;
  double reco_beam_PFP_trackScore;
  double reco_beam_PFP_emScore;
  double reco_beam_PFP_michelScore;
  double reco_beam_PFP_trackScore_collection;
  double reco_beam_PFP_emScore_collection;
  double reco_beam_PFP_michelScore_collection;

  int    reco_beam_allTrack_ID;
  bool   reco_beam_allTrack_beam_cuts, reco_beam_allTrack_flipped;
  double reco_beam_allTrack_len;
  double reco_beam_allTrack_startX, reco_beam_allTrack_startY, reco_beam_allTrack_startZ;
  double reco_beam_allTrack_endX, reco_beam_allTrack_endY, reco_beam_allTrack_endZ;
  double reco_beam_allTrack_trackDirX, reco_beam_allTrack_trackDirY, reco_beam_allTrack_trackDirZ;
  double reco_beam_allTrack_trackEndDirX, reco_beam_allTrack_trackEndDirY, reco_beam_allTrack_trackEndDirZ;
  std::vector< double > reco_beam_allTrack_resRange;
  std::vector< double > reco_beam_allTrack_calibrated_dEdX;
  double reco_beam_allTrack_Chi2_proton;
  int    reco_beam_allTrack_Chi2_ndof;

  /////////////////////////////////////////////////////
  //Info from the BI if using Real Data
  /////////////////////////////////////////////////////
  double data_BI_P;
  std::vector< int > data_BI_PDG_candidates;
  double data_BI_X, data_BI_Y, data_BI_Z;
  double data_BI_dirX, data_BI_dirY, data_BI_dirZ;
  int data_BI_nFibersP1, data_BI_nFibersP2, data_BI_nFibersP3;
  int data_BI_nTracks, data_BI_nMomenta;
  ////////////////////////////////////////////////////



  //EDIT: quality_reco_xxx
  bool quality_reco_view_0_hits_in_TPC5, quality_reco_view_1_hits_in_TPC5, quality_reco_view_2_hits_in_TPC5;
  ///BR-MS
  std::vector< double > quality_reco_view_0_wire, quality_reco_view_0_tick;
  std::vector< double > quality_reco_view_1_wire, quality_reco_view_1_tick;
  std::vector< double > quality_reco_view_2_wire, quality_reco_view_2_tick;
  std::vector< double > quality_reco_view_2_z;
  double quality_reco_view_0_max_segment, quality_reco_view_1_max_segment, quality_reco_view_2_max_segment;
  double quality_reco_view_0_wire_backtrack, quality_reco_view_1_wire_backtrack, quality_reco_view_2_wire_backtrack;

  double quality_reco_max_lateral, quality_reco_max_segment; 
  //////






  //Reco-level info of the reconstructed daughters coming out of the
  //reconstructed beam tracl
  //
  //
  //EDIT: daughter_xxx --> daughter_trk_xxx
  //
  //quality_reco_daughter_trk_byY_completeness...
  /*
  std::vector< double > reco_daughter_true_byE_completeness;

  //EDIT: truth --> true_byY_xxx
  std::vector< int > reco_daughter_true_byE_PDG;
  std::vector< int > reco_daughter_true_byE_ID;
  std::vector< int > reco_daughter_true_byE_origin;
  std::vector< int > reco_daughter_true_byE_parID;
  std::vector< int > reco_daughter_true_byE_parPDG;
  std::vector< std::string > reco_daughter_true_byE_process;
  std::vector< double > reco_daughter_true_byE_purity;

  std::vector< int > reco_daughter_true_byHits_PDG;
  std::vector< int > reco_daughter_true_byHits_ID;
  std::vector< int > reco_daughter_true_byHits_origin;
  std::vector< int > reco_daughter_true_byHits_parID;
  std::vector< int > reco_daughter_true_byHits_parPDG;
  std::vector< std::string > reco_daughter_true_byHits_process;
  std::vector< double > reco_daughter_true_byHits_purity;
  std::vector< size_t > reco_daughter_true_byHits_sharedHits, reco_daughter_true_byHits_emHits;

  std::vector< double > reco_daughter_true_byHits_len;
  std::vector< double > reco_daughter_true_byHits_startX;
  std::vector< double > reco_daughter_true_byHits_startY;
  std::vector< double > reco_daughter_true_byHits_startZ;
  std::vector< double > reco_daughter_true_byHits_endX;
  std::vector< double > reco_daughter_true_byHits_endY;
  std::vector< double > reco_daughter_true_byHits_endZ;
  std::vector< double > reco_daughter_true_byHits_startPx;
  std::vector< double > reco_daughter_true_byHits_startPy;
  std::vector< double > reco_daughter_true_byHits_startPz;
  std::vector< double > reco_daughter_true_byHits_startP;
  std::vector< double > reco_daughter_true_byHits_startE;

  //EDIT: reco_daughter_true_byXXX_isPrimary 
  bool reco_daughter_true_byE_isPrimary;
  */
  /// Add by hits?


  //Alternative Reco values
  //EDIT: track_score --> trkScore, etc.
  std::vector< int > reco_daughter_PFP_ID;
  std::vector< int > reco_daughter_PFP_nHits;
  std::vector< double > reco_daughter_PFP_trackScore;
  std::vector< double > reco_daughter_PFP_emScore;
  std::vector< double > reco_daughter_PFP_michelScore;
  std::vector< double > reco_daughter_PFP_trackScore_collection;
  std::vector< double > reco_daughter_PFP_emScore_collection;
  std::vector< double > reco_daughter_PFP_michelScore_collection;
  


  //EDIT: reco_daughter_PFP_true_byY_XXX
  std::vector< int > reco_daughter_PFP_true_byHits_PDG;
  std::vector< int > reco_daughter_PFP_true_byHits_ID;
  std::vector< int > reco_daughter_PFP_true_byHits_origin;
  std::vector< int > reco_daughter_PFP_true_byHits_parID;
  std::vector< int > reco_daughter_PFP_true_byHits_parPDG;
  std::vector< std::string > reco_daughter_PFP_true_byHits_process;
  std::vector< double > reco_daughter_PFP_true_byHits_purity;///EDIT: quality
  std::vector< size_t > reco_daughter_PFP_true_byHits_sharedHits, reco_daughter_PFP_true_byHits_emHits;

  std::vector< double > reco_daughter_PFP_true_byHits_len;
  std::vector< double > reco_daughter_PFP_true_byHits_startX;
  std::vector< double > reco_daughter_PFP_true_byHits_startY;
  std::vector< double > reco_daughter_PFP_true_byHits_startZ;
  std::vector< double > reco_daughter_PFP_true_byHits_endX;
  std::vector< double > reco_daughter_PFP_true_byHits_endY;
  std::vector< double > reco_daughter_PFP_true_byHits_endZ;

  std::vector< double > reco_daughter_PFP_true_byHits_startPx;
  std::vector< double > reco_daughter_PFP_true_byHits_startPy;
  std::vector< double > reco_daughter_PFP_true_byHits_startPz;
  std::vector< double > reco_daughter_PFP_true_byHits_startE;
  std::vector< double > reco_daughter_PFP_true_byHits_startP;

  std::vector< std::string > reco_daughter_PFP_true_byHits_endProcess;

  //////////////////////////////////////

  //EDIT: reco_daughter_allTrack_XXX
  std::vector< int > reco_daughter_allTrack_ID;
  std::vector< double > reco_daughter_allTrack_Theta;
  std::vector< double > reco_daughter_allTrack_Phi;
  std::vector< std::vector< double > > reco_daughter_allTrack_dQdX, reco_daughter_allTrack_dEdX, reco_daughter_allTrack_resRange;
  std::vector< std::vector< double > > reco_daughter_allTrack_dQdX_SCE, reco_daughter_allTrack_dEdX_SCE, reco_daughter_allTrack_resRange_SCE;
  std::vector< std::vector< double > > reco_daughter_allTrack_calibrated_dEdX, reco_daughter_allTrack_calibrated_dEdX_SCE;
  std::vector< double > reco_daughter_allTrack_Chi2_proton;
  std::vector< int >    reco_daughter_allTrack_Chi2_ndof;

  std::vector< double > reco_daughter_allTrack_startX, reco_daughter_allTrack_endX;
  std::vector< double > reco_daughter_allTrack_startY, reco_daughter_allTrack_endY;
  std::vector< double > reco_daughter_allTrack_startZ, reco_daughter_allTrack_endZ;
  std::vector< double > reco_daughter_allTrack_dR;
  std::vector< double > reco_daughter_allTrack_len;
  std::vector< double > reco_daughter_allTrack_to_vertex;
  //
  
  std::vector< int >    reco_daughter_allShower_ID;
  std::vector< double > reco_daughter_allShower_len;
  std::vector< double > reco_daughter_allShower_startX;
  std::vector< double > reco_daughter_allShower_startY;
  std::vector< double > reco_daughter_allShower_startZ;


  //EDIT: STANDARDIZE
  //
  //EDIT: reco_daughter_show_true_byHits_PDG
  /*
  std::vector< int > reco_daughter_shower_true_byHits_PDG;
  std::vector< int > reco_daughter_shower_true_byHits_ID;
  std::vector< int > reco_daughter_shower_true_byHits_origin;
  std::vector< int > reco_daughter_shower_true_byHits_parID;
  std::vector< int > reco_daughter_shower_true_byHits_parPDG;
  std::vector< std::string > reco_daughter_shower_true_byHits_process;
  std::vector< double > reco_daughter_shower_true_byHits_purity;
  std::vector< double > reco_daughter_shower_true_byHits_startPx;
  std::vector< double > reco_daughter_shower_true_byHits_startPy;
  std::vector< double > reco_daughter_shower_true_byHits_startPz;
  std::vector< double > reco_daughter_shower_true_byHits_startP;
  std::vector< std::string>  reco_daughter_shower_true_byHits_endProcess;
  //EDIT: reco_daughter_show_ID
  std::vector< int > reco_daughter_showerID;
  std::vector< int > reco_daughter_shower_true_byE_PDG;
  std::vector< int > reco_daughter_shower_true_byE_ID;
  std::vector< int > reco_daughter_shower_true_byE_origin;
  std::vector< int > reco_daughter_shower_true_byE_parID;
  std::vector< int > reco_daughter_shower_true_byE_parPDG;
  std::vector< double > reco_daughter_shower_true_byE_startPx;
  std::vector< double > reco_daughter_shower_true_byE_startPy;
  std::vector< double > reco_daughter_shower_true_byE_startPz;
  std::vector< double > reco_daughter_shower_true_byE_startP;
  std::vector< std::string>  reco_daughter_shower_true_byE_endProcess;
  std::vector< std::vector< double > > reco_daughter_dEdX, reco_daughter_dQdX, reco_daughter_resRange;
  */

  
  ///Reconstructed Daughter Info
  //  --- Tracks
  /*
  std::vector< int > reco_daughter_trackID;
  std::vector< double > reco_daughter_startX, reco_daughter_endX;
  std::vector< double > reco_daughter_startY, reco_daughter_endY;
  std::vector< double > reco_daughter_startZ, reco_daughter_endZ;
  std::vector< double > reco_daughter_deltaR;
  std::vector< double > reco_daughter_dR;
  std::vector< double > reco_daughter_to_vertex;
  std::vector< int >    reco_daughter_slice;
  std::vector< double > reco_daughter_len;
  std::vector< double > reco_daughter_trackScore;
  std::vector< double > reco_daughter_emScore;
  std::vector< double > reco_daughter_michelScore;
  std::vector< double > reco_daughter_Chi2_proton;
  std::vector< int >    reco_daughter_Chi2_ndof;

  std::vector< double > reco_daughter_momByRange_proton;
  std::vector< double > reco_daughter_momByRange_muon;
  */
  std::vector< double > reco_daughter_allTrack_momByRange_proton;
  std::vector< double > reco_daughter_allTrack_momByRange_muon;



  ///Reconstructed Daughter Info
  //  --- Showers
  /*
  std::vector< double > reco_daughter_shower_startX;
  std::vector< double > reco_daughter_shower_startY;
  std::vector< double > reco_daughter_shower_startZ;
  std::vector< double > reco_daughter_shower_len;
  std::vector< double > reco_daughter_shower_to_vertex;
  std::vector< std::vector< double > > reco_daughter_shower_dEdX, reco_daughter_shower_dQdX, reco_daughter_shower_resRange;
  std::vector< double > reco_daughter_shower_trackScore;
  std::vector< double > reco_daughter_shower_emScore;
  std::vector< double > reco_daughter_shower_michelScore;
  std::vector< double > reco_daughter_shower_Chi2_proton;
  std::vector< int >    reco_daughter_shower_Chi2_ndof;
  */

  //New hits info
  std::vector< double > reco_beam_spacePts_X, reco_beam_spacePts_Y, reco_beam_spacePts_Z;
  //reco_daughter_(trk/show)_spacePts_(X,Y,Z)
  std::vector< std::vector< double > > reco_daughter_spacePts_X, reco_daughter_spacePts_Y, reco_daughter_spacePts_Z;
  std::vector< std::vector< double > > reco_daughter_shower_spacePts_X, reco_daughter_shower_spacePts_Y, reco_daughter_shower_spacePts_Z;


  //Geant4Reweight stuff
  //G4ReweighterFactory RWFactory;
  //G4Reweighter * theRW;
  //G4ReweightParameterMaker ParMaker;


  ////New section -- mechanical class members
  std::map< int, TProfile* > templates;

  //FCL pars
  std::string fCalorimetryTag;
  std::string fTrackerTag;    
  std::string fHitTag;    
  std::string fShowerTag;     
  std::string fPFParticleTag; 
  std::string fGeneratorTag;
  std::string fBeamModuleLabel;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
  std::string dEdX_template_name;
  TFile dEdX_template_file;
  bool fVerbose;             
  int fNSliceCheck; 
  fhicl::ParameterSet BeamPars;
  fhicl::ParameterSet BeamCuts;
  protoana::ProtoDUNEBeamCuts beam_cuts;
  fhicl::ParameterSet CalibrationPars;
  protoana::ProtoDUNECalibration calibration;
  bool fSaveHits;
  bool fCheckCosmics;
};


pionana::PionAnalyzer::PionAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  ,
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),

  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fHitTag(p.get<std::string>("HitTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fBeamModuleLabel(p.get<std::string>("BeamModuleLabel")),
  fBeamlineUtils(p.get< fhicl::ParameterSet >("BeamlineUtils")),
  dEdX_template_name(p.get<std::string>("dEdX_template_name")),
  dEdX_template_file( dEdX_template_name.c_str(), "OPEN" ),
  fVerbose(p.get<bool>("Verbose")),
  fNSliceCheck( p.get< int >("NSliceCheck") ),
  BeamPars(p.get<fhicl::ParameterSet>("BeamPars")),
  BeamCuts(p.get<fhicl::ParameterSet>("BeamCuts")),
  CalibrationPars(p.get<fhicl::ParameterSet>("CalibrationPars")),
  fSaveHits( p.get<bool>( "SaveHits" ) ),
  fCheckCosmics( p.get<bool>( "CheckCosmics" ) )
{

  templates[ 211 ]  = (TProfile*)dEdX_template_file.Get( "dedx_range_pi"  );
  templates[ 321 ]  = (TProfile*)dEdX_template_file.Get( "dedx_range_ka"  );
  templates[ 13 ]   = (TProfile*)dEdX_template_file.Get( "dedx_range_mu"  );
  templates[ 2212 ] = (TProfile*)dEdX_template_file.Get( "dedx_range_pro" );

  calibration = protoana::ProtoDUNECalibration( CalibrationPars );
  beam_cuts = protoana::ProtoDUNEBeamCuts( BeamCuts );

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pionana::PionAnalyzer::analyze(art::Event const& evt)
{

  //reset containers
  reset();  


  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();

  if( !evt.isRealData() ) MC = 1;
  else MC = 0;

   
  // Get various utilities 
  protoana::ProtoDUNEPFParticleUtils                    pfpUtil;
  auto pfpVec = evt.getValidHandle< std::vector< recob::PFParticle > >( fPFParticleTag );   

  protoana::ProtoDUNETrackUtils                         trackUtil;
  protoana::ProtoDUNEShowerUtils                        showerUtil;
  art::ServiceHandle<cheat::BackTrackerService>         bt_serv;
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  const sim::ParticleList & plist = pi_serv->ParticleList(); 

  protoana::ProtoDUNETruthUtils                         truthUtil;
  art::ServiceHandle < geo::Geometry > fGeometryService;
  trkf::TrackMomentumCalculator track_p_calc;
  ////////////////////////////////////////
  


  // This gets the true beam particle that generated the event
  const simb::MCParticle* true_beam_particle = 0x0;
  if( !evt.isRealData() ){
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    true_beam_particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
    if( !true_beam_particle ){
      std::cout << "No true beam particle" << std::endl;
      return;
    }
  }
  ////////////////////////////
  

  // Getting the BI from the data events
  if( evt.isRealData() ){
    auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamModuleLabel);
    
    std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
    if( beamHandle.isValid()){
      art::fill_ptr_vector(beamVec, beamHandle);
    }

    const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); //Should just have one

    if( !fBeamlineUtils.IsGoodBeamlineTrigger( evt ) ){
      std::cout << "Failed quality check" << std::endl;
      return;
    }

    int nTracks = beamEvent.GetBeamTracks().size();
    std::vector< double > momenta = beamEvent.GetRecoBeamMomenta();
    int nMomenta = momenta.size();

/*
    if( !( nMomenta == 1 && nTracks == 1 ) ){
      std::cout << "Malformed tracks and momenta" << std::endl;
      return;
    }*/
    
    if( nMomenta > 0 )
      data_BI_P = momenta[0];

    if( nTracks > 0 ){
      data_BI_X = beamEvent.GetBeamTracks()[0].Trajectory().End().X();
      data_BI_Y = beamEvent.GetBeamTracks()[0].Trajectory().End().Y();
      data_BI_Z = beamEvent.GetBeamTracks()[0].Trajectory().End().Z();

      data_BI_dirX = beamEvent.GetBeamTracks()[0].Trajectory().EndDirection().X(); 
      data_BI_dirY = beamEvent.GetBeamTracks()[0].Trajectory().EndDirection().Y(); 
      data_BI_dirZ = beamEvent.GetBeamTracks()[0].Trajectory().EndDirection().Z(); 
    }

    data_BI_nTracks = nTracks;
    data_BI_nMomenta = nMomenta;

    std::vector< int > pdg_cands = fBeamlineUtils.GetPID( beamEvent, 1. );
    data_BI_PDG_candidates.insert( data_BI_PDG_candidates.end(), pdg_cands.begin(), pdg_cands.end() );

    data_BI_nFibersP1 = beamEvent.GetActiveFibers( "XBPF022697" ).size();
    data_BI_nFibersP2 = beamEvent.GetActiveFibers( "XBPF022701" ).size();
    data_BI_nFibersP3 = beamEvent.GetActiveFibers( "XBPF022702" ).size();
  }
  else{
    try{
      auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("generator");
      
      std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
      if( beamHandle.isValid()){
        art::fill_ptr_vector(beamVec, beamHandle);
      }

      const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); //Should just have one

      std::cout << "Got beam event" << std::endl;

      int nTracks = beamEvent.GetBeamTracks().size();
      std::cout << "Got " << nTracks << " Tracks" << std::endl;
      std::vector< double > momenta = beamEvent.GetRecoBeamMomenta();
      int nMomenta = momenta.size();
      std::cout << "Got " << nMomenta << " Momenta" << std::endl;

      if( nMomenta > 0 ){
        data_BI_P = momenta[0];
        std::cout << "reco P " << data_BI_P << std::endl;
      }

      if( nTracks > 0 ){
        data_BI_X = beamEvent.GetBeamTracks()[0].Trajectory().End().X();
        data_BI_Y = beamEvent.GetBeamTracks()[0].Trajectory().End().Y();
        data_BI_Z = beamEvent.GetBeamTracks()[0].Trajectory().End().Z();

        data_BI_dirX = beamEvent.GetBeamTracks()[0].Trajectory().EndDirection().X(); 
        data_BI_dirY = beamEvent.GetBeamTracks()[0].Trajectory().EndDirection().Y(); 
        data_BI_dirZ = beamEvent.GetBeamTracks()[0].Trajectory().EndDirection().Z(); 
      }

      data_BI_nTracks = nTracks;
      data_BI_nMomenta = nMomenta;

      /*
      std::vector< int > pdg_cands = fBeamlineUtils.GetPID( beamEvent, 1. );
      data_BI_PDG_candidates.insert( data_BI_PDG_candidates.end(), pdg_cands.begin(), pdg_cands.end() );
      */

      data_BI_nFibersP1 = beamEvent.GetActiveFibers( "XBPF022697" ).size();
      data_BI_nFibersP2 = beamEvent.GetActiveFibers( "XBPF022701" ).size();
      data_BI_nFibersP3 = beamEvent.GetActiveFibers( "XBPF022702" ).size();
    }
    catch( const cet::exception &e ){
      std::cout << "BeamEvent generator object not found, moving on" << std::endl;
    }
  }
  ////////////////////////////
  

  // Helper to get hits and the 4 associated CNN outputs
  // CNN Outputs: EM, Track, Michel, Empty
  // outputNames: track, em, none, michel
  anab::MVAReader<recob::Hit,4> hitResults(evt, /*fNNetModuleLabel*/ "emtrkmichelid:emtrkmichel" );

  auto allHits = evt.getValidHandle<std::vector<recob::Hit> >(fHitTag);

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fTrackerTag);
  art::FindManyP<recob::Hit> findHits(recoTracks,evt,fTrackerTag);

  auto recoShowers = evt.getValidHandle< std::vector< recob::Shower > >(fShowerTag);
  art::FindManyP<recob::Hit> findHitsFromShowers(recoShowers,evt,fShowerTag);

  // clang noticed we weren't using this variable
  //auto allPFPs = evt.getValidHandle< std::vector< recob::PFParticle > >( fPFParticleTag );

  std::map< int, std::vector< int > > trueToPFPs = truthUtil.GetMapMCToPFPs_ByHits( evt, fPFParticleTag, fHitTag );
  

  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);

  if(beamParticles.size() == 0){
    std::cerr << "We found no beam particles for this event... moving on" << std::endl;
    return;
  }
  //////////////////////////////////////////////////////////////////


  // Get the reconstructed PFParticle tagged as beam by Pandora
  const recob::PFParticle* particle = beamParticles.at(0);
  const simb::MCParticle* trueParticle = 0x0;

  if( !evt.isRealData() ){
    protoana::MCParticleSharedHits beam_match  = truthUtil.GetMCParticleByHits( *particle, evt, fPFParticleTag, fHitTag ); 
    if( beam_match.particle ){
      //Check that this is the correct true particle
      //if( beam_match.particle->TrackId() == true_beam_particle->TrackId() )
      //  reco_beam_true_byHits_matched = true;
      reco_beam_true_byHits_matched = ( beam_match.particle->TrackId() == true_beam_particle->TrackId() );
      reco_beam_true_byHits_PDG = beam_match.particle->PdgCode();
      reco_beam_true_byHits_ID = beam_match.particle->TrackId();
      std::cout << "Truth ID: " << reco_beam_true_byHits_ID << std::endl;

      reco_beam_true_byHits_process = beam_match.particle->Process();
      reco_beam_true_byHits_endProcess = beam_match.particle->EndProcess();
      reco_beam_true_byHits_origin = pi_serv->TrackIdToMCTruth_P(beam_match.particle->TrackId())->Origin();

      reco_beam_true_byHits_startPx = beam_match.particle->Px();
      reco_beam_true_byHits_startPy = beam_match.particle->Py();
      reco_beam_true_byHits_startPz = beam_match.particle->Pz();
      reco_beam_true_byHits_startP  = sqrt( reco_beam_true_byHits_startPx*reco_beam_true_byHits_startPx 
                                     + reco_beam_true_byHits_startPy*reco_beam_true_byHits_startPy 
                                     + reco_beam_true_byHits_startPz*reco_beam_true_byHits_startPz );
      reco_beam_true_byHits_startE = beam_match.particle->E();

      size_t np = beam_match.particle->NumberTrajectoryPoints();
      if( np > 1 ){
        reco_beam_true_byHits_endPx = beam_match.particle->Px( np - 2 );
        reco_beam_true_byHits_endPy = beam_match.particle->Py( np - 2 );
        reco_beam_true_byHits_endPz = beam_match.particle->Pz( np - 2 );
        reco_beam_true_byHits_endP  = sqrt( reco_beam_true_byHits_endPx*reco_beam_true_byHits_endPx 
                                     + reco_beam_true_byHits_endPy*reco_beam_true_byHits_endPy 
                                     + reco_beam_true_byHits_endPz*reco_beam_true_byHits_endPz );
        reco_beam_true_byHits_endE  = beam_match.particle->E( np - 2 );
      }

      auto list = truthUtil.GetMCParticleListByHits( *particle, evt, fPFParticleTag, fHitTag );
      double total = 0.;
      double matched_hits = 0.;
      for( size_t j = 0; j < list.size(); ++j ){
      //  std::cout << "Contrib " << j << " " << list[j].first->TrackId() << " " << list[j].second << std::endl;
        //std::cout << "Contrib " << j << " " << list[j].particle->TrackId() << " " << list[j].particle->PdgCode() 
        //           << " " << pi_serv->TrackIdToMCTruth_P(list[j].particle->TrackId())->Origin() 
        //           << " " << list[j].nSharedHits << " " << list[j].nSharedDeltaRayHits << std::endl;

        if( list[j].particle == beam_match.particle ){
           matched_hits = list[j].nSharedHits + list[j].nSharedDeltaRayHits;
        }

        total += list[j].nSharedHits + list[j].nSharedDeltaRayHits;
      }

      reco_beam_true_byHits_purity = ( matched_hits / total );

    }
      
    trueParticle = truthUtil.GetMCParticleFromPFParticle(*particle, evt, fPFParticleTag);
    if( trueParticle ){

      //Check that this is the correct true particle
      if( trueParticle->TrackId() == true_beam_particle->TrackId() ){
        reco_beam_true_byE_matched = true;
      }

      reco_beam_true_byE_PDG = trueParticle->PdgCode();
      reco_beam_true_byE_ID = trueParticle->TrackId();
      std::cout << "Truth ID: " << reco_beam_true_byE_ID << std::endl;

      reco_beam_true_byE_process = trueParticle->Process();
      reco_beam_true_byE_endProcess = trueParticle->EndProcess();
      reco_beam_true_byE_origin = pi_serv->TrackIdToMCTruth_P(trueParticle->TrackId())->Origin();

      reco_beam_true_byE_startPx = trueParticle->Px();
      reco_beam_true_byE_startPy = trueParticle->Py();
      reco_beam_true_byE_startPz = trueParticle->Pz();
      reco_beam_true_byE_startP  = sqrt( reco_beam_true_byE_startPx*reco_beam_true_byE_startPx 
                                     + reco_beam_true_byE_startPy*reco_beam_true_byE_startPy 
                                     + reco_beam_true_byE_startPz*reco_beam_true_byE_startPz );
      reco_beam_true_byE_startE = trueParticle->E();

      size_t np = trueParticle->NumberTrajectoryPoints();
      if( np > 1 ){
        reco_beam_true_byE_endPx = trueParticle->Px( np - 2 );
        reco_beam_true_byE_endPy = trueParticle->Py( np - 2 );
        reco_beam_true_byE_endPz = trueParticle->Pz( np - 2 );
        reco_beam_true_byE_endP  = sqrt( reco_beam_true_byE_endPx*reco_beam_true_byE_endPx 
                                     + reco_beam_true_byE_endPy*reco_beam_true_byE_endPy 
                                     + reco_beam_true_byE_endPz*reco_beam_true_byE_endPz );
        reco_beam_true_byE_endE  = trueParticle->E( np - 2 );
      }

    }

    //Some truth information
    true_beam_endProcess = true_beam_particle->EndProcess();
    
    true_beam_PDG         = true_beam_particle->PdgCode();
    true_beam_ID          = true_beam_particle->TrackId();
    std::cout << "Checking true beam par: " << true_beam_particle->TrackId() << " " 
              << pi_serv->TrackIdToMotherParticle_P( true_beam_particle->TrackId() )->TrackId() << std::endl;
    true_beam_endX = true_beam_particle->EndX();
    true_beam_endY = true_beam_particle->EndY();
    true_beam_endZ = true_beam_particle->EndZ();
    true_beam_startX     = true_beam_particle->Position(0).X();
    true_beam_startY     = true_beam_particle->Position(0).Y();
    true_beam_startZ     = true_beam_particle->Position(0).Z();

    true_beam_startPx    = true_beam_particle->Px();
    true_beam_startPy    = true_beam_particle->Py();
    true_beam_startPz    = true_beam_particle->Pz();
    true_beam_startP     = true_beam_particle->P();

    size_t true_np = true_beam_particle->NumberTrajectoryPoints();

    true_beam_endPx    = true_beam_particle->Px(true_np-2);
    true_beam_endPy    = true_beam_particle->Py(true_np-2);
    true_beam_endPz    = true_beam_particle->Pz(true_np-2);
    true_beam_endP     = true_beam_particle->P(true_np-2);

    true_beam_startDirX  = true_beam_startPx / true_beam_startP;
    true_beam_startDirY  = true_beam_startPy / true_beam_startP;
    true_beam_startDirZ  = true_beam_startPz / true_beam_startP;

    true_beam_nHits = truthUtil.GetMCParticleHits( *true_beam_particle, evt, fHitTag ).size();

    true_beam_reco_byHits_PFP_ID.push_back( std::vector< int >() );
    true_beam_reco_byHits_PFP_nHits.push_back( std::vector< int >() );
    true_beam_reco_byHits_allTrack_ID.push_back( std::vector< int >() );
    for( size_t i = 0; i < trueToPFPs[ true_beam_ID ].size(); ++i ){
      true_beam_reco_byHits_PFP_ID.back().push_back( trueToPFPs[ true_beam_ID ][i] );

      const recob::PFParticle * thePFP = &(pfpVec->at( trueToPFPs[ true_beam_ID ][i] ));
      true_beam_reco_byHits_PFP_nHits.back().push_back(
        pfpUtil.GetPFParticleHits_Ptrs( *thePFP, evt, fPFParticleTag ).size()
      );

      const recob::Track* pandora2Track = 0x0;
      
      try{ 
        pandora2Track = pfpUtil.GetPFParticleTrack( *thePFP, evt, fPFParticleTag, "pandora2Track" );
      }
      catch( const cet::exception &e ){
        std::cout << "pandora2Track object not found, moving on" << std::endl;
      }

      if( pandora2Track ){
        true_beam_reco_byHits_allTrack_ID.back().push_back( pandora2Track->ID() );
      }
      else{
        true_beam_reco_byHits_allTrack_ID.back().push_back( -1 );
      }
    }


    for( int i = 0; i < true_beam_particle->NumberDaughters(); ++i ){
      int daughterID = true_beam_particle->Daughter(i);

      std::cout << "Daughter " << i << " ID: " << daughterID << std::endl;
      auto part = plist[ daughterID ];
      int pid = part->PdgCode();
      true_beam_daughter_PDG.push_back(pid);
      true_beam_daughter_ID.push_back( part->TrackId() );      

      std::cout << "Checking true daughter par: " << part->TrackId() << " " << pi_serv->TrackIdToMotherParticle_P( part->TrackId() )->TrackId() << std::endl;
      std::cout << "Mother check: " << part->Mother() << " " << plist[ part->Mother() ]->TrackId() << std::endl;

      true_beam_daughter_len.push_back( part->Trajectory().TotalLength() );

      true_beam_daughter_startX.push_back( part->Position(0).X() );
      true_beam_daughter_startY.push_back( part->Position(0).Y() );
      true_beam_daughter_startZ.push_back( part->Position(0).Z() );

      true_beam_daughter_endX.push_back( part->EndX() );
      true_beam_daughter_endY.push_back( part->EndY() );
      true_beam_daughter_endZ.push_back( part->EndZ() );

      true_beam_daughter_startPx.push_back( part->Px() );
      true_beam_daughter_startPy.push_back( part->Py() );
      true_beam_daughter_startPz.push_back( part->Pz() );
      true_beam_daughter_startP.push_back( part->P() );

      true_beam_daughter_Process.push_back( part->Process() );
      true_beam_daughter_endProcess.push_back( part->EndProcess() );

      std::cout << "Proccess: " << part->Process() << std::endl; 
      std::cout << "PID: " << pid << std::endl;
      std::cout << "Start: " << part->Position(0).X() << " " << part->Position(0).Y() << " " << part->Position(0).Z() << std::endl;
      std::cout << "End: " << part->EndPosition().X() << " " << part->EndPosition().Y() << " " << part->EndPosition().Z() << std::endl;
      std::cout << "Len: " << part->Trajectory().TotalLength() << std::endl;

      if( part->Process().find( "Inelastic" ) != std::string::npos ){
        std::cout << "Inelastic" << std::endl;
        if( pid == 211  ) ++true_daughter_nPiPlus;
        if( pid == -211 ) ++true_daughter_nPiMinus;
        if( pid == 111  ) ++true_daughter_nPi0;
        if( pid == 2212 ) ++true_daughter_nProton;
        if( pid == 2112 ) ++true_daughter_nNeutron;
        if( pid > 2212  ) ++true_daughter_nNucleus; 
      }

      //Look for the gammas coming out of the pi0s
      if( pid == 111 ){
        //std::cout << "Found pi0. Looking at true daughters" << std::endl;
        for( int j = 0; j < part->NumberDaughters(); ++j ){
          int pi0_decay_daughter_ID = part->Daughter(j);
          auto pi0_decay_part = plist[ pi0_decay_daughter_ID ];
          true_beam_Pi0_decay_PDG.push_back( pi0_decay_part->PdgCode() );
          true_beam_Pi0_decay_ID.push_back( pi0_decay_part->TrackId() );
          true_beam_Pi0_decay_startP.push_back( pi0_decay_part->P() );
          true_beam_Pi0_decay_parID.push_back( pi0_decay_part->Mother() );

          true_beam_Pi0_decay_len.push_back( pi0_decay_part->Trajectory().TotalLength() );
          true_beam_Pi0_decay_nHits.push_back( truthUtil.GetMCParticleHits( *pi0_decay_part, evt, fHitTag ).size() );

          true_beam_Pi0_decay_reco_byHits_PFP_ID.push_back( std::vector<int>() );
          true_beam_Pi0_decay_reco_byHits_PFP_nHits.push_back( std::vector<int>() );
          true_beam_Pi0_decay_reco_byHits_PFP_trackScore.push_back( std::vector<double>() );

          true_beam_Pi0_decay_reco_byHits_allTrack_ID.push_back( std::vector<int>() );
          true_beam_Pi0_decay_reco_byHits_allTrack_startX.push_back( std::vector<double>() );
          true_beam_Pi0_decay_reco_byHits_allTrack_startY.push_back( std::vector<double>() );
          true_beam_Pi0_decay_reco_byHits_allTrack_startZ.push_back( std::vector<double>() );
          true_beam_Pi0_decay_reco_byHits_allTrack_len.push_back( std::vector<double>() );
          true_beam_Pi0_decay_reco_byHits_allTrack_endX.push_back( std::vector<double>() );
          true_beam_Pi0_decay_reco_byHits_allTrack_endY.push_back( std::vector<double>() );
          true_beam_Pi0_decay_reco_byHits_allTrack_endZ.push_back( std::vector<double>() );
          true_beam_Pi0_decay_reco_byHits_allShower_ID.push_back( std::vector<int>() );
          true_beam_Pi0_decay_reco_byHits_allShower_startX.push_back( std::vector<double>() );
          true_beam_Pi0_decay_reco_byHits_allShower_startY.push_back( std::vector<double>() );
          true_beam_Pi0_decay_reco_byHits_allShower_startZ.push_back( std::vector<double>() );
          true_beam_Pi0_decay_reco_byHits_allShower_len.push_back( std::vector<double>() );

          for( size_t k = 0; k < trueToPFPs[ pi0_decay_part->TrackId() ].size(); ++k ){
            true_beam_Pi0_decay_reco_byHits_PFP_ID.back().push_back( trueToPFPs[ pi0_decay_part->TrackId() ][k] );

            const recob::PFParticle * thePFP = &(pfpVec->at( trueToPFPs[ pi0_decay_part->TrackId() ][k] ));
            true_beam_Pi0_decay_reco_byHits_PFP_nHits.back().push_back(
              pfpUtil.GetPFParticleHits_Ptrs( *thePFP, evt, fPFParticleTag ).size()
            );

            cnnOutput2D theCNNResults = GetCNNOutputFromPFParticle( *thePFP, evt, hitResults, pfpUtil, fPFParticleTag );
            true_beam_Pi0_decay_reco_byHits_PFP_trackScore.back().push_back( ( ( theCNNResults.nHits > 0 ) ? ( theCNNResults.track / theCNNResults.nHits ) : -999. ) );

            const recob::Track* pandora2Track = 0x0; 
            try{
              pandora2Track = pfpUtil.GetPFParticleTrack( *thePFP, evt, fPFParticleTag, "pandora2Track" );
            }
            catch( const cet::exception &e ){
              std::cout << "pandora2Track object not found, moving on" << std::endl;
            }

            if( pandora2Track ){
              true_beam_Pi0_decay_reco_byHits_allTrack_ID.back().push_back( pandora2Track->ID() );
              true_beam_Pi0_decay_reco_byHits_allTrack_startX.back().push_back( pandora2Track->Trajectory().Start().X() );
              true_beam_Pi0_decay_reco_byHits_allTrack_startY.back().push_back( pandora2Track->Trajectory().Start().Y() );
              true_beam_Pi0_decay_reco_byHits_allTrack_startZ.back().push_back( pandora2Track->Trajectory().Start().Z() );
              true_beam_Pi0_decay_reco_byHits_allTrack_endX.back().push_back( pandora2Track->Trajectory().End().X() );
              true_beam_Pi0_decay_reco_byHits_allTrack_endY.back().push_back( pandora2Track->Trajectory().End().Y() );
              true_beam_Pi0_decay_reco_byHits_allTrack_endZ.back().push_back( pandora2Track->Trajectory().End().Z() );
              true_beam_Pi0_decay_reco_byHits_allTrack_len.back().push_back( pandora2Track->Length() );
            
            }
            else{
              true_beam_Pi0_decay_reco_byHits_allTrack_ID.back().push_back( -1 );
              true_beam_Pi0_decay_reco_byHits_allTrack_startX.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allTrack_startY.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allTrack_startZ.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allTrack_endX.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allTrack_endY.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allTrack_endZ.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allTrack_len.back().push_back( -999. );
            }

            const recob::Shower* pandora2Shower = 0x0;
            try{ 
              pandora2Shower = pfpUtil.GetPFParticleShower( *thePFP, evt, fPFParticleTag, "pandora2Shower" );
            }
            catch( const cet::exception &e ){
              std::cout << "pandora2Shower object not found, moving on" << std::endl;
            }

            if( pandora2Shower ){
              true_beam_Pi0_decay_reco_byHits_allShower_ID.back().push_back( pandora2Shower->ID() );
              true_beam_Pi0_decay_reco_byHits_allShower_startX.back().push_back( pandora2Shower->ShowerStart().X() );
              true_beam_Pi0_decay_reco_byHits_allShower_startY.back().push_back( pandora2Shower->ShowerStart().Y() );
              true_beam_Pi0_decay_reco_byHits_allShower_startZ.back().push_back( pandora2Shower->ShowerStart().Z() );
              true_beam_Pi0_decay_reco_byHits_allShower_len.back().push_back( pandora2Shower->Length() );
            }
            else{
              true_beam_Pi0_decay_reco_byHits_allShower_ID.back().push_back( -1 );
              true_beam_Pi0_decay_reco_byHits_allShower_startX.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allShower_startY.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allShower_startZ.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allShower_len.back().push_back( -999. );
            }

          }

        }
      }

      for( int j = 0; j < part->NumberDaughters(); ++j ){
        int grand_daughter_ID = part->Daughter(j);
        auto grand_daughter_part = plist[ grand_daughter_ID ];
        true_beam_grand_daughter_PDG.push_back( grand_daughter_part->PdgCode() );
        true_beam_grand_daughter_ID.push_back(  grand_daughter_part->TrackId() );
        true_beam_grand_daughter_parID.push_back(  part->TrackId() );
        true_beam_grand_daughter_nHits.push_back( truthUtil.GetMCParticleHits( *grand_daughter_part, evt, fHitTag ).size() );
        true_beam_grand_daughter_Process.push_back( grand_daughter_part->Process() );
        true_beam_grand_daughter_endProcess.push_back( grand_daughter_part->EndProcess() );
      }


      true_beam_daughter_reco_byHits_PFP_ID.push_back( std::vector<int>() );
      true_beam_daughter_reco_byHits_PFP_nHits.push_back( std::vector<int>() );
      true_beam_daughter_reco_byHits_PFP_trackScore.push_back( std::vector<double>() );

      true_beam_daughter_reco_byHits_allTrack_ID.push_back( std::vector<int>() );
      true_beam_daughter_reco_byHits_allTrack_startX.push_back( std::vector<double>() );
      true_beam_daughter_reco_byHits_allTrack_startY.push_back( std::vector<double>() );
      true_beam_daughter_reco_byHits_allTrack_startZ.push_back( std::vector<double>() );
      true_beam_daughter_reco_byHits_allTrack_len.push_back( std::vector<double>() );
      true_beam_daughter_reco_byHits_allTrack_endX.push_back( std::vector<double>() );
      true_beam_daughter_reco_byHits_allTrack_endY.push_back( std::vector<double>() );
      true_beam_daughter_reco_byHits_allTrack_endZ.push_back( std::vector<double>() );
      true_beam_daughter_reco_byHits_allShower_ID.push_back( std::vector<int>() );
      true_beam_daughter_reco_byHits_allShower_startX.push_back( std::vector<double>() );
      true_beam_daughter_reco_byHits_allShower_startY.push_back( std::vector<double>() );
      true_beam_daughter_reco_byHits_allShower_startZ.push_back( std::vector<double>() );
      true_beam_daughter_reco_byHits_allShower_len.push_back( std::vector<double>() );


      for( size_t j = 0; j < trueToPFPs[ part->TrackId() ].size(); ++j ){
        true_beam_daughter_reco_byHits_PFP_ID.back().push_back( trueToPFPs[ part->TrackId() ][j] );

        const recob::PFParticle * thePFP = &(pfpVec->at( trueToPFPs[ part->TrackId() ][j] ));
        true_beam_daughter_reco_byHits_PFP_nHits.back().push_back(
          pfpUtil.GetPFParticleHits_Ptrs( *thePFP, evt, fPFParticleTag ).size()
        );

        cnnOutput2D theCNNResults = GetCNNOutputFromPFParticle( *thePFP, evt, hitResults, pfpUtil, fPFParticleTag );
        true_beam_daughter_reco_byHits_PFP_trackScore.back().push_back( ( ( theCNNResults.nHits > 0 ) ? ( theCNNResults.track / theCNNResults.nHits ) : -999. ) );

        const recob::Track* pandora2Track = 0x0;
        try{
           pandora2Track = pfpUtil.GetPFParticleTrack( *thePFP, evt, fPFParticleTag, "pandora2Track" );
        }
        catch( const cet::exception &e ){
          std::cout << "pandora2Track object not found, moving on" << std::endl;
        }

        if( pandora2Track ){
          true_beam_daughter_reco_byHits_allTrack_ID.back().push_back( pandora2Track->ID() );
          true_beam_daughter_reco_byHits_allTrack_startX.back().push_back( pandora2Track->Trajectory().Start().X() );
          true_beam_daughter_reco_byHits_allTrack_startY.back().push_back( pandora2Track->Trajectory().Start().Y() );
          true_beam_daughter_reco_byHits_allTrack_startZ.back().push_back( pandora2Track->Trajectory().Start().Z() );
          true_beam_daughter_reco_byHits_allTrack_endX.back().push_back( pandora2Track->Trajectory().End().X() );
          true_beam_daughter_reco_byHits_allTrack_endY.back().push_back( pandora2Track->Trajectory().End().Y() );
          true_beam_daughter_reco_byHits_allTrack_endZ.back().push_back( pandora2Track->Trajectory().End().Z() );
          true_beam_daughter_reco_byHits_allTrack_len.back().push_back( pandora2Track->Length() );
        
        }
        else{
          true_beam_daughter_reco_byHits_allTrack_ID.back().push_back( -1 );
          true_beam_daughter_reco_byHits_allTrack_startX.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allTrack_startY.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allTrack_startZ.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allTrack_endX.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allTrack_endY.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allTrack_endZ.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allTrack_len.back().push_back( -999. );
        }

        const recob::Shower* pandora2Shower = 0x0; 
        try{ 
          pandora2Shower = pfpUtil.GetPFParticleShower( *thePFP, evt, fPFParticleTag, "pandora2Shower" );
        }
        catch( const cet::exception &e ){
          std::cout << "pandora2Shower object not found, moving on" << std::endl;
        }
        
        if( pandora2Shower ){
          true_beam_daughter_reco_byHits_allShower_ID.back().push_back( pandora2Shower->ID() );
          true_beam_daughter_reco_byHits_allShower_startX.back().push_back( pandora2Shower->ShowerStart().X() );
          true_beam_daughter_reco_byHits_allShower_startY.back().push_back( pandora2Shower->ShowerStart().Y() );
          true_beam_daughter_reco_byHits_allShower_startZ.back().push_back( pandora2Shower->ShowerStart().Z() );
          true_beam_daughter_reco_byHits_allShower_len.back().push_back( pandora2Shower->Length() );
        
        }
        else{
          true_beam_daughter_reco_byHits_allShower_ID.back().push_back( -1 );
          true_beam_daughter_reco_byHits_allShower_startX.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allShower_startY.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allShower_startZ.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allShower_len.back().push_back( -999. );
        }

      }

      true_beam_daughter_nHits.push_back( truthUtil.GetMCParticleHits( *part, evt, fHitTag ).size() );

    }
  } 

  reco_beam_PFP_ID = particle->Self();
  const std::vector< art::Ptr< recob::Hit > > beamPFP_hits = pfpUtil.GetPFParticleHits_Ptrs( *particle, evt, fPFParticleTag );
  reco_beam_PFP_nHits = beamPFP_hits.size();

  cnnOutput2D cnn = GetCNNOutputFromPFParticle( *particle, evt, hitResults, pfpUtil, fPFParticleTag );
  if( cnn.nHits > 0 ){
    reco_beam_PFP_trackScore = (cnn.track / cnn.nHits);
    reco_beam_PFP_emScore = (cnn.em / cnn.nHits);
    reco_beam_PFP_michelScore = (cnn.michel / cnn.nHits);
  }
  else{
    reco_beam_PFP_trackScore =  -999.;
    reco_beam_PFP_emScore = -999.;
    reco_beam_PFP_michelScore = -999.;
  }

  cnnOutput2D cnn_collection = GetCNNOutputFromPFParticleFromPlane( *particle, evt, hitResults, pfpUtil, fPFParticleTag, 2 );
  if( cnn_collection.nHits > 0 ){
    reco_beam_PFP_trackScore_collection = (cnn_collection.track / cnn_collection.nHits);
    reco_beam_PFP_emScore_collection = (cnn_collection.em / cnn_collection.nHits);
    reco_beam_PFP_michelScore_collection = (cnn_collection.michel / cnn_collection.nHits);
  }
  else{
    reco_beam_PFP_trackScore_collection =  -999.;
    reco_beam_PFP_emScore_collection = -999.;
    reco_beam_PFP_michelScore_collection = -999.;
  }


  

  // Determine if the beam particle is track-like or shower-like
  const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
  if( thisTrack ){
    // Now we can look for the interaction point of the particle if one exists, i.e where the particle
    // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
    // check this by making sure we get a sensible position
    const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    reco_beam_vtxX = interactionVtx.X();
    reco_beam_vtxY = interactionVtx.Y();
    reco_beam_vtxZ = interactionVtx.Z();
    ////////////////////////////////////////////
    

    std::cout << "Beam particle is track-like " << thisTrack->ID() << std::endl;
    reco_beam_type = 13;

    reco_beam_passes_beam_cuts = beam_cuts.IsBeamlike( *thisTrack, evt, "1" );
    std::cout << "Beam Cuts " << reco_beam_passes_beam_cuts << std::endl;


    reco_beam_trackID = thisTrack->ID();

    reco_beam_startX = thisTrack->Trajectory().Start().X();
    reco_beam_startY = thisTrack->Trajectory().Start().Y();
    reco_beam_startZ = thisTrack->Trajectory().Start().Z();
    reco_beam_endX = thisTrack->Trajectory().End().X();
    reco_beam_endY = thisTrack->Trajectory().End().Y();
    reco_beam_endZ = thisTrack->Trajectory().End().Z();

    auto startDir = thisTrack->StartDirection();
    auto endDir   = thisTrack->EndDirection();

    //try flipping
    if( reco_beam_startZ > reco_beam_endZ ){
      reco_beam_flipped = true;
      reco_beam_endX = thisTrack->Trajectory().Start().X();
      reco_beam_endY = thisTrack->Trajectory().Start().Y();
      reco_beam_endZ = thisTrack->Trajectory().Start().Z();
      reco_beam_startX = thisTrack->Trajectory().End().X();
      reco_beam_startY = thisTrack->Trajectory().End().Y();
      reco_beam_startZ = thisTrack->Trajectory().End().Z();
      
      reco_beam_trackDirX =  -1. * endDir.X(); 
      reco_beam_trackDirY =  -1. * endDir.Y(); 
      reco_beam_trackDirZ =  -1. * endDir.Z(); 

      reco_beam_trackEndDirX =  -1. * startDir.X(); 
      reco_beam_trackEndDirY =  -1. * startDir.Y(); 
      reco_beam_trackEndDirZ =  -1. * startDir.Z(); 
    }
    else{
      reco_beam_flipped = false;
      reco_beam_trackDirX    =  startDir.X(); 
      reco_beam_trackDirY    =  startDir.Y(); 
      reco_beam_trackDirZ    =  startDir.Z(); 
      reco_beam_trackEndDirX =  endDir.X(); 
      reco_beam_trackEndDirY =  endDir.Y(); 
      reco_beam_trackEndDirZ =  endDir.Z(); 
    }

    reco_beam_len  = thisTrack->Length();    
    ////////////////////////////////////////////////////////////////

    std::cout << "N Reco Traj Pts: " << thisTrack->NumberTrajectoryPoints() << std::endl;
    
    TVector3 start( reco_beam_startX, reco_beam_startY, reco_beam_startZ );
    TVector3 dir( reco_beam_trackDirX, reco_beam_trackDirY, reco_beam_trackDirZ );
    for( size_t i = 0; i < thisTrack->NumberTrajectoryPoints(); ++i ){
      auto pt = thisTrack->Trajectory().LocationAtPoint(i);
      if( ( pt.X() - -999. ) < 1.e-6 ) continue;

      TVector3 p( pt.X(), pt.Y(), pt.Z() );
      double dist = lateralDist( dir, start, p );

      if( dist > quality_reco_max_lateral ) quality_reco_max_lateral = dist;

      if( i < thisTrack->NumberTrajectoryPoints() - 1 ){
        auto next_pt = thisTrack->Trajectory().LocationAtPoint(i+1);
        if( ( next_pt.X() - -999. ) > 1.e-6 ){

          TVector3 next_p( next_pt.X(), next_pt.Y(), next_pt.Z() );
          double segment = ( next_p - p ).Mag();
          if( segment > quality_reco_max_segment ) quality_reco_max_segment = segment;
        }
      }
    }



///isRealData
    //Thin slice
    //
    std::map< const recob::Hit *, int > hitsToSlices;
    std::map< int, std::vector< const recob::Hit * > > slicesToHits;
    art::ServiceHandle<geo::Geometry> geom;
    
    double z0 = geom->Wire( geo::WireID(0, 1, 2, 0) ).GetCenter().Z();
    std::cout << "Z0: " << z0 << std::endl;
                                  //p, t, c 
    double pitch = geom->WirePitch( 2, 1, 0);
    std::cout << "Pitch: " << pitch << std::endl;

    size_t nWires = geom->Nwires( 2, 1, 0 );
    std::cout << "nWires: " << nWires << std::endl;

    //Looking at the hits in the beam track
    std::map< size_t, const recob::Hit * > trajPtsToHits = trackUtil.GetRecoHitsFromTrajPoints( *thisTrack, evt, fTrackerTag );
    std::cout << "Hits" << std::endl;
    double max_X = 0.;
    double max_Y = 0.;
    double max_Z = 0.;

    std::vector< int > view_0_TPC;
    std::vector< int > view_1_TPC;
    std::vector< int > view_2_TPC;

    for( auto it = trajPtsToHits.begin(); it != trajPtsToHits.end(); ++it ){

      const recob::Hit * theHit = it->second;
      size_t i = it->first;

      double x = thisTrack->Trajectory().LocationAtPoint(i).X();
      double y = thisTrack->Trajectory().LocationAtPoint(i).Y();
      double z = thisTrack->Trajectory().LocationAtPoint(i).Z();

      if( fSaveHits ){
        //saving all hit coordinates for beamtrack
        reco_beam_spacePts_X.push_back(thisTrack->Trajectory().LocationAtPoint(i).X());
        reco_beam_spacePts_Y.push_back(thisTrack->Trajectory().LocationAtPoint(i).Y());
        reco_beam_spacePts_Z.push_back(thisTrack->Trajectory().LocationAtPoint(i).Z());
      }

      int slice = std::floor( ( thisTrack->Trajectory().LocationAtPoint(i).Z() - z0 ) / pitch );
      hitsToSlices[ theHit ] = slice;
      slicesToHits[ slice ].push_back( theHit );
      std::cout << i << " Position: " <<  x << " " << y << " " << z << " Hit Wire: " << theHit->WireID().Wire << " " 
                << theHit->PeakTime()  << " " << theHit->StartTick() << " " << theHit->EndTick() << std::endl;


      if( thisTrack->Trajectory().LocationAtPoint(i).Z() > max_Z ){
        max_Z = thisTrack->Trajectory().LocationAtPoint(i).Z();
        max_X = thisTrack->Trajectory().LocationAtPoint(i).X();
        max_Y = thisTrack->Trajectory().LocationAtPoint(i).Y();
      }


      //std::cout << "View: " << theHit->View() << std::endl;
      switch( theHit->View() ){
        case 0: 
          if( theHit->WireID().TPC == 5 )
            quality_reco_view_0_hits_in_TPC5 = true;
          quality_reco_view_0_wire.push_back( theHit->WireID().Wire );
          quality_reco_view_0_tick.push_back( theHit->PeakTime() );
          view_0_TPC.push_back( theHit->WireID().TPC );
          break;
        case 1:
          if( theHit->WireID().TPC == 5 )
            quality_reco_view_1_hits_in_TPC5 = true;
          quality_reco_view_1_wire.push_back( theHit->WireID().Wire );
          quality_reco_view_1_tick.push_back( theHit->PeakTime() );
          view_1_TPC.push_back( theHit->WireID().TPC );
          break;
        case 2: 
          if( theHit->WireID().TPC == 5 )
            quality_reco_view_2_hits_in_TPC5 = true;
          quality_reco_view_2_wire.push_back( theHit->WireID().Wire );
          quality_reco_view_2_z.push_back( thisTrack->Trajectory().LocationAtPoint(i).Z() );
          quality_reco_view_2_tick.push_back( theHit->PeakTime() );
          view_2_TPC.push_back( theHit->WireID().TPC );
          break;
        default:
          break;
      }
    }

 //   std::cout << "View 0: " << quality_reco_view_0_wire.size() << std::endl;
    for( size_t i = 1; i < quality_reco_view_0_wire.size(); ++i ){
 //     std::cout << i << " " << i-1 << std::endl;
      double segment = sqrt( (quality_reco_view_0_wire[i] - quality_reco_view_0_wire[i-1])*(quality_reco_view_0_wire[i] - quality_reco_view_0_wire[i-1]) 
                           + (quality_reco_view_0_tick[i] - quality_reco_view_0_tick[i-1])*(quality_reco_view_0_tick[i] - quality_reco_view_0_tick[i-1]) );
      if( segment > quality_reco_view_0_max_segment ) quality_reco_view_0_max_segment = segment;                           

      if( quality_reco_view_0_wire[i] < quality_reco_view_0_wire[i-1] && ( view_0_TPC[i] != 5 && view_0_TPC[i-1] != 5)  ){
        quality_reco_view_0_wire_backtrack += (quality_reco_view_0_wire[i-1] - quality_reco_view_0_wire[i]);
      }
    }

 //   std::cout << "View 1: " << quality_reco_view_1_wire.size() << std::endl;
    for( size_t i = 1; i < quality_reco_view_1_wire.size(); ++i ){
 //     std::cout << i << " " << i-1 << std::endl;
      double segment = sqrt( (quality_reco_view_1_wire[i] - quality_reco_view_1_wire[i-1])*(quality_reco_view_1_wire[i] - quality_reco_view_1_wire[i-1]) 
                           + (quality_reco_view_1_tick[i] - quality_reco_view_1_tick[i-1])*(quality_reco_view_1_tick[i] - quality_reco_view_1_tick[i-1]) );
      if( segment > quality_reco_view_1_max_segment ) quality_reco_view_1_max_segment = segment;                           

      if( quality_reco_view_1_wire[i] > quality_reco_view_1_wire[i-1]  && ( view_1_TPC[i] != 5 && view_1_TPC[i-1] != 5)){
        quality_reco_view_1_wire_backtrack += (quality_reco_view_1_wire[i] - quality_reco_view_1_wire[i-1]);
      }
    }

//    std::cout << "View 2: " << quality_reco_view_2_wire.size() << std::endl;
    for( size_t i = 1; i < quality_reco_view_2_wire.size(); ++i ){
      //std::cout << i << " " << i-1 << std::endl;
      double segment = sqrt( (quality_reco_view_2_wire[i] - quality_reco_view_2_wire[i-1])*(quality_reco_view_2_wire[i] - quality_reco_view_2_wire[i-1]) 
                           + (quality_reco_view_2_tick[i] - quality_reco_view_2_tick[i-1])*(quality_reco_view_2_tick[i] - quality_reco_view_2_tick[i-1]) );
      if( segment > quality_reco_view_2_max_segment ) quality_reco_view_2_max_segment = segment;                           

      if( quality_reco_view_2_wire[i] < quality_reco_view_2_wire[i-1]  && ( view_2_TPC[i] != 5 && view_2_TPC[i-1] != 5)){
        quality_reco_view_2_wire_backtrack += (quality_reco_view_2_wire[i-1] - quality_reco_view_2_wire[i]);
      }
    }

    reco_beam_vertex_slice = slicesToHits.rbegin()->first;
    std::cout << "Vertex slice: " << reco_beam_vertex_slice << std::endl;


    //Go through the hits in the last slice, then backtrack to the IDs
    //std::vector< const recob::Hit * > vertex_hits = slicesToHits.rbegin()->second;
    std::vector< const recob::Hit * > vertex_hits;
    int n_slices = 0;
    auto itHits = slicesToHits.rbegin();
    std::cout << "SliceCheck: " << fNSliceCheck << std::endl;

    std::vector< int > temp_hits_slices;

    while( n_slices < fNSliceCheck && itHits != slicesToHits.rend() ){
      
      std::cout << n_slices << std::endl;

      std::vector< const recob::Hit * > temp_hits = itHits->second;
      vertex_hits.insert( vertex_hits.end(), temp_hits.begin(), temp_hits.end() );

      std::vector<int> hits_slices = std::vector<int>(temp_hits.size(), n_slices );
      temp_hits_slices.insert( temp_hits_slices.end(), hits_slices.begin(), hits_slices.end() );

      ++itHits; 
      ++n_slices;
    }
    ////////////////////////////////

    //Primary Track Calorimetry 
    std::vector< anab::Calorimetry> calo = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
    std::cout << "N Calos: " << calo.size() << std::endl;
    auto calo_dQdX = calo[0].dQdx();
    auto calo_dEdX = calo[0].dEdx();
    auto calo_range = calo[0].ResidualRange();
    auto TpIndices = calo[0].TpIndices();
    std::cout << TpIndices.size() << std::endl;
    for( size_t i = 0; i < calo_dQdX.size(); ++i ){
      reco_beam_dQdX.push_back( calo_dQdX[i] );
      reco_beam_dEdX.push_back( calo_dEdX[i] );
      reco_beam_resRange.push_back( calo_range[i] );
      reco_beam_TrkPitch.push_back( calo[0].TrkPitchVec()[i] );

      std::cout << TpIndices[i] << std::endl;
      recob::Hit theHit = (*allHits)[ TpIndices[i] ];
      reco_beam_calo_wire.push_back( theHit.WireID().Wire );
      reco_beam_calo_tick.push_back( theHit.PeakTime() );
      //auto theHit = trajPtsToHits[ TpIndices[i] ];
      //std::cout << i << " Pt: " << TpIndices[i] << " Hit: " << theHit->WireID().Wire << " " << theHit->View() << std::endl;
    }
    ////////////////////////////////////////////

    auto TpIndices1 = calo[1].TpIndices();
    auto TpIndices2 = calo[2].TpIndices();
    
    //New Calibration
    std::vector< float > new_dEdX = calibration.GetCalibratedCalorimetry(  *thisTrack, evt, fTrackerTag, fCalorimetryTag );
    std::cout << "n dEdX: " << reco_beam_dEdX.size() << " " << new_dEdX.size() << std::endl;
    for( size_t i = 0; i < new_dEdX.size(); ++i ){ reco_beam_calibrated_dEdX.push_back( new_dEdX[i] ); }
    ////////////////////////////////////////////

    std::pair< double, int > pid_chi2_ndof = trackUtil.Chi2PID( reco_beam_calibrated_dEdX, reco_beam_resRange, templates[ 2212 ] );
    reco_beam_Chi2_proton = pid_chi2_ndof.first; 
    reco_beam_Chi2_ndof = pid_chi2_ndof.second;
  
    std::cout << "Proton chi2: " << reco_beam_Chi2_proton << std::endl;

    //Doing thin slice
    if( reco_beam_calibrated_dEdX.size() && reco_beam_calibrated_dEdX.size() == reco_beam_TrkPitch.size() && reco_beam_calibrated_dEdX.size() == reco_beam_calo_wire.size() ){
      std::vector< calo_point > reco_beam_calo_points;
      for( size_t i = 0; i < reco_beam_calibrated_dEdX.size(); ++i ){
        reco_beam_calo_points.push_back(
          calo_point( reco_beam_calo_wire[i], reco_beam_TrkPitch[i], reco_beam_calibrated_dEdX[i] )
        );
      }

      //Sort
      std::sort( reco_beam_calo_points.begin(), reco_beam_calo_points.end(), [](calo_point a, calo_point b) {return ( a.wire < b.wire );} ); 


      //Get the initial Energy KE
      double mass = 0.;
      double init_KE = 0.;
      if( evt.isRealData() ){      
        mass = 139.57;

        init_KE =  sqrt( 1.e6*data_BI_P*data_BI_P + mass*mass ) - mass;
      }
      else{
        if( true_beam_PDG == 2212 ) mass = 938.27;
        else if( abs(true_beam_PDG) == 211 ) mass = 139.57;
        else if( abs(true_beam_PDG) == 11 ) mass = .511;
        else if( abs(true_beam_PDG) == 321 ) mass = 321;
        else if( abs(true_beam_PDG) == 13 )  mass = 105.66;      

        init_KE = sqrt( 1.e6 * true_beam_startP*true_beam_startP + mass*mass ) - mass;
      }

      reco_beam_incidentEnergies.push_back( init_KE );
      for( size_t i = 0; i < reco_beam_calo_points.size() - 1; ++i ){ //-1 to not count the last slice
        double this_energy = reco_beam_incidentEnergies.back() - ( reco_beam_calo_points[i].dEdX * reco_beam_calo_points[i].pitch );
        reco_beam_incidentEnergies.push_back( this_energy ); 
      }
      if( reco_beam_incidentEnergies.size() ) reco_beam_interactingEnergy = reco_beam_incidentEnergies.back();


    }

   

    if( !evt.isRealData() ){

      //Go through the true processes within the MCTrajectory
      const simb::MCTrajectory & true_beam_trajectory = true_beam_particle->Trajectory();
      auto true_beam_proc_map = true_beam_trajectory.TrajectoryProcesses();
      std::cout << "Processes: " << std::endl;

      for( auto itProc = true_beam_proc_map.begin(); itProc != true_beam_proc_map.end(); ++itProc ){
        int index = itProc->first;
        std::string process = true_beam_trajectory.KeyToProcess(itProc->second);
        std::cout << index << " " << process << std::endl;

        true_beam_processes.push_back( process );

        if( process == "hadElastic" ){

          ++true_beam_nElasticScatters;

          double process_X = true_beam_trajectory.X( index );
          double process_Y = true_beam_trajectory.Y( index );
          double process_Z = true_beam_trajectory.Z( index );

          double PX      = true_beam_trajectory.Px( index );
          double next_PX = true_beam_trajectory.Px( index + 1 );
          double PY      = true_beam_trajectory.Py( index );
          double next_PY = true_beam_trajectory.Py( index + 1 );
          double PZ      = true_beam_trajectory.Pz( index );
          double next_PZ = true_beam_trajectory.Pz( index + 1 );

          double total_P = sqrt( PX*PX + PY*PY + PZ*PZ );
          double total_next_P = sqrt( next_PX*next_PX + next_PY*next_PY + next_PZ*next_PZ );

          //Get the angle between the direction of this step and the next
          true_beam_elastic_costheta.push_back(
            ( ( PX * next_PX ) + ( PY * next_PY ) + ( PZ * next_PZ ) ) / ( total_P * total_next_P )
          );

          true_beam_elastic_X.push_back( process_X );
          true_beam_elastic_Y.push_back( process_Y );
          true_beam_elastic_Z.push_back( process_Z );

        }
      }



      for( size_t i = 0; i < vertex_hits.size(); ++i ){
        std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(vertex_hits[i]) );
        for( size_t j = 0; j < ides.size(); ++j ){
          reco_beam_vertex_hits_slices.push_back( temp_hits_slices[i] );
        }
      }

      //Also, get the distance between all of the IDEs in the last slice to the location of processes in the beam trajectory 

      for( auto itProc = true_beam_proc_map.begin(); itProc != true_beam_proc_map.end(); ++itProc ){
        double procX = true_beam_trajectory.X( itProc->first );
        double procY = true_beam_trajectory.Y( itProc->first );
        double procZ = true_beam_trajectory.Z( itProc->first );

        std::cout << std::endl << "Process: " << true_beam_trajectory.KeyToProcess(itProc->second) << procX << " " << procY << " " << procZ << std::endl; 

        std::vector< double > temp_dRs;

        int nIDEs = 0;

        for( size_t i = 0; i < vertex_hits.size(); ++i ){
          

          std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(vertex_hits[i]) );
          std::cout << "Hit: " << vertex_hits[i]->WireID().Wire << " " << vertex_hits[i]->WireID().TPC  << " " << vertex_hits[i]->WireID().Plane << std::endl;
          for( size_t j = 0; j < ides.size(); ++j ){
            std::cout << "\tIDE: " << ides[j]->trackID << " " << ides[j]->x << " " << ides[j]->y << " " << ides[j]->z << std::endl;
            temp_dRs.push_back( sqrt( std::pow( (ides[j]->x - procX), 2 ) +
                                      std::pow( (ides[j]->y - procY), 2 ) +
                                      std::pow( (ides[j]->z - procZ), 2 ) ) );
            
            ++nIDEs;
          }
        }
        reco_beam_vertex_dRs.push_back( temp_dRs );

      }

      if( true_beam_endProcess.find( "Inelastic" ) == std::string::npos ){
        true_beam_processes.push_back( true_beam_endProcess );

        double procX = true_beam_endX;
        double procY = true_beam_endY;
        double procZ = true_beam_endZ;

        std::vector< double > temp_dRs;

        int nIDEs = 0;

        for( size_t i = 0; i < vertex_hits.size(); ++i ){
          std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(vertex_hits[i]) );
          for( size_t j = 0; j < ides.size(); ++j ){
          temp_dRs.push_back( sqrt( std::pow( (ides[j]->x - procX), 2 ) +
                                    std::pow( (ides[j]->y - procY), 2 ) +
                                    std::pow( (ides[j]->z - procZ), 2 ) ) );
            
            ++nIDEs;

          }
        }
        reco_beam_vertex_dRs.push_back( temp_dRs );

      }


      double IDE_max_z = 0.; 
      const sim::IDE * max_IDE = 0x0;

      for( size_t i = 0; i < vertex_hits.size(); ++i ){
        std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(vertex_hits[i]) );
        for( size_t j = 0; j < ides.size(); ++j ){
          if( ides[j]->trackID == true_beam_ID ){
            //std::cout << "Found Track: " << ides[j]->trackID << " Hit: " << i << " ide: " <<  j <<  " " << ides[j] << std::endl;
            //std::cout << "\t" << ides[j]->x << " " << ides[j]->y << " " << ides[j]->z << " " << ides[j]->energy << std::endl;
            if( ides[j]->z > IDE_max_z ){
              IDE_max_z = ides[j]->z;
              max_IDE = ides[j];
            }
          }
        }
      }
      
      std::cout << "Looking at IDEs" << std::endl;
      auto view2_IDEs = bt_serv->TrackIdToSimIDEs_Ps( true_beam_ID, fGeometryService->View(2) );
      std::sort( view2_IDEs.begin(), view2_IDEs.end(), sort_IDEs );

      double new_total_dE = 0.;
      if( max_IDE ){
        true_beam_IDE_found_in_recoVtx = true;
        for( size_t i = 0; i < view2_IDEs.size(); ++i ){
          auto theIDE = view2_IDEs[i]; 

          if( theIDE->z > IDE_max_z )
            break;

          //std::cout << view2_IDEs[i]->z << " " << view2_IDEs[i]->energy << std::endl;
          new_total_dE += view2_IDEs[i]->energy;
        }
      }
      else{
        true_beam_IDE_found_in_recoVtx = false;
        for( size_t i = 0; i < view2_IDEs.size(); ++i ){
          new_total_dE += view2_IDEs[i]->energy;
        }
      }
      true_beam_IDE_totalDep = new_total_dE;
      std::cout << "New total: " << new_total_dE << std::endl;


      //Do the true xsec measurement
      double mass = 139.57; 
      if( true_beam_PDG == 2212 ) mass = 938.27;
      else if( abs(true_beam_PDG) == 211 ) mass = 139.57;
      else if( abs(true_beam_PDG) == -11 ) mass = .511;
      else if( abs(true_beam_PDG) == 321 ) mass = 321;
      else if( abs(true_beam_PDG) == 13 )  mass = 105.66;      

      double init_KE = sqrt( 1.e6 * true_beam_startP*true_beam_startP + mass*mass ) - mass;
      true_beam_incidentEnergies.push_back( init_KE );

      double slice_end = pitch;
      double slice_edep = 0.;
      for( size_t i = 0; i < view2_IDEs.size(); ++i ){

        auto theIDE = view2_IDEs[i];

        if( theIDE->z < 0. ) continue;

        if( theIDE->z > slice_end ){
          true_beam_incidentEnergies.push_back( true_beam_incidentEnergies.back() - slice_edep ); 
          slice_edep = 0.;
          slice_end += pitch;
        }

        slice_edep += theIDE->energy; 
      }

      //Remove the last. It's not considered an 'experiment'
      true_beam_incidentEnergies.pop_back();
      if( true_beam_incidentEnergies.size() ) true_beam_interactingEnergy = true_beam_incidentEnergies.back();

    }


    //Looking at reco daughters from the reco beam track
    /*
    const std::vector<const recob::Track*> trackDaughters = pfpUtil.GetPFParticleDaughterTracks(*particle,evt,fPFParticleTag,fTrackerTag);  
    const std::vector<const recob::Shower*> showerDaughters = pfpUtil.GetPFParticleDaughterShowers(*particle,evt,fPFParticleTag,fShowerTag);  
    std::cout << "Beam particle has " << trackDaughters.size() << " track-like daughters and " << showerDaughters.size() << " shower-like daughters." << std::endl;
    std::cout << std::endl;

    std::map< int, std::vector< size_t > > sliceDaughters;

    reco_beam_nTrackDaughters = trackDaughters.size();

    for( size_t i = 0; i < trackDaughters.size(); ++i ){
      auto daughterTrack = trackDaughters.at(i);
      
      std::cout << "Track daughter " << i << " has ID " << daughterTrack->ID() << " and len " << daughterTrack->Length() << std::endl; 


      reco_daughter_len.push_back( daughterTrack->Length() );
      reco_daughter_startX.push_back( daughterTrack->Trajectory().Start().X() );
      reco_daughter_startY.push_back( daughterTrack->Trajectory().Start().Y() );
      reco_daughter_startZ.push_back( daughterTrack->Trajectory().Start().Z() );
      reco_daughter_endX.push_back( daughterTrack->Trajectory().End().X() );
      reco_daughter_endY.push_back( daughterTrack->Trajectory().End().Y() );
      reco_daughter_endZ.push_back( daughterTrack->Trajectory().End().Z() );

      reco_daughter_momByRange_proton.push_back( track_p_calc.GetTrackMomentum( daughterTrack->Length(), 2212 ) );
      reco_daughter_momByRange_muon.push_back( track_p_calc.GetTrackMomentum(   daughterTrack->Length(), 13  ) );

      //Write out Hits of Daughter Track Particles
      if( fSaveHits ){
        reco_daughter_spacePts_X.push_back( std::vector< double >() );
        reco_daughter_spacePts_Y.push_back( std::vector< double >() );
        reco_daughter_spacePts_Z.push_back( std::vector< double >() );

        for(size_t j = 0; j < daughterTrack->NumberTrajectoryPoints(); ++j){

          reco_daughter_spacePts_X.back().push_back( daughterTrack->Trajectory().LocationAtPoint(j).X() );
          reco_daughter_spacePts_Y.back().push_back( daughterTrack->Trajectory().LocationAtPoint(j).Y() );
          reco_daughter_spacePts_Z.back().push_back( daughterTrack->Trajectory().LocationAtPoint(j).Z() );

          //std::cout << "Daughter | X = " << daughterTrack->Trajectory().LocationAtPoint(j).X() << " |Y = " << daughterTrack->Trajectory().LocationAtPoint(j).Y() << " |Z = " << daughterTrack->Trajectory().LocationAtPoint(j).Z() << std::endl;

        }
      }


      double delta_start = sqrt( ( reco_daughter_startX.back() - reco_beam_endX )*( reco_daughter_startX.back() - reco_beam_endX ) + 
                                 ( reco_daughter_startY.back() - reco_beam_endY )*( reco_daughter_startY.back() - reco_beam_endY ) +
                                 ( reco_daughter_startZ.back() - reco_beam_endZ )*( reco_daughter_startZ.back() - reco_beam_endZ ) );

      double delta_end =   sqrt( ( reco_daughter_endX.back() - reco_beam_endX )*( reco_daughter_endX.back() - reco_beam_endX ) + 
                                 ( reco_daughter_endY.back() - reco_beam_endY )*( reco_daughter_endY.back() - reco_beam_endY ) +
                                 ( reco_daughter_endZ.back() - reco_beam_endZ )*( reco_daughter_endZ.back() - reco_beam_endZ ) );

      reco_daughter_deltaR.push_back( ( delta_start < delta_end ) ? delta_start : delta_end );

      reco_daughter_trackID.push_back( daughterTrack->ID() );



      std::vector< anab::Calorimetry > dummy_calo = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackerTag, fCalorimetryTag);
      auto dummy_dQdx = dummy_calo[0].dQdx();
      auto dummy_dEdx = dummy_calo[0].dEdx();
      auto dummy_Range = dummy_calo[0].ResidualRange();
 
      reco_daughter_dQdX.push_back( std::vector<double>() );   
      reco_daughter_resRange.push_back( std::vector<double>() );
      reco_daughter_dEdX.push_back( std::vector<double>() );

      for( size_t j = 0; j < dummy_dQdx.size(); ++j ){
        reco_daughter_dQdX.back().push_back( dummy_dQdx[j] );
        reco_daughter_resRange.back().push_back( dummy_Range[j] );
        reco_daughter_dEdX.back().push_back( dummy_dEdx[j] );
      }

      std::pair< double, int > this_chi2_ndof = trackUtil.Chi2PID( reco_daughter_dEdX.back(), reco_daughter_resRange.back(), templates[ 2212 ] );
      reco_daughter_Chi2_proton.push_back( this_chi2_ndof.first );
      reco_daughter_Chi2_ndof.push_back( this_chi2_ndof.second );


      //Go through the hits of each daughter and check the CNN scores 
      //Look at the hits from the track:
      std::cout << "Getting hits from daughter" << std::endl;
      std::vector< art::Ptr< recob::Hit > > daughter_hits = findHits.at( daughterTrack->ID() );

      double track_total = 0.;
      double em_total = 0.;
      double michel_total = 0.;   
      for( size_t h = 0; h < daughter_hits.size(); ++h ){
        std::array<float,4> cnn_out = hitResults.getOutput( daughter_hits[h ] );
        track_total  += cnn_out[ hitResults.getIndex("track") ];
        em_total     += cnn_out[ hitResults.getIndex("em") ];
        michel_total += cnn_out[ hitResults.getIndex("michel") ];
      }

      if( daughter_hits.size() > 0 ){
        reco_daughter_trackScore.push_back( track_total / daughter_hits.size() );
        reco_daughter_emScore.push_back( em_total / daughter_hits.size() );
        reco_daughter_michelScore.push_back( michel_total / daughter_hits.size() );
      }
      else{
        reco_daughter_trackScore.push_back( -999. );
        reco_daughter_emScore.push_back( -999. );
        reco_daughter_michelScore.push_back( -999. );
      }




      //Match the daughters to a slice
      //First, check whether the start or end of the daughter track are closer
      double d_startX = daughterTrack->Trajectory().Start().X();
      double d_startY = daughterTrack->Trajectory().Start().Y();
      double d_startZ = daughterTrack->Trajectory().Start().Z();

      double d_endX = daughterTrack->Trajectory().End().X();
      double d_endY = daughterTrack->Trajectory().End().Y();
      double d_endZ = daughterTrack->Trajectory().End().Z();

      double to_start_of_daughter = sqrt(
        ( d_startX - max_X ) * ( d_startX - max_X ) + 
        ( d_startY - max_Y ) * ( d_startY - max_Y ) + 
        ( d_startZ - max_Z ) * ( d_startZ - max_Z )  
      );
      double to_end_of_daughter = sqrt(
        ( d_endX - max_X ) * ( d_endX - max_X ) + 
        ( d_endY - max_Y ) * ( d_endY - max_Y ) + 
        ( d_endZ - max_Z ) * ( d_endZ - max_Z )  
      );

      if ( to_end_of_daughter < to_start_of_daughter ){
        reco_daughter_to_vertex.push_back( to_end_of_daughter );
      }
      else{
        reco_daughter_to_vertex.push_back( to_start_of_daughter );
      }

      double dr_start = std::numeric_limits<double>::max();
      double dr_end = std::numeric_limits<double>::max();

      size_t min_start_index = 0;
      size_t min_end_index = 0;

      for( size_t j = 0; j < thisTrack->NumberTrajectoryPoints(); ++j ){
        double X = thisTrack->Trajectory().LocationAtPoint(j).X();
        double Y = thisTrack->Trajectory().LocationAtPoint(j).Y();
        double Z = thisTrack->Trajectory().LocationAtPoint(j).Z();

        double dr = sqrt(
          ( d_startX - X ) * ( d_startX - X ) + 
          ( d_startY - Y ) * ( d_startY - Y ) + 
          ( d_startZ - Z ) * ( d_startZ - Z )  
        );

        if( dr < dr_start ){
          dr_start = dr;    
          min_start_index = j;
        }

        dr = sqrt(
          ( d_endX - X ) * ( d_endX - X ) + 
          ( d_endY - Y ) * ( d_endY - Y ) + 
          ( d_endZ - Z ) * ( d_endZ - Z )  
        );

        if( dr < dr_end ){
          dr_end = dr;    
          min_end_index = j;
        }

      }
      
      size_t min_index = 0;
      if( dr_end < dr_start ){
        min_index = min_end_index; 

        std::cout << "dr between track and daughter: " << dr_end << std::endl;
        reco_daughter_dR.push_back( dr_end );
 
      }
      else{
        min_index = min_start_index; 
        std::cout << "dr between track and daughter: " << dr_start << std::endl;
        reco_daughter_dR.push_back( dr_start );
      }

      double min_z = thisTrack->Trajectory().LocationAtPoint(min_index).Z();
      int slice = std::floor( ( min_z - z0 ) / pitch );
      sliceDaughters[ slice ].push_back( i );
      reco_daughter_slice.push_back( slice );


      if( !evt.isRealData() ){
        reco_daughter_true_byE_completeness.push_back( truthUtil.GetCompleteness( *daughterTrack, evt, fTrackerTag, fHitTag ) );
        reco_daughter_true_byE_purity.push_back( truthUtil.GetPurity( *daughterTrack, evt, fTrackerTag#, fHitTag#) );

        //const simb::MCParticle * match = truthUtil.GetMCParticleByHits( *daughterTrack, evt, fTrackerTag, fHitTag );
        protoana::MCParticleSharedHits match = truthUtil.GetMCParticleByHits( *daughterTrack, evt, fTrackerTag, fHitTag );

        if( match.particle ){
          std::cout << std::endl << "Match: " << match.particle->PdgCode() << " " << match.particle->TrackId() << std::endl;
           
          reco_daughter_true_byHits_PDG.push_back( match.particle->PdgCode() );
          std::cout << "PDG" << std::endl;
          reco_daughter_true_byHits_ID.push_back( match.particle->TrackId() );
          std::cout << "ID" << std::endl;
          reco_daughter_true_byHits_parID.push_back( match.particle->Mother() );
          std::cout << "parID " << match.particle->Mother() << std::endl;
          //reco_daughter_true_byHits_parPDG.push_back( pi_serv->TrackIdToMotherParticle_P( match.particle->TrackId() )->PdgCode() );
          reco_daughter_true_byHits_parPDG.push_back(  
            ( (match.particle->Mother() > 0) ? plist[ match.particle->Mother() ]->PdgCode() : 0 )
          );
          std::cout << "parPDG" << std::endl;
          reco_daughter_true_byHits_process.push_back( match.particle->Process() );
          std::cout << "proc" << std::endl;
          reco_daughter_true_byHits_origin.push_back( 
            pi_serv->TrackIdToMCTruth_P(match.particle->TrackId())->Origin()
          );
          std::cout << "origin" << std::endl;
          reco_daughter_true_byHits_sharedHits.push_back( match.nSharedHits ); 
          std::cout << "shared" << std::endl;
          reco_daughter_true_byHits_emHits.push_back( match.nSharedDeltaRayHits ); 
          std::cout << "em" << std::endl;

          reco_daughter_true_byHits_len.push_back( match.particle->Trajectory().TotalLength() );
          reco_daughter_true_byHits_startX.push_back( match.particle->Position(0).X() );
          reco_daughter_true_byHits_startY.push_back( match.particle->Position(0).Y() );
          reco_daughter_true_byHits_startZ.push_back( match.particle->Position(0).Z() );

          reco_daughter_true_byHits_endX.push_back( match.particle->EndPosition().X() );
          reco_daughter_true_byHits_endY.push_back( match.particle->EndPosition().Y() );
          reco_daughter_true_byHits_endZ.push_back( match.particle->EndPosition().Z() );

          reco_daughter_true_byHits_startPx.push_back( match.particle->Px() );
          reco_daughter_true_byHits_startPy.push_back( match.particle->Py() );
          reco_daughter_true_byHits_startPz.push_back( match.particle->Pz() );
          reco_daughter_true_byHits_startP.push_back(
                          sqrt(match.particle->Px()*match.particle->Px() + 
                                  match.particle->Py()*match.particle->Py() + 
                                  match.particle->Pz()*match.particle->Pz()) );
          reco_daughter_true_byHits_startE.push_back( match.particle->E() );

          auto list = truthUtil.GetMCParticleListByHits( *daughterTrack, evt, fTrackerTag, fHitTag );
          for( size_t j = 0; j < list.size(); ++j ){
          //  std::cout << "Contrib " << j << " " << list[j].first->TrackId() << " " << list[j].second << std::endl;
            std::cout << "Contrib " << j << " " << list[j].particle->TrackId() << " " << list[j].particle->PdgCode() << " " << list[j].nSharedHits << " " << list[j].nSharedDeltaRayHits << std::endl;
          }
          double total = 0.;
          double matched_hits = 0.;
          for( size_t j = 0; j < list.size(); ++j ){
          //  std::cout << "Contrib " << j << " " << list[j].first->TrackId() << " " << list[j].second << std::endl;
            std::cout << "Contrib " << j << " " << list[j].particle->TrackId() << " " << list[j].particle->PdgCode() << " " << list[j].nSharedHits << " " << list[j].nSharedDeltaRayHits << std::endl;

            if( list[j].particle == match.particle ){
               matched_hits = list[j].nSharedHits + list[j].nSharedDeltaRayHits;
            }

            total += list[j].nSharedHits + list[j].nSharedDeltaRayHits;
          }

          reco_daughter_true_byHits_purity.push_back( matched_hits / total );
        }
        else{
          reco_daughter_true_byHits_PDG.push_back( -1 );
          reco_daughter_true_byHits_ID.push_back( -1 );
          reco_daughter_true_byHits_origin.push_back( -1 );
          reco_daughter_true_byHits_parID.push_back( -1 );
          reco_daughter_true_byHits_parPDG.push_back( -1 );
          reco_daughter_true_byHits_process.push_back( "empty" );
          reco_daughter_true_byHits_sharedHits.push_back( 0 ); 
          reco_daughter_true_byHits_emHits.push_back( 0 ); 

          reco_daughter_true_byHits_len.push_back( 0. );
          reco_daughter_true_byHits_startX.push_back( 0. );
          reco_daughter_true_byHits_startY.push_back( 0. );
          reco_daughter_true_byHits_startZ.push_back( 0. );
          reco_daughter_true_byHits_endX.push_back( 0. );
          reco_daughter_true_byHits_endY.push_back( 0. );
          reco_daughter_true_byHits_endZ.push_back( 0. );
          reco_daughter_true_byHits_purity.push_back( 0. );
          reco_daughter_true_byHits_startPx.push_back( 0. );
          reco_daughter_true_byHits_startPy.push_back( 0. );
          reco_daughter_true_byHits_startPz.push_back( 0. );
          reco_daughter_true_byHits_startP.push_back( 0. );
          reco_daughter_true_byHits_startE.push_back( 0. );
        }



        //For this reco daughter track, get the actual particle contributing to it
        const simb::MCParticle* trueDaughterParticle = truthUtil.GetMCParticleFromRecoTrack(*daughterTrack, evt, fTrackerTag);
        if( trueDaughterParticle ){
          reco_daughter_true_byE_PDG.push_back( trueDaughterParticle->PdgCode() );
          reco_daughter_true_byE_ID.push_back( trueDaughterParticle->TrackId() );
          if( trueDaughterParticle->TrackId() == true_beam_particle->TrackId() ){
            reco_daughter_true_byE_isPrimary = true;
          }
            
          reco_daughter_true_byE_origin.push_back( 
            pi_serv->TrackIdToMCTruth_P(trueDaughterParticle->TrackId())->Origin()
          );
          reco_daughter_true_byE_parID.push_back( trueDaughterParticle->Mother() );
          //reco_daughter_true_byE_parPDG.push_back( pi_serv->TrackIdToMotherParticle_P( trueDaughterParticle->TrackId() )->PdgCode() );
          reco_daughter_true_byE_parPDG.push_back(
            ( (trueDaughterParticle->Mother() > 0) ? plist[ trueDaughterParticle->Mother() ]->PdgCode() : 0 )
          );
          reco_daughter_true_byE_process.push_back( trueDaughterParticle->Process() );
        }
        else{
          reco_daughter_true_byE_PDG.push_back( -1 );
          reco_daughter_true_byE_ID.push_back( -1 );
          reco_daughter_true_byE_origin.push_back( -1 );
          reco_daughter_true_byE_parID.push_back( -1 );
          reco_daughter_true_byE_parPDG.push_back( -1 );
          reco_daughter_true_byE_process.push_back( "empty" );
        }
      }
    }

    ////Showers
    reco_beam_nShowerDaughters = showerDaughters.size();

    for( size_t i = 0; i < showerDaughters.size(); ++i ){
      auto daughterShowerFromRecoTrack = showerDaughters[i];

      std::cout << "Shower daughter " << i << " Starts at " << daughterShowerFromRecoTrack->ShowerStart().X() << " " << daughterShowerFromRecoTrack->ShowerStart().Y() << " " << daughterShowerFromRecoTrack->ShowerStart().Z() << std::endl;
      reco_daughter_showerID.push_back( daughterShowerFromRecoTrack->ID() );
      
      double d_startX = daughterShowerFromRecoTrack->ShowerStart().X();
      double d_startY = daughterShowerFromRecoTrack->ShowerStart().Y();
      double d_startZ = daughterShowerFromRecoTrack->ShowerStart().Z();

      reco_daughter_shower_startX.push_back( d_startX );
      reco_daughter_shower_startY.push_back( d_startY );
      reco_daughter_shower_startZ.push_back( d_startZ );

      
      double to_start_of_daughter = sqrt(
        ( d_startX - max_X ) * ( d_startX - max_X ) + 
        ( d_startY - max_Y ) * ( d_startY - max_Y ) + 
        ( d_startZ - max_Z ) * ( d_startZ - max_Z )  
      );

      reco_daughter_shower_to_vertex.push_back( to_start_of_daughter );


      reco_daughter_shower_len.push_back( daughterShowerFromRecoTrack->Length() );

      //Write out Hits of Daughter Shower Particles, need to associate hits with space points
      if( fSaveHits ){
        reco_daughter_shower_spacePts_X.push_back( std::vector< double >() );
        reco_daughter_shower_spacePts_Y.push_back( std::vector< double >() );
        reco_daughter_shower_spacePts_Z.push_back( std::vector< double >() );

        const std::vector< const recob::Hit*> sh_hits = showerUtil.GetRecoShowerHits(*daughterShowerFromRecoTrack, evt, fShowerTag);
        art::FindManyP< recob::SpacePoint > spFromShowerHits(sh_hits, evt, "pandora");
     
        for(size_t j = 0; j < sh_hits.size(); j++){

          std::vector< art::Ptr< recob::SpacePoint >> sp = spFromShowerHits.at(j);
          if( !sp.empty() ){
            //sp[0] apparently there can be more spacePoints for one shower Hit...?

            reco_daughter_shower_spacePts_X.back().push_back(sp[0]->XYZ()[0]); 
            reco_daughter_shower_spacePts_Y.back().push_back(sp[0]->XYZ()[1]);
            reco_daughter_shower_spacePts_Z.back().push_back(sp[0]->XYZ()[2]);
            
            //std::cout << "SHOWER Daughter | X = " << sp[0]->XYZ()[0] << " |Y = " << sp[0]->XYZ()[1] << " |Z = " << sp[0]->XYZ()[2] << std::endl;

          }
        }
      }


      std::vector< anab::Calorimetry > dummy_calo = showerUtil.GetRecoShowerCalorimetry(*daughterShowerFromRecoTrack, evt, fShowerTag, "pandoraShowercalo");

      reco_daughter_shower_dQdX.push_back( std::vector<double>() );   
      reco_daughter_shower_resRange.push_back( std::vector<double>() );
      reco_daughter_shower_dEdX.push_back( std::vector<double>() );

      if( dummy_calo.size() ){
        std::cout << "Plane: " << dummy_calo[0].PlaneID().toString() << std::endl;      
        auto dummy_dQdx = dummy_calo[0].dQdx();
        auto dummy_dEdx = dummy_calo[0].dEdx();
        auto dummy_Range = dummy_calo[0].ResidualRange();
 
        for( size_t j = 0; j < dummy_dQdx.size(); ++j ){
          reco_daughter_shower_dQdX.back().push_back( dummy_dQdx[j] );
          reco_daughter_shower_resRange.back().push_back( dummy_Range[j] );
          reco_daughter_shower_dEdX.back().push_back( dummy_dEdx[j] );
        }
      }

      std::pair< double, int > this_chi2_ndof = trackUtil.Chi2PID( reco_daughter_shower_dEdX.back(), reco_daughter_shower_resRange.back(), templates[ 2212 ] );
      reco_daughter_shower_Chi2_proton.push_back( this_chi2_ndof.first );
      reco_daughter_shower_Chi2_ndof.push_back( this_chi2_ndof.second );

      //Go through the hits of each daughter and check the CNN scores 
      //Look at the hits from the track:
      int shower_index = showerUtil.GetShowerIndex( *daughterShowerFromRecoTrack, evt, fShowerTag );
      std::vector< art::Ptr< recob::Hit > > daughter_hits = findHitsFromShowers.at( shower_index );
      std::cout << "Got " << daughter_hits.size() << " hits from shower " << daughterShowerFromRecoTrack->ID() << " " << shower_index << std::endl;

      double track_total = 0.;
      double em_total = 0.;
      double michel_total = 0.;   
      for( size_t h = 0; h < daughter_hits.size(); ++h ){
        std::array<float,4> cnn_out = hitResults.getOutput( daughter_hits[h ] );
        track_total  += cnn_out[ hitResults.getIndex("track") ];
        em_total     += cnn_out[ hitResults.getIndex("em") ];
        michel_total += cnn_out[ hitResults.getIndex("michel") ];
      }

      if( daughter_hits.size() > 0 ){
        reco_daughter_shower_trackScore.push_back( track_total / daughter_hits.size() );
        reco_daughter_shower_emScore.push_back( em_total / daughter_hits.size() );
        reco_daughter_shower_michelScore.push_back( michel_total / daughter_hits.size() );
      }
      else{
        reco_daughter_shower_trackScore.push_back( -999. );
        reco_daughter_shower_emScore.push_back( -999. );
        reco_daughter_shower_michelScore.push_back( -999. );
      }
      
      if( !evt.isRealData() ){
        //For this reco daughter track, get the actual particle contributing to it
        const simb::MCParticle* trueDaughterParticle = truthUtil.GetMCParticleFromRecoShower(*daughterShowerFromRecoTrack, evt, fShowerTag);
        if( trueDaughterParticle ){
          reco_daughter_shower_true_byE_PDG.push_back( trueDaughterParticle->PdgCode() );
          reco_daughter_shower_true_byE_ID.push_back( trueDaughterParticle->TrackId() );
          reco_daughter_shower_true_byE_origin.push_back( 
            pi_serv->TrackIdToMCTruth_P(trueDaughterParticle->TrackId())->Origin()
          );
          reco_daughter_shower_true_byE_parID.push_back( trueDaughterParticle->Mother() );
          //reco_daughter_shower_true_byE_parPDG.push_back( pi_serv->TrackIdToMotherParticle_P( trueDaughterParticle->TrackId() )->PdgCode() );
          reco_daughter_shower_true_byE_parPDG.push_back( 
            ( trueDaughterParticle->Mother() > 0 ? plist[ trueDaughterParticle->Mother() ]->PdgCode() : 0 )
          );
          reco_daughter_shower_true_byE_startPx.push_back( trueDaughterParticle->Px());
          reco_daughter_shower_true_byE_startPy.push_back( trueDaughterParticle->Py());
          reco_daughter_shower_true_byE_startPz.push_back( trueDaughterParticle->Pz());
          reco_daughter_shower_true_byE_startP.push_back( trueDaughterParticle->P());
          reco_daughter_shower_true_byE_endProcess.push_back( trueDaughterParticle->EndProcess());
        }
        else{
          reco_daughter_shower_true_byE_PDG.push_back( -1 );
          reco_daughter_shower_true_byE_ID.push_back( -1 );
          reco_daughter_shower_true_byE_origin.push_back( -1 );
          reco_daughter_shower_true_byE_parID.push_back( -1 );
          reco_daughter_shower_true_byE_parPDG.push_back( -1 );
          reco_daughter_shower_true_byE_startPx.push_back(-999.);
          reco_daughter_shower_true_byE_startPy.push_back(-999.);
          reco_daughter_shower_true_byE_startPz.push_back(-999.);
          reco_daughter_shower_true_byE_startP.push_back(-999.);
          reco_daughter_shower_true_byE_endProcess.push_back("empty");
        }

        const simb::MCParticle * match = truthUtil.GetMCParticleByHits( *daughterShowerFromRecoTrack, evt, fShowerTag, fHitTag ).particle;

        if( match ){
          std::cout << std::endl << "Match: " << match->PdgCode() << " " << match->TrackId() << std::endl;
           
          reco_daughter_shower_true_byHits_PDG.push_back( match->PdgCode() );
          reco_daughter_shower_true_byHits_ID.push_back( match->TrackId() );
          reco_daughter_shower_true_byHits_parID.push_back( match->Mother() );
          //reco_daughter_shower_true_byHits_parPDG.push_back( pi_serv->TrackIdToMotherParticle_P( match->TrackId() )->PdgCode() );
          reco_daughter_shower_true_byHits_parPDG.push_back( 
            ( ( match->Mother() > 0 ) ? plist[ match->Mother() ]->PdgCode() : 0 )
          );
          reco_daughter_shower_true_byHits_process.push_back( match->Process() );
          reco_daughter_shower_true_byHits_origin.push_back( 
            pi_serv->TrackIdToMCTruth_P(match->TrackId())->Origin()
          );
          reco_daughter_shower_true_byHits_startPx.push_back( match->Px());
          reco_daughter_shower_true_byHits_startPy.push_back( match->Py());
          reco_daughter_shower_true_byHits_startPz.push_back( match->Pz());
          reco_daughter_shower_true_byHits_startP.push_back( match->P());
          reco_daughter_shower_true_byHits_endProcess.push_back( match->EndProcess());
        }
        else{
          reco_daughter_shower_true_byHits_PDG.push_back( -1 );
          reco_daughter_shower_true_byHits_ID.push_back( -1 );
          reco_daughter_shower_true_byHits_origin.push_back( -1 );
          reco_daughter_shower_true_byHits_parID.push_back( -1 );
          reco_daughter_shower_true_byHits_parPDG.push_back( -1 );
          reco_daughter_shower_true_byHits_process.push_back( "empty" );
          reco_daughter_shower_true_byHits_startPx.push_back(-999.);
          reco_daughter_shower_true_byHits_startPy.push_back(-999.);
          reco_daughter_shower_true_byHits_startPz.push_back(-999.);
          reco_daughter_shower_true_byHits_startP.push_back(-999.);
          reco_daughter_shower_true_byHits_endProcess.push_back("empty");
        }

        auto list = truthUtil.GetMCParticleListByHits( *daughterShowerFromRecoTrack, evt, fShowerTag, fHitTag );
        double total = 0.;
        double matched_hits = 0.;
        for( size_t j = 0; j < list.size(); ++j ){
        //  std::cout << "Contrib " << j << " " << list[j].first->TrackId() << " " << list[j].second << std::endl;
          std::cout << "Contrib " << j << " " << list[j].particle->TrackId() << " " << list[j].particle->PdgCode() << " " << list[j].nSharedHits << " " << list[j].nSharedDeltaRayHits << std::endl;

          if( list[j].particle == match ){
             matched_hits = list[j].nSharedHits + list[j].nSharedDeltaRayHits;
          }

          total += list[j].nSharedHits + list[j].nSharedDeltaRayHits;
        }

        reco_daughter_shower_true_byHits_purity.push_back( matched_hits / total );
      }
    }
    */
    ////End Of Daughters


    // Alternative Reconstruction.
    //
    // Loop over all of the PFParticles associated as daughters.
    // Then, check the CNN score (later implement the GNN score)
    //
    // Also, get the forced-tracking (pandora2) and 
    // get calorimetry + other info
    
    for( size_t daughterID : particle->Daughters() ){
      const recob::PFParticle * daughterPFP = &(pfpVec->at( daughterID ));
      reco_daughter_PFP_ID.push_back( daughterID );

      const std::vector< art::Ptr< recob::Hit > > daughterPFP_hits = pfpUtil.GetPFParticleHits_Ptrs( *daughterPFP, evt, fPFParticleTag );
      std::cout << "Got " << daughterPFP_hits.size() << " hits from daughter " << daughterID << std::endl;

      reco_daughter_PFP_nHits.push_back( daughterPFP_hits.size() );

      double track_total = 0.;
      double em_total = 0.;
      double michel_total = 0.;   
      double none_total = 0.;
      for( size_t h = 0; h < daughterPFP_hits.size(); ++h ){
        std::array<float,4> cnn_out = hitResults.getOutput( daughterPFP_hits[h] );
        track_total  += cnn_out[ hitResults.getIndex("track") ];
        em_total     += cnn_out[ hitResults.getIndex("em") ];
        michel_total += cnn_out[ hitResults.getIndex("michel") ];
        none_total   += cnn_out[ hitResults.getIndex("none") ];

      }

      cnnOutput2D theCNNResults = GetCNNOutputFromPFParticle( *daughterPFP, evt, hitResults, pfpUtil, fPFParticleTag );

      std::cout << "Testing new CNN: " << std::endl;
      std::cout << track_total << " " << theCNNResults.track << std::endl;
      std::cout << em_total << " " << theCNNResults.em << std::endl;
      std::cout << michel_total << " " << theCNNResults.michel << std::endl;
      std::cout << none_total << " " << theCNNResults.none << std::endl;

      const std::vector< const recob::SpacePoint* > spVec = pfpUtil.GetPFParticleSpacePoints( *daughterPFP, evt, fPFParticleTag );
      std::cout << "Got " << spVec.size() << " SpacePoints" << std::endl;

      if( daughterPFP_hits.size() > 0 ){
        reco_daughter_PFP_trackScore.push_back( track_total / daughterPFP_hits.size() );
        reco_daughter_PFP_emScore.push_back( em_total / daughterPFP_hits.size() );
        reco_daughter_PFP_michelScore.push_back( michel_total / daughterPFP_hits.size() );
      }
      else{
        reco_daughter_PFP_trackScore.push_back( -999. );
        reco_daughter_PFP_emScore.push_back( -999. );
        reco_daughter_PFP_michelScore.push_back( -999. );
      }

      cnnOutput2D cnn_collection = GetCNNOutputFromPFParticleFromPlane( *daughterPFP, evt, hitResults, pfpUtil, fPFParticleTag, 2 );
      if( cnn_collection.nHits > 0 ){
        reco_daughter_PFP_trackScore_collection.push_back( cnn_collection.track / cnn_collection.nHits );
        reco_daughter_PFP_emScore_collection.push_back( cnn_collection.em / cnn_collection.nHits );
        reco_daughter_PFP_michelScore_collection.push_back( cnn_collection.michel / cnn_collection.nHits );
      }
      else{
        reco_daughter_PFP_trackScore_collection.push_back( -999. );
        reco_daughter_PFP_emScore_collection.push_back( -999. );
        reco_daughter_PFP_michelScore_collection.push_back( -999. );
      }
      


      if( !evt.isRealData() ){
        protoana::MCParticleSharedHits match = truthUtil.GetMCParticleByHits( *daughterPFP, evt, fPFParticleTag, fHitTag );

        if( match.particle ){
          std::cout << std::endl << "Match: " << match.particle->PdgCode() << " " << match.particle->TrackId() << std::endl;
           
          reco_daughter_PFP_true_byHits_PDG.push_back( match.particle->PdgCode() );
          reco_daughter_PFP_true_byHits_ID.push_back( match.particle->TrackId() );
          reco_daughter_PFP_true_byHits_parID.push_back( match.particle->Mother() );
          //reco_daughter_PFP_true_byHits_parPDG.push_back( pi_serv->TrackIdToMotherParticle_P( match.particle->TrackId() )->PdgCode() );
          reco_daughter_PFP_true_byHits_parPDG.push_back( 
            ( (match.particle->Mother() > 0) ? plist[ match.particle->Mother() ]->PdgCode() : 0 )
          );

          std::cout << "Checking par: " << match.particle->TrackId() << " " << pi_serv->TrackIdToMotherParticle_P( match.particle->TrackId() )->TrackId() << std::endl;

          reco_daughter_PFP_true_byHits_process.push_back( match.particle->Process() );
          reco_daughter_PFP_true_byHits_origin.push_back( 
            pi_serv->TrackIdToMCTruth_P(match.particle->TrackId())->Origin()
          );
          reco_daughter_PFP_true_byHits_sharedHits.push_back( match.nSharedHits ); 
          reco_daughter_PFP_true_byHits_emHits.push_back( match.nSharedDeltaRayHits ); 

          reco_daughter_PFP_true_byHits_len.push_back( match.particle->Trajectory().TotalLength() );
          reco_daughter_PFP_true_byHits_startX.push_back( match.particle->Position(0).X() );
          reco_daughter_PFP_true_byHits_startY.push_back( match.particle->Position(0).Y() );
          reco_daughter_PFP_true_byHits_startZ.push_back( match.particle->Position(0).Z() );

          reco_daughter_PFP_true_byHits_endX.push_back( match.particle->EndPosition().X() );
          reco_daughter_PFP_true_byHits_endY.push_back( match.particle->EndPosition().Y() );
          reco_daughter_PFP_true_byHits_endZ.push_back( match.particle->EndPosition().Z() );

          reco_daughter_PFP_true_byHits_startPx.push_back( match.particle->Px() );
          reco_daughter_PFP_true_byHits_startPy.push_back( match.particle->Py() );
          reco_daughter_PFP_true_byHits_startPz.push_back( match.particle->Pz() );
          reco_daughter_PFP_true_byHits_startE.push_back( match.particle->E() );
          reco_daughter_PFP_true_byHits_startP.push_back( 
                          sqrt(match.particle->Px()*match.particle->Px() + 
                                  match.particle->Py()*match.particle->Py() + 
                                  match.particle->Pz()*match.particle->Pz()) );
          reco_daughter_PFP_true_byHits_endProcess.push_back( match.particle->EndProcess());

          auto list = truthUtil.GetMCParticleListByHits( *daughterPFP, evt, fPFParticleTag, fHitTag );
          double total = 0.;
          double matched_hits = 0.;
          for( size_t j = 0; j < list.size(); ++j ){
          //  std::cout << "Contrib " << j << " " << list[j].first->TrackId() << " " << list[j].second << std::endl;
            //std::cout << "Contrib " << j << " " << list[j].particle->TrackId() << " " << list[j].particle->PdgCode() << " " << list[j].nSharedHits << " " << list[j].nSharedDeltaRayHits << std::endl;

            if( list[j].particle == match.particle ){
               matched_hits = list[j].nSharedHits + list[j].nSharedDeltaRayHits;
            }

            total += list[j].nSharedHits + list[j].nSharedDeltaRayHits;
          }

          reco_daughter_PFP_true_byHits_purity.push_back( matched_hits / total );

        }
        else{
          reco_daughter_PFP_true_byHits_PDG.push_back( -1 );
          reco_daughter_PFP_true_byHits_ID.push_back( -1 );
          reco_daughter_PFP_true_byHits_origin.push_back( -1 );
          reco_daughter_PFP_true_byHits_parID.push_back( -1 );
          reco_daughter_PFP_true_byHits_parPDG.push_back( -1 );
          reco_daughter_PFP_true_byHits_process.push_back( "empty" );
          reco_daughter_PFP_true_byHits_sharedHits.push_back( 0 ); 
          reco_daughter_PFP_true_byHits_emHits.push_back( 0 ); 

          reco_daughter_PFP_true_byHits_len.push_back( -999. );
          reco_daughter_PFP_true_byHits_startX.push_back( -999. );
          reco_daughter_PFP_true_byHits_startY.push_back( -999. );
          reco_daughter_PFP_true_byHits_startZ.push_back( -999. );
          reco_daughter_PFP_true_byHits_endX.push_back( -999. );
          reco_daughter_PFP_true_byHits_endY.push_back( -999. );
          reco_daughter_PFP_true_byHits_endZ.push_back( -999. );
          reco_daughter_PFP_true_byHits_startPx.push_back( -999. );
          reco_daughter_PFP_true_byHits_startPy.push_back( -999. );
          reco_daughter_PFP_true_byHits_startPz.push_back( -999. );
          reco_daughter_PFP_true_byHits_startP.push_back( -999. );
          reco_daughter_PFP_true_byHits_startE.push_back( -999. );
          reco_daughter_PFP_true_byHits_endProcess.push_back("empty");
          reco_daughter_PFP_true_byHits_purity.push_back( -999. );
        }
      }

      try{
        const recob::Track* pandora2Track = pfpUtil.GetPFParticleTrack( *daughterPFP, evt, fPFParticleTag, "pandora2Track" );
        std::cout << "pandora2 track: " << pandora2Track << std::endl;

        if( pandora2Track ){
          reco_daughter_allTrack_ID.push_back( pandora2Track->ID() );

          std::vector< anab::Calorimetry > dummy_calo = trackUtil.GetRecoTrackCalorimetry(*pandora2Track, evt, "pandora2Track", "pandora2calo");
          std::vector< anab::Calorimetry > dummy_caloSCE = trackUtil.GetRecoTrackCalorimetry(*pandora2Track, evt, "pandora2Track", "pandora2caloSCE");

          auto dummy_dEdx = dummy_calo[0].dEdx();
          auto dummy_dQdx = dummy_calo[0].dQdx();
          auto dummy_Range = dummy_calo[0].ResidualRange();

          auto dummy_dEdx_SCE = dummy_caloSCE[0].dEdx();
          auto dummy_dQdx_SCE = dummy_caloSCE[0].dQdx();
          auto dummy_Range_SCE = dummy_caloSCE[0].ResidualRange();

          std::vector< float > cali_dEdX = calibration.GetCalibratedCalorimetry(  *pandora2Track, evt, "pandora2Track", "pandora2calo" );
          std::vector< float > cali_dEdX_SCE = calibration.GetCalibratedCalorimetry(  *pandora2Track, evt, "pandora2Track", "pandora2caloSCE" );
 
          reco_daughter_allTrack_resRange.push_back( std::vector<double>() );
          reco_daughter_allTrack_dEdX.push_back( std::vector<double>() );
          reco_daughter_allTrack_dQdX.push_back( std::vector<double>() );

          for( size_t j = 0; j < dummy_dEdx.size(); ++j ){
            reco_daughter_allTrack_resRange.back().push_back( dummy_Range[j] );
            reco_daughter_allTrack_dEdX.back().push_back( dummy_dEdx[j] );
            reco_daughter_allTrack_dQdX.back().push_back( dummy_dQdx[j] );
          }

          reco_daughter_allTrack_resRange_SCE.push_back( std::vector<double>() );
          reco_daughter_allTrack_dEdX_SCE.push_back( std::vector<double>() );
          reco_daughter_allTrack_dQdX_SCE.push_back( std::vector<double>() );

          for( size_t j = 0; j < dummy_dEdx_SCE.size(); ++j ){
            reco_daughter_allTrack_resRange_SCE.back().push_back( dummy_Range_SCE[j] );
            reco_daughter_allTrack_dEdX_SCE.back().push_back( dummy_dEdx_SCE[j] );
            reco_daughter_allTrack_dQdX_SCE.back().push_back( dummy_dQdx_SCE[j] );
          }

          reco_daughter_allTrack_calibrated_dEdX.push_back( std::vector<double>() );
          reco_daughter_allTrack_calibrated_dEdX_SCE.push_back( std::vector<double>() );
          for( size_t j = 0; j < cali_dEdX.size(); ++j ){
            reco_daughter_allTrack_calibrated_dEdX.back().push_back( cali_dEdX[j] );
          }
          for( size_t j = 0; j < cali_dEdX_SCE.size(); ++j ){
            reco_daughter_allTrack_calibrated_dEdX_SCE.back().push_back( cali_dEdX_SCE[j] );
          } 

          std::pair< double, int > this_chi2_ndof = trackUtil.Chi2PID( reco_daughter_allTrack_calibrated_dEdX_SCE.back(), reco_daughter_allTrack_resRange.back(), templates[ 2212 ] );
          reco_daughter_allTrack_Chi2_proton.push_back( this_chi2_ndof.first );
          reco_daughter_allTrack_Chi2_ndof.push_back( this_chi2_ndof.second );
          
          reco_daughter_allTrack_Theta.push_back(  pandora2Track->Theta() );
          reco_daughter_allTrack_Phi.push_back(  pandora2Track->Phi() );

          reco_daughter_allTrack_len.push_back(    pandora2Track->Length() );
          reco_daughter_allTrack_startX.push_back( pandora2Track->Trajectory().Start().X() );
          reco_daughter_allTrack_startY.push_back( pandora2Track->Trajectory().Start().Y() );
          reco_daughter_allTrack_startZ.push_back( pandora2Track->Trajectory().Start().Z() );
          reco_daughter_allTrack_endX.push_back(   pandora2Track->Trajectory().End().X() );
          reco_daughter_allTrack_endY.push_back(   pandora2Track->Trajectory().End().Y() );
          reco_daughter_allTrack_endZ.push_back(   pandora2Track->Trajectory().End().Z() );

          reco_daughter_allTrack_momByRange_proton.push_back( track_p_calc.GetTrackMomentum( pandora2Track->Length(), 2212 ) );
          reco_daughter_allTrack_momByRange_muon.push_back(   track_p_calc.GetTrackMomentum( pandora2Track->Length(), 13  ) );

          //Match the daughters to a slice
          //First, check whether the start or end of the daughter track are closer
          double d_startX = pandora2Track->Trajectory().Start().X();
          double d_startY = pandora2Track->Trajectory().Start().Y();
          double d_startZ = pandora2Track->Trajectory().Start().Z();

          double d_endX = pandora2Track->Trajectory().End().X();
          double d_endY = pandora2Track->Trajectory().End().Y();
          double d_endZ = pandora2Track->Trajectory().End().Z();

          double to_start_of_daughter = sqrt(
            ( d_startX - max_X ) * ( d_startX - max_X ) + 
            ( d_startY - max_Y ) * ( d_startY - max_Y ) + 
            ( d_startZ - max_Z ) * ( d_startZ - max_Z )  
          );
          double to_end_of_daughter = sqrt(
            ( d_endX - max_X ) * ( d_endX - max_X ) + 
            ( d_endY - max_Y ) * ( d_endY - max_Y ) + 
            ( d_endZ - max_Z ) * ( d_endZ - max_Z )  
          );

          if ( to_end_of_daughter < to_start_of_daughter ){
            reco_daughter_allTrack_to_vertex.push_back( to_end_of_daughter );
          }
          else{
            reco_daughter_allTrack_to_vertex.push_back( to_start_of_daughter );
          }

          double dr_start = std::numeric_limits<double>::max();
          double dr_end = std::numeric_limits<double>::max();

          //size_t min_start_index = 0;
          //size_t min_end_index = 0;

          for( size_t j = 0; j < thisTrack->NumberTrajectoryPoints(); ++j ){
            double X = thisTrack->Trajectory().LocationAtPoint(j).X();
            double Y = thisTrack->Trajectory().LocationAtPoint(j).Y();
            double Z = thisTrack->Trajectory().LocationAtPoint(j).Z();

            double dr = sqrt(
              ( d_startX - X ) * ( d_startX - X ) + 
              ( d_startY - Y ) * ( d_startY - Y ) + 
              ( d_startZ - Z ) * ( d_startZ - Z )  
            );

            if( dr < dr_start ){
              dr_start = dr;    
              //min_start_index = j;
            }

            dr = sqrt(
              ( d_endX - X ) * ( d_endX - X ) + 
              ( d_endY - Y ) * ( d_endY - Y ) + 
              ( d_endZ - Z ) * ( d_endZ - Z )  
            );

            if( dr < dr_end ){
              dr_end = dr;    
              //min_end_index = j;
            }

          }
          
          //size_t min_index = 0;
          if( dr_end < dr_start ){
            //min_index = min_end_index; 

            std::cout << "dr between track and daughter: " << dr_end << std::endl;
            reco_daughter_allTrack_dR.push_back( dr_end );
 
          }
          else{
            //min_index = min_start_index; 
            std::cout << "dr between track and daughter: " << dr_start << std::endl;
            reco_daughter_allTrack_dR.push_back( dr_start );
          }

        }
        else{
          reco_daughter_allTrack_ID.push_back( -1 );
          reco_daughter_allTrack_resRange.push_back( std::vector<double>() );
          reco_daughter_allTrack_dEdX.push_back( std::vector<double>() );
          reco_daughter_allTrack_dQdX.push_back( std::vector<double>() );
          reco_daughter_allTrack_resRange_SCE.push_back( std::vector<double>() );
          reco_daughter_allTrack_dEdX_SCE.push_back( std::vector<double>() );
          reco_daughter_allTrack_dQdX_SCE.push_back( std::vector<double>() );

          reco_daughter_allTrack_calibrated_dEdX.push_back( std::vector<double>() );
          reco_daughter_allTrack_calibrated_dEdX_SCE.push_back( std::vector<double>() );


          reco_daughter_allTrack_Chi2_proton.push_back( -999. );
          reco_daughter_allTrack_Chi2_ndof.push_back( 0 );
          reco_daughter_allTrack_Theta.push_back(-999. );
          reco_daughter_allTrack_Phi.push_back(-999.);
          reco_daughter_allTrack_len.push_back( -999. );
          reco_daughter_allTrack_startX.push_back( -999. );
          reco_daughter_allTrack_startY.push_back( -999. );
          reco_daughter_allTrack_startZ.push_back( -999. );
          reco_daughter_allTrack_endX.push_back(   -999. );
          reco_daughter_allTrack_endY.push_back(   -999. );
          reco_daughter_allTrack_endZ.push_back(   -999. );
          reco_daughter_allTrack_to_vertex.push_back( -999. );
          reco_daughter_allTrack_dR.push_back( -1. );
          reco_daughter_allTrack_momByRange_proton.push_back(-999.);
          reco_daughter_allTrack_momByRange_muon.push_back(-999.);


        }
      }
      catch( const cet::exception &e ){
        std::cout << "pandora2Track object not found, moving on" << std::endl;
      }

        
      try{
        const recob::Shower* pandora2Shower = pfpUtil.GetPFParticleShower( *daughterPFP, evt, fPFParticleTag, "pandora2Shower" );
        std::cout << "pandora2 shower: " << pandora2Shower << std::endl;

        if( pandora2Shower ){
          reco_daughter_allShower_ID.push_back(     pandora2Shower->ID() );
          reco_daughter_allShower_len.push_back(    pandora2Shower->Length() );
          reco_daughter_allShower_startX.push_back( pandora2Shower->ShowerStart().X() );
          reco_daughter_allShower_startY.push_back( pandora2Shower->ShowerStart().Y() );
          reco_daughter_allShower_startZ.push_back( pandora2Shower->ShowerStart().Z() );
        }
        else{
          reco_daughter_allShower_ID.push_back(       -1  );
          reco_daughter_allShower_len.push_back(    -999. );
          reco_daughter_allShower_startX.push_back( -999. );
          reco_daughter_allShower_startY.push_back( -999. );
          reco_daughter_allShower_startZ.push_back( -999. );
        }

      }
      catch( const cet::exception &e ){
        std::cout << "pandora2Shower object not found, moving on" << std::endl;
      }
        
      
    }



    if( fCheckCosmics ){
      if( quality_reco_view_2_wire.size() ){

        std::map< int, std::pair< double, double > > UpperCosmicROILimits, LowerCosmicROILimits;
        std::map< int, int > wire_to_n_ticks;
        std::map< int, double > wire_to_avg_ticks;
        bool inTPC1 = false;
        for( size_t i = 0; i < quality_reco_view_2_wire.size(); ++i ){
          if( view_2_TPC[i] != 1 ) continue;
          else inTPC1 = true;
          wire_to_n_ticks[ int(quality_reco_view_2_wire[i]) ]++;
          wire_to_avg_ticks[ int(quality_reco_view_2_wire[i]) ] += quality_reco_view_2_tick[i];

        }

        if( inTPC1 ){

          double max_tick = 0.;
          double min_tick = 99999.;
          for( auto it = wire_to_n_ticks.begin(); it != wire_to_n_ticks.end(); ++it ){
            wire_to_avg_ticks[ it->first ] /= it->second;

            if( wire_to_avg_ticks[ it->first ] > max_tick ){
              max_tick = wire_to_avg_ticks[ it->first ];
            }
            if( wire_to_avg_ticks[ it->first ] < min_tick ){
              min_tick = wire_to_avg_ticks[ it->first ];
            }
          }

          min_tick -= 100;
          max_tick += 100;

          std::vector< double > these_wires, these_ticks;

          for( auto it = wire_to_avg_ticks.begin(); it != wire_to_avg_ticks.end(); ++it ){
            these_wires.push_back( it->first );
            these_ticks.push_back( it->second );

            std::cout << it->first << " " << it->second << std::endl;
          }
          TGraph gr_wire_ticks( these_wires.size(), &these_wires[0], &these_ticks[0] );




          std::cout << "Min tick: " << min_tick << " Max tick: " << max_tick << std::endl;
          std::cout << "Min wire: " << wire_to_avg_ticks.begin()->first << " Max wire: " << wire_to_avg_ticks.rbegin()->first << std::endl;
          

          //1st Get all the reco tracks -- or PFP?
          //for( size_t i = 0; i < recoTracks->size(); ++i ){
          for( size_t i = 0; i < pfpVec->size(); ++i ){

            //if( (*recoTracks)[i].ID() == thisTrack->ID() ) continue;
            if( &(*pfpVec)[i] == particle ) continue; // Check if the same pointer

            int nUpperCosmicROI = 0;
            int nLowerCosmicROI = 0;
            //std::cout << "Checking Track " << (*recoTracks)[i].ID() << std::endl;

            //auto planeHits = trackUtil.GetRecoTrackHitsFromPlane( (*recoTracks)[i], evt, fTrackerTag, 2 );
            auto planeHits = pfpUtil.GetPFParticleHitsFromPlane( (*pfpVec)[i], evt, fPFParticleTag, 2 );
            bool found_new = false;

            protoana::MCParticleSharedHits match = protoana::MCParticleSharedHits();
            if( !evt.isRealData() )
              match = truthUtil.GetMCParticleByHits( (*pfpVec)[i], evt, fPFParticleTag, fHitTag );


            for( size_t j = 0; j < planeHits.size(); ++j ){
              auto theHit = planeHits[j];
              if( theHit->WireID().TPC == 1 ){                    
                
                if( int(theHit->WireID().Wire) > wire_to_avg_ticks.begin()->first && int(theHit->WireID().Wire) < wire_to_avg_ticks.rbegin()->first &&
                    theHit->PeakTime() > min_tick && theHit->PeakTime() < max_tick ){

                  std::cout << "Checking " << theHit->WireID().Wire << " " << theHit->PeakTime() << std::endl;
                  std::cout << "\tBeam: " << gr_wire_ticks.Eval( theHit->WireID().Wire ) << std::endl;
                  
                  if( theHit->PeakTime() > gr_wire_ticks.Eval( theHit->WireID().Wire ) ){
                    ++nUpperCosmicROI; 
                  }
                  else{
                    ++nLowerCosmicROI; 
                  }

                }
                
                if( !found_new && !evt.isRealData() ){
                  if( match.particle ){
                    if( pi_serv->TrackIdToMCTruth_P(match.particle->TrackId())->Origin() == 2 ){
                      std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *theHit );
                      for( size_t j = 0; j < ides.size(); ++j ){
                        if( true_beam_ID == ides[j]->trackID ){
                          //cosmic_has_beam_IDE.push_back( (*recoTracks)[i].ID() );
                          cosmic_has_beam_IDE.push_back( i );
                          found_new = true;
                          break;
                        }
                      }
                    }
                  }
                }
              }
            }

            n_cosmics_with_beam_IDE = cosmic_has_beam_IDE.size();
            
            std::cout << "NHits in Upper ROI: " << nUpperCosmicROI << std::endl;
            std::cout << "NHits in Lower ROI: " << nLowerCosmicROI << std::endl;

            if( nLowerCosmicROI || nUpperCosmicROI ){
              reco_beam_cosmic_candidate_lower_hits.push_back( nLowerCosmicROI );
              reco_beam_cosmic_candidate_upper_hits.push_back( nUpperCosmicROI );
              reco_beam_cosmic_candidate_ID.push_back( i );
            }

          }
        }
      }

      if( !evt.isRealData() ){
        std::cout << "Checking beam for cosmic" << std::endl;
        auto planeHits = trackUtil.GetRecoTrackHitsFromPlane( *thisTrack, evt, fTrackerTag, 2 );
        for( size_t i = 0; i < planeHits.size(); ++i ){      
          auto theHit = planeHits[i];
          //std::cout << theHit->WireID().TPC << std::endl;
          if( theHit->WireID().TPC == 1 ){
            std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *theHit );
            for( size_t j = 0; j < ides.size(); ++j ){
              if( pi_serv->TrackIdToMCTruth_P( ides[j]->trackID )->Origin() == 2 ){
                beam_has_cosmic_IDE = true;
                break;
              }
            }
          }
          if( beam_has_cosmic_IDE ) break;
        }
      }
    }

  }
  else if( thisShower ){
    std::cout << "Beam particle is shower-like" << std::endl;
    reco_beam_type = 11;
    reco_beam_trackID = thisShower->ID();
    std::cout << thisShower->ShowerStart().X() << " " << thisShower->ShowerStart().Y() << " " << thisShower->ShowerStart().Z() << std::endl;
    std::cout << thisShower->Direction().X() << " " << thisShower->Direction().Y() << " " << thisShower->Direction().Z() << std::endl;
    std::cout << beam_cuts.IsBeamlike( *thisShower, evt, "1" ) << std::endl;
  }

  //Forced tracking for beam particle
  try{
    const recob::Track* pandora2Track = pfpUtil.GetPFParticleTrack( *particle, evt, fPFParticleTag, "pandora2Track" );
    std::cout << "pandora2 track: " << pandora2Track << std::endl;


    if( pandora2Track ){
      reco_beam_allTrack_ID = pandora2Track->ID();
      reco_beam_allTrack_beam_cuts = beam_cuts.IsBeamlike( *pandora2Track, evt, "1" );
      reco_beam_allTrack_startX = pandora2Track->Trajectory().Start().X();
      reco_beam_allTrack_startY = pandora2Track->Trajectory().Start().Y();
      reco_beam_allTrack_startZ = pandora2Track->Trajectory().Start().Z();
      reco_beam_allTrack_endX = pandora2Track->Trajectory().End().X();
      reco_beam_allTrack_endY = pandora2Track->Trajectory().End().Y();
      reco_beam_allTrack_endZ = pandora2Track->Trajectory().End().Z();

      auto startDir = pandora2Track->StartDirection();
      auto endDir   = pandora2Track->EndDirection();

      //try flipping
      if( reco_beam_allTrack_startZ > reco_beam_endZ ){
        reco_beam_allTrack_flipped = true;
        reco_beam_allTrack_endX = pandora2Track->Trajectory().Start().X();
        reco_beam_allTrack_endY = pandora2Track->Trajectory().Start().Y();
        reco_beam_allTrack_endZ = pandora2Track->Trajectory().Start().Z();
        reco_beam_allTrack_startX = pandora2Track->Trajectory().End().X();
        reco_beam_allTrack_startY = pandora2Track->Trajectory().End().Y();
        reco_beam_allTrack_startZ = pandora2Track->Trajectory().End().Z();
        
        reco_beam_allTrack_trackDirX =  -1. * endDir.X(); 
        reco_beam_allTrack_trackDirY =  -1. * endDir.Y(); 
        reco_beam_allTrack_trackDirZ =  -1. * endDir.Z(); 

        reco_beam_allTrack_trackEndDirX =  -1. * startDir.X(); 
        reco_beam_allTrack_trackEndDirY =  -1. * startDir.Y(); 
        reco_beam_allTrack_trackEndDirZ =  -1. * startDir.Z(); 
      }
      else{
        reco_beam_allTrack_flipped = false;
        reco_beam_allTrack_trackDirX    =  startDir.X(); 
        reco_beam_allTrack_trackDirY    =  startDir.Y(); 
        reco_beam_allTrack_trackDirZ    =  startDir.Z(); 
        reco_beam_allTrack_trackEndDirX =  endDir.X(); 
        reco_beam_allTrack_trackEndDirY =  endDir.Y(); 
        reco_beam_allTrack_trackEndDirZ =  endDir.Z(); 
      }

      reco_beam_allTrack_len  = pandora2Track->Length();    

      std::vector< anab::Calorimetry> calo = trackUtil.GetRecoTrackCalorimetry(*pandora2Track, evt, fTrackerTag, fCalorimetryTag);
      if( calo.size() ){
        auto calo_range = calo[0].ResidualRange();
        for( size_t i = 0; i < calo_range.size(); ++i ){
          reco_beam_allTrack_resRange.push_back( calo_range[i] );
        }

        //New Calibration
        std::vector< float > new_dEdX = calibration.GetCalibratedCalorimetry(  *pandora2Track, evt, fTrackerTag, fCalorimetryTag );
        for( size_t i = 0; i < new_dEdX.size(); ++i ){ reco_beam_allTrack_calibrated_dEdX.push_back( new_dEdX[i] ); }
        ////////////////////////////////////////////

        std::pair< double, int > pid_chi2_ndof = trackUtil.Chi2PID( reco_beam_allTrack_calibrated_dEdX, reco_beam_allTrack_resRange, templates[ 2212 ] );
        reco_beam_allTrack_Chi2_proton = pid_chi2_ndof.first; 
        reco_beam_allTrack_Chi2_ndof = pid_chi2_ndof.second;
      }
      else{
        reco_beam_allTrack_Chi2_proton = -999;
        reco_beam_allTrack_Chi2_ndof = -999;
      }
    }
    else{
      reco_beam_allTrack_ID = -999;
      reco_beam_allTrack_beam_cuts = -999;
      reco_beam_allTrack_startX = -999;
      reco_beam_allTrack_startY = -999;
      reco_beam_allTrack_startZ = -999;
      reco_beam_allTrack_endX = -999;
      reco_beam_allTrack_endY = -999;
      reco_beam_allTrack_endZ = -999;
      reco_beam_allTrack_Chi2_proton = -999;
      reco_beam_allTrack_Chi2_ndof = -999;
      reco_beam_allTrack_trackDirX    =  -999; 
      reco_beam_allTrack_trackDirY    =  -999; 
      reco_beam_allTrack_trackDirZ    =  -999; 
      reco_beam_allTrack_trackEndDirX =  -999; 
      reco_beam_allTrack_trackEndDirY =  -999; 
      reco_beam_allTrack_trackEndDirZ =  -999; 

    }
  }
  catch( const cet::exception &e ){
    std::cout << "beam pandora2Track object not found, moving on" << std::endl;
  }

  fTree->Fill();
}

void pionana::PionAnalyzer::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("beamana","beam analysis tree");

  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("MC", &MC);



  ///Reconstructed info
  fTree->Branch("reco_beam_type", &reco_beam_type);
  fTree->Branch("reco_beam_startX", &reco_beam_startX);
  fTree->Branch("reco_beam_startY", &reco_beam_startY);
  fTree->Branch("reco_beam_startZ", &reco_beam_startZ);
  fTree->Branch("reco_beam_endX", &reco_beam_endX);
  fTree->Branch("reco_beam_endY", &reco_beam_endY);
  fTree->Branch("reco_beam_endZ", &reco_beam_endZ);
  fTree->Branch("reco_beam_len", &reco_beam_len);
  fTree->Branch("reco_beam_trackDirX", &reco_beam_trackDirX);
  fTree->Branch("reco_beam_trackDirY", &reco_beam_trackDirY);
  fTree->Branch("reco_beam_trackDirZ", &reco_beam_trackDirZ);
  fTree->Branch("reco_beam_trackEndDirX", &reco_beam_trackEndDirX);
  fTree->Branch("reco_beam_trackEndDirY", &reco_beam_trackEndDirY);
  fTree->Branch("reco_beam_trackEndDirZ", &reco_beam_trackEndDirZ);
  fTree->Branch("reco_beam_vtxX", &reco_beam_vtxX);
  fTree->Branch("reco_beam_vtxY", &reco_beam_vtxY);
  fTree->Branch("reco_beam_vtxZ", &reco_beam_vtxZ);
  fTree->Branch("reco_beam_trackID", &reco_beam_trackID);
  fTree->Branch("reco_beam_dQdX", &reco_beam_dQdX);
  fTree->Branch("reco_beam_dEdX", &reco_beam_dEdX);
  fTree->Branch("reco_beam_calibrated_dEdX", &reco_beam_calibrated_dEdX);
  fTree->Branch("reco_beam_resRange", &reco_beam_resRange);
  fTree->Branch("reco_beam_TrkPitch", &reco_beam_TrkPitch);
  fTree->Branch("reco_beam_calo_wire", &reco_beam_calo_wire);
  fTree->Branch("reco_beam_calo_tick", &reco_beam_calo_tick);
  fTree->Branch("reco_beam_nTrackDaughters", &reco_beam_nTrackDaughters);
  fTree->Branch("reco_beam_nShowerDaughters", &reco_beam_nShowerDaughters);
  fTree->Branch("reco_beam_flipped", &reco_beam_flipped);
  fTree->Branch("reco_beam_passes_beam_cuts", &reco_beam_passes_beam_cuts);

  fTree->Branch("reco_beam_PFP_ID", &reco_beam_PFP_ID);
  fTree->Branch("reco_beam_PFP_nHits", &reco_beam_PFP_nHits);
  fTree->Branch("reco_beam_PFP_trackScore", &reco_beam_PFP_trackScore);
  fTree->Branch("reco_beam_PFP_emScore", &reco_beam_PFP_emScore);
  fTree->Branch("reco_beam_PFP_michelScore", &reco_beam_PFP_michelScore);
  fTree->Branch("reco_beam_PFP_trackScore_collection", &reco_beam_PFP_trackScore_collection);
  fTree->Branch("reco_beam_PFP_emScore_collection", &reco_beam_PFP_emScore_collection);
  fTree->Branch("reco_beam_PFP_michelScore_collection", &reco_beam_PFP_michelScore_collection);

  fTree->Branch("reco_beam_allTrack_ID",              &reco_beam_allTrack_ID);
  fTree->Branch("reco_beam_allTrack_beam_cuts",       &reco_beam_allTrack_beam_cuts);
  fTree->Branch("reco_beam_allTrack_flipped",         &reco_beam_allTrack_flipped);
  fTree->Branch("reco_beam_allTrack_len",             &reco_beam_allTrack_len);
  fTree->Branch("reco_beam_allTrack_startX",          &reco_beam_allTrack_startX);
  fTree->Branch("reco_beam_allTrack_startY",          &reco_beam_allTrack_startY);
  fTree->Branch("reco_beam_allTrack_startZ",          &reco_beam_allTrack_startZ);
  fTree->Branch("reco_beam_allTrack_endX",            &reco_beam_allTrack_endX);
  fTree->Branch("reco_beam_allTrack_endY",            &reco_beam_allTrack_endY);
  fTree->Branch("reco_beam_allTrack_endZ",            &reco_beam_allTrack_endZ);
  fTree->Branch("reco_beam_allTrack_trackDirX",       &reco_beam_allTrack_trackDirX);
  fTree->Branch("reco_beam_allTrack_trackDirY",       &reco_beam_allTrack_trackDirY);
  fTree->Branch("reco_beam_allTrack_trackDirZ",       &reco_beam_allTrack_trackDirZ);
  fTree->Branch("reco_beam_allTrack_trackEndDirX",    &reco_beam_allTrack_trackEndDirX);
  fTree->Branch("reco_beam_allTrack_trackEndDirY",    &reco_beam_allTrack_trackEndDirY);
  fTree->Branch("reco_beam_allTrack_trackEndDirZ",    &reco_beam_allTrack_trackEndDirZ);
  fTree->Branch("reco_beam_allTrack_resRange",        &reco_beam_allTrack_resRange);
  fTree->Branch("reco_beam_allTrack_calibrated_dEdX", &reco_beam_allTrack_calibrated_dEdX);
  fTree->Branch("reco_beam_allTrack_Chi2_proton",     &reco_beam_allTrack_Chi2_proton);
  fTree->Branch("reco_beam_allTrack_Chi2_ndof",       &reco_beam_allTrack_Chi2_ndof);


  //Reconstructed info -- daughters
  /*
  fTree->Branch("reco_daughter_trackID", &reco_daughter_trackID);
  fTree->Branch("reco_daughter_true_byE_completeness", &reco_daughter_true_byE_completeness);
  fTree->Branch("reco_daughter_true_byE_purity", &reco_daughter_true_byE_purity);
  fTree->Branch("reco_daughter_true_byE_PDG", &reco_daughter_true_byE_PDG);
  fTree->Branch("reco_daughter_true_byE_ID", &reco_daughter_true_byE_ID);
  fTree->Branch("reco_daughter_true_byE_origin", &reco_daughter_true_byE_origin);
  fTree->Branch("reco_daughter_true_byE_parID", &reco_daughter_true_byE_parID);
  fTree->Branch("reco_daughter_true_byE_parPDG", &reco_daughter_true_byE_parPDG);
  fTree->Branch("reco_daughter_true_byE_process", &reco_daughter_true_byE_process);

  fTree->Branch("reco_daughter_true_byHits_PDG", &reco_daughter_true_byHits_PDG);
  fTree->Branch("reco_daughter_true_byHits_ID", &reco_daughter_true_byHits_ID);
  fTree->Branch("reco_daughter_true_byHits_origin", &reco_daughter_true_byHits_origin);
  fTree->Branch("reco_daughter_true_byHits_parID", &reco_daughter_true_byHits_parID);
  fTree->Branch("reco_daughter_true_byHits_parPDG", &reco_daughter_true_byHits_parPDG);
  fTree->Branch("reco_daughter_true_byHits_process", &reco_daughter_true_byHits_process);
  fTree->Branch("reco_daughter_true_byHits_purity", &reco_daughter_true_byHits_purity);
  fTree->Branch("reco_daughter_true_byHits_sharedHits", &reco_daughter_true_byHits_sharedHits);
  fTree->Branch("reco_daughter_true_byHits_emHits", &reco_daughter_true_byHits_emHits);

  fTree->Branch("reco_daughter_true_byHits_len", &reco_daughter_true_byHits_len);
  fTree->Branch("reco_daughter_true_byHits_startX", &reco_daughter_true_byHits_startX);
  fTree->Branch("reco_daughter_true_byHits_startY", &reco_daughter_true_byHits_startY);
  fTree->Branch("reco_daughter_true_byHits_startZ", &reco_daughter_true_byHits_startZ);
  fTree->Branch("reco_daughter_true_byHits_endX", &reco_daughter_true_byHits_endX);
  fTree->Branch("reco_daughter_true_byHits_endY", &reco_daughter_true_byHits_endY);
  fTree->Branch("reco_daughter_true_byHits_endZ", &reco_daughter_true_byHits_endZ);

  fTree->Branch("reco_daughter_true_byHits_startPx", &reco_daughter_true_byHits_startPx);
  fTree->Branch("reco_daughter_true_byHits_startPy", &reco_daughter_true_byHits_startPy);
  fTree->Branch("reco_daughter_true_byHits_startPz", &reco_daughter_true_byHits_startPz);
  fTree->Branch("reco_daughter_true_byHits_startP", &reco_daughter_true_byHits_startP);
  fTree->Branch("reco_daughter_true_byHits_startE", &reco_daughter_true_byHits_startE);
  */
  //Alternative reco
  fTree->Branch("reco_daughter_PFP_true_byHits_PDG", &reco_daughter_PFP_true_byHits_PDG);
  fTree->Branch("reco_daughter_PFP_true_byHits_ID", &reco_daughter_PFP_true_byHits_ID);
  fTree->Branch("reco_daughter_PFP_true_byHits_origin", &reco_daughter_PFP_true_byHits_origin);
  fTree->Branch("reco_daughter_PFP_true_byHits_parID", &reco_daughter_PFP_true_byHits_parID);
  fTree->Branch("reco_daughter_PFP_true_byHits_parPDG", &reco_daughter_PFP_true_byHits_parPDG);
  fTree->Branch("reco_daughter_PFP_true_byHits_process", &reco_daughter_PFP_true_byHits_process);
  fTree->Branch("reco_daughter_PFP_true_byHits_sharedHits", &reco_daughter_PFP_true_byHits_sharedHits);
  fTree->Branch("reco_daughter_PFP_true_byHits_emHits", &reco_daughter_PFP_true_byHits_emHits);

  fTree->Branch("reco_daughter_PFP_true_byHits_len", &reco_daughter_PFP_true_byHits_len);
  fTree->Branch("reco_daughter_PFP_true_byHits_startX", &reco_daughter_PFP_true_byHits_startX);
  fTree->Branch("reco_daughter_PFP_true_byHits_startY", &reco_daughter_PFP_true_byHits_startY);
  fTree->Branch("reco_daughter_PFP_true_byHits_startZ", &reco_daughter_PFP_true_byHits_startZ);
  fTree->Branch("reco_daughter_PFP_true_byHits_endX", &reco_daughter_PFP_true_byHits_endX);
  fTree->Branch("reco_daughter_PFP_true_byHits_endY", &reco_daughter_PFP_true_byHits_endY);
  fTree->Branch("reco_daughter_PFP_true_byHits_endZ", &reco_daughter_PFP_true_byHits_endZ);

  fTree->Branch("reco_daughter_PFP_true_byHits_startPx", &reco_daughter_PFP_true_byHits_startPx);
  fTree->Branch("reco_daughter_PFP_true_byHits_startPy", &reco_daughter_PFP_true_byHits_startPy);
  fTree->Branch("reco_daughter_PFP_true_byHits_startPz", &reco_daughter_PFP_true_byHits_startPz);
  fTree->Branch("reco_daughter_PFP_true_byHits_startP", &reco_daughter_PFP_true_byHits_startP);
  fTree->Branch("reco_daughter_PFP_true_byHits_startE", &reco_daughter_PFP_true_byHits_startE);
  fTree->Branch("reco_daughter_PFP_true_byHits_endProcess", &reco_daughter_PFP_true_byHits_endProcess);
  fTree->Branch("reco_daughter_PFP_true_byHits_purity", &reco_daughter_PFP_true_byHits_purity);

  fTree->Branch("reco_daughter_allTrack_ID", &reco_daughter_allTrack_ID);
  fTree->Branch("reco_daughter_allTrack_dEdX", &reco_daughter_allTrack_dEdX);
  fTree->Branch("reco_daughter_allTrack_dQdX", &reco_daughter_allTrack_dQdX);
  fTree->Branch("reco_daughter_allTrack_resRange", &reco_daughter_allTrack_resRange);
  fTree->Branch("reco_daughter_allTrack_dQdX_SCE", &reco_daughter_allTrack_dQdX_SCE);
  fTree->Branch("reco_daughter_allTrack_dEdX_SCE", &reco_daughter_allTrack_dEdX);
  fTree->Branch("reco_daughter_allTrack_resRange_SCE", &reco_daughter_allTrack_resRange);

  fTree->Branch("reco_daughter_allTrack_calibrated_dEdX", &reco_daughter_allTrack_calibrated_dEdX);
  fTree->Branch("reco_daughter_allTrack_calibrated_dEdX_SCE", &reco_daughter_allTrack_calibrated_dEdX_SCE);

  fTree->Branch("reco_daughter_allTrack_Chi2_proton", &reco_daughter_allTrack_Chi2_proton);
  fTree->Branch("reco_daughter_allTrack_Chi2_ndof", &reco_daughter_allTrack_Chi2_ndof);

  fTree->Branch("reco_daughter_allTrack_Theta", &reco_daughter_allTrack_Theta);
  fTree->Branch("reco_daughter_allTrack_Phi", &reco_daughter_allTrack_Phi);

  fTree->Branch("reco_daughter_allTrack_len", &reco_daughter_allTrack_len);
  fTree->Branch("reco_daughter_allTrack_startX", &reco_daughter_allTrack_startX);
  fTree->Branch("reco_daughter_allTrack_startY", &reco_daughter_allTrack_startY);
  fTree->Branch("reco_daughter_allTrack_startZ", &reco_daughter_allTrack_startZ);
  fTree->Branch("reco_daughter_allTrack_endX", &reco_daughter_allTrack_endX);
  fTree->Branch("reco_daughter_allTrack_endY", &reco_daughter_allTrack_endY);
  fTree->Branch("reco_daughter_allTrack_endZ", &reco_daughter_allTrack_endZ);
  fTree->Branch("reco_daughter_allTrack_dR", &reco_daughter_allTrack_dR);
  fTree->Branch("reco_daughter_allTrack_to_vertex", &reco_daughter_allTrack_to_vertex);
  //////

  fTree->Branch("reco_daughter_allShower_ID", &reco_daughter_allShower_ID);
  fTree->Branch("reco_daughter_allShower_len", &reco_daughter_allShower_len);
  fTree->Branch("reco_daughter_allShower_startX", &reco_daughter_allShower_startX);
  fTree->Branch("reco_daughter_allShower_startY", &reco_daughter_allShower_startY);
  fTree->Branch("reco_daughter_allShower_startZ", &reco_daughter_allShower_startZ);


/*
  fTree->Branch("reco_daughter_shower_true_byE_PDG", &reco_daughter_shower_true_byE_PDG);
  fTree->Branch("reco_daughter_shower_true_byE_ID", &reco_daughter_shower_true_byE_ID);
  fTree->Branch("reco_daughter_shower_true_byE_origin", &reco_daughter_shower_true_byE_origin);
  fTree->Branch("reco_daughter_shower_true_byE_parID", &reco_daughter_shower_true_byE_parID);
  fTree->Branch("reco_daughter_shower_true_byE_parPDG", &reco_daughter_shower_true_byE_parPDG);
 
  fTree->Branch("reco_daughter_shower_true_byE_startPx", &reco_daughter_shower_true_byE_startPx);
  fTree->Branch("reco_daughter_shower_true_byE_startPy", &reco_daughter_shower_true_byE_startPy);
  fTree->Branch("reco_daughter_shower_true_byE_startPz", &reco_daughter_shower_true_byE_startPz);
  fTree->Branch("reco_daughter_shower_true_byE_startP", &reco_daughter_shower_true_byE_startP);
  fTree->Branch("reco_daughter_shower_true_byE_endProcess", &reco_daughter_shower_true_byE_endProcess);


  fTree->Branch("reco_daughter_shower_true_byHits_PDG", &reco_daughter_shower_true_byHits_PDG);
  fTree->Branch("reco_daughter_shower_true_byHits_ID", &reco_daughter_shower_true_byHits_ID);
  fTree->Branch("reco_daughter_shower_true_byHits_origin", &reco_daughter_shower_true_byHits_origin);
  fTree->Branch("reco_daughter_shower_true_byHits_parID", &reco_daughter_shower_true_byHits_parID);
  fTree->Branch("reco_daughter_shower_true_byHits_parPDG", &reco_daughter_shower_true_byHits_parPDG);
  fTree->Branch("reco_daughter_shower_true_byHits_process", &reco_daughter_shower_true_byHits_process);
  fTree->Branch("reco_daughter_shower_true_byHits_purity", &reco_daughter_shower_true_byHits_purity);

  fTree->Branch("reco_daughter_shower_true_byHits_startPx", &reco_daughter_shower_true_byHits_startPx);
  fTree->Branch("reco_daughter_shower_true_byHits_startPy", &reco_daughter_shower_true_byHits_startPy);
  fTree->Branch("reco_daughter_shower_true_byHits_startPz", &reco_daughter_shower_true_byHits_startPz);
  fTree->Branch("reco_daughter_shower_true_byHits_startP", &reco_daughter_shower_true_byHits_startP);
  fTree->Branch("reco_daughter_shower_true_byHits_endProcess", &reco_daughter_shower_true_byHits_endProcess);
  */


  ///Reconstructed info -- daughter
  /*
  fTree->Branch("reco_daughter_showerID", &reco_daughter_showerID);
  fTree->Branch("reco_daughter_dQdX", &reco_daughter_dQdX);
  fTree->Branch("reco_daughter_dEdX", &reco_daughter_dEdX);
  fTree->Branch("reco_daughter_resRange", &reco_daughter_resRange);
  fTree->Branch("reco_daughter_shower_dQdX", &reco_daughter_shower_dQdX);
  fTree->Branch("reco_daughter_shower_dEdX", &reco_daughter_shower_dEdX);
  fTree->Branch("reco_daughter_shower_resRange", &reco_daughter_shower_resRange);
  fTree->Branch("reco_daughter_len", &reco_daughter_len);
  fTree->Branch("reco_daughter_startX", &reco_daughter_startX);
  fTree->Branch("reco_daughter_startY", &reco_daughter_startY);
  fTree->Branch("reco_daughter_startZ", &reco_daughter_startZ);
  fTree->Branch("reco_daughter_endX", &reco_daughter_endX);
  fTree->Branch("reco_daughter_endY", &reco_daughter_endY);
  fTree->Branch("reco_daughter_endZ", &reco_daughter_endZ);
  fTree->Branch("reco_daughter_deltaR", &reco_daughter_deltaR);

  fTree->Branch("reco_daughter_dR", &reco_daughter_dR);
  fTree->Branch("reco_daughter_to_vertex", &reco_daughter_to_vertex);
  fTree->Branch("reco_daughter_slice", &reco_daughter_slice);

  fTree->Branch("reco_daughter_shower_to_vertex", &reco_daughter_shower_to_vertex);

  fTree->Branch("reco_daughter_shower_startX", &reco_daughter_shower_startX);
  fTree->Branch("reco_daughter_shower_startY", &reco_daughter_shower_startY);
  fTree->Branch("reco_daughter_shower_startZ", &reco_daughter_shower_startZ);
  fTree->Branch("reco_daughter_shower_len", &reco_daughter_shower_len);
  */


  fTree->Branch("reco_daughter_PFP_ID", &reco_daughter_PFP_ID);
  fTree->Branch("reco_daughter_PFP_nHits", &reco_daughter_PFP_nHits);
  fTree->Branch("reco_daughter_PFP_trackScore", &reco_daughter_PFP_trackScore);
  fTree->Branch("reco_daughter_PFP_emScore", &reco_daughter_PFP_emScore);
  fTree->Branch("reco_daughter_PFP_michelScore", &reco_daughter_PFP_michelScore);
  fTree->Branch("reco_daughter_PFP_trackScore_collection", &reco_daughter_PFP_trackScore_collection);
  fTree->Branch("reco_daughter_PFP_emScore_collection", &reco_daughter_PFP_emScore_collection);
  fTree->Branch("reco_daughter_PFP_michelScore_collection", &reco_daughter_PFP_michelScore_collection);




  fTree->Branch("true_beam_PDG", &true_beam_PDG);
  fTree->Branch("true_beam_ID", &true_beam_ID);
  fTree->Branch("true_beam_endProcess", &true_beam_endProcess);
  fTree->Branch("true_beam_endX", &true_beam_endX);
  fTree->Branch("true_beam_endY", &true_beam_endY);
  fTree->Branch("true_beam_endZ", &true_beam_endZ);
  fTree->Branch("true_beam_startX", &true_beam_startX);
  fTree->Branch("true_beam_startY", &true_beam_startY);
  fTree->Branch("true_beam_startZ", &true_beam_startZ);

  fTree->Branch("true_beam_startPx", &true_beam_startPx);
  fTree->Branch("true_beam_startPy", &true_beam_startPy);
  fTree->Branch("true_beam_startPz", &true_beam_startPz);
  fTree->Branch("true_beam_startP", &true_beam_startP);

  fTree->Branch("true_beam_endPx", &true_beam_endPx);
  fTree->Branch("true_beam_endPy", &true_beam_endPy);
  fTree->Branch("true_beam_endPz", &true_beam_endPz);
  fTree->Branch("true_beam_endP", &true_beam_endP);

  fTree->Branch("true_beam_startDirX", &true_beam_startDirX);
  fTree->Branch("true_beam_startDirY", &true_beam_startDirY);
  fTree->Branch("true_beam_startDirZ", &true_beam_startDirZ);

  fTree->Branch("true_beam_nElasticScatters", &true_beam_nElasticScatters);
  fTree->Branch("true_beam_elastic_costheta", &true_beam_elastic_costheta);
  fTree->Branch("true_beam_elastic_X", &true_beam_elastic_X);
  fTree->Branch("true_beam_elastic_Y", &true_beam_elastic_Y);
  fTree->Branch("true_beam_elastic_Z", &true_beam_elastic_Z);
  fTree->Branch("true_beam_IDE_totalDep",    &true_beam_IDE_totalDep);
  fTree->Branch("true_beam_IDE_found_in_recoVtx",    &true_beam_IDE_found_in_recoVtx);

  fTree->Branch("true_beam_nHits", &true_beam_nHits);
  fTree->Branch("true_beam_reco_byHits_PFP_ID", &true_beam_reco_byHits_PFP_ID);
  fTree->Branch("true_beam_reco_byHits_PFP_nHits", &true_beam_reco_byHits_PFP_nHits);
  fTree->Branch("true_beam_reco_byHits_allTrack_ID", &true_beam_reco_byHits_allTrack_ID);

  fTree->Branch("true_daughter_nPi0", &true_daughter_nPi0);
  fTree->Branch("true_daughter_nPiPlus", &true_daughter_nPiPlus);
  fTree->Branch("true_daughter_nProton", &true_daughter_nProton);
  fTree->Branch("true_daughter_nNeutron", &true_daughter_nNeutron);
  fTree->Branch("true_daughter_nPiMinus", &true_daughter_nPiMinus);
  fTree->Branch("true_daughter_nNucleus", &true_daughter_nNucleus);

  fTree->Branch("reco_beam_vertex_slice", &reco_beam_vertex_slice);


  fTree->Branch("reco_beam_vertex_dRs", &reco_beam_vertex_dRs);
  fTree->Branch("reco_beam_vertex_hits_slices", &reco_beam_vertex_hits_slices);

  fTree->Branch("true_beam_daughter_PDG", &true_beam_daughter_PDG);
  fTree->Branch("true_beam_daughter_ID", &true_beam_daughter_ID);
  fTree->Branch("true_beam_daughter_len", &true_beam_daughter_len);
  fTree->Branch("true_beam_daughter_startX", &true_beam_daughter_startX);
  fTree->Branch("true_beam_daughter_startY", &true_beam_daughter_startY);
  fTree->Branch("true_beam_daughter_startZ", &true_beam_daughter_startZ);
  fTree->Branch("true_beam_daughter_startPx", &true_beam_daughter_startPx);
  fTree->Branch("true_beam_daughter_startPy", &true_beam_daughter_startPy);
  fTree->Branch("true_beam_daughter_startPz", &true_beam_daughter_startPz);
  fTree->Branch("true_beam_daughter_startP", &true_beam_daughter_startP);
  fTree->Branch("true_beam_daughter_endX", &true_beam_daughter_endX);
  fTree->Branch("true_beam_daughter_endY", &true_beam_daughter_endY);
  fTree->Branch("true_beam_daughter_endZ", &true_beam_daughter_endZ);
  fTree->Branch("true_beam_daughter_Process", &true_beam_daughter_Process);
  fTree->Branch("true_beam_daughter_endProcess", &true_beam_daughter_endProcess);
  fTree->Branch("true_beam_daughter_nHits", &true_beam_daughter_nHits);

  fTree->Branch("true_beam_daughter_reco_byHits_PFP_ID", &true_beam_daughter_reco_byHits_PFP_ID);
  fTree->Branch("true_beam_daughter_reco_byHits_PFP_nHits", &true_beam_daughter_reco_byHits_PFP_nHits);
  fTree->Branch("true_beam_daughter_reco_byHits_PFP_trackScore", &true_beam_daughter_reco_byHits_PFP_trackScore);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_ID", &true_beam_daughter_reco_byHits_allTrack_ID);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_startX", &true_beam_daughter_reco_byHits_allTrack_startX);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_startY", &true_beam_daughter_reco_byHits_allTrack_startY);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_startZ", &true_beam_daughter_reco_byHits_allTrack_startZ);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_endX", &true_beam_daughter_reco_byHits_allTrack_endX);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_endY", &true_beam_daughter_reco_byHits_allTrack_endY);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_endZ", &true_beam_daughter_reco_byHits_allTrack_endZ);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_len", &true_beam_daughter_reco_byHits_allTrack_len);

  fTree->Branch("true_beam_daughter_reco_byHits_allShower_ID", &true_beam_daughter_reco_byHits_allShower_ID);
  fTree->Branch("true_beam_daughter_reco_byHits_allShower_startX", &true_beam_daughter_reco_byHits_allShower_startX);
  fTree->Branch("true_beam_daughter_reco_byHits_allShower_startY", &true_beam_daughter_reco_byHits_allShower_startY);
  fTree->Branch("true_beam_daughter_reco_byHits_allShower_startZ", &true_beam_daughter_reco_byHits_allShower_startZ);
  fTree->Branch("true_beam_daughter_reco_byHits_allShower_len", &true_beam_daughter_reco_byHits_allShower_len);

  fTree->Branch("true_beam_Pi0_decay_ID", &true_beam_Pi0_decay_ID);
  fTree->Branch("true_beam_Pi0_decay_parID", &true_beam_Pi0_decay_parID);
  fTree->Branch("true_beam_Pi0_decay_PDG", &true_beam_Pi0_decay_PDG);
  fTree->Branch("true_beam_Pi0_decay_startP", &true_beam_Pi0_decay_startP);
  fTree->Branch("true_beam_Pi0_decay_len", &true_beam_Pi0_decay_len);
  fTree->Branch("true_beam_Pi0_decay_nHits", &true_beam_Pi0_decay_nHits);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_PFP_ID", &true_beam_Pi0_decay_reco_byHits_PFP_ID);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_PFP_nHits", &true_beam_Pi0_decay_reco_byHits_PFP_nHits);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_PFP_trackScore", &true_beam_Pi0_decay_reco_byHits_PFP_trackScore);

  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_ID", &true_beam_Pi0_decay_reco_byHits_allTrack_ID);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_startX", &true_beam_Pi0_decay_reco_byHits_allTrack_startX);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_startY", &true_beam_Pi0_decay_reco_byHits_allTrack_startY);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_startZ", &true_beam_Pi0_decay_reco_byHits_allTrack_startZ);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_endX", &true_beam_Pi0_decay_reco_byHits_allTrack_endX);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_endY", &true_beam_Pi0_decay_reco_byHits_allTrack_endY);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_endZ", &true_beam_Pi0_decay_reco_byHits_allTrack_endZ);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_len", &true_beam_Pi0_decay_reco_byHits_allTrack_len);

  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allShower_ID", &true_beam_Pi0_decay_reco_byHits_allShower_ID);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allShower_startX", &true_beam_Pi0_decay_reco_byHits_allShower_startX);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allShower_startY", &true_beam_Pi0_decay_reco_byHits_allShower_startY);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allShower_startZ", &true_beam_Pi0_decay_reco_byHits_allShower_startZ);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allShower_len", &true_beam_Pi0_decay_reco_byHits_allShower_len);

  fTree->Branch("true_beam_grand_daughter_ID", &true_beam_grand_daughter_ID);
  fTree->Branch("true_beam_grand_daughter_parID", &true_beam_grand_daughter_parID);
  fTree->Branch("true_beam_grand_daughter_PDG", &true_beam_grand_daughter_PDG);
  fTree->Branch("true_beam_grand_daughter_nHits", &true_beam_grand_daughter_nHits);
  fTree->Branch("true_beam_grand_daughter_Process", &true_beam_grand_daughter_Process);
  fTree->Branch("true_beam_grand_daughter_endProcess", &true_beam_grand_daughter_endProcess);

  ////Matching reco to truth
  fTree->Branch("reco_beam_true_byE_endProcess", &reco_beam_true_byE_endProcess);
  fTree->Branch("reco_beam_true_byE_process", &reco_beam_true_byE_process);
  fTree->Branch("reco_beam_true_byE_origin", &reco_beam_true_byE_origin);
  fTree->Branch("reco_beam_true_byE_PDG", &reco_beam_true_byE_PDG);
  fTree->Branch("reco_beam_true_byE_ID", &reco_beam_true_byE_ID);

  fTree->Branch("reco_beam_true_byHits_endProcess", &reco_beam_true_byHits_endProcess);
  fTree->Branch("reco_beam_true_byHits_process", &reco_beam_true_byHits_process);
  fTree->Branch("reco_beam_true_byHits_origin", &reco_beam_true_byHits_origin);
  fTree->Branch("reco_beam_true_byHits_PDG", &reco_beam_true_byHits_PDG);
  fTree->Branch("reco_beam_true_byHits_ID", &reco_beam_true_byHits_ID);

  fTree->Branch("reco_beam_true_byE_matched", &reco_beam_true_byE_matched);
  fTree->Branch("reco_beam_true_byHits_matched", &reco_beam_true_byHits_matched);
  fTree->Branch("reco_beam_true_byHits_purity", &reco_beam_true_byHits_purity);

  fTree->Branch("true_beam_processes", &true_beam_processes);
  //fTree->Branch("reco_daughter_true_byE_isPrimary", &reco_daughter_true_byE_isPrimary);

  fTree->Branch("data_BI_P", &data_BI_P);
  fTree->Branch("data_BI_X", &data_BI_X);
  fTree->Branch("data_BI_Y", &data_BI_Y);
  fTree->Branch("data_BI_Z", &data_BI_Z);
  fTree->Branch("data_BI_dirX", &data_BI_dirX);
  fTree->Branch("data_BI_dirY", &data_BI_dirY);
  fTree->Branch("data_BI_dirZ", &data_BI_dirZ);

  fTree->Branch("data_BI_nFibersP1", &data_BI_nFibersP1);
  fTree->Branch("data_BI_nFibersP2", &data_BI_nFibersP2);
  fTree->Branch("data_BI_nFibersP3", &data_BI_nFibersP3);
  fTree->Branch("data_BI_PDG_candidates", &data_BI_PDG_candidates);
  fTree->Branch("data_BI_nTracks", &data_BI_nTracks);
  fTree->Branch("data_BI_nMomenta", &data_BI_nMomenta);


  fTree->Branch("quality_reco_view_0_hits_in_TPC5", &quality_reco_view_0_hits_in_TPC5);
  fTree->Branch("quality_reco_view_1_hits_in_TPC5", &quality_reco_view_1_hits_in_TPC5);
  fTree->Branch("quality_reco_view_2_hits_in_TPC5", &quality_reco_view_2_hits_in_TPC5);
  fTree->Branch("quality_reco_max_lateral", &quality_reco_max_lateral);
  fTree->Branch("quality_reco_max_segment", &quality_reco_max_segment);
  fTree->Branch("quality_reco_view_0_max_segment", &quality_reco_view_0_max_segment);
  fTree->Branch("quality_reco_view_1_max_segment", &quality_reco_view_1_max_segment);
  fTree->Branch("quality_reco_view_2_max_segment", &quality_reco_view_2_max_segment);

  fTree->Branch("quality_reco_view_0_wire_backtrack", &quality_reco_view_0_wire_backtrack);
  fTree->Branch("quality_reco_view_1_wire_backtrack", &quality_reco_view_1_wire_backtrack);
  fTree->Branch("quality_reco_view_2_wire_backtrack", &quality_reco_view_2_wire_backtrack);

  fTree->Branch("quality_reco_view_0_wire", &quality_reco_view_0_wire);
  fTree->Branch("quality_reco_view_1_wire", &quality_reco_view_1_wire);
  fTree->Branch("quality_reco_view_2_wire", &quality_reco_view_2_wire);

  fTree->Branch("quality_reco_view_2_z", &quality_reco_view_2_z);

  fTree->Branch("quality_reco_view_0_tick", &quality_reco_view_0_tick);
  fTree->Branch("quality_reco_view_1_tick", &quality_reco_view_1_tick);
  fTree->Branch("quality_reco_view_2_tick", &quality_reco_view_2_tick);

  fTree->Branch("reco_beam_Chi2_proton", &reco_beam_Chi2_proton);
  fTree->Branch("reco_beam_Chi2_ndof", &reco_beam_Chi2_ndof);

  fTree->Branch("reco_beam_cosmic_candidate_lower_hits", &reco_beam_cosmic_candidate_lower_hits);
  fTree->Branch("reco_beam_cosmic_candidate_upper_hits", &reco_beam_cosmic_candidate_upper_hits);
  fTree->Branch("reco_beam_cosmic_candidate_ID", &reco_beam_cosmic_candidate_ID);
  fTree->Branch("beam_has_cosmic_IDE", &beam_has_cosmic_IDE);
  fTree->Branch("cosmic_has_beam_IDE", &cosmic_has_beam_IDE);
  fTree->Branch("n_cosmics_with_beam_IDE", &n_cosmics_with_beam_IDE);

  /*
  fTree->Branch("reco_daughter_Chi2_proton", &reco_daughter_Chi2_proton);
  fTree->Branch("reco_daughter_Chi2_ndof", &reco_daughter_Chi2_ndof);
  fTree->Branch("reco_daughter_momByRange_proton", &reco_daughter_momByRange_proton);
  fTree->Branch("reco_daughter_momByRange_muon", &reco_daughter_momByRange_muon);
  fTree->Branch("reco_daughter_allTrack_momByRange_proton", &reco_daughter_allTrack_momByRange_proton);
  fTree->Branch("reco_daughter_allTrack_momByRange_muon", &reco_daughter_allTrack_momByRange_muon);

  fTree->Branch("reco_daughter_shower_Chi2_proton", &reco_daughter_shower_Chi2_proton);
  fTree->Branch("reco_daughter_shower_Chi2_ndof", &reco_daughter_shower_Chi2_ndof);

  fTree->Branch("reco_daughter_trackScore", &reco_daughter_trackScore);
  fTree->Branch("reco_daughter_emScore", &reco_daughter_emScore);
  fTree->Branch("reco_daughter_michelScore", &reco_daughter_michelScore);

  fTree->Branch("reco_daughter_shower_trackScore", &reco_daughter_shower_trackScore);
  fTree->Branch("reco_daughter_shower_emScore", &reco_daughter_shower_emScore);
  fTree->Branch("reco_daughter_shower_michelScore", &reco_daughter_shower_michelScore);
  */

  fTree->Branch("reco_beam_true_byE_endPx", &reco_beam_true_byE_endPx);
  fTree->Branch("reco_beam_true_byE_endPy", &reco_beam_true_byE_endPy);
  fTree->Branch("reco_beam_true_byE_endPz", &reco_beam_true_byE_endPz);
  fTree->Branch("reco_beam_true_byE_endE", &reco_beam_true_byE_endE);
  fTree->Branch("reco_beam_true_byE_endP", &reco_beam_true_byE_endP);

  fTree->Branch("reco_beam_true_byE_startPx", &reco_beam_true_byE_startPx);
  fTree->Branch("reco_beam_true_byE_startPy", &reco_beam_true_byE_startPy);
  fTree->Branch("reco_beam_true_byE_startPz", &reco_beam_true_byE_startPz);
  fTree->Branch("reco_beam_true_byE_startE", &reco_beam_true_byE_startE);
  fTree->Branch("reco_beam_true_byE_startP", &reco_beam_true_byE_startP);


  fTree->Branch("reco_beam_true_byHits_endPx", &reco_beam_true_byHits_endPx);
  fTree->Branch("reco_beam_true_byHits_endPy", &reco_beam_true_byHits_endPy);
  fTree->Branch("reco_beam_true_byHits_endPz", &reco_beam_true_byHits_endPz);
  fTree->Branch("reco_beam_true_byHits_endE", &reco_beam_true_byHits_endE);
  fTree->Branch("reco_beam_true_byHits_endP", &reco_beam_true_byHits_endP);

  fTree->Branch("reco_beam_true_byHits_startPx", &reco_beam_true_byHits_startPx);
  fTree->Branch("reco_beam_true_byHits_startPy", &reco_beam_true_byHits_startPy);
  fTree->Branch("reco_beam_true_byHits_startPz", &reco_beam_true_byHits_startPz);
  fTree->Branch("reco_beam_true_byHits_startE", &reco_beam_true_byHits_startE);
  fTree->Branch("reco_beam_true_byHits_startP", &reco_beam_true_byHits_startP);

  fTree->Branch("reco_beam_incidentEnergies", &reco_beam_incidentEnergies);
  fTree->Branch("reco_beam_interactingEnergy", &reco_beam_interactingEnergy);
  fTree->Branch("true_beam_incidentEnergies", &true_beam_incidentEnergies);
  fTree->Branch("true_beam_interactingEnergy", &true_beam_interactingEnergy);

  if( fSaveHits ){
    fTree->Branch( "reco_beam_spacePts_X", &reco_beam_spacePts_X );
    fTree->Branch( "reco_beam_spacePts_Y", &reco_beam_spacePts_Y );
    fTree->Branch( "reco_beam_spacePts_Z", &reco_beam_spacePts_Z );

    fTree->Branch( "reco_daughter_spacePts_X", &reco_daughter_spacePts_X );
    fTree->Branch( "reco_daughter_spacePts_Y", &reco_daughter_spacePts_Y );
    fTree->Branch( "reco_daughter_spacePts_Z", &reco_daughter_spacePts_Z );

    fTree->Branch( "reco_daughter_shower_spacePts_X", &reco_daughter_shower_spacePts_X );
    fTree->Branch( "reco_daughter_shower_spacePts_Y", &reco_daughter_shower_spacePts_Y );
    fTree->Branch( "reco_daughter_shower_spacePts_Z", &reco_daughter_shower_spacePts_Z );
  }

}

void pionana::PionAnalyzer::endJob()
{
  dEdX_template_file.Close();
}

double pionana::PionAnalyzer::lateralDist(TVector3 &n, TVector3 &x0, TVector3 &p){
  TVector3 x = ( (p - x0)*n )*n;
  return (x - (p - x0)).Mag();
}

void pionana::PionAnalyzer::reset()
{
  reco_beam_startX = -1;
  reco_beam_startY = -1;
  reco_beam_startZ = -1;
  reco_beam_endX = -1;
  reco_beam_endY = -1;
  reco_beam_endZ = -1;
  reco_beam_flipped = false;

  reco_beam_len = -1;
  reco_beam_type = -1;
  reco_beam_passes_beam_cuts = false;
  
  reco_beam_vertex_slice = std::numeric_limits<int>::max();
  reco_beam_vertex_dRs.clear();
  reco_beam_vertex_hits_slices.clear();

  true_daughter_nPi0 = 0;
  true_daughter_nPiPlus = 0;
  true_daughter_nPiMinus = 0;
  true_daughter_nProton = 0;
  true_daughter_nNeutron = 0;
  true_daughter_nNucleus = 0;

  reco_beam_true_byE_PDG = 0;
  reco_beam_true_byE_ID = 0;
  reco_beam_true_byHits_PDG = 0;
  reco_beam_true_byHits_ID = 0;

  true_beam_PDG = 0;
  true_beam_ID = 0;
  true_beam_endProcess ="";
  true_beam_endX = 0.;
  true_beam_endY = 0.;
  true_beam_endZ = 0.;
  true_beam_startX = 0.;
  true_beam_startY = 0.;
  true_beam_startZ = 0.;

  true_beam_startPx   = 0.; 
  true_beam_startPy   = 0.; 
  true_beam_startPz   = 0.; 
  true_beam_startP    = 0.; 

  true_beam_endPx   = 0.; 
  true_beam_endPy   = 0.; 
  true_beam_endPz   = 0.; 
  true_beam_endP    = 0.; 

  true_beam_startDirX = 0.; 
  true_beam_startDirY = 0.; 
  true_beam_startDirZ = 0.; 
  true_beam_nHits = -1;


  true_beam_processes.clear();
  true_beam_nElasticScatters = 0;
  true_beam_elastic_costheta.clear();
  true_beam_elastic_X.clear();
  true_beam_elastic_Y.clear();
  true_beam_elastic_Z.clear();
  true_beam_IDE_totalDep = 0.;
  true_beam_IDE_found_in_recoVtx = false;

  true_beam_reco_byHits_PFP_ID.clear();
  true_beam_reco_byHits_PFP_nHits.clear();
  true_beam_reco_byHits_allTrack_ID.clear();

  reco_beam_true_byE_endProcess ="";
  reco_beam_true_byE_process ="";
  reco_beam_true_byE_origin = -1;

  reco_beam_true_byE_endPx = 0.;
  reco_beam_true_byE_endPy = 0.;
  reco_beam_true_byE_endPz = 0.;
  reco_beam_true_byE_endE = 0.;
  reco_beam_true_byE_endP = 0.;

  reco_beam_true_byE_startPx = 0.;
  reco_beam_true_byE_startPy = 0.;
  reco_beam_true_byE_startPz = 0.;
  reco_beam_true_byE_startE = 0.;
  reco_beam_true_byE_startP = 0.;

  reco_beam_true_byHits_endProcess ="";
  reco_beam_true_byHits_process ="";
  reco_beam_true_byHits_origin = -1;

  reco_beam_true_byHits_endPx = 0.;
  reco_beam_true_byHits_endPy = 0.;
  reco_beam_true_byHits_endPz = 0.;
  reco_beam_true_byHits_endE = 0.;
  reco_beam_true_byHits_endP = 0.;

  reco_beam_true_byHits_startPx = 0.;
  reco_beam_true_byHits_startPy = 0.;
  reco_beam_true_byHits_startPz = 0.;
  reco_beam_true_byHits_startE = 0.;
  reco_beam_true_byHits_startP = 0.;

  reco_beam_true_byE_matched = false;
  reco_beam_true_byHits_matched = false;
  reco_beam_true_byHits_purity = 0.;


  //reco_daughter_true_byE_isPrimary = false;
  reco_beam_Chi2_proton = 999.;

  reco_beam_cosmic_candidate_lower_hits.clear();
  reco_beam_cosmic_candidate_upper_hits.clear();
  reco_beam_cosmic_candidate_ID.clear();
  beam_has_cosmic_IDE = false;
  cosmic_has_beam_IDE.clear();
  n_cosmics_with_beam_IDE = -1;


  data_BI_P = 0.;
  data_BI_X = 0.;
  data_BI_Y = 0.;
  data_BI_Z = 0.;
  data_BI_dirX = 0.;
  data_BI_dirY = 0.;
  data_BI_dirZ = 0.;
  data_BI_nFibersP1 = 0;
  data_BI_nFibersP2 = 0;
  data_BI_nFibersP3 = 0;
  data_BI_PDG_candidates.clear();
  data_BI_nTracks = -1;
  data_BI_nMomenta = -1;


  quality_reco_view_0_hits_in_TPC5 = false;
  quality_reco_view_1_hits_in_TPC5 = false;
  quality_reco_view_2_hits_in_TPC5 = false;
  quality_reco_max_lateral = -999.;
  quality_reco_max_segment = -999.;

  quality_reco_view_0_max_segment = -999.;
  quality_reco_view_1_max_segment = -999.;
  quality_reco_view_2_max_segment = -999.;

  quality_reco_view_0_wire.clear(); 
  quality_reco_view_1_wire.clear(); 
  quality_reco_view_2_wire.clear(); 

  quality_reco_view_2_z.clear(); 

  quality_reco_view_0_tick.clear(); 
  quality_reco_view_1_tick.clear(); 
  quality_reco_view_2_tick.clear(); 

  quality_reco_view_0_wire_backtrack = 0.;
  quality_reco_view_1_wire_backtrack = 0.;
  quality_reco_view_2_wire_backtrack = 0.;

  reco_beam_Chi2_ndof = -1;

  /*
  reco_daughter_Chi2_proton.clear();
  reco_daughter_Chi2_ndof.clear();
  reco_daughter_momByRange_proton.clear();
  reco_daughter_momByRange_muon.clear();
  reco_daughter_allTrack_momByRange_proton.clear();
  reco_daughter_allTrack_momByRange_muon.clear();

  reco_daughter_shower_Chi2_proton.clear();
  reco_daughter_shower_Chi2_ndof.clear();

  reco_daughter_trackScore.clear();
  reco_daughter_emScore.clear();
  reco_daughter_michelScore.clear();
  
  reco_daughter_shower_trackScore.clear();
  reco_daughter_shower_emScore.clear();
  reco_daughter_shower_michelScore.clear();
  */

  true_beam_daughter_PDG.clear();
  true_beam_daughter_len.clear();
  true_beam_daughter_startX.clear();
  true_beam_daughter_startY.clear();
  true_beam_daughter_startZ.clear();
  true_beam_daughter_startPx.clear();
  true_beam_daughter_startPy.clear();
  true_beam_daughter_startPz.clear();
  true_beam_daughter_startP.clear();
  true_beam_daughter_endX.clear();
  true_beam_daughter_endY.clear();
  true_beam_daughter_endZ.clear();
  true_beam_daughter_Process.clear();
  true_beam_daughter_endProcess.clear();
  true_beam_daughter_nHits.clear();

  true_beam_daughter_reco_byHits_PFP_ID.clear();
  true_beam_daughter_reco_byHits_PFP_nHits.clear();
  true_beam_daughter_reco_byHits_PFP_trackScore.clear();

  true_beam_daughter_reco_byHits_allTrack_ID.clear();
  true_beam_daughter_reco_byHits_allTrack_startX.clear();
  true_beam_daughter_reco_byHits_allTrack_startY.clear();
  true_beam_daughter_reco_byHits_allTrack_startZ.clear();
  true_beam_daughter_reco_byHits_allTrack_endX.clear();
  true_beam_daughter_reco_byHits_allTrack_endY.clear();
  true_beam_daughter_reco_byHits_allTrack_endZ.clear();
  true_beam_daughter_reco_byHits_allTrack_len.clear();

  true_beam_daughter_reco_byHits_allShower_ID.clear();
  true_beam_daughter_reco_byHits_allShower_startX.clear();
  true_beam_daughter_reco_byHits_allShower_startY.clear();
  true_beam_daughter_reco_byHits_allShower_startZ.clear();
  true_beam_daughter_reco_byHits_allShower_len.clear();

  true_beam_Pi0_decay_ID.clear();
  true_beam_Pi0_decay_parID.clear();
  true_beam_Pi0_decay_startP.clear();
  true_beam_Pi0_decay_PDG.clear();
  true_beam_Pi0_decay_len.clear();
  true_beam_Pi0_decay_nHits.clear();
  true_beam_Pi0_decay_reco_byHits_PFP_ID.clear();
  true_beam_Pi0_decay_reco_byHits_PFP_nHits.clear();
  true_beam_Pi0_decay_reco_byHits_PFP_trackScore.clear();

  true_beam_Pi0_decay_reco_byHits_allTrack_ID.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_startX.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_startY.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_startZ.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_endX.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_endY.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_endZ.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_len.clear();

  true_beam_Pi0_decay_reco_byHits_allShower_ID.clear();
  true_beam_Pi0_decay_reco_byHits_allShower_startX.clear();
  true_beam_Pi0_decay_reco_byHits_allShower_startY.clear();
  true_beam_Pi0_decay_reco_byHits_allShower_startZ.clear();
  true_beam_Pi0_decay_reco_byHits_allShower_len.clear();

  true_beam_grand_daughter_ID.clear();
  true_beam_grand_daughter_parID.clear();
  true_beam_grand_daughter_PDG.clear();
  true_beam_grand_daughter_nHits.clear();
  true_beam_grand_daughter_Process.clear();
  true_beam_grand_daughter_endProcess.clear();
  true_beam_daughter_ID.clear();

  reco_beam_nTrackDaughters = -1;
  reco_beam_nShowerDaughters = -1;

  reco_daughter_PFP_ID.clear();
  reco_daughter_PFP_nHits.clear();
  reco_daughter_PFP_trackScore.clear();
  reco_daughter_PFP_emScore.clear();
  reco_daughter_PFP_michelScore.clear();
  reco_daughter_PFP_trackScore_collection.clear();
  reco_daughter_PFP_emScore_collection.clear();
  reco_daughter_PFP_michelScore_collection.clear();

  reco_beam_PFP_ID = -999;
  reco_beam_PFP_nHits = -999;
  reco_beam_PFP_trackScore = -999;
  reco_beam_PFP_emScore = -999;
  reco_beam_PFP_michelScore = -999;
  reco_beam_PFP_trackScore_collection = -999;
  reco_beam_PFP_emScore_collection = -999;
  reco_beam_PFP_michelScore_collection = -999;

  reco_beam_allTrack_ID = -999;
  reco_beam_allTrack_beam_cuts = -999;
  reco_beam_allTrack_flipped = -999;
  reco_beam_allTrack_len = -999;
  reco_beam_allTrack_startX = -999;
  reco_beam_allTrack_startY = -999;
  reco_beam_allTrack_startZ = -999;
  reco_beam_allTrack_endX = -999;
  reco_beam_allTrack_endY = -999;
  reco_beam_allTrack_endZ = -999;
  reco_beam_allTrack_trackDirX = -999;
  reco_beam_allTrack_trackDirY = -999;
  reco_beam_allTrack_trackDirZ = -999;
  reco_beam_allTrack_trackEndDirX = -999;
  reco_beam_allTrack_trackEndDirY = -999;
  reco_beam_allTrack_trackEndDirZ = -999;
  reco_beam_allTrack_resRange.clear();
  reco_beam_allTrack_calibrated_dEdX.clear();
  reco_beam_allTrack_Chi2_proton = -999;
  reco_beam_allTrack_Chi2_ndof = -999;




  reco_beam_dQdX.clear();
  reco_beam_dEdX.clear();
  reco_beam_calibrated_dEdX.clear();
  reco_beam_vtxX = -1.;
  reco_beam_vtxY = -1.;
  reco_beam_vtxZ = -1.;
  /*
  reco_daughter_startX.clear();
  reco_daughter_startY.clear();
  reco_daughter_startZ.clear();
  reco_daughter_endX.clear();
  reco_daughter_endY.clear();
  reco_daughter_endZ.clear();
  reco_daughter_deltaR.clear();
  reco_daughter_dR.clear();
  reco_daughter_to_vertex.clear();
  reco_daughter_slice.clear();

  reco_daughter_shower_to_vertex.clear();

  reco_daughter_shower_startX.clear();
  reco_daughter_shower_startY.clear();
  reco_daughter_shower_startZ.clear();

  reco_daughter_shower_len.clear(); 
  */

  reco_beam_resRange.clear();
  reco_beam_TrkPitch.clear();
  reco_beam_calo_wire.clear();
  reco_beam_calo_tick.clear();

/*
  reco_daughter_dQdX.clear();
  reco_daughter_dEdX.clear();
  reco_daughter_resRange.clear();
  reco_daughter_shower_dQdX.clear();
  reco_daughter_shower_dEdX.clear();
  reco_daughter_shower_resRange.clear();
  reco_daughter_len.clear();
*/
  reco_beam_trackID = -1;

  reco_beam_incidentEnergies.clear();
  reco_beam_interactingEnergy = -999.;
  true_beam_incidentEnergies.clear();
  true_beam_interactingEnergy = -999.;

  /*
  reco_daughter_trackID.clear();
  reco_daughter_true_byE_completeness.clear();
  reco_daughter_true_byE_purity.clear();
  reco_daughter_true_byE_PDG.clear();
  reco_daughter_true_byE_ID.clear();
  reco_daughter_true_byE_origin.clear();
  reco_daughter_true_byE_parID.clear();
  reco_daughter_true_byE_parPDG.clear();
  reco_daughter_true_byE_process.clear();

  reco_daughter_true_byHits_PDG.clear();
  reco_daughter_true_byHits_ID.clear();
  reco_daughter_true_byHits_origin.clear();
  reco_daughter_true_byHits_parID.clear();
  reco_daughter_true_byHits_parPDG.clear();
  reco_daughter_true_byHits_process.clear();
  reco_daughter_true_byHits_sharedHits.clear();
  reco_daughter_true_byHits_emHits.clear();

  reco_daughter_true_byHits_len.clear();
  reco_daughter_true_byHits_startX.clear();
  reco_daughter_true_byHits_startY.clear();
  reco_daughter_true_byHits_startZ.clear();
  reco_daughter_true_byHits_endX.clear();
  reco_daughter_true_byHits_endY.clear();
  reco_daughter_true_byHits_endZ.clear();

  reco_daughter_true_byHits_startPx.clear();
  reco_daughter_true_byHits_startPy.clear();
  reco_daughter_true_byHits_startPz.clear();
  reco_daughter_true_byHits_startP.clear();
  reco_daughter_true_byHits_startE.clear();
  */

  //Alternative Reco
  reco_daughter_PFP_true_byHits_PDG.clear();
  reco_daughter_PFP_true_byHits_ID.clear();
  reco_daughter_PFP_true_byHits_origin.clear();
  reco_daughter_PFP_true_byHits_parID.clear();
  reco_daughter_PFP_true_byHits_parPDG.clear();
  reco_daughter_PFP_true_byHits_process.clear();
  reco_daughter_PFP_true_byHits_sharedHits.clear();
  reco_daughter_PFP_true_byHits_emHits.clear();

  reco_daughter_PFP_true_byHits_len.clear();
  reco_daughter_PFP_true_byHits_startX.clear();
  reco_daughter_PFP_true_byHits_startY.clear();
  reco_daughter_PFP_true_byHits_startZ.clear();
  reco_daughter_PFP_true_byHits_endX.clear();
  reco_daughter_PFP_true_byHits_endY.clear();
  reco_daughter_PFP_true_byHits_endZ.clear();

  reco_daughter_PFP_true_byHits_startPx.clear();
  reco_daughter_PFP_true_byHits_startPy.clear();
  reco_daughter_PFP_true_byHits_startPz.clear();
  reco_daughter_PFP_true_byHits_startP.clear();
  reco_daughter_PFP_true_byHits_startE.clear();
  reco_daughter_PFP_true_byHits_endProcess.clear();
  reco_daughter_PFP_true_byHits_purity.clear();

  reco_daughter_allTrack_ID.clear();
  reco_daughter_allTrack_dEdX.clear();
  reco_daughter_allTrack_dQdX.clear();
  reco_daughter_allTrack_resRange.clear();
  reco_daughter_allTrack_dEdX_SCE.clear();
  reco_daughter_allTrack_dQdX_SCE.clear();
  reco_daughter_allTrack_resRange_SCE.clear();

  reco_daughter_allTrack_calibrated_dEdX.clear();
  reco_daughter_allTrack_calibrated_dEdX_SCE.clear();

  reco_daughter_allTrack_Chi2_proton.clear();
  reco_daughter_allTrack_Chi2_ndof.clear();

  reco_daughter_allTrack_Theta.clear();
  reco_daughter_allTrack_Phi.clear();
  reco_daughter_allTrack_len.clear();
  reco_daughter_allTrack_startX.clear();
  reco_daughter_allTrack_startY.clear();
  reco_daughter_allTrack_startZ.clear();
  reco_daughter_allTrack_endX.clear();
  reco_daughter_allTrack_endY.clear();
  reco_daughter_allTrack_endZ.clear();
  reco_daughter_allTrack_dR.clear();
  reco_daughter_allTrack_to_vertex.clear();

  reco_daughter_allShower_ID.clear();
  reco_daughter_allShower_len.clear();
  reco_daughter_allShower_startX.clear();
  reco_daughter_allShower_startY.clear();
  reco_daughter_allShower_startZ.clear();

  ///////
  

  /*
  reco_daughter_shower_true_byHits_PDG.clear();
  reco_daughter_shower_true_byHits_ID.clear();
  reco_daughter_shower_true_byHits_origin.clear();
  reco_daughter_shower_true_byHits_parID.clear();
  reco_daughter_shower_true_byHits_parPDG.clear();
  reco_daughter_shower_true_byHits_process.clear();
  reco_daughter_shower_true_byHits_purity.clear();
  reco_daughter_true_byHits_purity.clear();
  reco_daughter_shower_true_byHits_startPx.clear();
  reco_daughter_shower_true_byHits_startPy.clear();
  reco_daughter_shower_true_byHits_startPz.clear();
  reco_daughter_shower_true_byHits_startP.clear();
  reco_daughter_shower_true_byHits_endProcess.clear();





  reco_daughter_showerID.clear();
  reco_daughter_shower_true_byE_PDG.clear();
  reco_daughter_shower_true_byE_ID.clear();
  reco_daughter_shower_true_byE_origin.clear();
  reco_daughter_shower_true_byE_parID.clear();
  reco_daughter_shower_true_byE_parPDG.clear();

  reco_daughter_shower_true_byE_startPx.clear();
  reco_daughter_shower_true_byE_startPy.clear();
  reco_daughter_shower_true_byE_startPz.clear();
  reco_daughter_shower_true_byE_startP.clear();
  reco_daughter_shower_true_byE_endProcess.clear();
  */




  //New Hits info
  reco_beam_spacePts_X.clear();
  reco_beam_spacePts_Y.clear();
  reco_beam_spacePts_Z.clear();

  reco_daughter_spacePts_X.clear();
  reco_daughter_spacePts_Y.clear();
  reco_daughter_spacePts_Z.clear();

  reco_daughter_shower_spacePts_X.clear();
  reco_daughter_shower_spacePts_Y.clear();
  reco_daughter_shower_spacePts_Z.clear();
  //

}

DEFINE_ART_MODULE(pionana::PionAnalyzer)
