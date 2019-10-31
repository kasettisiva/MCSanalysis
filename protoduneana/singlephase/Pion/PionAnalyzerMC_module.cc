////////////////////////////////////////////////////////////////////////
// Class:       PionAnalyzerMC
// Plugin Type: analyzer (art v3_00_00)
// File:        PionAnalyzerMC_module.cc
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

#include "dune/Protodune/singlephase/DataUtils/ProtoDUNETrackUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEShowerUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNETruthUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEBeamlineUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEBeamCuts.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

#include "lardata/ArtDataHelper/MVAReader.h"

#include "art_root_io/TFileService.h"
#include "TProfile.h"
#include "TFile.h"

// ROOT includes
#include "TTree.h"

namespace pionana {
  class PionAnalyzerMC;
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
}



class pionana::PionAnalyzerMC : public art::EDAnalyzer {
public:
  explicit PionAnalyzerMC(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PionAnalyzerMC(PionAnalyzerMC const&) = delete;
  PionAnalyzerMC(PionAnalyzerMC&&) = delete;
  PionAnalyzerMC& operator=(PionAnalyzerMC const&) = delete;
  PionAnalyzerMC& operator=(PionAnalyzerMC&&) = delete;

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




  /******************************/
  //Truth level info of the primary beam particle
  //that generated the event
  int true_beam_PDG;
  int true_beam_ID;
  std::string true_beam_EndProcess;
  double true_beam_EndVertex_X;
  double true_beam_EndVertex_Y;
  double true_beam_EndVertex_Z;
  double true_beam_Start_X;
  double true_beam_Start_Y;
  double true_beam_Start_Z;

  double true_beam_Start_DirX;
  double true_beam_Start_DirY;
  double true_beam_Start_DirZ;

  double true_beam_Start_Px;
  double true_beam_Start_Py;
  double true_beam_Start_Pz;
  double true_beam_Start_P;

  double true_beam_End_Px;
  double true_beam_End_Py;
  double true_beam_End_Pz;
  double true_beam_End_P;

    //Truth level info of the daughter MCParticles coming out of the 
    //true primary particle
  std::vector< int > true_beam_daughter_PDGs;
  std::vector< int > true_beam_daughter_IDs;
  std::vector< double > true_beam_daughter_lens;
  std::vector< std::string > true_beam_daughter_Process;

  std::vector< double > true_beam_daughter_startX, true_beam_daughter_startY, true_beam_daughter_startZ;
  std::vector< double > true_beam_daughter_startP, true_beam_daughter_startPX, true_beam_daughter_startPY, true_beam_daughter_startPZ;


  //Decay products from pi0s
  std::vector< int > true_beam_Pi0_decay_PDGs, true_beam_Pi0_decay_IDs;
  std::vector< int > true_beam_grand_daughter_PDGs, true_beam_grand_daughter_IDs, true_beam_grand_daughter_ParIDs;

  //How many of each true particle came out of the true primary beam particle?
  int nPiPlus_truth, nPiMinus_truth, nPi0_truth;
  int nProton_truth, nNeutron_truth, nNucleus_truth;
  /*****************************/

  //Matched to vertex/slice?
  int vertex_type = kUnmatched;
  int vertex_slice;
  int endZ_slice, startZ_slice;
  bool vertex_matched_to_beam, endZ_matched_to_beam, startZ_matched_to_beam;

  std::vector< int > vertex_IDs, endZ_IDs, startZ_IDs;
  std::vector< std::vector< double > > vertex_dRs, endZ_dRs, startZ_dRs;
  std::vector< double > centroid_dR;



  //Reconstructed track info
  double startX, startY, startZ;
  double endX, endY, endZ;
  double vtxX, vtxY, vtxZ; 
  double len;
  double trackDirX, trackDirY, trackDirZ;
  double trackEndDirX, trackEndDirY, trackEndDirZ;
  std::vector< double > dEdX, dQdX, resRange;
  int beamTrackID;

  std::string reco_beam_truth_EndProcess, alt_reco_beam_truth_EndProcess; //What process ended the reco beam particle
  std::string reco_beam_truth_Process, alt_reco_beam_truth_Process;    //What process created the reco beam particle
  
  int reco_beam_truth_PDG, alt_reco_beam_truth_PDG; //What is the PDG of the true MC particle contributing the most to 
                           //the reconstructed beam track
  int reco_beam_truth_ID, alt_reco_beam_truth_ID;
  
  bool reco_beam_good, alt_reco_beam_good; //Does the true particle contributing most to the 
                       //reconstructed beam track coincide with the actual
                       //beam particle that generated the event
  bool passes_beam_cuts;                       

  double alt_reco_beam_purity;                      
  
  bool flipped;
                       
  bool elastic_candidate, found_elastic_scatter;
  int  nElasticScatters;
  std::vector< double > elastic_costheta, elastic_X, elastic_Y, elastic_Z;
  int nReconstructed;
  std::vector< int > reconstructedIDs;
  bool daughter_is_primary;


  std::vector< std::string > processes;

  int reco_beam_truth_origin, alt_reco_beam_truth_origin; //What is the origin of the reconstructed beam track?

  double traj_cos;
  bool view_0_hits_in_TPC5, view_1_hits_in_TPC5, view_2_hits_in_TPC5;

  std::vector< double > view_0_wire, view_0_tick;
  std::vector< double > view_1_wire, view_1_tick;
  std::vector< double > view_2_wire, view_2_tick;
  std::vector< double > view_2_z;
  double view_0_max_segment, view_1_max_segment, view_2_max_segment;
  double view_0_wire_backtrack, view_1_wire_backtrack, view_2_wire_backtrack;

  double max_lateral, max_segment; 

  double reco_beam_Chi2_proton;
  int    reco_beam_Chi2_ndof;

  std::vector< double > reco_daughter_Chi2_proton;
  std::vector< int >    reco_daughter_Chi2_ndof;

  std::vector< double > reco_daughter_shower_Chi2_proton;
  std::vector< int >    reco_daughter_shower_Chi2_ndof;

  //Truth-level info of the reconstructed particles coming out of the 
  //reconstructed beam track
  std::vector< bool >   reco_beam_truth_daughter_good_reco;
  std::vector< int >    reco_beam_truth_daughter_true_PDGs;
  std::vector< int >    reco_beam_truth_daughter_true_IDs;
  std::vector< double > reco_beam_truth_daughter_true_lens;

  std::vector< double > reco_beam_truth_daughter_shower_true_lens;
  std::vector< bool >   reco_beam_truth_daughter_shower_good_reco;
  std::vector< int >    reco_beam_truth_daughter_shower_true_PDGs;
  std::vector< int >    reco_beam_truth_daughter_shower_true_IDs;

  double reco_beam_truth_End_Px,   alt_reco_beam_truth_End_Px;
  double reco_beam_truth_End_Py,   alt_reco_beam_truth_End_Py;
  double reco_beam_truth_End_Pz,   alt_reco_beam_truth_End_Pz;
  double reco_beam_truth_End_E,    alt_reco_beam_truth_End_E;
  double reco_beam_truth_End_P,    alt_reco_beam_truth_End_P;
                                   
  double reco_beam_truth_Start_Px, alt_reco_beam_truth_Start_Px;
  double reco_beam_truth_Start_Py, alt_reco_beam_truth_Start_Py;
  double reco_beam_truth_Start_Pz, alt_reco_beam_truth_Start_Pz;
  double reco_beam_truth_Start_E,  alt_reco_beam_truth_Start_E;
  double reco_beam_truth_Start_P,  alt_reco_beam_truth_Start_P;




  //Reco-level info of the reconstructed daughters coming out of the
  //reconstructed beam tracl
  std::vector< int > reco_daughter_trackID;
  std::vector< double > reco_daughter_completeness;
  std::vector< double > reco_daughter_purity;
  std::vector< int > reco_daughter_truth_PDG;
  std::vector< int > reco_daughter_truth_ID;
  std::vector< int > reco_daughter_truth_Origin;
  std::vector< int > reco_daughter_truth_ParID;
  std::vector< std::string > reco_daughter_truth_Process;
  std::vector< double > reco_daughter_truth_purity;

  std::vector< int > alt_reco_daughter_truth_PDG;
  std::vector< int > alt_reco_daughter_truth_ID;
  std::vector< int > alt_reco_daughter_truth_Origin;
  std::vector< int > alt_reco_daughter_truth_ParID;
  std::vector< std::string > alt_reco_daughter_truth_Process;
  std::vector< double > alt_reco_daughter_truth_purity;
  std::vector< size_t > alt_reco_daughter_truth_shared_hits, alt_reco_daughter_truth_em_hits;

  std::vector< double > alt_reco_daughter_truth_len;
  std::vector< double > alt_reco_daughter_truth_startX;
  std::vector< double > alt_reco_daughter_truth_startY;
  std::vector< double > alt_reco_daughter_truth_startZ;
  std::vector< double > alt_reco_daughter_truth_endX;
  std::vector< double > alt_reco_daughter_truth_endY;
  std::vector< double > alt_reco_daughter_truth_endZ;

  //Alternative Reco values
  std::vector< int > reco_daughter_PFP_ID;
  std::vector< double > reco_daughter_PFP_track_score;
  std::vector< double > reco_daughter_PFP_em_score;
  std::vector< double > reco_daughter_PFP_michel_score;
  std::vector< double > reco_daughter_PFP_none_score;


  std::vector< int > alt_reco_daughter_PFP_truth_PDG;
  std::vector< int > alt_reco_daughter_PFP_truth_ID;
  std::vector< int > alt_reco_daughter_PFP_truth_Origin;
  std::vector< int > alt_reco_daughter_PFP_truth_ParID;
  std::vector< std::string > alt_reco_daughter_PFP_truth_Process;
  std::vector< double > alt_reco_daughter_PFP_truth_purity;
  std::vector< size_t > alt_reco_daughter_PFP_truth_shared_hits, alt_reco_daughter_PFP_truth_em_hits;

  std::vector< double > alt_reco_daughter_PFP_truth_len;
  std::vector< double > alt_reco_daughter_PFP_truth_startX;
  std::vector< double > alt_reco_daughter_PFP_truth_startY;
  std::vector< double > alt_reco_daughter_PFP_truth_startZ;
  std::vector< double > alt_reco_daughter_PFP_truth_endX;
  std::vector< double > alt_reco_daughter_PFP_truth_endY;
  std::vector< double > alt_reco_daughter_PFP_truth_endZ;

  std::vector< int > reco_daughter_PMA_ID;
  std::vector< std::vector< double > > reco_daughter_PMA_dEdX, reco_daughter_PMA_resRange;
  std::vector< double > reco_daughter_PMA_Chi2_proton;
  std::vector< int >    reco_daughter_PMA_Chi2_ndof;

  std::vector< int > reco_daughter_pandora2_ID;
  std::vector< std::vector< double > > reco_daughter_pandora2_dEdX, reco_daughter_pandora2_resRange;
  std::vector< double > reco_daughter_pandora2_Chi2_proton;
  std::vector< int >    reco_daughter_pandora2_Chi2_ndof;

  std::vector< double > reco_daughter_PMA_startX, reco_daughter_PMA_endX;
  std::vector< double > reco_daughter_PMA_startY, reco_daughter_PMA_endY;
  std::vector< double > reco_daughter_PMA_startZ, reco_daughter_PMA_endZ;
  std::vector< double > reco_daughter_PMA_dR;
  std::vector< double > reco_daughter_PMA_to_vertex;

  std::vector< double > reco_daughter_pandora2_startX, reco_daughter_pandora2_endX;
  std::vector< double > reco_daughter_pandora2_startY, reco_daughter_pandora2_endY;
  std::vector< double > reco_daughter_pandora2_startZ, reco_daughter_pandora2_endZ;
  std::vector< double > reco_daughter_pandora2_dR;
  std::vector< double > reco_daughter_pandora2_to_vertex;
  //


  std::vector< int > alt_reco_daughter_shower_truth_PDG;
  std::vector< int > alt_reco_daughter_shower_truth_ID;
  std::vector< int > alt_reco_daughter_shower_truth_Origin;
  std::vector< int > alt_reco_daughter_shower_truth_ParID;
  std::vector< std::string > alt_reco_daughter_shower_truth_Process;
  std::vector< double > alt_reco_daughter_shower_truth_purity;

  std::vector< int > reco_daughter_showerID;
  std::vector< int > reco_daughter_shower_truth_PDG;
  std::vector< int > reco_daughter_shower_truth_ID;
  std::vector< int > reco_daughter_shower_truth_Origin;
  std::vector< int > reco_daughter_shower_truth_ParID;
  std::vector< std::vector< double > > reco_daughter_dEdX, reco_daughter_dQdX, reco_daughter_resRange;
  std::vector< std::vector< double > > reco_daughter_shower_dEdX, reco_daughter_shower_dQdX, reco_daughter_shower_resRange;
  std::vector< double > reco_daughter_startX, reco_daughter_endX;
  std::vector< double > reco_daughter_startY, reco_daughter_endY;
  std::vector< double > reco_daughter_startZ, reco_daughter_endZ;
  std::vector< double > reco_daughter_deltaR;
  std::vector< double > reco_daughter_dR;
  std::vector< double > reco_daughter_to_vertex;
  std::vector< double > reco_daughter_shower_to_vertex;
  std::vector< int > reco_daughter_slice;
  std::vector< double > reco_daughter_shower_startX;
  std::vector< double > reco_daughter_shower_startY;
  std::vector< double > reco_daughter_shower_startZ;
  std::vector< double > reco_daughter_shower_len;
  std::vector< double > reco_daughter_len;
  std::vector< double > reco_daughter_track_score;
  std::vector< double > reco_daughter_em_score;
  std::vector< double > reco_daughter_none_score;
  std::vector< double > reco_daughter_michel_score;

  std::vector< double > reco_daughter_shower_track_score;
  std::vector< double > reco_daughter_shower_em_score;
  std::vector< double > reco_daughter_shower_none_score;
  std::vector< double > reco_daughter_shower_michel_score;


  int nTrackDaughters, nShowerDaughters;

  int type;
  int nBeamParticles;

  std::map< int, TProfile* > templates;

  //FCL pars
  std::string fCalorimetryTag;
  std::string fTrackerTag;    
  std::string fHitTag;    
  std::string fShowerTag;     
  std::string fPFParticleTag; 
  std::string fGeneratorTag;
  std::string dEdX_template_name;
  TFile dEdX_template_file;
  bool fVerbose;             
  int fNSliceCheck; 
  fhicl::ParameterSet BeamPars;
  fhicl::ParameterSet BeamCuts;
  protoana::ProtoDUNEBeamCuts beam_cuts;
};


pionana::PionAnalyzerMC::PionAnalyzerMC(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  ,
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),

  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fHitTag(p.get<std::string>("HitTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  dEdX_template_name(p.get<std::string>("dEdX_template_name")),
  dEdX_template_file( dEdX_template_name.c_str(), "OPEN" ),
  fVerbose(p.get<bool>("Verbose")),
  fNSliceCheck( p.get< int >("NSliceCheck") ),
  BeamPars(p.get<fhicl::ParameterSet>("BeamPars")),
  BeamCuts(p.get<fhicl::ParameterSet>("BeamCuts"))
{

  templates[ 211 ]  = (TProfile*)dEdX_template_file.Get( "dedx_range_pi"  );
  templates[ 321 ]  = (TProfile*)dEdX_template_file.Get( "dedx_range_ka"  );
  templates[ 13 ]   = (TProfile*)dEdX_template_file.Get( "dedx_range_mu"  );
  templates[ 2212 ] = (TProfile*)dEdX_template_file.Get( "dedx_range_pro" );

  beam_cuts = protoana::ProtoDUNEBeamCuts( BeamCuts );

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pionana::PionAnalyzerMC::analyze(art::Event const& evt)
{

  //reset containers
  reset();  


  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();

   
  // Get various utilities 
  protoana::ProtoDUNEPFParticleUtils                    pfpUtil;
  protoana::ProtoDUNETrackUtils                         trackUtil;
  protoana::ProtoDUNEShowerUtils                        showerUtil;
  art::ServiceHandle<cheat::BackTrackerService>         bt_serv;
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  protoana::ProtoDUNETruthUtils                         truthUtil;
  ////////////////////////////////////////
  


  // This gets the true beam particle that generated the event
  auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
  const simb::MCParticle* true_beam_particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
  if( !true_beam_particle ){
    std::cout << "No true beam particle" << std::endl;
    return;
  }
  ////////////////////////////
  
  

  // Helper to get hits and the 4 associated CNN outputs
  // CNN Outputs: EM, Track, Michel, Empty
  // outputNames: track, em, none, michel
  anab::MVAReader<recob::Hit,4> hitResults(evt, /*fNNetModuleLabel*/ "emtrkmichelid:emtrkmichel" );


  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fTrackerTag);
  art::FindManyP<recob::Hit> findHits(recoTracks,evt,fTrackerTag);

  auto recoShowers = evt.getValidHandle< std::vector< recob::Shower > >(fShowerTag);
  art::FindManyP<recob::Hit> findHitsFromShowers(recoShowers,evt,fShowerTag);

  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  nBeamParticles = beamParticles.size();

  if(beamParticles.size() == 0){
    std::cerr << "We found no beam particles for this event... moving on" << std::endl;
    return;
  }
  //////////////////////////////////////////////////////////////////


  // Get the reconstructed PFParticle tagged as beam by Pandora
  const recob::PFParticle* particle = beamParticles.at(0);


  //Reco track vertex and direction/length
  const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);

  // Determine if the beam particle is track-like or shower-like
  const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
  const simb::MCParticle* trueParticle = 0x0;

  if(thisShower != 0x0){
    std::cout << "Beam particle is shower-like" << std::endl;
    type = 11;
    //This gets the MCParticle contributing the most to the shower
    trueParticle = truthUtil.GetMCParticleFromRecoShower(*thisShower, evt, fShowerTag);
  }
  if(thisTrack != 0x0){

    //This gets the MCParticle contributing the most to the track
    trueParticle = truthUtil.GetMCParticleFromRecoTrack(*thisTrack, evt, fTrackerTag);

    
    
    std::cout << "Trying " << pi_serv->TrackIdToParticle_P(trueParticle->TrackId())->TrackId() << " " 
              << pi_serv->TrackIdToParticle_P(-1*trueParticle->TrackId())->TrackId() << std::endl;


    std::cout << "Beam particle is track-like " << thisTrack->ID() << std::endl;
    type = 13;

    // Now we can look for the interaction point of the particle if one exists, i.e where the particle
    // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
    // check this by making sure we get a sensible position
    const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    vtxX = interactionVtx.X();
    vtxY = interactionVtx.Y();
    vtxZ = interactionVtx.Z();
    std::cout << "Secondary Vertex: " << vtxX << " " <<  vtxY << " " <<  vtxZ << std::endl;
    ////////////////////////////////////////////
    
    auto list = truthUtil.GetMCParticleListByHits( *thisTrack, evt, fTrackerTag, fHitTag );
    for( size_t j = 0; j < list.size(); ++j ){
      //std::cout << "Contrib " << j << " " << list[j].first->TrackId() << " " << list[j].second << std::endl;
      std::cout << "Contrib " << j << " " << list[j].particle->TrackId() << " " << list[j].nSharedHits << " " << list[j].nSharedDeltaRayHits << std::endl;
    }

    protoana::MCParticleSharedHits beam_match = truthUtil.GetMCParticleByHits( *thisTrack, evt, fTrackerTag, fHitTag );
    if( beam_match.particle ){
      //Check that this is the correct true particle
      if( beam_match.particle->TrackId() == true_beam_particle->TrackId() ){
        alt_reco_beam_good = true;
      }

      alt_reco_beam_truth_PDG = beam_match.particle->PdgCode();
      alt_reco_beam_truth_ID = beam_match.particle->TrackId();
      std::cout << "Truth ID: " << alt_reco_beam_truth_ID << std::endl;

      alt_reco_beam_truth_Process = beam_match.particle->Process();
      alt_reco_beam_truth_EndProcess = beam_match.particle->EndProcess();
      alt_reco_beam_truth_origin = pi_serv->TrackIdToMCTruth_P(beam_match.particle->TrackId())->Origin();

      alt_reco_beam_truth_Start_Px = beam_match.particle->Px();
      alt_reco_beam_truth_Start_Py = beam_match.particle->Py();
      alt_reco_beam_truth_Start_Pz = beam_match.particle->Pz();
      alt_reco_beam_truth_Start_P  = sqrt( alt_reco_beam_truth_Start_Px*alt_reco_beam_truth_Start_Px 
                                     + alt_reco_beam_truth_Start_Py*alt_reco_beam_truth_Start_Py 
                                     + alt_reco_beam_truth_Start_Pz*alt_reco_beam_truth_Start_Pz );
      alt_reco_beam_truth_Start_E = beam_match.particle->E();

      size_t np = beam_match.particle->NumberTrajectoryPoints();
      if( np > 1 ){
        alt_reco_beam_truth_End_Px = beam_match.particle->Px( np - 2 );
        alt_reco_beam_truth_End_Py = beam_match.particle->Py( np - 2 );
        alt_reco_beam_truth_End_Pz = beam_match.particle->Pz( np - 2 );
        alt_reco_beam_truth_End_P  = sqrt( alt_reco_beam_truth_End_Px*alt_reco_beam_truth_End_Px 
                                     + alt_reco_beam_truth_End_Py*alt_reco_beam_truth_End_Py 
                                     + alt_reco_beam_truth_End_Pz*alt_reco_beam_truth_End_Pz );
        alt_reco_beam_truth_End_E  = beam_match.particle->E( np - 2 );
      }

      auto list = truthUtil.GetMCParticleListByHits( *thisTrack, evt, fTrackerTag, fHitTag );
      double total = 0.;
      double matched_hits = 0.;
      for( size_t j = 0; j < list.size(); ++j ){
      //  std::cout << "Contrib " << j << " " << list[j].first->TrackId() << " " << list[j].second << std::endl;
        std::cout << "Contrib " << j << " " << list[j].particle->TrackId() << " " << list[j].particle->PdgCode() 
                   << " " << pi_serv->TrackIdToMCTruth_P(list[j].particle->TrackId())->Origin() 
                   << " " << list[j].nSharedHits << " " << list[j].nSharedDeltaRayHits << std::endl;

        if( list[j].particle == beam_match.particle ){
           matched_hits = list[j].nSharedHits + list[j].nSharedDeltaRayHits;
        }

        total += list[j].nSharedHits + list[j].nSharedDeltaRayHits;
      }

      alt_reco_beam_purity = ( matched_hits / total );


    }
    
  }
  //////////////////////////////////////////////////

  if( trueParticle ){

    //Check that this is the correct true particle
    if( trueParticle->TrackId() == true_beam_particle->TrackId() ){
      reco_beam_good = true;
    }

    reco_beam_truth_PDG = trueParticle->PdgCode();
    reco_beam_truth_ID = trueParticle->TrackId();
    std::cout << "Truth ID: " << reco_beam_truth_ID << std::endl;

    reco_beam_truth_Process = trueParticle->Process();
    reco_beam_truth_EndProcess = trueParticle->EndProcess();
    reco_beam_truth_origin = pi_serv->TrackIdToMCTruth_P(trueParticle->TrackId())->Origin();
    //What created the MCPraticle that created the track?
    /*
    std::cout << "Track-to-True origin: ";
    switch( reco_beam_truth_origin ){
      case simb::kCosmicRay: 
        std::cout << "Cosmic" << std::endl;
        break;
      case simb::kSingleParticle:
        std::cout << "Beam" << std::endl;
        break;
      default:
        std::cout << "Other" << std::endl;
    }
    */

    reco_beam_truth_Start_Px = trueParticle->Px();
    reco_beam_truth_Start_Py = trueParticle->Py();
    reco_beam_truth_Start_Pz = trueParticle->Pz();
    reco_beam_truth_Start_P  = sqrt( reco_beam_truth_Start_Px*reco_beam_truth_Start_Px 
                                   + reco_beam_truth_Start_Py*reco_beam_truth_Start_Py 
                                   + reco_beam_truth_Start_Pz*reco_beam_truth_Start_Pz );
    reco_beam_truth_Start_E = trueParticle->E();

    size_t np = trueParticle->NumberTrajectoryPoints();
    if( np > 1 ){
      reco_beam_truth_End_Px = trueParticle->Px( np - 2 );
      reco_beam_truth_End_Py = trueParticle->Py( np - 2 );
      reco_beam_truth_End_Pz = trueParticle->Pz( np - 2 );
      reco_beam_truth_End_P  = sqrt( reco_beam_truth_End_Px*reco_beam_truth_End_Px 
                                   + reco_beam_truth_End_Py*reco_beam_truth_End_Py 
                                   + reco_beam_truth_End_Pz*reco_beam_truth_End_Pz );
      reco_beam_truth_End_E  = trueParticle->E( np - 2 );
    }

  }
  //////////////////////////////////////////////////////////////
  //





  //Some truth information
  true_beam_EndProcess = true_beam_particle->EndProcess();
  
  true_beam_PDG         = true_beam_particle->PdgCode();
  true_beam_ID          = true_beam_particle->TrackId();
  true_beam_EndVertex_X = true_beam_particle->EndX();
  true_beam_EndVertex_Y = true_beam_particle->EndY();
  true_beam_EndVertex_Z = true_beam_particle->EndZ();
  true_beam_Start_X     = true_beam_particle->Position(0).X();
  true_beam_Start_Y     = true_beam_particle->Position(0).Y();
  true_beam_Start_Z     = true_beam_particle->Position(0).Z();

  true_beam_Start_Px    = true_beam_particle->Px();
  true_beam_Start_Py    = true_beam_particle->Py();
  true_beam_Start_Pz    = true_beam_particle->Pz();
  true_beam_Start_P     = true_beam_particle->P();

  size_t true_np = true_beam_particle->NumberTrajectoryPoints();

  true_beam_End_Px    = true_beam_particle->Px(true_np-2);
  true_beam_End_Py    = true_beam_particle->Py(true_np-2);
  true_beam_End_Pz    = true_beam_particle->Pz(true_np-2);
  true_beam_End_P     = true_beam_particle->P(true_np-2);

  true_beam_Start_DirX  = true_beam_Start_Px / true_beam_Start_P;
  true_beam_Start_DirY  = true_beam_Start_Py / true_beam_Start_P;
  true_beam_Start_DirZ  = true_beam_Start_Pz / true_beam_Start_P;
  

  //Go through and look at the reconstructed track->true IDs. Look for a match to the true
  //incident beam particle
  for( auto const & tr : *recoTracks ){
   const simb::MCParticle* true_p = truthUtil.GetMCParticleFromRecoTrack(tr, evt, fTrackerTag); 

   if(true_p){
     if( true_p->TrackId() == true_beam_particle->TrackId() ){
       ++nReconstructed;
       reconstructedIDs.push_back( tr.ID() );
     }
   }
  }



  //Look at true daughters coming out of the true beam particle
  //std::cout << "Has " << true_beam_particle->NumberDaughters() << " daughters" << std::endl;

  const sim::ParticleList & plist = pi_serv->ParticleList(); 
  for( int i = 0; i < true_beam_particle->NumberDaughters(); ++i ){
    int daughterID = true_beam_particle->Daughter(i);

    std::cout << "Daughter " << i << " ID: " << daughterID << std::endl;
    auto part = plist[ daughterID ];
    int pid = part->PdgCode();
    true_beam_daughter_PDGs.push_back(pid);
    true_beam_daughter_IDs.push_back( part->TrackId() );      
    true_beam_daughter_lens.push_back( part->Trajectory().TotalLength() );

    true_beam_daughter_startX.push_back( part->Position(0).X() );
    true_beam_daughter_startY.push_back( part->Position(0).Y() );
    true_beam_daughter_startZ.push_back( part->Position(0).Z() );

    true_beam_daughter_startPX.push_back( part->Px() );
    true_beam_daughter_startPY.push_back( part->Py() );
    true_beam_daughter_startPZ.push_back( part->Pz() );
    true_beam_daughter_startP.push_back( part->P() );

    true_beam_daughter_Process.push_back( part->Process() );

    std::cout << "Proccess: " << part->Process() << std::endl; 
    std::cout << "PID: " << pid << std::endl;
    std::cout << "Start: " << part->Position(0).X() << " " << part->Position(0).Y() << " " << part->Position(0).Z() << std::endl;
    std::cout << "End: " << part->EndPosition().X() << " " << part->EndPosition().Y() << " " << part->EndPosition().Z() << std::endl;
    std::cout << "Len: " << part->Trajectory().TotalLength() << std::endl;

    if( part->Process().find( "Inelastic" ) != std::string::npos ){
      std::cout << "Inelastic" << std::endl;
      if( pid == 211  ) ++nPiPlus_truth;
      if( pid == -211 ) ++nPiMinus_truth;
      if( pid == 111  ) ++nPi0_truth;
      if( pid == 2212 ) ++nProton_truth;
      if( pid == 2112 ) ++nNeutron_truth;
      if( pid > 2212  ) ++nNucleus_truth; 
    }

    //Look for the gammas coming out of the pi0s
    if( pid == 111 ){
      //std::cout << "Found pi0. Looking at true daughters" << std::endl;
      for( int j = 0; j < part->NumberDaughters(); ++j ){
        int pi0_decay_daughter_ID = part->Daughter(j);
        auto pi0_decay_part = plist[ pi0_decay_daughter_ID ];
        true_beam_Pi0_decay_PDGs.push_back( pi0_decay_part->PdgCode() );
        true_beam_Pi0_decay_IDs.push_back( pi0_decay_part->TrackId() );
      }
    }

    for( int j = 0; j < part->NumberDaughters(); ++j ){
      int grand_daughter_ID = part->Daughter(j);
      auto grand_daughter_part = plist[ grand_daughter_ID ];
      true_beam_grand_daughter_PDGs.push_back( grand_daughter_part->PdgCode() );
      true_beam_grand_daughter_IDs.push_back(  grand_daughter_part->TrackId() );
      true_beam_grand_daughter_ParIDs.push_back(  part->TrackId() );
    }

  }

  if( thisTrack ){

    passes_beam_cuts = beam_cuts.IsBeamlike( *thisTrack, evt, "1" );
    std::cout << "Beam Cuts " << passes_beam_cuts << std::endl;


    beamTrackID = thisTrack->ID();

    startX = thisTrack->Trajectory().Start().X();
    startY = thisTrack->Trajectory().Start().Y();
    startZ = thisTrack->Trajectory().Start().Z();
    endX = thisTrack->Trajectory().End().X();
    endY = thisTrack->Trajectory().End().Y();
    endZ = thisTrack->Trajectory().End().Z();

    auto startDir = thisTrack->StartDirection();
    auto endDir   = thisTrack->EndDirection();

    //try flipping
    if( startZ > endZ ){
      flipped = true;
      std::cout << "startZ > endZ: " << startZ << " " << endZ << std::endl;
      endX = thisTrack->Trajectory().Start().X();
      endY = thisTrack->Trajectory().Start().Y();
      endZ = thisTrack->Trajectory().Start().Z();
      startX = thisTrack->Trajectory().End().X();
      startY = thisTrack->Trajectory().End().Y();
      startZ = thisTrack->Trajectory().End().Z();
      
      trackDirX =  -1. * endDir.X(); 
      trackDirY =  -1. * endDir.Y(); 
      trackDirZ =  -1. * endDir.Z(); 

      trackEndDirX =  -1. * startDir.X(); 
      trackEndDirY =  -1. * startDir.Y(); 
      trackEndDirZ =  -1. * startDir.Z(); 
    }
    else{
      flipped = false;
      std::cout << "endZ > startZ: " << startZ << " " << endZ << std::endl;
      trackDirX    =  startDir.X(); 
      trackDirY    =  startDir.Y(); 
      trackDirZ    =  startDir.Z(); 
      trackEndDirX =  endDir.X(); 
      trackEndDirY =  endDir.Y(); 
      trackEndDirZ =  endDir.Z(); 
    }
    std::cout << "StartDir: " << startDir.Z() << std::endl;
    std::cout << "EndDir: "   << endDir.Z() << std::endl;

    std::cout << "trackDirX: " << trackDirX << std::endl;
    std::cout << "trackDirY: " << trackDirY << std::endl;
    std::cout << "trackDirZ: " << trackDirZ << std::endl;

    len  = thisTrack->Length();    
    std::cout << "Start: " << startX << " " << startY << " " << startZ << std::endl;
    std::cout << "End: " << endX << " " << endY << " " << endZ << std::endl;
    std::cout << "len: " << len << std::endl;
    ////////////////////////////////////////////////////////////////


    //Go through the true processes within the MCTrajectory
    const simb::MCTrajectory & true_beam_trajectory = true_beam_particle->Trajectory();
    auto true_beam_proc_map = true_beam_trajectory.TrajectoryProcesses();
    std::cout << "Processes: " << std::endl;

    std::vector< size_t > elastic_indices;

    for( auto itProc = true_beam_proc_map.begin(); itProc != true_beam_proc_map.end(); ++itProc ){
      int index = itProc->first;
      std::string process = true_beam_trajectory.KeyToProcess(itProc->second);
      std::cout << index << " " << process << std::endl;

      processes.push_back( process );

      if( process == "hadElastic" ){
        elastic_indices.push_back( index );

        ++nElasticScatters;

        double process_X = true_beam_trajectory.X( index );
        double process_Y = true_beam_trajectory.Y( index );
        double process_Z = true_beam_trajectory.Z( index );

        std::cout << "At "  
                  << process_X << " " 
                  << process_Y << " " 
                  << process_Z << std::endl;

        double delta = sqrt( 
          (process_X - endX)*(process_X - endX) +
          (process_Y - endY)*(process_Y - endY) +
          (process_Z - endZ)*(process_Z - endZ)
        );

        double PX      = true_beam_trajectory.Px( index );
        double next_PX = true_beam_trajectory.Px( index + 1 );
        double PY      = true_beam_trajectory.Py( index );
        double next_PY = true_beam_trajectory.Py( index + 1 );
        double PZ      = true_beam_trajectory.Pz( index );
        double next_PZ = true_beam_trajectory.Pz( index + 1 );

        double total_P = sqrt( PX*PX + PY*PY + PZ*PZ );
        double total_next_P = sqrt( next_PX*next_PX + next_PY*next_PY + next_PZ*next_PZ );

        //Get the angle between the direction of this step and the next
        elastic_costheta.push_back(
          ( ( PX * next_PX ) + ( PY * next_PY ) + ( PZ * next_PZ ) ) / ( total_P * total_next_P )
        );

        elastic_X.push_back( process_X );
        elastic_Y.push_back( process_Y );
        elastic_Z.push_back( process_Z );

        //for now: 10cm
        if( delta < 10. ){
          std::cout << "Found possible elastic scatter with vertex" << std::endl;
          elastic_candidate = true;
        }
      }
    }


    //Testing
    std::cout << "N Traj Pts: " << true_beam_trajectory.size() << std::endl;
    std::cout << true_beam_particle->TrackId() << std::endl;

    std::cout << "N Reco Traj Pts: " << thisTrack->NumberTrajectoryPoints() << std::endl;
    
    TVector3 start( startX, startY, startZ );
    TVector3 dir( trackDirX, trackDirY, trackDirZ );
    for( size_t i = 0; i < thisTrack->NumberTrajectoryPoints(); ++i ){
      auto pt = thisTrack->Trajectory().LocationAtPoint(i);
      if( ( pt.X() - -999. ) < 1.e-6 ) continue;

      TVector3 p( pt.X(), pt.Y(), pt.Z() );
      double dist = lateralDist( dir, start, p );

      if( dist > max_lateral ) max_lateral = dist;

      if( i < thisTrack->NumberTrajectoryPoints() - 1 ){
        auto next_pt = thisTrack->Trajectory().LocationAtPoint(i+1);
        if( ( next_pt.X() - -999. ) > 1.e-6 ){

          TVector3 next_p( next_pt.X(), next_pt.Y(), next_pt.Z() );
          double segment = ( next_p - p ).Mag();
          if( segment > max_segment ) max_segment = segment;
        }
      }
    }


    //Get the sim::IDEs for each true trajectory point
    std::map< size_t, std::vector< const sim::IDE * > > trueTrajPtsToSimIDEs = truthUtil.GetSimIDEs( *true_beam_particle );
    std::cout << "Traj Pts & Sim IDEs: " << std::endl;
    /*
    for( auto it = trueTrajPtsToSimIDEs.begin(); it != trueTrajPtsToSimIDEs.end(); ++it ){
      std::cout << it->first << " " << it->second.size() << std::endl;
    }
    */


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

    auto traj = thisTrack->Trajectory();
    size_t np = thisTrack->NumberTrajectoryPoints(); 
    traj_cos = ( ( traj.DirectionAtPoint(0).X() * traj.DirectionAtPoint(np-1).X() )
               + ( traj.DirectionAtPoint(0).Y() * traj.DirectionAtPoint(np-1).Y() )
               + ( traj.DirectionAtPoint(0).Z() * traj.DirectionAtPoint(np-1).Z() ) );
    
    traj_cos = traj_cos / sqrt(traj.DirectionAtPoint(0).Mag2() * traj.DirectionAtPoint(np-1).Mag2());    

    std::vector< int > view_0_TPC;
    std::vector< int > view_1_TPC;
    std::vector< int > view_2_TPC;

    for( auto it = trajPtsToHits.begin(); it != trajPtsToHits.end(); ++it ){

      const recob::Hit * theHit = it->second;
      size_t i = it->first;

      double x = thisTrack->Trajectory().LocationAtPoint(i).X();
      double y = thisTrack->Trajectory().LocationAtPoint(i).Y();
      double z = thisTrack->Trajectory().LocationAtPoint(i).Z();

      int slice = std::floor( ( thisTrack->Trajectory().LocationAtPoint(i).Z() - z0 ) / pitch );
      hitsToSlices[ theHit ] = slice;
      slicesToHits[ slice ].push_back( theHit );
      std::cout << "Position: " <<  x << " " << y << " " << z << " Hit Wire: " << theHit->WireID().Wire << " " 
                << theHit->PeakTime()  << " " << theHit->StartTick() << " " << theHit->EndTick() << std::endl;
      std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *theHit );
      for( size_t j = 0; j < ides.size(); ++j ){
        std::cout << "\tIDE: " << ides[j]->trackID << " " << ides[j]->x << " " << ides[j]->y << " " << ides[j]->z << std::endl;
      }


      if( thisTrack->Trajectory().LocationAtPoint(i).Z() > max_Z ){
        max_Z = thisTrack->Trajectory().LocationAtPoint(i).Z();
        max_X = thisTrack->Trajectory().LocationAtPoint(i).X();
        max_Y = thisTrack->Trajectory().LocationAtPoint(i).Y();
      }


      //std::cout << "View: " << theHit->View() << std::endl;
      switch( theHit->View() ){
        case 0: 
          if( theHit->WireID().TPC == 5 )
            view_0_hits_in_TPC5 = true;
          view_0_wire.push_back( theHit->WireID().Wire );
          view_0_tick.push_back( theHit->PeakTime() );
          view_0_TPC.push_back( theHit->WireID().TPC );
          break;
        case 1:
          if( theHit->WireID().TPC == 5 )
            view_1_hits_in_TPC5 = true;
          view_1_wire.push_back( theHit->WireID().Wire );
          view_1_tick.push_back( theHit->PeakTime() );
          view_1_TPC.push_back( theHit->WireID().TPC );
          break;
        case 2: 
          if( theHit->WireID().TPC == 5 )
            view_2_hits_in_TPC5 = true;
          view_2_wire.push_back( theHit->WireID().Wire );
          view_2_z.push_back( thisTrack->Trajectory().LocationAtPoint(i).Z() );
          view_2_tick.push_back( theHit->PeakTime() );
          view_2_TPC.push_back( theHit->WireID().TPC );
          break;
        default:
          break;
      }
    }

 //   std::cout << "View 0: " << view_0_wire.size() << std::endl;
    for( size_t i = 1; i < view_0_wire.size(); ++i ){
 //     std::cout << i << " " << i-1 << std::endl;
      double segment = sqrt( (view_0_wire[i] - view_0_wire[i-1])*(view_0_wire[i] - view_0_wire[i-1]) 
                           + (view_0_tick[i] - view_0_tick[i-1])*(view_0_tick[i] - view_0_tick[i-1]) );
      if( segment > view_0_max_segment ) view_0_max_segment = segment;                           

      if( view_0_wire[i] < view_0_wire[i-1] && ( view_0_TPC[i] != 5 && view_0_TPC[i-1] != 5)  ){
        view_0_wire_backtrack += (view_0_wire[i-1] - view_0_wire[i]);
      }
    }

 //   std::cout << "View 1: " << view_1_wire.size() << std::endl;
    for( size_t i = 1; i < view_1_wire.size(); ++i ){
 //     std::cout << i << " " << i-1 << std::endl;
      double segment = sqrt( (view_1_wire[i] - view_1_wire[i-1])*(view_1_wire[i] - view_1_wire[i-1]) 
                           + (view_1_tick[i] - view_1_tick[i-1])*(view_1_tick[i] - view_1_tick[i-1]) );
      if( segment > view_1_max_segment ) view_1_max_segment = segment;                           

      if( view_1_wire[i] > view_1_wire[i-1]  && ( view_1_TPC[i] != 5 && view_1_TPC[i-1] != 5)){
        view_1_wire_backtrack += (view_1_wire[i] - view_1_wire[i-1]);
      }
    }

//    std::cout << "View 2: " << view_2_wire.size() << std::endl;
    for( size_t i = 1; i < view_2_wire.size(); ++i ){
      //std::cout << i << " " << i-1 << std::endl;
      double segment = sqrt( (view_2_wire[i] - view_2_wire[i-1])*(view_2_wire[i] - view_2_wire[i-1]) 
                           + (view_2_tick[i] - view_2_tick[i-1])*(view_2_tick[i] - view_2_tick[i-1]) );
      if( segment > view_2_max_segment ) view_2_max_segment = segment;                           

      if( view_2_wire[i] < view_2_wire[i-1]  && ( view_2_TPC[i] != 5 && view_2_TPC[i-1] != 5)){
        view_2_wire_backtrack += (view_2_wire[i-1] - view_2_wire[i]);
      }
    }

    vertex_slice = slicesToHits.rbegin()->first;
    std::cout << "Vertex slice: " << vertex_slice << std::endl;

    endZ_slice = std::floor( ( endZ   - z0 ) / pitch );
    startZ_slice = std::floor( ( startZ   - z0 ) / pitch );


    std::cout << "StartZ: " << startZ << " " << std::floor( ( startZ - z0 ) / pitch ) << std::endl;
    std::cout << "EndZ: "   << endZ   << " " << std::floor( ( endZ   - z0 ) / pitch ) << std::endl;


    //Go through the hits in the last slice, then backtrack to the IDs
    std::vector< const recob::Hit * > vertex_hits = slicesToHits.rbegin()->second;
    for( size_t i = 0; i < vertex_hits.size(); ++i ){
      std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(vertex_hits[i]) );
      for( size_t j = 0; j < ides.size(); ++j ){
        vertex_IDs.push_back( ides[j]->trackID );
      }
    }

    vertex_matched_to_beam = ( std::find( vertex_IDs.begin(), vertex_IDs.end(), true_beam_ID ) != vertex_IDs.end() );

    std::vector< const recob::Hit * > endZ_hits = slicesToHits.rbegin()->second;
    for( size_t i = 0; i < endZ_hits.size(); ++i ){
      std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(endZ_hits[i]) );
      for( size_t j = 0; j < ides.size(); ++j ){
        endZ_IDs.push_back( ides[j]->trackID );
      }
    }
    endZ_matched_to_beam = ( std::find( endZ_IDs.begin(), endZ_IDs.end(), true_beam_ID ) != endZ_IDs.end() );

    std::vector< const recob::Hit * > startZ_hits = slicesToHits.rbegin()->second;
    for( size_t i = 0; i < startZ_hits.size(); ++i ){
      std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(startZ_hits[i]) );
      for( size_t j = 0; j < ides.size(); ++j ){
        startZ_IDs.push_back( ides[j]->trackID );
      }
    }
    startZ_matched_to_beam = ( std::find( startZ_IDs.begin(), startZ_IDs.end(), true_beam_ID ) != startZ_IDs.end() );

    //Also, get the distance between all of the IDEs in the last slice to the location of processes in the beam trajectory 

    for( auto itProc = true_beam_proc_map.begin(); itProc != true_beam_proc_map.end(); ++itProc ){
      double procX = true_beam_trajectory.X( itProc->first );
      double procY = true_beam_trajectory.Y( itProc->first );
      double procZ = true_beam_trajectory.Z( itProc->first );

      std::cout << std::endl << "Process: " << true_beam_trajectory.KeyToProcess(itProc->second) << procX << " " << procY << " " << procZ << std::endl; 

      std::vector< double > temp_dRs;

      double centroidX = 0.;
      double centroidY = 0.;
      double centroidZ = 0.;

      int nIDEs = 0;

      for( size_t i = 0; i < vertex_hits.size(); ++i ){
        

        std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(vertex_hits[i]) );
        std::cout << "Hit: " << vertex_hits[i]->WireID().Wire << " " << vertex_hits[i]->WireID().TPC  << " " << vertex_hits[i]->WireID().Plane << std::endl;
        for( size_t j = 0; j < ides.size(); ++j ){
          std::cout << "\tIDE: " << ides[j]->trackID << " " << ides[j]->x << " " << ides[j]->y << " " << ides[j]->z << std::endl;
          temp_dRs.push_back( sqrt( std::pow( (ides[j]->x - procX), 2 ) +
                                    std::pow( (ides[j]->y - procY), 2 ) +
                                    std::pow( (ides[j]->z - procZ), 2 ) ) );
          
          centroidX += ides[j]->x;
          centroidY += ides[j]->y;
          centroidZ += ides[j]->z;
          
          ++nIDEs;
        }
      }
      vertex_dRs.push_back( temp_dRs );

      centroidX = centroidX / nIDEs;
      centroidY = centroidY / nIDEs;
      centroidZ = centroidZ / nIDEs;

      centroid_dR.push_back( sqrt( std::pow( (centroidX - procX), 2 ) +
                                   std::pow( (centroidY - procY), 2 ) +
                                   std::pow( (centroidZ - procZ), 2 ) ) );

      temp_dRs.clear();
      for( size_t i = 0; i < endZ_hits.size(); ++i ){
        std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(endZ_hits[i]) );
        std::vector< double > temp_dRs;
        for( size_t j = 0; j < ides.size(); ++j ){
        temp_dRs.push_back( sqrt( std::pow( (ides[j]->x - procX), 2 ) +
                                  std::pow( (ides[j]->y - procY), 2 ) +
                                  std::pow( (ides[j]->z - procZ), 2 ) ) );
        }
      }
      endZ_dRs.push_back( temp_dRs );

      temp_dRs.clear();
      for( size_t i = 0; i < startZ_hits.size(); ++i ){
        std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(startZ_hits[i]) );
        for( size_t j = 0; j < ides.size(); ++j ){
        temp_dRs.push_back( sqrt( std::pow( (ides[j]->x - procX), 2 ) +
                                  std::pow( (ides[j]->y - procY), 2 ) +
                                  std::pow( (ides[j]->z - procZ), 2 ) ) );
        }
      }
      startZ_dRs.push_back( temp_dRs );

    }

    if( true_beam_EndProcess.find( "Inelastic" ) == std::string::npos ){
      processes.push_back( true_beam_EndProcess );

      double procX = true_beam_EndVertex_X;
      double procY = true_beam_EndVertex_Y;
      double procZ = true_beam_EndVertex_Z;

      std::vector< double > temp_dRs;

      double centroidX = 0.;
      double centroidY = 0.;
      double centroidZ = 0.;

      int nIDEs = 0;

      for( size_t i = 0; i < vertex_hits.size(); ++i ){
        std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(vertex_hits[i]) );
        for( size_t j = 0; j < ides.size(); ++j ){
        temp_dRs.push_back( sqrt( std::pow( (ides[j]->x - procX), 2 ) +
                                  std::pow( (ides[j]->y - procY), 2 ) +
                                  std::pow( (ides[j]->z - procZ), 2 ) ) );

          centroidX += ides[j]->x;
          centroidY += ides[j]->y;
          centroidZ += ides[j]->z;
          
          ++nIDEs;

        }
      }
      vertex_dRs.push_back( temp_dRs );

      centroidX = centroidX / nIDEs;
      centroidY = centroidY / nIDEs;
      centroidZ = centroidZ / nIDEs;

      centroid_dR.push_back( sqrt( std::pow( (centroidX - procX), 2 ) +
                                   std::pow( (centroidY - procY), 2 ) +
                                   std::pow( (centroidZ - procZ), 2 ) ) );


      temp_dRs.clear();
      for( size_t i = 0; i < endZ_hits.size(); ++i ){
        std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(endZ_hits[i]) );
        std::vector< double > temp_dRs;
        for( size_t j = 0; j < ides.size(); ++j ){
        temp_dRs.push_back( sqrt( std::pow( (ides[j]->x - procX), 2 ) +
                                  std::pow( (ides[j]->y - procY), 2 ) +
                                  std::pow( (ides[j]->z - procZ), 2 ) ) );
        }
      }
      endZ_dRs.push_back( temp_dRs );

      temp_dRs.clear();
      for( size_t i = 0; i < startZ_hits.size(); ++i ){
        std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps( *(startZ_hits[i]) );
        for( size_t j = 0; j < ides.size(); ++j ){
        temp_dRs.push_back( sqrt( std::pow( (ides[j]->x - procX), 2 ) +
                                  std::pow( (ides[j]->y - procY), 2 ) +
                                  std::pow( (ides[j]->z - procZ), 2 ) ) );
        }
      }
      startZ_dRs.push_back( temp_dRs );

    }

    //Primary Track Calorimetry 
    std::vector< anab::Calorimetry> calo = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
    std::cout << "N Calos: " << calo.size() << std::endl;
    auto calo_dQdX = calo[0].dQdx();
    auto calo_dEdX = calo[0].dEdx();
    auto calo_range = calo[0].ResidualRange();
    for( size_t i = 0; i < calo_dQdX.size(); ++i ){
      dQdX.push_back( calo_dQdX[i] );
      dEdX.push_back( calo_dEdX[i] );
      resRange.push_back( calo_range[i] );
    }
    ////////////////////////////////////////////



    //Looking at reco daughters from the reco beam track
    const std::vector<const recob::Track*> trackDaughters = pfpUtil.GetPFParticleDaughterTracks(*particle,evt,fPFParticleTag,fTrackerTag);  
    const std::vector<const recob::Shower*> showerDaughters = pfpUtil.GetPFParticleDaughterShowers(*particle,evt,fPFParticleTag,fShowerTag);  
    std::cout << "Beam particle has " << trackDaughters.size() << " track-like daughters and " << showerDaughters.size() << " shower-like daughters." << std::endl;
    std::cout << std::endl;

    std::map< int, std::vector< size_t > > sliceDaughters;

    nTrackDaughters = trackDaughters.size();

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

      double delta_start = sqrt( ( reco_daughter_startX.back() - endX )*( reco_daughter_startX.back() - endX ) + 
                                 ( reco_daughter_startY.back() - endY )*( reco_daughter_startY.back() - endY ) +
                                 ( reco_daughter_startZ.back() - endZ )*( reco_daughter_startZ.back() - endZ ) );

      double delta_end =   sqrt( ( reco_daughter_endX.back() - endX )*( reco_daughter_endX.back() - endX ) + 
                                 ( reco_daughter_endY.back() - endY )*( reco_daughter_endY.back() - endY ) +
                                 ( reco_daughter_endZ.back() - endZ )*( reco_daughter_endZ.back() - endZ ) );

      reco_daughter_deltaR.push_back( ( delta_start < delta_end ) ? delta_start : delta_end );

      reco_daughter_trackID.push_back( daughterTrack->ID() );

      reco_daughter_completeness.push_back( truthUtil.GetCompleteness( *daughterTrack, evt, fTrackerTag, fHitTag ) );
      reco_daughter_purity.push_back( truthUtil.GetPurity( *daughterTrack, evt, fTrackerTag/*, fHitTag*/) );

      //const simb::MCParticle * match = truthUtil.GetMCParticleByHits( *daughterTrack, evt, fTrackerTag, fHitTag );
      protoana::MCParticleSharedHits match = truthUtil.GetMCParticleByHits( *daughterTrack, evt, fTrackerTag, fHitTag );

      if( match.particle ){
        std::cout << std::endl << "Match: " << match.particle->PdgCode() << " " << match.particle->TrackId() << std::endl;
         
        alt_reco_daughter_truth_PDG.push_back( match.particle->PdgCode() );
        alt_reco_daughter_truth_ID.push_back( match.particle->TrackId() );
        alt_reco_daughter_truth_ParID.push_back( match.particle->Mother() );
        alt_reco_daughter_truth_Process.push_back( match.particle->Process() );
        alt_reco_daughter_truth_Origin.push_back( 
          pi_serv->TrackIdToMCTruth_P(match.particle->TrackId())->Origin()
        );
        alt_reco_daughter_truth_shared_hits.push_back( match.nSharedHits ); 
        alt_reco_daughter_truth_em_hits.push_back( match.nSharedDeltaRayHits ); 

        alt_reco_daughter_truth_len.push_back( match.particle->Trajectory().TotalLength() );
        alt_reco_daughter_truth_startX.push_back( match.particle->Position(0).X() );
        alt_reco_daughter_truth_startY.push_back( match.particle->Position(0).Y() );
        alt_reco_daughter_truth_startZ.push_back( match.particle->Position(0).Z() );

        alt_reco_daughter_truth_endX.push_back( match.particle->EndPosition().X() );
        alt_reco_daughter_truth_endY.push_back( match.particle->EndPosition().Y() );
        alt_reco_daughter_truth_endZ.push_back( match.particle->EndPosition().Z() );

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

        alt_reco_daughter_truth_purity.push_back( matched_hits / total );
      }
      else{
        alt_reco_daughter_truth_PDG.push_back( -1 );
        alt_reco_daughter_truth_ID.push_back( -1 );
        alt_reco_daughter_truth_Origin.push_back( -1 );
        alt_reco_daughter_truth_ParID.push_back( -1 );
        alt_reco_daughter_truth_Process.push_back( "empty" );
        alt_reco_daughter_truth_shared_hits.push_back( 0 ); 
        alt_reco_daughter_truth_em_hits.push_back( 0 ); 

        alt_reco_daughter_truth_len.push_back( 0. );
        alt_reco_daughter_truth_startX.push_back( 0. );
        alt_reco_daughter_truth_startY.push_back( 0. );
        alt_reco_daughter_truth_startZ.push_back( 0. );
        alt_reco_daughter_truth_endX.push_back( 0. );
        alt_reco_daughter_truth_endY.push_back( 0. );
        alt_reco_daughter_truth_endZ.push_back( 0. );
        alt_reco_daughter_truth_purity.push_back( 0. );
      }



      //For this reco daughter track, get the actual particle contributing to it
      const simb::MCParticle* trueDaughterParticle = truthUtil.GetMCParticleFromRecoTrack(*daughterTrack, evt, fTrackerTag);
      if( trueDaughterParticle ){
        reco_daughter_truth_PDG.push_back( trueDaughterParticle->PdgCode() );
        reco_daughter_truth_ID.push_back( trueDaughterParticle->TrackId() );
        if( trueDaughterParticle->TrackId() == true_beam_particle->TrackId() ){
          daughter_is_primary = true;
        }
          
        reco_daughter_truth_Origin.push_back( 
          pi_serv->TrackIdToMCTruth_P(trueDaughterParticle->TrackId())->Origin()
        );
        reco_daughter_truth_ParID.push_back( trueDaughterParticle->Mother() );
        reco_daughter_truth_Process.push_back( trueDaughterParticle->Process() );
      }
      else{
        reco_daughter_truth_PDG.push_back( -1 );
        reco_daughter_truth_ID.push_back( -1 );
        reco_daughter_truth_Origin.push_back( -1 );
        reco_daughter_truth_ParID.push_back( -1 );
        reco_daughter_truth_Process.push_back( "empty" );
      }

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
      double none_total = 0.;
      double michel_total = 0.;   
      for( size_t h = 0; h < daughter_hits.size(); ++h ){
        std::array<float,4> cnn_out = hitResults.getOutput( daughter_hits[h ] );
        track_total  += cnn_out[ hitResults.getIndex("track") ];
        em_total     += cnn_out[ hitResults.getIndex("em") ];
        none_total   += cnn_out[ hitResults.getIndex("none") ];
        michel_total += cnn_out[ hitResults.getIndex("michel") ];
      }

      if( daughter_hits.size() > 0 ){
        reco_daughter_track_score.push_back( track_total / daughter_hits.size() );
        reco_daughter_em_score.push_back( em_total / daughter_hits.size() );
        reco_daughter_none_score.push_back( none_total / daughter_hits.size() );
        reco_daughter_michel_score.push_back( michel_total / daughter_hits.size() );
      }
      else{
        reco_daughter_track_score.push_back( -999. );
        reco_daughter_em_score.push_back( -999. );
        reco_daughter_none_score.push_back( -999. );
        reco_daughter_michel_score.push_back( -999. );
      }

      ///Try to match the reconstructed daughter tracks from the reco'd/tagged beam track
      //   to the true daughter particles coming out of the true beam track
      bool found_daughter = false;
      const simb::MCParticle* daughterParticleFromRecoTrack = truthUtil.GetMCParticleFromRecoTrack(*daughterTrack, evt, fTrackerTag); 

      if( daughterParticleFromRecoTrack ){
        int loc = 0;
        for( size_t j = 0; j < true_beam_daughter_IDs.size(); ++j ){
          
          //do the checking
          if ( true_beam_daughter_IDs[j] == daughterParticleFromRecoTrack->TrackId() ){
            found_daughter = true;
            loc = j;
            break;
          }
        }

        reco_beam_truth_daughter_good_reco.push_back( found_daughter );
        if( found_daughter ){
          reco_beam_truth_daughter_true_PDGs.push_back( true_beam_daughter_PDGs[loc] ); 
          reco_beam_truth_daughter_true_IDs.push_back( true_beam_daughter_IDs[loc] );
          reco_beam_truth_daughter_true_lens.push_back( true_beam_daughter_lens[loc] );
        }
        else{

          reco_beam_truth_daughter_true_PDGs.push_back( daughterParticleFromRecoTrack->PdgCode() );
          reco_beam_truth_daughter_true_IDs.push_back( daughterParticleFromRecoTrack->TrackId() );
          reco_beam_truth_daughter_true_lens.push_back( daughterParticleFromRecoTrack->Trajectory().TotalLength() );
        }
      }
      else{
        reco_beam_truth_daughter_good_reco.push_back( false );
        reco_beam_truth_daughter_true_PDGs.push_back( -1 );
        reco_beam_truth_daughter_true_IDs.push_back( -1 );
        reco_beam_truth_daughter_true_lens.push_back( -1 );
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
        
      std::cout << std::endl;
    }

    for( auto it = sliceDaughters.begin(); it != sliceDaughters.end(); ++it ){
      std::cout << "Slice: " << it->first << " Has " << it->second.size() << " Daughters" << std::endl;
    }

    nShowerDaughters = showerDaughters.size();

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

      std::cout << "Looking at shower dedx: ";
      std::cout << daughterShowerFromRecoTrack->dEdx().size() << std::endl;

      //For this reco daughter track, get the actual particle contributing to it
      const simb::MCParticle* trueDaughterParticle = truthUtil.GetMCParticleFromRecoShower(*daughterShowerFromRecoTrack, evt, fShowerTag);
      if( trueDaughterParticle ){
        reco_daughter_shower_truth_PDG.push_back( trueDaughterParticle->PdgCode() );
        reco_daughter_shower_truth_ID.push_back( trueDaughterParticle->TrackId() );
        reco_daughter_shower_truth_Origin.push_back( 
          pi_serv->TrackIdToMCTruth_P(trueDaughterParticle->TrackId())->Origin()
        );
        reco_daughter_shower_truth_ParID.push_back( trueDaughterParticle->Mother() );

      }
      else{
        reco_daughter_shower_truth_PDG.push_back( -1 );
        reco_daughter_shower_truth_ID.push_back( -1 );
        reco_daughter_shower_truth_Origin.push_back( -1 );
        reco_daughter_shower_truth_ParID.push_back( -1 );

      }

      const simb::MCParticle * match = truthUtil.GetMCParticleByHits( *daughterShowerFromRecoTrack, evt, fShowerTag, fHitTag ).particle;

      if( match ){
        std::cout << std::endl << "Match: " << match->PdgCode() << " " << match->TrackId() << std::endl;
         
        alt_reco_daughter_shower_truth_PDG.push_back( match->PdgCode() );
        alt_reco_daughter_shower_truth_ID.push_back( match->TrackId() );
        alt_reco_daughter_shower_truth_ParID.push_back( match->Mother() );
        alt_reco_daughter_shower_truth_Process.push_back( match->Process() );
        alt_reco_daughter_shower_truth_Origin.push_back( 
          pi_serv->TrackIdToMCTruth_P(match->TrackId())->Origin()
        );
      }
      else{
        alt_reco_daughter_shower_truth_PDG.push_back( -1 );
        alt_reco_daughter_shower_truth_ID.push_back( -1 );
        alt_reco_daughter_shower_truth_Origin.push_back( -1 );
        alt_reco_daughter_shower_truth_ParID.push_back( -1 );
        alt_reco_daughter_shower_truth_Process.push_back( "empty" );
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

      alt_reco_daughter_shower_truth_purity.push_back( matched_hits / total );



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
      double none_total = 0.;
      double michel_total = 0.;   
      for( size_t h = 0; h < daughter_hits.size(); ++h ){
        std::array<float,4> cnn_out = hitResults.getOutput( daughter_hits[h ] );
        track_total  += cnn_out[ hitResults.getIndex("track") ];
        em_total     += cnn_out[ hitResults.getIndex("em") ];
        none_total   += cnn_out[ hitResults.getIndex("none") ];
        michel_total += cnn_out[ hitResults.getIndex("michel") ];


        //In addition, try to get the dQ from each hit
        //auto theHit = daughter_hits[ h ];
       // double dQ = theHit->Integral();
        //std::cout << "Hit: " << h << " has dQ " << dQ << std::endl;
        
      }

      if( daughter_hits.size() > 0 ){
        reco_daughter_shower_track_score.push_back( track_total / daughter_hits.size() );
        reco_daughter_shower_em_score.push_back( em_total / daughter_hits.size() );
        reco_daughter_shower_none_score.push_back( none_total / daughter_hits.size() );
        reco_daughter_shower_michel_score.push_back( michel_total / daughter_hits.size() );
      }
      else{
        reco_daughter_shower_track_score.push_back( -999. );
        reco_daughter_shower_em_score.push_back( -999. );
        reco_daughter_shower_none_score.push_back( -999. );
        reco_daughter_shower_michel_score.push_back( -999. );
      }
      



      bool found_daughter = false;
      const simb::MCParticle* daughterParticleFromRecoShower = truthUtil.GetMCParticleFromRecoShower(*daughterShowerFromRecoTrack, evt, fShowerTag); 
      std::cout << "Shower daughter " << i << " " << daughterParticleFromRecoShower << std::endl;

      if( daughterParticleFromRecoShower ){
        std::cout << "Shower daughter ID: " << daughterParticleFromRecoShower->TrackId() << std::endl;
        int loc = 0;
        for( size_t j = 0; j < true_beam_daughter_IDs.size(); ++j ){
          if ( true_beam_daughter_IDs[j] == daughterParticleFromRecoShower->TrackId() ){
            found_daughter = true;
            loc = j;
            break;
          }
        }

        reco_beam_truth_daughter_shower_good_reco.push_back( found_daughter );
        if( found_daughter ){
          reco_beam_truth_daughter_shower_true_PDGs.push_back( true_beam_daughter_PDGs[loc] ); 
          reco_beam_truth_daughter_shower_true_IDs.push_back( true_beam_daughter_IDs[loc] );
          reco_beam_truth_daughter_shower_true_lens.push_back( true_beam_daughter_lens[loc] );
        }
        else{
          reco_beam_truth_daughter_shower_true_PDGs.push_back( daughterParticleFromRecoShower->PdgCode() );
          reco_beam_truth_daughter_shower_true_IDs.push_back( daughterParticleFromRecoShower->TrackId() );
          reco_beam_truth_daughter_shower_true_lens.push_back( daughterParticleFromRecoShower->Trajectory().TotalLength() );
        }
      }

    }

    std::pair< double, int > pid_chi2_ndof = trackUtil.Chi2PID( dEdX, resRange, templates[ 2212 ] );
    reco_beam_Chi2_proton = pid_chi2_ndof.first; 
    reco_beam_Chi2_ndof = pid_chi2_ndof.second;
  
    std::cout << "Proton chi2: " << reco_beam_Chi2_proton << std::endl;


    //Looking at shower/track discrimination
    std::cout << "MVA" << std::endl;

    double track_total = 0.;
    double em_total = 0.;
    double none_total = 0.;
    double michel_total = 0.;   
    std::vector< art::Ptr< recob::Hit > > track_hits = findHits.at( thisTrack->ID() );
    for( size_t h = 0; h < track_hits.size(); ++h ){
      std::array<float,4> cnn_out = hitResults.getOutput( track_hits[h ] );
      track_total  += cnn_out[ hitResults.getIndex("track") ];
      em_total     += cnn_out[ hitResults.getIndex("em") ];
      none_total   += cnn_out[ hitResults.getIndex("none") ];
      michel_total += cnn_out[ hitResults.getIndex("michel") ];
    }
    std::cout << "track Total "  << track_total  << " " << track_hits.size() << std::endl;
    std::cout << "em Total "     << em_total     << " " << track_hits.size() << std::endl;
    std::cout << "none Total "   << none_total   << " " << track_hits.size() << std::endl;
    std::cout << "michel Total " << michel_total << " " << track_hits.size() << std::endl;



    // Alternative Reconstruction.
    //
    // Loop over all of the PFParticles associated as daughters.
    // Then, check the CNN score (later implement the GNN score)
    //
    // Also, get the forced-tracking (both pma and pandora2) and 
    // get calorimetry + other info
    
    auto pfpVec = evt.getValidHandle< std::vector< recob::PFParticle > >( fPFParticleTag );
    

    for( size_t daughterID : particle->Daughters() ){
      const recob::PFParticle * daughterPFP = &(pfpVec->at( daughterID ));
      reco_daughter_PFP_ID.push_back( daughterID );

      const std::vector< art::Ptr< recob::Hit > > daughterPFP_hits = pfpUtil.GetPFParticleHits_Ptrs( *daughterPFP, evt, fPFParticleTag );
      std::cout << "Got " << daughterPFP_hits.size() << " hits from daughter " << daughterID << std::endl;

      double track_total = 0.;
      double em_total = 0.;
      double none_total = 0.;
      double michel_total = 0.;   
      for( size_t h = 0; h < daughterPFP_hits.size(); ++h ){
        std::array<float,4> cnn_out = hitResults.getOutput( daughterPFP_hits[h] );
        track_total  += cnn_out[ hitResults.getIndex("track") ];
        em_total     += cnn_out[ hitResults.getIndex("em") ];
        none_total   += cnn_out[ hitResults.getIndex("none") ];
        michel_total += cnn_out[ hitResults.getIndex("michel") ];


        /*
        const recob::Hit * theHit = daughterPFP_hits[h].get();
        std::cout << theHit->View() << std::endl;
        */
      }

      const std::vector< const recob::SpacePoint* > spVec = pfpUtil.GetPFParticleSpacePoints( *daughterPFP, evt, fPFParticleTag );
      std::cout << "Got " << spVec.size() << " SpacePoints" << std::endl;

      if( daughterPFP_hits.size() > 0 ){
        reco_daughter_PFP_track_score.push_back( track_total / daughterPFP_hits.size() );
        reco_daughter_PFP_em_score.push_back( em_total / daughterPFP_hits.size() );
        reco_daughter_PFP_michel_score.push_back( michel_total / daughterPFP_hits.size() );
        reco_daughter_PFP_none_score.push_back( none_total / daughterPFP_hits.size() );
      }
      else{
        reco_daughter_PFP_track_score.push_back( -999. );
        reco_daughter_PFP_em_score.push_back( -999. );
        reco_daughter_PFP_none_score.push_back( -999. );
        reco_daughter_PFP_michel_score.push_back( -999. );
      }


      protoana::MCParticleSharedHits match = truthUtil.GetMCParticleByHits( *daughterPFP, evt, fPFParticleTag, fHitTag );

      if( match.particle ){
        std::cout << std::endl << "Match: " << match.particle->PdgCode() << " " << match.particle->TrackId() << std::endl;
         
        alt_reco_daughter_PFP_truth_PDG.push_back( match.particle->PdgCode() );
        alt_reco_daughter_PFP_truth_ID.push_back( match.particle->TrackId() );
        alt_reco_daughter_PFP_truth_ParID.push_back( match.particle->Mother() );
        alt_reco_daughter_PFP_truth_Process.push_back( match.particle->Process() );
        alt_reco_daughter_PFP_truth_Origin.push_back( 
          pi_serv->TrackIdToMCTruth_P(match.particle->TrackId())->Origin()
        );
        alt_reco_daughter_PFP_truth_shared_hits.push_back( match.nSharedHits ); 
        alt_reco_daughter_PFP_truth_em_hits.push_back( match.nSharedDeltaRayHits ); 

        alt_reco_daughter_PFP_truth_len.push_back( match.particle->Trajectory().TotalLength() );
        alt_reco_daughter_PFP_truth_startX.push_back( match.particle->Position(0).X() );
        alt_reco_daughter_PFP_truth_startY.push_back( match.particle->Position(0).Y() );
        alt_reco_daughter_PFP_truth_startZ.push_back( match.particle->Position(0).Z() );

        alt_reco_daughter_PFP_truth_endX.push_back( match.particle->EndPosition().X() );
        alt_reco_daughter_PFP_truth_endY.push_back( match.particle->EndPosition().Y() );
        alt_reco_daughter_PFP_truth_endZ.push_back( match.particle->EndPosition().Z() );

      }
      else{
        alt_reco_daughter_PFP_truth_PDG.push_back( -1 );
        alt_reco_daughter_PFP_truth_ID.push_back( -1 );
        alt_reco_daughter_PFP_truth_Origin.push_back( -1 );
        alt_reco_daughter_PFP_truth_ParID.push_back( -1 );
        alt_reco_daughter_PFP_truth_Process.push_back( "empty" );
        alt_reco_daughter_PFP_truth_shared_hits.push_back( 0 ); 
        alt_reco_daughter_PFP_truth_em_hits.push_back( 0 ); 

        alt_reco_daughter_PFP_truth_len.push_back( 0. );
        alt_reco_daughter_PFP_truth_startX.push_back( 0. );
        alt_reco_daughter_PFP_truth_startY.push_back( 0. );
        alt_reco_daughter_PFP_truth_startZ.push_back( 0. );
        alt_reco_daughter_PFP_truth_endX.push_back( 0. );
        alt_reco_daughter_PFP_truth_endY.push_back( 0. );
        alt_reco_daughter_PFP_truth_endZ.push_back( 0. );
      }


      //Now, go through the new reco (if they exist) and get the forced track objects to these.
      try{
        const recob::Track* pmaTrack = pfpUtil.GetPFParticleTrack( *daughterPFP, evt, fPFParticleTag, "pmtrajfit" );
        std::cout << "PMA track: " << pmaTrack << std::endl;

        if( pmaTrack ){ 
          reco_daughter_PMA_ID.push_back( pmaTrack->ID() );

          std::vector< anab::Calorimetry > dummy_calo = trackUtil.GetRecoTrackCalorimetry(*pmaTrack, evt, "pmtrajfit", "pmtrajfitcalo");
          auto dummy_dEdx = dummy_calo[0].dEdx();
          auto dummy_Range = dummy_calo[0].ResidualRange();
 
          reco_daughter_PMA_resRange.push_back( std::vector<double>() );
          reco_daughter_PMA_dEdX.push_back( std::vector<double>() );

          for( size_t j = 0; j < dummy_dEdx.size(); ++j ){
            reco_daughter_PMA_resRange.back().push_back( dummy_Range[j] );
            reco_daughter_PMA_dEdX.back().push_back( dummy_dEdx[j] );
          }

          std::pair< double, int > this_chi2_ndof = trackUtil.Chi2PID( reco_daughter_PMA_dEdX.back(), reco_daughter_PMA_resRange.back(), templates[ 2212 ] );
          reco_daughter_PMA_Chi2_proton.push_back( this_chi2_ndof.first );
          reco_daughter_PMA_Chi2_ndof.push_back( this_chi2_ndof.second );

          reco_daughter_PMA_startX.push_back( pmaTrack->Trajectory().Start().X() );
          reco_daughter_PMA_startY.push_back( pmaTrack->Trajectory().Start().Y() );
          reco_daughter_PMA_startZ.push_back( pmaTrack->Trajectory().Start().Z() );
          reco_daughter_PMA_endX.push_back(   pmaTrack->Trajectory().End().X() );
          reco_daughter_PMA_endY.push_back(   pmaTrack->Trajectory().End().Y() );
          reco_daughter_PMA_endZ.push_back(   pmaTrack->Trajectory().End().Z() );

          //Match the daughters to a slice
          //First, check whether the start or end of the daughter track are closer
          double d_startX = pmaTrack->Trajectory().Start().X();
          double d_startY = pmaTrack->Trajectory().Start().Y();
          double d_startZ = pmaTrack->Trajectory().Start().Z();

          double d_endX = pmaTrack->Trajectory().End().X();
          double d_endY = pmaTrack->Trajectory().End().Y();
          double d_endZ = pmaTrack->Trajectory().End().Z();

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
            reco_daughter_PMA_to_vertex.push_back( to_end_of_daughter );
          }
          else{
            reco_daughter_PMA_to_vertex.push_back( to_start_of_daughter );
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
           // min_index = min_end_index; 

            std::cout << "dr between track and daughter: " << dr_end << std::endl;
            reco_daughter_PMA_dR.push_back( dr_end );
 
          }
          else{
           // min_index = min_start_index; 
            std::cout << "dr between track and daughter: " << dr_start << std::endl;
            reco_daughter_PMA_dR.push_back( dr_start );
          }
        }
        else{
          reco_daughter_PMA_ID.push_back( -1 );
          reco_daughter_PMA_resRange.push_back( std::vector<double>() );
          reco_daughter_PMA_dEdX.push_back( std::vector<double>() );
          reco_daughter_PMA_Chi2_proton.push_back( -999. );
          reco_daughter_PMA_Chi2_ndof.push_back( 0 );
          reco_daughter_PMA_startX.push_back( -1. );
          reco_daughter_PMA_startY.push_back( -1. );
          reco_daughter_PMA_startZ.push_back( -1. );
          reco_daughter_PMA_endX.push_back(   -1. );
          reco_daughter_PMA_endY.push_back(   -1. );
          reco_daughter_PMA_endZ.push_back(   -1. );
          reco_daughter_PMA_to_vertex.push_back( -1. );
          reco_daughter_PMA_dR.push_back( -1. );
        }


      }
      catch( const cet::exception &e ){
        std::cout << "pmtrajfit object not found, moving on" << std::endl;
      }

      try{
        const recob::Track* pandora2Track = pfpUtil.GetPFParticleTrack( *daughterPFP, evt, fPFParticleTag, "pandora2Track" );
        std::cout << "pandora2 track: " << pandora2Track << std::endl;

        if( pandora2Track ){
          reco_daughter_pandora2_ID.push_back( pandora2Track->ID() );

          std::vector< anab::Calorimetry > dummy_calo = trackUtil.GetRecoTrackCalorimetry(*pandora2Track, evt, "pandora2Track", "pandora2calo");
          auto dummy_dEdx = dummy_calo[0].dEdx();
          auto dummy_Range = dummy_calo[0].ResidualRange();
 
          reco_daughter_pandora2_resRange.push_back( std::vector<double>() );
          reco_daughter_pandora2_dEdX.push_back( std::vector<double>() );

          for( size_t j = 0; j < dummy_dEdx.size(); ++j ){
            reco_daughter_pandora2_resRange.back().push_back( dummy_Range[j] );
            reco_daughter_pandora2_dEdX.back().push_back( dummy_dEdx[j] );
          }

          std::pair< double, int > this_chi2_ndof = trackUtil.Chi2PID( reco_daughter_pandora2_dEdX.back(), reco_daughter_pandora2_resRange.back(), templates[ 2212 ] );
          reco_daughter_pandora2_Chi2_proton.push_back( this_chi2_ndof.first );
          reco_daughter_pandora2_Chi2_ndof.push_back( this_chi2_ndof.second );

          reco_daughter_pandora2_startX.push_back( pandora2Track->Trajectory().Start().X() );
          reco_daughter_pandora2_startY.push_back( pandora2Track->Trajectory().Start().Y() );
          reco_daughter_pandora2_startZ.push_back( pandora2Track->Trajectory().Start().Z() );
          reco_daughter_pandora2_endX.push_back(   pandora2Track->Trajectory().End().X() );
          reco_daughter_pandora2_endY.push_back(   pandora2Track->Trajectory().End().Y() );
          reco_daughter_pandora2_endZ.push_back(   pandora2Track->Trajectory().End().Z() );

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
            reco_daughter_pandora2_to_vertex.push_back( to_end_of_daughter );
          }
          else{
            reco_daughter_pandora2_to_vertex.push_back( to_start_of_daughter );
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
            reco_daughter_pandora2_dR.push_back( dr_end );
 
          }
          else{
            //min_index = min_start_index; 
            std::cout << "dr between track and daughter: " << dr_start << std::endl;
            reco_daughter_pandora2_dR.push_back( dr_start );
          }

        }
        else{
          reco_daughter_pandora2_ID.push_back( -1 );
          reco_daughter_pandora2_resRange.push_back( std::vector<double>() );
          reco_daughter_pandora2_dEdX.push_back( std::vector<double>() );
          reco_daughter_pandora2_Chi2_proton.push_back( -999. );
          reco_daughter_pandora2_Chi2_ndof.push_back( 0 );
          reco_daughter_pandora2_startX.push_back( -1. );
          reco_daughter_pandora2_startY.push_back( -1. );
          reco_daughter_pandora2_startZ.push_back( -1. );
          reco_daughter_pandora2_endX.push_back(   -1. );
          reco_daughter_pandora2_endY.push_back(   -1. );
          reco_daughter_pandora2_endZ.push_back(   -1. );
          reco_daughter_pandora2_to_vertex.push_back( -1. );
          reco_daughter_pandora2_dR.push_back( -1. );
        }
      }
      catch( const cet::exception &e ){
        std::cout << "pandora2Track object not found, moving on" << std::endl;
      }

      try{
        const recob::Shower* pandora2Shower = pfpUtil.GetPFParticleShower( *daughterPFP, evt, fPFParticleTag, "pandora2Shower" );
        std::cout << "pandora2 shower: " << pandora2Shower << std::endl;
      }
      catch( const cet::exception &e ){
        std::cout << "pandora2Shower object not found, moving on" << std::endl;
      }

      
    }



  }
  else if( thisShower ){
    beamTrackID = thisShower->ID();
    std::cout << thisShower->ShowerStart().X() << " " << thisShower->ShowerStart().Y() << " " << thisShower->ShowerStart().Z() << std::endl;
    std::cout << thisShower->Direction().X() << " " << thisShower->Direction().Y() << " " << thisShower->Direction().Z() << std::endl;
    std::cout << beam_cuts.IsBeamlike( *thisShower, evt, "1" ) << std::endl;
  }



  fTree->Fill();
}

void pionana::PionAnalyzerMC::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("beamana","beam analysis tree");

  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("type", &type);
  fTree->Branch("MC", &MC);


  fTree->Branch("startX", &startX);
  fTree->Branch("startY", &startY);
  fTree->Branch("startZ", &startZ);
  fTree->Branch("endX", &endX);
  fTree->Branch("endY", &endY);
  fTree->Branch("endZ", &endZ);
  fTree->Branch("len", &len);
  fTree->Branch("trackDirX", &trackDirX);
  fTree->Branch("trackDirY", &trackDirY);
  fTree->Branch("trackDirZ", &trackDirZ);
  fTree->Branch("trackEndDirX", &trackEndDirX);
  fTree->Branch("trackEndDirY", &trackEndDirY);
  fTree->Branch("trackEndDirZ", &trackEndDirZ);
  fTree->Branch("vtxX", &vtxX);
  fTree->Branch("vtxY", &vtxY);
  fTree->Branch("vtxZ", &vtxZ);
  fTree->Branch("beamTrackID", &beamTrackID);
  fTree->Branch("dQdX", &dQdX);
  fTree->Branch("dEdX", &dEdX);
  fTree->Branch("resRange", &resRange);
  fTree->Branch("nBeamParticles", &nBeamParticles);


  fTree->Branch("reco_daughter_trackID", &reco_daughter_trackID);
  fTree->Branch("reco_daughter_completeness", &reco_daughter_completeness);
  fTree->Branch("reco_daughter_purity", &reco_daughter_purity);
  fTree->Branch("reco_daughter_truth_PDG", &reco_daughter_truth_PDG);
  fTree->Branch("reco_daughter_truth_ID", &reco_daughter_truth_ID);
  fTree->Branch("reco_daughter_truth_Origin", &reco_daughter_truth_Origin);
  fTree->Branch("reco_daughter_truth_ParID", &reco_daughter_truth_ParID);
  fTree->Branch("reco_daughter_truth_Process", &reco_daughter_truth_Process);

  fTree->Branch("alt_reco_daughter_truth_PDG", &alt_reco_daughter_truth_PDG);
  fTree->Branch("alt_reco_daughter_truth_ID", &alt_reco_daughter_truth_ID);
  fTree->Branch("alt_reco_daughter_truth_Origin", &alt_reco_daughter_truth_Origin);
  fTree->Branch("alt_reco_daughter_truth_ParID", &alt_reco_daughter_truth_ParID);
  fTree->Branch("alt_reco_daughter_truth_Process", &alt_reco_daughter_truth_Process);
  fTree->Branch("alt_reco_daughter_truth_purity", &alt_reco_daughter_truth_purity);
  fTree->Branch("alt_reco_daughter_truth_shared_hits", &alt_reco_daughter_truth_shared_hits);
  fTree->Branch("alt_reco_daughter_truth_em_hits", &alt_reco_daughter_truth_em_hits);

  fTree->Branch("alt_reco_daughter_truth_len", &alt_reco_daughter_truth_len);
  fTree->Branch("alt_reco_daughter_truth_startX", &alt_reco_daughter_truth_startX);
  fTree->Branch("alt_reco_daughter_truth_startY", &alt_reco_daughter_truth_startY);
  fTree->Branch("alt_reco_daughter_truth_startZ", &alt_reco_daughter_truth_startZ);
  fTree->Branch("alt_reco_daughter_truth_endX", &alt_reco_daughter_truth_endX);
  fTree->Branch("alt_reco_daughter_truth_endY", &alt_reco_daughter_truth_endY);
  fTree->Branch("alt_reco_daughter_truth_endZ", &alt_reco_daughter_truth_endZ);

  //Alternative reco
  fTree->Branch("alt_reco_daughter_PFP_truth_PDG", &alt_reco_daughter_PFP_truth_PDG);
  fTree->Branch("alt_reco_daughter_PFP_truth_ID", &alt_reco_daughter_PFP_truth_ID);
  fTree->Branch("alt_reco_daughter_PFP_truth_Origin", &alt_reco_daughter_PFP_truth_Origin);
  fTree->Branch("alt_reco_daughter_PFP_truth_ParID", &alt_reco_daughter_PFP_truth_ParID);
  fTree->Branch("alt_reco_daughter_PFP_truth_Process", &alt_reco_daughter_PFP_truth_Process);
  fTree->Branch("alt_reco_daughter_PFP_truth_shared_hits", &alt_reco_daughter_PFP_truth_shared_hits);
  fTree->Branch("alt_reco_daughter_PFP_truth_em_hits", &alt_reco_daughter_PFP_truth_em_hits);

  fTree->Branch("alt_reco_daughter_PFP_truth_len", &alt_reco_daughter_PFP_truth_len);
  fTree->Branch("alt_reco_daughter_PFP_truth_startX", &alt_reco_daughter_PFP_truth_startX);
  fTree->Branch("alt_reco_daughter_PFP_truth_startY", &alt_reco_daughter_PFP_truth_startY);
  fTree->Branch("alt_reco_daughter_PFP_truth_startZ", &alt_reco_daughter_PFP_truth_startZ);
  fTree->Branch("alt_reco_daughter_PFP_truth_endX", &alt_reco_daughter_PFP_truth_endX);
  fTree->Branch("alt_reco_daughter_PFP_truth_endY", &alt_reco_daughter_PFP_truth_endY);
  fTree->Branch("alt_reco_daughter_PFP_truth_endZ", &alt_reco_daughter_PFP_truth_endZ);

  fTree->Branch("reco_daughter_PMA_ID", &reco_daughter_PMA_ID);
  fTree->Branch("reco_daughter_PMA_dEdX", &reco_daughter_PMA_dEdX);
  fTree->Branch("reco_daughter_PMA_resRange", &reco_daughter_PMA_resRange);
  fTree->Branch("reco_daughter_PMA_Chi2_proton", &reco_daughter_PMA_Chi2_proton);
  fTree->Branch("reco_daughter_PMA_Chi2_ndof", &reco_daughter_PMA_Chi2_ndof);

  fTree->Branch("reco_daughter_pandora2_ID", &reco_daughter_pandora2_ID);
  fTree->Branch("reco_daughter_pandora2_dEdX", &reco_daughter_pandora2_dEdX);
  fTree->Branch("reco_daughter_pandora2_resRange", &reco_daughter_pandora2_resRange);
  fTree->Branch("reco_daughter_pandora2_Chi2_proton", &reco_daughter_pandora2_Chi2_proton);
  fTree->Branch("reco_daughter_pandora2_Chi2_ndof", &reco_daughter_pandora2_Chi2_ndof);

  fTree->Branch("reco_daughter_PMA_startX", &reco_daughter_PMA_startX);
  fTree->Branch("reco_daughter_PMA_startY", &reco_daughter_PMA_startY);
  fTree->Branch("reco_daughter_PMA_startZ", &reco_daughter_PMA_startZ);
  fTree->Branch("reco_daughter_PMA_endX", &reco_daughter_PMA_endX);
  fTree->Branch("reco_daughter_PMA_endY", &reco_daughter_PMA_endY);
  fTree->Branch("reco_daughter_PMA_endZ", &reco_daughter_PMA_endZ);
  fTree->Branch("reco_daughter_PMA_dR", &reco_daughter_PMA_dR);
  fTree->Branch("reco_daughter_PMA_to_vertex", &reco_daughter_PMA_to_vertex);

  fTree->Branch("reco_daughter_pandora2_startX", &reco_daughter_pandora2_startX);
  fTree->Branch("reco_daughter_pandora2_startY", &reco_daughter_pandora2_startY);
  fTree->Branch("reco_daughter_pandora2_startZ", &reco_daughter_pandora2_startZ);
  fTree->Branch("reco_daughter_pandora2_endX", &reco_daughter_pandora2_endX);
  fTree->Branch("reco_daughter_pandora2_endY", &reco_daughter_pandora2_endY);
  fTree->Branch("reco_daughter_pandora2_endZ", &reco_daughter_pandora2_endZ);
  fTree->Branch("reco_daughter_pandora2_dR", &reco_daughter_pandora2_dR);
  fTree->Branch("reco_daughter_pandora2_to_vertex", &reco_daughter_pandora2_to_vertex);
  //////



  fTree->Branch("reco_daughter_shower_truth_PDG", &reco_daughter_shower_truth_PDG);
  fTree->Branch("reco_daughter_shower_truth_ID", &reco_daughter_shower_truth_ID);
  fTree->Branch("reco_daughter_shower_truth_Origin", &reco_daughter_shower_truth_Origin);
  fTree->Branch("reco_daughter_shower_truth_ParID", &reco_daughter_shower_truth_ParID);

  fTree->Branch("alt_reco_daughter_shower_truth_PDG", &alt_reco_daughter_shower_truth_PDG);
  fTree->Branch("alt_reco_daughter_shower_truth_ID", &alt_reco_daughter_shower_truth_ID);
  fTree->Branch("alt_reco_daughter_shower_truth_Origin", &alt_reco_daughter_shower_truth_Origin);
  fTree->Branch("alt_reco_daughter_shower_truth_ParID", &alt_reco_daughter_shower_truth_ParID);
  fTree->Branch("alt_reco_daughter_shower_truth_Process", &alt_reco_daughter_shower_truth_Process);
  fTree->Branch("alt_reco_daughter_shower_truth_purity", &alt_reco_daughter_shower_truth_purity);

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
  fTree->Branch("nTrackDaughters", &nTrackDaughters);
  fTree->Branch("nShowerDaughters", &nShowerDaughters);


  fTree->Branch("reco_daughter_PFP_ID", &reco_daughter_PFP_ID);
  fTree->Branch("reco_daughter_PFP_track_score", &reco_daughter_PFP_track_score);
  fTree->Branch("reco_daughter_PFP_em_score", &reco_daughter_PFP_em_score);
  fTree->Branch("reco_daughter_PFP_michel_score", &reco_daughter_PFP_michel_score);
  fTree->Branch("reco_daughter_PFP_none_score", &reco_daughter_PFP_none_score);


  fTree->Branch("true_beam_PDG", &true_beam_PDG);
  fTree->Branch("true_beam_ID", &true_beam_ID);
  fTree->Branch("true_beam_EndProcess", &true_beam_EndProcess);
  fTree->Branch("true_beam_EndVertex_X", &true_beam_EndVertex_X);
  fTree->Branch("true_beam_EndVertex_Y", &true_beam_EndVertex_Y);
  fTree->Branch("true_beam_EndVertex_Z", &true_beam_EndVertex_Z);
  fTree->Branch("true_beam_Start_X", &true_beam_Start_X);
  fTree->Branch("true_beam_Start_Y", &true_beam_Start_Y);
  fTree->Branch("true_beam_Start_Z", &true_beam_Start_Z);

  fTree->Branch("true_beam_Start_Px", &true_beam_Start_Px);
  fTree->Branch("true_beam_Start_Py", &true_beam_Start_Py);
  fTree->Branch("true_beam_Start_Pz", &true_beam_Start_Pz);
  fTree->Branch("true_beam_Start_P", &true_beam_Start_P);

  fTree->Branch("true_beam_End_Px", &true_beam_End_Px);
  fTree->Branch("true_beam_End_Py", &true_beam_End_Py);
  fTree->Branch("true_beam_End_Pz", &true_beam_End_Pz);
  fTree->Branch("true_beam_End_P", &true_beam_End_P);

  fTree->Branch("true_beam_Start_DirX", &true_beam_Start_DirX);
  fTree->Branch("true_beam_Start_DirY", &true_beam_Start_DirY);
  fTree->Branch("true_beam_Start_DirZ", &true_beam_Start_DirZ);

  fTree->Branch("nPi0_truth", &nPi0_truth);
  fTree->Branch("nPiPlus_truth", &nPiPlus_truth);
  fTree->Branch("nProton_truth", &nProton_truth);
  fTree->Branch("nNeutron_truth", &nNeutron_truth);
  fTree->Branch("nPiMinus_truth", &nPiMinus_truth);
  fTree->Branch("nNucleus_truth", &nNucleus_truth);

  fTree->Branch("vertex_type", &vertex_type);
  fTree->Branch("vertex_slice", &vertex_slice);
  fTree->Branch("endZ_slice", &endZ_slice);
  fTree->Branch("startZ_slice", &startZ_slice);
  fTree->Branch("vertex_matched_to_beam", &vertex_matched_to_beam);
  fTree->Branch("endZ_matched_to_beam", &endZ_matched_to_beam);
  fTree->Branch("startZ_matched_to_beam", &startZ_matched_to_beam);

  fTree->Branch("vertex_IDs", &vertex_IDs);
  fTree->Branch("endZ_IDs", &endZ_IDs);
  fTree->Branch("startZ_IDs", &startZ_IDs);

  fTree->Branch("vertex_dRs", &vertex_dRs);
  fTree->Branch("endZ_dRs", &endZ_dRs);
  fTree->Branch("startZ_dRs", &startZ_dRs);

  fTree->Branch("centroid_dR", &centroid_dR);

  fTree->Branch("true_beam_daughter_PDGs", &true_beam_daughter_PDGs);
  fTree->Branch("true_beam_daughter_IDs", &true_beam_daughter_IDs);
  fTree->Branch("true_beam_daughter_lens", &true_beam_daughter_lens);
  fTree->Branch("true_beam_daughter_startX", &true_beam_daughter_startX);
  fTree->Branch("true_beam_daughter_startY", &true_beam_daughter_startY);
  fTree->Branch("true_beam_daughter_startZ", &true_beam_daughter_startZ);
  fTree->Branch("true_beam_daughter_startPX", &true_beam_daughter_startPX);
  fTree->Branch("true_beam_daughter_startPY", &true_beam_daughter_startPY);
  fTree->Branch("true_beam_daughter_startPZ", &true_beam_daughter_startPZ);
  fTree->Branch("true_beam_daughter_startP", &true_beam_daughter_startP);
  fTree->Branch("true_beam_daughter_Process", &true_beam_daughter_Process);

  fTree->Branch("true_beam_Pi0_decay_IDs", &true_beam_Pi0_decay_IDs);
  fTree->Branch("true_beam_Pi0_decay_PDGs", &true_beam_Pi0_decay_PDGs);
  fTree->Branch("true_beam_grand_daughter_IDs", &true_beam_grand_daughter_IDs);
  fTree->Branch("true_beam_grand_daughter_ParIDs", &true_beam_grand_daughter_ParIDs);
  fTree->Branch("true_beam_grand_daughter_PDGs", &true_beam_grand_daughter_PDGs);

  fTree->Branch("reco_beam_truth_EndProcess", &reco_beam_truth_EndProcess);
  fTree->Branch("reco_beam_truth_Process", &reco_beam_truth_Process);
  fTree->Branch("reco_beam_truth_origin", &reco_beam_truth_origin);
  fTree->Branch("reco_beam_truth_PDG", &reco_beam_truth_PDG);
  fTree->Branch("reco_beam_truth_ID", &reco_beam_truth_ID);

  fTree->Branch("alt_reco_beam_truth_EndProcess", &alt_reco_beam_truth_EndProcess);
  fTree->Branch("alt_reco_beam_truth_Process", &alt_reco_beam_truth_Process);
  fTree->Branch("alt_reco_beam_truth_origin", &alt_reco_beam_truth_origin);
  fTree->Branch("alt_reco_beam_truth_PDG", &alt_reco_beam_truth_PDG);
  fTree->Branch("alt_reco_beam_truth_ID", &alt_reco_beam_truth_ID);

  fTree->Branch("reco_beam_good", &reco_beam_good);
  fTree->Branch("alt_reco_beam_good", &alt_reco_beam_good);
  fTree->Branch("passes_beam_cuts", &passes_beam_cuts);
  fTree->Branch("alt_reco_beam_purity", &alt_reco_beam_purity);
  fTree->Branch("flipped", &flipped);
  fTree->Branch("elastic_candidate", &elastic_candidate);
  fTree->Branch("processes", &processes);
  fTree->Branch("found_elastic_scatter", &found_elastic_scatter);
  fTree->Branch("nElasticScatters", &nElasticScatters);
  fTree->Branch("elastic_costheta", &elastic_costheta);
  fTree->Branch("elastic_X", &elastic_X);
  fTree->Branch("elastic_Y", &elastic_Y);
  fTree->Branch("elastic_Z", &elastic_Z);
  fTree->Branch("nReconstructed", &nReconstructed);
  fTree->Branch("reconstructedIDs", &reconstructedIDs);
  fTree->Branch("daughter_is_primary", &daughter_is_primary);

  fTree->Branch("traj_cos", &traj_cos);
  fTree->Branch("view_0_hits_in_TPC5", &view_0_hits_in_TPC5);
  fTree->Branch("view_1_hits_in_TPC5", &view_1_hits_in_TPC5);
  fTree->Branch("view_2_hits_in_TPC5", &view_2_hits_in_TPC5);
  fTree->Branch("max_lateral", &max_lateral);
  fTree->Branch("max_segment", &max_segment);
  fTree->Branch("view_0_max_segment", &view_0_max_segment);
  fTree->Branch("view_1_max_segment", &view_1_max_segment);
  fTree->Branch("view_2_max_segment", &view_2_max_segment);

  fTree->Branch("view_0_wire_backtrack", &view_0_wire_backtrack);
  fTree->Branch("view_1_wire_backtrack", &view_1_wire_backtrack);
  fTree->Branch("view_2_wire_backtrack", &view_2_wire_backtrack);

  fTree->Branch("view_0_wire", &view_0_wire);
  fTree->Branch("view_1_wire", &view_1_wire);
  fTree->Branch("view_2_wire", &view_2_wire);

  fTree->Branch("view_2_z", &view_2_z);

  fTree->Branch("view_0_tick", &view_0_tick);
  fTree->Branch("view_1_tick", &view_1_tick);
  fTree->Branch("view_2_tick", &view_2_tick);

  fTree->Branch("reco_beam_Chi2_proton", &reco_beam_Chi2_proton);
  fTree->Branch("reco_beam_Chi2_ndof", &reco_beam_Chi2_ndof);

  fTree->Branch("reco_daughter_Chi2_proton", &reco_daughter_Chi2_proton);
  fTree->Branch("reco_daughter_Chi2_ndof", &reco_daughter_Chi2_ndof);

  fTree->Branch("reco_daughter_shower_Chi2_proton", &reco_daughter_shower_Chi2_proton);
  fTree->Branch("reco_daughter_shower_Chi2_ndof", &reco_daughter_shower_Chi2_ndof);

  fTree->Branch("reco_daughter_track_score", &reco_daughter_track_score);
  fTree->Branch("reco_daughter_em_score", &reco_daughter_em_score);
  fTree->Branch("reco_daughter_none_score", &reco_daughter_none_score);
  fTree->Branch("reco_daughter_michel_score", &reco_daughter_michel_score);

  fTree->Branch("reco_daughter_shower_track_score", &reco_daughter_shower_track_score);
  fTree->Branch("reco_daughter_shower_em_score", &reco_daughter_shower_em_score);
  fTree->Branch("reco_daughter_shower_none_score", &reco_daughter_shower_none_score);
  fTree->Branch("reco_daughter_shower_michel_score", &reco_daughter_shower_michel_score);

  fTree->Branch("reco_beam_truth_End_Px", &reco_beam_truth_End_Px);
  fTree->Branch("reco_beam_truth_End_Py", &reco_beam_truth_End_Py);
  fTree->Branch("reco_beam_truth_End_Pz", &reco_beam_truth_End_Pz);
  fTree->Branch("reco_beam_truth_End_E", &reco_beam_truth_End_E);
  fTree->Branch("reco_beam_truth_End_P", &reco_beam_truth_End_P);

  fTree->Branch("reco_beam_truth_Start_Px", &reco_beam_truth_Start_Px);
  fTree->Branch("reco_beam_truth_Start_Py", &reco_beam_truth_Start_Py);
  fTree->Branch("reco_beam_truth_Start_Pz", &reco_beam_truth_Start_Pz);
  fTree->Branch("reco_beam_truth_Start_E", &reco_beam_truth_Start_E);
  fTree->Branch("reco_beam_truth_Start_P", &reco_beam_truth_Start_P);


  fTree->Branch("alt_reco_beam_truth_End_Px", &alt_reco_beam_truth_End_Px);
  fTree->Branch("alt_reco_beam_truth_End_Py", &alt_reco_beam_truth_End_Py);
  fTree->Branch("alt_reco_beam_truth_End_Pz", &alt_reco_beam_truth_End_Pz);
  fTree->Branch("alt_reco_beam_truth_End_E", &alt_reco_beam_truth_End_E);
  fTree->Branch("alt_reco_beam_truth_End_P", &alt_reco_beam_truth_End_P);

  fTree->Branch("alt_reco_beam_truth_Start_Px", &alt_reco_beam_truth_Start_Px);
  fTree->Branch("alt_reco_beam_truth_Start_Py", &alt_reco_beam_truth_Start_Py);
  fTree->Branch("alt_reco_beam_truth_Start_Pz", &alt_reco_beam_truth_Start_Pz);
  fTree->Branch("alt_reco_beam_truth_Start_E", &alt_reco_beam_truth_Start_E);
  fTree->Branch("alt_reco_beam_truth_Start_P", &alt_reco_beam_truth_Start_P);

  fTree->Branch("reco_beam_truth_daughter_good_reco", &reco_beam_truth_daughter_good_reco);
  fTree->Branch("reco_beam_truth_daughter_true_PDGs", &reco_beam_truth_daughter_true_PDGs);
  fTree->Branch("reco_beam_truth_daughter_true_IDs", &reco_beam_truth_daughter_true_IDs);

  fTree->Branch("reco_beam_truth_daughter_shower_good_reco", &reco_beam_truth_daughter_shower_good_reco);
  fTree->Branch("reco_beam_truth_daughter_shower_true_PDGs", &reco_beam_truth_daughter_shower_true_PDGs);
  fTree->Branch("reco_beam_truth_daughter_shower_true_IDs", &reco_beam_truth_daughter_shower_true_IDs);

  fTree->Branch("reco_beam_truth_daughter_true_lens", &reco_beam_truth_daughter_true_lens);
  fTree->Branch("reco_beam_truth_daughter_shower_true_lens", &reco_beam_truth_daughter_shower_true_lens);
  fTree->Branch("reco_beam_truth_daughter_shower_good_reco", &reco_beam_truth_daughter_shower_good_reco);

}

void pionana::PionAnalyzerMC::endJob()
{
  dEdX_template_file.Close();
}

double pionana::PionAnalyzerMC::lateralDist(TVector3 &n, TVector3 &x0, TVector3 &p){
  TVector3 x = ( (p - x0)*n )*n;
  return (x - (p - x0)).Mag();
}

void pionana::PionAnalyzerMC::reset()
{
  startX = -1;
  startY = -1;
  startZ = -1;
  endX = -1;
  endY = -1;
  endZ = -1;

  len = -1;
  type = -1;
  nBeamParticles = 0;
  
  MC = 0;
  nProton_truth = 0;
  nNeutron_truth = 0;
  nNucleus_truth = 0;
  vertex_type = kUnmatched;
  vertex_slice = std::numeric_limits<int>::max();
  endZ_slice = std::numeric_limits<int>::max();
  startZ_slice = std::numeric_limits<int>::max();

  vertex_matched_to_beam = false; 
  endZ_matched_to_beam   = false; 
  startZ_matched_to_beam = false; 

  vertex_IDs.clear();
  endZ_IDs.clear();
  startZ_IDs.clear();

  vertex_dRs.clear();
  endZ_dRs.clear();
  startZ_dRs.clear();

  centroid_dR.clear();

  nPi0_truth = 0;
  nPiPlus_truth = 0;
  nPiMinus_truth = 0;
  reco_beam_truth_PDG = 0;
  reco_beam_truth_ID = 0;
  alt_reco_beam_truth_PDG = 0;
  alt_reco_beam_truth_ID = 0;

  true_beam_PDG = 0;
  true_beam_ID = 0;
  true_beam_EndProcess ="";
  true_beam_EndVertex_X = 0.;
  true_beam_EndVertex_Y = 0.;
  true_beam_EndVertex_Z = 0.;
  true_beam_Start_X = 0.;
  true_beam_Start_Y = 0.;
  true_beam_Start_Z = 0.;

  true_beam_Start_Px   = 0.; 
  true_beam_Start_Py   = 0.; 
  true_beam_Start_Pz   = 0.; 
  true_beam_Start_P    = 0.; 

  true_beam_End_Px   = 0.; 
  true_beam_End_Py   = 0.; 
  true_beam_End_Pz   = 0.; 
  true_beam_End_P    = 0.; 

  true_beam_Start_DirX = 0.; 
  true_beam_Start_DirY = 0.; 
  true_beam_Start_DirZ = 0.; 

  reco_beam_truth_EndProcess ="";
  reco_beam_truth_Process ="";
  reco_beam_truth_origin = -1;

  reco_beam_truth_End_Px = 0.;
  reco_beam_truth_End_Py = 0.;
  reco_beam_truth_End_Pz = 0.;
  reco_beam_truth_End_E = 0.;
  reco_beam_truth_End_P = 0.;

  reco_beam_truth_Start_Px = 0.;
  reco_beam_truth_Start_Py = 0.;
  reco_beam_truth_Start_Pz = 0.;
  reco_beam_truth_Start_E = 0.;
  reco_beam_truth_Start_P = 0.;

  alt_reco_beam_truth_EndProcess ="";
  alt_reco_beam_truth_Process ="";
  alt_reco_beam_truth_origin = -1;

  alt_reco_beam_truth_End_Px = 0.;
  alt_reco_beam_truth_End_Py = 0.;
  alt_reco_beam_truth_End_Pz = 0.;
  alt_reco_beam_truth_End_E = 0.;
  alt_reco_beam_truth_End_P = 0.;

  alt_reco_beam_truth_Start_Px = 0.;
  alt_reco_beam_truth_Start_Py = 0.;
  alt_reco_beam_truth_Start_Pz = 0.;
  alt_reco_beam_truth_Start_E = 0.;
  alt_reco_beam_truth_Start_P = 0.;

  reco_beam_good = false;
  alt_reco_beam_good = false;
  passes_beam_cuts = false;
  alt_reco_beam_purity = 0.;

  flipped = false;
  elastic_candidate = false;
  processes.clear();
  found_elastic_scatter = false;
  nElasticScatters = 0;
  elastic_costheta.clear();
  elastic_X.clear();
  elastic_Y.clear();
  elastic_Z.clear();
  nReconstructed = 0;
  reconstructedIDs.clear();
  daughter_is_primary = false;
  reco_beam_Chi2_proton = 999.;
  traj_cos = 999.;
  view_0_hits_in_TPC5 = false;
  view_1_hits_in_TPC5 = false;
  view_2_hits_in_TPC5 = false;
  max_lateral = -999.;
  max_segment = -999.;

  view_0_max_segment = -999.;
  view_1_max_segment = -999.;
  view_2_max_segment = -999.;

  view_0_wire.clear(); 
  view_1_wire.clear(); 
  view_2_wire.clear(); 

  view_2_z.clear(); 

  view_0_tick.clear(); 
  view_1_tick.clear(); 
  view_2_tick.clear(); 

  view_0_wire_backtrack = 0.;
  view_1_wire_backtrack = 0.;
  view_2_wire_backtrack = 0.;

  reco_beam_Chi2_ndof = -1;

  reco_daughter_Chi2_proton.clear();
  reco_daughter_Chi2_ndof.clear();

  reco_daughter_shower_Chi2_proton.clear();
  reco_daughter_shower_Chi2_ndof.clear();

  reco_daughter_track_score.clear();
  reco_daughter_em_score.clear();
  reco_daughter_none_score.clear();
  reco_daughter_michel_score.clear();
  
  reco_daughter_shower_track_score.clear();
  reco_daughter_shower_em_score.clear();
  reco_daughter_shower_none_score.clear();
  reco_daughter_shower_michel_score.clear();

  reco_beam_truth_daughter_good_reco.clear();
  reco_beam_truth_daughter_true_PDGs.clear();
  reco_beam_truth_daughter_true_IDs.clear();

  reco_beam_truth_daughter_shower_good_reco.clear();
  reco_beam_truth_daughter_shower_true_PDGs.clear();
  reco_beam_truth_daughter_shower_true_IDs.clear();

  true_beam_daughter_PDGs.clear();
  true_beam_daughter_lens.clear();
  true_beam_daughter_startX.clear();
  true_beam_daughter_startY.clear();
  true_beam_daughter_startZ.clear();
  true_beam_daughter_startPX.clear();
  true_beam_daughter_startPY.clear();
  true_beam_daughter_startPZ.clear();
  true_beam_daughter_startP.clear();
  true_beam_daughter_Process.clear();

  true_beam_Pi0_decay_IDs.clear();
  true_beam_Pi0_decay_PDGs.clear();
  true_beam_grand_daughter_IDs.clear();
  true_beam_grand_daughter_ParIDs.clear();
  true_beam_grand_daughter_PDGs.clear();
  reco_beam_truth_daughter_true_lens.clear();
  reco_beam_truth_daughter_shower_true_lens.clear();
  true_beam_daughter_IDs.clear();
  reco_beam_truth_daughter_shower_good_reco.clear();

  nTrackDaughters = -1;
  nShowerDaughters = -1;

  reco_daughter_PFP_ID.clear();
  reco_daughter_PFP_track_score.clear();
  reco_daughter_PFP_em_score.clear();
  reco_daughter_PFP_michel_score.clear();
  reco_daughter_PFP_none_score.clear();

  dQdX.clear();
  dEdX.clear();
  vtxX = -1.;
  vtxY = -1.;
  vtxZ = -1.;
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

  resRange.clear();

  reco_daughter_dQdX.clear();
  reco_daughter_dEdX.clear();
  reco_daughter_resRange.clear();
  reco_daughter_shower_dQdX.clear();
  reco_daughter_shower_dEdX.clear();
  reco_daughter_shower_resRange.clear();
  reco_daughter_len.clear();

  beamTrackID = -1;
  reco_daughter_trackID.clear();
  reco_daughter_completeness.clear();
  reco_daughter_purity.clear();
  reco_daughter_truth_PDG.clear();
  reco_daughter_truth_ID.clear();
  reco_daughter_truth_Origin.clear();
  reco_daughter_truth_ParID.clear();
  reco_daughter_truth_Process.clear();

  alt_reco_daughter_truth_PDG.clear();
  alt_reco_daughter_truth_ID.clear();
  alt_reco_daughter_truth_Origin.clear();
  alt_reco_daughter_truth_ParID.clear();
  alt_reco_daughter_truth_Process.clear();
  alt_reco_daughter_truth_shared_hits.clear();
  alt_reco_daughter_truth_em_hits.clear();

  alt_reco_daughter_truth_len.clear();
  alt_reco_daughter_truth_startX.clear();
  alt_reco_daughter_truth_startY.clear();
  alt_reco_daughter_truth_startZ.clear();
  alt_reco_daughter_truth_endX.clear();
  alt_reco_daughter_truth_endY.clear();
  alt_reco_daughter_truth_endZ.clear();

  //Alternative Reco
  alt_reco_daughter_PFP_truth_PDG.clear();
  alt_reco_daughter_PFP_truth_ID.clear();
  alt_reco_daughter_PFP_truth_Origin.clear();
  alt_reco_daughter_PFP_truth_ParID.clear();
  alt_reco_daughter_PFP_truth_Process.clear();
  alt_reco_daughter_PFP_truth_shared_hits.clear();
  alt_reco_daughter_PFP_truth_em_hits.clear();

  alt_reco_daughter_PFP_truth_len.clear();
  alt_reco_daughter_PFP_truth_startX.clear();
  alt_reco_daughter_PFP_truth_startY.clear();
  alt_reco_daughter_PFP_truth_startZ.clear();
  alt_reco_daughter_PFP_truth_endX.clear();
  alt_reco_daughter_PFP_truth_endY.clear();
  alt_reco_daughter_PFP_truth_endZ.clear();

  reco_daughter_PMA_ID.clear();
  reco_daughter_PMA_dEdX.clear();
  reco_daughter_PMA_resRange.clear();
  reco_daughter_PMA_Chi2_proton.clear();
  reco_daughter_PMA_Chi2_ndof.clear();

  reco_daughter_pandora2_ID.clear();
  reco_daughter_pandora2_dEdX.clear();
  reco_daughter_pandora2_resRange.clear();
  reco_daughter_pandora2_Chi2_proton.clear();
  reco_daughter_pandora2_Chi2_ndof.clear();

  reco_daughter_PMA_startX.clear();
  reco_daughter_PMA_startY.clear();
  reco_daughter_PMA_startZ.clear();
  reco_daughter_PMA_endX.clear();
  reco_daughter_PMA_endY.clear();
  reco_daughter_PMA_endZ.clear();
  reco_daughter_PMA_dR.clear();
  reco_daughter_PMA_to_vertex.clear();

  reco_daughter_pandora2_startX.clear();
  reco_daughter_pandora2_startY.clear();
  reco_daughter_pandora2_startZ.clear();
  reco_daughter_pandora2_endX.clear();
  reco_daughter_pandora2_endY.clear();
  reco_daughter_pandora2_endZ.clear();
  reco_daughter_pandora2_dR.clear();
  reco_daughter_pandora2_to_vertex.clear();
  ///////
  

  alt_reco_daughter_shower_truth_PDG.clear();
  alt_reco_daughter_shower_truth_ID.clear();
  alt_reco_daughter_shower_truth_Origin.clear();
  alt_reco_daughter_shower_truth_ParID.clear();
  alt_reco_daughter_shower_truth_Process.clear();
  alt_reco_daughter_shower_truth_purity.clear();
  alt_reco_daughter_truth_purity.clear();

  reco_daughter_showerID.clear();
  reco_daughter_shower_truth_PDG.clear();
  reco_daughter_shower_truth_ID.clear();
  reco_daughter_shower_truth_Origin.clear();
  reco_daughter_shower_truth_ParID.clear();

}

DEFINE_ART_MODULE(pionana::PionAnalyzerMC)
