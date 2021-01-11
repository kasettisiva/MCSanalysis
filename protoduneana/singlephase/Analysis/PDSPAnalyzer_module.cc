////////////////////////////////////////////////////////////////////////
// Class:       PDSPAnalyzer
// Plugin Type: analyzer (art v3_00_00)
// File:        PDSPAnalyzer_module.cc
// Written by Jake Calcutt (calcuttj@msu.edu -- Slack: @jakecalcutt)
// Reach out for questions/issues/bugs
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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamCuts.h"
#include "protoduneana/Utilities/G4ReweightUtils.h"
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
#include "geant4reweight/src/ReweightBase/G4MultiReweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightManager.hh"

#include "cetlib/search_path.h"
#include "cetlib/filesystem.h"

#include "art_root_io/TFileService.h"
#include "TProfile.h"
#include "TFile.h"

// ROOT includes
#include "TTree.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TPolyLine3D.h"
#include "Math/Vector3D.h"
#include "TGraph2D.h"
#include "TROOT.h"

namespace pduneana {
  class PDSPAnalyzer;

  bool sort_IDEs( const sim::IDE * i1, const sim::IDE * i2){
    return( i1->z < i2->z );
  }

  double distance2(double x, double y, double z, double * p);
  void line(double t, double * p, double & x, double & y, double & z);
  void SumDistance2(int &, double *, double & sum, double * par, int);
  TVector3 FitLine(const std::vector<TVector3> & input);

  std::map<int, std::vector<const sim::IDE*>> slice_IDEs(
      std::vector<const sim::IDE*> ides, double the_z0, double the_pitch,
      double true_endZ){

    std::map<int, std::vector<const sim::IDE*>> results;

    for (size_t i = 0; i < ides.size(); ++i) {
      int slice_num = std::floor(
          (ides[i]->z - (the_z0 - the_pitch/2.)) / the_pitch);
      /*
      std::cout << "IDE: " << i << " ID: " << ides[i]->trackID << " Edep: "
                << ides[i]->energy << " (X,Y,Z): " << "(" << ides[i]->x << ","
                << ides[i]->y<<","<<ides[i]->z << ") Z0: " << the_z0
                << " Slice: " << slice_num << std::endl;*/
      results[slice_num].push_back(ides[i]);
    }

    return results;
  }

  TFile * OpenFile(const std::string filename) {
    TFile * theFile = 0x0;
    mf::LogInfo("pduneana::OpenFile") << "Searching for " << filename;
    if (cet::file_exists(filename)) {
      mf::LogInfo("pduneana::OpenFile") << "File exists. Opening " << filename;
      theFile = new TFile(filename.c_str());
      if (!theFile ||theFile->IsZombie() || !theFile->IsOpen()) {
        delete theFile;
        theFile = 0x0;
        throw cet::exception("PDSPAnalyzer_module.cc") << "Could not open " << filename;
      }
    }
    else {
      mf::LogInfo("pduneana::OpenFile") << "File does not exist here. Searching FW_SEARCH_PATH";
      cet::search_path sp{"FW_SEARCH_PATH"};
      std::string found_filename;
      auto found = sp.find_file(filename, found_filename);
      if (!found) {
        throw cet::exception("PDSPAnalyzer_module.cc") << "Could not find " << filename;
      }

      mf::LogInfo("pduneana::OpenFile") << "Found file " << found_filename;
      theFile = new TFile(found_filename.c_str());
      if (!theFile ||theFile->IsZombie() || !theFile->IsOpen()) {
        delete theFile;
        theFile = 0x0;
        throw cet::exception("PDSPAnalyzer_module.cc") << "Could not open " << found_filename;
      }
    }
    return theFile;
  }

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
    calo_point(size_t w, double p, double dqdx, double dedx, double dq,
               double cali_dedx, double r, size_t index, double input_z, int t)
        : wire(w), pitch(p), dQdX(dqdx), dEdX(dedx), calibrated_dEdX(cali_dedx),
          res_range(r), hit_index(index), z(input_z), tpc(t) {};

    size_t wire;
    double pitch;
    double dQdX;
    double dEdX;
    double dQ;
    double calibrated_dEdX;
    double res_range;
    size_t hit_index;
    double z;
    int tpc;
  };

  cnnOutput2D GetCNNOutputFromPFParticle(
      const recob::PFParticle & part, const art::Event & evt,
      const anab::MVAReader<recob::Hit,4> & CNN_results,
      protoana::ProtoDUNEPFParticleUtils & pfpUtil,
      std::string fPFParticleTag) {

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

using protoana::G4ReweightUtils::CreateRWTraj;
using protoana::G4ReweightUtils::CreateNRWTrajs;

pduneana::cnnOutput2D::cnnOutput2D() : track(0), em(0), michel(0), none(0), nHits(0) { }

class pduneana::PDSPAnalyzer : public art::EDAnalyzer {
public:
  explicit PDSPAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPAnalyzer(PDSPAnalyzer const&) = delete;
  PDSPAnalyzer(PDSPAnalyzer&&) = delete;
  PDSPAnalyzer& operator=(PDSPAnalyzer const&) = delete;
  PDSPAnalyzer& operator=(PDSPAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void reset();
  double lateralDist( TVector3 & n, TVector3 & x0, TVector3 & p );

private:

  void BeamTrackInfo(const art::Event & evt,
                     const recob::Track * thisTrack,
                     detinfo::DetectorClocksData const& clockData);
  void BeamShowerInfo(const recob::Shower* thisShower);
  void BeamPFPInfo(const art::Event & evt,
                   const recob::PFParticle* particle,
                   anab::MVAReader<recob::Hit,4> * hitResults);
  void BeamForcedTrackInfo(const art::Event & evt,
                           const recob::PFParticle * particle);
  void TrueBeamInfo(const art::Event & evt,
                    const simb::MCParticle* true_beam_particle,
                    detinfo::DetectorClocksData const& clockData,
                    const sim::ParticleList & plist,
                    std::map<int, std::vector<int>> & trueToPFPs,
                    anab::MVAReader<recob::Hit,4> * hitResults);
  void BeamInstInfo(const art::Event & evt);
  void DaughterPFPInfo(const art::Event & evt,
                       const recob::PFParticle* particle,
                       detinfo::DetectorClocksData const& clockData,
                       anab::MVAReader<recob::Hit,4> * hitResults);
  void DaughterMatchInfo(const art::Event & evt,
                         const recob::PFParticle * daughterPFP,
                         detinfo::DetectorClocksData const& clockData);
  void BeamMatchInfo(const art::Event & evt, const recob::PFParticle * particle,
                     const simb::MCParticle * true_beam_particle,
                     detinfo::DetectorClocksData const& clockData);
  void BeamForcedTrackInfo();
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
  double true_beam_mass;
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
  double true_beam_endP, true_beam_endP2, true_beam_last_len;

  int  true_beam_nElasticScatters;
  int  true_beam_nHits;
  std::vector< double > true_beam_elastic_costheta, true_beam_elastic_X,
                        true_beam_elastic_Y, true_beam_elastic_Z,
                        true_beam_elastic_deltaE, true_beam_elastic_IDE_edep;

  double true_beam_IDE_totalDep;
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
  std::vector< double > true_beam_Pi0_decay_startP, true_beam_Pi0_decay_startPx, true_beam_Pi0_decay_startPy, true_beam_Pi0_decay_startPz;
  std::vector< double > true_beam_Pi0_decay_startX, true_beam_Pi0_decay_startY, true_beam_Pi0_decay_startZ;
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
  ////////////////////////


  //Reconstructed track info
  //EDIT: STANDARDIZE
  double reco_beam_startX, reco_beam_startY, reco_beam_startZ;
  double reco_beam_endX, reco_beam_endY, reco_beam_endZ;
  double reco_beam_len, reco_beam_alt_len;
  double reco_beam_vertex_michel_score;
  int reco_beam_vertex_nHits;

  //position from SCE corrected calo
  double reco_beam_calo_startX, reco_beam_calo_startY, reco_beam_calo_startZ;
  double reco_beam_calo_endX, reco_beam_calo_endY, reco_beam_calo_endZ;
  std::vector<double> reco_beam_calo_startDirX, reco_beam_calo_endDirX;
  std::vector<double> reco_beam_calo_startDirY, reco_beam_calo_endDirY;
  std::vector<double> reco_beam_calo_startDirZ, reco_beam_calo_endDirZ;

  double reco_beam_trackDirX, reco_beam_trackDirY, reco_beam_trackDirZ;
  double reco_beam_trackEndDirX, reco_beam_trackEndDirY, reco_beam_trackEndDirZ;

  std::vector<double> reco_beam_dEdX_SCE, reco_beam_dQdX_SCE, reco_beam_resRange_SCE, reco_beam_TrkPitch_SCE;
  std::vector<double> reco_beam_calibrated_dEdX_SCE, reco_beam_dQ;

  std::vector<double> reco_beam_dEdX_NoSCE, reco_beam_dQdX_NoSCE, reco_beam_resRange_NoSCE, reco_beam_TrkPitch_NoSCE;
  std::vector<double> reco_beam_calibrated_dEdX_NoSCE;

  std::vector<double> reco_beam_calo_wire, reco_beam_calo_tick, reco_beam_calo_wire_z;
  std::vector<int> reco_beam_calo_TPC;

  int reco_beam_trackID;
  bool reco_beam_flipped;

  //fix
  bool reco_beam_passes_beam_cuts;              

  int reco_beam_type;
  double reco_beam_Chi2_proton;
  int    reco_beam_Chi2_ndof;

  ////////////////////////

  //For all track info 
  std::vector<double> reco_track_startX, reco_track_startY, reco_track_startZ,
                      reco_track_endX, reco_track_endY, reco_track_endZ,
                      reco_track_michel_score;
  std::vector<int> reco_track_ID, reco_track_nHits;


  //GeantReweight stuff
  // -- Maybe think of new naming scheme?
  std::vector<double> g4rw_primary_weights;
  std::vector<double> g4rw_primary_plus_sigma_weight;
  std::vector<double> g4rw_primary_minus_sigma_weight;
  std::vector<std::string> g4rw_primary_var;

  std::vector<double> g4rw_alt_primary_plus_sigma_weight;
  std::vector<double> g4rw_alt_primary_minus_sigma_weight;

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
  std::vector< int >    true_beam_slices, true_beam_slices_found;
  std::vector< double > true_beam_slices_deltaE;
  double true_beam_interactingEnergy;
  double em_energy;
  std::vector<double> true_beam_traj_X;
  std::vector<double> true_beam_traj_Y;
  std::vector<double> true_beam_traj_Z;
  std::vector<double> true_beam_traj_KE;

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
  //Info from the BI
  /////////////////////////////////////////////////////
  double beam_inst_P;
  bool beam_inst_valid;
  std::vector<double> beam_inst_TOF;
  std::vector< int > beam_inst_PDG_candidates, beam_inst_TOF_Chan;
  double beam_inst_X, beam_inst_Y, beam_inst_Z;
  double beam_inst_dirX, beam_inst_dirY, beam_inst_dirZ;
  int beam_inst_nFibersP1, beam_inst_nFibersP2, beam_inst_nFibersP3;
  int beam_inst_nTracks, beam_inst_nMomenta;
  ////////////////////////////////////////////////////


  //Alternative Reco values
  //EDIT: track_score --> trkScore, etc.
  std::vector< int > reco_daughter_PFP_ID;
  std::vector<int> reco_daughter_PFP_nHits,
                   reco_daughter_PFP_nHits_collection;
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
  std::vector< double > reco_daughter_PFP_true_byHits_completeness;

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

  std::vector< int > reco_daughter_PFP_true_byE_PDG;
  std::vector< double > reco_daughter_PFP_true_byE_len;
  std::vector< double > reco_daughter_PFP_true_byE_completeness;
  std::vector< double > reco_daughter_PFP_true_byE_purity;



  //////////////////////////////////////

  //EDIT: reco_daughter_allTrack_XXX
  std::vector< int > reco_daughter_allTrack_ID;
  std::vector< double > reco_daughter_allTrack_Theta;
  std::vector< double > reco_daughter_allTrack_Phi;
  std::vector< std::vector< double > > reco_daughter_allTrack_dQdX_SCE, reco_daughter_allTrack_dEdX_SCE, reco_daughter_allTrack_resRange_SCE;
  std::vector< std::vector< double > > reco_daughter_allTrack_calibrated_dEdX_SCE;
  std::vector< double > reco_daughter_allTrack_Chi2_proton;
  std::vector< int >    reco_daughter_allTrack_Chi2_ndof;

  //New: calorimetry + chi2 for planes 0 and 1
  std::vector<std::vector<double>>
      reco_daughter_allTrack_calibrated_dEdX_SCE_plane0,
      reco_daughter_allTrack_calibrated_dEdX_SCE_plane1;

  std::vector<std::vector<double>>
      reco_daughter_allTrack_resRange_plane0,
      reco_daughter_allTrack_resRange_plane1;

  std::vector<double> reco_daughter_allTrack_Chi2_proton_plane0,
                      reco_daughter_allTrack_Chi2_proton_plane1;

  std::vector<int> reco_daughter_allTrack_Chi2_ndof_plane0,
                   reco_daughter_allTrack_Chi2_ndof_plane1;
  //////////////////////////////////////////////

  std::vector< double > reco_daughter_allTrack_startX, reco_daughter_allTrack_endX;
  std::vector< double > reco_daughter_allTrack_startY, reco_daughter_allTrack_endY;
  std::vector< double > reco_daughter_allTrack_startZ, reco_daughter_allTrack_endZ;
  std::vector< double > reco_daughter_allTrack_dR;
  std::vector< double > reco_daughter_allTrack_len, reco_daughter_allTrack_alt_len;
  std::vector< double > reco_daughter_allTrack_to_vertex;

  std::vector<double> reco_daughter_allTrack_vertex_michel_score;
  std::vector<int> reco_daughter_allTrack_vertex_nHits;
  //

  std::vector<int>    reco_daughter_allShower_ID;
  std::vector<double> reco_daughter_allShower_len,
                      reco_daughter_allShower_startX,
                      reco_daughter_allShower_startY,
                      reco_daughter_allShower_startZ,
                      reco_daughter_allShower_dirX,
                      reco_daughter_allShower_dirY,
                      reco_daughter_allShower_dirZ,
                      reco_daughter_allShower_energy;


  std::vector<double> reco_daughter_allTrack_momByRange_proton;
  std::vector<double> reco_daughter_allTrack_momByRange_muon;
  double reco_beam_momByRange_proton;
  double reco_beam_momByRange_muon;

  std::vector<double> reco_daughter_allTrack_momByRange_alt_proton;
  std::vector<double> reco_daughter_allTrack_momByRange_alt_muon;
  double reco_beam_momByRange_alt_proton;
  double reco_beam_momByRange_alt_muon;

  //New hits info
  std::vector< double > reco_beam_spacePts_X, reco_beam_spacePts_Y, reco_beam_spacePts_Z;
  std::vector< std::vector< double > > reco_daughter_spacePts_X, reco_daughter_spacePts_Y, reco_daughter_spacePts_Z;
  std::vector< std::vector< double > > reco_daughter_shower_spacePts_X, reco_daughter_shower_spacePts_Y, reco_daughter_shower_spacePts_Z;




  ////New section -- mechanical class members
  std::map< int, TProfile* > templates;

  //FCL pars
  std::string fCalorimetryTagSCE;
  std::string fPandora2CaloSCE;
  std::string fCalorimetryTagNoSCE;
  std::string fTrackerTag;
  std::string fHitTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  std::string fBeamModuleLabel;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
  std::string dEdX_template_name;
  TFile * dEdX_template_file;
  bool fVerbose;    
  fhicl::ParameterSet BeamPars;
  fhicl::ParameterSet BeamCuts;
  protoana::ProtoDUNEBeamCuts beam_cuts;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  protoana::ProtoDUNETrackUtils trackUtil;
  protoana::ProtoDUNEShowerUtils showerUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  art::ServiceHandle<geo::Geometry> geom;
  protoana::ProtoDUNECalibration calibration_SCE;
  protoana::ProtoDUNECalibration calibration_NoSCE;
  bool fSaveHits;
  bool fSkipMVA;
  bool fTrueToReco;
  bool fDoReweight;
  bool fDoProtReweight;
  bool fGetTrackMichel;
  bool fMCHasBI;

  double fZ0, fPitch;

  //Geant4Reweight stuff
  TFile * FracsFile;
  TFile * ProtFracsFile;
  std::vector<fhicl::ParameterSet> ParSet;
  G4ReweightParameterMaker ParMaker;
  G4MultiReweighter * MultiRW, * ProtMultiRW;
  G4ReweightManager * RWManager;
};


pduneana::PDSPAnalyzer::PDSPAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  ,
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),
  fCalorimetryTagSCE(p.get<std::string>("CalorimetryTagSCE")),
  fPandora2CaloSCE(p.get<std::string>("Pandora2CaloSCE")),
  fCalorimetryTagNoSCE(p.get<std::string>("CalorimetryTagNoSCE")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fHitTag(p.get<std::string>("HitTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fBeamModuleLabel(p.get<std::string>("BeamModuleLabel")),
  fBeamlineUtils(p.get< fhicl::ParameterSet >("BeamlineUtils")),
  dEdX_template_name(p.get<std::string>("dEdX_template_name")),
  //dEdX_template_file( dEdX_template_name.c_str(), "OPEN" ),
  fVerbose(p.get<bool>("Verbose")),
  BeamPars(p.get<fhicl::ParameterSet>("BeamPars")),
  BeamCuts(p.get<fhicl::ParameterSet>("BeamCuts")),
  calibration_SCE(p.get<fhicl::ParameterSet>("CalibrationParsSCE")),
  calibration_NoSCE(p.get<fhicl::ParameterSet>("CalibrationParsNoSCE")),
  fSaveHits( p.get<bool>( "SaveHits" ) ),
  fSkipMVA( p.get<bool>( "SkipMVA" ) ),
  fTrueToReco( p.get<bool>( "TrueToReco" ) ),
  fDoReweight(p.get<bool>("DoReweight")),
  fDoProtReweight(p.get<bool>("DoProtReweight")),
  fGetTrackMichel(p.get<bool>("GetTrackMichel")) {

  dEdX_template_file = OpenFile(dEdX_template_name);
  templates[ 211 ]  = (TProfile*)dEdX_template_file->Get( "dedx_range_pi"  );
  templates[ 321 ]  = (TProfile*)dEdX_template_file->Get( "dedx_range_ka"  );
  templates[ 13 ]   = (TProfile*)dEdX_template_file->Get( "dedx_range_mu"  );
  templates[ 2212 ] = (TProfile*)dEdX_template_file->Get( "dedx_range_pro" );

  beam_cuts = protoana::ProtoDUNEBeamCuts( BeamCuts );

  if (fDoReweight || fDoProtReweight) {
    RWManager = new G4ReweightManager({p.get<fhicl::ParameterSet>("Material")});
  }

  if (fDoReweight) {
    FracsFile =  new TFile((p.get< std::string >( "FracsFile" )).c_str(), "OPEN" );
    ParSet = p.get<std::vector<fhicl::ParameterSet>>("ParameterSet");
    ParMaker = G4ReweightParameterMaker(ParSet);
    MultiRW = new G4MultiReweighter(211, *FracsFile, ParSet,
                                    p.get<fhicl::ParameterSet>("Material"),
                                    RWManager);
  }
  if (fDoProtReweight) {
    ProtFracsFile =  new TFile((p.get<std::string>("ProtFracsFile")).c_str(),
                               "OPEN");
    ParSet = p.get<std::vector<fhicl::ParameterSet>>("ParameterSet");
    ParMaker = G4ReweightParameterMaker(ParSet);
    ProtMultiRW = new G4MultiReweighter(2212, *ProtFracsFile, ParSet,
                                    p.get<fhicl::ParameterSet>("Material"),
                                    RWManager);
  }

  fZ0 = geom->Wire( geo::WireID(0, 1, 2, 0) ).GetCenter().Z();
  fPitch = geom->WirePitch( 2, 1, 0);

}

void pduneana::PDSPAnalyzer::analyze(art::Event const & evt) {

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

  if (fGetTrackMichel && !fSkipMVA) {
    for (const recob::PFParticle & pfp : (*pfpVec)) {
      const recob::Track* tempTrack = pfpUtil.GetPFParticleTrack(pfp, evt,
                                                                 fPFParticleTag,
                                                                 fTrackerTag);
      if (tempTrack) {
        double startX = tempTrack->Start().X();
        double startY = tempTrack->Start().Y();
        double startZ = tempTrack->Start().Z();

        double endX = tempTrack->End().X();
        double endY = tempTrack->End().Y();
        double endZ = tempTrack->End().Z();

        double start[3] = {startX, startY, startZ};
        double end[3] = {endX, endY, endZ};
        int end_tpc = geom->FindTPCAtPosition(end).TPC;
        int start_tpc = geom->FindTPCAtPosition(start).TPC;

        if (!((end_tpc == 1 || end_tpc == 5) &&
              (start_tpc == 1 || start_tpc == 5)))
          continue;

        std::pair<double, int> vertex_michel_score =
            trackUtil.GetVertexMichelScore(*tempTrack, evt, fTrackerTag,
                                           fHitTag);

        reco_track_michel_score.push_back(vertex_michel_score.first);
        reco_track_nHits.push_back(vertex_michel_score.second);
        reco_track_ID.push_back(tempTrack->ID());
        reco_track_startX.push_back(startX);
        reco_track_startY.push_back(startY);
        reco_track_startZ.push_back(startZ);
        reco_track_endX.push_back(endX);
        reco_track_endY.push_back(endY);
        reco_track_endZ.push_back(endZ);
      }
    }
  }

  const sim::ParticleList & plist = pi_serv->ParticleList();

  art::ServiceHandle < geo::Geometry > fGeometryService;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp =  art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
  trkf::TrackMomentumCalculator track_p_calc;
  ////////////////////////////////////////


  //double z0 = geom->Wire( geo::WireID(0, 1, 2, 0) ).GetCenter().Z();
  //double pitch = geom->WirePitch( 2, 1, 0);

  if (fVerbose) {
    std::cout << "Z0: " << fZ0 << std::endl;
    std::cout << "Pitch: " << fPitch << std::endl;

    double z0_APA2 = geom->Wire(geo::WireID(0, 5, 2, 0)).GetCenter().Z();
    std::cout << "APA 2 Z0: " << z0_APA2 << std::endl;
  }

  // This gets the true beam particle that generated the event
  const simb::MCParticle* true_beam_particle = 0x0;
  if( !evt.isRealData() ){
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    true_beam_particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
    if( !true_beam_particle ){
      MF_LOG_WARNING("PDSPAnalyzer") << "No true beam particle" << std::endl;
      return;
    }
    if (fVerbose) {
      std::cout << "Got " << (*mcTruths)[0].NParticles() <<
                   " particles in mcTruth" << std::endl;
      for (int i = 0; i < (*mcTruths)[0].NParticles(); ++i) {
        simb::MCParticle part = (*mcTruths)[0].GetParticle(i);
        std::cout << part.Process() << " " << part.TrackId() << " " <<
                     part.PdgCode() << std::endl;

      }
    }
  }
  ////////////////////////////

  //Get the beam instrumentation from the event
  BeamInstInfo(evt);


  // Helper to get hits and the 4 associated CNN outputs
  // CNN Outputs: EM, Track, Michel, Empty
  // outputNames: track, em, none, michel
  anab::MVAReader<recob::Hit,4> * hitResults = 0x0;
  if (!fSkipMVA) { 
    hitResults = new anab::MVAReader<recob::Hit, 4>(evt, "emtrkmichelid:emtrkmichel" );
  }

  auto allHits = evt.getValidHandle<std::vector<recob::Hit> >(fHitTag);

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fTrackerTag);
  art::FindManyP<recob::Hit> findHits(recoTracks,evt,fTrackerTag);

  auto recoShowers = evt.getValidHandle< std::vector< recob::Shower > >(fShowerTag);
  art::FindManyP<recob::Hit> findHitsFromShowers(recoShowers,evt,fShowerTag);

  std::map< int, std::vector< int > > trueToPFPs;
  if( fTrueToReco ){
    trueToPFPs = truthUtil.GetMapMCToPFPs_ByHits( clockData, evt, fPFParticleTag, fHitTag );
  }


  ///Gets the beam pfparticle
  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);

  if(beamParticles.size() == 0){
    std::cout << "We found no beam particles for this event... moving on" << std::endl;
    return;
  }
  else {
    std::cout << "Found " << beamParticles.size() << " particles" << std::endl;
  }

  // Get the reconstructed PFParticle tagged as beam by Pandora
  const recob::PFParticle* particle = beamParticles.at(0);
  //////////////////////////////////////////////////////////////////
  
  
  //If MC, attempt to match to some MCParticle
  if( !evt.isRealData() ){
    BeamMatchInfo(evt, particle, true_beam_particle, clockData);
    TrueBeamInfo(evt, true_beam_particle, clockData, plist, trueToPFPs, hitResults);
  }

  BeamPFPInfo(evt, particle, hitResults);

  // Determine if the beam particle is track-like or shower-like
  const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
  if( thisTrack ){
    BeamTrackInfo(evt, thisTrack, clockData);
  }
  else if( thisShower ){
    BeamShowerInfo(thisShower);
  }

  DaughterPFPInfo(evt, particle, clockData, hitResults);

  BeamForcedTrackInfo(evt, particle);

  //New geant4reweight stuff
  if (!evt.isRealData() && fDoReweight) {
    if (fVerbose) std::cout << "Doing reweight" << std::endl;
    if (true_beam_PDG == 211) {
      G4ReweightTraj theTraj(true_beam_ID, true_beam_PDG, 0, event, {0,0});
      bool created = CreateRWTraj(*true_beam_particle, plist,
                                  fGeometryService, event, &theTraj);
      if (created && theTraj.GetNSteps() > 0) {
        g4rw_primary_weights.push_back(MultiRW->GetWeightFromNominal(theTraj));

        std::cout << g4rw_primary_weights.size() << std::endl;
        std::vector<double> weights_vec = MultiRW->GetWeightFromAll1DThrows(
            theTraj);
        std::cout << weights_vec.size() << std::endl;
        g4rw_primary_weights.insert(g4rw_primary_weights.end(),
                                    weights_vec.begin(), weights_vec.end());


        //g4rw_primary_plus_sigma_weight = pm_weights.first;
        //g4rw_primary_minus_sigma_weight = pm_weights.second;

        for (size_t i = 0; i < ParSet.size(); ++i) {
          std::pair<double, double> pm_weights =
              MultiRW->GetPlusMinusSigmaParWeight(theTraj, i);

          g4rw_primary_plus_sigma_weight.push_back(pm_weights.first);
          g4rw_primary_minus_sigma_weight.push_back(pm_weights.second);
          g4rw_primary_var.push_back(ParSet[i].get<std::string>("Name"));
        }

      }

      std::vector<G4ReweightTraj *> trajs = CreateNRWTrajs(
          *true_beam_particle, plist,
          fGeometryService, event);
      
      bool added = false;
      for (size_t i = 0; i < trajs.size(); ++i) {
        if (trajs[i]->GetNSteps() > 0) {
          //std::cout << i << " " << trajs[i]->GetNSteps() << std::endl;
          for (size_t j = 0; j < ParSet.size(); ++j) {
            std::pair<double, double> pm_weights =
                MultiRW->GetPlusMinusSigmaParWeight((*trajs[i]), j);
            //std::cout << "got weights" << std::endl;
            //std::cout << pm_weights.first << " " << pm_weights.second << std::endl;

            if (!added) {
              g4rw_alt_primary_plus_sigma_weight.push_back(pm_weights.first);
              g4rw_alt_primary_minus_sigma_weight.push_back(pm_weights.second);
            }
            else {
              g4rw_alt_primary_plus_sigma_weight[j] *= pm_weights.first;
              g4rw_alt_primary_minus_sigma_weight[j] *= pm_weights.second;
            }
          }
          added = true;
        }
      }
      
    }
  }
  if (!evt.isRealData() && fDoProtReweight && true_beam_PDG == 2212) {
    std::vector<G4ReweightTraj *> trajs = CreateNRWTrajs(
        *true_beam_particle, plist,
        fGeometryService, event);
    
    bool added = false;
    for (size_t i = 0; i < trajs.size(); ++i) {
      if (trajs[i]->GetNSteps() > 0) {
        //std::cout << i << " " << trajs[i]->GetNSteps() << std::endl;
        for (size_t j = 0; j < ParSet.size(); ++j) {
          std::pair<double, double> pm_weights =
              ProtMultiRW->GetPlusMinusSigmaParWeight((*trajs[i]), j);
          //std::cout << "alt got weights" << std::endl;
          //std::cout << pm_weights.first << " " << pm_weights.second <<
          //             std::endl;

          if (!added) {
            g4rw_alt_primary_plus_sigma_weight.push_back(pm_weights.first);
            g4rw_alt_primary_minus_sigma_weight.push_back(pm_weights.second);
          }
          else {
            g4rw_alt_primary_plus_sigma_weight[j] *= pm_weights.first;
            g4rw_alt_primary_minus_sigma_weight[j] *= pm_weights.second;
          }
        }
        added = true;
      }
    }
  }

  fTree->Fill();
}

void pduneana::PDSPAnalyzer::beginJob()
{

  gROOT->SetBatch(1);

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
  fTree->Branch("reco_beam_alt_len", &reco_beam_alt_len);
  fTree->Branch("reco_beam_calo_startX", &reco_beam_calo_startX);
  fTree->Branch("reco_beam_calo_startY", &reco_beam_calo_startY);
  fTree->Branch("reco_beam_calo_startZ", &reco_beam_calo_startZ);
  fTree->Branch("reco_beam_calo_endX", &reco_beam_calo_endX);
  fTree->Branch("reco_beam_calo_endY", &reco_beam_calo_endY);
  fTree->Branch("reco_beam_calo_endZ", &reco_beam_calo_endZ);
  fTree->Branch("reco_beam_calo_startDirX", &reco_beam_calo_startDirX);
  fTree->Branch("reco_beam_calo_startDirY", &reco_beam_calo_startDirY);
  fTree->Branch("reco_beam_calo_startDirZ", &reco_beam_calo_startDirZ);
  fTree->Branch("reco_beam_calo_endDirX", &reco_beam_calo_endDirX);
  fTree->Branch("reco_beam_calo_endDirY", &reco_beam_calo_endDirY);
  fTree->Branch("reco_beam_calo_endDirZ", &reco_beam_calo_endDirZ);

  fTree->Branch("reco_beam_trackDirX", &reco_beam_trackDirX);
  fTree->Branch("reco_beam_trackDirY", &reco_beam_trackDirY);
  fTree->Branch("reco_beam_trackDirZ", &reco_beam_trackDirZ);
  fTree->Branch("reco_beam_trackEndDirX", &reco_beam_trackEndDirX);
  fTree->Branch("reco_beam_trackEndDirY", &reco_beam_trackEndDirY);
  fTree->Branch("reco_beam_trackEndDirZ", &reco_beam_trackEndDirZ);
  fTree->Branch("reco_beam_vertex_nHits", &reco_beam_vertex_nHits);
  fTree->Branch("reco_beam_vertex_michel_score", &reco_beam_vertex_michel_score);
  fTree->Branch("reco_beam_trackID", &reco_beam_trackID);

  fTree->Branch("reco_beam_dQdX_SCE", &reco_beam_dQdX_SCE);
  fTree->Branch("reco_beam_dQ", &reco_beam_dQ);
  fTree->Branch("reco_beam_dEdX_SCE", &reco_beam_dEdX_SCE);
  fTree->Branch("reco_beam_calibrated_dEdX_SCE", &reco_beam_calibrated_dEdX_SCE);
  fTree->Branch("reco_beam_resRange_SCE", &reco_beam_resRange_SCE);
  fTree->Branch("reco_beam_TrkPitch_SCE", &reco_beam_TrkPitch_SCE);

  fTree->Branch("reco_beam_dQdX_NoSCE", &reco_beam_dQdX_NoSCE);
  fTree->Branch("reco_beam_dEdX_NoSCE", &reco_beam_dEdX_NoSCE);
  fTree->Branch("reco_beam_calibrated_dEdX_NoSCE", &reco_beam_calibrated_dEdX_NoSCE);
  fTree->Branch("reco_beam_resRange_NoSCE", &reco_beam_resRange_NoSCE);
  fTree->Branch("reco_beam_TrkPitch_NoSCE", &reco_beam_TrkPitch_NoSCE);

  fTree->Branch("reco_beam_calo_wire", &reco_beam_calo_wire);
  fTree->Branch("reco_beam_calo_wire_z", &reco_beam_calo_wire_z);
  fTree->Branch("reco_beam_calo_tick", &reco_beam_calo_tick);
  fTree->Branch("reco_beam_calo_TPC", &reco_beam_calo_TPC);

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

  fTree->Branch("reco_track_startX", &reco_track_startX);
  fTree->Branch("reco_track_startY", &reco_track_startY);
  fTree->Branch("reco_track_startZ", &reco_track_startZ);
  fTree->Branch("reco_track_endX", &reco_track_endX);
  fTree->Branch("reco_track_endY", &reco_track_endY);
  fTree->Branch("reco_track_endZ", &reco_track_endZ);
  fTree->Branch("reco_track_michel_score", &reco_track_michel_score);
  fTree->Branch("reco_track_ID", &reco_track_ID);
  fTree->Branch("reco_track_nHits", &reco_track_nHits);

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
  fTree->Branch("reco_daughter_PFP_true_byHits_completeness", &reco_daughter_PFP_true_byHits_completeness);
  fTree->Branch("reco_daughter_PFP_true_byE_PDG", &reco_daughter_PFP_true_byE_PDG);
  fTree->Branch("reco_daughter_PFP_true_byE_len", &reco_daughter_PFP_true_byE_len);
  fTree->Branch("reco_daughter_PFP_true_byE_completeness", &reco_daughter_PFP_true_byE_completeness);
  fTree->Branch("reco_daughter_PFP_true_byE_purity", &reco_daughter_PFP_true_byE_purity);

  fTree->Branch("reco_daughter_allTrack_ID", &reco_daughter_allTrack_ID);
  fTree->Branch("reco_daughter_allTrack_dQdX_SCE", &reco_daughter_allTrack_dQdX_SCE);
  fTree->Branch("reco_daughter_allTrack_dEdX_SCE", &reco_daughter_allTrack_dEdX_SCE);
  fTree->Branch("reco_daughter_allTrack_resRange_SCE", &reco_daughter_allTrack_resRange_SCE);
  fTree->Branch("reco_daughter_allTrack_calibrated_dEdX_SCE", &reco_daughter_allTrack_calibrated_dEdX_SCE);

  fTree->Branch("reco_daughter_allTrack_Chi2_proton", &reco_daughter_allTrack_Chi2_proton);
  fTree->Branch("reco_daughter_allTrack_Chi2_ndof", &reco_daughter_allTrack_Chi2_ndof);


  ///Calorimetry/chi2 planes 0 and 1
  fTree->Branch("reco_daughter_allTrack_Chi2_proton_plane0",
                &reco_daughter_allTrack_Chi2_proton_plane0);
  fTree->Branch("reco_daughter_allTrack_Chi2_proton_plane1",
                &reco_daughter_allTrack_Chi2_proton_plane1);

  fTree->Branch("reco_daughter_allTrack_Chi2_ndof_plane0",
                &reco_daughter_allTrack_Chi2_ndof_plane0);
  fTree->Branch("reco_daughter_allTrack_Chi2_ndof_plane1",
                &reco_daughter_allTrack_Chi2_ndof_plane1);

  fTree->Branch("reco_daughter_allTrack_calibrated_dEdX_SCE_plane0",
                &reco_daughter_allTrack_calibrated_dEdX_SCE_plane0);
  fTree->Branch("reco_daughter_allTrack_calibrated_dEdX_SCE_plane1",
                &reco_daughter_allTrack_calibrated_dEdX_SCE_plane1);
  fTree->Branch("reco_daughter_allTrack_resRange_plane0",
                &reco_daughter_allTrack_resRange_plane0);
  fTree->Branch("reco_daughter_allTrack_resRange_plane1",
                &reco_daughter_allTrack_resRange_plane1);
  ///////////////////////////////////

  fTree->Branch("reco_daughter_allTrack_Theta", &reco_daughter_allTrack_Theta);
  fTree->Branch("reco_daughter_allTrack_Phi", &reco_daughter_allTrack_Phi);

  fTree->Branch("reco_daughter_allTrack_len", &reco_daughter_allTrack_len);
  fTree->Branch("reco_daughter_allTrack_alt_len", &reco_daughter_allTrack_alt_len);
  fTree->Branch("reco_daughter_allTrack_startX", &reco_daughter_allTrack_startX);
  fTree->Branch("reco_daughter_allTrack_startY", &reco_daughter_allTrack_startY);
  fTree->Branch("reco_daughter_allTrack_startZ", &reco_daughter_allTrack_startZ);
  fTree->Branch("reco_daughter_allTrack_endX", &reco_daughter_allTrack_endX);
  fTree->Branch("reco_daughter_allTrack_endY", &reco_daughter_allTrack_endY);
  fTree->Branch("reco_daughter_allTrack_endZ", &reco_daughter_allTrack_endZ);
  fTree->Branch("reco_daughter_allTrack_dR", &reco_daughter_allTrack_dR);
  fTree->Branch("reco_daughter_allTrack_to_vertex", &reco_daughter_allTrack_to_vertex);

  fTree->Branch("reco_daughter_allTrack_vertex_michel_score",
                &reco_daughter_allTrack_vertex_michel_score);
  fTree->Branch("reco_daughter_allTrack_vertex_nHits",
                &reco_daughter_allTrack_vertex_nHits);
  //////

  fTree->Branch("reco_daughter_allShower_ID", &reco_daughter_allShower_ID);
  fTree->Branch("reco_daughter_allShower_len", &reco_daughter_allShower_len);
  fTree->Branch("reco_daughter_allShower_startX", &reco_daughter_allShower_startX);
  fTree->Branch("reco_daughter_allShower_startY", &reco_daughter_allShower_startY);
  fTree->Branch("reco_daughter_allShower_startZ", &reco_daughter_allShower_startZ);
  fTree->Branch("reco_daughter_allShower_dirX", &reco_daughter_allShower_dirX);
  fTree->Branch("reco_daughter_allShower_dirY", &reco_daughter_allShower_dirY);
  fTree->Branch("reco_daughter_allShower_dirZ", &reco_daughter_allShower_dirZ);
  fTree->Branch("reco_daughter_allShower_energy", &reco_daughter_allShower_energy);



  fTree->Branch("reco_daughter_PFP_ID", &reco_daughter_PFP_ID);
  fTree->Branch("reco_daughter_PFP_nHits", &reco_daughter_PFP_nHits);
  fTree->Branch("reco_daughter_PFP_nHits_collection",
                &reco_daughter_PFP_nHits_collection);
  fTree->Branch("reco_daughter_PFP_trackScore", &reco_daughter_PFP_trackScore);
  fTree->Branch("reco_daughter_PFP_emScore", &reco_daughter_PFP_emScore);
  fTree->Branch("reco_daughter_PFP_michelScore", &reco_daughter_PFP_michelScore);
  fTree->Branch("reco_daughter_PFP_trackScore_collection", &reco_daughter_PFP_trackScore_collection);
  fTree->Branch("reco_daughter_PFP_emScore_collection", &reco_daughter_PFP_emScore_collection);
  fTree->Branch("reco_daughter_PFP_michelScore_collection", &reco_daughter_PFP_michelScore_collection);




  fTree->Branch("true_beam_PDG", &true_beam_PDG);
  fTree->Branch("true_beam_mass", &true_beam_mass);
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
  fTree->Branch("true_beam_endP2", &true_beam_endP2);
  fTree->Branch("true_beam_last_len", &true_beam_last_len);

  fTree->Branch("true_beam_startDirX", &true_beam_startDirX);
  fTree->Branch("true_beam_startDirY", &true_beam_startDirY);
  fTree->Branch("true_beam_startDirZ", &true_beam_startDirZ);

  fTree->Branch("true_beam_nElasticScatters", &true_beam_nElasticScatters);
  fTree->Branch("true_beam_elastic_costheta", &true_beam_elastic_costheta);
  fTree->Branch("true_beam_elastic_X", &true_beam_elastic_X);
  fTree->Branch("true_beam_elastic_Y", &true_beam_elastic_Y);
  fTree->Branch("true_beam_elastic_Z", &true_beam_elastic_Z);
  fTree->Branch("true_beam_elastic_deltaE", &true_beam_elastic_deltaE);
  fTree->Branch("true_beam_elastic_IDE_edep", &true_beam_elastic_IDE_edep);
  fTree->Branch("true_beam_IDE_totalDep",    &true_beam_IDE_totalDep);

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
  fTree->Branch("true_beam_Pi0_decay_startPx", &true_beam_Pi0_decay_startPx);
  fTree->Branch("true_beam_Pi0_decay_startPy", &true_beam_Pi0_decay_startPy);
  fTree->Branch("true_beam_Pi0_decay_startPz", &true_beam_Pi0_decay_startPz);
  fTree->Branch("true_beam_Pi0_decay_startX", &true_beam_Pi0_decay_startX);
  fTree->Branch("true_beam_Pi0_decay_startY", &true_beam_Pi0_decay_startY);
  fTree->Branch("true_beam_Pi0_decay_startZ", &true_beam_Pi0_decay_startZ);

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

  fTree->Branch("beam_inst_P", &beam_inst_P);
  fTree->Branch("beam_inst_TOF", &beam_inst_TOF);
  fTree->Branch("beam_inst_TOF_Chan", &beam_inst_TOF_Chan);
  fTree->Branch("beam_inst_X", &beam_inst_X);
  fTree->Branch("beam_inst_Y", &beam_inst_Y);
  fTree->Branch("beam_inst_Z", &beam_inst_Z);
  fTree->Branch("beam_inst_dirX", &beam_inst_dirX);
  fTree->Branch("beam_inst_dirY", &beam_inst_dirY);
  fTree->Branch("beam_inst_dirZ", &beam_inst_dirZ);

  fTree->Branch("beam_inst_nFibersP1", &beam_inst_nFibersP1);
  fTree->Branch("beam_inst_nFibersP2", &beam_inst_nFibersP2);
  fTree->Branch("beam_inst_nFibersP3", &beam_inst_nFibersP3);
  fTree->Branch("beam_inst_PDG_candidates", &beam_inst_PDG_candidates);
  fTree->Branch("beam_inst_nTracks", &beam_inst_nTracks);
  fTree->Branch("beam_inst_nMomenta", &beam_inst_nMomenta);
  fTree->Branch("beam_inst_valid", &beam_inst_valid);


  fTree->Branch("reco_beam_Chi2_proton", &reco_beam_Chi2_proton);
  fTree->Branch("reco_beam_Chi2_ndof", &reco_beam_Chi2_ndof);


  fTree->Branch("reco_daughter_allTrack_momByRange_proton", &reco_daughter_allTrack_momByRange_proton);
  fTree->Branch("reco_daughter_allTrack_momByRange_muon", &reco_daughter_allTrack_momByRange_muon);
  fTree->Branch("reco_beam_momByRange_proton", &reco_beam_momByRange_proton);
  fTree->Branch("reco_beam_momByRange_muon", &reco_beam_momByRange_muon);

  fTree->Branch("reco_daughter_allTrack_momByRange_alt_proton", &reco_daughter_allTrack_momByRange_alt_proton);
  fTree->Branch("reco_daughter_allTrack_momByRange_alt_muon", &reco_daughter_allTrack_momByRange_alt_muon);
  fTree->Branch("reco_beam_momByRange_alt_proton", &reco_beam_momByRange_alt_proton);
  fTree->Branch("reco_beam_momByRange_alt_muon", &reco_beam_momByRange_alt_muon);

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
  fTree->Branch("true_beam_slices", &true_beam_slices);
  fTree->Branch("true_beam_slices_found", &true_beam_slices_found);
  fTree->Branch("true_beam_slices_deltaE", &true_beam_slices_deltaE);
  fTree->Branch("em_energy", &em_energy);
  fTree->Branch("true_beam_traj_X", &true_beam_traj_X);
  fTree->Branch("true_beam_traj_Y", &true_beam_traj_Y);
  fTree->Branch("true_beam_traj_Z", &true_beam_traj_Z);
  fTree->Branch("true_beam_traj_KE", &true_beam_traj_KE);

  fTree->Branch("g4rw_primary_weights", &g4rw_primary_weights);
  fTree->Branch("g4rw_primary_plus_sigma_weight", &g4rw_primary_plus_sigma_weight);
  fTree->Branch("g4rw_primary_minus_sigma_weight", &g4rw_primary_minus_sigma_weight);
  fTree->Branch("g4rw_primary_var", &g4rw_primary_var);

  fTree->Branch("g4rw_alt_primary_plus_sigma_weight",
                &g4rw_alt_primary_plus_sigma_weight);
  fTree->Branch("g4rw_alt_primary_minus_sigma_weight",
                &g4rw_alt_primary_minus_sigma_weight);

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

void pduneana::PDSPAnalyzer::endJob()
{
  dEdX_template_file->Close();
}

double pduneana::PDSPAnalyzer::lateralDist(TVector3 &n, TVector3 &x0, TVector3 &p){
  TVector3 x = ( (p - x0)*n )*n;
  return (x - (p - x0)).Mag();
}

void pduneana::PDSPAnalyzer::reset()
{
  reco_beam_startX = -999;
  reco_beam_startY = -999;
  reco_beam_startZ = -999;
  reco_beam_endX = -999;
  reco_beam_endY = -999;
  reco_beam_endZ = -999;
  reco_beam_flipped = false;
  reco_beam_trackEndDirX = -999;
  reco_beam_trackEndDirY = -999;
  reco_beam_trackEndDirZ = -999;
  reco_beam_trackDirX = -999;
  reco_beam_trackDirY = -999;
  reco_beam_trackDirZ = -999;

  reco_beam_len = -999;
  reco_beam_alt_len = -999;
  reco_beam_calo_startX = -999;
  reco_beam_calo_startY = -999;
  reco_beam_calo_startZ = -999;
  reco_beam_calo_endX = -999;
  reco_beam_calo_endY = -999;
  reco_beam_calo_endZ = -999;
  reco_beam_calo_startDirX.clear();
  reco_beam_calo_startDirY.clear();
  reco_beam_calo_startDirZ.clear();
  reco_beam_calo_endDirX.clear();
  reco_beam_calo_endDirY.clear();
  reco_beam_calo_endDirZ.clear();

  reco_track_startX.clear();
  reco_track_startY.clear();
  reco_track_startZ.clear();
  reco_track_endX.clear();
  reco_track_endY.clear();
  reco_track_endZ.clear();
  reco_track_michel_score.clear();
  reco_track_ID.clear();
  reco_track_nHits.clear();

  reco_beam_type = -999;
  reco_beam_passes_beam_cuts = false;

  reco_beam_vertex_slice = std::numeric_limits<int>::max();

  true_daughter_nPi0 = 0;
  true_daughter_nPiPlus = 0;
  true_daughter_nPiMinus = 0;
  true_daughter_nProton = 0;
  true_daughter_nNeutron = 0;
  true_daughter_nNucleus = 0;

  reco_beam_true_byE_PDG = -999;
  reco_beam_true_byE_ID = -999;
  reco_beam_true_byHits_PDG = -999;
  reco_beam_true_byHits_ID = -999;

  true_beam_PDG = -999;
  true_beam_mass = -999.;
  true_beam_ID = -999;
  true_beam_endProcess ="";
  true_beam_endX = -999.;
  true_beam_endY = -999.;
  true_beam_endZ = -999.;
  true_beam_startX = -999.;
  true_beam_startY = -999.;
  true_beam_startZ = -999.;

  true_beam_startPx   = -999.;
  true_beam_startPy   = -999.;
  true_beam_startPz   = -999.;
  true_beam_startP    = -999.;

  true_beam_endPx   = -999.;
  true_beam_endPy   = -999.;
  true_beam_endPz   = -999.;
  true_beam_endP    = -999.;
  true_beam_endP2    = -999.;
  true_beam_last_len    = -999.;

  true_beam_startDirX = -999.;
  true_beam_startDirY = -999.;
  true_beam_startDirZ = -999.;
  true_beam_nHits = -999;


  true_beam_processes.clear();
  true_beam_nElasticScatters = 0;
  true_beam_elastic_costheta.clear();
  true_beam_elastic_X.clear();
  true_beam_elastic_Y.clear();
  true_beam_elastic_Z.clear();
  true_beam_elastic_deltaE.clear();
  true_beam_elastic_IDE_edep.clear();
  true_beam_IDE_totalDep = -999.;

  true_beam_reco_byHits_PFP_ID.clear();
  true_beam_reco_byHits_PFP_nHits.clear();
  true_beam_reco_byHits_allTrack_ID.clear();

  reco_beam_true_byE_endProcess ="";
  reco_beam_true_byE_process ="";
  reco_beam_true_byE_origin = -999;

  reco_beam_true_byE_endPx = -999.;
  reco_beam_true_byE_endPy = -999.;
  reco_beam_true_byE_endPz = -999.;
  reco_beam_true_byE_endE = -999.;
  reco_beam_true_byE_endP = -999.;

  reco_beam_true_byE_startPx = -999.;
  reco_beam_true_byE_startPy = -999.;
  reco_beam_true_byE_startPz = -999.;
  reco_beam_true_byE_startE = -999.;
  reco_beam_true_byE_startP = -999.;

  reco_beam_true_byHits_endProcess ="";
  reco_beam_true_byHits_process ="";
  reco_beam_true_byHits_origin = -999;

  reco_beam_true_byHits_endPx = -999.;
  reco_beam_true_byHits_endPy = -999.;
  reco_beam_true_byHits_endPz = -999.;
  reco_beam_true_byHits_endE = -999.;
  reco_beam_true_byHits_endP = -999.;

  reco_beam_true_byHits_startPx = -999.;
  reco_beam_true_byHits_startPy = -999.;
  reco_beam_true_byHits_startPz = -999.;
  reco_beam_true_byHits_startE = -999.;
  reco_beam_true_byHits_startP = -999.;

  reco_beam_true_byE_matched = false;
  reco_beam_true_byHits_matched = false;
  reco_beam_true_byHits_purity = -999.;


  //reco_daughter_true_byE_isPrimary = false;
  reco_beam_Chi2_proton = 999.;



  beam_inst_P = -999.;
  beam_inst_X = -999.;
  beam_inst_Y = -999.;
  beam_inst_Z = -999.;
  beam_inst_dirX = -999.;
  beam_inst_dirY = -999.;
  beam_inst_dirZ = -999.;
  beam_inst_nFibersP1 = -999;
  beam_inst_nFibersP2 = -999;
  beam_inst_nFibersP3 = -999;
  beam_inst_PDG_candidates.clear();
  beam_inst_TOF.clear();
  beam_inst_TOF_Chan.clear();
  beam_inst_nTracks = -999;
  beam_inst_nMomenta = -999;
  beam_inst_valid = true;

  reco_beam_Chi2_ndof = -999;

  reco_daughter_allTrack_momByRange_proton.clear();
  reco_daughter_allTrack_momByRange_muon.clear();
  reco_beam_momByRange_proton = -999.;
  reco_beam_momByRange_muon = -999.;

  reco_daughter_allTrack_momByRange_alt_proton.clear();
  reco_daughter_allTrack_momByRange_alt_muon.clear();
  reco_beam_momByRange_alt_proton = -999.;
  reco_beam_momByRange_alt_muon = -999.;

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
  true_beam_Pi0_decay_startPx.clear();
  true_beam_Pi0_decay_startPy.clear();
  true_beam_Pi0_decay_startPz.clear();
  true_beam_Pi0_decay_startX.clear();
  true_beam_Pi0_decay_startY.clear();
  true_beam_Pi0_decay_startZ.clear();
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

  reco_daughter_PFP_ID.clear();
  reco_daughter_PFP_nHits.clear();
  reco_daughter_PFP_nHits_collection.clear();
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

  reco_beam_dQdX_NoSCE.clear();
  reco_beam_dEdX_NoSCE.clear();
  reco_beam_calibrated_dEdX_NoSCE.clear();
  reco_beam_resRange_NoSCE.clear();
  reco_beam_TrkPitch_NoSCE.clear();

  reco_beam_dQdX_SCE.clear();
  reco_beam_dQ.clear();
  reco_beam_dEdX_SCE.clear();
  reco_beam_calibrated_dEdX_SCE.clear();
  reco_beam_vertex_nHits = -999;
  reco_beam_vertex_michel_score = -999.;

  reco_beam_resRange_SCE.clear();
  reco_beam_TrkPitch_SCE.clear();
  reco_beam_calo_wire.clear();
  reco_beam_calo_wire_z.clear();
  reco_beam_calo_tick.clear();
  reco_beam_calo_TPC.clear();

  reco_beam_trackID = -999;

  reco_beam_incidentEnergies.clear();
  reco_beam_interactingEnergy = -999.;
  true_beam_incidentEnergies.clear();
  true_beam_slices.clear();
  true_beam_slices_found.clear();
  true_beam_slices_deltaE.clear();
  true_beam_interactingEnergy = -999.;
  em_energy = 0.;
  true_beam_traj_X.clear();
  true_beam_traj_Y.clear();
  true_beam_traj_Z.clear();
  true_beam_traj_KE.clear();

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
  reco_daughter_PFP_true_byHits_completeness.clear();

  reco_daughter_PFP_true_byE_PDG.clear();
  reco_daughter_PFP_true_byE_len.clear();
  reco_daughter_PFP_true_byE_completeness.clear();
  reco_daughter_PFP_true_byE_purity.clear();



  reco_daughter_allTrack_ID.clear();
  reco_daughter_allTrack_dEdX_SCE.clear();
  reco_daughter_allTrack_dQdX_SCE.clear();
  reco_daughter_allTrack_resRange_SCE.clear();


  //Calorimetry + chi2 for planes 0 and 1
  reco_daughter_allTrack_calibrated_dEdX_SCE_plane0.clear();
  reco_daughter_allTrack_calibrated_dEdX_SCE_plane1.clear();

  reco_daughter_allTrack_resRange_plane0.clear();
  reco_daughter_allTrack_resRange_plane1.clear();

  reco_daughter_allTrack_Chi2_proton_plane0.clear();
  reco_daughter_allTrack_Chi2_proton_plane1.clear();

  reco_daughter_allTrack_Chi2_ndof_plane0.clear();
  reco_daughter_allTrack_Chi2_ndof_plane1.clear();
  ///////////////////////////////////////////


  //reco_daughter_allTrack_calibrated_dEdX.clear();
  reco_daughter_allTrack_calibrated_dEdX_SCE.clear();

  reco_daughter_allTrack_Chi2_proton.clear();
  reco_daughter_allTrack_Chi2_ndof.clear();

  reco_daughter_allTrack_Theta.clear();
  reco_daughter_allTrack_Phi.clear();
  reco_daughter_allTrack_len.clear();
  reco_daughter_allTrack_alt_len.clear();
  reco_daughter_allTrack_startX.clear();
  reco_daughter_allTrack_startY.clear();
  reco_daughter_allTrack_startZ.clear();
  reco_daughter_allTrack_endX.clear();
  reco_daughter_allTrack_endY.clear();
  reco_daughter_allTrack_endZ.clear();
  reco_daughter_allTrack_dR.clear();
  reco_daughter_allTrack_to_vertex.clear();
  reco_daughter_allTrack_vertex_michel_score.clear();
  reco_daughter_allTrack_vertex_nHits.clear();

  reco_daughter_allShower_ID.clear();
  reco_daughter_allShower_len.clear();
  reco_daughter_allShower_startX.clear();
  reco_daughter_allShower_startY.clear();
  reco_daughter_allShower_startZ.clear();

  reco_daughter_allShower_dirX.clear();
  reco_daughter_allShower_dirY.clear();
  reco_daughter_allShower_dirZ.clear();
  reco_daughter_allShower_energy.clear();
  ///////


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

  g4rw_primary_weights.clear();
  g4rw_primary_plus_sigma_weight.clear();
  g4rw_primary_minus_sigma_weight.clear();
  g4rw_primary_var.clear();
  g4rw_alt_primary_plus_sigma_weight.clear();
  g4rw_alt_primary_minus_sigma_weight.clear();
}


void pduneana::PDSPAnalyzer::BeamPFPInfo(
    const art::Event & evt, const recob::PFParticle* particle,
    anab::MVAReader<recob::Hit,4> * hitResults) {
  //Get CNN output for the beam
  reco_beam_PFP_ID = particle->Self();
  const std::vector< art::Ptr< recob::Hit > > beamPFP_hits = pfpUtil.GetPFParticleHits_Ptrs( *particle, evt, fPFParticleTag );
  reco_beam_PFP_nHits = beamPFP_hits.size();
  if (!fSkipMVA) {
    cnnOutput2D cnn = GetCNNOutputFromPFParticle( *particle, evt, *hitResults, pfpUtil, fPFParticleTag );
    reco_beam_PFP_trackScore = (cnn.nHits > 0 ?
                                (cnn.track / cnn.nHits) : -999.);
    reco_beam_PFP_emScore = (cnn.nHits > 0 ?
                             (cnn.em / cnn.nHits) : -999.);
    reco_beam_PFP_michelScore = (cnn.nHits > 0 ?
                                 (cnn.michel / cnn.nHits) : -999.);

    cnnOutput2D cnn_collection = GetCNNOutputFromPFParticleFromPlane( *particle, evt, *hitResults, pfpUtil, fPFParticleTag, 2 );
    reco_beam_PFP_trackScore_collection = (cnn_collection.nHits > 0 ?
                                           cnn_collection.track / cnn_collection.nHits : -999.);
    reco_beam_PFP_emScore_collection = (cnn_collection.nHits > 0 ?
                                        cnn_collection.em / cnn_collection.nHits : -999.);
    reco_beam_PFP_michelScore_collection = (cnn_collection.nHits > 0 ?
                                            cnn_collection.michel / cnn_collection.nHits : -999.);
  }
  else{
    reco_beam_PFP_trackScore =  -999.;
    reco_beam_PFP_emScore = -999.;
    reco_beam_PFP_michelScore = -999.;
    reco_beam_PFP_trackScore_collection =  -999.;
    reco_beam_PFP_emScore_collection = -999.;
    reco_beam_PFP_michelScore_collection = -999.;
  }

}

void pduneana::PDSPAnalyzer::BeamTrackInfo(
    const art::Event & evt, const recob::Track * thisTrack,
    detinfo::DetectorClocksData const& clockData) {
  auto allHits = evt.getValidHandle<std::vector<recob::Hit> >(fHitTag);

  if (!fSkipMVA) {
    std::pair<double, int> vertex_michel_score =
        trackUtil.GetVertexMichelScore(*thisTrack, evt, fTrackerTag, fHitTag/*,
                                       0., -500., 500., 0., 500., 0.*/);
    reco_beam_vertex_nHits = vertex_michel_score.second;
    reco_beam_vertex_michel_score = vertex_michel_score.first;
  }
  else {
    reco_beam_vertex_nHits = -999;
    reco_beam_vertex_michel_score = -999.;
  }

  if (fVerbose) std::cout << "Beam particle is track-like " << thisTrack->ID() << std::endl;
  reco_beam_type = 13;

  reco_beam_passes_beam_cuts = beam_cuts.IsBeamlike( *thisTrack, evt, "1" );
  if (fVerbose) std::cout << "Beam Cuts " << reco_beam_passes_beam_cuts << std::endl;


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
  trkf::TrackMomentumCalculator track_p_calc;
  reco_beam_momByRange_proton = track_p_calc.GetTrackMomentum(
      thisTrack->Length(), 2212);
  reco_beam_momByRange_muon = track_p_calc.GetTrackMomentum(
      thisTrack->Length(), 13);
  ////////////////////////////////////////////////////////////////


  std::map< const recob::Hit *, int > hitsToSlices;
  std::map< int, std::vector< const recob::Hit * > > slicesToHits;

  //Looking at the hits in the beam track
  std::map< size_t, const recob::Hit * > trajPtsToHits = trackUtil.GetRecoHitsFromTrajPoints( *thisTrack, evt, fTrackerTag );
  //double max_X = 0.;
  //double max_Y = 0.;
  //double max_Z = 0.;

  for( auto it = trajPtsToHits.begin(); it != trajPtsToHits.end(); ++it ){

    const recob::Hit * theHit = it->second;
    size_t i = it->first;

    double x = thisTrack->Trajectory().LocationAtPoint(i).X();
    double y = thisTrack->Trajectory().LocationAtPoint(i).Y();
    double z = thisTrack->Trajectory().LocationAtPoint(i).Z();

    if( fSaveHits ){
      //saving all hit coordinates for beamtrack
      reco_beam_spacePts_X.push_back(x);
      reco_beam_spacePts_Y.push_back(y);
      reco_beam_spacePts_Z.push_back(z);
    }

    //This creates the slices for the thin slice method.
    int slice = std::floor((z - fZ0) / fPitch);
    hitsToSlices[theHit] = slice;
    slicesToHits[slice].push_back(theHit);

    //if (z > max_Z){
    //  max_Z = z;
    //  max_X = y;
    //  max_Y = x;
    //}
  }

  //Last point is the vertex slice
  reco_beam_vertex_slice = slicesToHits.rbegin()->first;

  //Primary Track Calorimetry
  auto calo = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag,
                                                fCalorimetryTagSCE);
  bool found_calo = false;
  size_t index = 0;
  for ( index = 0; index < calo.size(); ++index) {
    if (calo[index].PlaneID().Plane == 2) {
      found_calo = true;
      break; 
    }
  }

  if (found_calo) {
    reco_beam_momByRange_alt_proton = track_p_calc.GetTrackMomentum(
        calo[index].Range(), 2212);
    reco_beam_momByRange_alt_muon = track_p_calc.GetTrackMomentum(
        calo[index].Range(), 13);
    reco_beam_alt_len = calo[index].Range();

    auto calo_dQdX = calo[index].dQdx();
    auto calo_dEdX = calo[index].dEdx();
    auto calo_range = calo[index].ResidualRange();
    auto TpIndices = calo[index].TpIndices();

    //For Prod3 only
    if (fCalorimetryTagSCE == "pandoracali") {
      auto pandoracalo = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, "pandoracalo");
      size_t this_index = 0;
      for ( this_index = 0; this_index < pandoracalo.size(); ++this_index) {
        if (pandoracalo[this_index].PlaneID().Plane == 2) {
          break; 
        }
      }
      std::cout << this_index << std::endl;
      TpIndices = pandoracalo[this_index].TpIndices();
      std::cout << "pandoracalo hits " << pandoracalo[this_index].dQdx().size() << std::endl;
    }
    std::cout << calo_dQdX.size() << std::endl;
    std::cout << calo[index].PlaneID().Plane << std::endl;

    auto theXYZPoints = calo[index].XYZ();
    std::vector< size_t > calo_hit_indices;
    for( size_t i = 0; i < calo_dQdX.size(); ++i ){
      if (fVerbose) std::cout << i << std::endl;
      reco_beam_dQdX_SCE.push_back( calo_dQdX[i] );
      reco_beam_dEdX_SCE.push_back( calo_dEdX[i] );
      reco_beam_resRange_SCE.push_back( calo_range[i] );
      reco_beam_TrkPitch_SCE.push_back( calo[index].TrkPitchVec()[i] );

      const recob::Hit & theHit = (*allHits)[ TpIndices[i] ];
      reco_beam_dQ.push_back(theHit.Integral());
      reco_beam_calo_TPC.push_back(theHit.WireID().TPC);
      if (theHit.WireID().TPC == 1) {
        reco_beam_calo_wire.push_back( theHit.WireID().Wire );
      }
      else if (theHit.WireID().TPC == 5) {
        reco_beam_calo_wire.push_back( theHit.WireID().Wire + 479);
      }
      //Need other TPCs?
      else {
        reco_beam_calo_wire.push_back(theHit.WireID().Wire );
      }
      reco_beam_calo_tick.push_back( theHit.PeakTime() );
      calo_hit_indices.push_back( TpIndices[i] );

      reco_beam_calo_wire_z.push_back(
          geom->Wire(theHit.WireID()).GetCenter().Z());

      if (fVerbose)
        std::cout << theXYZPoints[i].X() << " " << theXYZPoints[i].Y() << " " <<
                     theXYZPoints[i].Z() << " " << theHit.WireID().Wire << " " <<
                     geom->Wire(theHit.WireID()).GetCenter().Z() << " " <<
                     theHit.WireID().TPC << " " << std::endl;
    }

    //Getting the SCE corrected start/end positions & directions
    std::sort(theXYZPoints.begin(), theXYZPoints.end(), [](auto a, auto b)
        {return (a.Z() < b.Z());});

    //std::cout << theXYZPoints.size() << std:;endl;
    if (theXYZPoints.size()) {
      reco_beam_calo_startX = theXYZPoints[0].X();
      reco_beam_calo_startY = theXYZPoints[0].Y();
      reco_beam_calo_startZ = theXYZPoints[0].Z();
      reco_beam_calo_endX = theXYZPoints.back().X();
      reco_beam_calo_endY = theXYZPoints.back().Y();
      reco_beam_calo_endZ = theXYZPoints.back().Z();

      TVector3 dir((theXYZPoints.back().X() - theXYZPoints[0].X()),
                   (theXYZPoints.back().Y() - theXYZPoints[0].Y()),
                   (theXYZPoints.back().Z() - theXYZPoints[0].Z()));
      reco_beam_calo_startDirX.push_back(dir.Unit().X());
      reco_beam_calo_endDirX.push_back(dir.Unit().X());
      reco_beam_calo_startDirY.push_back(dir.Unit().Y());
      reco_beam_calo_endDirY.push_back(dir.Unit().Y());
      reco_beam_calo_startDirZ.push_back(dir.Unit().Z());
      reco_beam_calo_endDirZ.push_back(dir.Unit().Z());
    }
    else {
      reco_beam_calo_startDirX.push_back(-999.);
      reco_beam_calo_endDirX.push_back(-999.);
      reco_beam_calo_startDirY.push_back(-999.);
      reco_beam_calo_endDirY.push_back(-999.);
      reco_beam_calo_startDirZ.push_back(-999.);
      reco_beam_calo_endDirZ.push_back(-999.);
    }

    if (theXYZPoints.size() > 1) {
      TVector3 start_p1(theXYZPoints[0].X(),
          theXYZPoints[0].Y(), theXYZPoints[0].Z());
      TVector3 start_p2(theXYZPoints[1].X(),
          theXYZPoints[1].Y(), theXYZPoints[1].Z());
      TVector3 start_diff = start_p2 - start_p1;

      reco_beam_calo_startDirX.push_back(start_diff.Unit().X());
      reco_beam_calo_startDirY.push_back(start_diff.Unit().Y());
      reco_beam_calo_startDirZ.push_back(start_diff.Unit().Z());

      size_t nPoints = theXYZPoints.size();
      TVector3 end_p1(theXYZPoints[nPoints - 2].X(),
          theXYZPoints[nPoints - 2].Y(), theXYZPoints[nPoints - 2].Z());
      TVector3 end_p2(theXYZPoints[nPoints - 1].X(),
          theXYZPoints[nPoints - 1].Y(), theXYZPoints[nPoints - 1].Z());
      TVector3 end_diff = end_p2 - end_p1;

      reco_beam_calo_endDirX.push_back(end_diff.Unit().X());
      reco_beam_calo_endDirY.push_back(end_diff.Unit().Y());
      reco_beam_calo_endDirZ.push_back(end_diff.Unit().Z());
    }
    else {
      reco_beam_calo_startDirX.push_back(-999.);
      reco_beam_calo_endDirX.push_back(-999.);
      reco_beam_calo_startDirY.push_back(-999.);
      reco_beam_calo_endDirY.push_back(-999.);
      reco_beam_calo_startDirZ.push_back(-999.);
      reco_beam_calo_endDirZ.push_back(-999.);
    }

    if (theXYZPoints.size() > 2) {
      std::vector<TVector3> input;
      for (size_t iP = 0; iP < 3; ++iP) {
        input.push_back(TVector3(theXYZPoints[iP].X(),
                                 theXYZPoints[iP].Y(),
                                 theXYZPoints[iP].Z()));
      }

      TVector3 startDiff = FitLine(input);
      reco_beam_calo_startDirX.push_back(startDiff.Unit().X());
      reco_beam_calo_startDirY.push_back(startDiff.Unit().Y());
      reco_beam_calo_startDirZ.push_back(startDiff.Unit().Z());

      std::vector<TVector3> end_input;
      size_t nPoints = theXYZPoints.size();
      for (size_t iP = nPoints - 3; iP < nPoints; ++iP) {
        end_input.push_back(TVector3(theXYZPoints[iP].X(),
                                     theXYZPoints[iP].Y(),
                                     theXYZPoints[iP].Z()));
      }

      TVector3 endDiff = FitLine(end_input);
      reco_beam_calo_endDirX.push_back(endDiff.Unit().X());
      reco_beam_calo_endDirY.push_back(endDiff.Unit().Y());
      reco_beam_calo_endDirZ.push_back(endDiff.Unit().Z());
    }
    else {
      reco_beam_calo_startDirX.push_back(-999.);
      reco_beam_calo_endDirX.push_back(-999.);
      reco_beam_calo_startDirY.push_back(-999.);
      reco_beam_calo_endDirY.push_back(-999.);
      reco_beam_calo_startDirZ.push_back(-999.);
      reco_beam_calo_endDirZ.push_back(-999.);
    }

    if (theXYZPoints.size() > 3) {
      std::vector<TVector3> input;
      for (size_t iP = 0; iP < 4; ++iP) {
        input.push_back(TVector3(theXYZPoints[iP].X(),
                                 theXYZPoints[iP].Y(),
                                 theXYZPoints[iP].Z()));
      }

      TVector3 startDiff = FitLine(input);
      reco_beam_calo_startDirX.push_back(startDiff.Unit().X());
      reco_beam_calo_startDirY.push_back(startDiff.Unit().Y());
      reco_beam_calo_startDirZ.push_back(startDiff.Unit().Z());

      std::vector<TVector3> end_input;
      size_t nPoints = theXYZPoints.size();
      for (size_t iP = nPoints - 4; iP < nPoints; ++iP) {
        end_input.push_back(TVector3(theXYZPoints[iP].X(),
                                     theXYZPoints[iP].Y(),
                                     theXYZPoints[iP].Z()));
      }

      TVector3 endDiff = FitLine(end_input);
      reco_beam_calo_endDirX.push_back(endDiff.Unit().X());
      reco_beam_calo_endDirY.push_back(endDiff.Unit().Y());
      reco_beam_calo_endDirZ.push_back(endDiff.Unit().Z());

    }
    else {
      reco_beam_calo_startDirX.push_back(-999.);
      reco_beam_calo_endDirX.push_back(-999.);
      reco_beam_calo_startDirY.push_back(-999.);
      reco_beam_calo_endDirY.push_back(-999.);
      reco_beam_calo_startDirZ.push_back(-999.);
      reco_beam_calo_endDirZ.push_back(-999.);
    }
    ////////////////////////////////////////////

    //New Calibration
    std::cout << "Getting reco beam calo" << std::endl;
    std::vector< float > new_dEdX = calibration_SCE.GetCalibratedCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTagSCE, 2, -1.);
    std::cout << new_dEdX.size() << " " << reco_beam_resRange_SCE.size() << std::endl;
    for( size_t i = 0; i < new_dEdX.size(); ++i ){ reco_beam_calibrated_dEdX_SCE.push_back( new_dEdX[i] ); }
    std::cout << "got calibrated dedx" << std::endl;
    ////////////////////////////////////////////

    std::pair< double, int > pid_chi2_ndof = trackUtil.Chi2PID( reco_beam_calibrated_dEdX_SCE, reco_beam_resRange_SCE, templates[ 2212 ] );
    std::cout << "got chi2" << std::endl;
    reco_beam_Chi2_proton = pid_chi2_ndof.first;
    reco_beam_Chi2_ndof = pid_chi2_ndof.second;

    if (fVerbose)
      std::cout << "Calo check: " << reco_beam_calibrated_dEdX_SCE.size() << " " <<
                   reco_beam_TrkPitch_SCE.size() << std::endl;
    std::vector< calo_point > reco_beam_calo_points;
    //Doing thin slice
    if (reco_beam_calibrated_dEdX_SCE.size() &&
        reco_beam_calibrated_dEdX_SCE.size() == reco_beam_TrkPitch_SCE.size() &&
        reco_beam_calibrated_dEdX_SCE.size() == reco_beam_calo_wire.size()) {

      for( size_t i = 0; i < reco_beam_calibrated_dEdX_SCE.size(); ++i ){
        reco_beam_calo_points.push_back(
          calo_point(reco_beam_calo_wire[i], reco_beam_TrkPitch_SCE[i],
                     reco_beam_dQdX_SCE[i], reco_beam_dEdX_SCE[i],
                     reco_beam_dQ[i],
                     reco_beam_calibrated_dEdX_SCE[i],
                     reco_beam_resRange_SCE[i], calo_hit_indices[i],
                     reco_beam_calo_wire_z[i], reco_beam_calo_TPC[i]));
      }

      //std::cout << "N Calo points: " << reco_beam_calo_points.size() << std::endl;
      //Sort
      std::sort( reco_beam_calo_points.begin(), reco_beam_calo_points.end(), [](calo_point a, calo_point b) {return ( a.z < b.z );} );

      //And also put these in the right order
      for( size_t i = 0; i < reco_beam_calo_points.size(); ++i ){
        calo_point thePoint = reco_beam_calo_points[i];
        reco_beam_calo_wire[i] = thePoint.wire;
        reco_beam_calibrated_dEdX_SCE[i] = thePoint.calibrated_dEdX;
        reco_beam_TrkPitch_SCE[i] = thePoint.pitch;
        calo_hit_indices[i] = thePoint.hit_index;
        reco_beam_calo_wire_z[i] = thePoint.z;
        reco_beam_calo_TPC[i] = thePoint.tpc;
        reco_beam_resRange_SCE[i] = thePoint.res_range;
        reco_beam_dQdX_SCE[i] = thePoint.dQdX;
        reco_beam_dQ[i] = thePoint.dQ;
        reco_beam_dEdX_SCE[i] = thePoint.dEdX;
      }

      //Get the initial Energy KE
      //double mass = 0.;
      double init_KE = 0.;
      //std::cout << "Has BI? " << fMCHasBI << " " << evt.isRealData() << std::endl;
      if (evt.isRealData() || fMCHasBI) {
        double mass = 139.57;

        init_KE =  sqrt( 1.e6*beam_inst_P*beam_inst_P + mass*mass ) - mass;
       // std::cout << "MC has BI: " << init_KE << std::endl;
      }
      else{
        init_KE = sqrt(1.e6*true_beam_startP*true_beam_startP +
                       true_beam_mass*true_beam_mass) - true_beam_mass;
      }

      reco_beam_incidentEnergies.push_back( init_KE );
      for( size_t i = 0; i < reco_beam_calo_points.size() - 1; ++i ){ //-1 to not count the last slice
        //use dedx * pitch or new hit calculation?
        double this_energy = reco_beam_incidentEnergies.back() - ( reco_beam_calo_points[i].calibrated_dEdX * reco_beam_calo_points[i].pitch );
        reco_beam_incidentEnergies.push_back( this_energy );
      }
      if( reco_beam_incidentEnergies.size() ) reco_beam_interactingEnergy = reco_beam_incidentEnergies.back();
    }
  }

  //Uncorrected
  auto calo_NoSCE = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt,
                                                      fTrackerTag,
                                                      fCalorimetryTagNoSCE);
  found_calo = false;
  index = 0;
  for (index = 0; index < calo_NoSCE.size(); ++index) {
    if (calo_NoSCE[index].PlaneID().Plane == 2) {
      found_calo = true;
      break; 
    }
  }
  std::cout << "Cali index " << index << std::endl;

  if (found_calo) {
    auto calo_dQdX = calo_NoSCE[index].dQdx();
    auto calo_dEdX = calo_NoSCE[index].dEdx();
    auto calo_range = calo_NoSCE[index].ResidualRange();
    auto TpIndices = calo_NoSCE[index].TpIndices();

    std::vector< size_t > calo_hit_indices;
    for( size_t i = 0; i < calo_dQdX.size(); ++i ){
      if (fVerbose) std::cout << i << std::endl;
      reco_beam_dQdX_NoSCE.push_back( calo_dQdX[i] );
      reco_beam_dEdX_NoSCE.push_back( calo_dEdX[i] );
      reco_beam_resRange_NoSCE.push_back( calo_range[i] );
      reco_beam_TrkPitch_NoSCE.push_back( calo_NoSCE[index].TrkPitchVec()[i] );
      calo_hit_indices.push_back( TpIndices[i] );
    }

    std::vector< float > new_dEdX = calibration_NoSCE.GetCalibratedCalorimetry(  *thisTrack, evt, fTrackerTag, fCalorimetryTagNoSCE, 2, -1.);
    for( size_t i = 0; i < new_dEdX.size(); ++i ){ reco_beam_calibrated_dEdX_NoSCE.push_back( new_dEdX[i] ); }
    std::cout << "got calibrated dedx" << std::endl;
    ////////////////////////////////////////////

    std::vector< calo_point > reco_beam_calo_points;
    std::cout << "wire: " << reco_beam_calo_wire.size() << std::endl;
    std::cout << "trkpitch: " << reco_beam_TrkPitch_NoSCE.size() << std::endl;
    std::cout << "dqdx: " << reco_beam_dQdX_NoSCE.size() << std::endl;
    std::cout << "dedx: " << reco_beam_dEdX_NoSCE.size() << std::endl;
    std::cout << "cal dedx: " << reco_beam_calibrated_dEdX_NoSCE.size() << std::endl;
    std::cout << "ind: " << calo_hit_indices.size() << std::endl;
    std::cout << "range: " << reco_beam_resRange_NoSCE.size() << std::endl;
    std::cout << "wire_z: " << reco_beam_calo_wire_z.size() << std::endl;
    std::cout << "TPC: " << reco_beam_calo_TPC.size() << std::endl;
    if (reco_beam_calibrated_dEdX_NoSCE.size() &&
        reco_beam_calibrated_dEdX_NoSCE.size() == reco_beam_TrkPitch_NoSCE.size() &&
        reco_beam_calibrated_dEdX_NoSCE.size() == reco_beam_calo_wire.size()) {

      for( size_t i = 0; i < reco_beam_calibrated_dEdX_NoSCE.size(); ++i ){
        reco_beam_calo_points.push_back(
          calo_point(reco_beam_calo_wire[i], reco_beam_TrkPitch_NoSCE[i],
                     reco_beam_dQdX_NoSCE[i], reco_beam_dEdX_NoSCE[i],
                     reco_beam_dQ[i],
                     reco_beam_calibrated_dEdX_NoSCE[i],
                     calo_hit_indices[i], reco_beam_resRange_NoSCE[i],
                     reco_beam_calo_wire_z[i], reco_beam_calo_TPC[i]));
      }

      //std::cout << "N Calo points: " << reco_beam_calo_points.size() << std::endl;
      //Sort
      std::sort( reco_beam_calo_points.begin(), reco_beam_calo_points.end(), [](calo_point a, calo_point b) {return ( a.z < b.z );} );

      //And also put these in the right order
      for( size_t i = 0; i < reco_beam_calo_points.size(); ++i ){
        calo_point thePoint = reco_beam_calo_points[i];
        reco_beam_calo_wire[i] = thePoint.wire;
        reco_beam_calibrated_dEdX_NoSCE[i] = thePoint.calibrated_dEdX;
        reco_beam_TrkPitch_NoSCE[i] = thePoint.pitch;
        calo_hit_indices[i] = thePoint.hit_index;
        reco_beam_calo_wire_z[i] = thePoint.z;
        reco_beam_calo_TPC[i] = thePoint.tpc;
        reco_beam_resRange_NoSCE[i] = thePoint.res_range;
        reco_beam_dQdX_NoSCE[i] = thePoint.dQdX;
        reco_beam_dEdX_NoSCE[i] = thePoint.dEdX;
      }
    }
  }
}

void pduneana::PDSPAnalyzer::BeamShowerInfo(const recob::Shower* thisShower) {
  reco_beam_type = 11;
  reco_beam_trackID = thisShower->ID();

  if (fVerbose) {
    std::cout << "Beam particle is shower-like" << std::endl;
    std::cout << thisShower->ShowerStart().X() << " " << thisShower->ShowerStart().Y() << " " << thisShower->ShowerStart().Z() << std::endl;
    std::cout << thisShower->Direction().X() << " " << thisShower->Direction().Y() << " " << thisShower->Direction().Z() << std::endl;
  }
}


void pduneana::PDSPAnalyzer::TrueBeamInfo(
    const art::Event & evt, const simb::MCParticle* true_beam_particle,
    detinfo::DetectorClocksData const& clockData,
    const sim::ParticleList & plist,
    std::map<int, std::vector<int>> & trueToPFPs,
    anab::MVAReader<recob::Hit,4> * hitResults) {
  auto pfpVec
      = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  true_beam_endProcess = true_beam_particle->EndProcess();

  true_beam_PDG         = true_beam_particle->PdgCode();
  true_beam_mass        = true_beam_particle->Mass()*1.e3;
  true_beam_ID          = true_beam_particle->TrackId();
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
  true_beam_endP2     = true_beam_particle->P(true_np-1);

  true_beam_last_len
      = sqrt(std::pow((true_beam_particle->Position(true_np-2).X() -
                       true_beam_particle->Position(true_np-1).X()), 2) +
             std::pow((true_beam_particle->Position(true_np-2).Y() -
                       true_beam_particle->Position(true_np-1).Y()), 2) +
             std::pow((true_beam_particle->Position(true_np-2).Z() -
                       true_beam_particle->Position(true_np-1).Z()), 2));

  true_beam_startDirX  = true_beam_startPx / true_beam_startP;
  true_beam_startDirY  = true_beam_startPy / true_beam_startP;
  true_beam_startDirZ  = true_beam_startPz / true_beam_startP;

  true_beam_nHits = truthUtil.GetMCParticleHits( clockData, *true_beam_particle, evt, fHitTag ).size();
  true_beam_reco_byHits_PFP_ID.push_back( std::vector< int >() );
  true_beam_reco_byHits_PFP_nHits.push_back( std::vector< int >() );
  true_beam_reco_byHits_allTrack_ID.push_back( std::vector< int >() );
  if( fTrueToReco ){
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
        MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Track object not found, moving on" << std::endl;
      }

      if( pandora2Track ){
        true_beam_reco_byHits_allTrack_ID.back().push_back( pandora2Track->ID() );
      }
      else{
        true_beam_reco_byHits_allTrack_ID.back().push_back( -1 );
      }
    }
  }

  //Truth thin slice info
  //Go through the true processes within the MCTrajectory
  const simb::MCTrajectory & true_beam_trajectory = true_beam_particle->Trajectory();
  auto true_beam_proc_map = true_beam_trajectory.TrajectoryProcesses();
  if (fVerbose) std::cout << "Processes: " << std::endl;

  for( auto itProc = true_beam_proc_map.begin(); itProc != true_beam_proc_map.end(); ++itProc ){
    int index = itProc->first;
    std::string process = true_beam_trajectory.KeyToProcess(itProc->second);
    if (fVerbose) std::cout << index << " " << process << std::endl;

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

      double total_E = sqrt(total_P*total_P*1.e6 +
                            true_beam_mass*true_beam_mass);
      double total_next_E = sqrt(total_next_P*total_next_P*1.e6 +
                                 true_beam_mass*true_beam_mass);

      true_beam_elastic_X.push_back( process_X );
      true_beam_elastic_Y.push_back( process_Y );
      true_beam_elastic_Z.push_back( process_Z );

      true_beam_elastic_deltaE.push_back(total_E - total_next_E);

      std::vector<const sim::IDE *> ides_between_points =
          truthUtil.GetSimIDEsBetweenPoints(
              *true_beam_particle, true_beam_trajectory.Position(index),
              true_beam_trajectory.Position(index +1));

      double total_edep = 0.;
      for (size_t i = 0; i < ides_between_points.size(); ++i) {
        total_edep += ides_between_points[i]->energy;
      }
      true_beam_elastic_IDE_edep.push_back(total_edep);

    }
  }
  if( true_beam_endProcess.find( "Inelastic" ) == std::string::npos ){
    true_beam_processes.push_back( true_beam_endProcess );
  }

  if (fVerbose) std::cout << "Looking at IDEs" << std::endl;

  auto view2_IDEs = bt_serv->TrackIdToSimIDEs_Ps( true_beam_ID, geo::View_t(2) );

  if (fVerbose) std::cout << "N view2 IDEs: " << view2_IDEs.size() << std::endl;

  true_beam_IDE_totalDep = 0.;
  for (size_t i = 0; i < view2_IDEs.size(); ++i) {
    const sim::IDE * ide = view2_IDEs[i];
    std::cout << " " << ide->trackID << " " << ide->z << " " << ide->energy <<
                 std::endl;

    true_beam_IDE_totalDep += ide->energy;
  }
  //Sort based on ide z-position
  std::sort( view2_IDEs.begin(), view2_IDEs.end(), sort_IDEs );
  std::cout << "Sorted" << std::endl;

  //This is an attempt to remove IDEs from things like delta-rays
  //that have a large gap in z to the previous
  size_t remove_index = 0;
  bool   do_remove = false;
  if( view2_IDEs.size() ){
    for( size_t i = 1; i < view2_IDEs.size()-1; ++i ){
      const sim::IDE * prev_IDE = view2_IDEs[i-1];
      const sim::IDE * this_IDE = view2_IDEs[i];

      if (this_IDE->trackID < 0.) {
        em_energy += this_IDE->energy;
      }


      if( this_IDE->trackID < 0 && ( this_IDE->z - prev_IDE->z ) > 5 ){
        remove_index = i;
        do_remove = true;
        break;   
      }
    }
  }

  if( do_remove ){
    view2_IDEs.erase( view2_IDEs.begin() + remove_index, view2_IDEs.end() );
  }

  double init_KE = sqrt(1.e6 * true_beam_startP*true_beam_startP +
                        true_beam_mass*true_beam_mass) - true_beam_mass;

  //slice up the view2_IDEs up by the wire pitch
  auto sliced_ides = slice_IDEs( view2_IDEs, fZ0, fPitch, true_beam_endZ);
  //Get the momentum at the start of the slices.
  //Get the first slice
  if (fVerbose) std::cout << "size: " << sliced_ides.size() << std::endl;
  if (sliced_ides.size()) {
    auto first_slice = sliced_ides.begin();
    if (fVerbose) std::cout << "Got first slice" << std::endl;

    //Check it has any IDEs
    auto theIDEs = first_slice->second;
    if (fVerbose) std::cout << "Got ides" << std::endl;
    if (theIDEs.size()) {
      //Get the first ide z position
      double ide_z = theIDEs[0]->z;

      //Go through the trajectory position
      //and check for the position that comes immediately before the
      //first ide
      for (size_t i = 1; i < true_beam_trajectory.size(); ++i) {
        double z0 = true_beam_trajectory.Z(i-1);
        double z1 = true_beam_trajectory.Z(i);

        double x0 = true_beam_trajectory.X(i-1);
        double y0 = true_beam_trajectory.Y(i-1);
        double x1 = true_beam_trajectory.X(i);
        double y1 = true_beam_trajectory.Y(i);
        geo::Point_t test_point{x0, y0, z0};
        geo::Point_t test_point_1{x1, y1, z1};
        const TGeoMaterial * mat = geom->Material(test_point);
        std::cout << "Point " << i-1 << " (" << x0 << ", " << y0 <<
                     ", " << z0 << ") Material " << mat->GetName() <<
                     std::endl;
        if (z0 < ide_z && z1 > ide_z) {
          const TGeoMaterial * mat1 = geom->Material(test_point_1);
          std::cout << "Point " << i << " (" << x1 << ", " << y1 <<
                       ", " << z1 << ") Material " << mat1->GetName() <<
                       std::endl;
          init_KE = 1.e3 * true_beam_trajectory.E(i-1) - true_beam_mass;
          if (fVerbose) {
            std::cout << "Found matching position" << z0 << " " << ide_z <<
                         " " << z1 << std::endl;
            std::cout << "init KE: " << init_KE << " " <<
                         (1.e3*true_beam_trajectory.E(i) - true_beam_mass) << std::endl;
          }
          break;
        }
      }
    }
  }

  //Go through the sliced up IDEs to create the thin targets 
  true_beam_incidentEnergies.push_back( init_KE );
  for( auto it = sliced_ides.begin(); it != sliced_ides.end(); ++it ){

    auto theIDEs = it->second;

    true_beam_slices.push_back( it->first );

    double deltaE = 0.;
    for( size_t i = 0; i < theIDEs.size(); ++i ){
      deltaE += theIDEs[i]->energy;
    }

    true_beam_slices_deltaE.push_back( deltaE );
    true_beam_incidentEnergies.push_back( true_beam_incidentEnergies.back() - deltaE );
  }
  true_beam_incidentEnergies.pop_back();
  if( true_beam_incidentEnergies.size() ) true_beam_interactingEnergy = true_beam_incidentEnergies.back();

  //Save the trajectory points
  for (size_t i = 0; i < true_beam_trajectory.size(); ++i) {
    true_beam_traj_X.push_back(true_beam_trajectory.X(i));
    true_beam_traj_Y.push_back(true_beam_trajectory.Y(i));
    true_beam_traj_Z.push_back(true_beam_trajectory.Z(i));
    
    true_beam_traj_KE.push_back(true_beam_trajectory.E(i)*1.e3 - true_beam_mass);
  }

  //Look through the daughters
  for( int i = 0; i < true_beam_particle->NumberDaughters(); ++i ){
    int daughterID = true_beam_particle->Daughter(i);

    if (fVerbose) std::cout << "Daughter " << i << " ID: " << daughterID << std::endl;
    auto part = plist[ daughterID ];
    int pid = part->PdgCode();

    std::string process = part->Process();

    if (process == "muIoni" || process == "hIoni")
      continue;

    true_beam_daughter_PDG.push_back(pid);
    true_beam_daughter_ID.push_back( part->TrackId() );

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

    if (fVerbose) {
      std::cout << "Process: " << part->Process() << std::endl;
      std::cout << "PID: " << pid << std::endl;
      std::cout << "Start: " << part->Position(0).X() << " " << part->Position(0).Y() << " " << part->Position(0).Z() << std::endl;
      std::cout << "End: " << part->EndPosition().X() << " " << part->EndPosition().Y() << " " << part->EndPosition().Z() << std::endl;
      std::cout << "Len: " << part->Trajectory().TotalLength() << std::endl;
    }

    if( part->Process().find( "Inelastic" ) != std::string::npos ){
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
        true_beam_Pi0_decay_startPx.push_back( pi0_decay_part->Px() );
        true_beam_Pi0_decay_startPy.push_back( pi0_decay_part->Py() );
        true_beam_Pi0_decay_startPz.push_back( pi0_decay_part->Pz() );
        true_beam_Pi0_decay_startX.push_back( pi0_decay_part->Position(0).X() );
        true_beam_Pi0_decay_startY.push_back( pi0_decay_part->Position(0).Y() );
        true_beam_Pi0_decay_startZ.push_back( pi0_decay_part->Position(0).Z() );
        true_beam_Pi0_decay_parID.push_back( pi0_decay_part->Mother() );

        true_beam_Pi0_decay_len.push_back( pi0_decay_part->Trajectory().TotalLength() );
        true_beam_Pi0_decay_nHits.push_back( truthUtil.GetMCParticleHits( clockData, *pi0_decay_part, evt, fHitTag ).size() );

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
        if( fTrueToReco ){
          for( size_t k = 0; k < trueToPFPs[ pi0_decay_part->TrackId() ].size(); ++k ){
            true_beam_Pi0_decay_reco_byHits_PFP_ID.back().push_back( trueToPFPs[ pi0_decay_part->TrackId() ][k] );

            const recob::PFParticle * thePFP = &(pfpVec->at( trueToPFPs[ pi0_decay_part->TrackId() ][k] ));
            true_beam_Pi0_decay_reco_byHits_PFP_nHits.back().push_back(
              pfpUtil.GetPFParticleHits_Ptrs( *thePFP, evt, fPFParticleTag ).size()
            );

            if (!fSkipMVA) {
              cnnOutput2D theCNNResults = GetCNNOutputFromPFParticle( *thePFP, evt, *hitResults, pfpUtil, fPFParticleTag );
              double track_score = ((theCNNResults.nHits > 0) ? 
                                    (theCNNResults.track / theCNNResults.nHits) :
                                    -999.);
              true_beam_Pi0_decay_reco_byHits_PFP_trackScore.back().push_back(track_score);
            }
            else {
              true_beam_Pi0_decay_reco_byHits_PFP_trackScore.back().push_back(-999.);
            }

            const recob::Track* pandora2Track = 0x0;
            try{
              pandora2Track = pfpUtil.GetPFParticleTrack( *thePFP, evt, fPFParticleTag, "pandora2Track" );
            }
            catch( const cet::exception &e ){
              MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Track object not found, moving on" << std::endl;
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
              true_beam_Pi0_decay_reco_byHits_allTrack_ID.back().push_back( -999 );
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
              MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Shower object not found, moving on" << std::endl;
            }

            if( pandora2Shower ){
              true_beam_Pi0_decay_reco_byHits_allShower_ID.back().push_back( pandora2Shower->ID() );
              true_beam_Pi0_decay_reco_byHits_allShower_startX.back().push_back( pandora2Shower->ShowerStart().X() );
              true_beam_Pi0_decay_reco_byHits_allShower_startY.back().push_back( pandora2Shower->ShowerStart().Y() );
              true_beam_Pi0_decay_reco_byHits_allShower_startZ.back().push_back( pandora2Shower->ShowerStart().Z() );
              true_beam_Pi0_decay_reco_byHits_allShower_len.back().push_back( pandora2Shower->Length() );
            }
            else{
              true_beam_Pi0_decay_reco_byHits_allShower_ID.back().push_back( -999 );
              true_beam_Pi0_decay_reco_byHits_allShower_startX.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allShower_startY.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allShower_startZ.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allShower_len.back().push_back( -999. );
            }

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
      true_beam_grand_daughter_nHits.push_back( truthUtil.GetMCParticleHits( clockData, *grand_daughter_part, evt, fHitTag ).size() );
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
    if( fTrueToReco ){
      for( size_t j = 0; j < trueToPFPs[ part->TrackId() ].size(); ++j ){
        true_beam_daughter_reco_byHits_PFP_ID.back().push_back( trueToPFPs[ part->TrackId() ][j] );

        const recob::PFParticle * thePFP = &(pfpVec->at( trueToPFPs[ part->TrackId() ][j] ));
        true_beam_daughter_reco_byHits_PFP_nHits.back().push_back(
          pfpUtil.GetPFParticleHits_Ptrs( *thePFP, evt, fPFParticleTag ).size()
        );

        if (!fSkipMVA) {
          cnnOutput2D theCNNResults = GetCNNOutputFromPFParticle( *thePFP, evt, *hitResults, pfpUtil, fPFParticleTag );
          double track_score = ((theCNNResults.nHits > 0) ? 
                                (theCNNResults.track / theCNNResults.nHits) :
                                -999.);
          true_beam_daughter_reco_byHits_PFP_trackScore.back().push_back(track_score);
        }
        else {
          true_beam_daughter_reco_byHits_PFP_trackScore.back().push_back(-999.);
        }

        //cnnOutput2D theCNNResults = GetCNNOutputFromPFParticle( *thePFP, evt, *hitResults, pfpUtil, fPFParticleTag );
        //true_beam_daughter_reco_byHits_PFP_trackScore.back().push_back( ( ( theCNNResults.nHits > 0 ) ? ( theCNNResults.track / theCNNResults.nHits ) : -999. ) );

        const recob::Track* pandora2Track = 0x0;
        try{
           pandora2Track = pfpUtil.GetPFParticleTrack( *thePFP, evt, fPFParticleTag, "pandora2Track" );
        }
        catch( const cet::exception &e ){
          MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Track object not found, moving on" << std::endl;
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
          true_beam_daughter_reco_byHits_allTrack_ID.back().push_back( -999 );
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
          MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Shower object not found, moving on" << std::endl;
        }
 
        if( pandora2Shower ){
          true_beam_daughter_reco_byHits_allShower_ID.back().push_back( pandora2Shower->ID() );
          true_beam_daughter_reco_byHits_allShower_startX.back().push_back( pandora2Shower->ShowerStart().X() );
          true_beam_daughter_reco_byHits_allShower_startY.back().push_back( pandora2Shower->ShowerStart().Y() );
          true_beam_daughter_reco_byHits_allShower_startZ.back().push_back( pandora2Shower->ShowerStart().Z() );
          true_beam_daughter_reco_byHits_allShower_len.back().push_back( pandora2Shower->Length() );
 
        }
        else{
          true_beam_daughter_reco_byHits_allShower_ID.back().push_back( -999 );
          true_beam_daughter_reco_byHits_allShower_startX.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allShower_startY.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allShower_startZ.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allShower_len.back().push_back( -999. );
        }

      }
    }
    true_beam_daughter_nHits.push_back( truthUtil.GetMCParticleHits( clockData, *part, evt, fHitTag ).size() );

  }
}

void pduneana::PDSPAnalyzer::BeamInstInfo(const art::Event & evt) {
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  if (evt.isRealData()) {
    auto beamHandle
        = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(
            fBeamModuleLabel);

    if (beamHandle.isValid()) {
      art::fill_ptr_vector(beamVec, beamHandle);
    }
  }
  else {
    try {
      auto beamHandle
          = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(
              "generator");

      if( beamHandle.isValid()){
        art::fill_ptr_vector(beamVec, beamHandle);
      }
      fMCHasBI = true;
    }
    catch (const cet::exception &e) {
      MF_LOG_WARNING("PDSPAnalyzer") << "BeamEvent generator object not " <<
                                        "found, moving on" << std::endl;
      fMCHasBI = false;
      return;
    }
  }

  const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); 
  if (evt.isRealData() && !fBeamlineUtils.IsGoodBeamlineTrigger(evt)) {
    MF_LOG_WARNING("PDSPAnalyzer") << "Failed beam quality check" << std::endl;
    beam_inst_valid = false;
    return;
  }

  int nTracks = beamEvent.GetBeamTracks().size();
  std::vector< double > momenta = beamEvent.GetRecoBeamMomenta();
  int nMomenta = momenta.size();

  if (fVerbose) {
    std::cout << "Got beam event" << std::endl;
    std::cout << "Got " << nTracks << " Tracks" << std::endl;
    std::cout << "Got " << nMomenta << " Momenta" << std::endl;
  }

  if( nMomenta > 0 ){
    beam_inst_P = momenta[0];
    if (fVerbose) std::cout << "reco P " << beam_inst_P << std::endl;
  }

  const std::vector<double> the_tofs = beamEvent.GetTOFs();
  const std::vector<int> the_chans = beamEvent.GetTOFChans();
  for (size_t iTOF = 0; iTOF < the_tofs.size(); ++iTOF) {
    beam_inst_TOF.push_back(the_tofs[iTOF]);
    beam_inst_TOF_Chan.push_back(the_chans[iTOF]);
  }

  if( nTracks > 0 ){
    beam_inst_X = beamEvent.GetBeamTracks()[0].Trajectory().End().X();
    beam_inst_Y = beamEvent.GetBeamTracks()[0].Trajectory().End().Y();
    beam_inst_Z = beamEvent.GetBeamTracks()[0].Trajectory().End().Z();

    beam_inst_dirX = beamEvent.GetBeamTracks()[0].Trajectory().EndDirection().X();
    beam_inst_dirY = beamEvent.GetBeamTracks()[0].Trajectory().EndDirection().Y();
    beam_inst_dirZ = beamEvent.GetBeamTracks()[0].Trajectory().EndDirection().Z();
  }

  beam_inst_nTracks = nTracks;
  beam_inst_nMomenta = nMomenta;

  if (evt.isRealData()) {
    std::vector< int > pdg_cands = fBeamlineUtils.GetPID( beamEvent, 1. );
    beam_inst_PDG_candidates.insert( beam_inst_PDG_candidates.end(), pdg_cands.begin(), pdg_cands.end() );
  }

  beam_inst_nFibersP1 = beamEvent.GetActiveFibers( "XBPF022697" ).size();
  beam_inst_nFibersP2 = beamEvent.GetActiveFibers( "XBPF022701" ).size();
  beam_inst_nFibersP3 = beamEvent.GetActiveFibers( "XBPF022702" ).size();
}

void pduneana::PDSPAnalyzer::DaughterPFPInfo(
    const art::Event & evt, const recob::PFParticle* particle,
    detinfo::DetectorClocksData const& clockData,
    anab::MVAReader<recob::Hit,4> * hitResults) {

  auto pfpVec
      = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  auto const detProp
      = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(
          evt, clockData);
  trkf::TrackMomentumCalculator track_p_calc;
  for (size_t daughterID : particle->Daughters()) {
    const recob::PFParticle * daughterPFP = &(pfpVec->at( daughterID ));
    reco_daughter_PFP_ID.push_back( daughterID );

    const std::vector< art::Ptr< recob::Hit > > daughterPFP_hits = pfpUtil.GetPFParticleHits_Ptrs( *daughterPFP, evt, fPFParticleTag );
    if (fVerbose) std::cout << "Got " << daughterPFP_hits.size() << " hits from daughter " << daughterID << std::endl;

    reco_daughter_PFP_nHits.push_back( daughterPFP_hits.size() );
    size_t nHits_coll = 0;
    for (size_t i = 0; i < daughterPFP_hits.size(); ++i) {
      if (daughterPFP_hits[i]->View() == 2) {
        ++nHits_coll;
      }
    }
    reco_daughter_PFP_nHits_collection.push_back(nHits_coll);

    if (!fSkipMVA) {
      cnnOutput2D theCNNResults = GetCNNOutputFromPFParticle(
          *daughterPFP, evt, *hitResults, pfpUtil, fPFParticleTag);
      double track_score = (theCNNResults.nHits > 0 ?
                            (theCNNResults.track / theCNNResults.nHits) :
                            -999.);
      double em_score = (theCNNResults.nHits > 0 ?
                            (theCNNResults.em / theCNNResults.nHits) :
                            -999.);
      double michel_score = (theCNNResults.nHits > 0 ?
                            (theCNNResults.michel / theCNNResults.nHits) :
                            -999.);
      reco_daughter_PFP_trackScore.push_back(track_score);
      reco_daughter_PFP_emScore.push_back(em_score);
      reco_daughter_PFP_michelScore.push_back(michel_score);

      cnnOutput2D cnn_collection = GetCNNOutputFromPFParticleFromPlane(
          *daughterPFP, evt, *hitResults, pfpUtil, fPFParticleTag, 2);
      double track_score_collection = (cnn_collection.nHits > 0 ?
                            (cnn_collection.track / cnn_collection.nHits) :
                            -999.);
      double em_score_collection = (cnn_collection.nHits > 0 ?
                            (cnn_collection.em / cnn_collection.nHits) :
                            -999.);
      double michel_score_collection = (cnn_collection.nHits > 0 ?
                            (cnn_collection.michel / cnn_collection.nHits) :
                            -999.);
      reco_daughter_PFP_trackScore_collection.push_back(track_score_collection);
      reco_daughter_PFP_emScore_collection.push_back(em_score_collection);
      reco_daughter_PFP_michelScore_collection.push_back(michel_score_collection);
    }
    else{
      reco_daughter_PFP_trackScore.push_back( -999. );
      reco_daughter_PFP_emScore.push_back( -999. );
      reco_daughter_PFP_michelScore.push_back( -999. );
      reco_daughter_PFP_trackScore_collection.push_back( -999. );
      reco_daughter_PFP_emScore_collection.push_back( -999. );
      reco_daughter_PFP_michelScore_collection.push_back( -999. );
    }

    if( !evt.isRealData() ){
      //Matching by hits
      DaughterMatchInfo(evt, daughterPFP, clockData);
    }

    try{
      const recob::Track* pandora2Track = pfpUtil.GetPFParticleTrack( *daughterPFP, evt, fPFParticleTag, "pandora2Track" );
      if (fVerbose) std::cout << "pandora2 track: " << pandora2Track << std::endl;

      if( pandora2Track ){
        reco_daughter_allTrack_ID.push_back( pandora2Track->ID() );

        std::cout << "Getting calo for " << pandora2Track->ID() << std::endl;

        auto dummy_caloSCE =
            trackUtil.GetRecoTrackCalorimetry(
                *pandora2Track, evt, "pandora2Track", fPandora2CaloSCE);
        std::cout << dummy_caloSCE.size() << std::endl;
        bool found_calo = false;
        size_t index = 0;
        for ( index = 0; index < dummy_caloSCE.size(); ++index) {
          if (dummy_caloSCE[index].PlaneID().Plane == 2) {
            found_calo = true;
            break; 
          }
        }
        std::cout << index << " " << found_calo << std::endl;

        reco_daughter_allTrack_resRange_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_dEdX_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_dQdX_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_calibrated_dEdX_SCE.push_back( std::vector<double>() );

        if (found_calo) {
          auto dummy_dEdx_SCE = dummy_caloSCE[index].dEdx();
          auto dummy_dQdx_SCE = dummy_caloSCE[index].dQdx();
          auto dummy_Range_SCE = dummy_caloSCE[index].ResidualRange();

          reco_daughter_allTrack_momByRange_alt_proton.push_back( track_p_calc.GetTrackMomentum( dummy_caloSCE[index].Range(), 2212 ) );
          reco_daughter_allTrack_momByRange_alt_muon.push_back(   track_p_calc.GetTrackMomentum( dummy_caloSCE[index].Range(), 13  ) );
          reco_daughter_allTrack_alt_len.push_back(    dummy_caloSCE[index].Range() );

          std::vector<float> cali_dEdX_SCE = calibration_SCE.GetCalibratedCalorimetry(*pandora2Track, evt, "pandora2Track", fPandora2CaloSCE, 2);

          for( size_t j = 0; j < dummy_dEdx_SCE.size(); ++j ){
            reco_daughter_allTrack_resRange_SCE.back().push_back( dummy_Range_SCE[j] );
            reco_daughter_allTrack_dEdX_SCE.back().push_back( dummy_dEdx_SCE[j] );
            reco_daughter_allTrack_dQdX_SCE.back().push_back( dummy_dQdx_SCE[j] );
          }

          for( size_t j = 0; j < cali_dEdX_SCE.size(); ++j ){
            reco_daughter_allTrack_calibrated_dEdX_SCE.back().push_back( cali_dEdX_SCE[j] );
          }

          std::pair<double, int> this_chi2_ndof = trackUtil.Chi2PID(
              reco_daughter_allTrack_calibrated_dEdX_SCE.back(),
              reco_daughter_allTrack_resRange_SCE.back(), templates[2212]);

          reco_daughter_allTrack_Chi2_proton.push_back(this_chi2_ndof.first);
          reco_daughter_allTrack_Chi2_ndof.push_back(this_chi2_ndof.second);
        }
        else {
          reco_daughter_allTrack_momByRange_alt_proton.push_back(-999.);
          reco_daughter_allTrack_momByRange_alt_muon.push_back(-999.);
          reco_daughter_allTrack_alt_len.push_back(-999.);
          reco_daughter_allTrack_Chi2_proton.push_back(-999.);
          reco_daughter_allTrack_Chi2_ndof.push_back(-999);
        }

        //Calorimetry + chi2 for planes 0 and 1
        size_t plane0_index = 0;
        bool found_plane0 = false;
        for ( plane0_index = 0; plane0_index < dummy_caloSCE.size(); ++plane0_index) {
          std::cout << dummy_caloSCE[plane0_index].PlaneID().Plane << std::endl;
          if (dummy_caloSCE[plane0_index].PlaneID().Plane == 0) {
            found_plane0 = true;
            break; 
          }
        }
        size_t plane1_index = 0;
        bool found_plane1 = false;
        for ( plane1_index = 0; plane1_index < dummy_caloSCE.size(); ++plane1_index) {
          std::cout << dummy_caloSCE[plane1_index].PlaneID().Plane << std::endl;
          if (dummy_caloSCE[plane1_index].PlaneID().Plane == 1) {
            found_plane1 = true;
            break; 
          }
        }
        if (fVerbose)
            std::cout << "Plane 0, 1 " << plane0_index << " " <<
                         plane1_index << std::endl;


        reco_daughter_allTrack_calibrated_dEdX_SCE_plane0.push_back(
            std::vector<double>());
        reco_daughter_allTrack_resRange_plane0.push_back(
            std::vector<double>());

        if (found_plane0) {
          auto resRange_plane0 = dummy_caloSCE[plane0_index].ResidualRange();
          for (size_t j = 0; j < resRange_plane0.size(); ++j) {
            reco_daughter_allTrack_resRange_plane0.back().push_back(
                resRange_plane0[j]);
          }

          std::vector<float> dEdX_plane0 = calibration_SCE.GetCalibratedCalorimetry(
              *pandora2Track, evt, "pandora2Track", fPandora2CaloSCE, 0);
          for (size_t j = 0; j < dEdX_plane0.size(); ++j) {
            reco_daughter_allTrack_calibrated_dEdX_SCE_plane0.back().push_back(
                dEdX_plane0[j]);
          }
          std::pair<double, int> plane0_chi2_ndof = trackUtil.Chi2PID(
              reco_daughter_allTrack_calibrated_dEdX_SCE_plane0.back(),
              reco_daughter_allTrack_resRange_plane0.back(), templates[2212]);
          reco_daughter_allTrack_Chi2_proton_plane0.push_back(
              plane0_chi2_ndof.first);
          reco_daughter_allTrack_Chi2_ndof_plane0.push_back(
              plane0_chi2_ndof.second);
        }
        else {
          reco_daughter_allTrack_Chi2_proton_plane0.push_back(
              -999.);
          reco_daughter_allTrack_Chi2_ndof_plane0.push_back(
              -999);
        }


        reco_daughter_allTrack_calibrated_dEdX_SCE_plane1.push_back(
            std::vector<double>());

        reco_daughter_allTrack_resRange_plane1.push_back(
            std::vector<double>());

        if (found_plane1) {
          auto resRange_plane1 = dummy_caloSCE[plane1_index].ResidualRange();
          std::vector<float> dEdX_plane1 = calibration_SCE.GetCalibratedCalorimetry(
              *pandora2Track, evt, "pandora2Track", fPandora2CaloSCE, 1);

          for (size_t j = 0; j < resRange_plane1.size(); ++j) {
            reco_daughter_allTrack_resRange_plane1.back().push_back(
                resRange_plane1[j]);
          }

          for (size_t j = 0; j < dEdX_plane1.size(); ++j) {
            reco_daughter_allTrack_calibrated_dEdX_SCE_plane1.back().push_back(
                dEdX_plane1[j]);
          }

          std::pair<double, int> plane1_chi2_ndof = trackUtil.Chi2PID(
              reco_daughter_allTrack_calibrated_dEdX_SCE_plane1.back(),
              reco_daughter_allTrack_resRange_plane1.back(), templates[2212]);
          reco_daughter_allTrack_Chi2_proton_plane1.push_back(
              plane1_chi2_ndof.first);
          reco_daughter_allTrack_Chi2_ndof_plane1.push_back(
              plane1_chi2_ndof.second);
        }
        else {
          reco_daughter_allTrack_Chi2_proton_plane1.push_back(
              -999.);
          reco_daughter_allTrack_Chi2_ndof_plane1.push_back(
              -999);
        }
        //////////////////////////////////////

        reco_daughter_allTrack_Theta.push_back(  pandora2Track->Theta() );
        reco_daughter_allTrack_Phi.push_back(  pandora2Track->Phi() );

        reco_daughter_allTrack_len.push_back(    pandora2Track->Length() );
        reco_daughter_allTrack_startX.push_back( pandora2Track->Trajectory().Start().X() );
        reco_daughter_allTrack_startY.push_back( pandora2Track->Trajectory().Start().Y() );
        reco_daughter_allTrack_startZ.push_back( pandora2Track->Trajectory().Start().Z() );
        reco_daughter_allTrack_endX.push_back(   pandora2Track->Trajectory().End().X() );
        reco_daughter_allTrack_endY.push_back(   pandora2Track->Trajectory().End().Y() );
        reco_daughter_allTrack_endZ.push_back(   pandora2Track->Trajectory().End().Z() );

        std::cout << "Getting michel" << std::endl;
        //Using new michel tagging
        if (!fSkipMVA) {
          std::pair<double, int> vertex_results =
              trackUtil.GetVertexMichelScore(
                  *pandora2Track, evt, "pandora2Track", fHitTag,
                  0., -500., 500., 0., 500., 0., false,
                  reco_beam_endX, reco_beam_endY, reco_beam_endZ);

          reco_daughter_allTrack_vertex_michel_score.push_back(
              vertex_results.first);
          reco_daughter_allTrack_vertex_nHits.push_back(
              vertex_results.second);
        }
        else {
          reco_daughter_allTrack_vertex_michel_score.push_back(
              -999.);
          reco_daughter_allTrack_vertex_nHits.push_back(
              -999);
        }

        if (fVerbose) std::cout << "pandora2Length " << pandora2Track->Length() << std::endl;
        reco_daughter_allTrack_momByRange_proton.push_back( track_p_calc.GetTrackMomentum( pandora2Track->Length(), 2212 ) );
        reco_daughter_allTrack_momByRange_muon.push_back(   track_p_calc.GetTrackMomentum( pandora2Track->Length(), 13  ) );

        //Match the daughters to a slice
        //First, check whether the start or end of the daughter track are closer
        /*
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
          reco_daughter_allTrack_dR.push_back( dr_end );
        }
        else{
          //min_index = min_start_index;
          reco_daughter_allTrack_dR.push_back( dr_start );
        }
        */
        reco_daughter_allTrack_dR.push_back(-999.);
        reco_daughter_allTrack_to_vertex.push_back(-999.);

      }
      else{
        reco_daughter_allTrack_ID.push_back( -999 );
        reco_daughter_allTrack_resRange_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_dEdX_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_dQdX_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_calibrated_dEdX_SCE.push_back(std::vector<double>());
        reco_daughter_allTrack_Chi2_proton.push_back( -999. );
        reco_daughter_allTrack_Chi2_ndof.push_back( -999 );

        //Calorimetry + chi2 for planes 0 and 1
        reco_daughter_allTrack_calibrated_dEdX_SCE_plane0.push_back(
            std::vector<double>());
        reco_daughter_allTrack_calibrated_dEdX_SCE_plane1.push_back(
            std::vector<double>());
        reco_daughter_allTrack_resRange_plane0.push_back(
            std::vector<double>());
        reco_daughter_allTrack_resRange_plane1.push_back(
            std::vector<double>());

        reco_daughter_allTrack_Chi2_proton_plane0.push_back( -999. );
        reco_daughter_allTrack_Chi2_ndof_plane0.push_back( -999 );
        reco_daughter_allTrack_Chi2_proton_plane1.push_back( -999. );
        reco_daughter_allTrack_Chi2_ndof_plane1.push_back( -999 );
        //////////////////////////////////

        reco_daughter_allTrack_Theta.push_back(-999. );
        reco_daughter_allTrack_Phi.push_back(-999.);
        reco_daughter_allTrack_len.push_back( -999. );
        reco_daughter_allTrack_alt_len.push_back(-999.);
        reco_daughter_allTrack_startX.push_back( -999. );
        reco_daughter_allTrack_startY.push_back( -999. );
        reco_daughter_allTrack_startZ.push_back( -999. );
        reco_daughter_allTrack_endX.push_back(   -999. );
        reco_daughter_allTrack_endY.push_back(   -999. );
        reco_daughter_allTrack_endZ.push_back(   -999. );
        reco_daughter_allTrack_to_vertex.push_back( -999. );
        reco_daughter_allTrack_dR.push_back( -999. );
        reco_daughter_allTrack_momByRange_proton.push_back(-999.);
        reco_daughter_allTrack_momByRange_muon.push_back(-999.);

        reco_daughter_allTrack_momByRange_alt_proton.push_back(-999.);
        reco_daughter_allTrack_momByRange_alt_muon.push_back(-999.);

        reco_daughter_allTrack_vertex_michel_score.push_back(-999.);
        reco_daughter_allTrack_vertex_nHits.push_back(-999);

      }
    }
    catch( const cet::exception &e ){
      MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Track object not found, moving on" << std::endl;
    }


    try{
      const recob::Shower* pandora2Shower = pfpUtil.GetPFParticleShower( *daughterPFP, evt, fPFParticleTag, "pandora2Shower" );
      if (fVerbose) std::cout << "pandora2 shower: " << pandora2Shower << std::endl;

      if( pandora2Shower ){
        reco_daughter_allShower_ID.push_back(     pandora2Shower->ID() );
        reco_daughter_allShower_len.push_back(    pandora2Shower->Length() );
        reco_daughter_allShower_startX.push_back( pandora2Shower->ShowerStart().X() );
        reco_daughter_allShower_startY.push_back( pandora2Shower->ShowerStart().Y() );
        reco_daughter_allShower_startZ.push_back( pandora2Shower->ShowerStart().Z() );

        reco_daughter_allShower_dirX.push_back( pandora2Shower->Direction().X() );
        reco_daughter_allShower_dirY.push_back( pandora2Shower->Direction().Y() );
        reco_daughter_allShower_dirZ.push_back( pandora2Shower->Direction().Z() );


        const std::vector<art::Ptr<recob::Hit>> hits =
            showerUtil.GetRecoShowerArtHits(
                *pandora2Shower, evt, "pandora2Shower");

        art::FindManyP<recob::SpacePoint> spFromHits(hits, evt, fHitTag);
        //double total_shower_energy = 0.;
        //need to get average y
        std::vector<double> x_vec, y_vec, z_vec;
        double total_y = 0.;
        int n_good_y = 0;
        std::vector<art::Ptr<recob::Hit>> good_hits;

        for (size_t iHit = 0; iHit < hits.size(); ++iHit) {
          auto theHit = hits[iHit];
          if (theHit->View() != 2) continue; //skip induction planes

          good_hits.push_back(theHit);

          double shower_hit_x = detProp.ConvertTicksToX(
              theHit->PeakTime(),
              theHit->WireID().Plane,
              theHit->WireID().TPC, 0);

          double shower_hit_z = geom->Wire(theHit->WireID()).GetCenter().Z();

          x_vec.push_back(shower_hit_x);
          z_vec.push_back(shower_hit_z);

          std::vector<art::Ptr<recob::SpacePoint>> sps = spFromHits.at(iHit);
          //std::cout << shower_hit_x << " " << shower_hit_z << " ";
          if (!sps.empty()) {
            y_vec.push_back(sps[0]->XYZ()[1]);
            total_y += y_vec.back();
            ++n_good_y;
            //std::cout << shower_hit_y_vec.back();
          }
          else {
           y_vec.push_back(-999.);
          }
          //std::cout << std::endl;
        }

        if (n_good_y < 1) {
          reco_daughter_allShower_energy.push_back(-999.);
        }
        else {
          double total_shower_energy = 0.;
          for (size_t iHit = 0; iHit < good_hits.size(); ++iHit) {
            auto theHit = good_hits[iHit];
            if (theHit->View() != 2) continue; //skip induction planes

            if (y_vec[iHit] < -100.)
              y_vec[iHit] = total_y / n_good_y;

            total_shower_energy += calibration_SCE.HitToEnergy(
                good_hits[iHit], x_vec[iHit], y_vec[iHit], z_vec[iHit]);
          }
          reco_daughter_allShower_energy.push_back(total_shower_energy);
        }
      }
      else{
        reco_daughter_allShower_ID.push_back(-999);
        reco_daughter_allShower_len.push_back(-999.);
        reco_daughter_allShower_startX.push_back(-999.);
        reco_daughter_allShower_startY.push_back(-999.);
        reco_daughter_allShower_startZ.push_back(-999.);
        reco_daughter_allShower_dirX.push_back(-999.);
        reco_daughter_allShower_dirY.push_back(-999.);
        reco_daughter_allShower_dirZ.push_back(-999.);
        reco_daughter_allShower_energy.push_back(-999.);
      }

    }
    catch( const cet::exception &e ){
      MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Shower object not found, moving on" << std::endl;
    }
  }
}

void pduneana::PDSPAnalyzer::BeamMatchInfo(
    const art::Event & evt, const recob::PFParticle * particle,
    const simb::MCParticle * true_beam_particle,
    detinfo::DetectorClocksData const& clockData) {
  protoana::MCParticleSharedHits beam_match  = truthUtil.GetMCParticleByHits( clockData, *particle, evt, fPFParticleTag, fHitTag );
  if( beam_match.particle ){
    reco_beam_true_byHits_matched = ( beam_match.particle->TrackId() == true_beam_particle->TrackId() );
    reco_beam_true_byHits_PDG = beam_match.particle->PdgCode();
    reco_beam_true_byHits_ID = beam_match.particle->TrackId();

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

    auto list = truthUtil.GetMCParticleListByHits( clockData, *particle, evt, fPFParticleTag, fHitTag );
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

  const simb::MCParticle* trueParticle = 0x0;
  trueParticle = truthUtil.GetMCParticleFromPFParticle(clockData, *particle, evt, fPFParticleTag);
  if( trueParticle ){

    //Check that this is the correct true particle
    if( trueParticle->TrackId() == true_beam_particle->TrackId() ){
      reco_beam_true_byE_matched = true;
    }

    reco_beam_true_byE_PDG = trueParticle->PdgCode();
    reco_beam_true_byE_ID = trueParticle->TrackId();

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
}

void pduneana::PDSPAnalyzer::DaughterMatchInfo(
    const art::Event & evt, const recob::PFParticle * daughterPFP,
    detinfo::DetectorClocksData const& clockData) {

  const sim::ParticleList & plist = pi_serv->ParticleList();
  protoana::MCParticleSharedHits match = truthUtil.GetMCParticleByHits( clockData, *daughterPFP, evt, fPFParticleTag, fHitTag );
  if( match.particle ){

    reco_daughter_PFP_true_byHits_PDG.push_back( match.particle->PdgCode() );
    reco_daughter_PFP_true_byHits_ID.push_back( match.particle->TrackId() );
    reco_daughter_PFP_true_byHits_parID.push_back( match.particle->Mother() );
    reco_daughter_PFP_true_byHits_parPDG.push_back(
      ( (match.particle->Mother() > 0) ? plist[ match.particle->Mother() ]->PdgCode() : 0 )
    );

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

    auto list = truthUtil.GetMCParticleListByHits( clockData, *daughterPFP, evt, fPFParticleTag, fHitTag );
    double total = 0.;
    double matched_hits = 0.;
    for( size_t j = 0; j < list.size(); ++j ){
      if( list[j].particle == match.particle ){
         matched_hits = list[j].nSharedHits + list[j].nSharedDeltaRayHits;
      }
      total += list[j].nSharedHits + list[j].nSharedDeltaRayHits;
    }

    reco_daughter_PFP_true_byHits_purity.push_back( matched_hits / total );

    double totalTruth = truthUtil.GetMCParticleHits( clockData, *match.particle, evt, fHitTag).size();
    double sharedHits = truthUtil.GetSharedHits( clockData, *match.particle, *daughterPFP, evt, fPFParticleTag).size();
    reco_daughter_PFP_true_byHits_completeness.push_back( sharedHits/totalTruth );
  }
  else{
    reco_daughter_PFP_true_byHits_PDG.push_back( -999 );
    reco_daughter_PFP_true_byHits_ID.push_back( -999 );
    reco_daughter_PFP_true_byHits_origin.push_back( -999 );
    reco_daughter_PFP_true_byHits_parID.push_back( -999 );
    reco_daughter_PFP_true_byHits_parPDG.push_back( -999 );
    reco_daughter_PFP_true_byHits_process.push_back( "empty" );
    reco_daughter_PFP_true_byHits_sharedHits.push_back( -999 );
    reco_daughter_PFP_true_byHits_emHits.push_back( -999 );

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
    reco_daughter_PFP_true_byHits_completeness.push_back( -999. );
  }

  //Matching by energy
  const simb::MCParticle* true_daughter_byE = truthUtil.GetMCParticleFromPFParticle(clockData, *daughterPFP, evt, fPFParticleTag);
  if( true_daughter_byE ){
    reco_daughter_PFP_true_byE_PDG.push_back( true_daughter_byE->PdgCode() );
    reco_daughter_PFP_true_byE_len.push_back( true_daughter_byE->Trajectory().TotalLength() );
    double purity = truthUtil.GetPurity( clockData, *daughterPFP, evt, fPFParticleTag);
    double completeness = truthUtil.GetCompleteness( clockData, *daughterPFP, evt, fPFParticleTag, fHitTag );
    reco_daughter_PFP_true_byE_purity.push_back( purity );
    reco_daughter_PFP_true_byE_completeness.push_back( completeness );
  }
  else {
    reco_daughter_PFP_true_byE_PDG.push_back( -999 );
    reco_daughter_PFP_true_byE_len.push_back( -999. );
    reco_daughter_PFP_true_byE_purity.push_back( -999. );
    reco_daughter_PFP_true_byE_completeness.push_back( -999. );
  }

}

void pduneana::PDSPAnalyzer::BeamForcedTrackInfo(
    const art::Event & evt, const recob::PFParticle * particle) {
  try{
    const recob::Track* pandora2Track = pfpUtil.GetPFParticleTrack( *particle, evt, fPFParticleTag, "pandora2Track" );
    if (fVerbose) std::cout << "pandora2 track: " << pandora2Track << std::endl;

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

      auto calo = trackUtil.GetRecoTrackCalorimetry(*pandora2Track, evt, fTrackerTag, fPandora2CaloSCE);
      size_t index = 0;
      bool found_plane = false;
      for (index = 0; index < calo.size(); ++index) {
        if (calo[index].PlaneID().Plane == 2) {
          found_plane = true;
          break; 
        }
      }

      if (found_plane) {
        auto calo_range = calo[index].ResidualRange();
        for( size_t i = 0; i < calo_range.size(); ++i ){
          reco_beam_allTrack_resRange.push_back( calo_range[i] );
        }

        //New Calibration
        std::vector< float > new_dEdX = calibration_SCE.GetCalibratedCalorimetry(*pandora2Track, evt, fTrackerTag, fPandora2CaloSCE, 2);
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
    MF_LOG_WARNING("PDSPAnalyzer") << "beam pandora2Track object not found, moving on" << std::endl;
  }
}

/////////////Below: To get the fitted lines for start/end directions
// define the parameteric line equation
void pduneana::line(double t, double *p,
                                 double &x, double &y, double &z) {
   // a parameteric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
}

// calculate distance line-point
double pduneana::distance2(double x,double y,double z, double *p) {
   // distance line point is D= | (xp-x0) cross  ux |
   // where ux is direction of line and x0 is a point in the line (like t = 0)
   ROOT::Math::XYZVector xp(x,y,z);
   ROOT::Math::XYZVector x0(p[0], p[2], 0.);
   ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1.);
   ROOT::Math::XYZVector u = (x1-x0).Unit();
   double d2 = ((xp-x0).Cross(u)) .Mag2();
   return d2;
}

// function to be minimized
void pduneana::SumDistance2(int &, double *, double & sum,
                                         double * par, int) {
   // the TGraph must be a global variable
   TGraph2D * gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
   assert(gr != 0);
   double * x = gr->GetX();
   double * y = gr->GetY();
   double * z = gr->GetZ();
   int npoints = gr->GetN();
   sum = 0;
   for (int i  = 0; i < npoints; ++i) {
      double d = distance2(x[i],y[i],z[i],par);
      sum += d;
   }
}

TVector3 pduneana::FitLine(const std::vector<TVector3> & input) {
  TGraph2D * gr = new TGraph2D();
  for (size_t i = 0; i < input.size(); ++i) {
    gr->SetPoint(i, input[i].X(), input[i].Y(), input[i].Z());
  }

  TVirtualFitter * min = TVirtualFitter::Fitter(0,4);
  min->SetObjectFit(gr);
  min->SetFCN(SumDistance2);

  double arglist[10];
  
  arglist[0] = -1;
  min->ExecuteCommand("SET PRINT",arglist,1);
  

  double pStart[4] = {1,1,1,1};
  min->SetParameter(0,"x0",pStart[0],0.01,0,0);
  min->SetParameter(1,"Ax",pStart[1],0.01,0,0);
  min->SetParameter(2,"y0",pStart[2],0.01,0,0);
  min->SetParameter(3,"Ay",pStart[3],0.01,0,0);

  arglist[0] = 1000; // number of function calls 
  arglist[1] = 0.001; // tolerance 
  min->ExecuteCommand("MIGRAD", arglist, 2);

  // get fit parameters
  double parFit[4];
  for (int i = 0; i < 4; ++i) {
     parFit[i] = min->GetParameter(i);
  }
  double startX1, startY1, startZ1;
  double startX2, startY2, startZ2;
  line(0, parFit, startX1, startY1, startZ1);
  line(1, parFit, startX2, startY2, startZ2);
  
  TVector3 diff(startX2 - startX1,
                startY2 - startY1,
                startZ2 - startZ1);
  delete gr;
  delete min;
  return diff;
}
DEFINE_ART_MODULE(pduneana::PDSPAnalyzer)
