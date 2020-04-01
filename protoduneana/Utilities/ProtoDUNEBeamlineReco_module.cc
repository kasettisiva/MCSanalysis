////////////////////////////////////////////////////////////////////////
//
// ProtoDUNEBeamlineReco does tracking and reconstruction from the beamline
// monitors
//
// 2018 Jake Calcutt, calcuttj@msu.edu 
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "ProtoDUNETruthUtils.h"

#include "lardataobj/RecoBase/Track.h"

#include "TTree.h"

namespace protoana{
  class ProtoDUNEBeamlineReco;
}

class protoana::ProtoDUNEBeamlineReco : public art::EDAnalyzer{
public:
  
  explicit ProtoDUNEBeamlineReco(fhicl::ParameterSet const& pset);
  virtual ~ProtoDUNEBeamlineReco();
  
  virtual void beginJob() override;
  virtual void endJob() override;

  void analyze(art::Event const & evt) override;
  void reconfigure(fhicl::ParameterSet const& pset);
  void reset();
  
private:
  
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
  bool fDoGlitches;
  int  fReferenceMomentum;

  TTree * fOutTree;

  double tof;
  int event, run;
  int chan;
  double momentum;
  int c0, c1;
  std::vector<int> glitches_p1, glitches_p2, glitches_p3;
  std::vector< int > glitches_h_upstream, glitches_v_upstream, glitches_h_downstream, glitches_v_downstream;
  int nMomenta;
  std::vector<double> momenta;
  std::vector<int> possible_pdg;
  int nTracks;
  std::vector<double> trackEndX, trackEndY, trackEndZ;


  std::vector<short> fibers_h_upstream,   fibers_v_upstream, 
                     fibers_h_downstream, fibers_v_downstream,
                     fibers_p1,           fibers_p2,
                     fibers_p3;

  unsigned long long GenTrigTS;
  bool perfectP;

  int true_PDG;
  double true_X, true_Y, true_Z;
};
  
//-----------------------------------------------------------------------
protoana::ProtoDUNEBeamlineReco::ProtoDUNEBeamlineReco(fhicl::ParameterSet const& pset):
  EDAnalyzer(pset),
  fBeamlineUtils(pset.get<fhicl::ParameterSet>("BeamlineUtils")),
  fDoGlitches( pset.get< bool >( "DoGlitches" ) ),
  fReferenceMomentum( pset.get< int >( "ReferenceMomentum" ) )
{
  this->reconfigure(pset);
}

//-----------------------------------------------------------------------
protoana::ProtoDUNEBeamlineReco::~ProtoDUNEBeamlineReco(){}

//-----------------------------------------------------------------------
void protoana::ProtoDUNEBeamlineReco::beginJob() {
  
  art::ServiceHandle<art::TFileService> tfs;

  fOutTree = tfs->make<TTree>("tree", "");
  fOutTree->Branch("event", &event);
  fOutTree->Branch("run", &run);
  fOutTree->Branch("glitches_p1", &glitches_p1);
  fOutTree->Branch("glitches_p2", &glitches_p2);
  fOutTree->Branch("glitches_p3", &glitches_p3);
  fOutTree->Branch("TOF", &tof);
  fOutTree->Branch("perfectP", &perfectP);
  fOutTree->Branch("Chan", &chan);
  fOutTree->Branch("Momentum", &momentum);
  fOutTree->Branch("C0",&c0);
  fOutTree->Branch("C1",&c1);
  fOutTree->Branch("nMomenta", &nMomenta);
  fOutTree->Branch("Momenta", &momenta);
  fOutTree->Branch("nTracks", &nTracks);
  fOutTree->Branch("trackEndX", &trackEndX);
  fOutTree->Branch("trackEndY", &trackEndY);
  fOutTree->Branch("trackEndZ", &trackEndZ);
  fOutTree->Branch("possible_pdg", &possible_pdg);

  //Tracking Monitors
  fOutTree->Branch("fibers_h_upstream", &fibers_h_upstream);
  fOutTree->Branch("fibers_v_upstream", &fibers_v_upstream);
  fOutTree->Branch("fibers_h_downstream", &fibers_h_downstream);
  fOutTree->Branch("fibers_v_downstream", &fibers_v_downstream);

  fOutTree->Branch("glitches_h_upstream", &glitches_h_upstream);
  fOutTree->Branch("glitches_v_upstream", &glitches_v_upstream);
  fOutTree->Branch("glitches_h_downstream", &glitches_h_downstream);
  fOutTree->Branch("glitches_v_downstream", &glitches_v_downstream);

  //Momentum Monitors
  fOutTree->Branch("fibers_p1", &fibers_p1);
  fOutTree->Branch("fibers_p2", &fibers_p2);
  fOutTree->Branch("fibers_p3", &fibers_p3);

  //Timestamp
  fOutTree->Branch("GenTrigTS", &GenTrigTS);

  //True info
  fOutTree->Branch("true_PDG", &true_PDG);
  fOutTree->Branch("true_X", &true_X);
  fOutTree->Branch("true_Y", &true_Y);
  fOutTree->Branch("true_Z", &true_Z);
}

//-----------------------------------------------------------------------
void protoana::ProtoDUNEBeamlineReco::reconfigure(fhicl::ParameterSet const& pset){
}

//-----------------------------------------------------------------------
void protoana::ProtoDUNEBeamlineReco::analyze(art::Event const & evt){

  reset();

  //std::cout << "Computed TOF "      << fBeamlineUtils.ComputeTOF( kProton, 2.0 ) << std::endl;
  //std::cout << "Computed Momentum " << fBeamlineUtils.ComputeMomentum( kProton, 130. ) << std::endl;

  c0 = -1; 
  c1 = -1;
  momentum = -1.;
  tof = -1.;
  chan = -1;
  GenTrigTS = 0;
  true_PDG = -1;
  true_X = 0.;
  true_Y = 0.;
  true_Z = 0.;

  event = evt.id().event();
  run   = evt.run();

  //Access the Beam Event
  /*
  auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("beamevent");
  
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  if( beamHandle.isValid()){
    art::fill_ptr_vector(beamVec, beamHandle);
  }

  const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); //Should just have one
  */
  const beam::ProtoDUNEBeamEvent & beamEvent = fBeamlineUtils.GetBeamEvent(evt);
  /////////////////////////////////////////////////////////////
  
  
  //Check the quality of the event
  std::cout << "Timing Trigger: " << beamEvent.GetTimingTrigger() << std::endl; 
  std::cout << "Is Matched: "     << beamEvent.CheckIsMatched() << std::endl << std::endl;

  if( !fBeamlineUtils.IsGoodBeamlineTrigger( evt ) ){
    std::cout << "Failed quality check" << std::endl;
    return;
  }

  std::cout << "Passed quality check!" << std::endl << std::endl;  
  /////////////////////////////////////////////////////////////
  

  //Access momentum
  const std::vector< double > & the_momenta = beamEvent.GetRecoBeamMomenta();
  std::cout << "Number of reconstructed momenta: " << the_momenta.size() << std::endl;

  if( the_momenta.size() > 0 ) 
    std::cout << "Measured Momentum: " << the_momenta[0] << std::endl;

  if( the_momenta.size()  == 1)
    momentum = the_momenta[0];

  momenta.insert( momenta.end(), the_momenta.begin(), the_momenta.end() );
  /*
  for( size_t i = 0; i < the_momenta.size(); ++i ){
    momenta.push_back( the_momenta[i] );
  }
  */
  nMomenta = momenta.size();
  ///////////////////////////////////////////////////////////// 


  std::cout << "Current: " << beamEvent.GetMagnetCurrent() << std::endl;

  //Access time of flight
  const std::vector< double > & the_tofs  = beamEvent.GetTOFs();
  const std::vector< int    > & the_chans = beamEvent.GetTOFChans();

  std::cout << "Number of measured TOF: " << the_tofs.size() << std::endl;
  std::cout << "First TOF: "              << beamEvent.GetTOF()         << std::endl;
  std::cout << "First TOF Channel: "      << beamEvent.GetTOFChan()     << std::endl << std::endl;

  std::cout << "All (TOF, Channels): " << std::endl;
  for( size_t i = 0; i < the_tofs.size(); ++i ){
    std::cout << "\t(" << the_tofs[i] << ", " << the_chans[i] << ")" << std::endl;
  }
  std::cout << std::endl;

  if( the_tofs.size() > 0){
    tof = the_tofs[0];
    chan = the_chans[0];
  }
  /////////////////////////////////////////////////////////////
  

  //Access Cerenkov info
  std::cout << "Cerenkov status, pressure:" << std::endl;
  std::cout << "C0: " << beamEvent.GetCKov0Status() << ", " << beamEvent.GetCKov0Pressure() << std::endl;
  std::cout << "C1: " << beamEvent.GetCKov1Status() << ", " << beamEvent.GetCKov1Pressure() << std::endl << std::endl;
  c0 = beamEvent.GetCKov0Status();
  c1 = beamEvent.GetCKov1Status();
  ///////////////////////////////////////////////////////////// 



  //Access PID
  std::vector< int > pids = fBeamlineUtils.GetPID( beamEvent, fReferenceMomentum );

  std::cout << "Possible particles" << std::endl;

  for( size_t i = 0; i < pids.size(); ++i ){ 
    std::cout << pids[i] << std::endl;
  }
  std::cout << std::endl;

  possible_pdg.insert(possible_pdg.end(), pids.begin(), pids.end());

  PossibleParticleCands candidates = fBeamlineUtils.GetPIDCandidates( beamEvent, fReferenceMomentum );
  std::cout << std::left << std::setw(10) << "electron " << candidates.electron << std::endl;
  std::cout << std::left << std::setw(10) << "muon "     << candidates.muon     << std::endl;
  std::cout << std::left << std::setw(10) << "pion "     << candidates.pion     << std::endl;
  std::cout << std::left << std::setw(10) << "kaon "     << candidates.kaon     << std::endl;
  std::cout << std::left << std::setw(10) << "proton "   << candidates.proton   << std::endl << std::endl;

  std::string candidates_string = fBeamlineUtils.GetPIDCandidates( beamEvent, fReferenceMomentum );
  std::cout << candidates_string << std::endl;
  ///////////////////////////////////////////////////////////// 
  
  //Tracking info
  nTracks = beamEvent.GetBeamTracks().size();
  if( nTracks == 1 ){
    std::cout << "X " << beamEvent.GetBeamTracks()[0].Trajectory().End().X() << std::endl;
    std::cout << "Y " << beamEvent.GetBeamTracks()[0].Trajectory().End().Y() << std::endl;
    std::cout << "Z " << beamEvent.GetBeamTracks()[0].Trajectory().End().Z() << std::endl;
  }
  for (size_t i = 0; i < beamEvent.GetBeamTracks().size(); ++i) {
    trackEndX.push_back(beamEvent.GetBeamTracks()[0].Trajectory().End().X());
    trackEndY.push_back(beamEvent.GetBeamTracks()[0].Trajectory().End().Y());
    trackEndZ.push_back(beamEvent.GetBeamTracks()[0].Trajectory().End().Z());
  }
  /////////////////////////////////////////////////////////////

  //Fibers
  std::cout << beamEvent.GetFiberTime( "XBPF022697" ) << std::endl;
  std::cout.precision(20);
  unsigned long test = (unsigned long)(beamEvent.GetT0().first*1e6) + (unsigned long)(beamEvent.GetT0().second/1.e3);
  std::cout << beamEvent.GetT0().first << " " << beamEvent.GetT0().second << std::endl;
  std::cout << test << std::endl;
  
  GenTrigTS = (unsigned long)(beamEvent.GetT0().first*1e6) + (unsigned long)(beamEvent.GetT0().second/1.e3);

  fibers_p1 = beamEvent.GetActiveFibers( "XBPF022697" );  
  fibers_p2 = beamEvent.GetActiveFibers( "XBPF022701" );  
  fibers_p3 = beamEvent.GetActiveFibers( "XBPF022702" );  

  
  


  fibers_h_upstream = beamEvent.GetActiveFibers( "XBPF022707" );  
  fibers_v_upstream = beamEvent.GetActiveFibers( "XBPF022708" );  
  fibers_h_downstream = beamEvent.GetActiveFibers( "XBPF022716" );  
  fibers_v_downstream = beamEvent.GetActiveFibers( "XBPF022717" );  

  if( fDoGlitches ){
    std::cout << "Doing glitches" << std::endl;
    std::array<short,192> glitches = beamEvent.GetFBM( "XBPF022697" ).glitch_mask;
    for( size_t i = 0; i < 192; ++i ){
      if( glitches[i] ) glitches_p1.push_back( i );
    }

    glitches = beamEvent.GetFBM( "XBPF022701" ).glitch_mask;
    for( size_t i = 0; i < 192; ++i ){
      if( glitches[i] ) glitches_p2.push_back( i );
    }

    glitches = beamEvent.GetFBM( "XBPF022702" ).glitch_mask;
    for( size_t i = 0; i < 192; ++i ){
      if( glitches[i] ) glitches_p3.push_back( i );
    }


    glitches = beamEvent.GetFBM( "XBPF022707" ).glitch_mask;
    for( size_t i = 0; i < 192; ++i ){
      if( glitches[i] ) glitches_h_upstream.push_back( i );
    }

    glitches = beamEvent.GetFBM( "XBPF022708" ).glitch_mask;
    for( size_t i = 0; i < 192; ++i ){
      if( glitches[i] ) glitches_v_upstream.push_back( i );
    }

    glitches = beamEvent.GetFBM( "XBPF022716" ).glitch_mask;
    for( size_t i = 0; i < 192; ++i ){
      if( glitches[i] ) glitches_h_downstream.push_back( i );
    }

    glitches = beamEvent.GetFBM( "XBPF022717" ).glitch_mask;
    for( size_t i = 0; i < 192; ++i ){
      if( glitches[i] ) glitches_v_downstream.push_back( i );
    }
  }
  /////////////////////////////////////////////////////////////   
  
  std::cout << "Perfect Momentum?" << fBeamlineUtils.HasPerfectBeamMomentum( evt ) << std::endl;
  perfectP = fBeamlineUtils.HasPerfectBeamMomentum( evt );
  
  if (!evt.isRealData()) {
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>("generator");
    protoana::ProtoDUNETruthUtils truthUtil;
    auto true_beam_particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
    if (!true_beam_particle) {
      std::cout << "No true beam particle" << std::endl;
    }
    else {
      std::cout << "True PDG: " << true_beam_particle->PdgCode() << std::endl;
      true_PDG = true_beam_particle->PdgCode();
      size_t nPoints = true_beam_particle->NumberTrajectoryPoints();
      for (size_t i = 0; i < nPoints; ++i) {
        if (true_beam_particle->Position(i).Z() > 0. &&
            true_beam_particle->Position(i).Z() < 1.) {
          true_X = true_beam_particle->Position(i).X();
          true_Y = true_beam_particle->Position(i).Y();
          true_Z = true_beam_particle->Position(i).Z();
          break;
        }
      }
    }
  }

  fOutTree->Fill();
}

void protoana::ProtoDUNEBeamlineReco::endJob() {}

void protoana::ProtoDUNEBeamlineReco::reset(){
  momenta.clear();

  fibers_p1.clear();
  fibers_p2.clear();
  fibers_p3.clear();

  fibers_h_upstream.clear();
  fibers_v_upstream.clear();
  fibers_h_downstream.clear();
  fibers_v_downstream.clear();
  glitches_p1.clear();
  glitches_p2.clear();
  glitches_p3.clear();

  glitches_h_upstream.clear();
  glitches_v_upstream.clear();
  glitches_h_downstream.clear();
  glitches_v_downstream.clear();

  perfectP = false;
  possible_pdg.clear();

  trackEndX.clear();
  trackEndY.clear();
  trackEndZ.clear();
}
 
DEFINE_ART_MODULE(protoana::ProtoDUNEBeamlineReco)
