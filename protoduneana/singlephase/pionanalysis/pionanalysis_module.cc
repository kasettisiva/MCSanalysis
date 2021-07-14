////////////////////////////////////////////////////////////////////////
// Class:       pionanalysis
// Plugin Type: analyzer (art v2_11_02)
// File:        pionanalysis_module.cc
// Pionanalysis module: Ajib Paudel 
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RawData/RDTimeStamp.h"

#include "lardata/ArtDataHelper/MVAReader.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamCuts.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TVector3.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include <map> 
#include <tuple>

//const int kMaxTracks  = 1000;
//const int kMaxHits = 10000;
double theta12(double x1, double x2, double y1, double y2, double z1, double z2,double x1p, double x2p, double y1p, double y2p, double z1p, double z2p){
  double numer=(x2-x1)*(x2p-x1p)+(y2-y1)*(y2p-y1p)+(z2-z1)*(z2p-z1p);
  double den1=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
  double den2=(x2p-x1p)*(x2p-x1p)+(y2p-y1p)*(y2p-y1p)+(z2p-z1p)*(z2p-z1p);
  return 180/3.14*acos(numer/sqrt(den1*den2));
} 

using namespace std;

namespace protoana {
  class pionanalysis;
}




class protoana::pionanalysis : public art::EDAnalyzer {
public:
  explicit pionanalysis(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  pionanalysis(pionanalysis const &) = delete;
  pionanalysis(pionanalysis &&) = delete;
  pionanalysis & operator = (pionanalysis const &) = delete;
  pionanalysis & operator = (pionanalysis &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void reset();

private:

  //const art::InputTag fSpacePointModuleLabel;
  const art::InputTag fTrackModuleLabel;
  //const art::InputTag fTimeDecoderModuleLabel;

  const art::InputTag fBeamModuleLabel;
  //std::string fBeamModuleLabel;


  TTree *fTree;
  // Run information
  int run;
  int subrun;
  int event;
  //int trigger;
  double evttime;
  int fNactivefembs[6];

  // beam information
  std::vector<double> beamPosx;
  std::vector<double> beamPosy;
  std::vector<double> beamPosz;
  
  std::vector<double> beamDirx;
  std::vector<double> beamDiry;
  std::vector<double> beamDirz;

  std::vector<double> beamMomentum;

  double tof;
  std::vector<double> tofs;
  std::vector<int> ch_tofs;
  short low_pressure_status;
  short high_pressure_status;
  double low_pressure;
  double high_pressure;

  //Beamline utils   
  protoana::ProtoDUNEBeamCuts beam_cuts; 
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
  protoana::ProtoDUNEDataUtils dataUtil;

  // fcl parameters for PFP particles
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fHitsModuleLabel;
  //std::string fGeneratorTag;
  
  //Beam Momentum
  fhicl::ParameterSet fBeamPars;
  //bool fUseCERNCalibSelection;
  bool fVerbose;


  // define parameters for primary tracks
  std::vector<double> primtrk_startx;
  std::vector<double> primtrk_starty;
  std::vector<double> primtrk_startz;

  std::vector<double> primtrk_endx;
  std::vector<double> primtrk_endy;
  std::vector<double> primtrk_endz;

  std::vector<double> primtrk_Dirx;
  std::vector<double> primtrk_Diry;
  std::vector<double> primtrk_Dirz;
  std::vector<double> primtrklen;
  std::vector<double> primtrkID;
  std::vector<int> primtrk_trktag;
  //hit wire and charge info
  std::vector<std::vector<int> > wireno_2;
  std::vector<std::vector<float> > peakTime_2;
  std::vector<std::vector<float> > dq_2;
  std::vector<std::vector<float> > trkhitx2;
  std::vector<std::vector<float> > trkhity2;
  std::vector<std::vector<float> > trkhitz2;

  
  //calo info
  std::vector< std::vector<double> > primtrk_dqdx;
  std::vector< std::vector<double> > primtrk_resrange;
  std::vector< std::vector<double> > primtrk_dedx;
  std::vector<double> primtrk_range;
  std::vector< std::vector<double> > primtrk_hitx;
  std::vector< std::vector<double> > primtrk_hity;
  std::vector< std::vector<double> > primtrk_hitz;
  std::vector< std::vector<double> > primtrk_pitch;

  double cosine_beam_primtrk;
 
  std::vector<int> pdg_code;
  std::vector<int> n_daughter;
  std::vector<int> n_beamparticle;
  std::vector<int> isPrimary;
  std::vector<int> pfp_self;
  //std::vector<int> pfp_parent;
  std::vector<int> pfp_daughter;


  ////Michel tagging
  std::vector< std::vector<int> > Mendhitssecondary;
  std::vector< std::vector<double> > Msecondarystartx;
  std::vector< std::vector<double> > Msecondaryendx;
  std::vector< std::vector<double> > Msecondarystarty;
  std::vector< std::vector<double> > Msecondaryendy;
  std::vector< std::vector<double> > Msecondarystartz;
  std::vector< std::vector<double> > Msecondaryendz;
  std::vector< std::vector<double> > MdQmichel;
  std::vector< std::vector<double> > MdQtrackend;
  std::vector< std::vector<double> > MdQtrackbegin;
  std::vector< std::vector<double> > Mprimsectheta;
  std::vector< std::vector<double> > Mtracklengthsecondary;
  std::vector< std::vector<int> > MtrackID;
  ////CNN for michel tagging
  std::vector< std::vector<double> > trackscore;
  std::vector< std::vector<double> > emscore;
  std::vector< std::vector<double> > michelscore;
  std::vector< std::vector<double> > nonescore;




};


protoana::pionanalysis::pionanalysis(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  //fSpacePointModuleLabel(p.get< art::InputTag >("SpacePointModuleLabel")),
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  beam_cuts(p.get<fhicl::ParameterSet>("BeamCuts")),
  fBeamlineUtils(p.get<fhicl::ParameterSet>("BeamlineUtils")),
  //fBeamModuleLabel(p.get<std::string>("BeamModuleLabel")),
  
  //fTimeDecoderModuleLabel(p.get< art::InputTag >("TimeDecoderModuleLabel")),

  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fHitsModuleLabel(p.get<std::string>("HitsModuleLabel")),
  //fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fBeamPars(p.get<fhicl::ParameterSet>("BeamPars")),
//fUseCERNCalibSelection(p.get<bool>("UseCERNCalibSelection")),
  fVerbose(p.get<bool>("Verbose"))
{
  //if (fSaveTrackInfo == false) fSaveCaloInfo = false;

}

void protoana::pionanalysis::analyze(art::Event const & evt)
{
  //reset containers
  protoana::pionanalysis::reset();  


  //Access the Beam Event ======================================================================//
  auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("beamevent");
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  if( beamHandle.isValid()){
    art::fill_ptr_vector(beamVec, beamHandle);
  }
  const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); //Should just have one
  /////////////////////////////////////////////////////////////
  
  //Check the quality of the event
  std::cout << "Timing Trigger: " << beamEvent.GetTimingTrigger() << std::endl; 
  std::cout << "Is Matched: "     << beamEvent.CheckIsMatched() << std::endl;

  if( !fBeamlineUtils.IsGoodBeamlineTrigger( evt ) ){
    std::cout << "Failed quality check!\n" << std::endl;
    return; //returns, does nothing when !=IsGoodBeamlineTrigger
  }

  std::cout << "Passed quality check!" << std::endl << std::endl;
  /////////////////////////////////////////////////////////////

  //Access momentum
  const std::vector< double > & momenta = beamEvent.GetRecoBeamMomenta();
  std::cout << "Number of reconstructed momenta: " << momenta.size() << std::endl;

  if( momenta.size() > 0 ) 
    std::cout << "Measured Momentum: " << momenta.at(0) << std::endl;
  ///////////////////////////////////////////////////////////// 

  //Access time of flight
  const std::vector< double > & the_tofs  = beamEvent.GetTOFs();
  const std::vector< int    > & the_chans = beamEvent.GetTOFChans();

  std::cout<<"run/subrun/event:"<<evt.run()<<"/"<<evt.subRun()<<"/"<<evt.id().event()<<std::endl;	
  std::cout << "Number of measured TOF: " << the_tofs.size() << std::endl;
  std::cout << "First TOF: "              << beamEvent.GetTOF()         << std::endl;
  std::cout << "First TOF Channel: "      << beamEvent.GetTOFChan()     << std::endl << std::endl;

  std::cout << "All (TOF, Channels): " << std::endl;
  for( size_t i = 0; i < the_tofs.size(); ++i ){
    std::cout << "\t(" << the_tofs[i] << ", " << the_chans[i] << ")" << std::endl;
  }
  std::cout << std::endl;
  /////////////////////////////////////////////////////////////
  
  //Access Cerenkov info
  std::cout << "Cerenkov status, pressure:" << std::endl;
  std::cout << "C0: " << beamEvent.GetCKov0Status() << ", " << beamEvent.GetCKov0Pressure() << std::endl;
  std::cout << "C1: " << beamEvent.GetCKov1Status() << ", " << beamEvent.GetCKov1Pressure() << std::endl << std::endl;

  //if (beamEvent.GetBITrigger() == 1) {
  //ckov0status = beamEvent.GetCKov0Status(); //low_pressure_status
  //ckov1status = beamEvent.GetCKov1Status(); //high_pressure_status

  //ckov0pressure=beamEvent.GetCKov0Pressure();
  //ckov1pressure=beamEvent.GetCKov1Pressure();
  //}

  if(beamEvent.GetCKov0Pressure() < beamEvent.GetCKov1Pressure() ){
    high_pressure_status = beamEvent.GetCKov1Status();
    low_pressure_status = beamEvent.GetCKov0Status();

    high_pressure=beamEvent.GetCKov1Pressure();
    low_pressure=beamEvent.GetCKov0Pressure(); 
  }
  else{
    high_pressure_status = beamEvent.GetCKov0Status();
    low_pressure_status = beamEvent.GetCKov1Status();

    high_pressure=beamEvent.GetCKov0Pressure();
    low_pressure=beamEvent.GetCKov1Pressure(); 
  }

  ///////////////////////////////////////////////////////////// 
  
  double nom_beam_mon=fBeamPars.get<double>("Momentum");
  std::cout<<"Selected nominal beam momentum is: "<<nom_beam_mon<<" GeV/c"<<std::endl;

  //Access PID
  //std::vector< int > pids = fBeamlineUtils.GetPID( beamEvent, nom_beam_mon);
  //for( size_t i = 0; i < pids.size(); ++i ){
  //std::cout << "pid["<<i<<"]: "<< pids[i] << std::endl;
  //}

  PossibleParticleCands candidates = fBeamlineUtils.GetPIDCandidates( beamEvent, nom_beam_mon);
  std::cout << "electron " << candidates.electron << std::endl;
  std::cout << "muon "     << candidates.muon     << std::endl;
  std::cout << "pion "     << candidates.pion     << std::endl;
  std::cout << "kaon "     << candidates.kaon     << std::endl;
  std::cout << "proton "   << candidates.proton   << std::endl;

  std::cout << std::endl;
  //if (candidates.proton==0) {
  //std::cout<<"no proton!!!!!!"<<std::endl;
  //}
  /////////////////////////////////////////////////////////////

  //Manual Event Selection -----------------------------------------------------------------------------------------------------------------------------------------------//
  /*bool IsProton=false;
  //old Tof: (!fUseCERNCalibSelection&&tof>170.)
  if (nom_beam_mon==1.) { 
  if ((beamEvent.GetTOF()>110.&&beamEvent.GetTOF()<160.)&&low_pressure_status==0) {
  IsProton=true;
  std::cout<<"-->> This is a proton candidate..."<<std::endl;
  }
  }
  */

  //HY::For Calo info
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
  else return;
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryTag);
  art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle, evt, "pandoraTrack");
  // std::vector<art::Ptr<recob::Hit>> pfpHits;//pfpHits definition
  art::Handle< std::vector<recob::Hit> > hitListHandle; // to get information about the hits
  std::vector<art::Ptr<recob::Hit>> hitlist;
  if(evt.getByLabel(fHitsModuleLabel, hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits
  // Implementation of required member function here.
  art::FindManyP<recob::Track> thass(hitListHandle, evt, fTrackModuleLabel); //to associate hit just trying
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  art::Timestamp ts = evt.time();
  //std::cout<<ts.timeHigh()<<" "<<ts.timeLow()<<std::endl;
  if (ts.timeHigh() == 0){
    //TTimeStamp tts(ts.timeLow());
    //evttime = tts.AsDouble();
    TTimeStamp ts2(ts.timeLow());
    evttime = ts2.AsDouble();
  }
  else{
    //TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    //evttime = tts.AsDouble();
    TTimeStamp ts2(ts.timeHigh(), ts.timeLow());
    evttime = ts2.AsDouble();
  }

  // Get number of active fembs
  if(!evt.isRealData()){
    for(int ik=0; ik<6; ik++)
      fNactivefembs[ik]=20;
  }
  else{
    for(int ik=0; ik<6; ik++)
      fNactivefembs[ik]=dataUtil.GetNActiveFembsForAPA(evt,ik);
  }
  //CNN
  anab::MVAReader<recob::Hit,4> hitResults(evt, "emtrkmichelid:emtrkmichel" );



  if (beamVec.size()&&candidates.pion==1) { //if beam pion
    if (beamEvent.GetTimingTrigger()==12) { //get beam timing trigger
      if (beamEvent.CheckIsMatched()){ //if CheckIsMatched
        //Get TOF info
        if (beamEvent.GetTOFChan() != -1) {//if TOFChan == -1, then there was not a successful match, if it's 0, 1, 2, or 3, then there was a good match corresponding to the different pair-wise combinations of the upstream and downstream channels
          tof = beamEvent.GetTOF();
  	  for( size_t ii = 0; ii < the_tofs.size(); ++ii ){
	    tofs.push_back(the_tofs[ii]);
	    ch_tofs.push_back(the_chans[ii]);
          } 
        }
        //Get beam particle trajectory info
        //auto & tracks = beaminfo[0]->GetBeamTracks();
        auto & tracks = beamEvent.GetBeamTracks();
	std::cout<<"ToF:"<<tof<<" [ns]"<<std::endl;
	std::cout<<"beam trk size:"<<tracks.size()<<std::endl;
        for (size_t i = 0; i<tracks.size(); ++i){
          beamPosx.push_back(tracks[i].End().X());
          beamPosy.push_back(tracks[i].End().Y());
          beamPosz.push_back(tracks[i].End().Z());
          beamDirx.push_back(tracks[i].StartDirection().X());
          beamDiry.push_back(tracks[i].StartDirection().Y());
          beamDirz.push_back(tracks[i].StartDirection().Z());

	  std::cout<<"run/subrun/evt:"<<run<<"/"<<subrun<<"/"<<event<<std::endl;	
	  std::cout<<"beamPosx/beamPosy/beamPosz:"<<tracks[i].End().X()<<"/"<<tracks[i].End().Y()<<"/"<<tracks[i].End().Z()<<std::endl;
	  std::cout<<"beamDirx/beamDiry/beamDirz:"<<tracks[i].StartDirection().X()<<"/"<<tracks[i].StartDirection().Y()<<"/"<<tracks[i].StartDirection().Z()<<std::endl;
	  //std::cout<<"beamDirx^2+beamDiry^2+beamDirz^2:"<<tracks[i].StartDirection().X()*tracks[i].StartDirection().X()+tracks[i].StartDirection().Y()*tracks[i].StartDirection().Y()+tracks[i].StartDirection().Z()*tracks[i].StartDirection().Z()<<std::endl;
        }
        //Get reconstructed beam momentum info
        //auto & beammom = beaminfo[0]->GetRecoBeamMomenta();
        //auto & beammom = beamEvent.GetRecoBeamMomenta();
        //auto & beammom = beamEvent.GetRecoBeamMomenta();
	//std::cout<<"==============================================================="<<std::endl;
	std::cout<<"beam mom size:"<<momenta.size()<<std::endl;
        for (size_t i = 0; i<momenta.size(); ++i){
          beamMomentum.push_back(momenta[i]);
	  std::cout<<"beam mom["<<i<<"]:"<<momenta[i]<<" [GeV]"<<std::endl;
        }
	//std::cout<<"==============================================================="<<std::endl;

        //put the track parameters in the cryo here ----------------------------------------------------------------------------//
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

	std::cout<<"\n*******************************************************"<<std::endl;
	std::cout<<"Moving on to the PFParticle section..."<<std::endl;	

  	// Get the PFParticle utility
  	protoana::ProtoDUNEPFParticleUtils pfpUtil;

        // Get the track utility
        protoana::ProtoDUNETrackUtils trackUtil;

  	// Get all of the PFParticles, by default from the "pandora" product
  	auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  	// We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility 
  	// to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These	
  	// are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those 
  	// PFParticles considered as primary
  	
        /// Use the pandora metadata to tell us if this is a beam particle or not
        //bool isBeamParticle=dataUtil.IsBeamParticle(fPFParticleTag, evt, );
        //bool IsBeamParticle(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;
        

  	//std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  	auto beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
        //HY:: From the new version, prepend recob::PFParticle* with "const"
  	if(beamParticles.size() == 0){
	  std::cerr << "We found no beam particles for this event... moving on" << std::endl;
	  return;
  	}

        //unsigned int nPrimPFPartices=pfpUtil.GetNumberPrimaryPFParticle(evt,fPFParticleTag); //count all the particles
	//std::cout<<"--> Number of all prim. PFPartices:"<<nPrimPFPartices<<std::endl;

	//int tmp_counter=0;
	n_beamparticle.push_back(beamParticles.size());
        std::cout<<"we have "<<beamParticles.size()<<" beam particle(s)"<<std::endl;
  	for(const recob::PFParticle* particle : beamParticles){
	  int nTrack=0;
	  int nShower=0;

	  // "particle" is the pointer to our beam particle. The recob::Track or recob::Shower object
	  // of this particle might be more helpful. These return null pointers if not track-like / shower-like
	  const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
	  //  art::Ptr<recob::Track> *track = thisTrack;
	  // art::Ptr<recob::Track> thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
	  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
	  auto recoTracks = evt.getValidHandle<std::vector<recob::Track>>(fTrackerTag);
	  art::FindManyP<recob::Hit> findHitsFromTracks(recoTracks,evt,fTrackerTag);
	  /////end of Michel hits loop
	  double beamstx=-30;
	  double beamendx=-30;
	  double beamsty=420;
	  double beamendy=420;
	  double beamstz=30;
	  double beamendz=100;    

	  if(thisTrack != 0x0){
	    if(!beam_cuts.IsBeamlike(*thisTrack, evt, "1")) return;
	    beamstx=thisTrack->Start().X();
	    beamsty=thisTrack->Start().Y();
	    beamstz=thisTrack->Start().Z();
	    beamendx=thisTrack->End().X();
	    beamendy=thisTrack->End().Y();
	    beamendz=thisTrack->End().Z();
	    std::cout<<"beamstx "<<beamstx<<std::endl;
	    //////Michel tagging here
	    std::vector<double> secondarystartx1;
	    std::vector<double> secondarystarty1;
	    std::vector<double> secondarystartz1;
	    std::vector<double> secondaryendx1;
	    std::vector<double> secondaryendy1;
	    std::vector<double> secondaryendz1;
	    std::vector<double> dQmichel1;
	    std::vector<double> dQtrackbegin1;
	    std::vector<double> dQtrackend1;
	    std::vector<double> tracklengthsecondary1;
	    std::vector<double> primsectheta1; 
	    std::vector<int> endhitssecondary1,trackID1;
	    std::vector<double> trks1;
	    std::vector<double> ems1;
	    std::vector<double> michels1;
	    std::vector<double> nones1;
	    size_t NTracks = tracklist.size();
	    for(size_t i=0;i<NTracks;i++){
	      art::Ptr<recob::Track> ptrack(trackListHandle, i);
	      const recob::Track& track = *ptrack;
	      auto pos = track.Vertex();
	      auto end = track.End();
	      int counter1=0;
	      double startx=pos.X();   
	      double starty=pos.Y();
	      double startz=pos.Z();
	      double endx=end.X();
	      double endy=end.Y();
	      double endz=end.Z();
	      if(track.Length()<5) continue;
	      //  if(TMath::Max(endy,starty)>520 || TMath::Min(endy, starty)<150 || TMath::Max(startx, endx)>20||TMath::Min(startx,endx)<-300||TMath::Max(startz,endz)<230) continue;
	      if(TMath::Max(endy,starty)>500 || TMath::Min(endy, starty)<200 || TMath::Max(startx, endx)>0||TMath::Min(startx,endx)<-200||TMath::Max(startz,endz)<30) continue;
	      //std::cout<<"event trackid "<<event<<" "<<track.ID()<<std::endl;
	      std::vector<int> wirenos;
	      std::vector<float> peakts,dqbuff1;
	      std::vector<float> dQstart,dQend;
	      std::vector<double> micheldq;
	      wirenos.clear();peakts.clear();dqbuff1.clear();
	      float peaktime=-1;
	      int wireno=-99999;
	      int tpcno=-1;
	      float zlast0=-99999;
	      float zlast=-99999;
	      std::vector<std::tuple<double,double,double,double,int,double>> buff_ZYXTWQ;
	      buff_ZYXTWQ.clear();
	      double thetavalue=theta12(beamstx,beamendx,beamsty,beamendy,beamstz,beamendz,startx,endx,starty,endy,startz,endz);
	      if(fmthm.isValid()){
		auto vhit=fmthm.at(i);
		auto vmeta=fmthm.data(i);
		for (size_t ii = 0; ii<vhit.size(); ++ii){ //loop over all meta data hit
		  bool fBadhit = false;
		  if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max())){
		    fBadhit = true;
		    //cout<<"fBadHit"<<fBadhit<<endl;
		    continue;
		  }
		  if (vmeta[ii]->Index()>=tracklist[i]->NumberTrajectoryPoints()){
		    throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<tracklist[i]->NumberTrajectoryPoints()<<" for track index "<<i<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
		  }
		  if (!tracklist[i]->HasValidPoint(vmeta[ii]->Index())){
		    fBadhit = true;
		    // cout<<"had valid point "<<fBadhit<<endl;
		    continue;
		  }
        
		  auto loc = tracklist[i]->LocationAtPoint(vmeta[ii]->Index());
		  if (fBadhit) continue; //HY::If BAD hit, skip this hit and go next
		  if (loc.Z()<-100) continue; //hit not on track
		  if(vhit[ii]->WireID().Plane==2){
		    buff_ZYXTWQ.push_back(std::make_tuple(loc.Z(),loc.Y(),loc.X(),vhit[ii]->PeakTime(),vhit[ii]->WireID().Wire,vhit[ii]->Integral()));
		    wirenos.push_back(vhit[ii]->WireID().Wire);
		    peakts.push_back(vhit[ii]->PeakTime());
		    zlast=loc.Z();
		    if(zlast>zlast0){
		      zlast0=zlast;
		      wireno=vhit[ii]->WireID().Wire;
		      peaktime=vhit[ii]->PeakTime();
		      tpcno=vhit[ii]->WireID().TPC;
		    }        
		  }//planenum 2
		}//loop over vhit
	      }//fmthm valid
	      //save start and end point of each track
	      //taking care of flipped start and end point
	      if(endz<startz){
		startx=end.X();   
		starty=end.Y();
		startz=end.Z();
		endx=pos.X();
		endy=pos.Y();
		endz=pos.Z();
	      }
	      double trk_score=0.0;
	      double em_score=0;
	      double michel_score=0;
	      double none_score=0;
	      for(size_t hitl=0;hitl<hitlist.size();hitl++){
		std::array<float,4> cnn_out=hitResults.getOutput(hitlist[hitl]);
		auto & tracks = thass.at(hitlist[hitl].key());
		// if (!tracks.empty() && tracks[0].key()!=ptrack.key() && tracklist[tracks[0].key()]->Length()>25) continue;
		if (!tracks.empty() && tracks[0].key()!=ptrack.key() && tracklist[tracks[0].key()]->Length()>25) continue;
		bool test=true;
		float peakth1=hitlist[hitl]->PeakTime();
		int wireh1=hitlist[hitl]->WireID().Wire;
		for(size_t m=0;m<wirenos.size();m++){
		  if(wireh1==wirenos[m] && peakth1==peakts[m]){
		    test=false;
		    break;
		  }
		}
		if(!test) continue;
		int planeid=hitlist[hitl]->WireID().Plane;
		int tpcid=hitlist[hitl]->WireID().TPC;
		if(abs(wireh1-wireno)<20 && abs(peakth1-peaktime)<150 && planeid==2 && tpcid==tpcno){
		  counter1++;
		  micheldq.push_back(hitlist[hitl]->Integral());
		  micheldq.push_back(hitlist[hitl]->Integral());
		  trk_score+=cnn_out[hitResults.getIndex("track")];
		  em_score+=cnn_out[hitResults.getIndex("em")];
		  michel_score+=cnn_out[hitResults.getIndex("michel")];
		  none_score+=cnn_out[hitResults.getIndex("none")];
		}
	      }//hitlist loop
	      if(buff_ZYXTWQ.size()<10) continue;
	      sort(buff_ZYXTWQ.begin(),buff_ZYXTWQ.end());
	      dQstart.clear(); dQend.clear();
	      int qi11=buff_ZYXTWQ.size();
	      for(int qi=5;qi<TMath::Min(15,qi11);qi++){
		dQstart.push_back(std::get<5>(buff_ZYXTWQ[qi]));
	      }
	     
	      for(int qi=qi11-5;qi<qi11;qi++){
		dQend.push_back(std::get<5>(buff_ZYXTWQ[qi]));
	      }
	      secondarystartx1.push_back(startx);
	      secondarystarty1.push_back(starty);
	      secondarystartz1.push_back(startz);
	      secondaryendx1.push_back(endx);
	      secondaryendy1.push_back(endy);
	      secondaryendz1.push_back(endz);
	      endhitssecondary1.push_back(counter1);
	      tracklengthsecondary1.push_back(track.Length());
	      trackID1.push_back(track.ID());
	      dQmichel1.push_back(TMath::Median(micheldq.size(),&micheldq[0]));
	      dQtrackbegin1.push_back(TMath::Median(dQstart.size(),&dQstart[0]));
	      dQtrackend1.push_back(TMath::Median(dQend.size(),&dQend[0]));
	      primsectheta1.push_back(thetavalue);
	      trks1.push_back(trk_score);
	      ems1.push_back(em_score);
	      michels1.push_back(michel_score);
	      nones1.push_back(none_score);
	    }//Ntracks
	    Msecondarystartx.push_back(secondarystartx1);
	    Msecondarystarty.push_back(secondarystarty1);
	    Msecondarystartz.push_back(secondarystartz1);
	    Msecondaryendx.push_back(secondaryendx1);
	    Msecondaryendy.push_back(secondaryendy1);
	    Msecondaryendz.push_back(secondaryendz1);
	    Mendhitssecondary.push_back(endhitssecondary1);
	    Mtracklengthsecondary.push_back(tracklengthsecondary1);
	    MtrackID.push_back(trackID1);
	    MdQmichel.push_back(dQmichel1);
	    MdQtrackbegin.push_back(dQtrackbegin1);
	    MdQtrackend.push_back(dQtrackend1);
	    Mprimsectheta.push_back(primsectheta1);
	    trackscore.push_back(trks1);
	    emscore.push_back(ems1);
	    michelscore.push_back(michels1);
	    nonescore.push_back(nones1);
	    //std::cout<<"clearing trees now "<<std::endl;
	    endhitssecondary1.clear();
	    secondarystartx1.clear();
	    secondaryendx1.clear();
	    secondarystarty1.clear();
	    secondaryendy1.clear();
	    secondarystartz1.clear();
	    secondaryendz1.clear();
	    dQmichel1.clear();
	    primsectheta1.clear();
	    dQtrackbegin1.clear();
	    dQtrackend1.clear();
	    tracklengthsecondary1.clear();
	    trks1.clear();
	    ems1.clear();
	    michels1.clear();
	    nones1.clear();
	  }
	  ///////////////End of Michel checking
	  /////defining secondary hits parameters

	  if(thisTrack != 0x0) {
	    if(!beam_cuts.IsBeamlike(*thisTrack, evt, "1")) return;
	    std::cout << "Beam particle is track-like" << std::endl;
	    nTrack++;
	    primtrk_trktag.push_back(1);
	   
	    //Adding hit coordinates info using space points
	    // pfpHits = findHitsFromTracks.at(thisTrack->ID());  
	    // art::FindManyP<recob::SpacePoint> spFromHits(pfpHits, evt, fHitsModuleLabel);
	    //hits and calorimetry loop
	    ////End of Hit Meta 
	    //HY::Get the Calorimetry(s) from track
	    std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
	    std::vector<double> tmp_primtrk_dqdx;	
	    std::vector<double> tmp_primtrk_resrange;	
	    std::vector<double> tmp_primtrk_dedx;	
	    std::vector<double> tmp_primtrk_hitx;	
	    std::vector<double> tmp_primtrk_hity;	
	    std::vector<double> tmp_primtrk_hitz;
	    std::vector<double> tmp_primtrk_pitch;
	    for (auto & calo : calovector){
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
		  //std::cout<<"dqdx="<<calo.dQdx()[ihit]<<"; resrange="<<calo.ResidualRange()[ihit]<<std::endl;
		  //std::cout<<"(X,Y,Z)="<<"("<<calo.XYZ()[ihit].X()<<","<<calo.XYZ()[ihit].Y()<<","<<calo.XYZ()[ihit].Z()<<")"<<std::endl;
		  //std::cout<<"(X,Y,Z)="<<"("<<primtrk_pos.X()<<","<<primtrk_pos.Y()<<","<<primtrk_pos.Z()<<")"<<std::endl;
		  //std::cout<<"(X,Y,Z)="<<"("<<tmp_primtrk_hitx[ihit]<<","<<tmp_primtrk_hity[ihit]<<","<<tmp_primtrk_hitz[ihit]<<")"<<std::endl;
		} //loop over hits
	      } //only collection plane
	    }
	    //primtrk_dqdx->push_back(tmp_primtrk_dqdx);
	    //primtrk_resrange->push_back(tmp_primtrk_resrange);
	    //primtrk_dedx->push_back(tmp_primtrk_dedx);
	    if (tmp_primtrk_hitz.size()!=0) { //prevent the zero vectors being push_back
	      primtrk_dqdx.push_back(tmp_primtrk_dqdx);
	      primtrk_resrange.push_back(tmp_primtrk_resrange);
	      primtrk_dedx.push_back(tmp_primtrk_dedx);
	      primtrk_hitx.push_back(tmp_primtrk_hitx);
	      primtrk_hity.push_back(tmp_primtrk_hity);
	      primtrk_hitz.push_back(tmp_primtrk_hitz);
	      primtrk_pitch.push_back(tmp_primtrk_pitch);
	    } //prevent the zero vectors being push_back

	    tmp_primtrk_dqdx.clear();
	    tmp_primtrk_resrange.clear();
	    tmp_primtrk_dedx.clear();
	    tmp_primtrk_hitx.clear();
	    tmp_primtrk_hity.clear();
	    tmp_primtrk_hitz.clear();
	    tmp_primtrk_pitch.clear();


	    //HY::Here comes the start/end position of primary track
	    // Find the particle vertex. We need the tracker tag here because we need to do a bit of
	    // additional work if the PFParticle is track-like to find the vertex. 
	    const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
	    const TVector3 vtx_end(thisTrack->Trajectory().End().X(), thisTrack->Trajectory().End().Y(), thisTrack->Trajectory().End().Z());
	    primtrk_startx.push_back(vtx.X());
	    primtrk_starty.push_back(vtx.Y());
	    primtrk_startz.push_back(vtx.Z());
	    std::cout<<"vtx_X:"<<vtx.X()<<" ; vtx_Y:"<<vtx.Y()<<" ; vtx_Z:"<<vtx.Z()<<std::endl;
                
	    primtrk_endx.push_back(vtx_end.X());
	    primtrk_endy.push_back(vtx_end.Y());
	    primtrk_endz.push_back(vtx_end.Z());
	    std::cout<<"vtx_end_X:"<<vtx_end.X()<<" ; vtx_end_Y:"<<vtx_end.Y()<<" ; vtx_end_Z:"<<vtx_end.Z()<<std::endl;

	    if (vtx.Z()==vtx_end.Z()) { //warning message if Z_start=Z_end
	      std::cout<<"WARNING!! StartZ and EndZ are the same!!"<<std::endl;	
	    }



	  }
	  if(thisShower != 0x0) { 
	    std::cout << "Beam particle is shower-like" << std::endl;
	    nShower++;
	    primtrk_trktag.push_back(-1);
	  }



	  //HY Add
	  pdg_code.push_back(particle->PdgCode());
	  n_daughter.push_back(particle->NumDaughters());
	  isPrimary.push_back(particle->IsPrimary());
	  pfp_self.push_back(particle->Self());
	  //pfp_parent.push_back(particle->Parent());
 
	  std::cout<<"pdg code:"<<particle->PdgCode()<<std::endl;
	  std::cout<<"IsPrimary:"<<particle->IsPrimary()<<std::endl;
	  std::cout<<"NumDaughters:"<<particle->NumDaughters()<<std::endl;
	  std::cout<<"Self:"<<particle->Self()<<std::endl;	
	  std::cout<<"Parent:"<<particle->Parent()<<std::endl;

	  if ((particle->NumDaughters())>0) {
	    for (int ii=0; ii<(particle->NumDaughters());++ii) {
	      std::cout<<"Daughter["<<ii<<"]:"<<particle->Daughter(ii)<<std::endl;
	      pfp_daughter.push_back(particle->Daughter(ii));
	    }
	  }
	  else {
	    pfp_daughter.push_back(-99);
	  }


	  //int pdg_code=pfpUtil.PdgCode(*particle,evt,fPFParticleTag,fShowerTag);
	  //std::cout<<"IsDaughters:"<<particle->Daughters()<<std::endl;	
	  //HY::Add parameters for the primary tracks here ---------------//
		
	  // Get track direction
	  if (thisTrack) {
	    //pdg_code=thisTrack->PdgCode(); 
	    //std::cout<<"pdg code:"<<pdg_code<<std::endl;
	
	    auto trackdir = thisTrack->StartDirection();
	    std::cout<<"run/subrun/event:"<<run<<"/"<<subrun<<"/"<<event<<std::endl;	
	    std::cout<<"trkDirx/trkDiry/trkDirz:"<<trackdir.X()<<"/"<<trackdir.Y()<<"/"<<trackdir.Z()<<std::endl;
	    primtrk_Dirx.push_back(trackdir.X());
	    primtrk_Diry.push_back(trackdir.Y());
	    primtrk_Dirz.push_back(trackdir.Z());

	    primtrklen.push_back(thisTrack->Length()); //track length
	    std::cout<<"trk length: "<<thisTrack->Length()<<" [cm]"<<std::endl;
	    primtrkID.push_back(thisTrack->ID());
	    std::cout<<"trk ID: "<<thisTrack->ID()<<""<<std::endl; //HY::Fix me::trk ID seems wrong 
	    //std::cout<<"trk ID: "<<thisTrack->TrackId()<<std::endl;
		    
	    //std::cout<<"trkDirx^2+trkDiry^2+trkDirz^2:"<<trackdir.X()*trackdir.X()+trackdir.Y()*trackdir.Y()+trackdir.Z()*trackdir.Z()<<std::endl;
	    //int nn=tracks.size()-1;
	    //cosine_beam_primtrk=tracks[nn].StartDirection().X()*trackdir.X()+tracks[nn].StartDirection().Y()*trackdir.Y()+tracks[nn].StartDirection().Z()*trackdir.Z();
	    if (tracks.size()){
	      cosine_beam_primtrk=tracks[0].StartDirection().X()*trackdir.X()+tracks[0].StartDirection().Y()*trackdir.Y()+tracks[0].StartDirection().Z()*trackdir.Z();
	    }
	    // fill a histogram of trackdir.X()*beamdir.X() + .....
	    // try to get calorimetry info of this track

	  }
	  // Now we can look for the interaction point of the particle if one exists, i.e where the particle
	  // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
	  // check this by making sure we get a sensible position
	  const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);


	  // Let's get the daughter PFParticles... we can do this simply without the utility
	  for(const int daughterID : particle->Daughters()){
	    // Daughter ID is the element of the original recoParticle vector
	    const recob::PFParticle *daughterParticle = &(recoParticles->at(daughterID));
	    std::cout << "Daughter " << daughterID << " has " << daughterParticle->NumDaughters() << " daughters" << std::endl;
	  }
 
	  // For actually studying the objects, it is easier to have the daughters in their track and shower forms.
	  // We can use the utility to get a vector of track-like and a vector of shower-like daughters
	  const std::vector<const recob::Track*> trackDaughters = pfpUtil.GetPFParticleDaughterTracks(*particle,evt,fPFParticleTag,fTrackerTag);  
	  const std::vector<const recob::Shower*> showerDaughters = pfpUtil.GetPFParticleDaughterShowers(*particle,evt,fPFParticleTag,fShowerTag);  
	  std::cout << "Beam particle has " << trackDaughters.size() << " track-like daughters and " << showerDaughters.size() << " shower-like daughters." << std::endl;

	  std::cout<<"# total Tracks:"<<nTrack<<std::endl;
	  std::cout<<"# total Showers:"<<nShower<<std::endl;
	  //tmp_counter++;
  	}
	//std::cout<<"tmp_counter:"<<tmp_counter<<"\n\n\n"<<std::endl;
	//'tmp_counter' should be the same as 'nBeamP'
	std::cout<<"*******************************************************"<<std::endl;

        //put the track parameters in the cryo here ----------------------------------------------------------------------------//


      } //if CheckIsMatched
    } //get beam timing trigger

  } //if beam pion

  fTree->Fill();
}

void protoana::pionanalysis::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("beamana","beam analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  //fTree->Branch("trigger",&trigger,"trigger/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("Nactivefembs",&fNactivefembs,"Nactivefembs[5]/I");

  fTree->Branch("beamPosx",&beamPosx);
  fTree->Branch("beamPosy",&beamPosy);
  fTree->Branch("beamPosz",&beamPosz);
  fTree->Branch("beamDirx",&beamDirx);
  fTree->Branch("beamDiry",&beamDiry);
  fTree->Branch("beamDirz",&beamDirz);
  fTree->Branch("beamMomentum",&beamMomentum);
  fTree->Branch("tof", &tof, "tof/D");
  fTree->Branch("tofs",&tofs);
  fTree->Branch("ch_tofs",&ch_tofs);
  fTree->Branch("low_pressure_status", &low_pressure_status, "low_pressure_status/S");
  fTree->Branch("high_pressure_status", &high_pressure_status, "high_pressure_status/S");
  fTree->Branch("low_pressure", &low_pressure, "low_pressure/D");
  fTree->Branch("high_pressure", &high_pressure, "high_pressure/D");

  fTree->Branch("cosine_beam_primtrk", &cosine_beam_primtrk, "cosine_beam_primtrk/D");
  fTree->Branch("primtrk_startx",&primtrk_startx);
  fTree->Branch("primtrk_starty",&primtrk_starty);
  fTree->Branch("primtrk_startz",&primtrk_startz);
  //hit wire info
  fTree->Branch("wireno_2",&wireno_2);
  fTree->Branch("peakTime_2",&peakTime_2);
  fTree->Branch("dq_2",&dq_2);
  fTree->Branch("trkhitx2",&trkhitx2);
  fTree->Branch("trkhity2",&trkhity2);
  fTree->Branch("trkhitz2",&trkhitz2);
  //hit wire info ends



  fTree->Branch("primtrk_endx",&primtrk_endx);
  fTree->Branch("primtrk_endy",&primtrk_endy);
  fTree->Branch("primtrk_endz",&primtrk_endz);

  fTree->Branch("primtrk_Dirx",&primtrk_Dirx);
  fTree->Branch("primtrk_Diry",&primtrk_Diry);
  fTree->Branch("primtrk_Dirz",&primtrk_Dirz);
  fTree->Branch("primtrklen",&primtrklen);
  fTree->Branch("primtrkID",&primtrkID);

  fTree->Branch("primtrk_dqdx",&primtrk_dqdx);
  fTree->Branch("primtrk_dedx",&primtrk_dedx);
  fTree->Branch("primtrk_resrange",&primtrk_resrange);
  fTree->Branch("primtrk_range",&primtrk_range);
  fTree->Branch("primtrk_hitx",&primtrk_hitx);
  fTree->Branch("primtrk_hity",&primtrk_hity);
  fTree->Branch("primtrk_hitz",&primtrk_hitz);
  fTree->Branch("primtrk_pitch",&primtrk_pitch);
  //stitched tracks info ends here
  fTree->Branch("pdg_code", &pdg_code);
  fTree->Branch("n_beamparticle", &n_beamparticle);
  fTree->Branch("n_daughter", &n_daughter);
  fTree->Branch("isPrimary", &isPrimary);
  fTree->Branch("pfp_self", &pfp_self);
  //fTree->Branch("pfp_parent", &pfp_parent);
  fTree->Branch("pfp_daughter", &pfp_daughter);
  fTree->Branch("primtrk_trktag", &primtrk_trktag);
  //Michel tagging
  fTree->Branch("Mendhitssecondary",&Mendhitssecondary);
  fTree->Branch("Msecondarystartx",&Msecondarystartx);
  fTree->Branch("Msecondarystarty",&Msecondarystarty);
  fTree->Branch("Msecondarystartz",&Msecondarystartz);
  fTree->Branch("Msecondaryendx",&Msecondaryendx);
  fTree->Branch("Msecondaryendy",&Msecondaryendy);
  fTree->Branch("Msecondaryendz",&Msecondaryendz);
  fTree->Branch("MdQmichel",&MdQmichel);
  fTree->Branch("Mprimsectheta",&Mprimsectheta);
  fTree->Branch("MdQtrackbegin",&MdQtrackbegin);
  fTree->Branch("MdQtrackend",&MdQtrackend);
  fTree->Branch("Mtracklengthsecondary",&Mtracklengthsecondary);
  fTree->Branch("MtrackID",&MtrackID);
  fTree->Branch("trackscore",&trackscore);
  fTree->Branch("emscore",&emscore);
  fTree->Branch("michelscore",&michelscore);
  fTree->Branch("nonescore",&nonescore);


}

void protoana::pionanalysis::endJob()
{

}

void protoana::pionanalysis::reset()
{
  for(int k=0; k < 5; k++)
    fNactivefembs[k] = -999;

  tof = -1;
  tofs.clear();
  ch_tofs.clear();

  low_pressure_status = -1;
  high_pressure_status = -1;
  low_pressure = -99;
  high_pressure = -99;

  beamPosx.clear();
  beamPosy.clear();
  beamPosz.clear();
  beamDirx.clear();
  beamDiry.clear();
  beamDirz.clear();
  beamMomentum.clear();

  primtrk_startx.clear();
  primtrk_starty.clear();
  primtrk_startz.clear();  

  primtrk_endx.clear();
  primtrk_endy.clear();
  primtrk_endz.clear();  

  primtrk_Dirx.clear();
  primtrk_Diry.clear();
  primtrk_Dirz.clear();
  primtrklen.clear();
  primtrkID.clear();
  primtrk_trktag.clear();

  //primtrk_dqdx->clear();
  //primtrk_resrange->clear();
  //primtrk_dedx->clear();
  primtrk_dqdx.clear();
  primtrk_resrange.clear();
  primtrk_dedx.clear();
  primtrk_range.clear();
  primtrk_hitx.clear();
  primtrk_hity.clear();
  primtrk_hitz.clear();
  primtrk_pitch.clear();
  //wire info
  wireno_2.clear();
  peakTime_2.clear();
  dq_2.clear();
  trkhitx2.clear();
  trkhity2.clear();
  trkhitz2.clear();
  //****************

  pdg_code.clear();
  n_beamparticle.clear();
  n_daughter.clear();

  isPrimary.clear();
  pfp_self.clear();
  //pfp_parent.clear();
  pfp_daughter.clear();

  cosine_beam_primtrk=-99;
  ////stitched track
  //New variables Michel
  Mendhitssecondary.clear();
  Msecondarystartx.clear();
  Msecondaryendx.clear();
  Msecondarystarty.clear();
  Msecondaryendy.clear();
  Msecondarystartz.clear();
  Msecondaryendz.clear();
  MdQmichel.clear();
  Mprimsectheta.clear();
  MdQtrackbegin.clear();
  MdQtrackend.clear();
  Mtracklengthsecondary.clear();
  MtrackID.clear();
  trackscore.clear();
  emscore.clear();
  michelscore.clear();
  nonescore.clear();
}

DEFINE_ART_MODULE(protoana::pionanalysis)
