////////////////////////////////////////////////////////////////////////
// Class:       pionanalysis
// Plugin Type: analyzer (art v2_11_02)
// File:        pionanalysis_module.cc
// Pionanalysis module: Contact Heng-Ye Liao (liao@phys.ksu.edu) 
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

#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
#include "protoduneana/Analysis/ProtoDUNETrackUtils.h"
#include "protoduneana/Analysis/ProtoDUNEShowerUtils.h"
#include "protoduneana/Analysis/ProtoDUNETruthUtils.h"
#include "protoduneana/Analysis/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Analysis/ProtoDUNEDataUtils.h"
#include "protoduneana/Analysis/ProtoDUNEBeamlineUtils.h"

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

//const int kMaxTracks  = 1000;
//const int kMaxHits = 10000;

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
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
  protoana::ProtoDUNEDataUtils dataUtil;

  // fcl parameters for PFP particles
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
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
  
  //carlo info
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
};


protoana::pionanalysis::pionanalysis(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  //fSpacePointModuleLabel(p.get< art::InputTag >("SpacePointModuleLabel")),
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fBeamlineUtils(p.get<fhicl::ParameterSet>("BeamlineUtils")),
  //fBeamModuleLabel(p.get<std::string>("BeamModuleLabel")),
  
  //fTimeDecoderModuleLabel(p.get< art::InputTag >("TimeDecoderModuleLabel")),

  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
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

  //HY::For Carlo info
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryTag);
  art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle, evt, "pandoraTrack");

  // Implementation of required member function here.
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
    		const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
    		if(thisTrack != 0x0) {
					std::cout << "Beam particle is track-like" << std::endl;
					nTrack++;
					primtrk_trktag.push_back(1);

                			//HY::Get the Calorimetry(s) from thisTrack
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

  } //if beam proton

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

  fTree->Branch("pdg_code", &pdg_code);
  fTree->Branch("n_beamparticle", &n_beamparticle);
  fTree->Branch("n_daughter", &n_daughter);
  fTree->Branch("isPrimary", &isPrimary);
  fTree->Branch("pfp_self", &pfp_self);
  //fTree->Branch("pfp_parent", &pfp_parent);
  fTree->Branch("pfp_daughter", &pfp_daughter);
  fTree->Branch("primtrk_trktag", &primtrk_trktag);

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

  pdg_code.clear();
  n_beamparticle.clear();
  n_daughter.clear();

  isPrimary.clear();
  pfp_self.clear();
  //pfp_parent.clear();
  pfp_daughter.clear();

  cosine_beam_primtrk=-99;

}

DEFINE_ART_MODULE(protoana::pionanalysis)
