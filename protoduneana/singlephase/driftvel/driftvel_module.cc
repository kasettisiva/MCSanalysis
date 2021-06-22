// Class:       driftvel
////////////////////////////////////////////////////////////////////////
// Plugin Type: analyzer (art v2_09_06)
// File:        driftvel_module.cc
//
// Generated at Wed Mar  7 09:17:09 2018 by Ajib Paudel using cetskelgen
// from cetlib version v3_01_03.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "lardataobj/RawData/RDTimeStamp.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"

#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>


const int kMaxTracks  = 300;
using namespace std;

namespace protoana{
  class driftvel : public art::EDAnalyzer {
  public:

    explicit driftvel(fhicl::ParameterSet const& pset);
    virtual ~driftvel();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void reset();
    
  private:
    ProtoDUNEDataUtils fDataUtils;
    TTree* fEventTree;
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;
    Double_t evttime; 
    Int_t    all_trks;
    int fNactivefembs[6];
    Float_t  trackthetaxz[kMaxTracks];
    Float_t  trackthetayz[kMaxTracks];
    Float_t trkstartx[kMaxTracks];
    Float_t trkstarty[kMaxTracks];
    Float_t trkstartz[kMaxTracks];
    Float_t trkendx[kMaxTracks];
    Float_t trkendy[kMaxTracks];
    Float_t trkendz[kMaxTracks];
    Float_t trklen[kMaxTracks];
    Int_t TrkID[kMaxTracks]; 
    Float_t  xprojectedlen[kMaxTracks];
    Int_t ntrkhits[kMaxTracks][3];
    Float_t hit_peakT[kMaxTracks][3][3000];
    Int_t hit_tpc[kMaxTracks][3][3000];
    Int_t hit_wire[kMaxTracks][3][3000];
    Int_t hit_channel[kMaxTracks][3][3000];
    Float_t trkhitx[kMaxTracks][3][3000];
    Float_t trkhity[kMaxTracks][3][3000];
    Float_t trkhitz[kMaxTracks][3][3000];
    Float_t trkdq_amp[kMaxTracks][3][3000];
    Float_t trkdq_int[kMaxTracks][3][3000];
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    bool  fSaveTrackInfo;
    // bool fSaveHitInfo;
    bool  fSaveCaloInfo;
  };

  //========================================================================
  driftvel::driftvel(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset),
    fDataUtils                  (pset.get<fhicl::ParameterSet>("DataUtils")),
    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","")         ), 
    fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
    fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
    fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false)),
    fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false))
   
    //  fSaveHitInfo              (pset.get< bool>("SaveHitInfo",false))

  {
    if (fSaveTrackInfo == false) fSaveCaloInfo = false;
  }
 
  //========================================================================
  driftvel::~driftvel(){
  }
  //========================================================================

  //========================================================================
  void driftvel::beginJob(){
    std::cout<<"job begin..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
    fEventTree->Branch("event", &event,"event/I");
    fEventTree->Branch("evttime",&evttime,"evttime/D");
    fEventTree->Branch("run", &run,"run/I");
    fEventTree->Branch("Nactivefembs",&fNactivefembs,"Nactivefembs[5]/I");
    fEventTree->Branch("subrun", &subrun,"surbrun/I");
    fEventTree->Branch("all_trks",&all_trks,"all_trks/I");
    fEventTree->Branch("ntrkhits",ntrkhits,"ntrkhits[all_trks][3]/I");
    fEventTree->Branch("xprojectedlen",xprojectedlen,"xprojectedlen[all_trks]/F");
    fEventTree->Branch("trackthetaxz",trackthetaxz,"trackthetaxz[all_trks]/F");
    fEventTree->Branch("trackthetayz",trackthetayz,"trackthetayz[all_trks]/F");
    fEventTree->Branch("trkstartx",trkstartx,"trkstartx[all_trks]/F");
    fEventTree->Branch("trkstarty",trkstarty,"trkstarty[all_trks]/F");
    fEventTree->Branch("trkstartz",trkstartz,"trkstartz[all_trks]/F");
    fEventTree->Branch("trkendx",trkendx,"trkendx[all_trks]/F");
    fEventTree->Branch("trkendy",trkendy,"trkendy[all_trks]/F");
    fEventTree->Branch("trkendz",trkendz,"trkendz[all_trks]/F");
    fEventTree->Branch("trklen",trklen,"trklen[all_trks]/F");
    fEventTree->Branch("TrkID",TrkID,"TrkID[all_trks]/I");
    fEventTree->Branch("hit_peakT",hit_peakT,"hit_peakT[all_trks][3][3000]/F");
    fEventTree->Branch("hit_tpc",hit_tpc,"hit_tpc[all_trks][3][3000]/I");
    fEventTree->Branch("hit_wire",hit_wire,"hit_wire[all_trks][3][3000]/I");
    fEventTree->Branch("hit_channel",hit_channel,"hit_channel[all_trks][3][3000]/I");
    fEventTree->Branch("trkhitx",trkhitx,"trkhitx[all_trks][3][3000]/F");
    fEventTree->Branch("trkhity",trkhity,"trkhity[all_trks][3][3000]/F");
    fEventTree->Branch("trkhitz",trkhitz,"trkhitz[all_trks][3][3000]/F");
    fEventTree->Branch("trkdq_int",trkdq_int,"trkdq_int[all_trks][3][3000]/F");
    fEventTree->Branch("trkdq_amp",trkdq_amp,"trkdq_amp[all_trks][3][3000]/F");
  }

  //========================================================================
  void driftvel::endJob(){     

  }

  //========================================================================
  void driftvel::beginRun(const art::Run&){
    mf::LogInfo("driftvel")<<"begin run..."<<std::endl;
  }
  //========================================================================

  //========================================================================

  //========================================================================

  void driftvel::analyze( const art::Event& evt){//analyze
    reset();  
   
    art::Handle< std::vector<recob::Track> > trackListHandle;
    art::Handle< std::vector<recob::PFParticle> > PFPListHandle; 
   
    std::vector<art::Ptr<recob::Track> > tracklist;
    std::vector<art::Ptr<recob::PFParticle> > pfplist;


    if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
    if(evt.getByLabel("pandora",PFPListHandle)) art::fill_ptr_vector(pfplist, PFPListHandle);
  
   
    art::Handle< std::vector<recob::Hit> > hitListHandle; // to get information about the hits
    std::vector<art::Ptr<recob::Hit>> hitlist;
    if(evt.getByLabel(fHitsModuleLabel, hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
       
    // art::FindManyP<recob::Hit> fmtht(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits

    art::FindManyP<anab::T0> trk_t0_assn_v(PFPListHandle, evt ,"pandora");
    art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle,evt,"pandoraTrack");
    art::FindManyP<anab::T0> fmT0(trackListHandle, evt ,"pmtrack");


    run = evt.run();
    subrun = evt.subRun();
    event = evt.id().event();
    art::Timestamp ts = evt.time();
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime=tts.AsDouble();



    // Get number of active fembs
    if(!evt.isRealData()){
      for(int k=0; k < 6; k++)
	fNactivefembs[k] = 20;
    }
    else{
      for(int k=0; k < 6; k++)
	fNactivefembs[k] = fDataUtils.GetNActiveFembsForAPA(evt, k);
    }


     
    all_trks=0;
    size_t NTracks = tracklist.size();
    for(size_t i=0; i<NTracks;++i){
      art::Ptr<recob::Track> ptrack(trackListHandle, i);
      /* if(fTrackModuleLabel=="pandoraTrack"){
	 std::vector<art::Ptr<recob::PFParticle>> pfps=pfp_trk_assn.at(i);
	 if(!pfps.size()) continue;
	 std::vector<art::Ptr<anab::T0>> t0s=trk_t0_assn_v.at(pfps[0].key());
	 if(!t0s.size()) continue;
	 //auto t0 = t0s.at(0);
	 // double t_zero=t0->Time();
	 }
      
	 if(fTrackModuleLabel=="pmtrack"){
	 std::vector<art::Ptr<anab::T0>> T0s=fmT0.at(i);
	 if(T0s.size()==0)
	 continue;
	 }
      */
      all_trks++;
      const recob::Track& track = *ptrack;
      if(track.Length()<300) continue;
      // TVector3 pos, dir_start, dir_end, end;
      auto pos = track.Vertex();
      auto dir_start = track.VertexDirection();
      //auto dir_end   = track.EndDirection();
      auto end = track.End();
      if(TMath::Abs(end.X()-pos.X())<330 ||TMath::Abs(end.X()-pos.X())>1500 ) continue;
      double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
      double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
      xprojectedlen[all_trks-1]=TMath::Abs(end.X()-pos.X());
      trackthetaxz[all_trks-1]=theta_xz;
      trackthetayz[all_trks-1]=theta_yz;
      trkstartx[all_trks-1]=pos.X();
      trkstarty[all_trks-1]=pos.Y();
      trkstartz[all_trks-1]=pos.Z();
      trkendx[all_trks-1]=end.X();
      trkendy[all_trks-1]=end.Y();
      trkendz[all_trks-1]=end.Z();
      trklen[all_trks-1]=track.Length();
      TrkID[all_trks-1]=track.ID();
      int planenum=999;
      int nhits_0=0;
      int nhits_1=0;
      int nhits_2=0;
      float xpos=-9999;
      float ypos=-9999;
      float zpos=-9999;
      // auto vhit=fmtht.at(i);
      cout<<"fmthm is valid "<<fmthm.isValid()<<endl;
      if(fmthm.isValid()){
	auto vhit=fmthm.at(i);
	auto vmeta=fmthm.data(i);
	cout<<"vhit vmetadata"<<endl;
	for (size_t ii = 0; ii<vhit.size(); ++ii){ //loop over all meta data hit
	  bool fBadhit = false;
	  if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max())){
	    fBadhit = true;
	    cout<<"fBadHit"<<fBadhit<<endl;
	    continue;
	  }
	  cout<<"first ii"<<ii<<endl;
	  if (vmeta[ii]->Index()>=tracklist[i]->NumberTrajectoryPoints()){
	    throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<tracklist[i]->NumberTrajectoryPoints()<<" for track index "<<i<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
	  }
	  cout<<"second ii"<<ii<<endl;
	  if (!tracklist[i]->HasValidPoint(vmeta[ii]->Index())){
	    fBadhit = true;
	    cout<<"had valid point "<<fBadhit<<endl;
	    continue;
	  }
	  cout<<"3rd ii"<<ii<<endl;
	  // TVector3 loc = tracklist[i]->LocationAtPoint(vmeta[ii]->Index());
	  auto loc = tracklist[i]->LocationAtPoint(vmeta[ii]->Index());
	  xpos=loc.X();

	  ypos=loc.Y();
	  zpos=loc.Z();
	  cout<<"x, y, z "<<xpos<<"  "<<ypos<<"  "<<zpos<<endl;
	  cout<<"BadHit"<<fBadhit<<endl;
	  if (fBadhit) continue; //HY::If BAD hit, skip this hit and go next
	  if (zpos<-100) continue; //hit not on track
	  planenum=vhit[ii]->WireID().Plane;
	  cout<<"plane num "<<planenum;
	  if(planenum==0){	 
	    nhits_0++;
	    cout<<"pl 0"<<endl;
	    hit_peakT[all_trks-1][planenum][nhits_0-1]=vhit[ii]->PeakTime();
	    int k=vhit[ii]->PeakTime();
	    cout<<"int k"<<k<<endl; 
	    cout<<"peakT"<<vhit[ii]->PeakTime();
	    hit_peakT[all_trks-1][planenum][nhits_0-1]=vhit[ii]->PeakTime();
	    hit_tpc[all_trks-1][planenum][nhits_0-1]=vhit[ii] ->WireID().TPC;
	    hit_wire[all_trks-1][planenum][nhits_0-1]=vhit[ii] ->WireID().Wire;
	    trkdq_int[all_trks-1][planenum][nhits_0-1]=vhit[ii] ->Integral();
	    trkdq_amp[all_trks-1][planenum][nhits_0-1]=vhit[ii] ->PeakAmplitude();
	    hit_channel[all_trks-1][planenum][nhits_0-1]=vhit[ii]->Channel();
	    cout<<vhit[ii]->WireID().Wire;
	    cout<<"peakT"<<vhit[ii]->PeakTime();
	    trkhitx[all_trks-1][planenum][nhits_0-1]=xpos;
	    trkhity[all_trks-1][planenum][nhits_0-1]=ypos;
	    trkhitz[all_trks-1][planenum][nhits_0-1]=zpos;
	  }//planenum 0
	  if(planenum==1){	 
	    nhits_1++;
	    cout<<"pl 1"<<endl;
	    hit_peakT[all_trks-1][planenum][nhits_1-1]=vhit[ii]->PeakTime();
	    hit_tpc[all_trks-1][planenum][nhits_1-1]=vhit[ii] ->WireID().TPC;
	    hit_wire[all_trks-1][planenum][nhits_1-1]=vhit[ii] ->WireID().Wire;
	    trkdq_int[all_trks-1][planenum][nhits_1-1]=vhit[ii] ->Integral();
	    trkdq_amp[all_trks-1][planenum][nhits_1-1]=vhit[ii] ->PeakAmplitude();
	    hit_channel[all_trks-1][planenum][nhits_1-1]=vhit[ii]->Channel();
	    trkhitx[all_trks-1][planenum][nhits_1-1]=xpos;
	    trkhity[all_trks-1][planenum][nhits_1-1]=ypos;
	    trkhitz[all_trks-1][planenum][nhits_1-1]=zpos;
	  }//planenum 1
	  if(planenum==2){	 
	    nhits_2++;
	    cout<<"pl 2"<<endl;
	    hit_peakT[all_trks-1][planenum][nhits_2-1]=vhit[ii]->PeakTime();
	    hit_tpc[all_trks-1][planenum][nhits_2-1]=vhit[ii] ->WireID().TPC;
	    hit_wire[all_trks-1][planenum][nhits_2-1]=vhit[ii] ->WireID().Wire;
	    trkdq_int[all_trks-1][planenum][nhits_2-1]=vhit[ii] ->Integral();
	    trkdq_amp[all_trks-1][planenum][nhits_2-1]=vhit[ii] ->PeakAmplitude();
	    hit_channel[all_trks-1][planenum][nhits_2-1]=vhit[ii]->Channel();
	    trkhitx[all_trks-1][planenum][nhits_2-1]=xpos;
	    trkhity[all_trks-1][planenum][nhits_2-1]=ypos;
	  
	    trkhitz[all_trks-1][planenum][nhits_2-1]=zpos;
	  }//planenum 2
	}//loop over vhit
      }//fmthm valid
      ntrkhits[all_trks-1][0]=nhits_0;
      ntrkhits[all_trks-1][1]=nhits_1;
      ntrkhits[all_trks-1][2]=nhits_2;
    } // loop over trks...
    
    fEventTree->Fill();
  } // end of analyze function
	   
  /////////////////// Defintion of reset function ///////////
  void driftvel::reset(){
    run = -9999;
    subrun = -9999;
    event = -9999;
    evttime = -9999;
    all_trks = -9999;
    for(int k=0; k < 6; k++)
      fNactivefembs[k] = -999;
    for(int i=0; i<kMaxTracks; i++){
      trackthetaxz[i]=-9999;
      trackthetayz[i]=-9999;
      trkstartx[i]=-9999;
      trkstarty[i]=-9999;
      trkstartz[i]=-9999;
      trkendx[i]=-9999;
      trkendy[i]=-9999;
      trkendz[i]=-9999;
      trklen[i]=-9999;
      TrkID[i]=-9999;
      xprojectedlen[i]=-9999;
      ntrkhits[i][0]=-9999;
      ntrkhits[i][1]=-9999;
      ntrkhits[i][2]=-9999;
      for(int j=0; j<3000; j++){
	hit_peakT[i][0][j]=-9999;
	hit_peakT[i][1][j]=-9999;
	hit_peakT[i][2][j]=-9999;
	hit_tpc[i][0][j]=-9999;
	hit_tpc[i][1][j]=-9999;
	hit_tpc[i][2][j]=-9999;
	hit_wire[i][0][j]=-9999;
	hit_wire[i][1][j]=-9999;
	hit_wire[i][2][j]=-9999;
	trkhitx[i][0][j]=-9999;
	trkhitx[i][1][j]=-9999;
	trkhitx[i][2][j]=-9999;
	trkhity[i][0][j]=-9999;
	trkhity[i][1][j]=-9999;
	trkhity[i][2][j]=-9999;
	trkhitz[i][0][j]=-9999;
	trkhitz[i][1][j]=-9999;
	trkhitz[i][2][j]=-9999;
	hit_channel[i][0][j]=-9999;
	hit_channel[i][1][j]=-9999;
	hit_channel[i][2][j]=-9999;
	trkdq_int[i][0][j]=-9999;
	trkdq_int[i][1][j]=-9999;
	trkdq_int[i][2][j]=-9999;
	trkdq_amp[i][0][j]=-9999;
	trkdq_amp[i][1][j]=-9999;
	trkdq_amp[i][2][j]=-9999;
      }
    }
  }
  //////////////////////// End of definition ///////////////	
	  
  DEFINE_ART_MODULE(driftvel)
}


