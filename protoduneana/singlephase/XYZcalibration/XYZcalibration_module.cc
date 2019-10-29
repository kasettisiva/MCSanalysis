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
#include "lardataobj/RecoBase/PFParticle.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"

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

const int kMaxTracks  = 200;

using namespace std;

namespace dune{

  class XYZcalibration : public art::EDAnalyzer {
  public:

    explicit XYZcalibration(fhicl::ParameterSet const& pset);
    virtual ~XYZcalibration();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void reset();
    
  private:
    TTree* fEventTree;
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;
    Double_t evttime; 
    Int_t    year_month_date;
    Int_t    hour_min_sec;
    Int_t    cross_trks;
    Int_t stopping_trks;
    Int_t    all_trks;
    Float_t  xprojectedlen[kMaxTracks];
    Float_t  trackthetaxz[kMaxTracks];
    Float_t  trackthetayz[kMaxTracks];
    Float_t trkstartx[kMaxTracks];
    Float_t trkstarty[kMaxTracks];
    Float_t trkstartz[kMaxTracks];
    Float_t trkendx[kMaxTracks];
    Float_t trkendy[kMaxTracks];
    Float_t trkendz[kMaxTracks];
    Float_t trklen[kMaxTracks];
    Int_t    TrkID[kMaxTracks]; 
    Float_t  trkstartcosxyz[kMaxTracks][3];
    Float_t  trkendcosxyz[kMaxTracks][3];
    Int_t    ntrkhits[kMaxTracks][3];
    Float_t  trkdqdx[kMaxTracks][3][3000];
    Float_t  trkdedx[kMaxTracks][3][3000];
    Float_t  trkresrange[kMaxTracks][3][3000];
    Float_t  trkhitx[kMaxTracks][3][3000];
    Float_t  trkhity[kMaxTracks][3][3000];
    Float_t  trkhitz[kMaxTracks][3][3000];
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    bool  fSaveCaloInfo;
    bool  fSaveTrackInfo;
  }; 

  //========================================================================
  XYZcalibration::XYZcalibration(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset),
    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","")         ),
    fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
    fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
    fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false)),
    fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false))
  {
    if (fSaveTrackInfo == false) fSaveCaloInfo = false;
  }
 
  //========================================================================
  XYZcalibration::~XYZcalibration(){
  }
  //========================================================================

  //========================================================================
  void XYZcalibration::beginJob(){
    std::cout<<"job begin..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;
    // deltaX = tfs->make<TH1D>("deltaX","Plot of deltaX",400,0,400);
    // deltaX->Sumw2();
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
    fEventTree->Branch("event", &event,"event/I");
    fEventTree->Branch("evttime",&evttime,"evttime/D");
    fEventTree->Branch("run", &run,"run/I");
    fEventTree->Branch("subrun", &subrun,"surbrun/I");
    fEventTree->Branch("year_month_date", &year_month_date,"year_month_date/I");
    fEventTree->Branch("hour_min_sec", &hour_min_sec,"hour_min_sec/I");
    fEventTree->Branch("cross_trks",&cross_trks,"cross_trks/I");
    fEventTree->Branch("stopping_trks",&stopping_trks,"stopping_trks/I");
    fEventTree->Branch("all_trks",&all_trks,"all_trks/I");
    fEventTree->Branch("xprojectedlen",xprojectedlen,"xprojectedlen[all_trks]/F");
    fEventTree->Branch("trackthetaxz",trackthetaxz,"trackthetaxz[cross_trks]/F");
    fEventTree->Branch("trackthetayz",trackthetayz,"trackthetayz[cross_trks]/F");
    fEventTree->Branch("trkstartx",trkstartx,"trkstartx[cross_trks]/F");
    fEventTree->Branch("trkstarty",trkstarty,"trkstarty[cross_trks]/F");
    fEventTree->Branch("trkstartz",trkstartz,"trkstartz[cross_trks]/F");
    fEventTree->Branch("trkendx",trkendx,"trkendx[cross_trks]/F");
    fEventTree->Branch("trkendy",trkendy,"trkendy[cross_trks]/F");
    fEventTree->Branch("trkendz",trkendz,"trkendz[cross_trks]/F");
    fEventTree->Branch("trklen",trklen,"trklen[cross_trks]/F");
    fEventTree->Branch("TrkID",TrkID,"TrkID[cross_trks]/I");
    fEventTree->Branch("trkstartcosxyz",trkstartcosxyz,"trkstartcosxyz[cross_trks][3]/F");
    fEventTree->Branch("trkendcosxyz",trkendcosxyz,"trkendcosxyz[cross_trks][3]/F");
    fEventTree->Branch("ntrkhits",ntrkhits,"ntrkhits[cross_trks][3]/I");
    fEventTree->Branch("trkdqdx",trkdqdx,"trkdqdx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkdedx",trkdedx,"trkdedx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkresrange",trkresrange,"trkresrange[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitx",trkhitx,"trkhitx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhity",trkhity,"trkhity[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitz",trkhitz,"trkhitz[cross_trks][3][3000]/F");
  }

  //========================================================================
  void XYZcalibration::endJob(){     

  }

  //========================================================================
  void XYZcalibration::beginRun(const art::Run&){
    mf::LogInfo("XYZcalibration")<<"begin run..."<<std::endl;
  }
  //========================================================================

  //========================================================================

  //========================================================================

  void XYZcalibration::analyze( const art::Event& evt){
    reset();  
     
    art::Handle< std::vector<recob::Track> > trackListHandle;
    art::Handle< std::vector<recob::PFParticle> > PFPListHandle; 
   
    std::vector<art::Ptr<recob::Track> > tracklist;
    std::vector<art::Ptr<recob::PFParticle> > pfplist;
  
    if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
    if(evt.getByLabel("pandora",PFPListHandle)) art::fill_ptr_vector(pfplist, PFPListHandle);
  
    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
    art::FindManyP<anab::T0> trk_t0_assn_v(PFPListHandle, evt ,"pandora");
    art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle,evt,"pandoraTrack");
    art::FindManyP<anab::T0> fmT0(trackListHandle, evt ,"pmtrack");

    run = evt.run();
    subrun = evt.subRun();
    event = evt.id().event();
    art::Timestamp ts = evt.time();
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime=tts.AsDouble();
     
    UInt_t year=0;
    UInt_t month=0;
    UInt_t day=0;
     
    year_month_date=tts.GetDate(kTRUE,0,&year,&month,&day);
     
    UInt_t hour=0;
    UInt_t min=0;
    UInt_t sec=0;
     
    hour_min_sec=tts.GetTime(kTRUE,0,&hour,&min,&sec);
  
    cross_trks=0;
    stopping_trks=0;
    all_trks=0;
  
    size_t NTracks = tracklist.size();
    for(size_t i=0; i<NTracks;++i){
      art::Ptr<recob::Track> ptrack(trackListHandle, i);
      if(fTrackModuleLabel=="pandoraTrack"){
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
      all_trks++;
      std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(i);     

      const recob::Track& track = *ptrack;
      auto pos = track.Vertex();
      auto dir_start = track.VertexDirection();
      auto dir_end   = track.EndDirection();
      auto end = track.End();
      double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
      double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
      //deltaX->Fill(X);
      cross_trks++;
      if(std::abs(pos.X())<340 && std::abs(end.X())<340 &&  pos.Y()>20 && pos.Y()<580 && end.Y()>20 && end.Y()<580 &&  pos.Z()>20 && pos.Z()<675 && end.Z()>20 && end.Z()<675) stopping_trks++;//this infact gives both ends confined, we have not counted the stopping muon tracks in this module
      trackthetaxz[cross_trks-1]=theta_xz;
      trackthetayz[cross_trks-1]=theta_yz;
      trkstartx[cross_trks-1]=pos.X();
      trkstarty[cross_trks-1]=pos.Y();
      trkstartz[cross_trks-1]=pos.Z();
      trkendx[cross_trks-1]=end.X();
      trkendy[cross_trks-1]=end.Y();
      trkendz[cross_trks-1]=end.Z();
      trklen[cross_trks-1]=track.Length();
      TrkID[cross_trks-1]=track.ID();
      trkstartcosxyz[cross_trks-1][0]=dir_start.X();
      trkstartcosxyz[cross_trks-1][1]=dir_start.Y();
      trkstartcosxyz[cross_trks-1][2]=dir_start.Z();
      trkendcosxyz[cross_trks-1][0]=dir_end.X();
      trkendcosxyz[cross_trks-1][1]=dir_end.Y();
      trkendcosxyz[cross_trks-1][2]=dir_end.Z();
				   
      for(size_t ical = 0; ical<calos.size(); ++ical){
	if(!calos[ical]) continue;
	if(!calos[ical]->PlaneID().isValid) continue;
	int planenum = calos[ical]->PlaneID().Plane;
	if(planenum<0||planenum>2) continue;
	const size_t NHits = calos[ical] -> dEdx().size();
	ntrkhits[cross_trks-1][planenum]=int(NHits);
	for(size_t iHit = 0; iHit < NHits; ++iHit){
	  const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
	  trkdqdx[cross_trks-1][planenum][iHit]=(calos[ical] -> dQdx())[iHit];
	  trkdedx[cross_trks-1][planenum][iHit]=(calos[ical] -> dEdx())[iHit];
	  trkresrange[cross_trks-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
	  trkhitx[cross_trks-1][planenum][iHit]=TrkPos.X();
	  trkhity[cross_trks-1][planenum][iHit]=TrkPos.Y();
	  trkhitz[cross_trks-1][planenum][iHit]=TrkPos.Z();
	} // loop over iHit..
      } // loop over ical 2 nd time...
    } // loop over trks...
    fEventTree->Fill();
  } // end of analyze function
	   
  /////////////////// Defintion of reset function ///////////
  void XYZcalibration::reset(){
    run = -9999;
    subrun = -9999;
    event = -9999;
    evttime = -9999;
    cross_trks = -9999;
    all_trks = -9999;
    year_month_date=-9999;
    hour_min_sec=-9999;
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
      trkstartcosxyz[i][0]=-9999;
      trkstartcosxyz[i][1]=-9999;
      trkstartcosxyz[i][2]=-9999; 
      trkendcosxyz[i][0]=-9999;
      trkendcosxyz[i][1]=-9999;
      trkendcosxyz[i][2]=-9999;
      ntrkhits[i][0] = -9999;
      ntrkhits[i][1] = -9999;
      ntrkhits[i][2] = -9999;
      for(int j=0; j<3; j++){
	for(int k=0; k<3000; k++){
	  trkdqdx[i][j][k]=-9999;
	  trkdedx[i][j][k]=-9999;
	  trkresrange[i][j][k]=-9999;
	  trkhitx[i][j][k]=-9999;
	  trkhity[i][j][k]=-9999;
	  trkhitz[i][j][k]=-9999;
	}
      }
    }
  }
  //////////////////////// End of definition ///////////////	
	  
  DEFINE_ART_MODULE(XYZcalibration)
}
