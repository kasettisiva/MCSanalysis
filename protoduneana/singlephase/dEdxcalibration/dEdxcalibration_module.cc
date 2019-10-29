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

const int kMaxTracks  = 200;

using namespace std;

namespace dune{

  class dEdxcalibration : public art::EDAnalyzer {
  public:

    explicit dEdxcalibration(fhicl::ParameterSet const& pset);
    virtual ~dEdxcalibration();

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
    Double_t Total_energy[kMaxTracks];
    Double_t EndE[kMaxTracks];
    Int_t pdg[kMaxTracks];
    Double_t evttime; 
    Int_t    year_month_date;
    Int_t    hour_min_sec;
    Int_t    cross_trks;
    Int_t true_stopping_muons;
    Int_t mc_stopping50;
    Int_t mc_stopping20;
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
    Int_t    TrackId_geant[kMaxTracks]; 
    Int_t origin[kMaxTracks];
    Float_t EndPointx[kMaxTracks];
    Float_t EndPointy[kMaxTracks];
    Float_t EndPointz[kMaxTracks];
    Float_t EndPointx_corrected[kMaxTracks];
    Float_t EndPointy_corrected[kMaxTracks];
    Float_t EndPointz_corrected[kMaxTracks];
    Float_t StartPointx[kMaxTracks];
    Float_t StartPointy[kMaxTracks];
    Float_t StartPointz[kMaxTracks];
    Float_t  trkstartcosxyz[kMaxTracks][3];
    Float_t  trkendcosxyz[kMaxTracks][3];
    Int_t    ntrkhits[kMaxTracks][3];
    Float_t  trkdqdx[kMaxTracks][3][3000];
    Float_t  trkdedx[kMaxTracks][3][3000];
    Float_t  trkresrange[kMaxTracks][3][3000];
    Float_t  trkhitx[kMaxTracks][3][3000];
    Float_t  trkhity[kMaxTracks][3][3000];
    Float_t  trkhitz[kMaxTracks][3][3000];
    Float_t  trkpitch[kMaxTracks][3][3000];
    Float_t peakT_max[kMaxTracks];
    Float_t peakT_min[kMaxTracks];
    Int_t daughter_electrons[kMaxTracks];
    Int_t daughter_positrons[kMaxTracks];
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    bool  fSaveCaloInfo;
    bool  fSaveTrackInfo;
  }; 

  //========================================================================
  dEdxcalibration::dEdxcalibration(fhicl::ParameterSet const& pset) :
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
  dEdxcalibration::~dEdxcalibration(){
  }
  //========================================================================

  //========================================================================
  void dEdxcalibration::beginJob(){
    std::cout<<"job begin..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
    fEventTree->Branch("event", &event,"event/I");
    fEventTree->Branch("evttime",&evttime,"evttime/D");
    fEventTree->Branch("run", &run,"run/I");
    fEventTree->Branch("subrun", &subrun,"surbrun/I");
    fEventTree->Branch("year_month_date", &year_month_date,"year_month_date/I");
    fEventTree->Branch("hour_min_sec", &hour_min_sec,"hour_min_sec/I");
    fEventTree->Branch("cross_trks",&cross_trks,"cross_trks/I");
    fEventTree->Branch("mc_stopping20",&mc_stopping20,"mc_stopping20/I");
    fEventTree->Branch("mc_stopping50",&mc_stopping50,"mc_stopping50/I");
    fEventTree->Branch("true_stopping_muons",&true_stopping_muons,"true_stopping_muons/I");
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
    fEventTree->Branch("Total_energy",Total_energy,"Total_energy[cross_trks]/D");
    fEventTree->Branch("EndE",EndE,"EndE[cross_trks]/D");
    fEventTree->Branch("peakT_max",peakT_max,"peakT_max[cross_trks]/F");
    fEventTree->Branch("peakT_min",peakT_min,"peakT_min[cross_trks]/F");
    fEventTree->Branch("pdg",pdg,"pdg[cross_trks]/I");
    fEventTree->Branch("origin",origin,"origin[cross_trks]/I");
    fEventTree->Branch("TrkID",TrkID,"TrkID[cross_trks]/I");
    fEventTree->Branch("TrackId_geant",TrackId_geant,"TrackId_geant[cross_trks]/I");
    fEventTree->Branch("daughter_electrons",daughter_electrons,"daughter_electrons[cross_trks]/I");
    fEventTree->Branch("daughter_positrons",daughter_positrons,"daughter_positrons[cross_trks]/I");
    fEventTree->Branch("EndPointx",EndPointx,"EndPointx[cross_trks]/F");
    fEventTree->Branch("EndPointy",EndPointy,"EndPointy[cross_trks]/F");
    fEventTree->Branch("EndPointz",EndPointz,"EndPointz[cross_trks]/F");
    fEventTree->Branch("EndPointx_corrected",EndPointx_corrected,"EndPointx_corrected[cross_trks]/F");
    fEventTree->Branch("EndPointy_corrected",EndPointy_corrected,"EndPointy_corrected[cross_trks]/F");
    fEventTree->Branch("EndPointz_corrected",EndPointz_corrected,"EndPointz_corrected[cross_trks]/F");
    fEventTree->Branch("StartPointx",StartPointx,"StartPointx[cross_trks]/F");
    fEventTree->Branch("StartPointy",StartPointy,"StartPointy[cross_trks]/F");
    fEventTree->Branch("StartPointz",StartPointz,"StartPointz[cross_trks]/F");
    fEventTree->Branch("trkstartcosxyz",trkstartcosxyz,"trkstartcosxyz[cross_trks][3]/F");
    fEventTree->Branch("trkendcosxyz",trkendcosxyz,"trkendcosxyz[cross_trks][3]/F");
    fEventTree->Branch("ntrkhits",ntrkhits,"ntrkhits[cross_trks][3]/I");
    fEventTree->Branch("trkdqdx",trkdqdx,"trkdqdx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkdedx",trkdedx,"trkdedx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkresrange",trkresrange,"trkresrange[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitx",trkhitx,"trkhitx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhity",trkhity,"trkhity[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitz",trkhitz,"trkhitz[cross_trks][3][3000]/F");
    fEventTree->Branch("trkpitch",trkpitch,"trkpitch[cross_trks][3][3000]/F");
  }

  //========================================================================
  void dEdxcalibration::endJob(){     

  }

  //========================================================================
  void dEdxcalibration::beginRun(const art::Run&){
    mf::LogInfo("dEdxcalibration")<<"begin run..."<<std::endl;
  }
  //========================================================================

  //========================================================================

  //========================================================================

  void dEdxcalibration::analyze( const art::Event& evt){
    reset();  
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;     
    

    art::Handle< std::vector<recob::Track> > trackListHandle;
    art::Handle< std::vector<recob::PFParticle> > PFPListHandle; 
   
    std::vector<art::Ptr<recob::Track> > tracklist;
    std::vector<art::Ptr<recob::PFParticle> > pfplist;


    if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
    if(evt.getByLabel("pandora",PFPListHandle)) art::fill_ptr_vector(pfplist, PFPListHandle);
  
    /******new lines*************************/
    art::Handle< std::vector<recob::Hit> > hitListHandle; // to get information about the hits
    std::vector<art::Ptr<recob::Hit>> hitlist;
    if(evt.getByLabel(fHitsModuleLabel, hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
       
    art::FindManyP<recob::Hit> fmtht(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits


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
    true_stopping_muons=0;
    mc_stopping20=0;
    mc_stopping50=0;
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
      auto pos = track.Vertex<TVector3>();
      auto dir_start = track.VertexDirection<TVector3>();
      auto dir_end   = track.EndDirection<TVector3>();
      auto end = track.End();
      double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
      double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());     
      float starty=pos.Y();
      float startz=pos.Z();
      float endy=end.Y();
      float endz=end.Z();
      float startcosx=dir_start.X();
      float startcosy=dir_start.Y();
      float startcosz=dir_start.Z();
      float endcosx=dir_end.X();
      float endcosy=dir_end.Y();
      float endcosz=dir_end.Z();
      float tracklength=track.Length();
      size_t count = 0;
      /********************************************************************************************************/
      /*********************storing hits for each tracks*******************************/


      std::vector<art::Ptr<recob::Hit>> allHits=fmtht.at(i); //storing hits for ith track
      int trackid=-1;
      art::ServiceHandle<cheat::BackTrackerService> bt_serv;
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      std::map<int,double> trkide;
      std::vector<float> hitpeakT;
      //selection based on hitpeakT
      for(size_t h1=0;h1<allHits.size();h1++){
	hitpeakT.push_back(allHits[h1]->PeakTime());
      }
      float max=-999999;
      float min=-999999;
      min=*min_element(hitpeakT.begin(),hitpeakT.end());
      max=*max_element(hitpeakT.begin(),hitpeakT.end());
      hitpeakT.clear();
      // if(min<250) continue;
      for(size_t h=0; h<allHits.size();h++){
	art::Ptr<recob::Hit> hit=allHits[h];
	std::vector<sim::TrackIDE> eveIDs = bt_serv->HitToTrackIDEs(hit);
	for(size_t e=0;e<eveIDs.size(); ++e){
	  trkide[eveIDs[e].trackID] += eveIDs[e].energy;
	}
      }
      double  maxe = -1;
      double tote = 0;
      for(std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
	tote += ii->second;
	if((ii->second)>maxe){
	  maxe = ii->second;
	  trackid = ii->first;
	 
	}
      }
      float EndE1=-99999;
      int pdg1=-99999;
      int mc_endx=-99999;
      int mc_endy=-99999;
      int mc_endz=-99999;
      const art::Ptr<simb::MCTruth> mc=pi_serv->TrackIdToMCTruth_P(trackid);
      const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(trackid);
      if(particle){
	EndE1 = particle->EndE();
	pdg1=particle->PdgCode();
	mc_endx=particle->EndPosition()[0];
	mc_endy=particle->EndPosition()[1];
	mc_endz=particle->EndPosition()[2];
      };
      if(std::abs(pdg1)==13 && std::abs(mc_endx)<310 && mc_endy>50 && mc_endy<550 && mc_endz>50 && mc_endz<645) mc_stopping50++;
      if(std::abs(pdg1)==13 && std::abs(mc_endx)<340 && mc_endy>20 && mc_endy<580 && mc_endz>20 && mc_endz<675) mc_stopping20++;
      if(std::abs(pdg1)==13 && EndE1<0.106) true_stopping_muons++;
      if(!particle) continue;

     



 
      if(!(((std::abs(pos.X())>340 ||  pos.Y()<20 || pos.Y()>580 || pos.Z()<20 || pos.Z()>675) && !(std::abs(end.X())>310 ||  end.Y()<50 || end.Y()>550 || end.Z()<50 || end.Z()>645))||(!(std::abs(pos.X())>310 ||  pos.Y()<50 || pos.Y()>550 || pos.Z()<50 || pos.Z()>645) && (std::abs(end.X())>340 ||  end.Y()<20 || end.Y()>580 || end.Z()<20 || end.Z()>675)))) continue;
      // if(!(((std::abs(pos.X())>340 ||  pos.Y()<20 || pos.Y()>580 || pos.Z()<20 || pos.Z()>675) && !(std::abs(end.X())>320 ||  end.Y()<40 || end.Y()>560 || end.Z()<40 || end.Z()>655))||(!(std::abs(pos.X())>320 ||  pos.Y()<40 || pos.Y()>560 || pos.Z()<40 || pos.Z()>655) && (std::abs(end.X())>340 ||  end.Y()<20 || end.Y()>580 || end.Z()<20 || end.Z()>675)))) continue;
      stopping_trks++;

      /**********************************broken tracks removal**************************************************************************************/
      for(size_t k=0;k<NTracks;++k){
	art::Ptr<recob::Track> ptrack_k(trackListHandle, k);
	const recob::Track& track_k = *ptrack_k;
	auto pos_k = track_k.Vertex();
	auto dir_pos_k= track_k.VertexDirection();
	auto dir_end_k   = track_k.EndDirection();
	auto end_k = track_k.End();
	if(k==i) continue;
	if((std::abs(((end_k.Y()-pos_k.Y())/(end_k.Z()-pos_k.Z()))*(endz-pos_k.Z())+pos_k.Y()-endy)<30||std::abs(((end_k.Y()-pos_k.Y())/(end_k.Z()-pos_k.Z()))*(startz-pos_k.Z())+pos_k.Y()-starty)<30)&&(std::abs(endcosx*dir_pos_k.X()+endcosy*dir_pos_k.Y()+endcosz*dir_pos_k.Z())>0.97||std::abs(startcosx*dir_pos_k.X()+startcosy*dir_pos_k.Y()+startcosz*dir_pos_k.Z())>0.97||std::abs(endcosx*dir_end_k.X()+endcosy*dir_end_k.Y()+endcosz*dir_end_k.Z())>0.97||std::abs(startcosx*dir_end_k.X()+startcosy*dir_end_k.Y()+startcosz*dir_end_k.Z())>0.97)) break;
	if((std::abs(((end_k.Y()-pos_k.Y())/(end_k.Z()-pos_k.Z()))*(endz-pos_k.Z())+pos_k.Y()-endy)<50||std::abs(((end_k.Y()-pos_k.Y())/(end_k.Z()-pos_k.Z()))*(startz-pos_k.Z())+pos_k.Y()-starty)<50)&&(std::abs(endcosx*dir_pos_k.X()+endcosy*dir_pos_k.Y()+endcosz*dir_pos_k.Z())>0.998||std::abs(startcosx*dir_pos_k.X()+startcosy*dir_pos_k.Y()+startcosz*dir_pos_k.Z())>0.998||std::abs(endcosx*dir_end_k.X()+endcosy*dir_end_k.Y()+endcosz*dir_end_k.Z())>0.998||std::abs(startcosx*dir_end_k.X()+startcosy*dir_end_k.Y()+startcosz*dir_end_k.Z())>0.998)) break;
	count++;
	 
      }
      if(!(count==NTracks-1)|| tracklength<100) continue;
      int geant_trkID=particle->TrackId();	   
      int n_daughters_electron=0;
      int n_daughters_positron=0;
      if(abs(pdg1)==13){
	const sim::ParticleList& plist = pi_serv->ParticleList();
	sim::ParticleList::const_iterator itPart=plist.begin(),pend=plist.end();
	for(size_t iPart = 0;(iPart<plist.size()) && (itPart!=pend); ++iPart){
	  const simb::MCParticle* pPart = (itPart++)->second;
	  int pdgcode = pPart->PdgCode();
	  // if(abs(pdgcode)==11) std::cout<<pdgcode<<"  "<<pPart->Process()<<"  "<<pPart->Mother()<<" "<<geant_trkID<<std::endl;
	  if(pPart->Mother()==geant_trkID && pdgcode==11)  n_daughters_electron++;
	  if(pPart->Mother()==geant_trkID && pdgcode==-11) n_daughters_positron++;
	}
      }
      /*************************************filling the values****************************/
      cross_trks++;

      auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
      EndPointx_corrected[cross_trks-1]=particle->EndPosition()[0]-SCE->GetPosOffsets(geo::Point_t(particle->EndPosition()[0],particle->EndPosition()[1],particle->EndPosition()[2])).X();
      EndPointy_corrected[cross_trks-1]=particle->EndPosition()[1]+SCE->GetPosOffsets(geo::Point_t(particle->EndPosition()[0],particle->EndPosition()[1],particle->EndPosition()[2])).Y();
      EndPointz_corrected[cross_trks-1]=particle->EndPosition()[2]+SCE->GetPosOffsets(geo::Point_t(particle->EndPosition()[0],particle->EndPosition()[1],particle->EndPosition()[2])).Z();
      origin[cross_trks-1]=mc->Origin();
      daughter_electrons[cross_trks-1]=n_daughters_electron++;
      daughter_positrons[cross_trks-1]=n_daughters_positron++;
      EndE[cross_trks-1]=particle->EndE();
      pdg[cross_trks-1]=particle->PdgCode();
      TrackId_geant[cross_trks-1]=particle->TrackId();
      EndPointx[cross_trks-1]=particle->EndPosition()[0];
      EndPointy[cross_trks-1]=particle->EndPosition()[1];
      EndPointz[cross_trks-1]=particle->EndPosition()[2];
      StartPointx[cross_trks-1]=particle->Vx();
      StartPointy[cross_trks-1]=particle->Vy();
      StartPointz[cross_trks-1]=particle->Vz();
      Total_energy[cross_trks-1]=tote;
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
      peakT_max[cross_trks-1]=max;
      peakT_min[cross_trks-1]=min;
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
	  trkpitch[cross_trks-1][planenum][iHit]=(calos[ical]->TrkPitchVec())[iHit];
	} // loop over iHit..
      } // loop over ical 2nd time...
    } // loop over trks...
    // stopping_muons->Fill(stopping_trks);
    fEventTree->Fill();
  } // end of analyze function
	   
  /////////////////// Defintion of reset function ///////////
  void dEdxcalibration::reset(){
    run = -99999;
    subrun = -99999;
    event = -99999;
    evttime = -99999;
    cross_trks = -99999;
    true_stopping_muons=-99999;
    mc_stopping20=-99999;
    mc_stopping50=-99999;
    all_trks = -99999;
    year_month_date=-99999;
    hour_min_sec=-99999;
    for(int i=0; i<kMaxTracks; i++){
      trackthetaxz[i]=-99999;
      trackthetayz[i]=-99999;
      trkstartx[i]=-99999;
      trkstarty[i]=-99999;
      trkstartz[i]=-99999;
      trkendx[i]=-99999;
      trkendy[i]=-99999;
      trkendz[i]=-99999;
      trklen[i]=-99999;
      TrkID[i]=-99999;
      TrackId_geant[i]=-99999;
      origin[i]=-99999;
      EndE[i]=-99999;
      pdg[i]=-99999;
      Total_energy[i]=-999999;
      daughter_electrons[i]=-99999;
      daughter_positrons[i]=-99999;
      EndPointx[i]=-99999;
      EndPointy[i]=-99999;
      EndPointz[i]=-99999;
      EndPointx_corrected[i]=-99999;
      EndPointy_corrected[i]=-99999;
      EndPointz_corrected[i]=-99999;
      StartPointx[i]=-99999;
      StartPointy[i]=-99999;
      StartPointz[i]=-99999;
      pdg[i]=-99999;
      EndE[i]=-99999;
      peakT_max[i]=-99999;
      peakT_min[i]=-99999;
      xprojectedlen[i]=-99999;
      trkstartcosxyz[i][0]=-99999;
      trkstartcosxyz[i][1]=-99999;
      trkstartcosxyz[i][2]=-99999; 
      trkendcosxyz[i][0]=-99999;
      trkendcosxyz[i][1]=-99999;
      trkendcosxyz[i][2]=-99999;
      ntrkhits[i][0] = -99999;
      ntrkhits[i][1] = -99999;
      ntrkhits[i][2] = -99999;
      for(int j=0; j<3; j++){
	for(int k=0; k<3000; k++){
	  trkdqdx[i][j][k]=-99999;
	  trkdedx[i][j][k]=-99999;
	  trkresrange[i][j][k]=-99999;
	  trkhitx[i][j][k]=-99999;
	  trkhity[i][j][k]=-99999;
	  trkhitz[i][j][k]=-99999;
	  trkpitch[i][j][k]=-99999;
	}
      }
    }
  }
  //////////////////////// End of definition ///////////////	
	  
  DEFINE_ART_MODULE(dEdxcalibration)
}
