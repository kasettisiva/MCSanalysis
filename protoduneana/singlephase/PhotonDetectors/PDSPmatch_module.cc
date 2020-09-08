////////////////////////////////////////////////////////////////////////
// Class:       PDSPmatch
// Plugin Type: analyzer (art v3_00_00)
// File:        PDSPmatch_module.cc
//
// Generated at Thu Feb 14 21:47:27 2019 by Bryan Joseph Ramson using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////

#include <numeric>
#include <sstream>
#include <string.h>
#include <bitset>
#include <vector>

#include "larana/OpticalDetector/OpFlashAlg.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
//#include "art/Framework/Services/Optional/TFileService.h"


#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "dune/Protodune/singlephase/CRT/data/CRTTrigger.h"
//#include "dune/Geometry/ProtoDUNESPCRTSorter.h"

#include "dune-raw-data/Overlays/CRTFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "canvas/Persistency/Common/FindManyP.h"

namespace pdsp {
  class PDSPmatch;
}


class pdsp::PDSPmatch : public art::EDAnalyzer {
public:
  explicit PDSPmatch(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPmatch(PDSPmatch const&) = delete;
  PDSPmatch(PDSPmatch&&) = delete;
  PDSPmatch& operator=(PDSPmatch const&) = delete;
  PDSPmatch& operator=(PDSPmatch&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  const art::InputTag fTimeLabel; //Label of timing products
  const art::InputTag fCTBLabel; //Label of pdsctbdata products
  const art::InputTag fCRTLabel; //Label of crt products
  const art::InputTag fOpHitLabel; //Label of ophits after rawdecoding
  const art::InputTag fPandoLabel; //Label for Pandora Tracks
  const art::InputTag fPFParListLabel; // Label for PF Particle containers
  
  const int64_t fCRTCTBOffset; //CRT offset to global trigger in 20mus ticks
  const uint64_t fCRTWindow; //coincidence window for CRT triggers in 20mus ticks
  const bool fMC; //MC flag
  
  TTree *fTree;
  int run;
  int subrun;
  int event;
  int64_t timing_time;
  int64_t fCTB_time;
  uint32_t fCTBChan;
  
  //  TH1F *CTBCRTMap[32];
  //TH1F *CTBPairPE[16];
  //TH1F *CosTrigger[16];
  //TH1F *CRTPairPE[16];
  //TH2F *PandaTrackProUp[16];
  //TH2F *PandaTrackProDown[16];
  
  //TH1F *PMTrack[16];
  //TH1F *CRTPairPE[16];
  //TH1F *trig

  std::vector<uint32_t> fCRTChan;
  std::vector<int64_t> fPDS_time, fPando_time, fCRT_time;
  
  std::vector<double> fOpChan, fPE;//, ft0, fx0, fy0, fz0;
  
  std::vector<double> fTrkStartx_Pando, fTrkStarty_Pando, fTrkStartz_Pando, fTrkEndx_Pando, fTrkEndy_Pando, fTrkEndz_Pando;
  std::vector<double> fTrkTheta_Pando, fTrkPhi_Pando, fTrkUSPos1x_Pando, fTrkUSPos1y_Pando, fTrkUSPos2x_Pando, fTrkUSPos2y_Pando, fTrkDSPosx_Pando, fTrkDSPosy_Pando; 
  std::vector<double> fTrigCos_Pando;
  std::vector<uint32_t> fTrkProUS_Pando, fTrkProDS_Pando;
  std::vector<double> fTrkUSPixelProx_Pando, fTrkUSPixelProy_Pando; 
  std::vector<double> fTrkDSPixelProx_Pando, fTrkDSPixelProy_Pando; 
  
  std::vector<TVector3> PandoTracks, PMTracks;
  // Declare member data here.
  
  TVector3 CTBPixelCtr[32];
};


pdsp::PDSPmatch::PDSPmatch(fhicl::ParameterSet const& p)
  : 
  EDAnalyzer(p), 
  fTimeLabel(p.get<art::InputTag>("TimingLabel")),
  fCTBLabel(p.get<art::InputTag>("CTBLabel")),
  fCRTLabel(p.get<art::InputTag>("CRTLabel")),
  fOpHitLabel(p.get<art::InputTag>("OpHitLabel")),
  fPandoLabel(p.get<art::InputTag>("PandoraLabel")),
  fPFParListLabel(p.get<art::InputTag>("PFParListLabel")),
  fCRTCTBOffset(p.get<int64_t>("CRTCTBOffset")),
  fCRTWindow(p.get<uint64_t>("CRTWindow")),
  fMC(p.get<bool>("MC"))
  //: EDAnalyzer{p}  // ,
  //More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  consumes<std::vector<raw::RDTimeStamp>>(fTimeLabel);
  consumes<std::vector<raw::ctb::pdspctb>>(fCTBLabel);
  consumes<std::vector<CRT::Trigger>>(fCRTLabel);
  consumes<std::vector<recob::OpHit>>(fOpHitLabel);
  consumes<std::vector<recob::Track>>(fPandoLabel);
  consumes<std::vector<recob::PFParticle>>(fPFParListLabel);
}

void pdsp::PDSPmatch::beginJob(){
  art::ServiceHandle<art::TFileService> tFileService;
  CTBPixelCtr[0].SetXYZ(-111.494,580.249,-972.204);
  CTBPixelCtr[1].SetXYZ(58.5057,580.249,-972.204);
  CTBPixelCtr[2].SetXYZ(307.506,580.249,-267.204);
  CTBPixelCtr[3].SetXYZ(477.506,580.249,-267.204);
  CTBPixelCtr[4].SetXYZ(477.506,410.249,-267.204);
  CTBPixelCtr[5].SetXYZ(477.506,235.249,-267.204);
  CTBPixelCtr[6].SetXYZ(477.506,65.2494,-267.204);
  CTBPixelCtr[7].SetXYZ(307.506,65.2494,-267.204);
  CTBPixelCtr[8].SetXYZ(58.5057,65.2494,-972.204);
  CTBPixelCtr[9].SetXYZ(-111.494,65.2494,-972.204);
  CTBPixelCtr[10].SetXYZ(-111.494,235.249,-972.204);
  CTBPixelCtr[11].SetXYZ(-111.494,410.249,-972.204);
  CTBPixelCtr[12].SetXYZ(58.5057,410.249,-972.204);
  CTBPixelCtr[13].SetXYZ(307.506,410.249,-267.204);
  CTBPixelCtr[14].SetXYZ(307.506,235.249,-267.204);
  CTBPixelCtr[15].SetXYZ(58.5057,235.249,-972.204);
  CTBPixelCtr[16].SetXYZ(293.506,580.249,1085.55);
  CTBPixelCtr[17].SetXYZ(123.506,580.249,1085.55);
  CTBPixelCtr[18].SetXYZ(-51.4993,580.249,1085.55);
  CTBPixelCtr[19].SetXYZ(-221.494,580.249,1085.55);
  CTBPixelCtr[20].SetXYZ(-221.494,410.249,1085.55);
  CTBPixelCtr[21].SetXYZ(-221.494,235.249,1085.55);
  CTBPixelCtr[22].SetXYZ(-221.494,65.2494,1085.55);
  CTBPixelCtr[23].SetXYZ(-51.4943,65.2494,1085.55);
  CTBPixelCtr[24].SetXYZ(123.506,65.2494,1085.55);
  CTBPixelCtr[25].SetXYZ(293.506,65.2494,1085.55);
  CTBPixelCtr[26].SetXYZ(293.506,235.249,1085.55);
  CTBPixelCtr[27].SetXYZ(293.506,410.249,1085.55);
  CTBPixelCtr[28].SetXYZ(123.506,410.249,1085.55);
  CTBPixelCtr[29].SetXYZ(-51.4993,410.249,1085.55);
  CTBPixelCtr[30].SetXYZ(-51.4943,235.249,1085.55);
  CTBPixelCtr[31].SetXYZ(123.506,235.249,1085.55);
  
  for(int h=0;h<32;h++){
    if(h==0 || h==1 || h==11 || h==12 || h ==10 || h==15 || h==9 || h==8) CTBPixelCtr[h].SetY(CTBPixelCtr[h].Y()-46.0); //US Jura/Beam Left
    if(h==2 || h==3 || h==13 || h==4 || h ==14 || h==5 || h==7 || h==6) CTBPixelCtr[h].SetY(CTBPixelCtr[h].Y()-47.5); //US Saleve/Beam Right
    if(h==16 || h==17 || h==27 || h==28 || h ==26 || h==31 || h==25 || h==24) CTBPixelCtr[h].SetY(CTBPixelCtr[h].Y()-151.0); //DS Saleve/Beam Right
    if(h==18 || h==19 || h==29 || h==20 || h ==21 || h==30 || h==22 || h==23) CTBPixelCtr[h].SetY(CTBPixelCtr[h].Y()-150.0); //DS Jura/Beam Left
  }								    
  art::ServiceHandle < geo::Geometry > geom;
 
  fTree = tFileService->make<TTree>("ProtoDUNE_Evt_Match","ProtoDUNE Matched Event Tree");
  
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("timing_time",&timing_time,"timing_time/L");
  
  fTree->Branch("fCRT_time",&fCRT_time);
  fTree->Branch("fCRTChan", &fCRTChan);
  fTree->Branch("fCTB_time",&fCTB_time,"fCTB_time/L");
  fTree->Branch("fCTBChan", &fCTBChan,"fCTBChan/i");
  
  fTree->Branch("fPDS_time",&fPDS_time);
  fTree->Branch("fOpChan", &fOpChan);
  fTree->Branch("fPE",&fPE);
 
  fTree->Branch("fPando_time",&fPando_time);
  fTree->Branch("fTrkStartx_Pando",&fTrkStartx_Pando);
  fTree->Branch("fTrkStarty_Pando",&fTrkStarty_Pando);
  fTree->Branch("fTrkStartz_Pando",&fTrkStartz_Pando);
  fTree->Branch("fTrkEndx_Pando",&fTrkEndx_Pando);
  fTree->Branch("fTrkEndy_Pando",&fTrkEndy_Pando);
  fTree->Branch("fTrkEndz_Pando",&fTrkEndz_Pando);
  

}

void pdsp::PDSPmatch::analyze(art::Event const& e){
  run = e.run();
  subrun = e.subRun();
  event = e.id().event();
  
  fCRT_time.clear();
  fCRTChan.clear();
    
  fPDS_time.clear();
  fOpChan.clear();
  fPE.clear();
  
  fPando_time.clear();
  fTrkStartx_Pando.clear();
  fTrkStarty_Pando.clear();
  fTrkStartz_Pando.clear();
  fTrkEndx_Pando.clear();
  fTrkEndy_Pando.clear();
  fTrkEndz_Pando.clear();
  fTrkTheta_Pando.clear();
  fTrkPhi_Pando.clear();
  
  int trigger = -1;
  const auto timeStamps = e.getValidHandle<std::vector<raw::RDTimeStamp>>(fTimeLabel);
  if(timeStamps.isValid() && timeStamps->size() == 1){
    const raw::RDTimeStamp& timeStamp = timeStamps->at(0);
    trigger = timeStamp.GetFlags();
    timing_time = timeStamp.GetTimeStamp();
  }
  else {
    mf::LogWarning("Empty Event TimeStamp") << "Invalid Event TimeStamp. Skipping. \n";
    return;
  }
    
  if (trigger !=13) return;
  
  const auto crtHandle = e.getValidHandle<std::vector<CRT::Trigger>>(fCRTLabel);
  if(crtHandle->empty()){
    mf::LogWarning("Empty CRT Fragment") << "Empty CRT fragments for this event. Skipping. \n";
    return;
  }
    
  const auto OpHitHandle = e.getValidHandle<std::vector<recob::OpHit>>(fOpHitLabel);
  if(OpHitHandle->empty()){
    mf::LogWarning("Empty OpHit Object") << "Empty OpHit Object. Error in retrieval. Skipping. \n";
    return;
  }
    
  const auto ctbHandle = e.getValidHandle<std::vector<raw::ctb::pdspctb>>(fCTBLabel);
  if(ctbHandle->empty()){
    mf::LogWarning("Empty CTB Fragment") << "Empty CTB fragment for this event. Skipping. \n";
    return;
  }
  
  if(ctbHandle->size() > 1){
    mf::LogWarning("Multiple CTB Triggers") << "Found " << ctbHandle->size() << " CTB data products. Skipping. \n";
    return;
  }
  
  const auto& ctbStatus = ctbHandle->front();
  const auto statuses = ctbStatus.GetChStatusAfterHLTs();
  
  bool pat = false;
  for(const auto& status: statuses){
    pat = false;
    
    fCRT_time.clear();
    fCRTChan.clear();
    
    fPDS_time.clear();
    fOpChan.clear();
    fPE.clear();
    
    fPando_time.clear();
    fTrkStartx_Pando.clear();
    fTrkStarty_Pando.clear();
    fTrkStartz_Pando.clear();
    fTrkEndx_Pando.clear();
    fTrkEndy_Pando.clear();
    fTrkEndz_Pando.clear();

    if(status.crt != 0) {
      fCTB_time = status.timestamp;
      fCTBChan = status.crt;
      if(!crtHandle->empty()){ 
	std::vector<CRT::Trigger> inWindow(crtHandle->size());
	const int64_t ctbetime = status.timestamp;
	const auto newEnd = std::copy_if(crtHandle->begin(), crtHandle->end(), inWindow.begin(),[this,ctbetime](const auto& trigger){
	    return (uint64_t)(labs(ctbetime - ((int64_t)(trigger.Timestamp()) - fCRTCTBOffset))) < fCRTWindow; });
	inWindow.resize(std::distance(inWindow.begin(), newEnd));   
	for(size_t k=0;k<inWindow.size();++k){
	  fCRT_time.push_back(inWindow[k].Timestamp());
	  fCRTChan.push_back(inWindow[k].Channel());
	}
      }
      for(const auto& OpHit: *OpHitHandle){
	if(OpHit.PE() < 10.0 || OpHit.PE() > 1000.00) continue;
	fPDS_time.push_back((OpHit.PeakTime()/3));
	fOpChan.push_back(OpHit.OpChannel());
	fPE.push_back(OpHit.PE());
      }
      
      art::Handle<std::vector<recob::Track>> PandoTrkHandle;
      std::vector<art::Ptr<recob::Track>> PandoTrk;
      if(e.getByLabel(fPandoLabel,PandoTrkHandle)) art::fill_ptr_vector(PandoTrk,PandoTrkHandle);
      else {
	mf::LogWarning("Empty PandoTrk Fragment") << "Empty PandoTrk Vector for this event. Skipping. \n";
	return;
      }
      
      art::Handle<std::vector<recob::PFParticle>> PFParListHandle;
      if(!(e.getByLabel(fPFParListLabel,PFParListHandle))){;
	mf::LogWarning("Empty PFParticle Vector") << "Empty PFParticle Vector for this event. Skipping. \n";
	return;
      }
      art::FindManyP<recob::PFParticle> PFPar(PandoTrkHandle,e,fPandoLabel);
      art::FindManyP<anab::T0> PFT0(PFParListHandle,e,fPFParListLabel);

      for(size_t p = 0;p<PandoTrk.size();++p){
	auto & Trk = PandoTrk[p];
	if(!((Trk->Vertex().Z() < 40) || (Trk->End().Z() < 40))) continue;
	if(!((Trk->Vertex().Z() > 660) || (Trk->End().Z() > 660))) continue;
	fTrkStartx_Pando.push_back(Trk->Vertex().X());
	fTrkStarty_Pando.push_back(Trk->Vertex().Y());
	fTrkStartz_Pando.push_back(Trk->Vertex().Z());
	fTrkEndx_Pando.push_back(Trk->End().X());
	fTrkEndy_Pando.push_back(Trk->End().Y());
	fTrkEndz_Pando.push_back(Trk->End().Z());
	double t0temp = 0;
	auto &PFPS = PFPar.at(Trk.key());
	if(!PFPS.empty()){
	  auto &T0S = PFT0.at(PFPS[0].key());
	  if(!T0S.empty()){
	    t0temp = T0S[0]->Time();
	  }
	}
	fPando_time.push_back(t0temp);
      }
      if(PandoTrk.size() > 0) pat = true;

      if(pat) fTree->Fill();
    }
  }
}
 DEFINE_ART_MODULE(pdsp::PDSPmatch)
