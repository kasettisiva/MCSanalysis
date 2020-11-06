////////////////////////////////////////////////////////////////////////
// Class:       CRTMatchTrackCaloAna
// Plugin Type: analyzer (art v3_02_06)
// File:        CRTMatchTrackCaloAna_module.cc
// This code is an adaptation of Tingjun Yang's T0 code with adjustments so it pulls a crt matched track and the TPC dE/dx information
// Generated at Tue Jul 30 21:43:10 2019 by Tingjun Yang using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "protoduneana/Utilities/ProtoDUNECalibration.h"
#include "TTree.h"

#include <vector>

namespace pdsp {
  class CRTMatchTrackCaloAna;
}


class pdsp::CRTMatchTrackCaloAna : public art::EDAnalyzer {
public:
  explicit CRTMatchTrackCaloAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTMatchTrackCaloAna(CRTMatchTrackCaloAna const&) = delete;
  CRTMatchTrackCaloAna(CRTMatchTrackCaloAna&&) = delete;
  CRTMatchTrackCaloAna& operator=(CRTMatchTrackCaloAna const&) = delete;
  CRTMatchTrackCaloAna& operator=(CRTMatchTrackCaloAna&&) = delete;
  void beginJob() override;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  protoana::ProtoDUNECalibration calibration;

  // Declare member data here.
  TTree *crtCalo;
  int run;
  int event;
  std::vector<int> trackid;
  std::vector<double> t0crt2;
  std::vector<double> t0pandora;
  std::vector<double> t0anodep;
  std::vector<double> t0truth;
  std::vector<double> trackstartx;
  std::vector<double> trackstarty;
  std::vector<double> trackstartz;
  std::vector<double> trackendx;
  std::vector<double> trackendy;
  std::vector<double> trackendz;
  std::vector<double> trackstartx_sce;
  std::vector<double> trackstarty_sce;
  std::vector<double> trackstartz_sce;
  std::vector<double> trackendx_sce;
  std::vector<double> trackendy_sce;
  std::vector<double> trackendz_sce;

  std::vector<double> crt2x0;
  std::vector<double> crt2y0;
  std::vector<double> crt2z0;
  std::vector<double> crt2x1;
  std::vector<double> crt2y1;
  std::vector<double> crt2z1;
  std::vector<double> trkhitx;
  std::vector<double> trkhity;
  std::vector<double> trkhitz;
  std::vector<double> WireID;
  std::vector<double> TPCID;
    std::vector<double>  trkpitch;
    std::vector<double>  trkdqdx;
    std::vector<double>  trkdedx;
    std::vector<double>  trkdedx_cali;
    std::vector<double>  trkresrange;
  double rdtimestamp_evt;
  std::vector<double> rdtimestamp_digits;
  fhicl::ParameterSet CalibrationPars;
    std::string fCalorimetryModuleLabel;
};


pdsp::CRTMatchTrackCaloAna::CRTMatchTrackCaloAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},  // ,
  // More initializers here.
  CalibrationPars(p.get<fhicl::ParameterSet>("CalorimetryParameters")),
    fCalorimetryModuleLabel (p.get<std::string>("CalorimetryModuleLabel"))
{
 calibration = protoana::ProtoDUNECalibration( CalibrationPars);
}

void pdsp::CRTMatchTrackCaloAna::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  run = e.run();
  event = e.id().event();

  rdtimestamp_evt = -1;
  rdtimestamp_digits.clear();

  //Services
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  //Space charge service
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  
   // Reconstruciton information
  art::Handle < std::vector < recob::Track > > trackListHandle;
  std::vector < art::Ptr < recob::Track > > trackList;
  if (e.getByLabel("pandoraTrack", trackListHandle)) {
    art::fill_ptr_vector(trackList, trackListHandle);
  }
  else return;

  //Get hits associated with track
  art::FindManyP < recob::Hit > hitsFromTrack(trackListHandle, e, "pandoraTrack");

  //Get PFParticles
  art::Handle< std::vector<recob::PFParticle> > pfpListHandle;
  e.getByLabel("pandora", pfpListHandle);
  
  //Get Track-PFParticle association
  art::FindManyP<recob::PFParticle> fmpfp(trackListHandle, e, "pandoraTrack");

  //Get Pandora T0-PFParticle association
  art::FindManyP<anab::T0> fmt0pandora(pfpListHandle, e, "pandora");

  //Get Anode Piercing T0
  art::FindManyP<anab::T0> fmt0anodepiercer(pfpListHandle, e, "anodepiercerst0");
  
  //Get 2-CRT T0
  art::FindManyP<anab::T0> fmt0crt2(trackListHandle, e, "crtreco");
  art::FindManyP<anab::CosmicTag> fmctcrt2(trackListHandle, e, "crtreco");
    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, e, "pandoracalo");


    art::Handle< std::vector<recob::Hit> > hitListHandle; // to get information about the hits
    std::vector<art::Ptr<recob::Hit>> hitlist;
    if(e.getByLabel("hitpdune", hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);


  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e, clockData);

  for (auto const& track : trackList){


	trackid.clear();
	t0crt2.clear();
	trackstartx.clear();
	trackstarty.clear();
	trackstartz.clear();
	trackendx.clear();
	trackendy.clear();
	trackendz.clear();
	trackstartx_sce.clear();
	trackstarty_sce.clear();
	trackstartz_sce.clear();
	trackendx_sce.clear();
	trackendy_sce.clear();
	trackendz_sce.clear();

	crt2x0.clear();
	crt2y0.clear();
	crt2z0.clear();
	crt2x1.clear();
	crt2y1.clear();
	crt2z1.clear();
	trkpitch.clear();
	WireID.clear();
	TPCID.clear();
	trkdedx_cali.clear();
	trkdedx.clear();
	trkpitch.clear();
	trkhitx.clear(); trkhity.clear(); trkhitz.clear();
    int this_trackid = track.key();
    double this_t0crt2 = -DBL_MAX;
    double this_t0pandora = -DBL_MAX;
    double this_t0anodep = -DBL_MAX;
    double this_t0truth = -DBL_MAX;
    double this_trackstartx = -DBL_MAX;
    double this_trackstarty = -DBL_MAX;
    double this_trackstartz = -DBL_MAX;
    double this_trackendx = -DBL_MAX;
    double this_trackendy = -DBL_MAX;
    double this_trackendz = -DBL_MAX;
    double this_trackstartx_sce = -DBL_MAX;
    double this_trackstarty_sce = -DBL_MAX;
    double this_trackstartz_sce = -DBL_MAX;
    double this_trackendx_sce = -DBL_MAX;
    double this_trackendy_sce = -DBL_MAX;
    double this_trackendz_sce = -DBL_MAX;
    double this_crt2x0 = -DBL_MAX;
    double this_crt2y0 = -DBL_MAX;
    double this_crt2z0 = -DBL_MAX;
    double this_crt2x1 = -DBL_MAX;
    double this_crt2y1 = -DBL_MAX;
    double this_crt2z1 = -DBL_MAX;
    if (!fmctcrt2.isValid()) continue;
    if (fmt0crt2.isValid()){
      auto const& vt0crt2 = fmt0crt2.at(track.key());
      if (!vt0crt2.empty()) this_t0crt2 = vt0crt2[0]->Time();
    }


    if (fmctcrt2.isValid()){
      auto const& vctcrt2 = fmctcrt2.at(track.key());
      if (!vctcrt2.empty()){
        this_crt2x0 = vctcrt2[0]->EndPoint1()[0];
        this_crt2y0 = vctcrt2[0]->EndPoint1()[1];
        this_crt2z0 = vctcrt2[0]->EndPoint1()[2];
        this_crt2x1 = vctcrt2[0]->EndPoint2()[0];
        this_crt2y1 = vctcrt2[0]->EndPoint2()[1];
        this_crt2z1 = vctcrt2[0]->EndPoint2()[2];
      }
    }



    if (this_t0crt2 > -DBL_MAX){

      auto const & allHits = hitsFromTrack.at(track.key());

      
      this_trackstartx = track->Vertex().X();
      this_trackstarty = track->Vertex().Y();
      this_trackstartz = track->Vertex().Z();
      this_trackendx = track->End().X();
      this_trackendy = track->End().Y();
      this_trackendz = track->End().Z();

      if (std::abs(this_t0pandora+DBL_MAX)<1e-10){
        //no pandora t0 found, correct for t0
        double ticksOffset = 0;
        if (this_t0crt2 > -DBL_MAX) ticksOffset = this_t0crt2/500.+detProp.GetXTicksOffset(allHits[0]->WireID());
        else if (this_t0anodep > -DBL_MAX) ticksOffset = this_t0anodep/500.+detProp.GetXTicksOffset(allHits[0]->WireID());
        double xOffset = detProp.ConvertTicksToX(ticksOffset,allHits[0]->WireID());
        this_trackstartx -= xOffset;
        this_trackendx -= xOffset;
      }

      auto const & posOffsets = SCE->GetCalPosOffsets(geo::Point_t(this_trackstartx, this_trackstarty, this_trackstartz), allHits[0]->WireID().TPC);
      std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(track.key());


	std::vector< float > calis = calibration.GetCalibratedCalorimetry(*track, e, "pandoraTrack", fCalorimetryModuleLabel, 2);
      double trkHitX, trkHitY, trkHitZ;
      for(size_t ical = 0; ical<calos.size(); ++ical){
	if (calos[ical]->PlaneID().Plane!=2) continue;
	
	const size_t NHits=calos[ical]->dEdx().size();
	for(size_t iHit = 0; iHit < NHits; ++iHit){
		   
	  const auto& TrkPos = (calos[ical] -> XYZ()[iHit]);

          trkdqdx.push_back(calos[ical]->dQdx()[iHit]);
	  size_t tpIndex=(calos[ical]->TpIndices()[iHit]);
	  //trkhitIntegral=hitlist[tpIndex]->Integral();
	  trkdedx_cali.push_back(calis[iHit]);
	  trkdedx.push_back(calos[ical]->dEdx()[iHit]);
	  trkHitX=TrkPos.X();
	  trkHitY=TrkPos.Y();
	  trkHitZ=TrkPos.Z();
	  //std::cout<<"Trk Hit Z:"<<iHit<<','<<trkHitZ<<std::endl;
          trkpitch.push_back(calos[ical]->TrkPitchVec()[iHit]);
          trkhitx.push_back(trkHitX);
	  trkhity.push_back(trkHitY);
	  trkhitz.push_back(trkHitZ);
	  WireID.push_back(hitlist[tpIndex]->WireID().Wire);
	  TPCID.push_back(hitlist[tpIndex]->WireID().TPC); 
	}
 
      } 
	

      this_trackstartx_sce = this_trackstartx - posOffsets.X();
      this_trackstarty_sce = this_trackstarty + posOffsets.Y();
      this_trackstartz_sce = this_trackstartz + posOffsets.Z();
      this_trackendx_sce = this_trackendx - posOffsets.X();
      this_trackendy_sce = this_trackendy + posOffsets.Y();
      this_trackendz_sce = this_trackendz + posOffsets.Z();

      trackid.push_back(this_trackid);
      t0crt2.push_back(this_t0crt2);
      t0truth.push_back(this_t0truth);
      trackstartx.push_back(this_trackstartx);
      trackstarty.push_back(this_trackstarty);
      trackstartz.push_back(this_trackstartz);
      trackendx.push_back(this_trackendx);
      trackendy.push_back(this_trackendy);
      trackendz.push_back(this_trackendz);
      trackstartx_sce.push_back(this_trackstartx_sce);
      trackstarty_sce.push_back(this_trackstarty_sce);
      trackstartz_sce.push_back(this_trackstartz_sce);
      trackendx_sce.push_back(this_trackendx_sce);
      trackendy_sce.push_back(this_trackendy_sce);
      trackendz_sce.push_back(this_trackendz_sce);
      crt2x0.push_back(this_crt2x0);
      crt2y0.push_back(this_crt2y0);
      crt2z0.push_back(this_crt2z0);
      crt2x1.push_back(this_crt2x1);
      crt2y1.push_back(this_crt2y1);
      crt2z1.push_back(this_crt2z1);
      if (!trackid.empty()) crtCalo->Fill();
    }
  }

 
}

void pdsp::CRTMatchTrackCaloAna::beginJob() {
  art::ServiceHandle<art::TFileService> fileServiceHandle;
  crtCalo = fileServiceHandle->make<TTree>("crtCalo", "t0 info");
  crtCalo->Branch("run", &run, "run/I");
  crtCalo->Branch("event", &event, "event/I");
  crtCalo->Branch("trackid", &trackid);
  crtCalo->Branch("t0crt2", &t0crt2);
  crtCalo->Branch("t0truth", &t0truth);
  crtCalo->Branch("trackstartx",&trackstartx);
  crtCalo->Branch("trackstarty",&trackstarty);
  crtCalo->Branch("trackstartz",&trackstartz);
  crtCalo->Branch("trackendx",&trackendx);
  crtCalo->Branch("trackendy",&trackendy);
  crtCalo->Branch("trackendz",&trackendz);
  crtCalo->Branch("trackstartx_sce",&trackstartx_sce);
  crtCalo->Branch("trackstarty_sce",&trackstarty_sce);
  crtCalo->Branch("trackstartz_sce",&trackstartz_sce);
  crtCalo->Branch("trackendx_sce",&trackendx_sce);
  crtCalo->Branch("trackendy_sce",&trackendy_sce);
  crtCalo->Branch("trackendz_sce",&trackendz_sce);
  crtCalo->Branch("crt2x0",&crt2x0);
  crtCalo->Branch("crt2y0",&crt2y0);
  crtCalo->Branch("crt2z0",&crt2z0);
  crtCalo->Branch("crt2x1",&crt2x1);
  crtCalo->Branch("crt2y1",&crt2y1);
  crtCalo->Branch("crt2z1",&crt2z1);


  crtCalo->Branch("trkdedx_cali",&trkdedx_cali);
  crtCalo->Branch("trkdedx",&trkdedx);
  crtCalo->Branch("trkpitch",&trkpitch);
  crtCalo->Branch("trkhitx",&trkhitx);
  crtCalo->Branch("trkhity",&trkhity);
  crtCalo->Branch("trkhitz",&trkhitz);
  crtCalo->Branch("TPCID",&TPCID);
  crtCalo->Branch("WireID",&WireID);
}
DEFINE_ART_MODULE(pdsp::CRTMatchTrackCaloAna)
