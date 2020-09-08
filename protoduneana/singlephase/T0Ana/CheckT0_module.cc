////////////////////////////////////////////////////////////////////////
// Class:       CheckT0
// Plugin Type: analyzer (art v3_02_06)
// File:        CheckT0_module.cc
//
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

#include "lardata/Utilities/AssociationUtil.h"
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

#include "TTree.h"

#include <vector>

namespace pdsp {
  class CheckT0;
}


class pdsp::CheckT0 : public art::EDAnalyzer {
public:
  explicit CheckT0(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CheckT0(CheckT0 const&) = delete;
  CheckT0(CheckT0&&) = delete;
  CheckT0& operator=(CheckT0 const&) = delete;
  CheckT0& operator=(CheckT0&&) = delete;
  void beginJob() override;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  TTree *t0tree;
  int run;
  int event;
  std::vector<int> trackid;
  std::vector<double> t0crt2;
  std::vector<double> t0crt1;
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
  std::vector<double> crt1x;
  std::vector<double> crt1y;
  std::vector<double> crt1z;
  std::vector<double> crt2x0;
  std::vector<double> crt2y0;
  std::vector<double> crt2z0;
  std::vector<double> crt2x1;
  std::vector<double> crt2y1;
  std::vector<double> crt2z1;
  double rdtimestamp_evt;
  std::vector<double> rdtimestamp_digits;

};


pdsp::CheckT0::CheckT0(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

}

void pdsp::CheckT0::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  run = e.run();
  event = e.id().event();
  trackid.clear();
  t0crt2.clear();
  t0crt1.clear();
  t0pandora.clear();
  t0anodep.clear();
  t0truth.clear();
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
  crt1x.clear();
  crt1y.clear();
  crt1z.clear();
  crt2x0.clear();
  crt2y0.clear();
  crt2z0.clear();
  crt2x1.clear();
  crt2y1.clear();
  crt2z1.clear();
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

  //Get 1-CRT T0
  art::FindManyP<anab::T0> fmt0crt1(trackListHandle, e, "crttag");
  art::FindManyP<anab::CosmicTag> fmctcrt1(trackListHandle, e, "crttag");
  
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e, clockData);

  for (auto const& track : trackList){
    int this_trackid = track.key();
    double this_t0crt2 = -DBL_MAX;
    double this_t0crt1 = -DBL_MAX;
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
    double this_crt1x = -DBL_MAX;
    double this_crt1y = -DBL_MAX;
    double this_crt1z = -DBL_MAX;
    double this_crt2x0 = -DBL_MAX;
    double this_crt2y0 = -DBL_MAX;
    double this_crt2z0 = -DBL_MAX;
    double this_crt2x1 = -DBL_MAX;
    double this_crt2y1 = -DBL_MAX;
    double this_crt2z1 = -DBL_MAX;
    
    if (fmt0crt2.isValid()){
      auto const& vt0crt2 = fmt0crt2.at(track.key());
      if (!vt0crt2.empty()) this_t0crt2 = vt0crt2[0]->Time();
    }

    if (fmt0crt1.isValid()){
      auto const& vt0crt1 = fmt0crt1.at(track.key());
      if (!vt0crt1.empty()) this_t0crt1 = vt0crt1[0]->Time();
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

    if (fmctcrt1.isValid()){
      auto const& vctcrt1 = fmctcrt1.at(track.key());
      if (!vctcrt1.empty()){
        this_crt1x = vctcrt1[0]->EndPoint1()[0];
        this_crt1y = vctcrt1[0]->EndPoint1()[1];
        this_crt1z = vctcrt1[0]->EndPoint1()[2];
      }
    }

    if (fmpfp.isValid()){
      auto const &pfps = fmpfp.at(track.key());
      if(!pfps.empty()){
        //Find T0 for PFParticle
        if (fmt0pandora.isValid()){
          auto const &t0s = fmt0pandora.at(pfps[0].key());
          if(!t0s.empty()){
            this_t0pandora = t0s[0]->Time();
          }
        }
        if (fmt0anodepiercer.isValid()){
          auto const &t0aps = fmt0anodepiercer.at(pfps[0].key());
          if (!t0aps.empty()){
            this_t0anodep = t0aps[0]->Time();
          }
        }
      }
    }
    if (this_t0crt2 > -DBL_MAX ||
        this_t0crt1 > -DBL_MAX ||
        this_t0pandora > -DBL_MAX ||
        this_t0anodep > -DBL_MAX){

      auto const & allHits = hitsFromTrack.at(track.key());

      if (!e.isRealData()){
        // Find true particle for reconstructed track
        int TrackID = 0;
        std::map<int,double> trkide;
        for(size_t h = 0; h < allHits.size(); ++h){
          art::Ptr<recob::Hit> hit = allHits[h];
          std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
          for(size_t e = 0; e < TrackIDs.size(); ++e){
            trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
          }	    
        }
        // Work out which IDE despoited the most charge in the hit if there was more than one.
        double maxe = -1;
        double tote = 0;
        for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
          tote += ii->second;
          if ((ii->second)>maxe){
            maxe = ii->second;
            TrackID = ii->first;
          }
        }
        // Now have trackID, so get PdG code and T0 etc.
        const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackID);
        if (particle){
          this_t0truth = particle->T();
        }
      }
      
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
        else if (this_t0crt1 > -DBL_MAX) ticksOffset = this_t0crt1/500.+detProp.GetXTicksOffset(allHits[0]->WireID());
        else if (this_t0anodep > -DBL_MAX) ticksOffset = this_t0anodep/500.+detProp.GetXTicksOffset(allHits[0]->WireID());
        double xOffset = detProp.ConvertTicksToX(ticksOffset,allHits[0]->WireID());
        this_trackstartx -= xOffset;
        this_trackendx -= xOffset;
      }

      auto const & posOffsets = SCE->GetCalPosOffsets(geo::Point_t(this_trackstartx, this_trackstarty, this_trackstartz), allHits[0]->WireID().TPC);
      this_trackstartx_sce = this_trackstartx - posOffsets.X();
      this_trackstarty_sce = this_trackstarty + posOffsets.Y();
      this_trackstartz_sce = this_trackstartz + posOffsets.Z();
      this_trackendx_sce = this_trackendx - posOffsets.X();
      this_trackendy_sce = this_trackendy + posOffsets.Y();
      this_trackendz_sce = this_trackendz + posOffsets.Z();

      trackid.push_back(this_trackid);
      t0crt2.push_back(this_t0crt2);
      t0crt1.push_back(this_t0crt1);
      t0pandora.push_back(this_t0pandora);
      t0anodep.push_back(this_t0anodep);
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
      crt1x.push_back(this_crt1x);
      crt1y.push_back(this_crt1y);
      crt1z.push_back(this_crt1z);
      crt2x0.push_back(this_crt2x0);
      crt2y0.push_back(this_crt2y0);
      crt2z0.push_back(this_crt2z0);
      crt2x1.push_back(this_crt2x1);
      crt2y1.push_back(this_crt2y1);
      crt2z1.push_back(this_crt2z1);
    }
  }

  //RDTimeStamps for the event
  art::Handle < std::vector < raw::RDTimeStamp > > rdts_evt;
  std::vector < art::Ptr < raw::RDTimeStamp > > rdtsevtList;
  if (e.getByLabel("timingrawdecoder:daq", rdts_evt)) {
    art::fill_ptr_vector(rdtsevtList, rdts_evt);
  }
  
  if (!rdtsevtList.empty()){
    rdtimestamp_evt = rdtsevtList[0]->GetTimeStamp();
  }
  
  //RDTimeStamps for each raw digit (channel)
  art::Handle < std::vector < raw::RDTimeStamp > > rdts_digit;
  std::vector < art::Ptr < raw::RDTimeStamp > > rdtsdigitList;
  if (e.getByLabel("tpcrawdecoder:daq", rdts_digit)) {
    art::fill_ptr_vector(rdtsdigitList, rdts_digit);
  }
  
  for (const auto & rdts : rdtsdigitList){
    rdtimestamp_digits.push_back(rdts->GetTimeStamp());
  }

  if (!trackid.empty()) t0tree->Fill();
}

void pdsp::CheckT0::beginJob() {
  art::ServiceHandle<art::TFileService> fileServiceHandle;
  t0tree = fileServiceHandle->make<TTree>("t0tree", "t0 info");
  t0tree->Branch("run", &run, "run/I");
  t0tree->Branch("event", &event, "event/I");
  t0tree->Branch("trackid", &trackid);
  t0tree->Branch("t0crt2", &t0crt2);
  t0tree->Branch("t0crt1", &t0crt1);
  t0tree->Branch("t0pandora",&t0pandora);
  t0tree->Branch("t0anodep",&t0anodep);
  t0tree->Branch("t0truth", &t0truth);
  t0tree->Branch("trackstartx",&trackstartx);
  t0tree->Branch("trackstarty",&trackstarty);
  t0tree->Branch("trackstartz",&trackstartz);
  t0tree->Branch("trackendx",&trackendx);
  t0tree->Branch("trackendy",&trackendy);
  t0tree->Branch("trackendz",&trackendz);
  t0tree->Branch("trackstartx_sce",&trackstartx_sce);
  t0tree->Branch("trackstarty_sce",&trackstarty_sce);
  t0tree->Branch("trackstartz_sce",&trackstartz_sce);
  t0tree->Branch("trackendx_sce",&trackendx_sce);
  t0tree->Branch("trackendy_sce",&trackendy_sce);
  t0tree->Branch("trackendz_sce",&trackendz_sce);
  t0tree->Branch("crt1x",&crt1x);
  t0tree->Branch("crt1y",&crt1y);
  t0tree->Branch("crt1z",&crt1z);
  t0tree->Branch("crt2x0",&crt2x0);
  t0tree->Branch("crt2y0",&crt2y0);
  t0tree->Branch("crt2z0",&crt2z0);
  t0tree->Branch("crt2x1",&crt2x1);
  t0tree->Branch("crt2y1",&crt2y1);
  t0tree->Branch("crt2z1",&crt2z1);
  t0tree->Branch("rdtimestamp_evt", &rdtimestamp_evt, "rdtimestamp_evt/D");
  t0tree->Branch("rdtimestamp_digits", &rdtimestamp_digits);
}
DEFINE_ART_MODULE(pdsp::CheckT0)
