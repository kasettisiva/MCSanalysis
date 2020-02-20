////////////////////////////////////////////////////////////////////////
// Class:       MDMAna
// Plugin Type: analyzer (art v3_04_00)
// File:        MDMAna_module.cc
//
// Generated at Wed Feb 19 22:04:32 2020 by Tingjun Yang using cetskelgen
// from cetlib version v3_09_00.
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
#include "art_root_io/TFileService.h"

#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "TTree.h"

#include <vector>

namespace pdsp {
  class MDMAna;
}


class pdsp::MDMAna : public art::EDAnalyzer {
public:
  explicit MDMAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MDMAna(MDMAna const&) = delete;
  MDMAna(MDMAna&&) = delete;
  MDMAna& operator=(MDMAna const&) = delete;
  MDMAna& operator=(MDMAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.
  TTree *mdmtree;
  int run;
  int subrun;
  int event;
  std::vector<int> trackid;
  std::vector<double> tracklen;
  std::vector<int> trackorg;
  std::vector<int> trackpdg;

};


pdsp::MDMAna::MDMAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdsp::MDMAna::analyze(art::Event const& e)
{
  run = e.run();
  subrun = e.subRun();
  event = e.id().event();
  trackid.clear();
  tracklen.clear();
  trackorg.clear();
  trackpdg.clear();

  //Services
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

   // Reconstruciton information
  art::Handle < std::vector < recob::Track > > trackListHandle;
  std::vector < art::Ptr < recob::Track > > trackList;
  if (e.getByLabel("pandoraTrack", trackListHandle)) {
    art::fill_ptr_vector(trackList, trackListHandle);
  }
  else return;

  //Get hits associated with track
  art::FindManyP < recob::Hit > hitsFromTrack(trackListHandle, e, "pandoraTrack");

  for (auto const& track : trackList){
    int this_trackid = track.key();
    double this_tracklen = track->Length();
    int this_trackorg = -1;
    int this_trackpdg = 0;

    auto const & allHits = hitsFromTrack.at(track.key());

    if (!e.isRealData()){
      // Find true particle for reconstructed track
      int TrackID = 0;
      std::map<int,double> trkide;
      for(size_t h = 0; h < allHits.size(); ++h){
        art::Ptr<recob::Hit> hit = allHits[h];
        std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(hit);
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
      if (TrackID){
        const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackID);
        if (particle){
          this_trackpdg = particle->PdgCode();
          auto & truth = pi_serv->ParticleToMCTruth_P(particle);
          this_trackorg = truth->Origin();
        }
      }
    }
    trackid.push_back(this_trackid);
    tracklen.push_back(this_tracklen);
    trackorg.push_back(this_trackorg);
    trackpdg.push_back(this_trackpdg);
  }
  if (!trackid.empty()) mdmtree->Fill();
}

void pdsp::MDMAna::beginJob()
{
  art::ServiceHandle<art::TFileService> fileServiceHandle;
  mdmtree = fileServiceHandle->make<TTree>("mdmtree", "MDM info");
  mdmtree->Branch("run", &run, "run/I");
  mdmtree->Branch("subrun", &subrun, "subrun/I");
  mdmtree->Branch("event", &event, "event/I");
  mdmtree->Branch("trackid", &trackid);
  mdmtree->Branch("tracklen", &tracklen);
  mdmtree->Branch("trackorg", &trackorg);
  mdmtree->Branch("trackpdg", &trackpdg);

}

DEFINE_ART_MODULE(pdsp::MDMAna)
