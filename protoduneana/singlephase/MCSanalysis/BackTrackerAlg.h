#ifndef BACKTRACKERALG_H
#define BACKTRACKERALG_H

#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"

namespace lsu
{
  class BackTrackerAlg
  {
   public:
    //Give a collection of hits, get an MCParticle
    const simb::MCParticle getMCParticle(const std::vector<art::Ptr<recob::Hit> >& hits, art::Event const & e);
    const simb::MCParticle getMCParticle(const std::vector<art::Ptr<recob::Hit> >& hits, art::Event const & e, bool& isNeg);


    //Give a track, get an MCParticle
    // Sometimes a delta ray will create hits and an MCParticle is not created.  
    // If a track's highest charge purity track ID is one of these delta rays, 
    // a negative track id is returned, corresponding to the MCParticle that created
    // this delta ray.  In this event, that MCParticle will be returned, and isNeg will be set.
    const simb::MCParticle getMCParticle(art::Ptr<recob::Track> Track,const art::Event & e, bool& isNeg);
    const simb::MCParticle getMCParticle(art::Ptr<recob::Track> Track,const art::Event & e);

    // Tag is for "pandoraTrack" or "pmtrack" etc.
    const simb::MCParticle getMCParticle(art::Ptr<recob::Track> Track,const art::Event & e, std::string tag);
    const simb::MCParticle getMCParticle(art::Ptr<recob::Track> Track,const art::Event & e, std::string tag, bool& isNeg);

   private:
    std::string _defaultTag = "pandoraTrack";

    // Helper Functions
    void getHCPTrackID(std::set<int> setOfTrackIDs,int &highestChargePurityTrackID,double &highestChargePurity,art::ServiceHandle<cheat::BackTrackerService> bt,const std::vector<art::Ptr<recob::Hit> > hits, const detinfo::DetectorClocksData &clockData);
    void finishGetting(int highestChargePurityTrackID,std::vector<art::Ptr<simb::MCParticle> > mcparticleList, art::Ptr<simb::MCParticle> &mcparticle);
    
  };
}
#endif
