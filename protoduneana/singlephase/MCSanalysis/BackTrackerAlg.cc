#include "BackTrackerAlg.h"

const simb::MCParticle lsu::BackTrackerAlg::getMCParticle(const std::vector<art::Ptr<recob::Hit> >& hits, art::Event const & e, bool& isNeg) {

  art::ServiceHandle<cheat::BackTrackerService> bt;
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  art::Ptr<simb::MCParticle> mcparticle;
  
  //Make MCParticle List
  art::Handle<std::vector<simb::MCParticle> > mcparticleListHandle;
  std::vector<art::Ptr<simb::MCParticle> > mcparticleList;
  if(e.getByLabel("largeant",mcparticleListHandle))
    art::fill_ptr_vector(mcparticleList,mcparticleListHandle);

  // Get Set of TrackIds and figure out which has the highest charge purity
  std::set<int> setOfTrackIDs = bt->GetSetOfTrackIds(clockData, hits); //all trackIDs for hits in this track
  int highestChargePurityTrackID = -999999;
  double highestChargePurity = -999999;
  getHCPTrackID(setOfTrackIDs,highestChargePurityTrackID,highestChargePurity,bt,hits, clockData);
  // At this point, highestChargePurityTrackID is set, so is highestChargePurity

  if(highestChargePurityTrackID >= 0)
    isNeg = false;
  else
    isNeg = true;

  mcparticle = mcparticleList.at(0);
  if(highestChargePurity >= 0.)
    finishGetting(highestChargePurityTrackID,mcparticleList,mcparticle);
  else
    std::cout << "Error: Couldn't find True Track ID for cluster " << std::endl << std::endl;

  return *mcparticle;

} // getMCParticle(hits, e, isNeg)

const simb::MCParticle lsu::BackTrackerAlg::getMCParticle(const std::vector<art::Ptr<recob::Hit> > & hits, art::Event const & e) {

  bool isNeg = false;
  const simb::MCParticle particle = getMCParticle(hits, e, isNeg);
  return particle;

}

// TODO: Make this the default one, not the other one.
const simb::MCParticle lsu::BackTrackerAlg::getMCParticle(art::Ptr<recob::Track> Track, const art::Event & e, std::string tag, bool& isNeg)
{
  std::vector<art::Ptr<recob::Track> > trackList;
  trackList.push_back(Track);
  const auto trackHitList = art::FindManyP<recob::Hit>(trackList,e,tag);
  const auto trackHits = trackHitList.at(0);
  const simb::MCParticle particle = getMCParticle(trackHits,e,isNeg);
  return particle;
}

const simb::MCParticle lsu::BackTrackerAlg::getMCParticle(art::Ptr<recob::Track> Track, const art::Event & e, bool& isNeg) {

  const simb::MCParticle particle = getMCParticle(Track, e, _defaultTag, isNeg);
  return particle;

} // getMCParticle(Track, e, isNeg)

const simb::MCParticle lsu::BackTrackerAlg::getMCParticle(art::Ptr<recob::Track> Track, const art::Event & e, std::string tag)
{
  bool isNeg = false;
  const simb::MCParticle particle = getMCParticle(Track, e, tag, isNeg);
  if(isNeg)
    std::cout << "WARNING: Track " << Track->ID() 
	      << " in event " << e.event() 
	      << " returns a negative track id that is not accounted for, please consider using an alternative getMCParticle function." 
	      << std::endl;
  return particle;
} // getMCParticle(Track, e, tag)

const simb::MCParticle lsu::BackTrackerAlg::getMCParticle(art::Ptr<recob::Track> Track, const art::Event & e) {

  const simb::MCParticle particle = getMCParticle(Track, e, _defaultTag);
  return particle;

} // getMCParticle(Track, e)

// Private Helper Functions
void lsu::BackTrackerAlg::getHCPTrackID(std::set<int> setOfTrackIDs,int &highestChargePurityTrackID,double &highestChargePurity,art::ServiceHandle<cheat::BackTrackerService> bt, const std::vector<art::Ptr<recob::Hit> > hits, const detinfo::DetectorClocksData &clockData)
{
  std::set<int> thisTrackIDSet;
  for(const auto & trackID : setOfTrackIDs)
    {
      thisTrackIDSet.clear();
      thisTrackIDSet.insert(trackID);
      // what fraction of the track charge came from this particle
      float chargePurity = bt->HitChargeCollectionPurity(clockData, thisTrackIDSet, hits);
      //we want the particle with the largest fraction of the charge
      if(chargePurity > highestChargePurity)
	{
	  highestChargePurity = chargePurity;
	  highestChargePurityTrackID = trackID;
	}
    } // trackID : setOfTrackIDs
}

void lsu::BackTrackerAlg::finishGetting(int highestChargePurityTrackID,std::vector<art::Ptr<simb::MCParticle> > mcparticleList,art::Ptr<simb::MCParticle> &mcparticle)
{
  // Find the MCParticle with this trackID
  for(auto mcpart: mcparticleList)
    {
      if(mcpart->TrackId() == highestChargePurityTrackID)
	{
	  mcparticle = mcpart;
	  break;
	}
    }
  if(mcparticle.isNull())
    {
      std::cout << "Warning: Couldn't find MCParticle for TrackID: " << highestChargePurityTrackID << " trying to find: " << abs(highestChargePurityTrackID) << std::endl;
      for(auto mcpart: mcparticleList)
	{
	  if(mcpart->TrackId() == abs(highestChargePurityTrackID))
	    {
	      mcparticle = mcpart;
	      break;
	    }
	}
      if(mcparticle.isNull())
	{
	  std::string message = "Couldn't find MCParticle for TrackID: ";
	  message.append(std::to_string(highestChargePurityTrackID));
	  throw cet::exception("MCparticleNotFound",message);
	}
    } // if can't find MCParticle for highest Track ID
}
