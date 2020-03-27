#ifndef CONE_H
#define CONE_H

#include <iostream>
#include <vector>
#include <string>

// art and larsoft includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larsim/MCCheater/BackTrackerService.h"

namespace pizero {

// Class to store information of hits that fall within it.
class Cone {
 private:
  const art::Event* m_evt;
  TVector3 m_start;
  TVector3 m_direction;
  double m_length;
  double m_ang = 15.0/360 * 2*3.1415926; // Opening angle of cone in radians.
  double m_width = 1; // Has to be overwritten.
  double m_startoffset = 10;

  std::vector<const recob::Hit*> m_hits;
  std::vector<const recob::SpacePoint*> m_spacepoints;

 public:
  Cone(const art::Event& evt, const TVector3& newStart, const TVector3& newDir,
       const double newLen):
       m_evt(&evt), m_start(newStart), m_direction(newDir.Unit()), m_length(newLen) {
    // Move the start back by a constant to catch all hits.
    m_start = m_start - m_direction.Unit()*m_startoffset;
    m_length += m_startoffset;
    // Determine the width of the cone at the end.
    m_width = m_length * tan(m_ang);

    // Get all spacepoints in the event.
    std::string spsLabel = "pandora";
    art::Handle<std::vector<recob::SpacePoint>> spHandle;
    if (!m_evt->getByLabel(spsLabel, spHandle)) { return; }
    // Get associations between the spacepoints and hits.
    const art::FindManyP<recob::Hit> sphits(spHandle, *m_evt, spsLabel);

    // Loop through spacepoints to test whether they're in the cone.
    for(const recob::SpacePoint& sp : *spHandle) {
      const TVector3 pos(sp.XYZ());
      // 0 <= projection <= cone length
      const double proj = (pos-m_start).Dot(m_direction);
      if(proj < 0 || proj > m_length) continue;
      // Cone radius at projection.
      const double cradius = (proj/m_length)*m_width/2;
      // rejection < cone radius
      const double rej = ((pos-m_start) - proj*m_direction).Mag();
      if(rej > cradius) continue;

      // This spacepoint is inside the cone. Put it among the results.
      for(const art::Ptr<recob::Hit>& hit : sphits.at(sp.ID())) {
        m_hits.push_back(&*hit);
      }
      m_spacepoints.push_back(&sp);
    }
  }

  // Member getters
  TVector3 start() const { return m_start; }
  TVector3 end() const { return m_start + m_direction*m_length; }
  TVector3 direction() const { return m_direction; }
  double length() const { return m_length; }
  double width() const { return m_width; }
  std::vector<const recob::Hit*> hits() const { return m_hits; }
  std::vector<const recob::SpacePoint*> spacepoints() const { return m_spacepoints; }

  double energy(calo::CalorimetryAlg caloAlg) const {
    protoana::ProtoDUNEShowerUtils shUtil;
    return shUtil.EstimateEnergyFromHitCharge(m_hits, caloAlg)[2];
  }
  double completeness(const simb::MCParticle& mcpart,
                      const std::string hitLabel = "hitpdune") const;
  double purity(const simb::MCParticle& mcpart) const;
};

double Cone::completeness(const simb::MCParticle& mcpart,
                          const std::string hitLabel) const {
  // Get all hits in the event.
  art::Handle<std::vector<recob::Hit>> hitHandle;
  if(!m_evt->getByLabel(hitLabel, hitHandle)) {
    std::cout << "pizero::Cone::completeness: could not find hits in event.\n";
    return 0;
  }

  // Check for each hit whether it came from the MCParticle and if so,
  // whether it fell inside the cone.
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  double sharedHits = 0;
  double MChits = 0;
  for(const recob::Hit& hit : *hitHandle) {
    for(const int trackId : bt_serv->HitToTrackIds(hit)) {
      if(std::abs(trackId) == std::abs(mcpart.TrackId())) {
        // Hit is in MCParticle.
        ++MChits;
        // Check if hit in cone.
        for(const recob::Hit* chit : m_hits) {
          // Compare the time and place of hit for lack of ID
          if(chit->PeakTime() - hit.PeakTime() < 1e-5 &&
             chit->Channel() == hit.Channel()) {
            ++sharedHits;
            // std::cout << "Found shared hit!\n";
            break;
          }
        }
        break;
      } // if hit is from MCP
    } // for trackID in hit
  } // for hit in event

  return MChits==0? 0: sharedHits/MChits;
} // Cone::completeness

double Cone::purity(const simb::MCParticle& mcpart) const {
  unsigned sharedHits = 0;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  for(const recob::Hit* hit : m_hits) {
    for(const int trackId : bt_serv->HitToTrackIds(*hit)) {
      if(std::abs(trackId) == std::abs(mcpart.TrackId())) {
        ++sharedHits;
        break;
      }
    }
  }

  return (double)sharedHits / m_hits.size();
}


} // namespace pizero

#endif // CONE_H
