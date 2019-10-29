#ifndef SHOWER_PROCESS_H
#define SHOWER_PROCESS_H

#include <iostream>
#include <vector>
#include <string>

// Larsoft includes
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "protoduneana/Analysis/ProtoDUNETruthUtils.h"

namespace pizero {

// Class to store all relevant objects and variables of a shower process.
class ShowerProcess {
 private:
  // Shower process objects.
  std::vector<const simb::MCParticle*> m_mcparts;
  std::vector<const recob::Shower*> m_showers;
  std::vector<const recob::Track*> m_tracks;

  // Event information.
  const art::Event* m_evt = 0x0;
  std::string m_showerLabel;
  std::string m_trackLabel;
  std::string m_simulationLabel;

  // MCParticle and shower finders through the truth utilities.
  void find_mcparticles();
  void find_showers();
  void find_tracks();

  // Truth utilities for convenience.
  protoana::ProtoDUNETruthUtils truthUtils;

 public:
  // Constructors from MCParticle, shower and both.
  ShowerProcess(const simb::MCParticle& mcpart, const art::Event& evt,
                std::string showerLabel = "pandoraShower",
                std::string trackLabel = "pandoraTrack",
                std::string simulationLabel = "largeant"):
                m_evt(&evt),
                m_showerLabel(showerLabel),
                m_trackLabel(trackLabel),
                m_simulationLabel(simulationLabel) {
    m_mcparts.push_back(&mcpart);
    find_showers();
    find_tracks();

    // std::cout << "ShowerProcess found " << m_mcparts.size() << " MCParticles, "
    //           << m_showers.size() << " showers and "
    //           << m_tracks.size() << " tracks.\n";
  }
  ShowerProcess(const recob::Shower& shower, const art::Event& evt,
                std::string showerLabel = "pandoraShower",
                std::string trackLabel = "pandoraTrack",
                std::string simulationLabel = "largeant"):
                m_evt(&evt),
                m_showerLabel(showerLabel),
                m_trackLabel(trackLabel),
                m_simulationLabel(simulationLabel) {
    m_showers.push_back(&shower);
    find_mcparticles();
    find_tracks();

    // std::cout << "ShowerProcess found " << m_mcparts.size() << " MCParticles, "
    //           << m_showers.size() << " showers and "
    //           << m_tracks.size() << " tracks.\n";
  }
  ShowerProcess(const simb::MCParticle& mcpart, const recob::Shower& shower,
                const art::Event& evt,
                std::string showerLabel = "pandoraShower",
                std::string trackLabel = "pandoraTrack",
                std::string simulationLabel = "largeant"):
                m_evt(&evt),
                m_showerLabel(showerLabel),
                m_trackLabel(trackLabel),
                m_simulationLabel(simulationLabel) {
    m_mcparts.push_back(&mcpart);
    m_showers.push_back(&shower);
    find_tracks();

    // std::cout << "ShowerProcess found " << m_mcparts.size() << " MCParticles, "
    //           << m_showers.size() << " showers and "
    //           << m_tracks.size() << " tracks.\n";
  }

  // Related object getters. The elements are guaranteed to be descending in energy.
  const simb::MCParticle* mcparticle() const {
    return m_mcparts.size()!=0? m_mcparts[0]: 0x0;
  }
  const recob::Shower* shower() const {
    return m_showers.size()!=0? m_showers[0]: 0x0;
  }
  const recob::Track* track() const {
    return m_tracks.size()!=0? m_tracks[0]: 0x0;
  }
  // Object vector getters.
  std::vector<const simb::MCParticle*> mcparticles() const { return m_mcparts; }
  std::vector<const recob::Shower*> showers() const { return m_showers; }
  std::vector<const recob::Track*> tracks() const { return m_tracks; }
}; // class ShowerProcess

// Find MCParticle based on the biggest shower via truth utilities.
void ShowerProcess::find_mcparticles() {
  if(m_evt->isRealData() || m_showers.size() == 0) return;
  m_mcparts.push_back(truthUtils.GetMCParticleFromReco(*m_showers[0], *m_evt, m_showerLabel));
}

// Find showers based on the given MCParticle.
void ShowerProcess::find_showers() {
  if(m_mcparts.size() == 0) return;

  // Find all showers this MCParticle contributed to.
  std::vector<std::pair<const recob::Shower*, double>> showers =
    truthUtils.GetRecoShowerListFromMCParticle(*m_mcparts[0], *m_evt, m_showerLabel);
  // Only showers to which the MCParticle was the primary contributor are saved.
  for(const std::pair<const recob::Shower*, double>& sh : showers) {
    const simb::MCParticle* sh_part =
      truthUtils.GetMCParticleFromReco(*sh.first, *m_evt, m_showerLabel);
    if(std::abs(sh_part->TrackId()) == std::abs(m_mcparts[0]->TrackId())) {
      m_showers.push_back(sh.first);
    }
  }
}

// Find tracks based on the given MCParticle or shower if no MC is present.
void ShowerProcess::find_tracks() {
  // For now return if no MC present.
  if(m_mcparts.size() == 0) return;

  // Find all tracks this MCParticle contributed to.
  std::vector<std::pair<const recob::Track*, double>> tracks =
    truthUtils.GetRecoTrackListFromMCParticle(*m_mcparts[0], *m_evt, m_trackLabel);
  // Only tracks to which the MCParticle was the primary contributor are saved.
  for(const std::pair<const recob::Track*, double>& tr : tracks) {
    const simb::MCParticle* tr_part =
      truthUtils.GetMCParticleFromReco(*tr.first, *m_evt, m_trackLabel);
    if(std::abs(tr_part->TrackId()) == std::abs(m_mcparts[0]->TrackId())) {
      m_tracks.push_back(tr.first);
    }
  }
}

} // namespace pizero

#endif // SHOWERPROCESS_H
