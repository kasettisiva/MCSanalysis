#ifndef SHOWER_PROCESS_H
#define SHOWER_PROCESS_H

#include <iostream>
#include <vector>
#include <string>

// Larsoft includes
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "dune/Protodune/singlephase/DataUtils/ProtoDUNETruthUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEShowerUtils.h"
#include "Cone.h"

#include "TFile.h"
#include "TH3F.h"

namespace pizero {

// Class to store all relevant objects and variables of a shower process.
class ShowerProcess {
 private:
  // Shower process objects.
  std::vector<const simb::MCParticle*> m_mcparts;
  std::vector<const recob::Shower*> m_showers;
  std::vector<const recob::Track*> m_tracks;
  // double cone_length = 100;
  // double cone_width = 20;
  Cone* m_cone = 0x0;

  // Event information.
  const art::Event* m_evt = 0x0;
  std::string m_showerLabel;
  std::string m_trackLabel;
  std::string m_simulationLabel;

  // MCParticle and reconstructed object finders.
  void find_mcparticles();
  void find_showers();
  void find_tracks();
  void fill_cone();

  // Utilities for convenience.
  protoana::ProtoDUNETruthUtils truthUtils;
  protoana::ProtoDUNEShowerUtils shUtils;

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
    fill_cone();
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
    fill_cone();
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
    fill_cone();
  }
  ~ShowerProcess() { if(m_cone != 0x0) delete m_cone; } // TODO: make unique_ptr.

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
  const Cone* cone() const {
    return m_cone;
  }
  // Object vector getters.
  std::vector<const simb::MCParticle*> mcparticles() const { return m_mcparts; }
  std::vector<const recob::Shower*> showers() const { return m_showers; }
  std::vector<const recob::Track*> tracks() const { return m_tracks; }
  // Calorimetry.
  double energy(const std::string& caloLabel, const double calib,
    const double norm, const std::string& scefile, const std::string& yzfile,
    const std::string& xfile) const;
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

  using weightedTrackPair = std::pair<const recob::Track*, double>;
  std::vector<weightedTrackPair> outVec;

  // Get all reconstructed tracks
  auto allRecoTracks = m_evt->getValidHandle<std::vector<recob::Track> >(m_trackLabel);

  // We need the association between the tracks and the hits
  const art::FindManyP<recob::Hit> findTrackHits(allRecoTracks, *m_evt, m_trackLabel);

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  std::unordered_map<int, double> recoTrack;

  // Record the energy contribution to the MCParticle of every relevant reco track
  double part_E = 0;
  for (recob::Track const & track : *allRecoTracks) {
    // Record the contribution of all MCParticles to this track.
    std::unordered_map<int, double> mcpmap;
    // Loop through the track's hits.
    for (art::Ptr<recob::Hit> const & hptr : findTrackHits.at(track.ID())) {
      // Loop over hit IDEs to find our MCParticle's part in it
      for (const sim::IDE* ide : bt_serv->HitToSimIDEs_Ps(hptr)) {
        mcpmap[abs(ide->trackID)] += ide->energy;
        if (abs(ide->trackID) == mcparticle()->TrackId()) {
          recoTrack[track.ID()] += ide->energy;
        }
        part_E += ide->energy;
      } // sim::IDE*
    } // art::Ptr<recob::Hit>

    // Remove track if the shower's MCParticle was not the main contributor.
    if(recoTrack.size() > 0) {
      int maxMCPID = mcpmap[0];
      if(mcpmap.size() > 1) {
        maxMCPID = std::max_element(mcpmap.begin(), mcpmap.end(),
          [](const std::pair<int, double>& a, const std::pair<int, double>& b){
            return a.second < b.second;
          })->first;
      }
      if(maxMCPID != mcparticle()->TrackId()) {
        recoTrack.erase(track.ID());
      }
    }
  } // const recob::Track&
  // for(auto const& p : recoTrack) {
  //   std::cout << "Track with ID " << p.first << " has energy " << p.second << '\n';
  // }

  // Fill and sort the output vector
  for (std::pair<int, double> const& p : recoTrack)
  {
    auto const trackIt = std::find_if(allRecoTracks->begin(), allRecoTracks->end(),
                               [&](recob::Track tr){ return tr.ID() == p.first; });
    outVec.push_back(std::make_pair(&*trackIt, p.second));
  }
  std::sort(outVec.begin(), outVec.end(),
    [](weightedTrackPair a, weightedTrackPair b){ return a.second > b.second;});

  // Normalise the vector weights
  if (part_E < 1e-5) { part_E = 1; } // Protect against zero division
  for (weightedTrackPair& p : outVec)
  {
    p.second /= part_E;
    m_tracks.push_back(p.first);
  }

  // // Find all tracks this MCParticle contributed to.
  // std::vector<std::pair<const recob::Track*, double>> tracks =
  //   truthUtils.GetRecoTrackListFromMCParticle(*m_mcparts[0], *m_evt, m_trackLabel);
  // // Only tracks to which the MCParticle was the primary contributor are saved.
  // for(const std::pair<const recob::Track*, double>& tr : tracks) {
  //   const simb::MCParticle* tr_part =
  //     truthUtils.GetMCParticleFromReco(*tr.first, *m_evt, m_trackLabel);
  //   if(std::abs(tr_part->TrackId()) == std::abs(m_mcparts[0]->TrackId())) {
  //     m_tracks.push_back(tr.first);
  //   }
  // }
}

// Fill the cone object based on the showers and tracks in the class.
void ShowerProcess::fill_cone() {
  // TODO: Find the start and direction of first reconstructed object.
  TVector3 cstart, cdir;
  if(m_showers.size() != 0) {
    cstart = m_showers[0]->ShowerStart();
    cdir = m_showers[0]->Direction();
    // Make cone object
    m_cone = new Cone(*m_evt, cstart, cdir, m_showers[0]->Length());
  }
  // else if(m_tracks.size() != 0) {
  //   const recob::tracking::Point_t tstart = m_tracks[0]->Start();
  //   const recob::tracking::Vector_t tdir = m_tracks[0]->StartDirection();
  //   cstart = TVector3(tstart.X(), tstart.Y(), tstart.Z());
  //   cdir = TVector3(tdir.X(), tdir.Y(), tdir.Z());
  // }

}

double ShowerProcess::energy(const std::string& caloLabel, const double calib,
  const double norm, const std::string& scefile, const std::string& yzfile,
  const std::string& xfile) const {
  double totE = 0;

  // Check whether there is a shower associated with this object.
  if(shower() == 0x0) return totE;

  // Constants from DUNE docDB 15974 by A Paudel.
  const double Wion = 23.6e-6; // ArgoNeuT-determined parameter at E0 kV/cm
  const double recob = 0.6417; // Recombination factor from Aaron's talk.

  // Get dQ/dx YZ and X correction factors.
  TFile* yzf = new TFile(yzfile.c_str());
  TH2F* yzcorrposhist =(TH2F*)yzf->Get("correction_dqdx_ZvsY_positiveX_hist_2");
  TH2F* yzcorrneghist =(TH2F*)yzf->Get("correction_dqdx_ZvsY_negativeX_hist_2");
  TFile* xf = new TFile(xfile.c_str());
  TH1F* xcorrhist = (TH1F*)xf->Get("dqdx_X_correction_hist_2");

  // Get the calorimetry objects from the event using the shower utilities.
  std::vector<anab::Calorimetry> calovec =
    shUtils.GetRecoShowerCalorimetry(*shower(), *m_evt, m_showerLabel, caloLabel);

  // Loop over all plane calorimetry objects, but only use the collection plane.
  for(const anab::Calorimetry& calo : calovec) {
    if(calo.PlaneID().Plane != 2) continue;

    // Loop over all hits in the collection plane calorimetry object.
    for(unsigned i = 0; i < calo.dQdx().size(); ++i) {

      double dQdx = calo.dQdx()[i];
      const double pitch = calo.TrkPitchVec()[i];
      const geo::Point_t& pos = calo.XYZ()[i];
      const double x = pos.X();
      const double y = pos.Y();
      const double z = pos.Z();
      // Can't apply location-dependent corrections if there's no location available.
      if(x==0 && y==0 && z==0) {
        totE += dQdx*pitch*norm/calib*Wion/recob * 1e-3;
        continue;
      }
      // Spatial corrections.
      double yzcorr;
      if(x >= 0) {
        yzcorr = yzcorrposhist->GetBinContent(yzcorrposhist->FindBin(z,y));
      } else {
        yzcorr = yzcorrneghist->GetBinContent(yzcorrneghist->FindBin(z,y));
      }
      // Correct dQ/dx.
      // std::cout << "Bare energy (GeV): " << dQdx*pitch*norm/calib*Wion/recob * 1e-3 << '\n';
      const double xcorr = xcorrhist->GetBinContent(xcorrhist->FindBin(x));
      // std::cout << "Hit energy (GeV): " << xcorr*yzcorr*dQdx*pitch*norm/calib*Wion/recob * 1e-3 << "\n\n";
      totE += xcorr*yzcorr*dQdx*pitch*norm/calib*Wion/recob * 1e-3; // Convert from MeV to GeV.
    } // for hits in calo
  } // for calo objects

  // std::cout << "Total energy: " << totE << '\n';
  return totE;
}

} // namespace pizero

#endif // SHOWERPROCESS_H
