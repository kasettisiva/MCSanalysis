#ifndef PIZERO_PROCESS_H
#define PIZERO_PROCESS_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "TTree.h"

// Larsoft includes
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEShowerUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNETruthUtils.h"
#include "protoduneana/singlephase/EMTaskForce/ShowerProcess.h"

namespace pizero {

// Function to find the closest distance between a line and a point.
double ClosestDistanceToPoint(const TVector3& p,
                  const TVector3& lineStart, const TVector3& lineDir) {
  // Difference vector between line start and point.
  const TVector3 diff = p - lineStart;
  // Get line projection onto difference vector.
  const TVector3 proj = lineDir * diff.Dot(lineDir) * (1.0 / lineDir.Mag2());
  // The difference between the difference and projection vectors is
  // called the rejection vector and returns the closest distance.
  return (diff - proj).Mag();
}

// Function to find the point in the middle of the shortest distance between two lines.
TVector3 ClosestPoint(const TVector3& a, const TVector3& adir,
                      const TVector3& b, const TVector3& bdir) {
  TVector3 result(-99999., -99999., -99999.);
  // Deal with both vectors being parallel.
  if (adir != bdir) {
    const TVector3 n1 = adir.Cross(bdir.Cross(adir));
    const TVector3 n2 = bdir.Cross(adir.Cross(bdir));
    const TVector3 c1 = a + adir * ((b - a).Dot(n2) / adir.Dot(n2));
    const TVector3 c2 = b + bdir * ((a - b).Dot(n1) / bdir.Dot(n1));

    result = c1 + 0.5 * (c2 - c1);
  }

  return result;
}
TVector3 ClosestPoint(const recob::Shower* showerA,
                      const recob::Shower* showerB) {
  const TVector3 candidate =
    ClosestPoint(showerA->ShowerStart(), showerA->Direction(),
                 showerB->ShowerStart(), showerB->Direction());

  return candidate;
}

// Function to find the shortest distance between two lines.
double ClosestDistance(const TVector3& a, const TVector3& adir,
                         const TVector3& b, const TVector3& bdir) {
  // Find unit vector perpendicular to both a and b.
  const TVector3 unit = adir.Cross(bdir) * (1.0 / adir.Cross(bdir).Mag());
  // Project difference vector on perpendicular vector.
  return abs((a - b).Dot(unit));
}
double ClosestDistance(const recob::Shower* showerA,
                         const recob::Shower* showerB) {
  // Check whether the intersection is behind both showers.
  const TVector3 candidate =
    ClosestPoint(showerA->ShowerStart(), showerA->Direction(),
                 showerB->ShowerStart(), showerB->Direction());
  const TVector3 startPlusDirA = showerA->ShowerStart() + showerA->Direction();
  const TVector3 startPlusDirB = showerB->ShowerStart() + showerB->Direction();
  const TVector3 startMinusDirA = showerA->ShowerStart() - showerA->Direction();
  const TVector3 startMinusDirB = showerB->ShowerStart() - showerB->Direction();
  if((startPlusDirA-candidate).Mag() < (startMinusDirA-candidate).Mag() ||
     (startPlusDirB-candidate).Mag() < (startMinusDirB-candidate).Mag()) {
    return -1;
  }

  return ClosestDistance(showerA->ShowerStart(), showerA->Direction(),
                                 showerB->ShowerStart(), showerB->Direction());
}

// Class to hold all relevant objects and information of a pi0 decay and
// reconstruction process.
class PiZeroProcess {
 private:
  // Objects associated with the process.
  std::unique_ptr<ShowerProcess> _shProcess1;
  std::unique_ptr<ShowerProcess> _shProcess2;
  // const recob::Shower* _shower1 = 0x0;
  // const recob::Shower* _shower2 = 0x0;

  const simb::MCParticle* _pi0 = 0x0;
  // const simb::MCParticle* _photon1 = 0x0;
  // const simb::MCParticle* _photon2 = 0x0;

  const art::Event& _evt;
  const std::string _showerLabel;

  // Utility objects.
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  protoana::ProtoDUNETruthUtils truthUtils;

 public:
  PiZeroProcess(const simb::MCParticle& mcp, const art::Event& evt,
                std::string showerLabel);
  PiZeroProcess(const recob::Shower& shower, const art::Event& evt,
                std::string showerLabel);
  PiZeroProcess(const recob::Shower& showerA, const recob::Shower& showerB,
                const art::Event& evt, std::string showerLabel);
  ~PiZeroProcess(){};

  bool _haveMCInfo;
  bool _haveRecoInfo;

  // Reset function to start with a clean slate.
  void reset();

  // Const access to private members.
  const simb::MCParticle* pi0() const { return _pi0; }
  const ShowerProcess* showerProcess1() const {
    return &*_shProcess1;
  }
  const ShowerProcess* showerProcess2() const {
    return &*_shProcess2;
  }
  const simb::MCParticle* photon1() const {
    return _shProcess1? _shProcess1->mcparticle(): 0x0;
  }
  const simb::MCParticle* photon2() const {
    return _shProcess2? _shProcess2->mcparticle(): 0x0;
  }
  const recob::Shower* shower1() const {
    return _shProcess1? _shProcess1->shower(): 0x0;
  }
  const recob::Shower* shower2() const {
    return _shProcess2? _shProcess2->shower(): 0x0;
  }
  const recob::Track* track1() const {
    return _shProcess1? _shProcess1->track(): 0x0;
  }
  const recob::Track* track2() const {
    return _shProcess2? _shProcess2->track(): 0x0;
  }
  const Cone* cone1() const {
    return _shProcess1? _shProcess1->cone(): 0x0;
  }
  const Cone* cone2() const {
    return _shProcess2? _shProcess2->cone(): 0x0;
  }
  const art::Event* evt() const { return &_evt; }
  std::string showerLabel() const { return _showerLabel; }

  // Functions to check whether data have been set (disallow partial reco).
  bool allMCSet() const;
  bool allRecoSet() const;

  // Function to set the pi0 and photons.
  void setPi0(const simb::MCParticle& newPi0);

  // Function to create the branches for all variables in a TTree.
  void createBranches(TTree* tree);

  // Function to fill MC info from the MCParticles.
  void fillMCInfo();

  // Function to automatically fill shower info from the found showers.
  void fillRecoInfo(calo::CalorimetryAlg caloAlg);

}; // class PiZeroProcess


// Constructor to set internal data from an MCParticle.
PiZeroProcess::PiZeroProcess(const simb::MCParticle& mcp, const art::Event& evt,
                std::string showerLabel): _evt(evt), _showerLabel(showerLabel) {
  if(mcp.PdgCode() == 111) {
    // Got a pi0.
    setPi0(mcp);
  } else if(mcp.PdgCode() == 22) {
    // Got a photon. Find the mother and run through the above process again.
    if(mcp.Mother() == 0) return;
    const simb::MCParticle* mcpmom = pi_serv->TrackIdToParticle_P(mcp.Mother());
    if(mcpmom->PdgCode() == 111) setPi0(*mcpmom);
  }

  // Quit if the system did not manage to find appropriate MCParticles.
  if(!allMCSet()) return;
}

// Constructor to set internal data from a shower.
PiZeroProcess::PiZeroProcess(const recob::Shower& shower, const art::Event& evt,
              std::string showerLabel): _evt(evt), _showerLabel(showerLabel) {
  // Get the event's showers.
  art::Handle<std::vector<recob::Shower>> showerHandle;
  if (!_evt.getByLabel(_showerLabel, showerHandle)) return;
  if(showerHandle->size() < 2) return;

  // For each shower, find the one which is closest to intersecting.
  double min_dist = 100000;
  const recob::Shower* shMatch = 0x0;
  for(const recob::Shower& potMatch : *showerHandle) {
    if(&potMatch == &shower) continue;

    const double this_dist = pizero::ClosestDistance(&shower, &potMatch);
    if(this_dist < min_dist) {
      min_dist = this_dist;
      shMatch = &potMatch;
    }
  }
  TVector3 reco_pi0_location = pizero::ClosestPoint(&shower, shMatch);

  // Find closest MC pi0 match to the reconstructed pi0 vertex.
  double reco_mismatch = 100000;
  const simb::MCParticle* min_pi0 = 0;
  for(auto const& p : pi_serv->ParticleList()) {
    const simb::MCParticle* true_pi0 = p.second;
    // std::cout << "Particle PDG: " << true_pi0->PdgCode() << '\n';
    if(true_pi0->PdgCode() != 111) continue;

    const double this_dist = (true_pi0->Position().Vect() - reco_pi0_location).Mag();
    if(this_dist < reco_mismatch) {
      reco_mismatch = this_dist;
      min_pi0 = true_pi0;
    }
  }
  if(min_pi0 == 0) return;

  // Fill internal data fields with the found pointers.
  setPi0(*min_pi0);
  // Check whether we have all necessary pointers.
  if(!allMCSet()) {
    _pi0 = 0x0;
    _shProcess1 = 0x0;
    _shProcess2 = 0x0;
    // _photon1 = 0x0;
    // _photon2 = 0x0;
    // _shower1 = 0x0;
    // _shower2 = 0x0;
    return;
  }
  // // Use truth utilities to see which shower matches the first photon best.
  // for(auto const& p : truthUtils.GetRecoShowerListFromMCParticle(*_photon1, _evt, _showerLabel)) {
  //   if(p.first->Length() - shower.Length() < 1e-5) {
  //     _shower1 = &shower;
  //     _shower2 = shMatch;
  //     break;
  //   } else if(p.first->Length() - shMatch->Length() < 1e-5) {
  //     _shower1 = shMatch;
  //     _shower2 = &shower;
  //     break;
  //   }
  // }
}

// Functions to check whether data have been set (disallow partial reco).
bool PiZeroProcess::allMCSet() const {
  return (pi0() != 0x0) && (photon1() != 0x0) && (photon2() != 0x0);
}
bool PiZeroProcess::allRecoSet() const {
  return (shower1() != 0x0) && (shower2() != 0x0);
}

// Function to set the pi0 and photons.
void PiZeroProcess::setPi0(const simb::MCParticle& newPi0) {
  _pi0 = &newPi0;
  if (_pi0->NumberDaughters() == 2) { // Skip rare decays for now.
    // Find the associated photons through particle inventory.
    const simb::MCParticle* photon1 = pi_serv->TrackIdToParticle_P(_pi0->Daughter(0));
    const simb::MCParticle* photon2 = pi_serv->TrackIdToParticle_P(_pi0->Daughter(1));
    // Set photon1 to be the more energetic one.
    if(photon1->E() < photon2->E()) std::swap(photon1, photon2);
    // Fill the showerProcesses.
    _shProcess1 = std::make_unique<ShowerProcess>(*photon1, _evt);
    _shProcess2 = std::make_unique<ShowerProcess>(*photon2, _evt);
  }
}

} // namespace pizero

#endif
