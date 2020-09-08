#ifndef G4REWEIGHTUTILS_h
#define G4REWEIGHTUTILS_h

#include <vector>
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"

namespace protoana {
 namespace G4ReweightUtils {

  //Function to create reweightable object
  bool CreateRWTraj(const simb::MCParticle & part,
                    const sim::ParticleList & plist,
                    art::ServiceHandle < geo::Geometry > geo_serv, int event,
                    G4ReweightTraj * theTraj); 

  std::vector<G4ReweightTraj *> CreateNRWTrajs(
      const simb::MCParticle & part,
      const sim::ParticleList & plist,
      art::ServiceHandle < geo::Geometry > geo_serv, int event,
      bool fVerbose=false);
 }
}

#endif
