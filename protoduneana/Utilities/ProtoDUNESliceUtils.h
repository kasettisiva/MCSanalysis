#ifndef PROTODUNE_SLICE_UTILS_H
#define PROTODUNE_SLICE_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNESliceUtils
//  - Class to help analysers access useful slice information
// 
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////

#include <map>

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"

#include "art/Framework/Principal/Event.h"

namespace protoana {

  class ProtoDUNESliceUtils {

  public:

    ProtoDUNESliceUtils();
    ~ProtoDUNESliceUtils();

    // Get the hits from the slice (using the slice object)
    const std::vector<const recob::Hit*> GetRecoSliceHits(const recob::Slice &slice, art::Event const &evt, const std::string sliceModule) const;
    const std::vector<const recob::Hit*> GetRecoSliceHits(unsigned int sliceID, art::Event const &evt, const std::string sliceModule) const;

    // Number of slice hits
    unsigned int GetNumberRecoSliceHits(const recob::Slice &slice, art::Event const &evt, const std::string sliceModule) const;
    unsigned int GetNumberRecoSliceHits(const unsigned int sliceID, art::Event const &evt, const std::string sliceModule) const;

    // A map of all hits in each slice
    const std::map<unsigned int, std::vector<const recob::Hit*>> GetRecoSliceHitMap(art::Event const &evt, const std::string sliceModule) const;
  
  private:


  };

}

#endif

