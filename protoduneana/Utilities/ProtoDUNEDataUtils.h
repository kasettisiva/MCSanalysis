#ifndef PROTODUNE_DATA_UTILS_H
#define PROTODUNE_DATA_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNEDataUtils
//  - Class to help analysers access useful beam data 
//    information
// 
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include <set>
#include "RtypesCore.h"
#include <stdint.h>

namespace protoana {

  class ProtoDUNEDataUtils {

  public:

    ProtoDUNEDataUtils(fhicl::ParameterSet const &pset);
    ~ProtoDUNEDataUtils();

    void reconfigure(fhicl::ParameterSet const &pset);

    /**
     * Returns true if the ProtoDUNE trigger says this is a beam trigger
     */
    bool IsBeamTrigger(art::Event const & evt) const;

    /// Get number of active fembs in an APA
    int GetNActiveFembsForAPA(art::Event const & evt, int apa) const;

    /// Check for consistency of timestamp values for a set of APAs.  True if consistent, false if there are mismatches
    bool CheckTimeStampConsistencyForAPAs(art::Event const & evt, std::set<int> apas, 
					  ULong64_t &timestamp, ULong64_t &timestamp2,
					  int &apainconsist) const;

  private:

    art::InputTag fTimingTag;
    art::InputTag fRawDigitTag;
    art::InputTag fRawDigitTimeStampTag;

  };

}

#endif

