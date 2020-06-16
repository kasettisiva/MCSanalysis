#ifndef PROTODUNE_SHOWER_UTILS_H
#define PROTODUNE_SHOWER_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNEShowerUtils
//  - Class to help analysers access useful shower information
// 
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "art/Framework/Principal/Event.h"

namespace protoana {

  class ProtoDUNEShowerUtils {

  public:

    ProtoDUNEShowerUtils();
    ~ProtoDUNEShowerUtils();

    /// Get the hits from a given reco shower
    const std::vector<const recob::Hit*> GetRecoShowerHits(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const;
    const std::vector<art::Ptr<recob::Hit>> GetRecoShowerArtHits(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const;

    /// Get the number of hits from a given reco shower
    unsigned int GetNumberRecoShowerHits(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const;
    /// Get the associated PCAxis object (from a principal component analysis)
    std::vector<const recob::PCAxis*> GetRecoShowerPCAxis(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const;
    // Estimate the energy of a shower from its hits
    std::vector<double> EstimateEnergyFromHitCharge(const std::vector<const recob::Hit*> &hits, calo::CalorimetryAlg caloAlg);
    /// If the shower.ID() isn't filled we must find the actual shower index ourselves
    int GetShowerIndex(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const;
    /// Get shower calo info
    std::vector<anab::Calorimetry> GetRecoShowerCalorimetry(const recob::Shower &shower, art::Event const &evt, const std::string showerModule, const std::string caloModule) const;
  private:


  };

}

#endif

