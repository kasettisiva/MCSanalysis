////////////////////////////////////////////////////////////////////////
// Class:       ClusterCounter2
// Module Type: analyzer
// File:        ClusterCounter2_module.cc
//
// Access clusters, fill ROOT tree with a simple info from clusters.
// Robert Sulej
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Cluster.h"

#include "TTree.h"

namespace tutorial {

class ClusterCounter2 : public art::EDAnalyzer {
public:
  explicit ClusterCounter2(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  ClusterCounter2(ClusterCounter2 const &) = delete;
  ClusterCounter2(ClusterCounter2 &&) = delete;
  ClusterCounter2 & operator = (ClusterCounter2 const &) = delete;
  ClusterCounter2 & operator = (ClusterCounter2 &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.  
  void beginJob() override;
  void endJob() override;

private:
  size_t fEvNumber;

  TTree *fEventTree;
  size_t fNClusters;

  TTree *fClusterTree;
  size_t fNHits;

  // ******* fcl parameters *******
  art::InputTag fClusterModuleLabel;
  size_t fMinSize;
};

ClusterCounter2::ClusterCounter2(fhicl::ParameterSet const & p) : EDAnalyzer(p)
{
    fClusterModuleLabel = p.get< std::string >("ClusterModuleLabel");
    fMinSize = p.get< size_t >("MinSize");
}

void ClusterCounter2::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs; // TTree's are created in the memory managed by ROOT (you don't delete them)

    fEventTree = tfs->make<TTree>("EventTree", "event by event info");
    fEventTree->Branch("event", &fEvNumber, "fEvNumber/I");
    fEventTree->Branch("nclusters", &fNClusters, "fNClusters/I");

    fClusterTree = tfs->make<TTree>("ClusterTree", "cluster by cluster info");
    fClusterTree->Branch("event", &fEvNumber, "fEvNumber/I");
    fClusterTree->Branch("nhits", &fNHits, "fNHits/I");
}

void ClusterCounter2::analyze(art::Event const & evt)
{
    fEvNumber = evt.id().event();
    mf::LogVerbatim("ClusterCounter2") << "ClusterCounter2 module on event #" << fEvNumber;

    // use auto to make the line shorter when you remember all art types,
    // the full type being de-referenced here is: art::ValidHandle< std::vector<recob::Cluster> >
    auto const & clusters = *evt.getValidHandle< std::vector<recob::Cluster> >(fClusterModuleLabel);

    fNClusters = 0;

    // if you are old c++ granpa (or need cluster index, which may indeed happen!):
    // for (size_t i = 0; i < clusters.size(); ++i)
    // or:
    for (auto const & clu : clusters) // loop over recob::Cluster's stored in the std::vector 
    {
        fNHits = clu.NHits();
        fClusterTree->Fill();

        if (fNHits >= fMinSize) { ++fNClusters; }
    }
    fEventTree->Fill();
}

void ClusterCounter2::endJob()
{
    mf::LogVerbatim("ClusterCounter2") << "ClusterCounter2 finished job";
}

} // tutorial namespace

DEFINE_ART_MODULE(tutorial::ClusterCounter2)
