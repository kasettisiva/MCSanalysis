////////////////////////////////////////////////////////////////////////
// Class:       ClusterCounter
// Module Type: analyzer
// File:        ClusterCounter_module.cc
//
// Just an empty module, outputs the event number.
// Robert Sulej
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace tutorial {

class ClusterCounter : public art::EDAnalyzer {
public:
  explicit ClusterCounter(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  ClusterCounter(ClusterCounter const &) = delete;
  ClusterCounter(ClusterCounter &&) = delete;
  ClusterCounter & operator = (ClusterCounter const &) = delete;
  ClusterCounter & operator = (ClusterCounter &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:
};

ClusterCounter::ClusterCounter(fhicl::ParameterSet const & p) : EDAnalyzer(p) {}

void ClusterCounter::analyze(art::Event const & e)
{
    std::cout << "ClusterCounter module on event #" << e.id().event() << std::endl;
}

} // tutorial namespace

DEFINE_ART_MODULE(tutorial::ClusterCounter)
