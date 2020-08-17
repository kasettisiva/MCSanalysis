////////////////////////////////////////////////////////////////////////
// Class:       TrueEdep
// Plugin Type: producer (art v3_05_01)
// File:        TrueEdep_module.cc
//
// Generated at Wed Jul  1 23:49:17 2020 by Tingjun Yang using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include <memory>

namespace pdsp {
  class TrueEdep;
}


class pdsp::TrueEdep : public art::EDProducer {
public:
  explicit TrueEdep(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrueEdep(TrueEdep const&) = delete;
  TrueEdep(TrueEdep&&) = delete;
  TrueEdep& operator=(TrueEdep const&) = delete;
  TrueEdep& operator=(TrueEdep&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  double nelectrons;

};


pdsp::TrueEdep::TrueEdep(fhicl::ParameterSet const& p)
  : EDProducer{p},
  nelectrons(p.get< double >("Nelectrons"))
{
  produces< std::vector<sim::SimEnergyDeposit> >();
}

void pdsp::TrueEdep::produce(art::Event& e)
{
  std::unique_ptr<std::vector<sim::SimEnergyDeposit> > edepos(new std::vector<sim::SimEnergyDeposit>);

  geo::Point_t start={-355,300,120};
  geo::Point_t end={-355,300,120};
  edepos->emplace_back(sim::SimEnergyDeposit(0,nelectrons,0,10,start,end));
  e.put(std::move(edepos));
}

DEFINE_ART_MODULE(pdsp::TrueEdep)
