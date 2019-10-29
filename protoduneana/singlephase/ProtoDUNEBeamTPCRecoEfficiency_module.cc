// leigh.howard.whitehead@cern.ch
// A simple module to determine the efficiency of reconstructing
// a beam particle using Pandora in events with a beam trigger

#include <iostream>
#include <utility>
#include <set>

#include "art/Framework/Core/EDAnalyzer.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 

#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "protoduneana/Analysis/ProtoDUNEPFParticleUtils.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"

namespace protoana{

  class ProtoDUNEBeamTPCRecoEfficiency : public art::EDAnalyzer {
  public:
    explicit ProtoDUNEBeamTPCRecoEfficiency(fhicl::ParameterSet const & pset);
    virtual ~ProtoDUNEBeamTPCRecoEfficiency() {};
    void analyze(art::Event const &evt) override;
    virtual void endJob() override;
  private:
    std::string fParticleLabel;
    unsigned int fBeamTriggers;
    unsigned int fBeamParticles;

    ProtoDUNEDataUtils fDataUtils;
  };

  ProtoDUNEBeamTPCRecoEfficiency::ProtoDUNEBeamTPCRecoEfficiency(fhicl::ParameterSet const & pset): 
    EDAnalyzer(pset),
    fDataUtils(pset.get<fhicl::ParameterSet>("DataUtils"))
  {
    fParticleLabel = pset.get<std::string>("ParticleLabel");
    fBeamTriggers = 0;
    fBeamParticles = 0;
  }

  void ProtoDUNEBeamTPCRecoEfficiency::analyze(art::Event const &evt) {

    ProtoDUNEPFParticleUtils pfpUtil;

    // Is this event from a beam trigger?
    if(fDataUtils.IsBeamTrigger(evt)){
      ++fBeamTriggers;
    }

    // Do we have a reconstructed beam slice from Pandora?
    if(pfpUtil.GetBeamSlice(evt,fParticleLabel) != 9999){
      ++fBeamParticles;
    }

  }

  void ProtoDUNEBeamTPCRecoEfficiency::endJob(){
    std::cout << "Beam triggered particle reconstruction efficiency = " << fBeamParticles/static_cast<float>(fBeamTriggers) 
              << " (" << fBeamParticles << "/" << fBeamTriggers << ")" << std::endl;
  }

  DEFINE_ART_MODULE(ProtoDUNEBeamTPCRecoEfficiency)

}
