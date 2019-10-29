// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

namespace beam{
  class ReadBeamInfo : public art::EDAnalyzer{
    public:

      struct Config {
        // save some typing:
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // one Atom for each parameter
        fhicl::Atom<art::InputTag> BeamLabel {
          Name("BeamLabel"),
          Comment("tag of the input data product with the detector simulation information")
        };
      }; // Config
      using Parameters = art::EDAnalyzer::Table<Config>;
      
//      explicit ReadBeamInfo(fhicl::ParameterSet const& pset);
      explicit ReadBeamInfo(Parameters const& pset);
      virtual ~ReadBeamInfo();

      // The analysis routine, called once per event. 
      void analyze (const art::Event& event) override;

    private:
      
      art::InputTag fBeamLabel; ///< The name of the beam event producer

  }; // class ReadBeamInfo

} // End the namespace.

//beam::ReadBeamInfo::ReadBeamInfo(fhicl::ParameterSet const& pset): art::EDAnalyzer(pset){
beam::ReadBeamInfo::ReadBeamInfo(ReadBeamInfo::Parameters const& pset): art::EDAnalyzer(pset){

//  fBeamLabel = pset.get<art::InputTag>("BeamLabel");
  fBeamLabel = pset().BeamLabel();

}

beam::ReadBeamInfo::~ReadBeamInfo(){

}

void beam::ReadBeamInfo::analyze(const art::Event& evt){

  // Just for MC!
  if(evt.isRealData()) return;

  // ProtoDUNE beam generator information
  art::Handle< std::vector<simb::MCTruth> > beamTruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > beamTruth;

  if(evt.getByLabel(fBeamLabel,beamTruthListHandle)){
    art::fill_ptr_vector(beamTruth,beamTruthListHandle);
  }
  else{
    // If we requested this info but it doesn't exist, don't use it.
    mf::LogError("ReadBeamInfo") << "Requested protoDUNE beam generator information with name " << fBeamLabel << " but none exists";
    return;
  }

  // There is only one MC truth, so use it to get the number of primaries
  Int_t nPrimaries = beamTruth[0]->NParticles();
 
  for(Int_t iPartp = 0; iPartp < nPrimaries; ++iPartp){
    const simb::MCParticle& partp(beamTruth[0]->GetParticle(iPartp));

    // Access the beam particle variables
    TVector3 beamVtx(partp.Vx(), partp.Vy(), partp.Vz());
    TVector3 beamMom(partp.Px(), partp.Py(), partp.Pz());
    TVector3 beamDir = beamMom.Unit();
    double beamTime = partp.T();
    double beamMomentum = partp.P();
    double beamEnergy = partp.E();
    int beamPDG = partp.PdgCode();

    // There are two types of beam particles...
    // The ones we trigger on, "GoodParticles" have Process() == "primary"
    // The background ones have Process() == "primaryBackground"
    // This is how we identify which ones are the good particle.
    bool beamIsGood = partp.Process()=="primary";

    std::string type = "good";
    if(!beamIsGood){
      type = "background";
    }

    std::cout << " -- Read " << type << " beam particle with pdg code " << beamPDG << " from event " << evt.event() << std::endl;

    if(beamIsGood){
      std::cout << " - Vtx    = (" << beamVtx.X() << ", " << beamVtx.Y() << ", " << beamVtx.Z() << ")" << std::endl; 
      std::cout << " - Time   = " << beamTime << std::endl;
      std::cout << " - Dir    = (" << beamDir.X() << ", " << beamDir.Y() << ", " << beamDir.Z() << ")" << std::endl; 
      std::cout << " - Mom    = (" << beamMom.X() << ", " << beamMom.Y() << ", " << beamMom.Z() << ")" << std::endl; 
      std::cout << " - |Mom|  = " << beamMomentum << std::endl;
      std::cout << " - Energy = " << beamEnergy << std::endl;
      std::cout << " - PDG    = " << beamPDG << std::endl;
    }

  } // End loop over the beam particles

}

namespace beam{

  DEFINE_ART_MODULE(ReadBeamInfo)

}

