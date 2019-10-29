#ifndef ProtoDUNEBeamCuts_h
#define ProtoDUNEBeamCuts_h

#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "art/Framework/Principal/Event.h"

#include <map>

namespace protoana{
 
  struct BeamVals{
    double X;
    double Y;
    double DirX;
    double DirY;
    double DirZ;
    bool Valid;
  };

  class ProtoDUNEBeamCuts{
   
    public:
      ProtoDUNEBeamCuts(){};
      ProtoDUNEBeamCuts(const fhicl::ParameterSet & pset);
      bool IsBeamlike( const recob::Track &, const art::Event &, std::string );
      bool IsBeamlike( const recob::Shower &, const art::Event &, std::string );

    private:

      BeamVals GetMCBeam( const art::Event & );
      BeamVals GetDataBeam( const art::Event & );

      std::map< std::string, fhicl::ParameterSet > DataCuts;
      std::map< std::string, fhicl::ParameterSet > MCCuts;

      std::vector< std::string > valid_momenta = {
        ".3", ".5", "1", "2", "3", "6", "7"
      };
   
  };
}

#endif
