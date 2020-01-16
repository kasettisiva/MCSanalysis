////////////////////////////////////////////////////////////////////////
//
// Filter to selects events with a given beam particle.
// Owen Goodwin
// owen.goodwin@postgrad.manchester.ac.uk
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"




// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
// C++ Includes
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>



namespace protoana{
  class ProtoDUNETruthBeamFilter;
}

class protoana::ProtoDUNETruthBeamFilter : public art::EDFilter{
public:
  
  explicit ProtoDUNETruthBeamFilter(fhicl::ParameterSet const& pset);
  virtual ~ProtoDUNETruthBeamFilter();
  
  void beginJob();
  bool filter(art::Event& evt); 
  void reconfigure(fhicl::ParameterSet const& pset);
  void endJob(); 
  
private:
  
  bool fDebug;                  /// Controls the verbosity
  
  std::vector<int> fPDG;      /// List of particle PDGs we want to keep 

 
 
  TH1D* fSelectedEvents;
  TH1D* fTotalEvents;


  
};
  
//-----------------------------------------------------------------------
protoana::ProtoDUNETruthBeamFilter::ProtoDUNETruthBeamFilter(fhicl::ParameterSet const& pset) 
: EDFilter(pset)
{
  this->reconfigure(pset);

}

//-----------------------------------------------------------------------
protoana::ProtoDUNETruthBeamFilter::~ProtoDUNETruthBeamFilter(){}

//-----------------------------------------------------------------------
void protoana::ProtoDUNETruthBeamFilter::beginJob() {

	art::ServiceHandle<art::TFileService> tfs;
    fSelectedEvents = tfs->make<TH1D>("fSelectedEvents", "Number of Selected Events", 3, 0, 3); //counts the number of selected events 
    fTotalEvents = tfs->make<TH1D>("fTotalEvents", "Total Events", 3, 0, 3); //counts the initial number of events in the unfiltered root input file
  
}

//-----------------------------------------------------------------------
void protoana::ProtoDUNETruthBeamFilter::reconfigure(fhicl::ParameterSet const& pset){

  
  fPDG             = pset.get< std::vector<int> >("PDG");
  fDebug             = pset.get<bool>("IsVerbose"); 

}

//-----------------------------------------------------------------------
bool protoana::ProtoDUNETruthBeamFilter::filter(art::Event& evt){



  if(fDebug) std::cout << "Reading Event" << std::endl;
    

  bool found = false;  
  // Access truth info
  art::Handle< std::vector<simb::MCTruth> > mctruthhandle;
  std::vector< art::Ptr<simb::MCTruth> > mclist;
  if(evt.getByLabel("generator",mctruthhandle)){
    art::fill_ptr_vector(mclist,mctruthhandle);
  }
  else{
    mf::LogError("protoana::ProtoDUNETruthBeamParticle::analyze") << "Requested protoDUNE beam generator information does not exist!";
    return found;
  }
  fTotalEvents->Fill(1); //count total events

  art::Ptr<simb::MCTruth> mctruth;
  mctruth = mclist[0];
  

  for(int i = 0; i < mctruth->NParticles(); ++i){
    const simb::MCParticle& part(mctruth->GetParticle(i));
  
	  bool isgoodparticle   = (part.Process() == "primary" && mctruth->Origin()==4); //primary particle and beam origin
	  if(!isgoodparticle) continue;
    // Select partciles with the chosen pdg
    for(unsigned int j = 0; j < fPDG.size(); j++){
		  if(part.PdgCode() == fPDG[j]){
	  	  found = true;
	  	  if(fDebug) std::cout << "this is a selected event with PDG "<< part.PdgCode() << std::endl;
	  	  fSelectedEvents->Fill(1); //count selected events
	 	    break;
	    }
    }
    
  }   
  
  return found;
  
}


void protoana::ProtoDUNETruthBeamFilter::endJob() {}
 
DEFINE_ART_MODULE(protoana::ProtoDUNETruthBeamFilter)

