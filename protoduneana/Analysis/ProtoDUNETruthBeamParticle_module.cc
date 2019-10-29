////////////////////////////////////////////////////////////////////////
//
// Example how to read and save the beam truth particle.
// Georgios Christodoulou 
// Georgios.Christodoulou@cern.ch
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
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

// C++ Includes
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

const int NMAXBEAMPARTICLES = 500;

namespace protoana{
  class ProtoDUNETruthBeamParticle;
}

class protoana::ProtoDUNETruthBeamParticle : public art::EDAnalyzer{
public:
  
  explicit ProtoDUNETruthBeamParticle(fhicl::ParameterSet const& pset);
  virtual ~ProtoDUNETruthBeamParticle();
  
  void beginJob();
  void analyze(const art::Event& evt);
  void reconfigure(fhicl::ParameterSet const& pset);
  void endJob();
  
private:
  
  int nbeamparticles;
  int run[NMAXBEAMPARTICLES];
  int subrun[NMAXBEAMPARTICLES];
  int event[NMAXBEAMPARTICLES];
  int pdg[NMAXBEAMPARTICLES];
  int goodparticle[NMAXBEAMPARTICLES];
  float mom[NMAXBEAMPARTICLES];
  float ene[NMAXBEAMPARTICLES];
  float startpos[NMAXBEAMPARTICLES][4];
  float p[NMAXBEAMPARTICLES][3];

  TTree *truth;
  
  // Parameters from datacard
  std::vector<int> p_pdg;
  bool p_savegoodparticle; 
  bool selectall;
  
};
  
//-----------------------------------------------------------------------
protoana::ProtoDUNETruthBeamParticle::ProtoDUNETruthBeamParticle(fhicl::ParameterSet const& pset) : EDAnalyzer(pset){

  this->reconfigure(pset);
  
}

//-----------------------------------------------------------------------
protoana::ProtoDUNETruthBeamParticle::~ProtoDUNETruthBeamParticle(){}

//-----------------------------------------------------------------------
void protoana::ProtoDUNETruthBeamParticle::beginJob() {
  
  art::ServiceHandle<art::TFileService> tfs;
  
  // Define output tree
  truth = tfs->make<TTree>("truth", "Beam Particle Truth Tree");
  truth->Branch("nbeamparticles",   &nbeamparticles,       "nbeamparticles/I");
  truth->Branch("run",              run,                   "run[nbeamparticles]/I");
  truth->Branch("subrun",           subrun,                "subrun[nbeamparticles]/I");
  truth->Branch("event",            event,                 "event[nbeamparticles]/I");
  truth->Branch("pdg",              pdg,                   "pdg[nbeamparticles]/I");
  truth->Branch("goodparticle",     goodparticle,          "goodparticle[nbeamparticles]/I");
  truth->Branch("mom",              mom,                   "mom[nbeamparticles]/F");
  truth->Branch("ene",              ene,                   "ene[nbeamparticles]/F");
  truth->Branch("startpos",         startpos,              "startpos[nbeamparticles][4]/F");
  truth->Branch("p",                p,                     "p[nbeamparticles][3]/F");
  
}

//-----------------------------------------------------------------------
void protoana::ProtoDUNETruthBeamParticle::reconfigure(fhicl::ParameterSet const& pset){

  p_savegoodparticle  = pset.get<bool>("SaveOnlyGoodParticle");
  p_pdg               = pset.get< std::vector<int> >("Pdg");
  if(p_pdg[0] == 0)
    selectall = true;
  else
    selectall = false;

}

//-----------------------------------------------------------------------
void protoana::ProtoDUNETruthBeamParticle::analyze(const art::Event& evt){
    
  // Access truth info
  art::Handle< std::vector<simb::MCTruth> > mctruthhandle;
  std::vector< art::Ptr<simb::MCTruth> > mclist;
  if(evt.getByLabel("generator",mctruthhandle)){
    art::fill_ptr_vector(mclist,mctruthhandle);
  }
  else{
    mf::LogError("protoana::ProtoDUNETruthBeamParticle::analyze") << "Requested protoDUNE beam generator information does not exist!";
    return;
  }
  
  art::Ptr<simb::MCTruth> mctruth;
  mctruth = mclist[0];
  
  nbeamparticles = 0;
  for(int i = 0; i < mctruth->NParticles(); ++i){
    const simb::MCParticle& part(mctruth->GetParticle(i));

    // Select partciles with the chosen pdg
    if(!selectall){
      bool found = false;
      for(unsigned int j = 0; j < p_pdg.size(); j++){
	if(part.PdgCode() == p_pdg[j]){
	  found = true;
	  break;
	}
      }
      if(!found) continue;
    }   
    
    // Option to select only good beam particles
    bool isgoodparticle   = (part.Process() == "primary");
    if(p_savegoodparticle && !isgoodparticle) continue;

    nbeamparticles++;

    run[i]            = evt.run();
    subrun[i]         = evt.subRun();
    event[i]          = evt.id().event();
    goodparticle[i]   = (part.Process() == "primary");
    startpos[i][0]    = part.Vx();
    startpos[i][1]    = part.Vy();
    startpos[i][2]    = part.Vz();
    startpos[i][3]    = part.T();
    p[i][0]           = part.Px();
    p[i][1]           = part.Py();
    p[i][2]           = part.Pz();
    mom[i]            = part.P();
    ene[i]            = part.E();
    pdg[i]            = part.PdgCode();
  }
  
  if(nbeamparticles > 0)
    truth->Fill();
  
}

void protoana::ProtoDUNETruthBeamParticle::endJob() {}
 
DEFINE_ART_MODULE(protoana::ProtoDUNETruthBeamParticle)
