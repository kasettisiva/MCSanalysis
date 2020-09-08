////////////////////////////////////////////////////////////////////////
// Class:       PDSPmatchmc
// Plugin Type: analyzer (art v3_01_02)
// File:        PDSPmatchmc_module.cc
//
// Generated at Wed May  1 16:53:19 2019 by Bryan Ramson using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include <numeric>
#include <sstream>
#include <string.h>
#include <bitset>
#include <vector>

#include "larana/OpticalDetector/OpFlashAlg.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//#include "art/Framework/Services/Optional/TFileService.h"
#include "art_root_io/TFileService.h"

#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "dune/Protodune/singlephase/CRT/data/CRTTrigger.h"
//#include "dune/Geometry/ProtoDUNESPCRTSorter.h"                                                                                                                                                    

#include "dune-raw-data/Overlays/CRTFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "nusimdata/SimulationBase/MCTruth.h"
//#include "nusimdata/SimulationBase/MCParticle.h"
//#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace pdsp {
  class PDSPmatchmc;
}


class pdsp::PDSPmatchmc : public art::EDAnalyzer {
public:
  explicit PDSPmatchmc(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPmatchmc(PDSPmatchmc const&) = delete;
  PDSPmatchmc(PDSPmatchmc&&) = delete;
  PDSPmatchmc& operator=(PDSPmatchmc const&) = delete;
  PDSPmatchmc& operator=(PDSPmatchmc&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  // Declare member data here.

  const art::InputTag fOpHitLabel;//MC light hits
  //const art::InputTag fPFParListLabel;//Label for MC PF Particle containers
  const art::InputTag fMCTruthLabel;//Label for MC Truth
  const art::InputTag fGEANTLabel;//GEANT Label for truth particles

  TTree *fTree;
  
  std::vector<double> fTrueID, fTrue_time, fTrueStartx, fTrueStarty, fTrueStartz, fTruePx, fTruePy, fTruePz;

  std::vector<int64_t> fPDS_time;
  std::vector<double> fOpChan, fPE;
};

pdsp::PDSPmatchmc::PDSPmatchmc(fhicl::ParameterSet const& p)
  :
  EDAnalyzer(p),
  fOpHitLabel(p.get<art::InputTag>("OpHitLabel_MC")),
  fMCTruthLabel(p.get<art::InputTag>("MCTruthLabel_MC")),
  fGEANTLabel(p.get<art::InputTag>("GEANTLabel_MC"))
  //: EDAnalyzer{p},
  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  consumes<std::vector<recob::OpHit>>(fOpHitLabel);
  consumes<std::vector<simb::MCTruth>>(fMCTruthLabel);
  consumes<std::vector<simb::MCParticle>>(fGEANTLabel);
}

void pdsp::PDSPmatchmc::analyze(art::Event const& e)
{
  art::ServiceHandle< cheat::PhotonBackTrackerService > pbt;

  bool op = false;
  // Implementation of required member function here.
  fTrueID.clear();
  fTrue_time.clear();
  fTrueStartx.clear();
  fTrueStarty.clear();
  fTrueStartz.clear();
  fTruePx.clear();
  fTruePy.clear();
  fTruePz.clear();
  
  fPDS_time.clear();
  fOpChan.clear();
  fPE.clear();

  const auto MClistHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);
  art::Ptr<simb::MCTruth> mctruth(MClistHandle,0);
  if(mctruth->NParticles() == 0){
    mf::LogError("PDSPmatchmc") << "No MCTruth Particles! Skipping.";
    return;
  }

  const auto OpHitHandle = e.getValidHandle<std::vector<recob::OpHit>>(fOpHitLabel);
 
  const simb::MCParticle& part(mctruth->GetParticle(0));
  fTrueID.push_back(part.PdgCode());
  fTrue_time.push_back(part.T());
  fTrueStartx.push_back(part.Vx());
  fTrueStarty.push_back(part.Vy());
  fTrueStartz.push_back(part.Vz());
  fTruePx.push_back(part.Px());
  fTruePy.push_back(part.Py());
  fTruePz.push_back(part.Pz());
  
  for(const auto& OpHit: *OpHitHandle){
    fPDS_time.push_back((OpHit.PeakTime()));
    fOpChan.push_back(OpHit.OpChannel());
    fPE.push_back(OpHit.PE());
  }
  if(fPDS_time.size() > 0) op=true;
  
  if(op) fTree->Fill();
  
}


void pdsp::PDSPmatchmc::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tFileService;
  
  fTree = tFileService->make<TTree>("ProtoDUNE_Evt_Match","ProtoDUNE Matched MC Event Tree");
  
  fTree->Branch("fTrueID",&fTrueID);
  fTree->Branch("fTrue_time",&fTrue_time);
  fTree->Branch("fTrueStartx",&fTrueStartx);
  fTree->Branch("fTrueStarty",&fTrueStarty);
  fTree->Branch("fTrueStartz",&fTrueStartz);
  fTree->Branch("fTruePx",&fTruePx);
  fTree->Branch("fTruePy",&fTruePy);
  fTree->Branch("fTruePz",&fTruePz);
   
  fTree->Branch("fPDS_time",&fPDS_time);
  fTree->Branch("fOpChan",&fOpChan);
  fTree->Branch("fPE",&fPE);
}

DEFINE_ART_MODULE(pdsp::PDSPmatchmc)
