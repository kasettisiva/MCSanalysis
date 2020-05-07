////////////////////////////////////////////////////////////////////////
// Class:       G4RWExampleAnalyzer
// Plugin Type: analyzer (art v3_05_00)
// File:        G4RWExampleAnalyzer_module.cc
//
// Generated at Wed Apr 22 10:59:17 2020 by Jacob Calcutt using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "TTree.h"

#include "geant4reweight/src/ReweightBase/G4ReweighterFactory.hh"
#include "geant4reweight/src/ReweightBase/G4Reweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/src/PropBase/G4ReweightParameterMaker.hh"
#include "geant4reweight/src/ReweightBase/G4MultiReweighter.hh"
#include "G4ReweightUtils.h"

#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"


namespace protoana {
  class G4RWExampleAnalyzer;
}

using protoana::G4ReweightUtils::CreateRWTraj;

class protoana::G4RWExampleAnalyzer : public art::EDAnalyzer {
public:
  explicit G4RWExampleAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  G4RWExampleAnalyzer(G4RWExampleAnalyzer const&) = delete;
  G4RWExampleAnalyzer(G4RWExampleAnalyzer&&) = delete;
  G4RWExampleAnalyzer& operator=(G4RWExampleAnalyzer const&) = delete;
  G4RWExampleAnalyzer& operator=(G4RWExampleAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  //Optional
  void beginJob() override;
  void reset();

private:

  // Declare member data here.
  TTree *fTree;
  int run, subrun, event;
  int true_beam_PDG;
  int true_beam_ID;
  double true_beam_len;
  int true_beam_nElasticScatters;
  std::vector<double> g4rw_primary_plus_sigma_weight;
  std::vector<double> g4rw_primary_minus_sigma_weight;
  std::vector<double> g4rw_primary_weights;
  std::vector<std::string> g4rw_primary_var;
  double g4rw_primary_singular_weight;

  
  std::string fGeneratorTag;
  //Geant4Reweight stuff
  int RW_PDG;
  TFile FracsFile, XSecFile;
  std::vector<fhicl::ParameterSet> ParSet;
  G4ReweightParameterMaker ParMaker;
  G4MultiReweighter MultiRW;
  G4ReweighterFactory RWFactory;
  G4Reweighter * theRW;

};


protoana::G4RWExampleAnalyzer::G4RWExampleAnalyzer(
    fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
      fGeneratorTag(p.get<std::string>("GeneratorTag")),
      RW_PDG(p.get<int>("RW_PDG")),
      FracsFile( (p.get< std::string >( "FracsFile" )).c_str(), "OPEN" ),
      XSecFile( (p.get< std::string >( "XSecFile" )).c_str(), "OPEN"),
      ParSet(p.get<std::vector<fhicl::ParameterSet>>("ParameterSet")),
      ParMaker(ParSet, RW_PDG),
      MultiRW(RW_PDG, XSecFile, FracsFile, ParSet) {

  theRW = RWFactory.BuildReweighter(RW_PDG, &XSecFile, &FracsFile,
                                    ParMaker.GetFSHists(),
                                    ParMaker.GetElasticHist()/*, true*/ );
}

void protoana::G4RWExampleAnalyzer::analyze(art::Event const& e) {

  reset();

  protoana::ProtoDUNETruthUtils truthUtil;
  art::ServiceHandle < geo::Geometry > fGeometryService;
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  const sim::ParticleList & plist = pi_serv->ParticleList(); 


  // This gets the true beam particle that generated the event
  const simb::MCParticle* true_beam_particle = 0x0;
  auto mcTruths = e.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
  true_beam_particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0], e);
  if (!true_beam_particle) {
    MF_LOG_WARNING("PionAnalyzer") << "No true beam particle" << std::endl;
    return;
  }
  ////////////////////////////

  true_beam_PDG = true_beam_particle->PdgCode();
  true_beam_ID = true_beam_particle->TrackId();
  true_beam_len = true_beam_particle->Trajectory().TotalLength();

  const simb::MCTrajectory & true_beam_trajectory =
      true_beam_particle->Trajectory();
  auto true_beam_proc_map = true_beam_trajectory.TrajectoryProcesses();
 
  for (auto itProc = true_beam_proc_map.begin();
       itProc != true_beam_proc_map.end(); ++itProc) {
    //int index = itProc->first;
    std::string process = true_beam_trajectory.KeyToProcess(itProc->second);

    if (process == "hadElastic") {
      ++true_beam_nElasticScatters;
    }
  }

  event = e.id().event();
  run = e.run();
  subrun = e.subRun();

  std::cout << "Doing reweight" << std::endl;
  if (true_beam_PDG == RW_PDG) {
    G4ReweightTraj theTraj(true_beam_ID, true_beam_PDG, 0, event, {0,0});
    bool created = CreateRWTraj(*true_beam_particle, plist,
                                fGeometryService, event, &theTraj);
    if (created && theTraj.GetNSteps()) {

      g4rw_primary_singular_weight = MultiRW.GetWeightFromNominal(theTraj);
      //the following method achieves the same result
      //g4rw_primary_singular_weight = theRW->GetWeight(&theTraj);
      
      std::vector<double> weights_vec = MultiRW.GetWeightFromAll1DThrows(
          theTraj);
      g4rw_primary_weights.insert(g4rw_primary_weights.end(),
                                  weights_vec.begin(), weights_vec.end());

      for (size_t i = 0; i < ParSet.size(); ++i) {
        std::pair<double, double> pm_weights =
            MultiRW.GetPlusMinusSigmaParWeight(theTraj, i);

        g4rw_primary_plus_sigma_weight.push_back(pm_weights.first);
        g4rw_primary_minus_sigma_weight.push_back(pm_weights.second);
        g4rw_primary_var.push_back(ParSet[i].get<std::string>("Name"));
      }

    }
  }

  fTree->Fill();
}

void protoana::G4RWExampleAnalyzer::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree","output tree");

  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("true_beam_ID", &true_beam_ID);
  fTree->Branch("true_beam_PDG", &true_beam_PDG);
  fTree->Branch("true_beam_len", &true_beam_len);
  fTree->Branch("true_beam_nElasticScatters", &true_beam_nElasticScatters);

  fTree->Branch("g4rw_primary_weights", &g4rw_primary_weights);
  fTree->Branch("g4rw_primary_singular_weight", &g4rw_primary_singular_weight);
  fTree->Branch("g4rw_primary_plus_sigma_weight", &g4rw_primary_plus_sigma_weight);
  fTree->Branch("g4rw_primary_minus_sigma_weight", &g4rw_primary_minus_sigma_weight);
  fTree->Branch("g4rw_primary_var", &g4rw_primary_var);
}

void protoana::G4RWExampleAnalyzer::reset() {
  true_beam_PDG = -1;
  true_beam_ID = -1;
  true_beam_len = -1.;
  true_beam_nElasticScatters = 0.;
  
  g4rw_primary_weights.clear();
  g4rw_primary_singular_weight = 1.;
  g4rw_primary_plus_sigma_weight.clear();
  g4rw_primary_minus_sigma_weight.clear();
  g4rw_primary_var.clear();
}
DEFINE_ART_MODULE(protoana::G4RWExampleAnalyzer)
