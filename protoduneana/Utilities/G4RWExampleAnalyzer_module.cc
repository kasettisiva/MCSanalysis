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
#include "geant4reweight/src/ReweightBase/G4ReweightStep.hh"
#include "geant4reweight/src/PropBase/G4ReweightParameterMaker.hh"
#include "geant4reweight/src/ReweightBase/G4MultiReweighter.hh"

#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace protoana {
  class G4RWExampleAnalyzer;
}


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

  //Function to create reweightable object
  bool CreateRWTraj(const simb::MCParticle & part,
                    const sim::ParticleList & plist,
                    art::ServiceHandle < geo::Geometry > geo_serv, int event,
                    G4ReweightTraj * theTraj);

  // Declare member data here.
  TTree *fTree;
  int run, subrun, event;
  int true_beam_PDG;
  int true_beam_ID;
  std::vector<double> g4rw_primary_plus_sigma_weight;
  std::vector<double> g4rw_primary_minus_sigma_weight;
  std::vector<double> g4rw_primary_weights;
  std::vector<std::string> g4rw_primary_var;

  
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
      ParMaker(ParSet),
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
  event = e.id().event();
  run = e.run();
  subrun = e.subRun();

  std::cout << "Doing reweight" << std::endl;
  if (true_beam_PDG == 211) {
    G4ReweightTraj theTraj(true_beam_ID, true_beam_PDG, 0, event, {0,0});
    bool created = CreateRWTraj(*true_beam_particle, plist,
                                fGeometryService, event, &theTraj);
    if (created) {
      g4rw_primary_weights.push_back(MultiRW.GetWeightFromNominal(theTraj));
      
      std::vector<double> weights_vec = MultiRW.GetWeightFromAll1DThrows(
          theTraj);
      g4rw_primary_weights.insert(g4rw_primary_weights.end(),
                                  weights_vec.begin(), weights_vec.end());


      //g4rw_primary_plus_sigma_weight = pm_weights.first;
      //g4rw_primary_minus_sigma_weight = pm_weights.second;

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

bool protoana::G4RWExampleAnalyzer::CreateRWTraj(
    const simb::MCParticle & part, const sim::ParticleList & plist,
    art::ServiceHandle < geo::Geometry > geo_serv, int event,
    G4ReweightTraj * theTraj) {

  //Loop over daughters
  for (int i = 0; i < part.NumberDaughters(); ++i) {
    int d_index = part.Daughter(i);
    auto d_part = plist[d_index];
    
    int d_PDG = d_part->PdgCode();
    int d_ID = d_part->TrackId();

    theTraj->AddChild(new G4ReweightTraj(d_ID, d_PDG, part.TrackId(),
                      event, {0,0}));
  }

  //Create process map
  auto procs = part.Trajectory().TrajectoryProcesses();
  std::map<size_t, std::string> proc_map;
  for (auto it = procs.begin(); it != procs.end(); ++it) {
    proc_map[it->first] = part.Trajectory().KeyToProcess(it->second);
  }

  std::vector<double> traj_X, traj_Y, traj_Z;
  std::vector<double> traj_PX, traj_PY, traj_PZ;
  std::vector<size_t> elastic_indices;

  bool found_last = false;
  //G4ReweightTraj theTraj(part.TrackId(), part.PdgCode(), 0, event, {0,0});
  for (size_t i = 0; i < part.NumberTrajectoryPoints(); ++i) {
    double x = part.Position(i).X();
    double y = part.Position(i).Y();
    double z = part.Position(i).Z();
    
    geo::Point_t test_point{x, y, z};
    const TGeoMaterial * test_material = geo_serv->Material(test_point);

    if (!strcmp(test_material->GetName(), "LAr")) {
      traj_X.push_back(x);
      traj_Y.push_back(y);
      traj_Z.push_back(z);

      traj_PX.push_back(part.Px(i));
      traj_PY.push_back(part.Py(i));
      traj_PZ.push_back(part.Pz(i));

      auto itProc = proc_map.find(i);
      if (itProc != proc_map.end() && itProc->second == "hadElastic") {
        elastic_indices.push_back(i);
      }
    }

    if (i == part.NumberTrajectoryPoints() - 1)
      found_last = true;
  }

  double mass = 0.;

  switch (abs(part.PdgCode())) {
    case 211: {
      mass = 139.57;
      break;
    }
    case 2212: {
      mass = 938.28;
      break;
    }
    default: {
      return false;
      break;
    }
  }

  for (size_t i = 1; i < traj_X.size(); ++i) {
    std::string proc = "default";
    if (found_last && i == traj_X.size() - 1) {
      proc = part.EndProcess();
    }
    else if (std::find(elastic_indices.begin(), elastic_indices.end(), i) !=
             elastic_indices.end()){
      proc = "hadElastic";
    }

    double dX = traj_X[i] - traj_X[i-1];
    double dY = traj_Y[i] - traj_Y[i-1];
    double dZ = traj_Z[i] - traj_Z[i-1];

    double len = sqrt(dX*dX + dY*dY + dZ*dZ);

    double preStepP[3] = {traj_PX[i-1]*1.e3, 
                          traj_PY[i-1]*1.e3, 
                          traj_PZ[i-1]*1.e3};

    double postStepP[3] = {traj_PX[i]*1.e3, 
                           traj_PY[i]*1.e3, 
                           traj_PZ[i]*1.e3};
    if (i == 1) {
      double p_squared = preStepP[0]*preStepP[0] + preStepP[1]*preStepP[1] +
                         preStepP[2]*preStepP[2];
      theTraj->SetEnergy(sqrt(p_squared + mass*mass));
    }

    G4ReweightStep * step = new G4ReweightStep(part.TrackId(), part.PdgCode(),
                                               0, event, preStepP, postStepP,
                                               len, proc);
    theTraj->AddStep(step);
  }

  return true;
}

void protoana::G4RWExampleAnalyzer::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("beamana","beam analysis tree");

  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("true_beam_ID", &true_beam_ID);
  fTree->Branch("true_beam_PDG", &true_beam_PDG);

  fTree->Branch("g4rw_primary_weights", &g4rw_primary_weights);
  fTree->Branch("g4rw_primary_plus_sigma_weight", &g4rw_primary_plus_sigma_weight);
  fTree->Branch("g4rw_primary_minus_sigma_weight", &g4rw_primary_minus_sigma_weight);
  fTree->Branch("g4rw_primary_var", &g4rw_primary_var);
}

void protoana::G4RWExampleAnalyzer::reset() {
  true_beam_PDG = -1;
  true_beam_ID = -1;
  
  g4rw_primary_weights.clear();
  g4rw_primary_plus_sigma_weight.clear();
  g4rw_primary_minus_sigma_weight.clear();
  g4rw_primary_var.clear();
}
DEFINE_ART_MODULE(protoana::G4RWExampleAnalyzer)
