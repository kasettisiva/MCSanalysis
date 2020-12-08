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
#include "geant4reweight/src/ReweightBase/G4ReweightManager.hh"
#include "G4ReweightUtils.h"

#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"


namespace protoana {
  class G4RWExampleAnalyzer;
}

using protoana::G4ReweightUtils::CreateRWTraj;
using protoana::G4ReweightUtils::CreateNRWTrajs;
using protoana::G4ReweightUtils::GetNTrajPMSigmaWeights;
using protoana::G4ReweightUtils::GetNTrajWeightFromSetPars;

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
  double true_beam_len, true_beam_startP, true_beam_endP;
  double true_beam_startX, true_beam_startY, true_beam_startZ,
         true_beam_endX, true_beam_endY, true_beam_endZ;
  int true_beam_nTrajPts;
  std::vector<double> true_beam_traj_dx;
  std::vector<double> true_beam_dEdX;
  std::vector<double> reco_beam_len,
                      reco_beam_startX, reco_beam_startY, reco_beam_startZ,
                      reco_beam_endX, reco_beam_endY, reco_beam_endZ;
  std::vector<int> reco_beam_ID, reco_beam_nDaughters;
  int true_beam_nPi0_daughter, true_beam_nPiPlus_daughter,
      true_beam_nPiMinus_daughter, true_beam_nProton_daughter,
      true_beam_nNeutron_daughter;
  std::string true_beam_endProcess;
  int true_beam_nElasticScatters;
  std::vector<double> g4rw_primary_plus_sigma_weight;
  std::vector<double> g4rw_primary_minus_sigma_weight;
  std::vector<double> g4rw_primary_weights;
  std::vector<std::string> g4rw_primary_var;
  std::vector<double> g4rw_alt_primary_plus_sigma_weight;
  std::vector<double> g4rw_alt_primary_minus_sigma_weight;
  std::vector<double> g4rw_test_primary_plus_sigma_weight;
  std::vector<double> g4rw_test_primary_minus_sigma_weight;
  /*std::vector<double> g4rw_full_primary_plus_sigma_weight;
  std::vector<double> g4rw_full_primary_minus_sigma_weight;*/
  double g4rw_primary_singular_weight;
  std::vector<double> g4rw_set_weights;

  
  std::string fGeneratorTag, fPFParticleTag, fTrackerTag;
  //Geant4Reweight stuff
  int RW_PDG;
  TFile FracsFile;
  TFile ProtFracsFile;
  std::vector<fhicl::ParameterSet> ParSet;
  G4ReweightParameterMaker ParMaker;
  G4ReweightManager RWManager;
  G4MultiReweighter MultiRW;
  //G4MultiReweighter ProtMultiRW;
  //G4ReweighterFactory RWFactory;
  //G4Reweighter * theRW;

};


protoana::G4RWExampleAnalyzer::G4RWExampleAnalyzer(
    fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
      fGeneratorTag(p.get<std::string>("GeneratorTag")),
      fPFParticleTag(p.get<std::string>("PFParticleTag")),
      fTrackerTag(p.get<std::string>("TrackerTag")),
      RW_PDG(p.get<int>("RW_PDG")),
      FracsFile( (p.get< std::string >( "FracsFile" )).c_str(), "OPEN" ),
      //ProtFracsFile( (p.get< std::string >( "ProtFracsFile" )).c_str(), "OPEN" ),
      ParSet(p.get<std::vector<fhicl::ParameterSet>>("ParameterSet")),
      ParMaker(ParSet, RW_PDG),
      RWManager({p.get<fhicl::ParameterSet>("Material")}),
      MultiRW(RW_PDG, FracsFile, ParSet,
              p.get<fhicl::ParameterSet>("Material"),
              &RWManager) {//,
      //ProtMultiRW(2212, ProtXSecFile, ProtFracsFile, ParSet) {

  //theRW = RWFactory.BuildReweighter(RW_PDG, &XSecFile, &FracsFile,
  //                                  ParMaker.GetFSHists(),
  //                                  ParMaker.GetElasticHist()/*, true*/ );
  std::cout << "done" << std::endl;
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
  true_beam_startX = true_beam_particle->Position(0).X();
  true_beam_startY = true_beam_particle->Position(0).Y();
  true_beam_startZ = true_beam_particle->Position(0).Z();
  true_beam_endX = true_beam_particle->EndX();
  true_beam_endY = true_beam_particle->EndY();
  true_beam_endZ = true_beam_particle->EndZ();
  true_beam_startP = true_beam_particle->P();
  true_beam_endP = true_beam_particle->P(
      (true_beam_particle->NumberTrajectoryPoints() - 2));
  true_beam_endProcess = true_beam_particle->EndProcess();

  for (int i = 0; i < true_beam_particle->NumberDaughters(); ++i) {
    int daughterID = true_beam_particle->Daughter(i);
    auto part = plist[ daughterID ];
    int pdg = part->PdgCode();

    if (part->Process().find("Inelastic") != std::string::npos) {
      if (pdg == 211) {
        ++true_beam_nPiPlus_daughter;
      }
      else if (pdg == -211) {
        ++true_beam_nPiMinus_daughter;
      }
      else if (pdg == 111) {
        ++true_beam_nPi0_daughter;
      }
      else if (pdg == 2212) {
        ++true_beam_nProton_daughter;
      }
      else if (pdg == 2112) {
        ++true_beam_nNeutron_daughter;
      }
    }

  }

  std::cout << "bting" << std::endl;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  auto view2_IDEs = bt_serv->TrackIdToSimIDEs_Ps(true_beam_ID, geo::View_t(2));
  if (view2_IDEs.size() > 1) {
    std::sort(view2_IDEs.begin(), view2_IDEs.end(),
              [](const sim::IDE * i1, const sim::IDE * i2)
                  {return (i1->z < i2->z);});
    for (size_t i = 0; i < view2_IDEs.size()-1; ++i) {
      const sim::IDE * i1 = view2_IDEs[i]; 
      const sim::IDE * i2 = view2_IDEs[i+1]; 

      double edep = i1->energy;
      double len = sqrt(std::pow((i1->x - i2->x), 2) + 
                        std::pow((i1->y - i2->y), 2) + 
                        std::pow((i1->z - i2->z), 2));
      if (len > 1.e-5)
        true_beam_dEdX.push_back((edep/len));
    }
  }
  std::cout << "done bt" << std::endl;

  try {
    protoana::ProtoDUNEPFParticleUtils pfpUtil;
    auto pfpVec =
        e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
    for (const recob::PFParticle & pfp : (*pfpVec)) {
      const recob::Track * thisTrack =
          pfpUtil.GetPFParticleTrack(pfp, e, fPFParticleTag, fTrackerTag);

      reco_beam_nDaughters.push_back(pfp.NumDaughters());
      if (thisTrack) {
        reco_beam_len.push_back(thisTrack->Length());
        reco_beam_ID.push_back(thisTrack->ID());
        reco_beam_startX.push_back(thisTrack->Start().X());
        reco_beam_startY.push_back(thisTrack->Start().Y());
        reco_beam_startZ.push_back(thisTrack->Start().Z());
        reco_beam_endX.push_back(thisTrack->End().X());
        reco_beam_endY.push_back(thisTrack->End().Y());
        reco_beam_endZ.push_back(thisTrack->End().Z());
      }
      else {
        reco_beam_len.push_back(-1.);
        reco_beam_ID.push_back(-1);
        reco_beam_startX.push_back(-1.);
        reco_beam_startY.push_back(-1.);
        reco_beam_startZ.push_back(-1.);
        reco_beam_endX.push_back(-1.);
        reco_beam_endY.push_back(-1.);
        reco_beam_endZ.push_back(-1.);
      }

    }
  }
  catch (const std::exception & e) {
    std::cout << "No reco info. Moving on" << std::endl;
  }

  const simb::MCTrajectory & true_beam_trajectory =
      true_beam_particle->Trajectory();
  true_beam_nTrajPts = true_beam_trajectory.size();
  for (int i = 1; i < true_beam_nTrajPts; ++i) {
    double x0 = true_beam_trajectory.X(i-1);
    double x1 = true_beam_trajectory.X(i);
    double y0 = true_beam_trajectory.Y(i-1);
    double y1 = true_beam_trajectory.Y(i);
    double z0 = true_beam_trajectory.Z(i-1);
    double z1 = true_beam_trajectory.Z(i);

    true_beam_traj_dx.push_back(sqrt((x1 - x0)*(x1 - x0) +
                                     (y1 - y0)*(y1 - y0) +
                                     (z1 - z0)*(z1 - z0)));
  }
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

  if (true_beam_PDG == RW_PDG) {
  std::cout << "Doing reweight" << std::endl;
    G4ReweightTraj theTraj(true_beam_ID, true_beam_PDG, 0, event, {0,0});
    /*
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

      //For testing with Heng-Ye's parameters
      if (ParSet.size() == 2) {
        for (size_t i = 0; i < 20; ++i) {
          for (size_t j = 0; j < 20; ++j) {
            std::vector<double> input_values = {(.1 + i*.1), (.1 + j*.1)};
            bool set_values = MultiRW.SetAllParameterValues(input_values);
            if (!set_values) continue;

            g4rw_set_weights.push_back(
                MultiRW.GetWeightFromSetParameters(theTraj));
          }
        }
      }


    }
    */

    std::vector<G4ReweightTraj *> trajs = CreateNRWTrajs(
        *true_beam_particle, plist,
        fGeometryService, event, true);
    
    bool added = false;
    for (size_t i = 0; i < trajs.size(); ++i) {
      if (trajs[i]->GetNSteps() > 0) {
        //std::cout << i << " " << trajs[i]->GetNSteps() << std::endl;
        for (size_t j = 0; j < ParSet.size(); ++j) {
          std::pair<double, double> pm_weights =
              MultiRW.GetPlusMinusSigmaParWeight((*trajs[i]), j);
          //std::cout << "got weights" << std::endl;
          //std::cout << pm_weights.first << " " << pm_weights.second << std::endl;

          if (!added) {
            g4rw_alt_primary_plus_sigma_weight.push_back(pm_weights.first);
            g4rw_alt_primary_minus_sigma_weight.push_back(pm_weights.second);
          }
          else {
            g4rw_alt_primary_plus_sigma_weight[j] *= pm_weights.first;
            g4rw_alt_primary_minus_sigma_weight[j] *= pm_weights.second;
          }
        }
        added = true;
      }
    }
    
    for (size_t i = 0; i < ParSet.size(); ++i) {
      std::pair<double, double> temp_weights = GetNTrajPMSigmaWeights(trajs, MultiRW, i);
      g4rw_test_primary_plus_sigma_weight.push_back(temp_weights.first);
      g4rw_test_primary_minus_sigma_weight.push_back(temp_weights.second);
    }

    std::cout << "Testing setting pars" << std::endl;
    std::vector<double> input(ParSet.size(), 1.);
    bool set_values = MultiRW.SetAllParameterValues(input);
    if (set_values) {
      std::cout << GetNTrajWeightFromSetPars(trajs, MultiRW);
    }

    //Heng-Ye's weights
    //for (ii reac loop)
    //for (jj elast loop)
    //  bool set_values = MultiRW.SetAllParameterValues(input_values);
    //  if (!set_values) continue;
    //  double temp_w = GetNTrajWeightFromSetPars(trajs, MultiRW);
    //  g4rw_set_weights.push_back(temp_w);

    /*
    std::vector<int> to_create = {true_beam_ID};
    std::vector<std::vector<G4ReweightTraj *>> created;
    while (to_create.size()) {
      auto part = plist[to_create[0]];
      std::vector<G4ReweightTraj *> temp_trajs =
          CreateNRWTrajs(*part, plist, fGeometryService,
                         event, true);
      std::cout << "size: " << temp_trajs.size() << std::endl;
      if (temp_trajs.size()) {
        auto last_traj = temp_trajs.back();
        std::cout << "created " << last_traj->GetTrackID() << " " <<
                     last_traj->GetPDG() << std::endl;
        for (size_t i = 0; i < last_traj->GetNChilds(); ++i) {
          if ((last_traj->GetChild(i)->GetPDG() == 2212) ||
              (last_traj->GetChild(i)->GetPDG() == 2112) ||
              (abs(last_traj->GetChild(i)->GetPDG()) == 211) ) {
            to_create.push_back(last_traj->GetChild(i)->GetTrackID());
            std::cout << "Adding daughter " << to_create.back() << std::endl;
          }
        }

        if (temp_trajs[0]->GetPDG() == RW_PDG) {
          created.push_back(temp_trajs);
        }
      }
      to_create.erase(to_create.begin());
    }

    std::cout << "Created " << created.size() << " reweightable pi+" << std::endl;

    bool new_added = false;
    for (size_t i = 0; i < created.size(); ++i) {
      std::vector<G4ReweightTraj *> temp_trajs = created[i];
      std::cout << i << " n trajs: " << temp_trajs.size() << std::endl;
      for (size_t j = 0; j < temp_trajs.size(); ++j) {
        G4ReweightTraj * this_traj = temp_trajs[j];
        if (this_traj->GetNSteps() > 0) {
          for (size_t k = 0; k < ParSet.size(); ++k) {
            std::pair<double, double> pm_weights =
                MultiRW.GetPlusMinusSigmaParWeight((*this_traj), k);

            if (!new_added) {
              g4rw_full_primary_plus_sigma_weight.push_back(pm_weights.first);
              g4rw_full_primary_minus_sigma_weight.push_back(pm_weights.second);
            }
            else {
              g4rw_full_primary_plus_sigma_weight[k] *= pm_weights.first;
              g4rw_full_primary_minus_sigma_weight[k] *= pm_weights.second;
            }
          }
          new_added = true;
        }
      }
    }*/
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
  fTree->Branch("true_beam_nTrajPts", &true_beam_nTrajPts);
  fTree->Branch("true_beam_traj_dx", &true_beam_traj_dx);
  fTree->Branch("true_beam_startP", &true_beam_startP);
  fTree->Branch("true_beam_startX", &true_beam_startX);
  fTree->Branch("true_beam_startY", &true_beam_startY);
  fTree->Branch("true_beam_startZ", &true_beam_startZ);
  fTree->Branch("true_beam_endX", &true_beam_endX);
  fTree->Branch("true_beam_endY", &true_beam_endY);
  fTree->Branch("true_beam_endZ", &true_beam_endZ);
  fTree->Branch("true_beam_endP", &true_beam_endP);
  fTree->Branch("true_beam_dEdX", &true_beam_dEdX);
  fTree->Branch("true_beam_nPi0_daughter", &true_beam_nPi0_daughter);
  fTree->Branch("true_beam_nPiPlus_daughter", &true_beam_nPiPlus_daughter);
  fTree->Branch("true_beam_nPiMinus_daughter", &true_beam_nPiMinus_daughter);
  fTree->Branch("true_beam_nProton_daughter", &true_beam_nProton_daughter);
  fTree->Branch("true_beam_nNeutron_daughter", &true_beam_nNeutron_daughter);
  fTree->Branch("true_beam_endProcess", &true_beam_endProcess);
  fTree->Branch("true_beam_nElasticScatters", &true_beam_nElasticScatters);

  fTree->Branch("reco_beam_len", &reco_beam_len);
  fTree->Branch("reco_beam_startX", &reco_beam_startX);
  fTree->Branch("reco_beam_startY", &reco_beam_startY);
  fTree->Branch("reco_beam_startZ", &reco_beam_startZ);
  fTree->Branch("reco_beam_endX", &reco_beam_endX);
  fTree->Branch("reco_beam_endY", &reco_beam_endY);
  fTree->Branch("reco_beam_endZ", &reco_beam_endZ);
  fTree->Branch("reco_beam_ID", &reco_beam_ID);
  fTree->Branch("reco_beam_nDaughters", &reco_beam_nDaughters);

  fTree->Branch("g4rw_primary_weights", &g4rw_primary_weights);
  fTree->Branch("g4rw_primary_singular_weight", &g4rw_primary_singular_weight);
  fTree->Branch("g4rw_primary_plus_sigma_weight", &g4rw_primary_plus_sigma_weight);
  fTree->Branch("g4rw_primary_minus_sigma_weight", &g4rw_primary_minus_sigma_weight);
  fTree->Branch("g4rw_primary_var", &g4rw_primary_var);
  fTree->Branch("g4rw_alt_primary_plus_sigma_weight",
                &g4rw_alt_primary_plus_sigma_weight);
  fTree->Branch("g4rw_alt_primary_minus_sigma_weight",
                &g4rw_alt_primary_minus_sigma_weight);
  fTree->Branch("g4rw_test_primary_plus_sigma_weight",
                &g4rw_test_primary_plus_sigma_weight);
  fTree->Branch("g4rw_test_primary_minus_sigma_weight",
                &g4rw_test_primary_minus_sigma_weight);
  /*
  fTree->Branch("g4rw_full_primary_plus_sigma_weight",
                &g4rw_full_primary_plus_sigma_weight);
  fTree->Branch("g4rw_full_primary_minus_sigma_weight",
                &g4rw_full_primary_minus_sigma_weight);*/
  fTree->Branch("g4rw_set_weights", &g4rw_set_weights);
}

void protoana::G4RWExampleAnalyzer::reset() {
  true_beam_PDG = -1;
  true_beam_ID = -1;
  true_beam_len = -1.;
  true_beam_nTrajPts = -1;
  true_beam_traj_dx.clear();
  true_beam_startP = -1.;
  true_beam_startX = -1.;
  true_beam_startY = -1.;
  true_beam_startZ = -1.;
  true_beam_endX = -1.;
  true_beam_endY = -1.;
  true_beam_endZ = -1.;
  true_beam_endP = -1.;
  true_beam_endProcess = "";
  true_beam_nElasticScatters = 0.;
  true_beam_nPi0_daughter = 0;
  true_beam_nPiPlus_daughter = 0;
  true_beam_nPiMinus_daughter = 0;
  true_beam_nProton_daughter = 0;
  true_beam_nNeutron_daughter = 0;
  true_beam_dEdX.clear();
  
  reco_beam_len.clear();
  reco_beam_ID.clear();
  reco_beam_nDaughters.clear();
  reco_beam_startX.clear();
  reco_beam_startY.clear();
  reco_beam_startZ.clear();
  reco_beam_endX.clear();
  reco_beam_endY.clear();
  reco_beam_endZ.clear();

  g4rw_primary_weights.clear();
  g4rw_primary_singular_weight = 1.;
  g4rw_primary_plus_sigma_weight.clear();
  g4rw_primary_minus_sigma_weight.clear();
  g4rw_primary_var.clear();
  g4rw_alt_primary_plus_sigma_weight.clear();
  g4rw_alt_primary_minus_sigma_weight.clear();
  g4rw_test_primary_plus_sigma_weight.clear();
  g4rw_test_primary_minus_sigma_weight.clear();
  /*g4rw_full_primary_plus_sigma_weight.clear();
  g4rw_full_primary_minus_sigma_weight.clear();*/
  g4rw_set_weights.clear();
}
DEFINE_ART_MODULE(protoana::G4RWExampleAnalyzer)
