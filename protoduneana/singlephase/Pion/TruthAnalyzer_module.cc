////////////////////////////////////////////////////////////////////////
// Class:       TruthAnalyzer
// Plugin Type: analyzer (art v3_05_01)
// File:        TruthAnalyzer_module.cc
//
// Generated at Tue Nov 24 17:23:39 2020 by Jacob Calcutt using cetskelgen
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
#include "larcore/Geometry/Geometry.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace pionana {
  class TruthAnalyzer;
}


class pionana::TruthAnalyzer : public art::EDAnalyzer {
public:
  explicit TruthAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TruthAnalyzer(TruthAnalyzer const&) = delete;
  TruthAnalyzer(TruthAnalyzer&&) = delete;
  TruthAnalyzer& operator=(TruthAnalyzer const&) = delete;
  TruthAnalyzer& operator=(TruthAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  art::InputTag fGeneratorTag;

};


pionana::TruthAnalyzer::TruthAnalyzer(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
      fGeneratorTag(p.get<art::InputTag>("GeneratorTag")) {
}

void pionana::TruthAnalyzer::analyze(art::Event const& e)
{
  protoana::ProtoDUNETruthUtils truthUtil;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  auto mcTruths = e.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
  art::ServiceHandle<geo::Geometry> geom;
  const simb::MCParticle* true_beam_particle
      = truthUtil.GetGeantGoodParticle((*mcTruths)[0],e);
  const simb::MCTrajectory & true_beam_trajectory = true_beam_particle->Trajectory();


  double init_KE = 0.;
  //Get the mass for the beam particle
  double mass = 139.57;
  //int true_beam_ID = true_beam_particle->TrackId();
  int true_beam_PDG = true_beam_particle->PdgCode();
  if( true_beam_PDG == 2212 ) mass = 938.27;
  else if( abs(true_beam_PDG) == 211 ) mass = 139.57;
  else if( abs(true_beam_PDG) == 11 ) mass = .511;
  else if( abs(true_beam_PDG) == 321 ) mass = 321;
  else if( abs(true_beam_PDG) == 13 )  mass = 105.66;

  //auto view2_IDEs = bt_serv->TrackIdToSimIDEs_Ps( true_beam_ID, geo::View_t(2) );
  //std::sort(view2_IDEs.begin(), view2_IDEs.end(),
  //          [](const sim::IDE * i1, const sim::IDE * i2) {return (i1->z < i2->z);});
  //if (view2_IDEs.size()) {
  //  double ide_z = view2_IDEs[0]->z;

  auto true_beam_procs = true_beam_trajectory.TrajectoryProcesses();
  std::map<size_t, unsigned char> proc_map(true_beam_procs.begin(),
                                           true_beam_procs.end());


  for (size_t i = 0; i < true_beam_trajectory.size(); ++i) {
    double z0 = true_beam_trajectory.Z(i);
    //double z1 = true_beam_trajectory.Z(i);

    double x0 = true_beam_trajectory.X(i);
    double y0 = true_beam_trajectory.Y(i);
    //double x1 = true_beam_trajectory.X(i);
    //double y1 = true_beam_trajectory.Y(i);
    geo::Point_t test_point{x0, y0, z0};
    //geo::Point_t test_point_1{x1, y1, z1};
    const TGeoMaterial * mat = geom->Material(test_point);
    std::cout << "Point " << i << " (" << x0 << ", " << y0 <<
                 ", " << z0 << ") Material " << mat->GetName() << " " <<
                 (proc_map.find(i) != proc_map.end() ?
                  true_beam_trajectory.KeyToProcess(proc_map[i]) :
                  "") <<
                 std::endl;
    //if (z0 < ide_z && z1 > ide_z) {
      //const TGeoMaterial * mat1 = geom->Material(test_point_1);
      //std::cout << "Point " << i << " (" << x1 << ", " << y1 <<
      //             ", " << z1 << ") Material " << mat1->GetName() <<
      //             std::endl;
      init_KE = 1.e3 * true_beam_trajectory.E(i) - mass;
      //std::cout << "Found matching position" << z0 << " " << ide_z <<
      //             " " << z1 << std::endl;
      std::cout << "init KE: " << init_KE << " " <<
                   (1.e3*true_beam_trajectory.E(i) - mass) << std::endl;
    //  break;
    //}
  }
  //}
}

DEFINE_ART_MODULE(pionana::TruthAnalyzer)
