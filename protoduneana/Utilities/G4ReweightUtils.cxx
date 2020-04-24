#include "G4ReweightUtils.h"
#include "geant4reweight/src/ReweightBase/G4ReweightStep.hh"

bool protoana::G4ReweightUtils::CreateRWTraj(
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
