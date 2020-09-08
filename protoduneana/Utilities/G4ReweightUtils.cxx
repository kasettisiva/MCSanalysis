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

std::vector<G4ReweightTraj *> protoana::G4ReweightUtils::CreateNRWTrajs(
    const simb::MCParticle & part, const sim::ParticleList & plist,
    art::ServiceHandle < geo::Geometry > geo_serv, int event, bool fVerbose) {
  std::vector<G4ReweightTraj *> results;


  //Create process map
  auto procs = part.Trajectory().TrajectoryProcesses();
  std::map<size_t, std::string> proc_map;
  for (auto it = procs.begin(); it != procs.end(); ++it) {
    proc_map[it->first] = part.Trajectory().KeyToProcess(it->second);
  }

  std::vector<double> traj_X, traj_Y, traj_Z;
  //std::vector<double> traj_PX, traj_PY, traj_PZ;
  //std::vector<size_t> elastic_indices;

  std::vector<std::pair<size_t, size_t>> ranges;

  //bool found_last = false;
  bool found_LAr = false;
  size_t start = 0, end = 0;
  //G4ReweightTraj theTraj(part.TrackId(), part.PdgCode(), 0, event, {0,0});
  for (size_t i = 0; i < part.NumberTrajectoryPoints(); ++i) {
    double x = part.Position(i).X();
    double y = part.Position(i).Y();
    double z = part.Position(i).Z();

    geo::Point_t test_point{x, y, z};
    const TGeoMaterial * test_material = geo_serv->Material(test_point);

    if (!strcmp(test_material->GetName(), "LAr")) {
      if (fVerbose) {
        std::cout << i << " " << "LAr: " << test_material->GetDensity() << " " <<
                     test_material->GetA() << " " << test_material->GetZ() <<
                     " " << x << " " << y << " " << z << 
                     std::endl;
      }

      if (!found_LAr) {
        found_LAr = true;
        start = i;
      }

      //traj_PX.push_back(part.Px(i));
      //traj_PY.push_back(part.Py(i));
      //traj_PZ.push_back(part.Pz(i));

      //auto itProc = proc_map.find(i);
      //if (itProc != proc_map.end() && itProc->second == "hadElastic") {
      //  elastic_indices.push_back(i);
      //}
    }
    else {
      if (fVerbose) {
        std::cout << i << " " << test_material->GetName() << " " <<
                     test_material->GetDensity() << " " <<
                     test_material->GetA() << " " << test_material->GetZ() <<
                     " " << x << " " << y << " " << z << 
                     std::endl;
      }
      if (found_LAr) {
        found_LAr = false;
        end = i;
        ranges.push_back({start, end});
      }
    }

    //if (i == part.NumberTrajectoryPoints() - 1)
    //  found_last = true;
  }
  if (found_LAr) {
    //size_t np = part.NumberTrajectoryPoints();
    ranges.push_back({start, part.NumberTrajectoryPoints() - 1});
    //double x = part.Position(np - 1).X();
    //double y = part.Position(np - 1).Y();
    //double z = part.Position(np - 1).Z();
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
      return results;
      break;
    }
  }

  for (size_t i = 0; i < ranges.size(); ++i) {
    //std::cout << ranges[i].first << " " << ranges[i].second << std::endl;
    G4ReweightTraj * theTraj = new G4ReweightTraj(i, part.PdgCode(), 0, event, {0,0});
    
    for (size_t j = ranges[i].first; j < ranges[i].second; ++j) {
      double dx = part.Position(j+1).X() - part.Position(j).X();
      double dy = part.Position(j+1).Y() - part.Position(j).Y();
      double dz = part.Position(j+1).Z() - part.Position(j).Z();

      double len = sqrt(dx*dx + dy*dy + dz*dz);
  
      double preStepP[3] = {part.Px(j)*1.e3,
                            part.Py(j)*1.e3,
                            part.Pz(j)*1.e3};
  
      double postStepP[3] = {part.Px(j + 1)*1.e3,
                             part.Py(j + 1)*1.e3,
                             part.Pz(j + 1)*1.e3};
      if (j == ranges[i].first) {
        double p_squared = preStepP[0]*preStepP[0] + preStepP[1]*preStepP[1] +
                           preStepP[2]*preStepP[2];
        theTraj->SetEnergy(sqrt(p_squared + mass*mass));
      }

      //prochere
      auto itProc = proc_map.find(j);
      std::string proc = "default";
      if (itProc != proc_map.end() &&
          j != (part.NumberTrajectoryPoints() - 2)) {
        proc = itProc->second;
      }
      //- 2 because the last element is the end of the last step
      else if (j == (part.NumberTrajectoryPoints() - 2)) {
        proc = part.EndProcess();
      }
      //std::cout << j << " Proc: " << proc << std::endl;
      G4ReweightStep * step = new G4ReweightStep(i, part.PdgCode(),
                                                 0, event, preStepP, postStepP,
                                                 len, proc);
      theTraj->AddStep(step);
    }

    results.push_back(theTraj);
  }

  if (results.size()) {
    //Loop over daughters
    for (int i = 0; i < part.NumberDaughters(); ++i) {
      int d_index = part.Daughter(i);
      auto d_part = plist[d_index];

      int d_PDG = d_part->PdgCode();
      int d_ID = d_part->TrackId();

      results.back()->AddChild(new G4ReweightTraj(d_ID, d_PDG,
                               results.size() - 1, event, {0,0}));
    }
  }
  return results;
  /*

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
  */
}
