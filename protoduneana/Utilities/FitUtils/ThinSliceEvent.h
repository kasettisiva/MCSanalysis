#ifndef THINSLICEEVENT_hh
#define THINSLICEEVENT_hh
namespace protoana {
class ThinSliceEvent {
 public:
  ThinSliceEvent(int event, int subrun, int run)
    : event_ID(event), subrun_ID(subrun), run_ID(run) {
    sample_ID = -999;
    selection_ID = -999;
    true_beam_interactingEnergy = -999;
    reco_beam_interactingEnergy = -999;
    true_beam_endP = -999;
    true_beam_startP = -999;
    true_beam_mass = -999;
    reco_beam_endZ = -999;
    beam_inst_P = -999;
    pdg = -999;
    reco_beam_incidentEnergies = std::vector<double>();
    true_beam_incidentEnergies = std::vector<double>();
    true_beam_traj_Z = std::vector<double>();
    true_beam_traj_KE = std::vector<double>();
    true_beam_slices = std::vector<int>();
    calibrated_dQdX = std::vector<double>();
    beam_EField = std::vector<double>();
    track_pitch = std::vector<double>();
    g4rw_weights = std::map<std::string, std::vector<double>>();
    reco_daughter_track_thetas = std::vector<double>();
    reco_daughter_track_scores = std::vector<double>();
    reco_daughter_track_dQdX = std::vector<std::vector<double>>();
    reco_daughter_track_res_range = std::vector<std::vector<double>>();
    reco_daughter_efield = std::vector<std::vector<double>>();
    has_pi0_shower = false;
  };

  /*
  int GetEventID() const {
    return event_ID;
  };

  int GetSubrunID() const {
    return subrun_ID;
  };

  int GetRunID() const {
    return run_ID;
  };*/

  int GetSampleID() const {
    return sample_ID;
  };
  void SetSampleID(int s) {
    sample_ID = s;
  };

  int GetSelectionID() const {
    return selection_ID;
  };
  void SetSelectionID(int s) {
    selection_ID = s;
  };

  bool GetHasPi0Shower() const {
    return has_pi0_shower;
  };
  void SetHasPi0Shower(bool s) {
    has_pi0_shower = s;
  };

  double GetTrueInteractingEnergy() const {
    return true_beam_interactingEnergy;
  };
  void SetTrueInteractingEnergy(double e) {
    true_beam_interactingEnergy = e;
  };

  double GetRecoInteractingEnergy() const {
    return reco_beam_interactingEnergy;
  };
  void SetRecoInteractingEnergy(double e) {
    reco_beam_interactingEnergy = e;
  };

  double GetTrueEndP() const {
    return true_beam_endP;
  };
  void SetTrueEndP(double p) {
    true_beam_endP = p;
  };

  double GetRecoEndZ() const {
    return reco_beam_endZ;
  };
  void SetRecoEndZ(double p) {
    reco_beam_endZ = p;
  };

  double GetTrueStartP() const {
    return true_beam_startP;
  };
  void SetTrueStartP(double p) {
    true_beam_startP = p;
  };

  double GetTrueMass() const {
    return true_beam_mass;
  };
  void SetTrueMass(double m) {
    true_beam_mass = m;
  };

  const std::vector<double> & GetRecoIncidentEnergies() const {
    return reco_beam_incidentEnergies;
  };
  void SetRecoIncidentEnergies(std::vector<double> v) {
    reco_beam_incidentEnergies = v;
  };

  const std::vector<double> & GetTrueIncidentEnergies() const {
    return true_beam_incidentEnergies;
  };
  void SetTrueIncidentEnergies(std::vector<double> v) {
    true_beam_incidentEnergies = v;
  };

  const std::vector<double> & GetTrueTrajZ() const {
    return true_beam_traj_Z;
  };
  void SetTrueTrajZ(std::vector<double> v) {
    true_beam_traj_Z = v;
  };

  const std::vector<double> & GetTrueTrajKE() const {
    return true_beam_traj_KE;
  };
  void SetTrueTrajKE(std::vector<double> v) {
    true_beam_traj_KE = v;
  };

  const std::vector<double> & GetRecoDaughterTrackThetas() const {
    return reco_daughter_track_thetas;
  };
  void SetRecoDaughterTrackThetas(std::vector<double> v) {
    reco_daughter_track_thetas = v;
  };

  const std::vector<double> & GetRecoDaughterTrackScores() const {
    return reco_daughter_track_scores;
  };
  void SetRecoDaughterTrackScores(std::vector<double> v) {
    reco_daughter_track_scores = v;
  };

  const std::vector<std::vector<double>>
      & GetRecoDaughterTrackResRanges() const {
    return reco_daughter_track_res_range;
  };
  void AddRecoDaughterTrackResRange(std::vector<double> v) {
    reco_daughter_track_res_range.push_back(v);
  };

  const std::vector<std::vector<double>>
      & GetRecoDaughterTrackdQdXs() const {
    return reco_daughter_track_dQdX;
  };
  void AddRecoDaughterTrackdQdX(std::vector<double> v) {
    reco_daughter_track_dQdX.push_back(v);
  };

  const std::vector<std::vector<double>>
      & GetRecoDaughterEFields() const {
    return reco_daughter_efield;
  };
  void AddRecoDaughterEField(std::vector<double> v) {
    reco_daughter_efield.push_back(v);
  };

  const std::vector<int> & GetTrueSlices() const {
    return true_beam_slices;
  };
  void SetTrueSlices(std::vector<int> v) {
    true_beam_slices = v;
  };

  const std::vector<double> & GetdQdXCalibrated() const {
    return calibrated_dQdX;
  };
  void SetdQdXCalibrated(std::vector<double> v) {
    calibrated_dQdX = v;
  };

  const std::vector<double> & GetEField() const {
    return beam_EField;
  };
  void SetEField(std::vector<double> v) {
    beam_EField = v;
  };

  const std::vector<double> & GetTrackPitch() const {
    return track_pitch;
  };
  void SetTrackPitch(std::vector<double> v) {
    track_pitch = v;
  };

  void SetBeamInstP(double p) {
    beam_inst_P = p;
  };
  double GetBeamInstP() const {
    return beam_inst_P;
  };

  void SetPDG(int p) {
    pdg = p;
  };
  int GetPDG() const {
    return pdg;
  };

  void MakeG4RWBranch(const std::string & br, const std::vector<double> & ws) {
    g4rw_weights[br] = ws;
  };
  double GetG4RWWeight(const std::string & br, size_t i) const {
    if (g4rw_weights.at(br).size() == 0) return 1.;
    return g4rw_weights.at(br).at(i); 
  };
  const std::map<std::string, std::vector<double>> & GetG4RWWeightMap() const {
    return g4rw_weights; 
  };
  const std::vector<double> & GetG4RWBranch(const std::string & br) const {
    return g4rw_weights.at(br);
  };
  bool HasG4RWBranch(const std::string & br) const {
    return (g4rw_weights.find(br) != g4rw_weights.end());
  };

  int GetEventID() const {return event_ID;};
  int GetSubrunID() const {return subrun_ID;};
  int GetRunID() const {return run_ID;};

 private:
  int event_ID, subrun_ID, run_ID;
  int sample_ID;
  int selection_ID;
  int pdg;
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  double true_beam_endP, true_beam_mass;
  double reco_beam_endZ, true_beam_startP;
  double beam_inst_P;
  bool has_pi0_shower;
  std::vector<double> reco_beam_incidentEnergies,
                      true_beam_incidentEnergies,
                      true_beam_traj_Z,
                      true_beam_traj_KE,
                      reco_daughter_track_thetas,
                      reco_daughter_track_scores;
  std::vector<std::vector<double>> reco_daughter_track_dQdX,
                                   reco_daughter_track_res_range,
                                   reco_daughter_efield;

  std::vector<int> true_beam_slices;
  std::vector<double> calibrated_dQdX, beam_EField,
                      track_pitch;
  std::map<std::string, std::vector<double>> g4rw_weights;
};
}
#endif
