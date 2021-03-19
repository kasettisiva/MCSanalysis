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
    reco_beam_incidentEnergies = std::vector<double>();
    true_beam_incidentEnergies = std::vector<double>();
    true_beam_traj_Z = std::vector<double>();
    true_beam_traj_KE = std::vector<double>();
    true_beam_slices = std::vector<int>();
    calibrated_dQdX = std::vector<double>();
    beam_EField = std::vector<double>();
    track_pitch = std::vector<double>();
  };

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
 private:
  int event_ID, subrun_ID, run_ID;
  int sample_ID;
  int selection_ID;
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  double true_beam_endP, true_beam_mass;
  double reco_beam_endZ, true_beam_startP;
  double beam_inst_P;
  std::vector<double> reco_beam_incidentEnergies,
                      true_beam_incidentEnergies,
                      true_beam_traj_Z,
                      true_beam_traj_KE;
  std::vector<int> true_beam_slices;
  std::vector<double> calibrated_dQdX, beam_EField,
                      track_pitch;

};
}
#endif
