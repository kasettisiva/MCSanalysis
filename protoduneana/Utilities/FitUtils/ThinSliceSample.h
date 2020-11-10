#include "TH1D.h"

#include <map>
#include <sstream>

namespace protoana {
std::string PreciseToString(const double val, const int n = 2);/*{
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << val;
  return out.str();
};*/


#ifndef THINSLICESAMPLE_hh
#define THINSLICESAMPLE_hh
class ThinSliceSample {
 public:
  ThinSliceSample(std::string name, int flux_type,
                  const std::map<int, std::string> & selections,
                  const std::vector<double> & incident_bins,
                  const std::vector<double> & selected_bins,
                  bool is_signal = false, std::pair<double, double> range = {0., 0.});
  ~ThinSliceSample(){};

  void SetFactor(double f) {fFactor = f;};

  const std::map<int, TH1D> & GetSelectionHists() const {
    return fSelectionHists;
  };

  TH1D & GetSelectionHist(int id) {
    return fSelectionHists.at(id);
  };

  TH1D & GetIncidentHist() {
    return fIncidentHist;
  };

  TH1D & GetRebinnedIncidentHist() {
    return fIncidentHistRebinned;
  };

  TH1D & GetRebinnedSelectionHist(int id) {
    return fSelectionHistsRebinned.at(id);
  };

  const std::string & GetName() const {
    return fSampleName;
  };

  const int & GetFluxType() const {
    return fFluxType;
  };

  const double & GetNominalFlux() const {
    return fNominalFlux;
  };

  void AddFlux(double val = 1.) {
    fNominalFlux += val;
  };

  void FillIncidentHist(const std::vector<double> & vals) {
    for (size_t i = 0; i < vals.size(); ++i) {
      fIncidentHist.Fill(vals.at(i));
    }
  };

  void FillSelectionHist(int id, double val) {
    if (fSelectionHists.find(id) != fSelectionHists.end()) {
      fSelectionHists.at(id).Fill(val);
    }
  };

  void ScaleHists(double val) {
    fIncidentHist.Scale(val);
    for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
      it->second.Scale(val);
    }
  };

  void SetDataMCScale(double val) {
    fDataMCScale = val;
    ScaleHists(fDataMCScale);
  };

  void SetFactorAndScale(double val) {
    ResetFactor();
    fFactor = val;
    ScaleHists(val);
  };

  void ResetFactor() {
    ScaleHists(1./fFactor);
    fFactor = 1.;
  };

  bool CheckIsSignal() {return fIsSignal;};
  bool CheckInSignalRange(double val) {return ((fRange.first < val) &&
                                               (val <= fRange.second));};
  const std::pair<double, double> & GetRange() const {return fRange;};

  void MakeRebinnedHists();
  void RefillRebinnedHists();

 private:
  double fFactor = 1.;
  std::map<int, TH1D> fSelectionHists;
  TH1D fIncidentHist;
  std::map<int, TH1D> fSelectionHistsRebinned;
  TH1D fIncidentHistRebinned;
  std::string fSampleName;
  int fFluxType;
  double fNominalFlux = 0.;
  double fDataMCScale = 1.;
  bool fIsSignal;
  std::pair<double, double> fRange;
  bool fMadeRebinned = false;
};

}
#endif
