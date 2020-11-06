#ifndef THINSLICESAMPLE_hh
#define THINSLICESAMPLE_hh
#include "TH1D.h"

#include <map>

namespace protoana {

class ThinSliceSample {
 public:
  ThinSliceSample(std::string name,
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

  const std::string & GetName() const {
    return fSampleName;
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
  }

  bool CheckIsSignal() {return fIsSignal;};
  bool CheckInSignalRange(double val) {return ((fRange.first < val) &&
                                               (val <= fRange.second));};
  const std::pair<double, double> & GetRange() const {return fRange;};

 private:
  double fFactor = 1.;
  std::map<int, TH1D> fSelectionHists;
  TH1D fIncidentHist;
  std::string fSampleName;
  bool fIsSignal;
  std::pair<double, double> fRange;
};

}
#endif
