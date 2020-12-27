#ifndef THINSLICEDATASET_hh
#define THINSLICEDATASET_hh

#include <vector>
#include <map>

#include "fhiclcpp/ParameterSet.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

namespace protoana {
class ThinSliceDataSet {
 public:
  ThinSliceDataSet() {};
  ThinSliceDataSet(const std::vector<double> & incident_bins,
                   const std::vector<fhicl::ParameterSet> & selections);
  ~ThinSliceDataSet() {};

  const std::map<int, std::string> & GetSelectionNames() const {
    return fSelectionNames;
  };

  std::string & GetSelectionName(int id) {
    return fSelectionNames.at(id);
  }

  std::map<int, TH1 *> & GetSelectionHists() {
    return fSelectionHists;
  };

  std::map<int, TH1 *> & GetRebinnedSelectionHists() {
    return fSelectionHistsRebinned;
  };

  TH1 * GetSelectionHist(int id) {
    return fSelectionHists.at(id);
  };

  TH1D & GetIncidentHist() {
    return fIncidentHist;
  };

  TH1D & GetRebinnedIncidentHist() {
    return fIncidentHistRebinned;
  };

  TH1 * GetRebinnedSelectionHist(int id) {
    return fSelectionHistsRebinned.at(id);
  };

  void FillIncidentHist(const std::vector<double> & vals) {
    for (size_t i = 0; i < vals.size(); ++i) {
      fIncidentHist.Fill(vals.at(i));
    }
  };

  void FillSelectionHist(int id, double val) {
    if (fSelectionHists.find(id) != fSelectionHists.end()) {
      fSelectionHists.at(id)->Fill(val);
    }
  };

  template <size_t N> void FillSelectionHist(int id, const double (& vals)[N]) {
    if (N < 1 || N > 3) {
      std::string message = "Error: trying to fill hists with too many values";
      message += std::to_string(N);
      throw std::runtime_error(message);
    }

    if (fSelectionHists.find(id) != fSelectionHists.end()) {
      if (N == 1) {
        fSelectionHists.at(id)->Fill(vals[0]);
      }
      else if (N == 2) {
        ((TH2D*)fSelectionHists.at(id))->Fill(vals[0], vals[1]);
      }
      else if (N == 3) {
        ((TH3D*)fSelectionHists.at(id))->Fill(vals[0], vals[1], vals[2]);
      }
    }
  }

  void MakeRebinnedHists();

 private:
  void Rebin1D(TH1 * sel_hist, TH1 * rebinned);
  void Rebin2D(TH1 * sel_hist, TH1 * rebinned);
  void Rebin3D(TH1 * sel_hist, TH1 * rebinned);
  std::map<int, TH1 *> fSelectionHists;
  TH1D fIncidentHist;
  std::map<int, TH1 *> fSelectionHistsRebinned;
  TH1D fIncidentHistRebinned;
  bool fMadeRebinned = false;
  std::map<int, std::string> fSelectionNames;

};
}
#endif
