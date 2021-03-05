#ifndef THINSLICESAMPLE_hh
#define THINSLICESAMPLE_hh

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include <map>
#include <sstream>
#include <stdexcept>

#include "fhiclcpp/ParameterSet.h"

namespace protoana {
std::string PreciseToString(const double val, const int n = 2);


//Need to split up samples by its true incident energy (from the beam line)

class ThinSliceSample {
 public:
  ThinSliceSample(std::string name, int flux_type,
                  const std::vector<fhicl::ParameterSet> & selections,
                  const std::vector<double> & incident_bins,
                  const std::vector<double> & true_incident_bins,
                  size_t beam_energy_bin,
                  bool is_signal = false, std::pair<double, double> range = {0., 0.});

  ~ThinSliceSample(){};

  void SetFactor(double f) {fFactor = f;};

  const std::map<int, TH1 *> & GetSelectionHists() const {
    return fSelectionHists;
  };

  const std::map<int, TH1 *> & GetRebinnedSelectionHists() const {
    return fSelectionHistsRebinned;
  };

/*
  const std::map<int, TH1 *> & GetBestFitSelectionHists() const {
    SetFactorAndScale(fBestFitFactor);
    return fSelectionHists;
    //return fBestFitSelectionHists;
  };

  const std::map<int, TH1 *> & GetBestFitRebinnedSelectionHists() const {
    SetFactorAndScale(fBestFitFactor);
    return fSelectionHistsRebinned;
    //return fBestFitSelectionHistsRebinned;
  };*/

  TH1 * GetSelectionHist(int id) {
    return fSelectionHists.at(id);
  };

  TH1D & GetIncidentHist() {
    return fIncidentHist;
  };

  TH1D & GetTrueIncidentHist() {
    return fTrueIncidentHist;
  };

  TH1D & GetRebinnedIncidentHist() {
    return fIncidentHistRebinned;
  };

  TH1 * GetRebinnedSelectionHist(int id) {
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

  const double & GetVariedFlux() const {
    return fVariedFlux;
  };

  void AddFlux(double val = 1.) {
    fNominalFlux += val;
    fVariedFlux += val;
  };

  void FillIncidentHist(const std::vector<double> & vals) {
    for (size_t i = 0; i < vals.size(); ++i) {
      fIncidentHist.Fill(vals.at(i));
    }
  };

  void AddIncidentEnergies(const std::vector<double> & vals) {
    for (auto v : vals)
      fIncidentEnergies.push_back({v, 1.});
  };

  void AddESliceEnergies(const std::pair<double, double> & vals) {
    fESliceEnergies.push_back(vals);
  };

  void FillTrueIncidentHist(const std::vector<double> & vals) {
    for (size_t i = 0; i < vals.size(); ++i) {
      fTrueIncidentHist.Fill(vals.at(i));
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

  void FillHistFromIncidentEnergies(TH1D & hist) {
    for (auto vals : fIncidentEnergies) {
      hist.Fill(vals.first, fFactor/*vals.second*/);
    }
  };

  void FillESliceHist(TH1D & hist) {
    for (auto vals : fESliceEnergies) {
      int last_bin = hist.FindBin(vals.first);
      int first_bin = hist.FindBin(vals.second);
      for (int i = first_bin; i <= last_bin; ++i)
        hist.AddBinContent(i, fFactor);
    }
  };

  void ScaleHists(double val) {
    fIncidentHist.Scale(val);
    for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
      it->second->Scale(val);
    }

    fTrueIncidentHist.Scale(val);
  };

  void ScaleIncidentEnergies(double val) {
    for (auto it = fIncidentEnergies.begin();
         it != fIncidentEnergies.end(); ++it) {
      it->second *= val;
    }
  };

  void SetDataMCScale(double val) {
    fDataMCScale = val;
    ScaleHists(fDataMCScale);
    //ScaleIncidentEnergies(fDataMCScale);
    fNominalFlux *= val;
    fVariedFlux *= val;
  };

  void SetFactorAndScale(double val) {
    ResetFactor();
    fFactor = val;
    //fNominalFlux *= val;
    fVariedFlux *= val;
    ScaleHists(val);
    //ScaleIncidentEnergies(val);
  };

  void ResetFactor() {
    ScaleHists(1./fFactor);
    //ScaleIncidentEnergies(1./fFactor);
    //fNominalFlux *= (1./fFactor);
    fVariedFlux *= (1./fFactor);
    fFactor = 1.;
  };

  void SetFactorToBestFit() {
    SetFactorAndScale(fBestFitFactor);
  };

  double GetBestFitFactor() {
    return fBestFitFactor;
  };

  void SetBestFit() {
    if (fBestFitIsSet) {
      return;
    }

    fBestFitFactor = fFactor;
    fBestFitIsSet = true;
  /*
    for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
      std::string name = it->second->GetName();
      name += "BestFit";
      fBestFitSelectionHists[it->first] = (TH1*)it->second->Clone(name.c_str());

      name = fSelectionHistsRebinned[it->first]->GetName();
      name += "BestFit";
      fBestFitSelectionHistsRebinned[it->first]
          = (TH1*)fSelectionHistsRebinned[it->first]->Clone(name.c_str());
    }*/
  };

  bool CheckIsSignal() {return fIsSignal;};
  bool CheckInSignalRange(double val) {return ((fRange.first < val) &&
                                               (val <= fRange.second));};
  double RangeLowEnd() {return fRange.first;};
  double RangeHighEnd() {return fRange.second;};

  const std::pair<double, double> & GetRange() const {return fRange;};

  void RefillRebinnedHists();
  void MakeRebinnedHists();

 private:
  double fFactor = 1., fBestFitFactor = 1.;
  bool fBestFitIsSet = false;
  std::string fSampleName;
  int fFluxType;
  double fNominalFlux = 0.;
  double fVariedFlux = 0.;
  double fDataMCScale = 1.;
  bool fIsSignal;
  std::pair<double, double> fRange;

  void Rebin1D(TH1 * sel_hist, TH1 * rebinned);
  void Rebin2D(TH1 * sel_hist, TH1 * rebinned);
  void Rebin3D(TH1 * sel_hist, TH1 * rebinned);
  std::map<int, TH1 *> fSelectionHists;
  //std::map<int, TH1 *> fBestFitSelectionHists;
  TH1D fIncidentHist;
  TH1D fTrueIncidentHist;
  std::map<int, TH1 *> fSelectionHistsRebinned;
  //std::map<int, TH1 *> fBestFitSelectionHistsRebinned;
  TH1D fIncidentHistRebinned;
  bool fMadeRebinned = false;

  std::vector<std::pair<double, double>> fIncidentEnergies;
  std::vector<std::pair<double, double>> fESliceEnergies;

};

}
#endif
