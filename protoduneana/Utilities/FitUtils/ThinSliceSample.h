#ifndef THINSLICESAMPLE_hh
#define THINSLICESAMPLE_hh

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TSpline.h"
#include "TDirectory.h"
#include "TCanvas.h"

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

  const std::map<int, std::vector<TH1 *>> &
      GetShifts(std::string syst_name) const {
    return fSystematicShifts.at(syst_name); 
  };

  const std::map<int, std::vector<TSpline3 *>> &
      GetSplines(std::string syst_name) const {
    return fSystematicSplines.at(syst_name);
  };

  double GetSplineWeight(std::string syst_name, double par_val,
                         int selection_ID, double val) const {
    int bin = fSelectionHists.at(selection_ID)->FindBin(val);
    return fSystematicSplines.at(syst_name).at(selection_ID).at(bin-1)->Eval(par_val);
  };

  void AddSystematicShift(TH1 * hist, std::string syst_name,
                          int selection_ID) {
    fSystematicShifts[syst_name][selection_ID].push_back(hist);
  };

  void SaveSystematics(std::string syst_name, TDirectory * dir) {
    dir->cd();
    for (auto it = fSystematicShifts[syst_name].begin();
         it != fSystematicShifts[syst_name].end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        it->second[i]->Write();
      }
    } 
    for (auto it = fSystematicSplines[syst_name].begin();
         it != fSystematicSplines[syst_name].end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        TCanvas c(it->second[i]->GetName(), "");
        it->second[i]->Draw();
        c.Write(it->second[i]->GetName());
      }
    }
  };

  const std::map<int, TH1 *> & GetRebinnedSelectionHists() const {
    return fSelectionHistsRebinned;
  };

  TH1 * GetSelectionHist(int id) {
    return fSelectionHists.at(id);
  };

  /*
  TH1D & GetIncidentHist() {
    return fIncidentHist;
  };*/

  TH1D & GetTrueIncidentHist() {
    return fTrueIncidentHist;
  };

/*
  TH1D & GetRebinnedIncidentHist() {
    return fIncidentHistRebinned;
  };*/

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

  void AddVariedFlux(double val = 1.) {
    fVariedFlux += val;
  }

  void FillSystematicShift(std::string syst_name,
                           int selection_ID,
                           const std::vector<double> & vals) {
    if (vals.size() !=
        fSystematicShifts[syst_name][selection_ID].size()) {
      std::string message = "ThinSliceSample: Input systematic shift values and number of shift hists differ"; 
      throw std::runtime_error(message);
    }

    for (size_t i = 0; i < vals.size(); ++i) {
      fSystematicShifts[syst_name][selection_ID][i]->Fill(vals[i]);
    }
  };

  void FillSystematicShift(std::string syst_name,
                           int selection_ID,
                           const std::vector<double> & vals,
                           const std::vector<double> & weights) {
    if (vals.size() !=
        fSystematicShifts[syst_name][selection_ID].size()) {
      std::string message = "ThinSliceSample: Input systematic shift values and number of shift hists differ"; 
      throw std::runtime_error(message);
    }

    for (size_t i = 0; i < vals.size(); ++i) {
      fSystematicShifts[syst_name][selection_ID][i]->Fill(vals[i], weights[i]);
    }
  };

  void SetSystematicVals(std::string syst_name, std::vector<double> & vals) {
    fSystematicVals[syst_name] = vals;
  };

  void MakeSystematicSplines(std::string syst_name/*, bool do_insert = true*/) {
    for (auto it = fSystematicShifts.begin(); it != fSystematicShifts.end();
         ++it) {
      for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
        int selection_ID = it2->first;
        std::vector<TH1*> hists = it2->second;
        TH1D * selection_hist = (TH1D*)fSelectionHists[selection_ID];
        for (int i = 1; i <= selection_hist->GetNbinsX(); ++i) {
          std::vector<double> vars;
          for (size_t j = 0; j < hists.size(); ++j) {
            if (selection_hist->GetBinContent(i) < 1.e-5) {
              vars.push_back(1.);
            }
            else {
              vars.push_back(
                  hists[j]->GetBinContent(i)/selection_hist->GetBinContent(i));
            }
          }
          //if (do_insert)
            vars.insert(vars.begin() + vars.size()/2, 1.);
          std::string spline_name = selection_hist->GetName();
          spline_name += "_" + syst_name + "_Spline_" + std::to_string(i);
          fSystematicSplines[syst_name][selection_ID].push_back(
              new TSpline3(spline_name.c_str(), &fSystematicVals[syst_name][0],
                           &vars[0], vars.size()));
        }
      }
    }
  };

/*
  void FillIncidentHist(const std::vector<double> & vals) {
    for (size_t i = 0; i < vals.size(); ++i) {
      fIncidentHist.Fill(vals.at(i));
    }
  };*/

  void AddIncidentEnergies(const std::vector<double> & vals, double weight = 1.) {
    for (auto v : vals)
      fIncidentEnergies.push_back({v, weight});
  };

  void AddESliceEnergies(const std::pair<double, double> & vals, double weight = 1.) {
    fESliceEnergies.push_back({vals, weight});
  };

  void FillTrueIncidentHist(const std::vector<double> & vals, double weight = 1.) {
    for (size_t i = 0; i < vals.size(); ++i) {
      fTrueIncidentHist.Fill(vals.at(i), weight);
    }
  };

  void FillSelectionHist(int id, double val, double weight = 1.) {
    if (fSelectionHists.find(id) != fSelectionHists.end()) {
      fSelectionHists.at(id)->Fill(val, weight);
    }
  };

  template <size_t N> void FillSelectionHist(int id, const double (& vals)[N],
                                             double weight = 1.) {
    if (N < 1 || N > 3) {
      std::string message = "Error: trying to fill hists with too many values";
      message += std::to_string(N);
      throw std::runtime_error(message);
    }

    if (fSelectionHists.find(id) != fSelectionHists.end()) {
      if (N == 1) {
        fSelectionHists.at(id)->Fill(vals[0], weight);
      }
      else if (N == 2) {
        ((TH2D*)fSelectionHists.at(id))->Fill(vals[0], vals[1], weight);
      }
      else if (N == 3) {
        ((TH3D*)fSelectionHists.at(id))->Fill(vals[0], vals[1], vals[2], weight);
      }
    }
  }

  void FillHistFromIncidentEnergies(TH1D & hist) {
    for (auto vals : fIncidentEnergies) {
      //hist.Fill(vals.first, fFactor/*vals.second*/);
      hist.Fill(vals.first, fFactor*vals.second);
    }
  };

  void FillESliceHist(TH1D & hist) {
    for (auto e : fESliceEnergies) {
      std::pair<double, double> vals = e.first;
      double w = e.second;
      int last_bin = hist.FindBin(vals.first);
      int first_bin = hist.FindBin(vals.second);
      for (int i = first_bin; i <= last_bin; ++i)
        //hist.AddBinContent(i, fFactor);
        hist.AddBinContent(i, fFactor*w);
    }
  };

  void ScaleHists(double val) {
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

  void ScaleESliceEnergies(double val) {
    for (auto it = fESliceEnergies.begin();
         it != fESliceEnergies.end(); ++it) {
      it->second *= val;
    }
  };

  void Reset() {
    fVariedFlux = 0.;
    for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
      it->second->Reset();
    }

    fTrueIncidentHist.Reset();

    fIncidentEnergies.clear();
    fESliceEnergies.clear();
    fFactor = 1.;
  };

  void ScaleVariedFlux(double val) {
    fVariedFlux *= val;
  };

  void ScaleToDataMC() {
    ScaleHists(fDataMCScale);
    ScaleIncidentEnergies(fDataMCScale);
    ScaleESliceEnergies(fDataMCScale);
    ScaleVariedFlux(fDataMCScale);
  }

  void SetDataMCScale(double val) {
    fDataMCScale = val;
    ScaleHists(fDataMCScale);
    fNominalFlux *= val;
    fVariedFlux *= val;
  };

  void SetFactorAndScale(double val) {
    ResetFactor();
    fFactor = val;
    fVariedFlux *= val;
    ScaleHists(val);
  };

  void ResetFactor() {
    ScaleHists(1./fFactor);
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
  //TH1D fIncidentHist;
  TH1D fTrueIncidentHist;
  std::map<int, TH1 *> fSelectionHistsRebinned;
  //TH1D fIncidentHistRebinned;

  bool fMadeRebinned = false;

  //Consider changing to a class itself
  //       syst name             selection id   
  // string in 1st map: syst name
  // int in 2nd map: selection ID
  // index of vector: bin of corresponding selection hist
  std::map<std::string, std::map<int, std::vector<TSpline3 *>>>
      fSystematicSplines;
  std::map<std::string, std::map<int, std::vector<TH1 *>>>
      fSystematicShifts;
  std::map<std::string, std::vector<double>> fSystematicVals;


  std::vector<std::pair<double, double>> fIncidentEnergies;
  std::vector<std::pair<std::pair<double, double>, double>> fESliceEnergies;

};

}
#endif
