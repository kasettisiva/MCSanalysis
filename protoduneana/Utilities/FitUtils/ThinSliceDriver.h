#ifndef THINSLICEDRIVER_hh
#define THINSLICEDRIVER_hh

#include <map>
#include <vector>

#include "TTree.h"
#include "TH1D.h"

#include "ThinSliceSample.h"

namespace protoana {
class ThinSliceDriver {
 public:
  ThinSliceDriver(const std::string & analysis);
  virtual ~ThinSliceDriver();
  virtual void BuildDataHists(
      TTree * tree, TH1D & incident_hist,
      std::map<int, TH1 *> & selected_hists) = 0;
  virtual void BuildMCSamples(
      TTree * tree,
      std::map<int, std::vector<ThinSliceSample>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::map<int, double> & nominal_fluxes,
      std::map<int, std::vector<double>> & fluxes_by_sample) = 0;
  //virtual std::pair<double, size_t> CalculateChi2();
  std::string GetAnalysis() {return fAnalysis;};
 protected:
  std::string fAnalysis;
 private:
};
}
#endif
