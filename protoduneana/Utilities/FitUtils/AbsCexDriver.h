#ifndef ABSCEXDRIVER_hh
#define ABSCEXDRIVER_hh

#include "ThinSliceDriver.h"
namespace protoana {
class AbsCexDriver : public ThinSliceDriver {
 public:
  AbsCexDriver(const std::string & analysis);
  virtual ~AbsCexDriver();
  void BuildDataHists(
      TTree * tree, TH1D & incident_hist,
      std::map<int, TH1 *> & selected_hists) override;
  void BuildMCSamples(
      TTree * tree,
      std::map<int, std::vector<ThinSliceSample>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::map<int, double> & nominal_fluxes,
      std::map<int, std::vector<double>> & fluxes_by_sample) override;
  //std::pair<double, size_t> CalculateChi2() override;
};
}
#endif
