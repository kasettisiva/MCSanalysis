#ifndef PROTODUNESELECTIONUTILS_h
#define PROTODUNESELECTIONUTILS_h

#include <iostream>
#include <vector>
#include <string>

#include <TH1.h>

namespace protoana{
  namespace ProtoDUNESelectionUtils{

    // Fill signal and background histograms. The event selection is done at this stage.
    TH1* FillMCBackgroundHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, int channel, int topo);
    TH1* FillMCSignalHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, int channel, int topo, double minval, double maxval);

    // Fill the number of truth beam pions
    TH1* FillMCTruthSignalHistogram_Pions(std::string filename, std::string treename, std::vector<double> truthBins, int channel);

    // Fill data
    TH1* FillDataHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, int channel);

  }
}

#endif
