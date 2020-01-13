#ifndef PROTODUNESELECTIONUTILS_h
#define PROTODUNESELECTIONUTILS_h

#include <iostream>
#include <vector>
#include <string>

#include <TH1.h>

namespace protoana{
  namespace ProtoDUNESelectionUtils{

    TH1* FillMCBackgroundHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, int channel, int topo);
    TH1* FillMCSignalHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, int channel, int topo, double minval, double maxval);
    //TH1* FillMCTruthSignalHistogram_Pions(std::string filename, std::string treename, int channel, double minval, double maxval);
    TH1* FillDataHistogram_Pions(std::string filename, std::string treename, std::vector<double> recoBins, int channel);

  }
}

#endif
