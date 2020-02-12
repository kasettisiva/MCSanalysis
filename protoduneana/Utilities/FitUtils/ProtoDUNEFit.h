#ifndef PROTOANA_PROTODUNEFIT_H
#define PROTOANA_PROTODUNEFIT_H

// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"

#include "RooStats/HistFactory/Measurement.h"

// C/C++ standard libraries
#include <string>
#include <vector>
#include <map>

namespace protoana{
  class ProtoDUNEFit{
  public:

    ProtoDUNEFit();
    ProtoDUNEFit(std::string configPath);
    virtual ~ProtoDUNEFit();

    // Create histograms, fit, and save
    void BuildWorkspace(TString Outputfile, int analysis=-1);

    // Apply systematics to sample
    bool ApplySystematicToSample(RooStats::HistFactory::Sample& sample, TH1* histo, std::vector<TH1*> systvec, bool hasnormfactor, bool isnorm);

  private:
    // Configure input from fcl file
    bool Configure(std::string configPath);

    // Read the input root file and create the histograms
    bool FillHistogramVectors_Pions();

    // Build samples and channels
    void AddSamplesAndChannelsToMeasurement(RooStats::HistFactory::Measurement& meas);
    void AddIncidentSamplesAndChannelsToMeasurement(RooStats::HistFactory::Measurement& meas);
    void AddSidebandSamplesAndChannelsToMeasurement(RooStats::HistFactory::Measurement& meas);

    std::string _RecoTreeName, _Minimizer, _TruthTreeName;

    std::vector<double> _RecoBinning, _TruthBinning;

    std::vector<int> _SignalTopology, _BackgroundTopology, _IncidentTopology;

    std::vector<std::string> _DataFileNames, _MCFileNames, _DataControlSampleFiles, _MCControlSampleFiles, _SystFileNames, _SystToConsider, _SystType, _BackgroundTopologyName, _SignalTopologyName, _ChannelNames, _IncidentMCFileNames, _IncidentTopologyName;

    int _FitStrategy, _NToys;
    double _IgnoreStatisticalErrorBelow, _IgnoreSystematicErrorBelow;
    bool _EnableMinosError, _DoAsimovFit, _EnableStatisticalError, _EnableSystematicError, _NormalisedSystematic;

    std::vector<TH1*> _bkghistos, _sighistos, _truthsighistos, _datahistos;
    std::vector<TH1*> _incbkghistos, _incsighistos, _incdatahistos;
    std::vector<TH1*> _sidhistos, _siddatahistos;

    std::vector<TGraphAsymmErrors*> _efficiencyGraphs;
  };
}

#endif
