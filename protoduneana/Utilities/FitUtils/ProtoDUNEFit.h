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
#include "TCanvas.h"

#include "RooStats/HistFactory/Measurement.h"
#include <RooFitResult.h>

// C/C++ standard libraries
#include <string>
#include <vector>
#include <map>

namespace protoana{

  enum HistType {
    kSignal,
    kBackground,
    kIncident
  };

  class ProtoDUNEFit{
  public:

    ProtoDUNEFit();
    ProtoDUNEFit(std::string configPath);
    virtual ~ProtoDUNEFit();

    // Create histograms, fit, and save
    void BuildWorkspace(TString Outputfile, int analysis=-1);

    // Apply systematics to sample
    bool ApplySystematicToSample(RooStats::HistFactory::Sample& sample, TH1* histo, std::vector<TH1*> systvec, bool hasnormfactor, bool isnorm);

    bool BuildSystThenApplyToSample(RooStats::HistFactory::Sample& sample,
                                    TH1* histo, bool hasnormfactor,
                                    bool isnorm,
                                    protoana::HistType this_histType,
                                    size_t iChan, size_t iTopo,
                                    size_t iTruthBin=999);
    bool BuildSignalSystThenApplyToSample(
        RooStats::HistFactory::Sample& sample, TH1* histo,
        bool hasnormfactor, bool isnorm,
        size_t iChan, size_t iTopo, size_t iTruthBin);

    bool BuildBackgroundSystThenApplyToSample(
        RooStats::HistFactory::Sample& sample, TH1* histo,
        bool hasnormfactor, bool isnorm,
        size_t iChan, size_t iTopo);

    bool BuildIncidentSystThenApplyToSample(
        RooStats::HistFactory::Sample& sample, TH1* histo,
        bool hasnormfactor, bool isnorm, size_t iTopo);

    std::vector<std::vector<std::pair<TH1*, TH1*>>>
        BuildIncidentSignalSyst(size_t iTopo);

    bool ApplyBuiltSystToSample(TH1 * histo, TH1 * high_hist, TH1 * low_hist,
                                RooStats::HistFactory::Sample& sample,
                                std::string syst_name, std::string syst_type,
                                bool hasnormfactor);
    std::vector<TH1 *> DrawXSecs(RooFitResult *fitresult = 0x0);
    std::vector<TH2 *> DrawSmearingMatrix(RooFitResult *fitresult, bool doPostFit = false);


  private:
    // Configure input from fcl file
    bool Configure(std::string configPath);

    // Read the input root file and create the histograms
    bool FillHistogramVectors_Pions();

    // Scale MC to Data
    void ScaleMCToData(bool data_is_mc = false);

    // Fix Muon content
    void ScaleMuonContent();

    // Build samples and channels
    void AddSamplesAndChannelsToMeasurement(RooStats::HistFactory::Measurement& meas);
    void AddIncidentSamplesAndChannelsToMeasurement(RooStats::HistFactory::Measurement& meas);
    void AddSidebandSamplesAndChannelsToMeasurement(RooStats::HistFactory::Measurement& meas);

    void DecorateEfficiency(TGraphAsymmErrors * eff,
                            std::string x_title = "E_{true} at vertex [MeV]");

    std::string _RecoTreeName, _Minimizer, _TruthTreeName;

    std::vector<double> _RecoBinning, _TruthBinning;

    std::vector<int> _SignalTopology, _BackgroundTopology, _IncidentTopology;

    std::vector<std::string> _DataFileNames, _MCFileNames,
                             _DataControlSampleFiles, _MCControlSampleFiles,
                             _SystFileNames, _SystToConsider, _SystType,
                             _BackgroundTopologyName, _SignalTopologyName,
                             _ChannelNames, _IncidentMCFileNames,
                             _IncidentDataFileNames, _IncidentTopologyName;

    int _FitStrategy, _NToys;
    double _IgnoreStatisticalErrorBelow, _IgnoreSystematicErrorBelow;
    bool _EnableMinosError, _DoAsimovFit, _EnableStatisticalError;
    bool _EnableSystematicError, _NormalisedSystematic, _FitInReco;
    bool _UseComputedSysts;

    std::vector<TH1*> _bkghistos, _sighistos, _truthsighistos, _datahistos;
    std::vector<TH1*> _syst_hists;
    std::vector<size_t> _bkg_chan_index, _sig_chan_index;
    std::vector<size_t> _bkg_topo_index, _sig_topo_index, _sig_truth_index,
                        _inc_sig_topo_index, _inc_bkg_topo_index;
    std::vector<TH1*> _incbkghistos, _incsighistos, _incdatahistos;
    std::vector<TH1*> _sidhistos, _siddatahistos;

    std::vector<TGraphAsymmErrors*> _efficiencyGraphs;

    TGraphAsymmErrors * _incidentEfficiency;
    TH1 * _incidentEfficiencyNum;
    TH1 * _incidentEfficiencyDenom;
    std::vector< TH1 * > _interactingEfficiencyDenoms;
    std::vector< TH1 * > _interactingEfficiencyNums;
    std::vector< TGraphAsymmErrors * > _interactingEfficiencies;
 
    bool _AddIncidentToMeasurement, _DoNegativeReco, _DoScaleMCToData;
    bool _DoScaleMuonContent;
    bool _AddBackgroundFactors, _AddIncidentBackgroundFactors;
    bool _DistinguishIncidentSignal, _OnlyDrawXSecs;
    bool _DataIsMC;
    double _EndZCut, _WirePitch;

    double _ScaleFactor = 1.;
    double _IncidentScaleFactor = 1.;

    double _MuonScaleFactor = 1., _PionScaleFactor = 1.;
  };
}

#endif
