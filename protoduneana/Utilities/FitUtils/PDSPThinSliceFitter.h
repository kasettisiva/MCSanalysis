#ifndef PDSPTHINSLICEFITTER_hh
#define PDSPTHINSLICEFITTER_hh

#include <string>
#include <vector>
#include <map>

#include "TTree.h"
#include "TFile.h"
#include "THStack.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#include "ThinSliceSample.h"

namespace protoana {

class PDSPThinSliceFitter {
 public:
  PDSPThinSliceFitter(std::string fcl_file, std::string output_file);
  void BuildMCSamples();
  void SaveMCSamples();
  void BuildAndSaveNominalStacks();
  void BuildAndSavePostFitStacks();
  void GetNominalFluxes();
  void BuildDataHists();
  void InitializeMCSamples();
  void CompareDataMC();
  void ScaleMCToData();
  void RunFitAndSave();
  ~PDSPThinSliceFitter();

 private:
  void Configure(std::string fcl_file);
  std::pair<double, size_t> CalculateChi2();
  void DefineFitFunction();
  void MakeMinimizer();

  std::map<int, std::vector<ThinSliceSample>> fSamples;
  std::map<int, bool> fIsSignalSample;
  TFile fMCFile;
  TTree * fMCTree;
  TFile fDataFile;
  TTree * fDataTree;
  TFile fOutputFile;
  ROOT::Math::Functor fFitFunction;
  std::unique_ptr<ROOT::Math::Minimizer> fMinimizer;

  std::map<int, TH1D> fSelectedDataHists;
  TH1D fIncidentDataHist;

  THStack * fNominalIncidentMCStack;
  THStack * fPostFitIncidentMCStack;
  std::map<int, THStack *> fNominalSelectedMCStacks;
  std::map<int, THStack *> fPostFitSelectedMCStacks;

  std::map<int, double> fNominalFluxes;
  std::map<int, std::vector<double>> fFluxesBySample;
  double fDataFlux;
  double fMCDataScale = 1.;

  std::map<int, std::vector<double>> fSignalParameters;
  std::map<int, std::vector<std::string>> fSignalParameterNames;
  size_t fTotalSignalParameters;

  std::map<int, double> fFluxParameters;
  std::map<int, std::string> fFluxParameterNames;
  size_t fTotalFluxParameters;

  //Configurable members
  std::string fMCFileName;
  std::string fDataFileName;
  std::string fTreeName;
  std::map<int, std::string> fSelectionIDs;
  std::vector<fhicl::ParameterSet> fSampleSets;
  std::map<int, std::string> fFluxTypes;
  int fMaxCalls;
  double fTolerance, fLowerLimit, fUpperLimit;
  
  std::vector<double> fSelectedRecoBins;
  std::vector<double> fIncidentRecoBins;
  //////////////////////////
};

}
#endif
