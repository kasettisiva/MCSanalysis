#ifndef PDSPTHINSLICEFITTER_hh
#define PDSPTHINSLICEFITTER_hh

#include <string>
#include <vector>
#include <map>

#include "TTree.h"
#include "TFile.h"
#include "THStack.h"

#include "ThinSliceSample.h"

namespace protoana {

class PDSPThinSliceFitter {
 public:
  PDSPThinSliceFitter(std::string fcl_file, std::string output_file);
  void BuildMCSamples();
  void SaveMCSamples();
  void BuildAndSaveNominalStacks();
  void BuildDataHists();
  void InitializeMCSamples();
  void CompareDataMC();
  ~PDSPThinSliceFitter();

 private:
  void Configure(std::string fcl_file);

  std::map<int, std::vector<ThinSliceSample>> fSamples;
  std::map<int, bool> fIsSignalSample;
  TFile fMCFile;
  TTree * fMCTree;
  TFile fDataFile;
  TTree * fDataTree;
  TFile fOutputFile;

  std::map<int, TH1D> fSelectedDataHists;
  TH1D fIncidentDataHist;

  THStack * fNominalIncidentMCStack;
  std::map<int, THStack *> fNominalSelectedMCStacks;

  //Configurable members
  std::string fMCFileName;
  std::string fDataFileName;
  std::string fTreeName;
  std::map<int, std::string> fSelectionIDs;
  std::vector<fhicl::ParameterSet> fSampleSets;
  
  std::vector<double> fSelectedRecoBins;
  std::vector<double> fIncidentRecoBins;
  //////////////////////////
};

}
#endif
