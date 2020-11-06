#ifndef PDSPTHINSLICEFITTER_hh
#define PDSPTHINSLICEFITTER_hh

#include <string>
#include <vector>
#include <map>

#include "TTree.h"
#include "TFile.h"

#include "ThinSliceSample.h"

namespace protoana {

class PDSPThinSliceFitter {
 public:
  PDSPThinSliceFitter(std::string fcl_file, std::string output_file);
  void BuildMCSamples();
  void SaveMCSamples();
  void BuildAndSaveStacks();
  void InitializeMCSamples();
  ~PDSPThinSliceFitter();

 private:
  void Configure(std::string fcl_file);

  std::map<int, std::vector<ThinSliceSample>> fSamples;
  std::map<int, bool> fIsSignalSample;
  TFile fMCFile;
  TTree * fMCTree;
  TFile fOutputFile;

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
