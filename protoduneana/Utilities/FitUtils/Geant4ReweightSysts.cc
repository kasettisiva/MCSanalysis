#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "cetlib/filepath_maker.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"

#include "ProtoDUNESelectionUtils.h"

#include "TFile.h"

//Necessary defs
std::string output_file_name; 
std::string fcl_file_name;
//////////////////////////////////////////

//Functions
bool parseArgs(int argc, char ** argv);
fhicl::ParameterSet makeParameterSet(std::string file);
void saveHists(std::vector<TH1 *> & vec, std::string addToName);
void makeFracErrorHist(TH1 * variation, TH1 * nominal);
/////////////////////////////////////////

typedef std::vector<std::string> string_v;
typedef std::vector<int> int_v;
typedef std::vector<double> double_v;

struct systConfig {

  systConfig(const fhicl::ParameterSet & pset)
    : MCFileNames(pset.get<string_v>("MCFileNames")),
      BackgroundTopologyName(
          pset.get<string_v>("BackgroundTopologyName")),
      SignalTopologyName(
          pset.get<string_v>("SignalTopologyName")),
      ChannelNames(pset.get<string_v>("ChannelNames")),
      IncidentMCFileNames(
          pset.get<string_v>("IncidentMCFileNames")),
      IncidentTopologyName(
          pset.get<string_v>("IncidentTopologyName")),
      RecoBinning(pset.get<double_v>("RecoBinning")),
      TruthBinning(pset.get<double_v>("TruthBinning")),
      RecoTreeName(pset.get<std::string>("RecoTreeName")),
      TruthTreeName(pset.get<std::string>("TruthTreeName")),
      SignalTopology(pset.get<int_v>("SignalTopology")),
      BackgroundTopology(pset.get<int_v>("BackgroundTopology")),
      IncidentTopology(pset.get<int_v>("IncidentTopology")),
      SystToConsider(pset.get<std::string>("SystToConsider")),
      EndZCut(pset.get<double>("EndZCut")) {};

  string_v MCFileNames;
  string_v BackgroundTopologyName;
  string_v SignalTopologyName;
  string_v ChannelNames;
  string_v IncidentMCFileNames;
  string_v IncidentTopologyName;
  double_v RecoBinning;
  double_v TruthBinning;
  std::string RecoTreeName;
  std::string TruthTreeName;

  int_v SignalTopology;
  int_v BackgroundTopology;
  int_v IncidentTopology;

  std::string SystToConsider;

  double EndZCut;
};


int main(int argc, char ** argv) {

  if (!parseArgs(argc, argv)) return 0; 

  fhicl::ParameterSet pset = makeParameterSet(fcl_file_name);

  systConfig cfg(pset);
  double tmin = cfg.TruthBinning[0];
  double tmax = cfg.TruthBinning.back();

  std::vector<TH1 *> bg_syst_hists;
  std::vector<TH1 *> bg_nom_hists;
  std::vector<TH1 *> sig_syst_hists;
  std::vector<TH1 *> sig_nom_hists;

  for (size_t i = 0; i < cfg.ChannelNames.size(); ++i) {
  
    //Backgrounds
    for (size_t j = 0; j < cfg.BackgroundTopology.size(); ++j) {
      int topo = cfg.BackgroundTopology[j];

      //Get nominal
      bg_nom_hists.push_back(
          protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(
              cfg.MCFileNames[i], cfg.RecoTreeName, cfg.RecoBinning,
              cfg.ChannelNames[i], cfg.BackgroundTopologyName[j], topo,
              cfg.EndZCut, tmin, tmax, false/*cfg.DoNegativeReco*/, 0));

      //Get Variations
      for (auto const k : {-1, 1}) {
        bg_syst_hists.push_back(
            protoana::ProtoDUNESelectionUtils::FillMCBackgroundHistogram_Pions(
                cfg.MCFileNames[i], cfg.RecoTreeName, cfg.RecoBinning,
                cfg.ChannelNames[i], cfg.BackgroundTopologyName[j], topo,
                cfg.EndZCut, tmin, tmax, false/*cfg.DoNegativeReco*/, /*doSyst=*/k));

        //Turn into frac errors
        makeFracErrorHist(bg_syst_hists.back(), bg_nom_hists.back());
      }
    }

    //Signals
    for (size_t j = 0; j < cfg.SignalTopology.size(); ++j) {
      int topo = cfg.SignalTopology[j];
      for (size_t k = 1; k < cfg.TruthBinning.size(); ++k) {
        
        //Get Nominal
        sig_nom_hists.push_back(
            protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(
                cfg.MCFileNames[i], cfg.RecoTreeName, cfg.RecoBinning,
                cfg.ChannelNames[i], cfg.SignalTopologyName[j], topo,
                cfg.TruthBinning[k-1], cfg.TruthBinning[k], cfg.EndZCut,
                false/*cfg.DoNegativeReco*/));

        //Get Variations
        for (auto const m : {-1, 1}) {
          sig_syst_hists.push_back(
              protoana::ProtoDUNESelectionUtils::FillMCSignalHistogram_Pions(
                  cfg.MCFileNames[i], cfg.RecoTreeName, cfg.RecoBinning,
                  cfg.ChannelNames[i], cfg.SignalTopologyName[j], topo,
                  cfg.TruthBinning[k-1], cfg.TruthBinning[k], cfg.EndZCut,
                  false/*cfg.DoNegativeReco*/, m));

          //Turn into frac errors
          makeFracErrorHist(sig_syst_hists.back(), sig_nom_hists.back());
        }
      }
    }
  }


  TFile output(output_file_name.c_str(), "RECREATE");
  output.cd();
  saveHists(bg_syst_hists, cfg.SystToConsider);
  saveHists(sig_syst_hists, cfg.SystToConsider);
  output.Close();
  return 0;
}

void saveHists(std::vector<TH1 *> & vec, std::string addToName) {
  for (auto const * h : vec) {
    std::string name = h->GetName();
    h->Write((name + "_" + addToName).c_str());
  }
}

void makeFracErrorHist(TH1 * variation, TH1 * nominal) {
  // Scale by -1 to subtract
  nominal->Scale(-1.);

  // Subtract from var
  variation->Add(nominal);

  // Scale back to normal
  nominal->Scale(-1.);

  // Divide to get fractional error
  variation->Divide(nominal);  
}

bool parseArgs(int argc, char ** argv) {

  bool found_fcl_file = false;
  bool found_output_file = false;

  for (int i = 1; i < argc; ++i) {
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      std::cout << "Usage: ./Geant4ReweightSysts -c <fitter_config>.fcl " <<
                   " -o <output_file>.root" << std::endl;
      std::cout << std::endl;

      return false;
    }

    else if (strcmp(argv[i], "-c") == 0) {
      fcl_file_name = argv[i+1];
      found_fcl_file = true;
    }

    else if (strcmp(argv[i], "-o") == 0) {
      output_file_name = argv[i+1];
      found_output_file = true;
    }

  }

  if (!found_fcl_file || !found_output_file) {
    std::cout << "Error: Must provide fcl file with options '-c' and '-o'" <<
                 std::endl;
    return false;
  }

  return true;
}

fhicl::ParameterSet makeParameterSet(std::string file) {
  fhicl::ParameterSet pset;

  // Configuration file lookup policy.
  char const* fhicl_env = getenv("FHICL_FILE_PATH");
  std::string search_path;

  if (fhicl_env == nullptr) {
    std::cerr << "Expected environment variable FHICL_FILE_PATH is missing" <<
                 " or empty: using \".\"\n";
    search_path = ".";
  }
  else {
    search_path = std::string{fhicl_env};
  }

  cet::filepath_first_absolute_or_lookup_with_dot lookupPolicy{search_path};
  fhicl::make_ParameterSet(file, lookupPolicy, pset);
  return pset;
}
