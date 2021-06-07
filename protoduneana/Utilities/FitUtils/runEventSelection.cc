#include <iostream>
#include "cetlib_except/exception.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"
#include <ROOT/RDataFrame.hxx>

class testing {
 private:
  int num; 
 public:
  testing(int n) : num(n){}
  int operator()() const {
    std::cout << num << std::endl;
    return num;
  }
};

int main(int argc, char ** argv){

  std::string fcl_file;
  std::string output_file;
  // Options to run
  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     fcl_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-o")) {
      output_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << "Usage: runPDSPThinSliceFit -c fclfile.fcl " << 
                    "-o outputfile.root " << std::endl;
      return 1;
    }
  }

  char const* fhicl_env = getenv("FHICL_FILE_PATH");
  std::string search_path;
  if (fhicl_env == nullptr) {
    std::cerr << "Expected environment variable FHICL_FILE_PATH is missing or empty: using \".\"\n";
    search_path = ".";
  }
  else {
    search_path = std::string{fhicl_env};
  }

  cet::filepath_first_absolute_or_lookup_with_dot lookupPolicy{search_path};
  fhicl::ParameterSet pset;
  //= fhicl::ParameterSet::make(fcl_file, lookupPolicy);
  fhicl::make_ParameterSet(fcl_file, lookupPolicy, pset);
  std::string tree_name = pset.get<std::string>("TreeName"); 

  ROOT::RDataFrame frame(tree_name, pset.get<std::string>("MCFile"));
  testing testing1(1);
  testing testing2(2);
  auto mc = frame.Define("testing1", testing1)
                 .Define("testing2", testing2);
  std::cout << "Defined" << std::endl;

  mc.Snapshot("new_tree", "testing.root");

  return 0;
}
