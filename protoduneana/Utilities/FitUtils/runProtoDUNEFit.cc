#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "cetlib/filepath_maker.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"

#include "ProtoDUNEFit.h"

#include "TFile.h"

int main(int argc, char ** argv){

  std::string fhcfile;
  TString Outputfile;
  int analysis = -1;
  // Options to run
  for(int iArg = 1; iArg < argc; iArg++){
    if((!strcasecmp(argv[iArg],"-f")) || (!strcasecmp(argv[iArg],"-fhc"))) {
     fhcfile = argv[++iArg];
    }
    if((!strcasecmp(argv[iArg],"-o")) || (!strcasecmp(argv[iArg],"-output"))) {
      Outputfile = TString(argv[++iArg]);
    }
    if((!strcasecmp(argv[iArg],"-h")) || (!strcasecmp(argv[iArg],"-help"))) {
      std::cout << "Usage: runProtoDUNEFit -f fclfile.fcl -o outputfile.root [-pions]" << std::endl;
      return 1;
    }
    if((!strcasecmp(argv[iArg],"-pions")) || (!strcasecmp(argv[iArg],"--pions"))) {
      analysis = 1;
    }
  }

  protoana::ProtoDUNEFit *pfit = new protoana::ProtoDUNEFit(fhcfile);
  pfit->BuildWorkspace(Outputfile, analysis);

  return 0;
}
