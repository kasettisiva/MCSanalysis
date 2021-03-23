#ifndef THINSLICESYSTEMATIC_hh
#define THINSLICESYSTEMATIC_hh
#include "fhiclcpp/ParameterSet.h"
namespace protoana {
class ThinSliceSystematic {
 public:
  ThinSliceSystematic(const fhicl::ParameterSet & pset)
    : fType(pset.get<std::string>("Name")),
      fVal(pset.get<double>("Central")),
      fCentral(pset.get<double>("Central")),
      fUpperLimit(pset.get<double>("UpperLimit")),
      fLowerLimit(pset.get<double>("LowerLimit")),
      fSigma(pset.get<double>("Sigma")),
      fOptions(pset.get<fhicl::ParameterSet>("Options")) {
    fName = "par_" + fType + "_syst";
  };

  ThinSliceSystematic(){};

  ~ThinSliceSystematic(){};

  const double GetValue() const {
    return fVal;
  };

  const double GetCentral() const {
    return fCentral;
  };

  const double GetUpperLimit() const {
    return fUpperLimit;
  };

  const double GetLowerLimit() const {
    return fLowerLimit;
  };

  const std::string GetName() const {
    return fName;
  };

  void SetValue(double v) {
    fVal = v;
  };

  void SetToNSigma(int n = 1) {
    SetValue(n*fSigma);
  }

  void SetToCentral() {
    SetValue(fCentral); 
  }

  template <typename T> const T GetOption(std::string name) const {
    return fOptions.get<T>(name);
  }

 private:

  std::string fType;
  double fVal;
  double fCentral;
  double fUpperLimit; 
  double fLowerLimit;
  double fSigma;

  fhicl::ParameterSet fOptions;

  std::string fName;
};
}
#endif
