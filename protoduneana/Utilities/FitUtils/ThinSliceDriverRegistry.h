#ifndef THINSLICEDRIVERREGISTRY_h
#define THINSLICEDRIVERREGISTRY_h

#include <map>
#include "fhiclcpp/ParameterSet.h"

namespace protoana {

class ThinSliceDriver;
class BaseThinSliceDriverFactory;

class ThinSliceDriverRegistry {
 public:
  static ThinSliceDriverRegistry * Instance();
  ~ThinSliceDriverRegistry();
  void AddFactory(std::string name, BaseThinSliceDriverFactory * factory);
  void PrintAvailableDrivers() const;
  ThinSliceDriver * GetDriver(
      const std::string & name,
      const fhicl::ParameterSet & extra_options);

 private:
  ThinSliceDriverRegistry();
  static ThinSliceDriverRegistry * fInstance;
  std::map<std::string, BaseThinSliceDriverFactory *> fFactories;
  
};
}
#endif
