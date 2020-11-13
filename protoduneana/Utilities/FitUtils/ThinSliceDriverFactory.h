#ifndef THINSLICEDRIVERFACTORY_H
#define THINSLICEDRIVERFACTORY_H

#include "ThinSliceDriverRegistry.h"
#include "ThinSliceDriver.h"

#include <string>

namespace protoana {
class BaseThinSliceDriverFactory {
 public:
  virtual ThinSliceDriver * Instantiate(
      const fhicl::ParameterSet & extra_options) = 0;
};

template <typename T> class ThinSliceDriverFactory
    : public BaseThinSliceDriverFactory {
 public:

  ThinSliceDriverFactory(const std::string name) {
    ThinSliceDriverRegistry::Instance()->AddFactory(name, this);
  }

  virtual ThinSliceDriver * Instantiate(
      const fhicl::ParameterSet & extra_options) {
    return new T(extra_options);
  }
};
}

#define DECLARE_THINSLICEDRIVER_FACTORY(driver) \
  const ThinSliceDriverFactory<driver>& driver##Factory = ThinSliceDriverFactory<driver>(#driver)

// support for drivers  defined within a namespace
// a bit tricky due to cpp macro expansion and the use of "::"
// use  DECLARE_THINSLICEDRIVER_FACTORY_NS( myns::MyDriver, myns, driverbase )  // without trailing ";"
#define DECLARE_THINSLICEDRIVER_FACTORY_NS( driver, nsname, driverbase )  \
  namespace nsname { \
    const ThinSliceDriverFactory<driver>& driverbase##Factory = ThinSliceDriverFactory<driver>(#driver); \
  }

#endif
