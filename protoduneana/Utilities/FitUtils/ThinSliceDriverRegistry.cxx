#include "ThinSliceDriverRegistry.h"
#include "ThinSliceDriver.h"
#include "ThinSliceDriverFactory.h"

#include <iostream>
#include <stdexcept>

protoana::ThinSliceDriverRegistry * protoana::ThinSliceDriverRegistry::fInstance = 0;

protoana::ThinSliceDriverRegistry * protoana::ThinSliceDriverRegistry::Instance() {
  if(!fInstance) {
    static protoana::ThinSliceDriverRegistry * manager_ptr = 0 ;
    if (!manager_ptr) {
      manager_ptr = new protoana::ThinSliceDriverRegistry;
    }
    protoana::ThinSliceDriverRegistry & manager = * manager_ptr;
    fInstance = & manager;
  }
  return fInstance;
}

protoana::ThinSliceDriverRegistry::ThinSliceDriverRegistry() {}

protoana::ThinSliceDriverRegistry::~ThinSliceDriverRegistry() {
  //Clean();
}

void protoana::ThinSliceDriverRegistry::PrintAvailableDrivers() const {
  std::cout << "####ThinSliceDriverRegistry####" << std::endl;
  std::cout << "Available Drivers:" << std::endl;
  for (auto it = fFactories.begin(); it != fFactories.end(); ++it) {
    std::cout << "Driver: " << it->first << std::endl;
  }
  std::cout << "###############################" << std::endl << std::endl;
}

void protoana::ThinSliceDriverRegistry::AddFactory(
    std::string name, protoana::BaseThinSliceDriverFactory * factory) {
  fFactories[name] = factory;
}

protoana::ThinSliceDriver * protoana::ThinSliceDriverRegistry::GetDriver(
    const std::string & name, const std::string & analysis) {
  if (fFactories.find(name) == fFactories.end()) {
    std::string message = "Could not find ThinSliceDriver of type: " +
                          name;
    throw std::runtime_error(message);
  }
  else {
    return fFactories[name]->Instantiate(analysis);
  }
}
