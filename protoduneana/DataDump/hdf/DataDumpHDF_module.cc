////////////////////////////////////////////////////////////////////////
// Class:       DataDumpHDF
// Plugin Type: analyzer (art v3_03_01)
// File:        DataDumpHDF_module.cc
//
// Generated at Tue Nov 19 13:37:52 2019 by Tingjun Yang using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"

#include "TTree.h"
#include "TTimeStamp.h"

#include <fstream>
#include <vector>

#include "hep_hpc/hdf5/File.hpp"
#include "hep_hpc/hdf5/Ntuple.hpp"

using wire_nt_t = hep_hpc::hdf5::Ntuple<unsigned int, float>;

inline std::array<unsigned int, 3>
get_eid(art::Event const& e)
{
  return { e.run(), e.subRun(), e.id().event() };
}

namespace pdune {
  class DataDumpHDF;
}

class pdune::DataDumpHDF : public art::EDAnalyzer {
public:
  explicit DataDumpHDF(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DataDumpHDF(DataDumpHDF const&) = delete;
  DataDumpHDF(DataDumpHDF&&) = delete;
  DataDumpHDF& operator=(DataDumpHDF const&) = delete;
  DataDumpHDF& operator=(DataDumpHDF&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

//geo::GeometryCore const* fGeom;
  //TTree *fTree;
  unsigned int run;
  unsigned int subrun;
  unsigned int event;
  double evttime;
//  std::vector<unsigned short> channel;
//  std::vector<unsigned short> tick;
//  std::vector<float> adc;
};


pdune::DataDumpHDF::DataDumpHDF(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdune::DataDumpHDF::analyze(art::Event const& e)
{

//fGeom = &*(art::ServiceHandle<geo::Geometry>());

  run = e.run();
  subrun = e.subRun();
  event = e.id().event();

  auto event_id = get_eid(e);

  art::Timestamp ts = e.time();
  if (ts.timeHigh() == 0){
    TTimeStamp tts(ts.timeLow());
    evttime = tts.AsDouble();
  }
  else{
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime = tts.AsDouble();
  }
//  channel.clear();
//  tick.clear();
//  adc.clear();

  //std::ofstream outfile (Form("r%de%d.txt",run,event));

  hep_hpc::hdf5::File hdffile(Form("r%de%d.h5",run,event), H5F_ACC_TRUNC);

  //wire_nt_t wiresigs(hdffile, "wiresigs", {{"eid",4}});
  wire_nt_t wiresigs(hdffile, "wiresigs", {{"eid",3}, "adc"});

  auto const& wires =
    e.getValidHandle<std::vector<recob::Wire> >("caldata:dataprep");
   std::vector<recob::Wire> const& wireVector(*wires);
/*
  for (auto & wire : * wires){
    int channel = wire.Channel();
    std::cout<<"Channel = "<<channel<<std::endl;
    if (channel<2080 || channel > 2559) continue;
    //for (auto & adc : wire.Signal()){
    for (size_t i = 0; i < wire.Signal().size(); ++i){
      //outfile<<channel<<" "<<i<<" "<<wire.Signal()[i]<<std::endl;
//      channel.push_back(channel_no);
//      tick.push_back(i);
//      adc.push_back(wire.Signal()[i]);
      //wiresigs.insert(run, subrun, event, wire.Signal()[i]);
      wiresigs.insert(event_id.data(), channel, i, wire.Signal()[i]);
    }
  }
*/

   // Note: the following code has two assumptions:
   //    1. channels as a whole should not be missing; 
   //    2. missing ticks are only the last few.

  for (int i=2080; i<=2559; i++){
    std::cout<<"Channel = "<<i<<std::endl;
    int nticks = wireVector.at(i).Signal().size();
    for (int j = 0; j < nticks; j++){
	float adc = wireVector.at(i).Signal()[j];
	wiresigs.insert(event_id.data(), adc);
    }
    if (nticks <6000){
        for (int j = nticks; j < 6000; j++){
            float adc = 0.;
	    wiresigs.insert(event_id.data(), adc);
        }
    }
  }

  
  //fTree->Fill();
  //outfile.close();
}

DEFINE_ART_MODULE(pdune::DataDumpHDF)
