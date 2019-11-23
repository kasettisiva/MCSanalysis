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

using wire_nt_t = hep_hpc::hdf5::Ntuple<float>;
using evt_nt_t = hep_hpc::hdf5::Ntuple<unsigned int, double>;

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
  void analyze(art::Event const& e) noexcept;
  void beginJob();
  virtual ~DataDumpHDF() noexcept;


protected:

//geo::GeometryCore const* fGeom;
  //TTree *fTree;
  double evttime;
//  std::vector<unsigned short> channel;
//  std::vector<unsigned short> tick;
//  std::vector<float> adc;
  hep_hpc::hdf5::File hdffile;

};


pdune::DataDumpHDF::DataDumpHDF(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}
pdune::DataDumpHDF::~DataDumpHDF() noexcept
{
}


void pdune::DataDumpHDF::analyze(art::Event const& e) noexcept
{

//fGeom = &*(art::ServiceHandle<geo::Geometry>());

  hdffile = hep_hpc::hdf5::File(Form("r%d_e%d.h5",e.run(),e.event()), H5F_ACC_TRUNC);
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

  //wire_nt_t wiresigs(hdffile, "wiresigs", {{"eid",4}});
  wire_nt_t wiresigs(hdffile, "wiresigs", {"adc"});
  evt_nt_t evtids(hdffile, "evtids", {{"eid",3}, "evttime"});
  evtids.insert(event_id.data(), evttime);

  auto const& wires =
    e.getValidHandle<std::vector<recob::Wire> >("caldata:dataprep");
   //std::vector<recob::Wire> const& wireVector(*wires);

   int fill_channel = 2080;
   for (auto & wire : * wires){
	   int channel = wire.Channel();
	   if (channel<2080 || channel > 2559) continue;
	   std::cout<<"Channel = "<<channel<<std::endl;
	   if (channel != fill_channel){
		   for  (int j = 0; j < 6000; j++)
			   wiresigs.insert(0.);
	   }
	   else{
		   int nticks = wire.Signal().size();
		   for (int j = 0; j < nticks; j++){
			   float adc = wire.Signal()[j];
			   wiresigs.insert(adc);
		   }
		   if (nticks <6000){
			   for (int j = nticks; j < 6000; j++)
				   wiresigs.insert(0.);
		   }
	   }
	   fill_channel++;
   }
   std::cout<<"event_time: "<<evttime<<std::endl;

}

void pdune::DataDumpHDF::beginJob()
{
}

DEFINE_ART_MODULE(pdune::DataDumpHDF)
