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
#include <algorithm>

#include "hep_hpc/hdf5/File.hpp"
#include "hep_hpc/hdf5/Ntuple.hpp"

using wire_nt_t = hep_hpc::hdf5::Ntuple<float>;
using evt_nt_t = hep_hpc::hdf5::Ntuple<unsigned int, double>;

inline std::array<unsigned int, 3>
get_eid(art::Event const& e)
{
  return { e.run(), e.subRun(), e.id().event() };
}

bool chIncrease(art::Ptr<recob::Wire> const & w1,
                art::Ptr<recob::Wire> const & w2){
  return (w1->Channel()<w2->Channel());
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

  art::Handle < std::vector < recob::Wire > > wireListHandle;
  std::vector < art::Ptr < recob::Wire > > wires;
  if (e.getByLabel("caldata:dataprep", wireListHandle)) {
    art::fill_ptr_vector(wires, wireListHandle);
  }

  std::sort(wires.begin(), wires.end(), chIncrease);
//  auto const& wires =
//    e.getValidHandle<std::vector<recob::Wire> >("caldata:dataprep");
//   //std::vector<recob::Wire> const& wireVector(*wires);

  for (auto & wire : wires){
    int channel = wire->Channel();
    if (!((channel>=2080 && channel < 2560)||
          (channel>=7200 && channel < 7680)||
          (channel>=12320 && channel < 12800))) continue;
    if (channel%1000==0) std::cout<<"Channel = "<<channel<<std::endl;

    const recob::Wire::RegionsOfInterest_t& signalROI = wire->SignalROI();
    int lasttick = 0;
    for(const auto& range : signalROI.get_ranges()){
      const auto& waveform = range.data();
      // ROI start time
      raw::TDCtick_t roiFirstBinTick = range.begin_index();
      for (int i = lasttick; i<roiFirstBinTick; ++i){
        wiresigs.insert(0.);
      }
      for(size_t idx = 0; idx < waveform.size(); idx++){
        wiresigs.insert(waveform[idx]);
        ++lasttick;
      }
    }
    for (int i = lasttick; i<6000; ++i){
      wiresigs.insert(0.);
    }
  }
//    int nticks = wire->Signal().size();
//    for (int j = 0; j < nticks; j++){
//      float adc = wire->Signal()[j];
//      wiresigs.insert(adc);
//    }
//    if (nticks <6000){
//      for (int j = nticks; j < 6000; j++)
//        wiresigs.insert(0.);
//    }

  std::cout<<"event_time: "<<evttime<<std::endl;
  
}

void pdune::DataDumpHDF::beginJob()
{
}

DEFINE_ART_MODULE(pdune::DataDumpHDF)
