////////////////////////////////////////////////////////////////////////
// Class:       DataDump
// Plugin Type: analyzer (art v3_03_01)
// File:        DataDump_module.cc
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

namespace pdune {
  class DataDump;
}


class pdune::DataDump : public art::EDAnalyzer {
public:
  explicit DataDump(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DataDump(DataDump const&) = delete;
  DataDump(DataDump&&) = delete;
  DataDump& operator=(DataDump const&) = delete;
  DataDump& operator=(DataDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

//geo::GeometryCore const* fGeom;
  TTree *fTree;
  unsigned int run;
  unsigned int event;
  double evttime;
  std::vector<unsigned short> channel;
  std::vector<unsigned short> tick;
  std::vector<float> adc;
};


pdune::DataDump::DataDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdune::DataDump::analyze(art::Event const& e)
{

//fGeom = &*(art::ServiceHandle<geo::Geometry>());

  run = e.run();
  event = e.id().event();

  art::Timestamp ts = e.time();
  if (ts.timeHigh() == 0){
    TTimeStamp tts(ts.timeLow());
    evttime = tts.AsDouble();
  }
  else{
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime = tts.AsDouble();
  }
  channel.clear();
  tick.clear();
  adc.clear();

  //std::ofstream outfile (Form("r%de%d.txt",run,event));

  auto const& wires =
    e.getValidHandle<std::vector<recob::Wire> >("caldata:dataprep");

  for (auto & wire : * wires){
    int channel_no = wire.Channel();
    std::cout<<"Channel = "<<channel_no<<std::endl;
    if (channel_no<2080 || channel_no > 2559) continue;
    //for (auto & adc : wire.Signal()){
    for (size_t i = 0; i < wire.Signal().size(); ++i){
      //outfile<<channel<<" "<<i<<" "<<wire.Signal()[i]<<std::endl;
      channel.push_back(channel_no);
      tick.push_back(i);
      adc.push_back(wire.Signal()[i]);
    }
  }
  fTree->Fill();
  //outfile.close();
}

void pdune::DataDump::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("spt","space point tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("channel", &channel);
  fTree->Branch("tick",&tick);
  fTree->Branch("adc",&adc);
}
DEFINE_ART_MODULE(pdune::DataDump)
