#include <random>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft libraries
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// DUNETPC specific includes
//#include "dune/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"
//#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"
#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/SimChannelExtractService.h"

#include "c2numpy.h"

using std::string;
using std::cout;
using std::endl;
using std::ofstream;

namespace pdune {
  class RawWaveformDump;
}

class pdune::RawWaveformDump : public art::EDAnalyzer {

public:

  explicit RawWaveformDump(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  RawWaveformDump(RawWaveformDump const&) = delete;
  RawWaveformDump(RawWaveformDump&&) = delete;
  RawWaveformDump& operator=(RawWaveformDump const&) = delete;
  RawWaveformDump& operator=(RawWaveformDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void reconfigure(fhicl::ParameterSet const & p);

  void beginJob() override;
  void endJob() override;

private:

  std::string fDumpWaveformsFileName;

  std::string	fSimulationProducerLabel; ///< The name of the producer that tracked simulated particles through the detector
  std::string	fDigitModuleLabel; ///< module that made digits

  std::string fSelectGenLabel;
  std::string fSelectProcID;
  int fSelectPDGCode;
  std::string fPlaneToDump;
  double fMinParticleEnergyGeV;
  double fMinEnergyDepositedMeV;
  int fMinNumberOfElectrons;
  int fMaxNumberOfElectrons;
  art::ServiceHandle<geo::Geometry> fgeom;
  art::ServiceHandle<cheat::ParticleInventoryService> PIS;
  //art::ServiceHandle<SimChannelExtractService> m_pscx;
  detinfo::DetectorClocks const * fClks;

  std::default_random_engine rndm_engine;

  c2numpy_writer npywriter;
};

//-----------------------------------------------------------------------
struct genFinder{
  private:
    typedef std::pair<int, std::string> track_id_to_string;
    std::vector<track_id_to_string> track_id_map;
    std::set<std::string> generator_names;
    bool isSorted=false;

  public:
    void sort_now(){
      std::sort(this->track_id_map.begin(), this->track_id_map.end(), [](const auto &a, const auto &b){return (a.first < b.first) ; } );
      isSorted=true;
    }
    void add(const int& track_id, const std::string& gname){
      this->track_id_map.push_back(std::make_pair(track_id, gname));
      generator_names.emplace(gname);
      isSorted=false;
    }
    bool has_gen(std::string gname){
      return static_cast<bool>(generator_names.count(gname));
    };
    std::string get_gen(int tid){
      if( !isSorted ){
	this->sort_now();
      }
      return std::lower_bound(track_id_map.begin(), track_id_map.end(), tid,[](const auto &a, const auto &b){return (a.first < b) ; } )->second;
    };

};
genFinder* gf = new genFinder();

//-----------------------------------------------------------------------
pdune::RawWaveformDump::RawWaveformDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fClks(lar::providerFrom<detinfo::DetectorClocksService>())
{
  this->reconfigure(p);
}

//-----------------------------------------------------------------------
void pdune::RawWaveformDump::reconfigure(fhicl::ParameterSet const & p)
{
  fDumpWaveformsFileName = p.get<std::string>("DumpWaveformsFileName","dumpwaveforms");

  fSimulationProducerLabel = p.get<std::string>("SimulationProducerLabel", "largeant");
  fDigitModuleLabel = p.get<std::string>("DigitModuleLabel", "daq");

  fSelectGenLabel = p.get<std::string>("SelectGenLabel","ANY");
  fSelectProcID   = p.get<std::string>("SelectProcID",  "ANY");
  fSelectPDGCode  = p.get<int>("SelectPDGCode", 0);

  fPlaneToDump = p.get<std::string>("PlaneToDump", "U");
  fMinParticleEnergyGeV = p.get<double>("MinParticleEnergyGeV",0.);
  fMinEnergyDepositedMeV = p.get<double>("MinEnergyDepositedMeV",0.);
  fMinNumberOfElectrons = p.get<int>("MinNumberOfElectrons",1000);
  fMaxNumberOfElectrons = p.get<int>("MaxNumberOfElectrons",100000);

  return;
}

//-----------------------------------------------------------------------
void pdune::RawWaveformDump::beginJob()
{
    std::random_device rndm_device;	// this will give us our seed
    rndm_engine.seed(rndm_device());

    c2numpy_init(&npywriter, fDumpWaveformsFileName, 50000);
    c2numpy_addcolumn(&npywriter, "evt",   C2NUMPY_UINT32);
    c2numpy_addcolumn(&npywriter, "gen",   (c2numpy_type)((int)C2NUMPY_STRING + 6));
    c2numpy_addcolumn(&npywriter, "trkid", C2NUMPY_INT32);
    c2numpy_addcolumn(&npywriter, "pdg",   C2NUMPY_INT32);
    c2numpy_addcolumn(&npywriter, "edepo", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&npywriter, "nelec", C2NUMPY_UINT32);
    c2numpy_addcolumn(&npywriter, "procid",(c2numpy_type)((int)C2NUMPY_STRING + 7));
    c2numpy_addcolumn(&npywriter, "chan",  C2NUMPY_UINT32);
    c2numpy_addcolumn(&npywriter, "view",  (c2numpy_type)((int)C2NUMPY_STRING + 1) );
    c2numpy_addcolumn(&npywriter, "stck_i", C2NUMPY_UINT16);
    c2numpy_addcolumn(&npywriter, "stck_f", C2NUMPY_UINT16);
    for(unsigned int i=0; i<6000; i++){
      std::ostringstream name;
      name << "tck_" << i;
      c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_INT16);
    } 
}

//-----------------------------------------------------------------------
void pdune::RawWaveformDump::endJob()
{
    c2numpy_close(&npywriter);
}

//-----------------------------------------------------------------------
void pdune::RawWaveformDump::analyze(art::Event const& evt)
{
  cout << "Event "
       << " " << evt.id().run()
       << " " << evt.id().subRun()
       << " " << evt.id().event() << endl;

  // ... Read in the digit List object(s). 
  art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
  evt.getByLabel(fDigitModuleLabel, digitVecHandle);

  if (!digitVecHandle->size())  return;

  // ... Use the handle to get a particular (0th) element of collection.
  art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);      
  unsigned int dataSize = digitVec0->Samples(); //size of raw data vectors
  if (dataSize!=6000){
    std::cout << "!!!!! Bad dataSize: " << dataSize << std::endl;
    return;
  } 

  std::vector<short> rawadc(dataSize);  // vector to hold uncompressed adc values later

  // ... Build a map from channel number -> rawdigitVec
  std::map< raw::ChannelID_t, art::Ptr<raw::RawDigit> > rawdigitMap;
  raw::ChannelID_t chnum = raw::InvalidChannelID;	// channel number
  for ( size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter ) {
    art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
    chnum = digitVec->Channel();
    if(chnum==raw::InvalidChannelID)continue;
    rawdigitMap[chnum]=digitVec;
  }

  // ... Read in MC particle list
  art::Handle< std::vector<simb::MCParticle> > particleHandle;
  if (!evt.getByLabel(fSimulationProducerLabel, particleHandle)){
    throw cet::exception("AnalysisExample")
      << " No simb::MCParticle objects in this event - "
      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  // ... Read in sim channel list
  auto simChannelHandle = evt.getValidHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);

  if (!simChannelHandle->size())  return;
  unsigned int simchansize = simChannelHandle->size();

  // ... Create a map of track IDs to generator labels
  //Get a list of generator names.
  std::vector< art::Handle< std::vector< simb::MCTruth > > > mcHandles;
  evt.getManyByType(mcHandles);
  std::vector< std::pair<int, std::string>> track_id_to_label;

  for( auto const& mcHandle : mcHandles ){
    const std::string& sModuleLabel = mcHandle.provenance()->moduleLabel();
    art::FindManyP<simb::MCParticle> findMCParts(mcHandle, evt, "largeant");
    std::vector<art::Ptr<simb::MCParticle> > mcParts = findMCParts.at(0);
    for( const art::Ptr<simb::MCParticle> ptr : mcParts){
      int track_id = ptr->TrackId();
      gf->add(track_id, sModuleLabel);
    }
  }

  // ... Loop over particles
  int sigchancount=0;
  for ( auto const& particle : (*particleHandle) ){

    int pdgcode = particle.PdgCode();
    int trackid = particle.TrackId();
    string procid = particle.Process();
    
    // .. get trackid of primary that ultimately created the particle that made current trackid.
    int eve_id = PIS->TrackIdToEveTrackId(trackid);
    std::string genlab = gf->get_gen(eve_id);
    
    // .. make sure particle comes from primary neutrino interaction vertex
    //if ( particle.Process() != "primary"  ||  pdgcode != 11 || particle.Mother() !=0 || particle.E() < 0.001)continue;
    if (fSelectGenLabel!="ANY"){
      if(genlab!=fSelectGenLabel)continue;
    }
    if (fSelectProcID!="ANY"){
      if(procid!=fSelectProcID)continue;
    }
    if (fSelectPDGCode!=0){
      if(pdgcode!=fSelectPDGCode)continue;
    }

    if ( particle.E() < fMinParticleEnergyGeV)continue;

    genlab.resize(6,' ');
    procid.resize(7,' ');

    bool foundsignal4P=false;
    std::vector<raw::ChannelID_t> sigchannels;
    std::vector<std::pair<unsigned int,unsigned int>>sigticks;
    std::vector<double> edepvec;
    std::vector<int> numelvec;
    sigchannels.reserve(simchansize);
    sigticks.reserve(simchansize);
    edepvec.reserve(simchansize);
    numelvec.reserve(simchansize);

    // ... Loop over simChannels
    for ( auto const& channel : (*simChannelHandle) ){
      bool foundsignal4CH=false;

      // .. get simChannel channel number
      const raw::ChannelID_t ch1 = channel.Channel();
      if(ch1==raw::InvalidChannelID)continue;
      if(geo::PlaneGeo::ViewName(fgeom->View(ch1))!=fPlaneToDump[0])continue;
      
      // .. now find RawDigit with this channel number
      auto search = rawdigitMap.find( ch1 );
      if ( search == rawdigitMap.end() ) continue;
      const raw::ChannelID_t ch2 = (*search).first;
      //art::Ptr<raw::RawDigit> rawdig = (*search).second;
      if(ch1!=ch2){
        cout << "!!!!!!!!!!!!!!!!!!!!! ch1 is not equal to ch2" << endl;
	exit(1);
      }

      // Create vector that holds the floating ADC count for each tick.
      //std::vector<AdcSignal> fChargeWork;
      //m_pscx->extract(&channel, fChargeWork);

      // ... Loop over all ticks with ionization energy deposited
      unsigned int tdcmin=dataSize-1;
      unsigned int tdcmax=0;
      int numel=0;
      double edep=0.,edep1=0.;
      auto const& timeSlices = channel.TDCIDEMap();
      for ( auto const& timeSlice : timeSlices ){

        auto const& energyDeposits = timeSlice.second;
	auto const tpctime = timeSlice.first;
        unsigned int tdctick = static_cast<unsigned int>(fClks->TPCTDC2Tick(double(tpctime)));
	if(tdctick!=tpctime)std::cout << "tpctime: " << tpctime << ", tdctick: " << tdctick << std::endl;
	if(tdctick<0||tdctick>(dataSize-1))continue;

	// ... Loop over all energy depositions in this tick
	for ( auto const& energyDeposit : energyDeposits ){

          edep+=energyDeposit.energy;
	  
	  // .. check to see if ide came from the primary particle
          if ( energyDeposit.trackID != trackid ) continue;
	  if (tdctick<tdcmin)tdcmin=tdctick;
	  if (tdctick>tdcmax)tdcmax=tdctick;
	  foundsignal4P=true;
          foundsignal4CH=true;
          edep1+=energyDeposit.energy;
	  numel+=energyDeposit.numElectrons;
        }
      }
      if(foundsignal4CH){
        if(numel>=fMinNumberOfElectrons && edep1>=fMinEnergyDepositedMeV){
	  if(fMaxNumberOfElectrons>=0 && numel>=fMaxNumberOfElectrons){
	    continue;
          } else {
	    sigchannels.push_back(ch2);
	    sigticks.emplace_back(tdcmin,tdcmax);
	    edepvec.push_back(edep1);
	    numelvec.push_back(numel);
	  }
	}
      }
    }
    if(foundsignal4P && sigchannels.size()>0){
      std::uniform_int_distribution<int> rndm_dist(0,sigchannels.size()-1);
      int i = rndm_dist(rndm_engine); // randomly select one channel with a signal from this particle
      chnum=sigchannels[i];
      unsigned int stck1 = sigticks[i].first;
      unsigned int stck2 = sigticks[i].second;
      double edepo=edepvec[i];
      int numelec=numelvec[i];
      auto search = rawdigitMap.find(chnum);
      if ( search == rawdigitMap.end() ) continue;
      art::Ptr<raw::RawDigit> rawdig = (*search).second;
      raw::Uncompress(rawdig->ADCs(), rawadc, rawdig->GetPedestal(), rawdig->Compression());
      c2numpy_uint32(&npywriter, evt.id().event());
      c2numpy_string(&npywriter, genlab.c_str());
      c2numpy_int32(&npywriter, trackid);
      c2numpy_int32(&npywriter, pdgcode);
      c2numpy_float32(&npywriter, edepo);
      c2numpy_uint32(&npywriter, numelec);
      c2numpy_string(&npywriter, procid.c_str());
      c2numpy_uint32(&npywriter, chnum);
      c2numpy_string(&npywriter, geo::PlaneGeo::ViewName(fgeom->View(chnum)).c_str());
      c2numpy_uint16(&npywriter, stck1);
      c2numpy_uint16(&npywriter, stck2);
      for ( unsigned int itck=0; itck<dataSize; ++itck ){
        rawadc[itck] -= rawdig->GetPedestal();
        c2numpy_int16(&npywriter, rawadc[itck]);
      }
      sigchancount++;
    }
  }
}
DEFINE_ART_MODULE(pdune::RawWaveformDump)
