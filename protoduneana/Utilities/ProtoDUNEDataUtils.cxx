#include "ProtoDUNEDataUtils.h"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RawData/RawDigit.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

protoana::ProtoDUNEDataUtils::ProtoDUNEDataUtils(fhicl::ParameterSet const& p){
  this->reconfigure(p);
}

protoana::ProtoDUNEDataUtils::~ProtoDUNEDataUtils(){

}

void protoana::ProtoDUNEDataUtils::reconfigure(fhicl::ParameterSet const& p){
  fTimingTag            = p.get<art::InputTag>("TimingTag");
  fRawDigitTag          = p.get<art::InputTag>("RawDigitTag");
  fRawDigitTimeStampTag = p.get<art::InputTag>("RawDigitTimeStampTag");
}

// Access the trigger information to see if this is a beam trigger
bool protoana::ProtoDUNEDataUtils::IsBeamTrigger(art::Event const & evt) const{

  bool isBeam = false;

  // Accessing the trigger information as done in DataPrepModule
  // The information is stored in the time stamps
  art::Handle<std::vector<raw::RDTimeStamp>> timeStamps;
  evt.getByLabel(fTimingTag,timeStamps);

  // Return false if we have no time stamps
  if(!timeStamps.isValid()) return isBeam;
  // We should only have one RDTimeStamp
  if(timeStamps->size() > 1) return isBeam;
  
  // Access the trigger information. Beam trigger flag = 0xc
  const raw::RDTimeStamp& timeStamp = timeStamps->at(0);
  isBeam = (timeStamp.GetFlags() == 0xc);
  
  return isBeam;
}


// ----------------------------------------------------------------------------
int protoana::ProtoDUNEDataUtils::GetNActiveFembsForAPA(art::Event const & evt, int apa) const {


  // Get pd channel map
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;

  // set only saves unique elements
  std::set<int> apaset;

  // Get raw digits time stamps
  art::Handle< std::vector<raw::RDTimeStamp> > RawdigitTSListHandle;
  std::vector<art::Ptr<raw::RDTimeStamp> > digitTSlist;
  

  // Get raw digits
  art::Handle< std::vector<raw::RawDigit> > RawdigitListHandle;
  std::vector<art::Ptr<raw::RawDigit> > digitlist;

  if (evt.getByLabel(fRawDigitTag, RawdigitListHandle)){
  	
    art::fill_ptr_vector(digitlist, RawdigitListHandle);  

    for(auto const & dptr : digitlist) {
      const raw::RawDigit& digit = *dptr;
      // Get the channel number for this digit
      uint32_t chan = digit.Channel();
    
      int iapa = channelMap->APAFromOfflineChannel(chan);
      if(iapa != apa) continue;
      // Get the channel FEMB and WIB
      int WIB = channelMap->WIBFromOfflineChannel(chan); // 0-4
      int FEMB = channelMap->FEMBFromOfflineChannel(chan); // 1-4
      //int FEMBchan = channelMap->FEMBChannelFromOfflineChannel(chan);
      int iFEMB = ((WIB*4)+(FEMB-1)); //index of the FEMB 0-19

      apaset.insert(iFEMB);
    }
  }

  else{ // if raw digits have been dropped use RDTimeStamps instead
    evt.getByLabel(fRawDigitTimeStampTag, RawdigitTSListHandle);
	
    art::fill_ptr_vector(digitTSlist, RawdigitTSListHandle);  
	
    for(auto const & dptr : digitTSlist) {
  	
      const raw::RDTimeStamp & digit = *dptr;
    
      // Get the channel number for this digit
      uint16_t chan = digit.GetFlags();
    
      int iapa = channelMap->APAFromOfflineChannel(chan);
      if(iapa != apa) continue;
      // Get the channel FEMB and WIB
      int WIB = channelMap->WIBFromOfflineChannel(chan); // 0-4
      int FEMB = channelMap->FEMBFromOfflineChannel(chan); // 1-4
      //int FEMBchan = channelMap->FEMBChannelFromOfflineChannel(chan);
      int iFEMB = ((WIB*4)+(FEMB-1)); //index of the FEMB 0-19

      apaset.insert(iFEMB);
    }
  }
  return (apaset.size());
}


// ----------------------------------------------------------------------------
// returns true if timestamps are consistent, and fills the timestamp number in the timestamp argument.  timestamp2
//   should be the same as timestamp if this method returns true, and apainconsist should be zero.
// returns false if inconsistent, and the timestamp argument is one of the timestamps present, and timestamp2 is
//   another one if it disagrees.  apainconsist is the apa number of the first clashing channel

bool protoana::ProtoDUNEDataUtils::CheckTimeStampConsistencyForAPAs(art::Event const & evt, std::set<int> apas, 
								    ULong64_t &timestamp, ULong64_t &timestamp2,
								    int &apainconsist) const
{
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;
  art::Handle< std::vector<raw::RDTimeStamp> > RawTSListHandle;
  std::vector<art::Ptr<raw::RDTimeStamp> > TSlist;
  evt.getByLabel(fRawDigitTimeStampTag, RawTSListHandle);	
  art::fill_ptr_vector(TSlist, RawTSListHandle);

  timestamp = 0;
  timestamp2 = 0;
  apainconsist = 0;
  bool tsi = false;

  for(auto const & tsptr : TSlist) 
    {
      const raw::RDTimeStamp & rdts = *tsptr;
      uint16_t chan = rdts.GetFlags();
      int iapa = channelMap->APAFromOfflineChannel(chan);
      if (apas.find(iapa) != apas.end())
	{
	  ULong64_t ts = rdts.GetTimeStamp();
	  if (tsi)
	    {
	      if (ts != timestamp)
		{
		  apainconsist = iapa;
		  timestamp2 = ts;
		  return false;
		}
	    }
	  else
	    {
	      tsi = true;
	      timestamp = ts;
	      timestamp2 = ts;
	    }
	}
    }
  //std::cout << "All timestamps equal at: " << timestamp << " for APAs: ";
  //for (auto const &iter : apas)
  // {
  //  std::cout << iter << " "; 
  //}
  //std::cout << std::endl;
  return true;
}
