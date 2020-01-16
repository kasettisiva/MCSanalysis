// A very simple module to filter out the events with beam triggers
// leigh.howard.whitehead@cern.ch

#include <iostream>
#include <utility>
#include <set>

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 

#include "protoduneana/Utilities/ProtoDUNEDataUtils.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"

namespace filt{

  class ProtoDUNETriggerFilter : public art::EDFilter {
  public:
    explicit ProtoDUNETriggerFilter(fhicl::ParameterSet const & pset);
    virtual ~ProtoDUNETriggerFilter() {};
    virtual bool filter(art::Event& e);

  private:

    unsigned int fLogLevel;
    std::vector<unsigned int> fTimingFlagSelectList;
    std::vector<unsigned int> fTimingFlagDeselectList;
    std::string fTimingLabel;
    std::string fTimingInstance;
    std::string fTriggerLabel;
    std::string fTriggerInstance;

  };

  ProtoDUNETriggerFilter::ProtoDUNETriggerFilter::ProtoDUNETriggerFilter(fhicl::ParameterSet const & pset)
: EDFilter(pset) {
    using std::cout;
    using std::endl;
    const std::string myname = "ProtoDUNETriggerFilter::ctor: ";
    fLogLevel = pset.get<unsigned int>("LogLevel");
    std::vector<unsigned int> defaulttriglist;
    defaulttriglist.push_back(0xc);
    fTimingFlagSelectList = pset.get<std::vector<unsigned int> >("TimingFlagSelectList",defaulttriglist);
    std::vector<unsigned int> emptylist;
    fTimingFlagDeselectList = pset.get<std::vector<unsigned int> >("TimingFlagDeselectList",emptylist);

    fTimingLabel = pset.get<std::string>("TimingLabel","timingrawdecoder");
    fTimingInstance = pset.get<std::string>("TimingInstance","daq");
    fTriggerLabel = pset.get<std::string>("TriggerLabel","ctbrawdecoder");
    fTriggerInstance = pset.get<std::string>("TriggerInstance","daq");
    if ( fLogLevel >= 1 ) {
      cout << myname << "                LogLevel: " << fLogLevel << endl;
      cout << myname << "    TimingFlagSelectList: [";
      bool first = true;
      for ( unsigned int flg : fTimingFlagSelectList ) {
        if ( first ) first = false;
        else cout << ", ";
        cout << flg;
      }
      cout << "]" << endl;
      cout << myname << "  TimingFlagDeselectList: [";
      first = true;
      for ( unsigned int flg : fTimingFlagDeselectList ) {
        if ( first ) first = false;
        else cout << ", ";
        cout << flg;
      }
      cout << "]" << endl;
      cout << myname << "             TimingLabel: " << fTimingLabel << endl;
      cout << myname << "          TimingInstance: " << fTimingInstance << endl;
      cout << myname << "            TriggerLabel: " << fTriggerLabel << endl;
      cout << myname << "         TriggerInstance: " << fTriggerInstance << endl;
    }
  }

  bool ProtoDUNETriggerFilter::filter(art::Event & evt) {
    using std::cout;
    using std::endl;
    const std::string myname = "ProtoDUNETriggerFilter::filter: ";

    bool keep = true;

    bool checkTriggerFlag = fTimingFlagSelectList.size() || fTimingFlagDeselectList.size();

    std::string stinfo = "Trigger check disabled.";
    if ( keep && checkTriggerFlag ) {
      // Fetch the trigger and timing clock.
      art::Handle<std::vector<raw::RDTimeStamp>> htims;
      evt.getByLabel(fTimingLabel, fTimingInstance, htims);
      //art::Handle<std::vector<raw::ctb::pdspctb> > hctb;
      //evt.getByLabel(fTriggerLabel, fTriggerInstance, hctb);

      if ( ! htims.isValid() ) {
        std::cout << myname << "WARNING: Timing clocks product not found." << std::endl;
        if ( fLogLevel >=2 ) stinfo = "Timing clocks product not found.";
      } else if (  htims->size() != 1 ) {
        std::cout << myname << "WARNING: Unexpected timing clocks size: " << htims->size() << std::endl;
        if ( fLogLevel >=2 ) stinfo = "Unexpected timing clocks size.";
        for ( unsigned int itim=0; itim<htims->size() && itim<50; ++itim ) {
          std::cout << myname << "  " << htims->at(itim).GetTimeStamp() << std::endl;
        }
      } else {
        const raw::RDTimeStamp& tim = htims->at(0);

        // See https://twiki.cern.ch/twiki/bin/view/CENF/TimingSystemAdvancedOp#Reference_info
        unsigned int trigFlag = tim.GetFlags();
        if ( fLogLevel >=2 ) stinfo = "Trigger flag: " + std::to_string(trigFlag);

        // If TimingFlagSelectList has entries, the trigger flag must be there.
        if ( fTimingFlagSelectList.size() ) {
          keep = false;
          for ( unsigned int flg : fTimingFlagSelectList ) {
            if ( keep ) break;
            if ( flg == trigFlag) keep = true;
          }
        }
  
        // The trigger flag must not be in TimingFlagDeselectList.
        for ( unsigned int flg : fTimingFlagDeselectList ) {
          if ( ! keep ) break;
          if ( flg == trigFlag ) keep = false;
        }
  
      }
    }

    if ( fLogLevel >=2 ) std::cout << myname << (keep ? "Keep" : "Reject") << "ing event " << evt.event()
                                   << ". "  << stinfo << endl;
    return keep;

  }

  DEFINE_ART_MODULE(ProtoDUNETriggerFilter)

}
