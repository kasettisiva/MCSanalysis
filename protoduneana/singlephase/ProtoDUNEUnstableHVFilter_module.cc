////////////////////////////////////////////////////////////////////////
//
// ProtoDUNEUnstableHVFilter class
// 
// author: Owen Goodwin
// email: owen.goodwin@manchester.ac.uk
//
// - A filter to reject events between given start and end times of unstable HV periods. Using raw decoder timestamp.
// 
//    - All dates and times should be in unix time (UTC)
//
////////////////////////////////////////////////////////////////////////
// C++
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
// ROOT
#include "TMath.h"
#include "TTimeStamp.h"

#include "TH1.h"
#include "TFile.h"
/// Framework 
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNETimeStamp.h"
///filters for events, etc
namespace filter {
        class ProtoDUNEUnstableHVFilter : public art::EDFilter  {
        public:
                explicit ProtoDUNEUnstableHVFilter(fhicl::ParameterSet const& ); 
                virtual ~ProtoDUNEUnstableHVFilter();      
                uint64_t GetRawDecoderInfo(art::Event & e);
                bool filter(art::Event& evt);
                void beginJob();
                void endJob();
                void reconfigure(fhicl::ParameterSet const& p);
        private:

                std::vector<std::pair<UInt_t,UInt_t>> fTimeRanges; 
                UInt_t fTimeRangeLow;
                UInt_t fTimeRangeHigh;
                bool fDebug;
                long long RDTSTime;
                double RDTSTimeSec;
                int RDTSTrigger;
                art::Handle< std::vector<raw::RDTimeStamp> > RDTimeStampHandle;


                TH1D* fSelectedEvents;
                TH1D* fTotalEvents;


        }; //class ProtoDUNEUnstableHVFilter
}
void filter::ProtoDUNEUnstableHVFilter::beginJob() { 

    art::ServiceHandle<art::TFileService> tfs;
    fSelectedEvents = tfs->make<TH1D>("fSelectedEvents", "Number of Selected Events", 3, 0, 3); //counts the number of selected events 
    fTotalEvents = tfs->make<TH1D>("fTotalEvents", "Total Events", 3, 0, 3); //counts the initial number of events in the unfiltered root input file
}
void filter::ProtoDUNEUnstableHVFilter::endJob() { }


void filter::ProtoDUNEUnstableHVFilter::reconfigure(fhicl::ParameterSet const& p) {
        fTimeRanges  = p.get<std::vector <std::pair<UInt_t,UInt_t >>>("TimeRanges");
        fDebug = p.get<int>("Debug");

}  


filter::ProtoDUNEUnstableHVFilter::ProtoDUNEUnstableHVFilter(fhicl::ParameterSet const& pset)
: EDFilter(pset) {
        this->reconfigure(pset);
}
filter::ProtoDUNEUnstableHVFilter::~ProtoDUNEUnstableHVFilter() { }


uint64_t filter::ProtoDUNEUnstableHVFilter::GetRawDecoderInfo(art::Event & e){
    MF_LOG_INFO("BeamEvent") << "\n";
    MF_LOG_INFO("BeamEvent") << "Getting Raw Decoder Info" << "\n";
    e.getByLabel("timingrawdecoder","daq",RDTimeStampHandle);
    MF_LOG_INFO("BeamEvent") << "RDTS valid? " << RDTimeStampHandle.isValid() << "\n";
    for (auto const & RDTS : *RDTimeStampHandle){
        MF_LOG_INFO("BeamEvent") << "High: " << RDTS.GetTimeStamp_High() << "\n";
        MF_LOG_INFO("BeamEvent") << "Low: " << RDTS.GetTimeStamp_Low() << "\n";

        uint64_t high = RDTS.GetTimeStamp_High();
        uint64_t low  = RDTS.GetTimeStamp_Low();

        high = high << 32;
        uint64_t joined = (high | low);

        MF_LOG_INFO("BeamEvent") << "Raw Decoder Timestamp: " << joined << "\n";

        RDTSTime = joined;
        RDTSTrigger = RDTS.GetFlags();
        MF_LOG_INFO("BeamEvent") << "Trigger: " << RDTSTrigger << "\n";

        //Separates seconds portion of the ticks 
        //From the nanoseconds
        long long RDTSTickSec = (RDTSTime * 2) / (int)(TMath::Power(10,8));
        RDTSTickSec = RDTSTickSec * (int)(TMath::Power(10,8)) / 2;
        //long long RDTSTickNano = RDTSTime - RDTSTickSec;

        //Units are 20 nanoseconds ticks
        RDTSTimeSec  = 20.e-9 * RDTSTickSec;
        //RDTSTimeNano = 20.    * RDTSTickNano;


  }
  return RDTSTimeSec;
}

bool filter::ProtoDUNEUnstableHVFilter::filter(art::Event &evt) {   


        fTotalEvents->Fill(1);



        if(!evt.isRealData()){
            fSelectedEvents->Fill(1); 
            return true;   //Filter is designed for Data only. Don't want to filter on MC
            }



        const std::string myname = "ProtoDUNEUnstableHVFilter::filter: ";
        bool keep = true;
        TTimeStamp * evtTTS;
        evtTTS = new TTimeStamp(GetRawDecoderInfo(evt));
        // if (evtTime.timeHigh() == 0) { evtTTS = new TTimeStamp(evtTime.timeLow()); }
        // else { evtTTS = new TTimeStamp(evtTime.timeHigh(), evtTime.timeLow()); }
        if (fDebug) std::cout << "Event time:  " << evtTTS -> AsString() << std::endl;
        // Requested time range lower end

        for (auto TimeRange : fTimeRanges){ //loop through unstable hv time ranges

            fTimeRangeLow=TimeRange.first;
            fTimeRangeHigh=TimeRange.second;


                        // Check that input time is in correct format
            
            
            // Check that input times are in correct format
            if (fTimeRangeHigh < 1000000000 || fTimeRangeLow < 1000000000 ) {
                    std::cout << "Warning: please provide time in POSIX foramt, event time "
                              << "filter returning false." << std::endl; 
                    return false;

            }



            if (fTimeRangeHigh < fTimeRangeLow ) {
                    std::cout << "Warning: Lower limit bigger than lower limit "
                              << "filter returning false." << std::endl; 
                    return false;

            }
            //no idea what this is supposed to be doing!! will reject any times before 10am
            // 
            // if (fTimeRangeHigh > 0 && fTimeRangeHigh < 100000) {
            //         std::cout << "Warning: please provide time in format HHMMSS, event time "
            //                   << "filter returning false.2" << std::endl; 
            //         return false;
            // }
            // if (fTimeRangeLow > 0 && fTimeRangeLow < 100000) {
            //         std::cout << "Warning: please provide time in format HHMMSS, event time "
            //                   << "filter returning false.3" << std::endl; 
            //         return false;
            // }
            // Event time
            //art::Timestamp evtTime = evt.time();


                // all the checking he dies then ask if its in the correct range

            TTimeStamp * ttsLow(nullptr); 
            
        
            ttsLow = new TTimeStamp(fTimeRangeLow); 
            // Requested time range higher end
            TTimeStamp * ttsHigh(nullptr);

            ttsHigh = new TTimeStamp(fTimeRangeHigh); 
  
            // Filter decision
            if(fDebug){
                std::cout << "Lower Limit:  " << ttsLow -> AsString() << std::endl;
                std::cout << "Upper Limit:  " << ttsHigh -> AsString() << std::endl;
            }
    
                
            if (evtTTS -> GetSec() > ttsLow -> GetSec() && 
                evtTTS -> GetSec() < ttsHigh -> GetSec()) { keep=false; }
                    

    }
    if ( fDebug ) std::cout << myname << (keep ? "Keep" : "Reject") << "ing event." << std::endl;
    if (keep==true) fSelectedEvents->Fill(1); //count total events
    //keep=false; //for testing
    return keep;    //need this to stop the get to end error make sure to check this is working as intended
}
 DEFINE_ART_MODULE(filter::ProtoDUNEUnstableHVFilter) 
