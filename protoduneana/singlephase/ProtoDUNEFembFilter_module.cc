// A module to filter out events with inactive FEMBs
// owen.goodwin@postgrad.manchester.ac.uk

#include <iostream>
#include <utility>
#include <set>

#include "TH1.h"
#include "TFile.h"

#include "art/Framework/Core/EDAnalyzer.h" 
#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art_root_io/TFileService.h"

#include "protoduneana/Analysis/ProtoDUNEDataUtils.h"
#include "lardataobj/RawData/RDTimeStamp.h"




namespace filt{

  class ProtoDUNEFembFilter : public art::EDFilter {
  public:
    explicit ProtoDUNEFembFilter(fhicl::ParameterSet const & pset);
    virtual ~ProtoDUNEFembFilter() {};
    virtual bool filter(art::Event& e);
    void beginJob();

  private:
    protoana::ProtoDUNEDataUtils fDataUtils;

    unsigned int fLogLevel;
    bool fRequireBeamsideFembsOnly;
    bool fRequireBeamsideTimestampConsistencyOnly;
    
    TH1D* fSelectedEvents;
    TH1D* fTotalEvents;

  };

  ProtoDUNEFembFilter::ProtoDUNEFembFilter::ProtoDUNEFembFilter(fhicl::ParameterSet const & pset):
    EDFilter(pset), fDataUtils(pset.get<fhicl::ParameterSet>("DataUtils"))
  {

    fLogLevel = pset.get<unsigned int>("LogLevel");
    fRequireBeamsideFembsOnly = pset.get<bool>("RequireBeamsideFembsOnly");
    fRequireBeamsideTimestampConsistencyOnly = pset.get<bool>("RequireBeamsideTimestampConsistencyOnly");
    if ( fLogLevel >= 1 ) {
      std::cout  << "                LogLevel: " << fLogLevel << std::endl;
      if(fRequireBeamsideFembsOnly){
        std::cout  << "Filtering events with inactive FEMBs on beamside APAs"<< std::endl;
      }    
      else{
        std::cout  << "Filtering events with any inactive FEMBs"<< std::endl;
      }
      if(fRequireBeamsideTimestampConsistencyOnly){
        std::cout  << "Filtering events with inconsistent timestmaps on beamside APAs"<< std::endl;
      }    
      else{
        std::cout  << "Filtering events with any inconsistent timestamps"<< std::endl;
      }
    }
    
  }

  void ProtoDUNEFembFilter::beginJob() {

    art::ServiceHandle<art::TFileService> tfs;
      fSelectedEvents = tfs->make<TH1D>("fSelectedEvents", "Number of Selected Events", 3, 0, 3); //counts the number of selected events 
      fTotalEvents = tfs->make<TH1D>("fTotalEvents", "Total Events", 3, 0, 3); //counts the initial number of events in the unfiltered root input file
    
  }

  bool ProtoDUNEFembFilter::filter(art::Event & evt) {


    fTotalEvents->Fill(1); //count total events

    
    if(!evt.isRealData()){
        fSelectedEvents->Fill(1); 
        return true;   //Filter is designed for Data only. Don't want to filter on MC
      }




    std::vector<int> BeamsideAPAs;
    std::vector<int> AllAPAs;
    std::vector<int> checkedAPAs;
    std::vector<int> TScheckedAPAs;

    // Add some elements to myIntVector
    BeamsideAPAs.push_back(0);
    BeamsideAPAs.push_back(2);
    BeamsideAPAs.push_back(4);

    AllAPAs.push_back(0);
    AllAPAs.push_back(1);
    AllAPAs.push_back(2);
    AllAPAs.push_back(3);
    AllAPAs.push_back(4);
    AllAPAs.push_back(5);

    const std::string myname = "ProtoDUNEFembFilter::filter: ";

    if(fRequireBeamsideFembsOnly){
      checkedAPAs=BeamsideAPAs;
    }
    else{
      checkedAPAs=AllAPAs;
    }
    if(fRequireBeamsideTimestampConsistencyOnly){
      TScheckedAPAs=BeamsideAPAs;
    }
    else{
      TScheckedAPAs=AllAPAs;
    }

    // make a set out of these for faster lookup by the timestamp checker

    std::set<int> checkedAPAset;
    for (size_t i=0; i < checkedAPAs.size(); ++i)
      {
	checkedAPAset.emplace(TScheckedAPAs.at(i));
      }
    
    bool keep = true;
    // Helper utility functions

    fTotalEvents->Fill(1); //count total events
    for (auto APA = checkedAPAs.begin(); APA != checkedAPAs.end(); ++APA){ //loop through beam side APAs
      //std::cout<<"APA:"<<*APA<<std::endl;
      //std::cout<<fDataUtils.GetNActiveFembsForAPA(evt, *APA)<<std::endl;
      if (fDataUtils.GetNActiveFembsForAPA(evt, *APA)!=20){ //check if APA has all 20 fembs active

        if (fLogLevel >=2) std::cout<<"Missing FEMBs on APA: "<<*APA<<std::endl; 
        keep=false; //if not remove event
      }      
    }

    // check timestamp consistency

    ULong64_t timestamp=0;
    ULong64_t timestamp2=0;
    int apainconsist=0;
    if (!fDataUtils.CheckTimeStampConsistencyForAPAs(evt, checkedAPAset, timestamp, timestamp2, apainconsist ))
      {
	keep = false;
        if (fLogLevel >=2) std::cout<<"ProtoDUNEFembFilter Timestamp mismatch: " << timestamp << " vs " << timestamp2 << " on TPC set " << apainconsist  << std::endl; 
      }

    if ( fLogLevel >=2 ) std::cout << myname << (keep ? "Keep" : "Reject") << "ing event." << std::endl;
    if (keep==true) fSelectedEvents->Fill(1); //count total events
    
    return keep;

  }

  DEFINE_ART_MODULE(ProtoDUNEFembFilter)

}
