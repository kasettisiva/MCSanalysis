// ProtoFlash_module
// Module that tries to associate T0s from reconstructed tracks and light flashes
// leigh.howard.whitehead@cern.ch
//

#ifndef ProtoFlash_H
#define ProtoFlash_H 1

// ROOT includes
#include "TH1D.h"
#include "TH2.h"
#include "TTree.h"

// C++ includes
#include <map>
#include <vector>
#include <iostream>

// LArSoft includes
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "dune/OpticalDetector/OpFlashSort.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

namespace protoana {
 
  class ProtoFlash : public art::EDAnalyzer{
  public:
 
    // Standard constructor and destructor for an ART module.
    ProtoFlash(const fhicl::ParameterSet&);
    virtual ~ProtoFlash();

    // This method is called once, at the start of the job. In this
    // example, it will define the histogram we'll write.
    void beginJob();

    // The analyzer routine, called once per event. 
    void analyze (const art::Event&); 

    // Extract true particle times and add them to a map
    void ExtractTrueTimes(const std::string &moduleName, std::map<int,double> &timeMap, const art::Event &evt);

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fOpFlashModuleLabel;       // Input tag for OpFlash collection
    std::string fOpHitModuleLabel;         // Input tag for OpHit collection
    std::string fSignalLabel;              // Input tag for the signal generator label
    std::string fCosmicLabel;              // Input tag for the cosmic generator label
    std::string fGeantLabel;               // Input tag for GEANT
    std::string fParticleLabel;            // Input tag for particle reconstruction
   
    TH1D *hFlashToTrueTime;
    TH1D *hFlashPurity;
    TH2D *hPurityTimeDiff;
    TH1D *hFlashToRecoTime;
    TH1D *hFlashToRecoWide;
    TH1D *hFlashToRecoWider;
    TH1D *hFlashTimes;
    TH1D *hHitTimes;
    TH1D *hNFlash;
    TH1D *hNHitPerFlash;
  };

} 

#endif // FlashMatchAna_H

namespace protoana {

  //-----------------------------------------------------------------------
  // Constructor
  ProtoFlash::ProtoFlash(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {

    // Indicate that the Input Module comes from .fcl
    fOpFlashModuleLabel = pset.get<std::string>("OpFlashModuleLabel");
    fOpHitModuleLabel   = pset.get<std::string>("OpHitModuleLabel");
    fSignalLabel        = pset.get<std::string>("SignalLabel");
    fCosmicLabel        = pset.get<std::string>("CosmicLabel");
    fGeantLabel         = pset.get<std::string>("GeantLabel");
    fParticleLabel      = pset.get<std::string>("ParticleLabel");

    art::ServiceHandle< art::TFileService > tfs;

    // Make a few plots
    hFlashToTrueTime    = tfs->make<TH1D>("hFlashToTrueTime","",50,-5,5);
    hFlashPurity     = tfs->make<TH1D>("hFlashPurity","",50,0.499,1.001);
    hFlashToRecoTime = tfs->make<TH1D>("hFlashToRecoTime","",50,-5,5);
    hFlashToRecoWide = tfs->make<TH1D>("hFlashToRecoWide","",100,-10,10);
    hFlashToRecoWider = tfs->make<TH1D>("hFlashToRecoWider","",100,-50,50);
    hNFlash          = tfs->make<TH1D>("hNFlash","",50,0,400);
    hNHitPerFlash    = tfs->make<TH1D>("hNHitPerFlash","",50,0,100);
    hFlashTimes      = tfs->make<TH1D>("hFlashTimes","",100,-4000,4000);
    hHitTimes        = tfs->make<TH1D>("hHitTimes","",100,-4000,4000);
    hPurityTimeDiff  = tfs->make<TH2D>("hPurityTimeDiff","",25,-5,5,25,0.499,1.001);
  }

  //-----------------------------------------------------------------------
  // Destructor
  ProtoFlash::~ProtoFlash() 
  {}
   
  //-----------------------------------------------------------------------
  void ProtoFlash::beginJob()
  {}

  //-----------------------------------------------------------------------
  void ProtoFlash::analyze(const art::Event& evt) 
  {

    // Make sure we can use this on data and MC
    bool isMC = !(evt.isRealData());
//    isMC = false;

    // Get flashes from event
    art::Handle< std::vector< recob::OpFlash > > FlashHandle;
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    if (evt.getByLabel(fOpFlashModuleLabel, FlashHandle)) {
      art::fill_ptr_vector(flashlist, FlashHandle);
      std::sort(flashlist.begin(), flashlist.end(), recob::OpFlashPtrSortByPE);
    }

    // Get assosciations between flashes and hits
    art::FindManyP< recob::OpHit > Assns(FlashHandle, evt, fOpFlashModuleLabel);

    // We want to store a map of track time to track id
    // This will later allow us to check how much of a flash comes from the best matched time
    std::map<int,double> trueTimeID;
    if(isMC){
      // Check for generator inputs - particle gun or protoDUNE beam
      if(fSignalLabel!=""){
        ExtractTrueTimes(fSignalLabel,trueTimeID,evt);
      }

      // Check for cosmic generator inputs
      if(fCosmicLabel!=""){
        ExtractTrueTimes(fCosmicLabel,trueTimeID,evt);
      }
    }

    //////////////////////////////////////
    // Access all the Flash Information //
    //////////////////////////////////////

    hNFlash->Fill(FlashHandle->size());

    // Use the clocks service to make sure to account for the offset between true times and the electronics clocks
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    // Store a map of flash times
    std::map<int,double> flashMap;

    // For every OpFlash in the vector
    for(unsigned int i = 0; i < FlashHandle->size(); ++i)
    {

      // Get OpFlash and associated hits
      art::Ptr< recob::OpFlash > TheFlashPtr(FlashHandle, i);
      recob::OpFlash TheFlash = *TheFlashPtr;
      art::FindManyP< recob::OpHit > Assns(FlashHandle, evt, fOpFlashModuleLabel);
      std::vector< art::Ptr<recob::OpHit> > matchedHits = Assns.at(i);

      // Account for the time offset in the TPC
      flashMap.insert(std::make_pair(i,TheFlash.Time() - clockData.TriggerOffsetTPC()));

      hNHitPerFlash->Fill(matchedHits.size());
      hFlashTimes->Fill(TheFlash.Time());

      for(auto const h : matchedHits){
        hHitTimes->Fill(h->PeakTime());
      }

      // Truth level code
      if(isMC){
        art::ServiceHandle<cheat::PhotonBackTrackerService> pbt;
    
        // Get the best match in time, then check how much of the light comes from
        // the best matched particle      
        int bestMatch = 0;
        double bestPurity = 0;
        double minTimeDiff = 1e20;
            
        // Iterate over the map
        for(auto &m : trueTimeID){
            
          double dist = fabs(m.second - TheFlash.Time());
          if(dist < minTimeDiff){
            // Does this track give us purity > 0.5?
            std::set<int> trackIDs;
            trackIDs.emplace(m.first);
            double purity = pbt->OpHitCollectionPurity(trackIDs,matchedHits);
            if(purity < 0.5) continue;
  
            // If purity is good enough, this is a success
            minTimeDiff = dist;
            bestMatch = m.first;
            bestPurity = purity;
          }
        }
        if(minTimeDiff > 1000){
  //        std::cout << "No sensible true match to flash " << i << std::endl;
        }
  
        if(bestPurity >= 0.5) std::cout << "Best true time matched particle " << bestMatch << " has delta T = " << minTimeDiff << " and gives a purity " << bestPurity << std::endl;
  
        hFlashToTrueTime->Fill(trueTimeID[bestMatch] - TheFlash.Time());
        hPurityTimeDiff->Fill(trueTimeID[bestMatch] - TheFlash.Time(),bestPurity);
        hFlashPurity->Fill(bestPurity);
      }
    } // End loop over flashes

    // Get reconstruction information from particles with a measured T0
    // Try finding some particles
    art::ValidHandle< std::vector<recob::PFParticle> > particleHandle
          = evt.getValidHandle<std::vector<recob::PFParticle> >(fParticleLabel);

    protoana::ProtoDUNEPFParticleUtils pfpUtil;

    for (size_t p = 0; p != particleHandle->size(); ++p){

      const auto thisParticle = (*particleHandle)[p];

      // Only consider primaries so we don't include deltas etc
      if(!(thisParticle.IsPrimary())) continue;

      // Did this particle have an associated T0?
      std::vector<anab::T0> t0s = pfpUtil.GetPFParticleT0(thisParticle,evt,fParticleLabel);
      if(t0s.size() == 0) continue;

      // Pandora gives us times in ns
      double recoT0 = t0s[0].Time() / 1000.; 
   
      int bestRecoMatch = 0;
      double minRecoTimeDiff = 1e20;

      // Now loop over the flashes
      for(auto &m : flashMap){
        double dist = fabs(recoT0 - m.second); 
        if(dist < minRecoTimeDiff){
          minRecoTimeDiff = dist;
          bestRecoMatch = m.first;
        }
      } 

      std::cout << "Particle " << p << " has a TPC T0 = " << recoT0 << std::endl;

      if(minRecoTimeDiff > 100){
        std::cout << "No sensible reco match to flash " << p << std::endl;
      }
      std::cout << "Best matched flash " << bestRecoMatch << " has delta T = " << minRecoTimeDiff << std::endl;

      hFlashToRecoTime->Fill(recoT0 - flashMap[bestRecoMatch]);
      hFlashToRecoWide->Fill(recoT0 - flashMap[bestRecoMatch]);
      hFlashToRecoWider->Fill(recoT0 - flashMap[bestRecoMatch]);
    }

  }

  void ProtoFlash::ExtractTrueTimes(const std::string &moduleName, std::map<int,double> &timeMap, const art::Event &evt){


    auto MCTruthsHandle = evt.getValidHandle<std::vector<simb::MCTruth> >(moduleName);
    art::Ptr<simb::MCTruth> thisMCTruth(MCTruthsHandle, 0);
    if (thisMCTruth->NParticles() == 0) {
      mf::LogError("FlashMatchAna") << "No Cosmic MCTruth Particles";
    }
    else std::cout << "MCTruth from " << moduleName  << " has " << thisMCTruth->NParticles() << " particles " << std::endl;

    // Get all the track ids associated with the cosmic events.
    art::FindManyP<simb::MCParticle> geantAssns(MCTruthsHandle,evt,fGeantLabel);

    for ( size_t i = 0; i < geantAssns.size(); i++) {
      auto parts = geantAssns.at(i);
      int counter = 0;
      for (auto part = parts.begin(); part != parts.end(); part++) {
        ++counter;
        timeMap.insert(std::make_pair((*part)->TrackId(),(*part)->T()/1000.));
      }
      std::cout << "Found " << counter << " GEANT particles associated to " << moduleName << " MCTruth " << i << std::endl;
    }

  }

} // namespace protoana

namespace protoana {
  DEFINE_ART_MODULE(ProtoFlash)
}
