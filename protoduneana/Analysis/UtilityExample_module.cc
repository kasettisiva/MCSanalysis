////////////////////////////////////////////////////////////////////////
// Class:       UtilityExample
// Plugin Type: analyzer (art v2_07_03)
// File:        UtilityExample_module.cc
//
// Generated at Mon Sep  4 06:55:33 2017 by Leigh Whitehead using cetskelgen
// from cetlib version v3_00_01.
//
// This module is designed to show some usage examples of the analysis
// tools that I have been producing for protoDUNE. The aim has been to
// simplify the associations between the different objects to make
// some of the low-level art features more transparent
//
// The code is split into a few sections with each focusing on a different
// type of initial object
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "protoduneana/Analysis/ProtoDUNETrackUtils.h"
#include "protoduneana/Analysis/ProtoDUNEShowerUtils.h"
#include "protoduneana/Analysis/ProtoDUNETruthUtils.h"
#include "protoduneana/Analysis/ProtoDUNEPFParticleUtils.h"

namespace protoana {
  class UtilityExample;
}


class protoana::UtilityExample : public art::EDAnalyzer {
public:

  explicit UtilityExample(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UtilityExample(UtilityExample const &) = delete;
  UtilityExample(UtilityExample &&) = delete;
  UtilityExample & operator = (UtilityExample const &) = delete;
  UtilityExample & operator = (UtilityExample &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // fcl parameters
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  bool fVerbose;

};


protoana::UtilityExample::UtilityExample(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fVerbose(p.get<bool>("Verbose"))
{

}

void protoana::UtilityExample::beginJob()
{

}

void protoana::UtilityExample::analyze(art::Event const & evt)
{

  // We must have MC for this module to make sense
  if(evt.isRealData()) return;

  // ------------------------- Pure Truth Information -------------------------
  // One of the slightly annoying things about the protoDUNE simulation is that
  // LArG4 loses the information about which one of the true particles is the 
  // one we triggered on - also called the Good Particle. This is the particle
  // from the beam that we are trying to detect.

  // Get the truth utility to help us out
  protoana::ProtoDUNETruthUtils truthUtil;

  // Get the generator MCTruth objects and find the GEANT track id of the good particle. The generator MCTruth
  // object has the information about which particle is the good one. The GetGeantGoodParticle(...) function
  // takes care of giving you back the Geant4 version of this particle 

  // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
  // simulation has fGeneratorTag = "generator"
  auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
  // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
  // of these, so we pass the first element into the function to get the good particle
  const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
  if(geantGoodParticle != 0x0){
    std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() << std::endl;
  }
 
  // ------------------------- Reco Track Information -------------------------
  // In this section of code we will look at reconstructed tracks. The main
  // reconstruction we use is Pandora, and we will look at all of the tracks
  // reconstructed by this algorithm.

  // Get the track utility
  protoana::ProtoDUNETrackUtils trackUtil;
  // A few variables to tell us a bit about our reconstructed tracks
  unsigned int nTracksWithTruth = 0;
  unsigned int nTracksWithT0    = 0;
  unsigned int nTracksWithTag   = 0;
  unsigned int nTracksWithCalo  = 0;
  // Get the reconstructed tracks, by default from the "pandoraTrack" module
  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fTrackerTag);

  // Loop over the reconstructed tracks
  for(unsigned int t = 0; t < recoTracks->size(); ++t){

    // Take the t-th element for our current track
    const recob::Track thisTrack = (*recoTracks)[t];
    // We can then use the truth utility to give us the true track matching out reconstructed track
    const simb::MCParticle *trueMatch = truthUtil.GetMCParticleFromRecoTrack(thisTrack,evt,fTrackerTag);
    bool hasTruth = (trueMatch != 0x0);
 
    // Some of the tracks, primarily those that cross the cathode, will have a reconstructed T0.
    // We can access these through the track utility
    std::vector<anab::T0> trackT0 = trackUtil.GetRecoTrackT0(thisTrack,evt,fTrackerTag);
    bool hasT0 = (trackT0.size() != 0);

    // Similarly, some tracks can have a cosmic tag (though not from the pandora reconstruction)
    std::vector<anab::CosmicTag> trackCosmic = trackUtil.GetRecoTrackCosmicTag(thisTrack,evt,fTrackerTag);
    bool hasTag = (trackCosmic.size() != 0);

    // All tracks have some calorimetry information that stores things such as dE/dx etc.
    std::vector<anab::Calorimetry> trackCalo = trackUtil.GetRecoTrackCalorimetry(thisTrack,evt,fTrackerTag,fCalorimetryTag);
    bool hasCalo = (trackCalo.size() != 0);

    // Note that the above three functions return std::vectors of objects, but they will all either have
    // size 0 or 1. 

    if(hasTruth) ++nTracksWithTruth;
    if(hasT0)    ++nTracksWithT0;
    if(hasTag)   ++nTracksWithTag;
    if(hasCalo)  ++nTracksWithCalo;
  } // End loop over reconstructed tracks

  // Report back what we have managed to find out about the tracks
  std::cout << "Found " << recoTracks->size() << " reconstructed tracks:" << std::endl;
  std::cout << " - " << nTracksWithTruth << " successfully associated to the truth information " << std::endl;
  std::cout << " - " << nTracksWithT0    << " have a reconstructed T0" << std::endl;
  std::cout << " - " << nTracksWithTag   << " have a cosmic tag" << std::endl;
  std::cout << " - " << nTracksWithCalo  << " have calorimetry info" << std::endl;

  // ------------------------- Reco Shower Information ------------------------
  // There isn't a whole lot that we can do with shower objects. Here I just
  // get the average number of hits from all of the showers.

  // Get the shower utility
  protoana::ProtoDUNEShowerUtils showerUtil;

  // Get the reconstructed showers, by default from the "pandoraShower" module
  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(fShowerTag);

  float totalHits = 0.0;

  for(unsigned int s = 0; s < recoShowers->size(); ++s){
    // Get the shower
    const recob::Shower thisShower = (*recoShowers)[s];
    // Get the number of hits
    totalHits += showerUtil.GetNumberRecoShowerHits(thisShower,evt,fShowerTag);
  }

  std::cout << "Found " << recoShowers->size() << " showers with an average of " << totalHits/(float)recoShowers->size() << " hits" << std::endl;

  
  // ------------------------- PFParticle Information -------------------------
  // PFParticles are "particle flow" output objects that represent an 
  // individual particle. These are the highest level objects and contain
  // links to their daughter particles and are the recommended starting
  // point to use the output from the pandora reconstruction

  // Get the PFParticle utility
  protoana::ProtoDUNEPFParticleUtils pfpUtil;

  // Get all of the PFParticles, by default from the "pandora" product
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  // Get some information about the particles
  unsigned int nParticlesPrimary   = 0;
  unsigned int nParticlesWithT0    = 0;
  unsigned int nParticlesWithTag   = 0;
  for(unsigned int p = 0; p < recoParticles->size(); ++p){
    // Get the PFParticle
    const recob::PFParticle thisParticle = (*recoParticles)[p];
 
    // Only consider primary particles here
    if(!thisParticle.IsPrimary()) continue;
 
    ++nParticlesPrimary;

    // Does this particle have a T0?
    if(pfpUtil.GetPFParticleT0(thisParticle,evt,fPFParticleTag).size() != 0) ++nParticlesWithT0;

    // Does this particle have a cosmic tag?
    if(pfpUtil.GetPFParticleCosmicTag(thisParticle,evt,fPFParticleTag).size() != 0) ++nParticlesWithTag;

  }
  // Print out a summary
  std::cout << "Found " << nParticlesPrimary << " reconstructed primary particles:" << std::endl;
  std::cout << " - Total number of particles = " << recoParticles->size() << std::endl;
  std::cout << " - " << nParticlesWithT0    << " have a reconstructed T0" << std::endl;
  std::cout << " - " << nParticlesWithTag   << " have a cosmic tag" << std::endl;

  // Pandora divides up the detector into slices representing spacial and temporal regions
  // of the event. We can use metadata from pandora to build a slice map if we want to.
  // This consists of the slice number from pandora and a vector of primary PFParticles 
  // in that slice
  std::map<unsigned int, std::vector<const recob::PFParticle*>> sliceMap;
  // Get the slice map from the PFParticle utility
  sliceMap = pfpUtil.GetPFParticleSliceMap(evt,fPFParticleTag);
  unsigned int nMultiSlice = 0;
  // Let's count the number of slices containing multiple primary PFParticles
  std::cout << "Found " << sliceMap.size() << " slices with PFParticles" << std::endl;
  for(auto slice : sliceMap){
    if(slice.second.size() > 1) ++nMultiSlice;
  }
  std::cout << " - " << nMultiSlice << " have at least one primary PFParticle" << std::endl;

  // Pandora will also tag which slice contains the beam particle. The PFParticle utility
  // lets us find out which one
  unsigned short beamSlice = pfpUtil.GetBeamSlice(evt,fPFParticleTag);
  std::cout << "- Beam slice = " << beamSlice << std::endl;
  if(beamSlice != 9999){
    // We have found the beam slice. Let's get the particles contained in the slice
    std::vector<const recob::PFParticle*> beamSlicePrimaries = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
    std::cout << " - Found the beam slice! " << beamSlicePrimaries.size() << " beam particles in slice " << beamSlice << std::endl;
    // Loop over the particles in the beam slice
    for (const recob::PFParticle* prim : beamSlicePrimaries){
      // Use the PFParticle utility to get the vertex of the PFParticle
      // NB: this is the most upstream position for beam particles, and the highest point for cosmics
      const TVector3 vtx = pfpUtil.GetPFParticleVertex(*prim,evt,fPFParticleTag,fTrackerTag);
      std::cout << "Beam particle vertex: " << std::endl;
      vtx.Print();
      // If this is a track-like particle then we can get a secondary interaction point from the
      // downstream end of the track.
      if(pfpUtil.IsPFParticleTracklike(*prim,evt,fPFParticleTag,fTrackerTag)){
        const TVector3 sec = pfpUtil.GetPFParticleSecondaryVertex(*prim,evt,fPFParticleTag,fTrackerTag);
        std::cout << "Beam particle interaction vertex: " << std::endl;
        sec.Print();
        // The following(complicated) statement gets the recob::Track from the recob::PFParticle and then
        // gets the number of recob::Hit objects associated to the recob::Track
        std::cout << "PFParticle track has " << trackUtil.GetNumberRecoTrackHits(*(pfpUtil.GetPFParticleTrack(*prim,evt,fPFParticleTag,fTrackerTag)),evt,fTrackerTag) << " hits" << std::endl;
      }
      // Showers do not have a secondary vertex
      else{
        // The following(complicated) statement gets the recob::Shower from the recob::PFParticle and then
        // gets the number of recob::Hit objects associated to the recob::Shower
        std::cout << "PFParticle shower has " << showerUtil.GetNumberRecoShowerHits(*(pfpUtil.GetPFParticleShower(*prim,evt,fPFParticleTag,fShowerTag)),evt,fShowerTag) << " hits" << std::endl;
      }
      // We can also ask the PFParticle how many recob::Hit and recob::SpacePoint objects it contains
      std::cout << "Beam particle has " << pfpUtil.GetNumberPFParticleHits(*prim,evt,fPFParticleTag) << " hits and " 
                << pfpUtil.GetNumberPFParticleSpacePoints(*prim,evt,fPFParticleTag) << " space points" << std::endl;
      // Try to get the daughter tracks and showers
      const std::vector<const recob::Track*> daughterTracks = pfpUtil.GetPFParticleDaughterTracks(*prim,evt,fPFParticleTag,fTrackerTag);
      const std::vector<const recob::Shower*> daughterShowers = pfpUtil.GetPFParticleDaughterShowers(*prim,evt,fPFParticleTag,fShowerTag);
      std::cout << "Beam particle has " << daughterTracks.size() << " daughter tracks and " << daughterShowers.size() << " daughter showers." << std::endl;
    }

    // Look for unassociated hits in the beam slice
    const std::vector<const recob::Hit*> unassocHits = pfpUtil.GetPFParticleSliceUnassociatedHits(*(beamSlicePrimaries[0]),evt,fPFParticleTag);
    std::cout << "The besm slice has " << unassocHits.size() << " hits not associated to any PFParticle " << std::endl;
  }

  // The first pass of pandora reconstruction identifies clear cosmic particles. We can get a vector of these from the PFParticle utility
  std::cout << "There are " << pfpUtil.GetClearCosmicPFParticles(evt,fPFParticleTag).size() << " clear cosmic muons in this event" << std::endl;

}

void protoana::UtilityExample::endJob()
{

}

DEFINE_ART_MODULE(protoana::UtilityExample)

