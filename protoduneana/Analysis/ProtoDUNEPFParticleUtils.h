#ifndef PROTODUNE_PFPARTICLE_UTILS_H
#define PROTODUNE_PFPARTICLE_UTILS_H

///////////////////////////////////////////////////////////////////
// ProtoDUNEPFParticleUtils
//  - Class to help analysers access useful PFParticle information
// 
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////////

#include <map>
#include <string>

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "art/Framework/Principal/Event.h"

#include "TVector3.h"

namespace protoana {

  class ProtoDUNEPFParticleUtils {

  public:

    ProtoDUNEPFParticleUtils();
    ~ProtoDUNEPFParticleUtils();

    /// Get the number of primary PFParticles
    unsigned int GetNumberPrimaryPFParticle(art::Event const &evt, const std::string particleLabel) const;

    /// Get a map of slice index to the primary PFParticles within it
    const std::map<unsigned int,std::vector<const recob::PFParticle*>> GetPFParticleSliceMap(art::Event const &evt, const std::string particleLabel) const;

    /// Get a map of slice index to all of the PFParticles within it
    const std::map<unsigned int,std::vector<const recob::PFParticle*>> GetAllPFParticleSliceMap(art::Event const &evt, const std::string particleLabel) const;

    /// Get the Primary PFParticles from a given slice. Returns an empty vector if the slice number is not valid
    const std::vector<const recob::PFParticle*> GetPFParticlesFromSlice(const unsigned short slice, art::Event const &evt, const std::string particleLabel) const;

    /// Get all of the PFParticles from a given slice. Returns an empty vector if the slice number is not valid
    const std::vector<const recob::PFParticle*> GetAllPFParticlesFromSlice(const unsigned short slice, art::Event const &evt, const std::string particleLabel) const;

    /// Try to get the slice tagged as beam. Returns 9999 if no beam slice was found
    unsigned short GetBeamSlice(art::Event const &evt, const std::string particleLabel) const;

    /// Return the pointers for the PFParticles in the beam slice. Returns an empty vector is no beam slice was found
    const std::vector<const recob::PFParticle*> GetPFParticlesFromBeamSlice(art::Event const &evt, const std::string particleLabel) const;

    /// Get the reconstructed slice associated with a particle
    const recob::Slice* GetPFParticleSlice(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get the reconstructed slice number associated with a particle
    unsigned short GetPFParticleSliceIndex(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// For a given PFParticle, return all hits from the slice it comes from
    const std::vector<const recob::Hit*> GetPFParticleSliceHits(const recob::PFParticle &particlei, art::Event const &evt, const std::string particleLabel) const;

    /// For a given PFParticle find its slice and return all those hits not associated to any PFParticle
    const std::vector<const recob::Hit*> GetPFParticleSliceUnassociatedHits(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get the cosmic tag(s) from a given PFParticle
    std::vector<anab::CosmicTag> GetPFParticleCosmicTag(const recob::PFParticle &particle, art::Event const &evt, std::string particleLabel) const;

    /// Get the T0(s) from a given PFParticle
    std::vector<anab::T0> GetPFParticleT0(const recob::PFParticle &particle, art::Event const &evt, std::string particleLabel) const;
    std::vector<anab::T0> GetPFParticleT0(const recob::PFParticle &particle, art::Event const &evt, std::string particleLabel, std::string t0Label) const;

    /// Access the BDT output used to decide if a slice is beam-like or cosmic-like
    float GetBeamCosmicScore(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Use the pandora metadata to tell us if this is a beam particle or not
    bool IsBeamParticle(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Pandora tags and removes clear cosmics before slicing, so check if this particle is a clear cosmic
    bool IsClearCosmic(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get all of the clear cosmic ray particles
    const std::vector<const recob::PFParticle*> GetClearCosmicPFParticles(art::Event const &evt, const std::string particleLabel) const;

    /// Get the metadata associated to a PFParticle from pandora
    const std::map<std::string,float> GetPFParticleMetaData(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Function to find the interaction vertex of a primary PFParticle
    const TVector3 GetPFParticleVertex(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string trackLabel) const;

    /// Function to find the secondary interaction vertex of a primary PFParticle
    const TVector3 GetPFParticleSecondaryVertex(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string trackLabel) const;

    /// Is the particle track-like?
    bool IsPFParticleTracklike(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string trackLabel) const;
    // Old deprecated version
    bool IsPFParticleTracklike(const recob::PFParticle &particle) const;

    /// Is the particle track-like?
    bool IsPFParticleShowerlike(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string showerLabel) const;
    // Old deprecated version
    bool IsPFParticleShowerlike(const recob::PFParticle &particle) const;

    /// Get the track associated to this particle. Returns a null pointer if not found.
    const recob::Track* GetPFParticleTrack(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string trackLabel) const;

    /// Get the shower associated to this particle. Returns a null pointer if not found.
    const recob::Shower* GetPFParticleShower(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string showerLabel) const;

    // Get the SpacePoints associated to the PFParticle
    const std::vector<const recob::SpacePoint*> GetPFParticleSpacePoints(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get the number of space points
    unsigned int GetNumberPFParticleSpacePoints(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const; 

    /// Get the clusters associated to the PFParticle
    const std::vector<const recob::Cluster*> GetPFParticleClusters(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get the number of clusters associated to the PFParticle
    unsigned int GetNumberPFParticleClusters(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get the hits associated to the PFParticle
    const std::vector<const recob::Hit*> GetPFParticleHits(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;
    const std::vector< art::Ptr< recob::Hit > > GetPFParticleHits_Ptrs(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get the number of hits
    unsigned int GetNumberPFParticleHits(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const; 

    /// Get the total hit charge for each view
    const std::vector<double> GetPFParticleHitsCharge(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get the earliest hit peak time
    const double GetPFParticleEarliestHitPeakTime(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get the daughter tracks from the PFParticle
    const std::vector<const recob::Track*> GetPFParticleDaughterTracks(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string trackLabel) const;

    /// Get the daughter showers from the PFParticle
    const std::vector<const recob::Shower*> GetPFParticleDaughterShowers(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string showerLabel) const;

  private:

    // Some functions that don't need to be exposed to the outside
    
    /// Helper to get the slice map and avoid code repetition
    const std::map<unsigned int,std::vector<const recob::PFParticle*>> SliceMapHelper(art::Event const &evt, const std::string particleLabel, bool primaryOnly) const;

    /// Look for entries in the meta data
    bool FindBoolInMetaData(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string entry) const;
    float FindFloatInMetaData(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string entry) const;

  };

}

#endif

