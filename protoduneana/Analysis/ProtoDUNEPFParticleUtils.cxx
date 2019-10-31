#include "protoduneana/Analysis/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Analysis/ProtoDUNESliceUtils.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

protoana::ProtoDUNEPFParticleUtils::ProtoDUNEPFParticleUtils(){

}

protoana::ProtoDUNEPFParticleUtils::~ProtoDUNEPFParticleUtils(){

}

// Get the number of primary PFParticles
unsigned int protoana::ProtoDUNEPFParticleUtils::GetNumberPrimaryPFParticle(art::Event const &evt, const std::string particleLabel) const{

  // Get the particles
  auto pfParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);

  unsigned int nPrimary = 0;
  for(unsigned int p = 0; p < pfParticles->size(); ++p){
    if(pfParticles->at(p).IsPrimary()){
      ++nPrimary;
    }
  }

  return nPrimary;
}

// Helper to get slice maps and avoid duplicate code
const std::map<unsigned int,std::vector<const recob::PFParticle*>> protoana::ProtoDUNEPFParticleUtils::SliceMapHelper(art::Event const &evt, const std::string particleLabel, bool primaryOnly) const{

  // Get the particles
  auto pfParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);

  std::map<unsigned int, std::vector<const recob::PFParticle*>> sliceMap;

  for(unsigned int p = 0; p < pfParticles->size(); ++p){
    const recob::PFParticle* particle = &(pfParticles->at(p));

    //  Only the primary particles have the slice association
    if(primaryOnly && !particle->IsPrimary()) continue;

    unsigned int thisSlice = GetPFParticleSliceIndex(*particle,evt,particleLabel);

    if(thisSlice != 9999){
      sliceMap[thisSlice].push_back(particle);
    }
  }

  return sliceMap;

}

// Return a map of primary particles grouped by their reconstructed slice. Useful for finding slices with multiple particles
const std::map<unsigned int,std::vector<const recob::PFParticle*>> protoana::ProtoDUNEPFParticleUtils::GetPFParticleSliceMap(art::Event const &evt, const std::string particleLabel) const{

  return SliceMapHelper(evt,particleLabel,true);

}

// Return a map of all particles grouped by their reconstructed slice. 
const std::map<unsigned int,std::vector<const recob::PFParticle*>> protoana::ProtoDUNEPFParticleUtils::GetAllPFParticleSliceMap(art::Event const &evt, const std::string particleLabel) const{

  return SliceMapHelper(evt,particleLabel,false);

}

// Get the cosmic tag(s) from a given PFParticle
std::vector<anab::CosmicTag> protoana::ProtoDUNEPFParticleUtils::GetPFParticleCosmicTag(const recob::PFParticle &particle, art::Event const &evt, std::string particleLabel) const{

  try{
    auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle> >(particleLabel);

    unsigned int pIndex = particle.Self();

    std::vector<anab::CosmicTag> pTags;

    const art::FindManyP<anab::CosmicTag> findCosmicTags(recoParticles,evt,particleLabel);
    for(unsigned int p = 0; p < findCosmicTags.at(pIndex).size(); ++p){
      pTags.push_back((*(findCosmicTags.at(pIndex)[p])));
    }

    return pTags;
  }
  catch(...){
    return std::vector<anab::CosmicTag>();
  }
}

// Get the T0(s) from a given PFParticle
std::vector<anab::T0> protoana::ProtoDUNEPFParticleUtils::GetPFParticleT0(const recob::PFParticle &particle, art::Event const &evt, std::string particleLabel) const{

  return GetPFParticleT0(particle,evt,particleLabel,particleLabel);

}

std::vector<anab::T0> protoana::ProtoDUNEPFParticleUtils::GetPFParticleT0(const recob::PFParticle &particle, art::Event const &evt, std::string particleLabel, std::string t0Label) const{

  try{
    auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle> >(particleLabel);

    unsigned int pIndex = particle.Self();

    std::vector<anab::T0> pT0s;

    const art::FindManyP<anab::T0> findParticleT0s(recoParticles,evt,t0Label);
    for(unsigned int p = 0; p < findParticleT0s.at(pIndex).size(); ++p){
      pT0s.push_back((*(findParticleT0s.at(pIndex)[p])));
    }

    return pT0s;
  }
  catch(...){
    return std::vector<anab::T0>();
  }

}

// Try to get the slice tagged as beam
unsigned short protoana::ProtoDUNEPFParticleUtils::GetBeamSlice(art::Event const &evt, const std::string particleLabel) const{

  const std::map<unsigned int, std::vector<const recob::PFParticle*>> sliceMap = GetPFParticleSliceMap(evt,particleLabel);

  for(auto slice : sliceMap){
    for(auto particle : slice.second){
      if(IsBeamParticle(*particle,evt,particleLabel)){
        return slice.first;
      }
    }
  }

  return 9999;

}

// Returns pointers for the primary PFParticles in a slice
const std::vector<const recob::PFParticle*> protoana::ProtoDUNEPFParticleUtils::GetPFParticlesFromSlice(const unsigned short slice, art::Event const &evt, const std::string particleLabel) const{

  const std::map<unsigned int, std::vector<const recob::PFParticle*>> sliceMap = GetPFParticleSliceMap(evt,particleLabel);
  
  if(sliceMap.find(slice) != sliceMap.end()){
    return sliceMap.at(slice);
  }
  else{
    return std::vector<const recob::PFParticle*>();
  }

}

// Returns pointers for all PFParticles in a slice
const std::vector<const recob::PFParticle*> protoana::ProtoDUNEPFParticleUtils::GetAllPFParticlesFromSlice(const unsigned short slice, art::Event const &evt, const std::string particleLabel) const{

  const std::map<unsigned int, std::vector<const recob::PFParticle*>> sliceMap = GetAllPFParticleSliceMap(evt,particleLabel);

  if(sliceMap.find(slice) != sliceMap.end()){
    return sliceMap.at(slice);
  }
  else{
    return std::vector<const recob::PFParticle*>();
  }

}

// Return the pointers for the PFParticles in the beam slice
const std::vector<const recob::PFParticle*> protoana::ProtoDUNEPFParticleUtils::GetPFParticlesFromBeamSlice(art::Event const &evt, const std::string particleLabel) const{

  unsigned short beamSlice = GetBeamSlice(evt,particleLabel);
  return GetPFParticlesFromSlice(beamSlice,evt,particleLabel);
}

// Access the BDT output used to decide if a slice is beam-like or cosmic-like
float protoana::ProtoDUNEPFParticleUtils::GetBeamCosmicScore(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  return FindFloatInMetaData(particle,evt,particleLabel,"TestBeamScore");

}

// Use the pandora metadata to tell us if this is a beam particle or not
bool protoana::ProtoDUNEPFParticleUtils::IsBeamParticle(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  return FindBoolInMetaData(particle,evt,particleLabel,"IsTestBeam");

}

bool protoana::ProtoDUNEPFParticleUtils::FindBoolInMetaData(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string entry) const{

  std::map<std::string,float> mdMap = GetPFParticleMetaData(particle,evt,particleLabel);
  
  if(mdMap.find(entry) != mdMap.end()){
    return true;
  }
  else{
    return false;
  }

}

float protoana::ProtoDUNEPFParticleUtils::FindFloatInMetaData(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string entry) const{

  std::map<std::string,float> mdMap = GetPFParticleMetaData(particle,evt,particleLabel);
  
  if(mdMap.find(entry) != mdMap.end()){
    return mdMap.at(entry);
  }
  else{
    return -999.;
  }

}


// Get the reconstructed slice associated with a particle
const recob::Slice* protoana::ProtoDUNEPFParticleUtils::GetPFParticleSlice(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  // Perhaps we should use the associations to do this? 
  auto pfParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);
  const art::FindOneP<recob::Slice> findSlice(pfParticles,evt,particleLabel);

  const recob::Slice* slice = findSlice.at(particle.Self()).get();

  return slice;
}

// Get the reconstructed slice associated with a particle
unsigned short protoana::ProtoDUNEPFParticleUtils::GetPFParticleSliceIndex(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  // Try to use slices if we can
  try{
    const recob::Slice* slice = GetPFParticleSlice(particle,evt,particleLabel);
    return slice->ID();
  }
  // Otherwise fall back on metadata
  catch(...){
    std::map<std::string,float> mdMap = GetPFParticleMetaData(particle,evt,particleLabel);
    std::string search = "SliceIndex";
    if(mdMap.find(search) != mdMap.end()){
      return static_cast<unsigned short>(mdMap.at(search));
    }
    else{
//    std::cerr << "Object has no slice index... returning 9999" << std::endl;
      return 9999;
    }
  }

}

// Get all hits from a PFParticle's slice
const std::vector<const recob::Hit*> protoana::ProtoDUNEPFParticleUtils::GetPFParticleSliceHits(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  // Use the slice utility to help us out 
  protoana::ProtoDUNESliceUtils sliceUtil;
  
  unsigned short sliceIndex = GetPFParticleSliceIndex(particle,evt,particleLabel);
  const std::vector<const recob::Hit*> hits = sliceUtil.GetRecoSliceHits(sliceIndex,evt,particleLabel);

  return hits;
}

// Get hits from a PFParticle's slice that are not associated to any PFParticle
const std::vector<const recob::Hit*> protoana::ProtoDUNEPFParticleUtils::GetPFParticleSliceUnassociatedHits(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  unsigned short sliceNumber = GetPFParticleSliceIndex(particle,evt,particleLabel);
  const std::vector<const recob::Hit*> hits = GetPFParticleSliceHits(particle,evt,particleLabel);

  std::vector<bool> hitUsed(hits.size(),false);
  std::vector<const recob::Hit*> unassocHits;

  // We need to iterate through the PFParticles in this slice and get their hits
  const std::vector<const recob::PFParticle*> slicePFPs = GetAllPFParticlesFromSlice(sliceNumber,evt,particleLabel);
 
  for(auto const p : slicePFPs){

    // We are in the correct slice, so get the hits
    const std::vector<const recob::Hit*> pfpHits = GetPFParticleHits(*p,evt,particleLabel);

    std::cout << "Hit vector sizes: " << hits.size() << ", " << pfpHits.size() << std::endl;

    for(const recob::Hit* h : pfpHits){
      // Loop through the hit collection from the slice
      for(unsigned int o = 0; o < hits.size(); ++o){
        if(fabs(h->Integral() - hits[o]->Integral()) < 1e-5 && fabs(h->PeakTime() - hits[o]->PeakTime()) < 1e-5){
          hitUsed[o] = true;
          break;
        }
      }
    }
  }

  // Now we can fill our unassociated hits vector
  for(unsigned int i = 0; i < hits.size(); ++i){
    if(!hitUsed[i]) unassocHits.push_back(hits[i]);
  }
  std::cout << " Number of unassociated hits = " << unassocHits.size() << " of " << hits.size() << std::endl;

  return unassocHits;
}

const std::map<std::string,float> protoana::ProtoDUNEPFParticleUtils::GetPFParticleMetaData(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const {

  // Get the particles
  auto pfParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);
  // And their meta data
  const art::FindManyP<larpandoraobj::PFParticleMetadata> findMetaData(pfParticles,evt,particleLabel);

  const larpandoraobj::PFParticleMetadata metaData = *((findMetaData.at(particle.Self())).at(0));

  return metaData.GetPropertiesMap();
}

// Pandora tags and removes clear cosmics before slicing, so check if this particle is a clear cosmic
bool protoana::ProtoDUNEPFParticleUtils::IsClearCosmic(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  return FindBoolInMetaData(particle,evt,particleLabel,"IsClearCosmic");

}

// Get all of the clear cosmic ray particles
const std::vector<const recob::PFParticle*> protoana::ProtoDUNEPFParticleUtils::GetClearCosmicPFParticles(art::Event const &evt, const std::string particleLabel) const{

  // Get the particles
  auto pfParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);

  std::vector<const recob::PFParticle*> cosmicParticles;

  for(unsigned int p = 0; p < pfParticles->size(); ++p){
    const recob::PFParticle* particle = &(pfParticles->at(p));

    //  Only consider primary particles
    if(!particle->IsPrimary()) continue;

    if(IsClearCosmic(*particle,evt,particleLabel)){
      cosmicParticles.push_back(particle);
    }

  }

  return cosmicParticles;

}

// Function to find the interaction vertex of a primary PFParticle
const TVector3 protoana::ProtoDUNEPFParticleUtils::GetPFParticleVertex(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string trackLabel) const{

  // Pandora produces associations between PFParticles and recob::Vertex objects
  auto pfParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);
  const art::FindManyP<recob::Vertex> findVertices(pfParticles,evt,particleLabel);
  const std::vector<art::Ptr<recob::Vertex>> vertices = findVertices.at(particle.Self());

  // What happens next depends on the type of event.
  // Not primary    -> just use the pfparticle vertex
  // Shower objects -> just use the pfparticle vertex
  // Cosmics        -> use track start point
  // Beam           -> use track start point

//  std::cout << "PFParticle daughters " << particle.NumDaughters() << std::endl;

  // Non-primary particle or shower-like primary particle
  if(!particle.IsPrimary() || !IsPFParticleTracklike(particle,evt,particleLabel,trackLabel)){
    if(vertices.size() != 0){
      const recob::Vertex* vtx = (vertices.at(0)).get();
      return TVector3(vtx->position().X(),vtx->position().Y(),vtx->position().Z());
    }
    else{
      std::cerr << "Non track-like PFParticle has no vertex?! Return default vector" << std::endl;
      return TVector3();
    }
  }
  else{
    // Cosmic or track-like beam primary particle
  
    const recob::Track* track = GetPFParticleTrack(particle,evt,particleLabel,trackLabel);

    if(track != 0x0){
      const TVector3 start(track->Trajectory().Start().X(),track->Trajectory().Start().Y(),track->Trajectory().Start().Z());
      const TVector3 end(track->Trajectory().End().X(),track->Trajectory().End().Y(),track->Trajectory().End().Z());
      // Return the most upstream point as some cases where the track is reversed...
      if(IsBeamParticle(particle,evt,particleLabel)){ 
        if(start.Z() < end.Z()) return start;
        else return end;
      }
      // Return the highest point for cosmics
      else{
        if(start.Y() > end.Y()) return start;
        else return end;
      }
    }
    else{
      std::cerr << "No track found for track-like PFParticle?! Return default vector" << std::endl;
      return TVector3();
    }

  }

}

// Function to find the secondary interaction vertex of a primary PFParticle
const TVector3 protoana::ProtoDUNEPFParticleUtils::GetPFParticleSecondaryVertex(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string trackLabel) const{

  // In this case we want to find the end of the track-like PFParticle
  // To do this, we need to access things via the track

  // Showers have no secondary vertex
  if(!IsPFParticleTracklike(particle,evt,particleLabel,trackLabel)){
    std::cerr << "Shower-like PFParticles have no secondary vertex. Returning default TVector3" << std::endl;
    return TVector3();
  }

  const recob::Track* track = GetPFParticleTrack(particle,evt,particleLabel,trackLabel);

  // Interaction vertex is the downstream end of the track, or the bottom for cosmics
  if(track != 0x0){
    const TVector3 start(track->Trajectory().Start().X(),track->Trajectory().Start().Y(),track->Trajectory().Start().Z());
    const TVector3 end(track->Trajectory().End().X(),track->Trajectory().End().Y(),track->Trajectory().End().Z());

      // Return the most downstream point as some cases where the track is reversed...
      if(IsBeamParticle(particle,evt,particleLabel)){
        if(start.Z() > end.Z()) return start;
        else return end;
      }
      // Return the lowest point for cosmics
      else{
        if(start.Y() < end.Y()) return start;
        else return end;
      }

  }
  else{
    std::cerr << "This track-like PFParticle has no associated track. Returning default TVector3" << std::endl;
    return TVector3();
  }

}

// Is the particle track-like?
bool protoana::ProtoDUNEPFParticleUtils::IsPFParticleTracklike(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string trackLabel) const{

  return (GetPFParticleTrack(particle,evt,particleLabel,trackLabel) != 0x0);

}

// This is the old and deprecated version of the function 
bool protoana::ProtoDUNEPFParticleUtils::IsPFParticleTracklike(const recob::PFParticle &particle) const{
  std::cerr << "This function is deprecated, please use the version with four arguments." << std::endl;
  if(abs(particle.PdgCode()) == 11){
    return false;
  }
  else{
    return true;
  }
}



bool protoana::ProtoDUNEPFParticleUtils::IsPFParticleShowerlike(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string showerLabel) const{

  return (GetPFParticleShower(particle,evt,particleLabel,showerLabel) != 0x0);

}

// This is the old and deprecated version of the function
bool protoana::ProtoDUNEPFParticleUtils::IsPFParticleShowerlike(const recob::PFParticle &particle) const{
  return !IsPFParticleTracklike(particle);
}

// Get the track associated to this particle. Returns a null pointer if not found.
const recob::Track* protoana::ProtoDUNEPFParticleUtils::GetPFParticleTrack(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string trackLabel) const{

  // Pandora produces associations between PFParticles and recob::Track objects
  auto particles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);
  const art::FindManyP<recob::Track> findTracks(particles,evt,trackLabel);
  const std::vector<art::Ptr<recob::Track>> pfpTracks = findTracks.at(particle.Self());

  // Check that the track exists
  if(pfpTracks.size() != 0){
    const recob::Track* track = (pfpTracks.at(0)).get();
    return track;
  }
  else{
//    std::cerr << "No track found, returning null pointer" << std::endl;
    return nullptr;
  }

}

// Get the shower associated to this particle. Returns a null pointer if not found.
const recob::Shower* protoana::ProtoDUNEPFParticleUtils::GetPFParticleShower(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel, const std::string showerLabel) const{

  // Pandora produces associations between PFParticles and recob::Track objects
  auto particles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);
  const art::FindManyP<recob::Shower> findShowers(particles,evt,showerLabel);
  const std::vector<art::Ptr<recob::Shower>> pfpShowers = findShowers.at(particle.Self());

  // Check that the shower exists
  if(pfpShowers.size() != 0){
    const recob::Shower* shw = (pfpShowers.at(0)).get();
    return shw;
  }
  else{
//    std::cerr << "No shower found, returning null pointer" << std::endl;
    return nullptr;
  }

}

// Get the space points associated to the PFParticle
const std::vector<const recob::SpacePoint*> protoana::ProtoDUNEPFParticleUtils::GetPFParticleSpacePoints(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  // Get the particles and their associations
  auto particles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);
  const art::FindManyP<recob::SpacePoint> findSpacePoints(particles,evt,particleLabel);
  const std::vector<art::Ptr<recob::SpacePoint>> pfpSpacePoints = findSpacePoints.at(particle.Self());

  // We don't want the art::Ptr so we need to get rid of it
  std::vector<const recob::SpacePoint*> sp;
  for(auto pointer : pfpSpacePoints){
    sp.push_back(pointer.get());
  }  

  return sp;
}

// Get the number of space points
unsigned int protoana::ProtoDUNEPFParticleUtils::GetNumberPFParticleSpacePoints(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  return GetPFParticleSpacePoints(particle,evt,particleLabel).size();

}

// Get the clusters associated to the PFParticle
const std::vector<const recob::Cluster*> protoana::ProtoDUNEPFParticleUtils::GetPFParticleClusters(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  // Get the particles and their associations
  auto particles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);
  const art::FindManyP<recob::Cluster> findClusters(particles,evt,particleLabel);
  const std::vector<art::Ptr<recob::Cluster>> pfpClusters = findClusters.at(particle.Self());

  // We don't want the art::Ptr so we need to get rid of it
  std::vector<const recob::Cluster*> clusters;
  for(auto pointer : pfpClusters){
    clusters.push_back(pointer.get());
  }  

  return clusters;
}

// Get the number of clusters associated to the PFParticle
unsigned int protoana::ProtoDUNEPFParticleUtils::GetNumberPFParticleClusters(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  return GetPFParticleClusters(particle,evt,particleLabel).size();

}

// Get the hits associated to the PFParticle
const std::vector<const recob::Hit*> protoana::ProtoDUNEPFParticleUtils::GetPFParticleHits(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  const std::vector<const recob::Cluster*> pfpClusters = GetPFParticleClusters(particle,evt,particleLabel);
  auto allClusters = evt.getValidHandle<std::vector<recob::Cluster>>(particleLabel);
  const art::FindManyP<recob::Hit> findHits(allClusters,evt,particleLabel);

  std::vector<const recob::Hit*> pfpHits;
 
  // Store all of the hits in a single vector 
  for(auto cluster : pfpClusters){
    const std::vector<art::Ptr<recob::Hit>> clusterHits = findHits.at(cluster->ID());
//    std::cout << "Cluster has " << clusterHits.size() << " hits" << std::endl;
    for(auto hit : clusterHits){
      pfpHits.push_back(hit.get());
    }
  }

  return pfpHits;
}

const std::vector< art::Ptr< recob::Hit > > protoana::ProtoDUNEPFParticleUtils::GetPFParticleHits_Ptrs(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  const std::vector<const recob::Cluster*> pfpClusters = GetPFParticleClusters(particle,evt,particleLabel);
  auto allClusters = evt.getValidHandle<std::vector<recob::Cluster>>(particleLabel);
  const art::FindManyP<recob::Hit> findHits(allClusters,evt,particleLabel);

  std::vector< art::Ptr< recob::Hit > > pfpHits;
 
  // Store all of the hits in a single vector 
  for(auto cluster : pfpClusters){
    const std::vector<art::Ptr<recob::Hit>> clusterHits = findHits.at(cluster->ID());
    for(auto hit : clusterHits){
      pfpHits.push_back(hit);
    }
  }

  return pfpHits;
}

// Get the number of hits
unsigned int protoana::ProtoDUNEPFParticleUtils::GetNumberPFParticleHits(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  return GetPFParticleHits(particle,evt,particleLabel).size();

}

// Get the total charge from all PFParticle hits
const std::vector<double> protoana::ProtoDUNEPFParticleUtils::GetPFParticleHitsCharge(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const {

  const std::vector<const recob::Hit*> hitvector = GetPFParticleHits(particle, evt, particleLabel);

  double hitcharge[3] = {0.0, 0.0, 0.0};
  for(auto hit : hitvector){
    int hit_plane = (int)hit->WireID().Plane;
    if(hit_plane < 0) continue;
    if(hit_plane > 2) continue;
    //hitcharge[hit_plane] += hit->SummedADC();
    hitcharge[hit_plane] += hit->Integral();
  }
  
  std::vector<double> hitchargevec;
  for(int i = 0; i < 3; i++)
    hitchargevec.push_back(hitcharge[i]);

  return hitchargevec;
}

// Get the earliest hit peak time
const double protoana::ProtoDUNEPFParticleUtils::GetPFParticleEarliestHitPeakTime(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const {

  const std::vector<const recob::Hit*> hitvector = GetPFParticleHits(particle, evt, particleLabel);
  
  double earliesthittime = 999999.0;
  for(auto hit : hitvector){
    if(hit->PeakTime() < earliesthittime)
      earliesthittime = (double)hit->PeakTime();
  }
  
  return earliesthittime;
}

// Get the daughter tracks from the PFParticle
const std::vector<const recob::Track*> protoana::ProtoDUNEPFParticleUtils::GetPFParticleDaughterTracks(const recob::PFParticle &particle, art::Event const &evt, 
                                                                                                       const std::string particleLabel, const std::string trackLabel) const{

  // Get the PFParticles
  auto particles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);

  std::vector<const recob::Track*> daughterTracks;

  // Loop over the daughters
  for(size_t daughterID : particle.Daughters()){
    const recob::Track* track = GetPFParticleTrack(particles->at(daughterID),evt,particleLabel,trackLabel);
    if(track != 0x0){
      daughterTracks.push_back(track);
    }  
  }

  return daughterTracks;

}

// Get the daughter showers from the PFParticle
const std::vector<const recob::Shower*> protoana::ProtoDUNEPFParticleUtils::GetPFParticleDaughterShowers(const recob::PFParticle &particle, art::Event const &evt, 
                                                                                                         const std::string particleLabel, const std::string showerLabel) const{

  // Get the PFParticles
  auto particles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);

  std::vector<const recob::Shower*> daughterShowers;

  // Loop over the daughters
  for(size_t daughterID : particle.Daughters()){
    const recob::Shower* shower = GetPFParticleShower(particles->at(daughterID),evt,particleLabel,showerLabel);
    if(shower != 0x0){
      daughterShowers.push_back(shower);
    }
  }

  return daughterShowers;

}




