#include "protoduneana/Analysis/ProtoDUNESliceUtils.h"

#include "lardataobj/RecoBase/Hit.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

protoana::ProtoDUNESliceUtils::ProtoDUNESliceUtils(){

}

protoana::ProtoDUNESliceUtils::~ProtoDUNESliceUtils(){

}

// Get the hits from a given reco slice
const std::vector<const recob::Hit*> protoana::ProtoDUNESliceUtils::GetRecoSliceHits(const recob::Slice &slice, art::Event const &evt, const std::string sliceModule) const{

  return GetRecoSliceHits(slice.ID(),evt,sliceModule);  

}

// Get the reco hits but using the slice id instead
const std::vector<const recob::Hit*> protoana::ProtoDUNESliceUtils::GetRecoSliceHits(const unsigned int sliceID, art::Event const &evt, const std::string sliceModule) const{

  auto recoSlices = evt.getValidHandle<std::vector<recob::Slice> >(sliceModule);
  art::FindManyP<recob::Hit> findHits(recoSlices,evt,sliceModule);
  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(sliceID);

  std::vector<const recob::Hit*> sliceHits;

  for(const art::Ptr<recob::Hit> hit : inputHits){

    sliceHits.push_back(hit.get());

  }

  return sliceHits;

}

// Get the number of hits from a given reco slice
unsigned int protoana::ProtoDUNESliceUtils::GetNumberRecoSliceHits(const recob::Slice &slice, art::Event const &evt, const std::string sliceModule) const{

  return GetRecoSliceHits(slice,evt,sliceModule).size();

}

// Get the number of hits from a given slice id
unsigned int protoana::ProtoDUNESliceUtils::GetNumberRecoSliceHits(const unsigned int sliceID, art::Event const &evt, const std::string sliceModule) const{

  return GetRecoSliceHits(sliceID,evt,sliceModule).size();

}

// Get a map of a slice number and all hits in the slice
const std::map<unsigned int, std::vector<const recob::Hit*>> protoana::ProtoDUNESliceUtils::GetRecoSliceHitMap(art::Event const &evt, const std::string sliceModule) const{

  auto recoSlices = evt.getValidHandle<std::vector<recob::Slice> >(sliceModule);
  std::map<unsigned int, std::vector<const recob::Hit*>> hitMap;

  for(auto const slice : *recoSlices){

    const std::vector<const recob::Hit*> constvec = GetRecoSliceHits(slice.ID(),evt,sliceModule);
    for(auto const h : constvec){
      hitMap[slice.ID()].push_back(h);
    }

  }

  return hitMap;

}

