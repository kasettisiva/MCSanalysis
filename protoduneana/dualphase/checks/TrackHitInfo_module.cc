////////////////////////////////////////////////////////////////////////
// Class:       CosmicsdQdx
// Plugin Type: analyzer (art v3_05_01)
// File:        CosmicsdQdx_module.cc
//
// Generate simple summary tree from track dQ/dx measurements for each
// collection view to check track calorimetry and effecit CRP gains
// Follow the same logic as in TrackCalorimetryAlg
//
// Generated at Tue May 26 16:28:38 2020 by Vyacheslav Galymov using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <utility>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" 
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/TrackUtils.h" 

//
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"


//
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::pair;

//
namespace pddpana {
  class TrackHitInfo;
}


//
class pddpana::TrackHitInfo : public art::EDAnalyzer {
public:
  explicit TrackHitInfo(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackHitInfo(TrackHitInfo const&) = delete;
  TrackHitInfo(TrackHitInfo&&) = delete;
  TrackHitInfo& operator=(TrackHitInfo const&) = delete;
  TrackHitInfo& operator=(TrackHitInfo&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  int      fLogLevel;
  string   fTrackModuleLabel;
  // track utils
  protoana::ProtoDUNETrackUtils trackUtil;

  // detector geometry
  const geo::Geometry* fGeom;
  // detector properties
  //const detinfo::DetectorProperties* fDetprop;
};


pddpana::TrackHitInfo::TrackHitInfo(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
  fLogLevel( p.get< int >("LogLevel") ),
  fTrackModuleLabel( p.get< std::string  >("TrackModuleLabel") )
  {
    fGeom    = &*art::ServiceHandle<geo::Geometry>();
    //fDetprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  }

//
void pddpana::TrackHitInfo::analyze(art::Event const& e)
{
  const string myname = "pddpana::TrackHitInfo::analyze: ";

  // get hits ValidHandle< std::vector<recob::Hits> >
  auto Tracks  = e.getValidHandle<vector<recob::Track>>(fTrackModuleLabel);
  
  if( fLogLevel >= 2 ){
    cout<<myname<<"The event contains "<< Tracks->size() <<" tracks\n";
  }

  unsigned EventId = e.id().event();
  
  unsigned countbad = 0;

  // loop over tracks
  for (unsigned itrk = 0; itrk < Tracks->size(); ++itrk) {
    const recob::Track& track = Tracks->at(itrk);
    unsigned TrackId   = itrk;
    unsigned TrajPnts  = track.NumberTrajectoryPoints();
    float    TrackLen  = track.Length();
          
    if( fLogLevel >= 3 ){
      cout<<myname<<"Track ID "<<TrackId<<" has "
	  <<TrajPnts<<" points and is "<<TrackLen<<" cm long"
	  << "\n  start at: ( " << track.Vertex().X()
	  << " ; " << track.Vertex().Y()
	  << " ; " << track.Vertex().Z()
	  << "\n  end at:   ( " << track.End().X()
	  << " ; " << track.End().Y()
	  << " ; " << track.End().Z()<<" )"<<endl;
    }

    vector<unsigned> hitsTpcId( fGeom->Nplanes() );

    // loop over the planes 
    for(size_t i_plane = 0; i_plane<fGeom->Nplanes(); i_plane++) {
      
      // get hits in this plane
      auto hits = trackUtil.GetRecoTrackHitsFromPlane( track, e, fTrackModuleLabel, i_plane );
      if( fLogLevel >= 3 ){
	cout<<myname<<"Hits in plane "<<i_plane<<" "<<hits.size()<<endl;
      }

      vector<unsigned> plane_hits_tpcid( hits.size() );
      // loop over hits
      for(size_t i_hit = 0; i_hit<hits.size(); i_hit++ ){
	plane_hits_tpcid[i_hit] = hits[i_hit]->WireID().TPC;
      }
      
      int this_tpcid = -1;
      auto start = plane_hits_tpcid.begin();
      auto end   = plane_hits_tpcid.end();
      if( std::equal(start + 1, end, start)) {
	this_tpcid = (int)(*start);
      } //
      else {
	if( fLogLevel>= 2){
	  cout<<myname<<"hits for the same plane have mixed TPC IDs. Skipping...\n";
	}
      }
      if( this_tpcid < 0 ){
	hitsTpcId.clear();
	break;
      }
      else{
	hitsTpcId[i_plane] = (unsigned)this_tpcid;
      }
    }// end plane loop
    
    if( hitsTpcId.empty() ) continue;
        
    //if( fLogLevel >= 1 ){
    auto start = hitsTpcId.begin();
    auto end   = hitsTpcId.end();
    if( !std::equal(start + 1, end, start)) {
      cout<<myname<<"mismatch in TPC ID for the initial hits: ";
      cout<<" event ID "<<EventId<<"; hit TPC IDs ";
      for( auto const &v: hitsTpcId ){ cout<<v<<" ";}
      cout<<endl;
      countbad++;
    } //

  }// end track loop
  
  //
  cout<<"Found "<<countbad<<" CRP mismatched tracks in total "<<Tracks->size()<<endl;
  
} // end analyze()

//
void pddpana::TrackHitInfo::beginJob()
{;}

//
void pddpana::TrackHitInfo::endJob()
{;}

DEFINE_ART_MODULE(pddpana::TrackHitInfo)
