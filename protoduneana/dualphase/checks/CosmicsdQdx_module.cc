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

// ROOT
#include "TTree.h"
//#include "TMath.h"
//#include "TH1F.h"
//#include "TFile.h"

//
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::pair;

//
namespace pddpana {
  class CosmicsdQdx;
}

//
namespace {
  //from calo::TrackCalorimetryAlg
  class dist_projected{
  public:
    dist_projected(recob::Hit const* h, geo::Geometry const* g):
      hit(h), geom(g){}
    bool operator() (std::pair<geo::WireID,float> i, std::pair<geo::WireID,float> j)
    {
      float dw_i = ((int)(i.first.Wire) - (int)(hit->WireID().Wire))*geom->WirePitch(i.first.Plane);
      float dw_j = ((int)(j.first.Wire) - (int)(hit->WireID().Wire))*geom->WirePitch(j.first.Plane);
      float dt_i = i.second - hit->PeakTime();
      float dt_j = j.second - hit->PeakTime();
      return (std::sqrt(dw_i*dw_i + dt_i*dt_i) < std::sqrt(dw_j*dw_j + dt_j*dt_j));
    }
  private:
    recob::Hit const* hit;
    geo::Geometry const* geom;
  };
}


class pddpana::CosmicsdQdx : public art::EDAnalyzer {
public:
  explicit CosmicsdQdx(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicsdQdx(CosmicsdQdx const&) = delete;
  CosmicsdQdx(CosmicsdQdx&&) = delete;
  CosmicsdQdx& operator=(CosmicsdQdx const&) = delete;
  CosmicsdQdx& operator=(CosmicsdQdx&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  int      fLogLevel;
  string   fTrackModuleLabel;
  float    fTrackMinLen;
  float    fTrackDriftCut;
  float    fTrackWallCut;
  unsigned fMaxHitMultiplicity;

  unsigned fDrift;
  
  // summary tree
  TTree   *fTree;

  //
  unsigned fEventNum;
  unsigned fTrackId;
  unsigned fTrajPoints;
  float    fTrackLen;
  float    fDriftOffset;
  //float    fThetaXZ;
  
  //
  unsigned fPlane;
  float    fHitAdcSum;
  float    fHitIntegral;
  float    fHitPeak;
  float    fHitTime;
  float    fPitch;
  float    fdQdx;
  float    fX;
  float    fY;
  float    fZ;
 
  // track utils
  protoana::ProtoDUNETrackUtils trackUtil;

  // detector geometry
  const geo::Geometry* fGeom;
  // detector properties
  const detinfo::DetectorProperties* fDetprop;

  bool checkCutsAndGetT0( const recob::Track& track, float &T0 );
};


pddpana::CosmicsdQdx::CosmicsdQdx(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
  fLogLevel( p.get< int >("LogLevel") ),
  fTrackModuleLabel( p.get< std::string  >("TrackModuleLabel") ),
  fTrackMinLen( p.get< float  >("TrackMinLen") ),
  fTrackDriftCut( p.get< float  >("TrackDriftCut") ),
  fTrackWallCut( p.get< float  >("TrackWallCut") ),
  fMaxHitMultiplicity( p.get< float  >("MaxHitMultiplicity") )
  {
    fGeom    = &*art::ServiceHandle<geo::Geometry>();
    fDetprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    //auto const &tpc = ;
    fDrift = (std::abs( fGeom->TPC(0).DetectDriftDirection() ));
    if( fDrift != 1 ){
      throw cet::exception("CosmicsdQdx") << "Drift direction "<<fDrift
					  << " is not support\n";
    }
  }

//
void pddpana::CosmicsdQdx::analyze(art::Event const& e)
{
  const string myname = "pddpana::CosmicsdQdx::analyze: ";

  // get hits ValidHandle< std::vector<recob::Hits> >
  auto Tracks  = e.getValidHandle<vector<recob::Track>>(fTrackModuleLabel);
  
  if( fLogLevel >= 2 ){
    cout<<myname<<"The event contains "<< Tracks->size() <<" tracks\n";
  }

  fEventNum = e.id().event();

  // loop over tracks
  for (unsigned itrk = 0; itrk < Tracks->size(); ++itrk) {
    const recob::Track& track = Tracks->at(itrk);
    fTrackId     = itrk;
    fTrajPoints  = track.NumberTrajectoryPoints();
    fTrackLen    = track.Length();
    fDriftOffset = 0;
    if( !checkCutsAndGetT0( track, fDriftOffset) ){
	if( fLogLevel >= 3 ){
	  cout<<myname<<"Skipping track ID "<<fTrackId<<" has "
	      <<fTrajPoints<<" points and is "<<fTrackLen<<" cm long\n";
	}	continue;
    }
      
      
    if( fLogLevel >= 3 ){
      cout<<myname<<"Track ID "<<fTrackId<<" has "
	  <<fTrajPoints<<" points and is "<<fTrackLen<<" cm long"
	  << "\n  start at: ( " << track.Vertex().X()
	  << " ; " << track.Vertex().Y()
	  << " ; " << track.Vertex().Z()
	  << "\n  end at:   ( " << track.End().X()
	  << " ; " << track.End().Y()
	  << " ; " << track.End().Z()<<" )"<<endl;
    }
      
    //loop over the planes 
    for(size_t i_plane=0; i_plane<fGeom->Nplanes(); i_plane++) {
      fPlane = i_plane;
      
      // get hits in this plane
      auto hits = trackUtil.GetRecoTrackHitsFromPlane( track, e, fTrackModuleLabel, i_plane );
      if( fLogLevel >= 3 ){
	cout<<myname<<"Hits in plane "<<i_plane<<" "<<hits.size()<<endl;
      }
      
      //project down the track into wire/tick space for this plane
      vector< unsigned > traj_points_idx;
      vector< pair<geo::WireID,float> > traj_points_in_plane; //(track.NumberTrajectoryPoints());
      for(size_t i_trjpt=0; i_trjpt<track.NumberTrajectoryPoints(); i_trjpt++){
	auto pnt        = track.LocationAtPoint(i_trjpt);
	if( !track.HasValidPoint( i_trjpt ) ){
	  if( fLogLevel>=3 ){
	    cerr<<"track point "<<i_trjpt<<" has position ("
		<<pnt.X()<<", "<<pnt.Y()<<", "<<pnt.Z()<<") \n"
		<<" and is not a valid point "<<track.HasValidPoint( i_trjpt )<<endl;
	  }
	  continue;
	}

	double x_pos    = pnt.X();
	float tick      = fDetprop->ConvertXToTicks(x_pos,(int)i_plane,0,0);
	auto tpcid      = fGeom->PositionToTPCID(pnt);
	if( tpcid.TPC >= fGeom->NTPC() ){
	  if( fLogLevel>=3 ){
	    cerr<<myname<<"tpc no "<<tpcid.TPC<<" is not valid \n"
		<<" track id "<<fTrackId<<" at ("
		<<pnt.X()<<", "<<pnt.Y()<<", "<<pnt.Z()<<") \n"
		<<" has valid point "<<track.HasValidPoint( i_trjpt )<<endl;
	  }
	  continue;
	}
	geo::WireID wid = fGeom->NearestWireID(pnt, geo::PlaneID(tpcid, i_plane));
	//traj_points_in_plane[i_trjpt] = std::make_pair(wid, tick);
	traj_points_in_plane.push_back( std::make_pair(wid, tick) );
	traj_points_idx.push_back( i_trjpt );
      }
      
      // from calo::TrackCalorimetryAlg::AnalyzeHit
      for(auto const &hit: hits ){
	//skip high mulitplicity hits
	if(hit->Multiplicity() > (int)fMaxHitMultiplicity) continue;
	//
	size_t traj_iter = std::distance( traj_points_in_plane.begin(),
					  std::min_element( traj_points_in_plane.begin(), 
							    traj_points_in_plane.end(),
							    dist_projected(hit, fGeom) ) );
	size_t traj_pnt_idx = traj_points_idx[traj_iter];
	try{
	  fPitch = lar::util::TrackPitchInView(track, fGeom->View(hit->WireID().Plane), traj_pnt_idx);
	}
	catch(cet::exception &e){
	  if( fLogLevel>=3 ){
	    cout<<e;
	  }
	  continue;
	}
	fHitAdcSum   = hit->SummedADC();
	fHitIntegral = hit->Integral();
	fHitPeak     = hit->PeakAmplitude();
	fHitTime     = hit->PeakTime();
	fdQdx        = fHitAdcSum / fPitch;
	auto pnt     = track.LocationAtPoint(traj_iter);
	fX = pnt.X();
	fY = pnt.Y();
	fZ = pnt.Z();

	fTree->Fill();
      }// loop over track plane hits

    }// end plane loop
    
  }// end track loop

} // end analyze()

//
void pddpana::CosmicsdQdx::beginJob()
{
  // init summary tree
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("cranaTree","Check cosmics dQdx");
  fTree->Branch("EventNum", &fEventNum, "EventNum/i");
  fTree->Branch("TrackId", &fTrackId, "TrackId/i");
  fTree->Branch("TrajPoints", &fTrajPoints, "TrajPoints/i");
  fTree->Branch("TrackLen",   &fTrackLen,   "TrackLen/F");
  fTree->Branch("DriftOffset",  &fDriftOffset,  "DriftOffset/F");
  
  fTree->Branch("Plane", &fPlane, "Plane/i");
  fTree->Branch("HitAdcSum", &fHitAdcSum, "HitAdcSum/F");  
  fTree->Branch("HitIntegral", &fHitIntegral, "HitIntegral/F");
  fTree->Branch("HitPeak", &fHitPeak, "HitPeak/F");
  fTree->Branch("HitTime", &fHitTime, "HitTime/F");
  fTree->Branch("Pitch", &fPitch, "Pitch/F");
  fTree->Branch("dQdx", &fdQdx, "dQdx/F");
  fTree->Branch("X", &fX, "X/F");
  fTree->Branch("Y", &fY, "Y/F");
  fTree->Branch("Z", &fZ, "Z/F");
}

//
void pddpana::CosmicsdQdx::endJob()
{}

//
bool pddpana::CosmicsdQdx::checkCutsAndGetT0( const recob::Track& track, float &T0 )
{
  const string myname = "pddpana::CosmicsdQdx::checkCutsAndGetT0: ";
  
  //
  T0 = 0; // track distance in cm from the top

  if( track.Length() < fTrackMinLen ){
    if( fLogLevel>=3 ){
      cout<<myname<<"Failed to pass length cut "<<track.Length()<<endl;
    }
    return false;
  }
  
  auto tstart = track.Vertex();
  auto tend   = track.End();
  
  // cosmics are downward going
  if( tstart.X() > tend.X() ){
    tstart = tend;
  }

  // find TPC volume
  auto const tpc = fGeom->PositionToTPC( tstart );

  // find track start in local TPC coordinates
  auto pos       = tpc.toLocalCoords( tstart );
  
  // check distance from anode 
  // assuming drift along X right now
  T0 = tpc.HalfSizeX() - pos.X();
  if( T0 < fTrackDriftCut ){ 
    if( fLogLevel>=3 ){
      cout<<myname<<"Failed to pass T0 cut "<<T0<<endl;
    }
    return false;
  }

  // check distance from the borders
  if( (tpc.HalfSizeY() - std::abs(pos.Y())) < fTrackWallCut || 
      (tpc.HalfSizeZ() - std::abs(pos.Z())) < fTrackWallCut ){
    
    if( fLogLevel>=3 ){
      cout<<myname<<"Failed to pass wall cut ("
	  <<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<")"<<endl;
    }
    
    return false;
  }

  // float x0 = track.Vertex().X();
  // float y0 = track.Vertex().Y();
  // float z0 = track.Vertex().Z();
  // float x1 = track.End().X();
  // float y1 = track.End().Y();
  // float z1 = track.End().Z();
  return true;
}

DEFINE_ART_MODULE(pddpana::CosmicsdQdx)
