///////////////////////////////////////////////////////////////////////                    
// Class:       DEdx                                                                        
// Module Type: analyzer                                                                    
// File:        DEdx_module.cc                  
//                                            
// This is a simple module that makes plot of dE/dx using 3D reconstructed tracks.
// It cab be run before each MC production on the sample of generated muon tracks.
// It can be also used on the data with reconstructed 3D tracks tagged as cosmic muons. 
//
// Authors: Casandra Hazel Morris (University of Houston) and Dorota Stefan (CERN/NCBJ)
//                                                                                                                                                     
////////////////////////////////////////////////////////////////////////                    

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardata/Utilities/DatabaseUtil.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TH1.h"
#include "TTree.h"

namespace proto
{
  class DEdx;
}

class proto::DEdx : public art::EDAnalyzer {
public:
  explicit DEdx(fhicl::ParameterSet const & p);     

  //plugins should not be copied or assigned.         
  DEdx(DEdx const &) = delete;
  DEdx(DEdx &&) = delete;
  DEdx & operator = (DEdx const &) = delete;
  DEdx & operator = (DEdx &&) = delete;

  // Required functions.                                          
  void analyze(art::Event const & e) override;
  void beginJob() override;
  void endJob() override;

  void reconfigure(fhicl::ParameterSet const& p) ;

private:

  // Declare member data here.                      

  // void CountdEdx(const std::vector < art::Ptr< recob::Hit > > & hits, const std::vector < recob::TrackHitMeta const* > & data); // MeV/cm
  void CountdEdx(detinfo::DetectorClocksData const& clockData,
                 detinfo::DetectorPropertiesData const& detProp,
                 const std::vector < art::Ptr< recob::Hit > > & hits, const std::vector< recob::TrackHitMeta const* > & data);
  void ResetVars();

  bool fCosmics;

  TTree *fTree;
  TTree *fTreere;
  int fRun;
  int fEvent;
  int fChosenView;
  size_t fNumberOfTracks;
  size_t fNumberOfTaggedTracks;

  double fdQdx;
  double fdEdx;
  double fdQ;
  double fdx;

  // Module labels to get data products  
  std::string fHitModuleLabel;
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;
  calo::CalorimetryAlg fCalorimetryAlg;

};

proto::DEdx::DEdx(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fCalorimetryAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg"))
  // More initializers here.                                
{
  reconfigure(p);
}

void proto::DEdx::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTreere = tfs->make<TTree>("entries","entries tree");
  fTree = tfs->make<TTree>("dedx","dedx tree");

  fTreere->Branch("fRun", &fRun, "fRun/I");
  fTreere->Branch("fEvent", &fEvent, "fEvent/I");
  fTreere->Branch("fChosenView", &fChosenView, "fChosenView/I");
  fTreere->Branch("fNumberOfTracks", &fNumberOfTracks, "fNumberOfTracks/I");
  fTreere->Branch("fNumberOfTaggedTracks", &fNumberOfTaggedTracks, "fNumberOfTaggedTracks/I");

  fTree->Branch("fdQdx", &fdQdx, "fdQdx/D");
  fTree->Branch("fdEdx", &fdEdx, "fdEdx/D");
  fTree->Branch("fdQ", &fdQ, "fdQ/D");
  fTree->Branch("fdx", &fdx, "fdx/D");
}

void proto::DEdx::endJob()
{
}

void proto::DEdx::reconfigure(fhicl::ParameterSet const & p)
{
  fChosenView = p.get<int>("ChosenView");
  fCosmics = p.get<bool>("Cosmics");
  fHitModuleLabel = p.get< std::string >("HitModuleLabel");
  fTrackModuleLabel = p.get< std::string >("TrackModuleLabel");
}

void proto::DEdx::analyze(art::Event const & e)
{
  // Implementation of required member function here.     
  ResetVars();
  fRun = e.run();
  fEvent = e.id().event();

  // tracks                                                                                             
  auto trkHandle = e.getValidHandle< std::vector<recob::Track> >(fTrackModuleLabel);
  fNumberOfTracks = trkHandle->size();

  const art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trkHandle, e, fTrackModuleLabel);

  // Find the tagged tracks as cosmic muons
  const art::FindManyP<anab::CosmicTag> ct(trkHandle, e, fTrackModuleLabel);

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e, clockData);
  if (fCosmics)
  {
    if (ct.isValid())
    {	
     	// loop over tracks
     	 for (size_t t = 0; t < trkHandle->size(); ++t)
     	 { 
		if (ct.at(t).size())
		{ 
			fNumberOfTaggedTracks++;
					                        
			auto vhit = fmthm.at(t);
			auto vmeta = fmthm.data(t);
	
                        CountdEdx(clockData, detProp, vhit, vmeta);
		}
      	 }
    }
  }	
  else if (fmthm.isValid())
  {
	
      // loop over tracks
      for (size_t t = 0; t < trkHandle->size(); ++t)
      {                           
	auto vhit = fmthm.at(t);
	auto vmeta = fmthm.data(t);
	
        CountdEdx(clockData, detProp, vhit, vmeta);
      }
  }

  fTreere->Fill();
}

void proto::DEdx::CountdEdx(detinfo::DetectorClocksData const& clockData,
                            detinfo::DetectorPropertiesData const& detProp,
                            const std::vector < art::Ptr< recob::Hit > > & hits,
                            const std::vector< recob::TrackHitMeta const* > & data) // MeV/cm
{
  for (size_t h = 0; h < hits.size(); ++h)
    {
	unsigned short plane = hits[h]->WireID().Plane;
      	unsigned short time = hits[h]->PeakTime();

	if (plane == fChosenView)
	{
      		double dqadc = hits[h]->Integral();
      		if (!std::isnormal(dqadc) || (dqadc < 0)) dqadc = 0.0;

      		double t0 = 0;

     		fdQdx = 0.0;
      		fdQ = dqadc;
      		fdx = data[h]->Dx();
      		if ((fdx > 0) && (fdQ > 0))
        	{
        	 fdQdx = fdQ/fdx;
                 fdEdx = fCalorimetryAlg.dEdx_AREA(clockData, detProp, fdQdx, time, plane, t0);
         	 if (fdEdx > 35) fdEdx = 35;
        	}
      		else if ((fdx == 0) && (fdQ > 0))
       	 	{
//	  		std::cout << " the charge is positive but dx is equal zero " << std::endl;
        	}

      		fTree->Fill();
	}
    }
}

void proto::DEdx::ResetVars()
{
  fdQdx = 0.0;
  fdEdx = 0.0;
  fdQ = 0.0;
  fdx = 0.0;
  fNumberOfTracks = 0;
  fNumberOfTaggedTracks = 0;
}


DEFINE_ART_MODULE(proto::DEdx)
