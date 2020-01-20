//////////////////////////////////////////////////////////////////////////
// Class:       RecoHitsEff
// Module Type: analyzer
// File:        RecoHitsEff_module.cc
// Author:      D.Stefan and R.Sulej
//
// Check of hit efficiency reconstruction. 
// It helps to tests various disambiguation algorithms. 
//
//////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "TTree.h"

namespace pdune
{
	class RecoHitsEff;
}

class pdune::RecoHitsEff : public art::EDAnalyzer {
public:

  struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> HitModuleLabel {
          Name("HitModuleLabel"), Comment("tag of hits producer")
      };
      
      fhicl::Atom<art::InputTag> HitDCheatLabel {
      		Name("HitDCheatLabel"), Comment("tag of dcheat hits producer")
      };
      
      fhicl::Atom<int> Plane {
      		Name("Plane"), Comment("check hit disambiguation in selected plane only")
      };
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit RecoHitsEff(Parameters const& config);

  RecoHitsEff(RecoHitsEff const &) = delete;
  RecoHitsEff(RecoHitsEff &&) = delete;
  RecoHitsEff & operator = (RecoHitsEff const &) = delete;
  RecoHitsEff & operator = (RecoHitsEff &&) = delete;

  void analyze(art::Event const & e) override;
	void beginJob() override;
	
private:

  void ResetVars();

  TTree *fTree;
  
  int fRun, fEvent;
	int fNhits; int fNdcheat;
	int fMatching;

	art::InputTag fHitModuleLabel;
	art::InputTag fHitDCheatLabel; 
	int fPlane; 
};

pdune::RecoHitsEff::RecoHitsEff(Parameters const& config) : EDAnalyzer(config),
	fHitModuleLabel(config().HitModuleLabel()),
	fHitDCheatLabel(config().HitDCheatLabel()),
	fPlane(config().Plane())
{
}

void pdune::RecoHitsEff::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;

	fTree = tfs->make<TTree>("events", "summary tree");
	fTree->Branch("fRun", &fRun, "fRun/I");
	fTree->Branch("fEvent", &fEvent, "fEvent/I");
	fTree->Branch("fPlane", &fPlane, "fPlane/I");
	fTree->Branch("fNhits", &fNhits, "fNhits/I");
	fTree->Branch("fNdcheat", &fNdcheat, "fNdcheat/I");
	fTree->Branch("fMatching", &fMatching, "fMatching/I");	
}

void pdune::RecoHitsEff::analyze(art::Event const & evt)
{
  ResetVars();
  
  fRun = evt.run();
  fEvent = evt.id().event();
  
  auto const & hitList = *evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);

  std::vector< size_t > hits;
  // get and sort hits
  int plane = 0; 
  for (size_t k = 0; k < hitList.size(); ++k)
  {
  	plane = hitList[k].WireID().Plane;
  	if (plane != fPlane) continue;
  	fNhits++;
  	
  	hits.push_back(k);
  }
  
  auto const & dcheatList = *evt.getValidHandle< std::vector<recob::Hit> >(fHitDCheatLabel);
  
  std::vector< size_t > hitsdc;
  // get and short hits_dc
  for (size_t k = 0; k < dcheatList.size(); ++k)
  {
  	plane = dcheatList[k].WireID().Plane;
  	if (plane != fPlane) continue;
  	fNdcheat++;
  	
  	hitsdc.push_back(k);
  }
  
	std::cout << " fNhits " << fNhits << std::endl;
	std::cout << " fNdcheat " << fNdcheat << std::endl;

	for (auto const h: hits)
	{
		for (auto const hdc: hitsdc)
		{
			if ((hitList[h].PeakTime() == dcheatList[hdc].PeakTime()) && (hitList[h].WireID() == dcheatList[hdc].WireID()))
			{
				fMatching++;
			}
		}
	}
	
	std::cout << " fMatching " << fMatching << std::endl;

  fTree->Fill();
}

void pdune::RecoHitsEff::ResetVars()
{
	fRun = 0;
	fEvent = 0;
	fNhits = 0;
	fNdcheat = 0;
	fMatching = 0;
}

DEFINE_ART_MODULE(pdune::RecoHitsEff)
