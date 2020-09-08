//////////////////////////////////////////////////////////////////////////
// Class:       RecoStats
// Module Type: analyzer
// File:        RecoStatistics_module.cc
// Author:      D.Stefan and R.Sulej
//
// Various tests of the reco properties specific to ProtoDUNE SP. Started
// with reco charge per TPC, for the purpose of optimizing the decision
// on the locations of the instrumented TPC's. 
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

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "TH2.h"
#include "TTree.h"

namespace pdune
{
	class RecoStats;
}

class pdune::RecoStats : public art::EDAnalyzer {
public:

  struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> EmTrkModuleLabel {
          Name("EmTrkModuleLabel"), Comment("tag of hits classifier")
      };

      fhicl::Table<calo::CalorimetryAlg::Config> CalorimetryAlg {
          Name("CalorimetryAlg"), Comment("calo alg to calibrate hit adc to energy")
      };
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit RecoStats(Parameters const& config);

  RecoStats(RecoStats const &) = delete;
  RecoStats(RecoStats &&) = delete;
  RecoStats & operator = (RecoStats const &) = delete;
  RecoStats & operator = (RecoStats &&) = delete;

  void analyze(art::Event const & e) override;
  void beginJob() override;
  void endJob() override;

private:
  template <size_t N> bool fillChargeHistos(const art::Event & evt,
                                        detinfo::DetectorClocksData const& clockData,
                                        detinfo::DetectorPropertiesData const& detProp);
  void fillChargeHisto(detinfo::DetectorClocksData const& clockData,
                       detinfo::DetectorPropertiesData const& detProp,
                       TH2F* histo, const std::vector<art::Ptr<recob::Hit>> & hits);
  double getHitMeV(detinfo::DetectorClocksData const& clockData,
                   detinfo::DetectorPropertiesData const& detProp,
                   const recob::Hit & hit) const;
  void ResetVars();

  double fElectronsToGeV;
  std::unordered_map<int, std::pair<int, int>> mapTPCtoXY;

  TH2F* fEmChargePerTpcHist;
  TH2F* fHadChargePerTpcHist;
  TH2F* fChargePerTpcHist;

  TTree *fTree;
  
  int fRun, fEvent;

  float fEkGen, fEkDep; // <- to be filled

  // ******** FHiCL parameters: ********
  art::InputTag fEmTrkModuleLabel;

  calo::CalorimetryAlg fCalorimetryAlg;
};

pdune::RecoStats::RecoStats(Parameters const& config) : EDAnalyzer(config),
    fEmTrkModuleLabel(config().EmTrkModuleLabel()),
    fCalorimetryAlg(config().CalorimetryAlg())
{
	art::ServiceHandle<sim::LArG4Parameters> larParameters;
    fElectronsToGeV = 1./larParameters->GeVToElectrons();

    mapTPCtoXY[0]  = std::make_pair(-1, -1); // <- in case something falls into not active tpc
    mapTPCtoXY[1]  = std::make_pair(0, 1);   // <- make 2D histo bins looking like in evd ortho
    mapTPCtoXY[2]  = std::make_pair(0, 0);
    mapTPCtoXY[3]  = std::make_pair(-1, -1);
    mapTPCtoXY[4]  = std::make_pair(-1, -1);
    mapTPCtoXY[5]  = std::make_pair(1, 1);
    mapTPCtoXY[6]  = std::make_pair(1, 0);
    mapTPCtoXY[7]  = std::make_pair(-1, -1);
    mapTPCtoXY[8]  = std::make_pair(-1, -1);
    mapTPCtoXY[9]  = std::make_pair(2, 1);
    mapTPCtoXY[10] = std::make_pair(2, 0);
    mapTPCtoXY[11] = std::make_pair(-1, -1);
}

void pdune::RecoStats::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;
	
	fEmChargePerTpcHist = tfs->make<TH2F>("EmChargePerTPC", "EM charge reconstucted", 3, 0., 3., 2, 0., 2.);
	fHadChargePerTpcHist = tfs->make<TH2F>("HadChargePerTPC", "hadronic charge reconstucted", 3, 0., 3., 2, 0., 2.);
	fChargePerTpcHist = tfs->make<TH2F>("ChargePerTPC", "total charge reconstucted", 3, 0., 3., 2, 0., 2.);

	fTree = tfs->make<TTree>("events", "summary tree");
	fTree->Branch("fRun", &fRun, "fRun/I");
	fTree->Branch("fEvent", &fEvent, "fEvent/I");
	fTree->Branch("fEkGen", &fEkGen, "fEkGen/F");
	fTree->Branch("fEkDep", &fEkDep, "fEkDep/F");
}

void pdune::RecoStats::endJob()
{
}

template <size_t N>
bool pdune::RecoStats::fillChargeHistos(const art::Event & evt,
                                        detinfo::DetectorClocksData const& clockData,
                                        detinfo::DetectorPropertiesData const& detProp)
{
    float trackLikeThr = 0.63;  // should move to fcl params
    geo::View_t view = geo::kZ; // can move to fcl params

    auto cluResults = anab::MVAReader< recob::Cluster, N >::create(evt, fEmTrkModuleLabel);
    if (!cluResults) { return false; }

    int trkLikeIdx = cluResults->getIndex("track");
    int emLikeIdx = cluResults->getIndex("em");
    if ((trkLikeIdx < 0) || (emLikeIdx < 0)) { return false; }

    const art::FindManyP< recob::Hit > hitsFromClusters(cluResults->dataHandle(), evt, cluResults->dataTag());
    const auto & cnnOuts = cluResults->outputs();
    const auto & clusters = cluResults->items();

    for (size_t i = 0; i < clusters.size(); ++i)
    {
        if (clusters[i].View() == view) { continue; } // take just one view

        double trackLike, trkOrEm = cnnOuts[i][trkLikeIdx] + cnnOuts[i][emLikeIdx];
        if (trkOrEm > 0) { trackLike = cnnOuts[i][trkLikeIdx] / trkOrEm; }
        else { trackLike = 0; }

        if (trackLike > trackLikeThr)
        {
          fillChargeHisto(clockData, detProp, fHadChargePerTpcHist, hitsFromClusters.at(i));
        }
        else
        {
            fillChargeHisto(clockData, detProp, fEmChargePerTpcHist, hitsFromClusters.at(i));
        }
        fillChargeHisto(clockData, detProp, fChargePerTpcHist, hitsFromClusters.at(i));
    }
    return true;
}

void pdune::RecoStats::fillChargeHisto(detinfo::DetectorClocksData const& clockData,
                                       detinfo::DetectorPropertiesData const& detProp,
                                       TH2F* histo, const std::vector<art::Ptr<recob::Hit>> & hits)
{
    for (const auto & h : hits)
    {
        auto xy = mapTPCtoXY[ h->WireID().TPC ];
        histo->Fill(xy.first, xy.second, getHitMeV(clockData, detProp, *h));
    }
}

double pdune::RecoStats::getHitMeV(detinfo::DetectorClocksData const& clockData,
                                   detinfo::DetectorPropertiesData const& detProp,
                                   const recob::Hit & hit) const
{
	double adc = hit.Integral();
	if (!std::isnormal(adc) || (adc < 0)) adc = 0.0;

	unsigned short plane = hit.WireID().Plane;
	double tdrift = hit.PeakTime();
	double dqel = fCalorimetryAlg.ElectronsFromADCArea(adc, plane);

        double correllifetime = fCalorimetryAlg.LifetimeCorrection(clockData, detProp, tdrift, 0); // for the moment T0 = 0 (beam particles)
	double dq = dqel * correllifetime * fElectronsToGeV * 1000;
	if (!std::isnormal(dq) || (dq < 0)) dq = 0.0;

	return dq; 
}

void pdune::RecoStats::analyze(art::Event const & evt)
{
  ResetVars();
  
  fRun = evt.run();
  fEvent = evt.id().event();

  // try to dig out 4- or 3-output MVA data product
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
  if (!fillChargeHistos<4>(evt, clockData, detProp) &&
      !fillChargeHistos<3>(evt, clockData, detProp))
  {
    throw cet::exception("PMAlgTrackMaker") << "No EM/track MVA data products." << std::endl;
  }

  fTree->Fill();
}

void pdune::RecoStats::ResetVars()
{
	fRun = 0;
	fEvent = 0;
	fEkGen = 0;
	fEkDep = 0;
}

DEFINE_ART_MODULE(pdune::RecoStats)
