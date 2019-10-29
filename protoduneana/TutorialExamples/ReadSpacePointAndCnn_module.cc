////////////////////////////////////////////////////////////////////////
// Class:       ReadSpacePointAndCnn
// Module Type: analyzer
// File:        ReadSpacePointAndCnn_module.cc
//
// Example prepared for the computing tutorial, saves SpacePoint
// coordinates in the ROOT tree, together with EM/track classification
// obtained from CNN. Output file to be displayed with SWAN example.
// Robert Sulej
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/ArtDataHelper/MVAReader.h"
#define MVA_LENGTH 4

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include "TTree.h"

namespace tutorial {

class ReadSpacePointAndCnn : public art::EDAnalyzer {
public:

  struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> HitsModuleLabel {
          Name("HitsModuleLabel"), Comment("tag of hits producer")
      };

      fhicl::Atom<art::InputTag> SpacePointModuleLabel {
          Name("SpacePointModuleLabel"), Comment("tag of spacepoint producer")
      };

      fhicl::Atom<art::InputTag> CnnModuleLabel {
          Name("CnnModuleLabel"), Comment("tag of CNN module for EM/track id")
      };
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit ReadSpacePointAndCnn(Parameters const& config);

  ReadSpacePointAndCnn(ReadSpacePointAndCnn const &) = delete;
  ReadSpacePointAndCnn(ReadSpacePointAndCnn &&) = delete;
  ReadSpacePointAndCnn & operator = (ReadSpacePointAndCnn const &) = delete;
  ReadSpacePointAndCnn & operator = (ReadSpacePointAndCnn &&) = delete;

  void analyze(art::Event const & e) override;

  void beginJob() override;
  void endJob() override;

private:
  void clear();

  size_t fEvNumber;

  TTree *fEventTree;

  // just some summary info
  size_t fNHits[3];
  size_t fNPoints;

  // let's sort hits by planes already now, it really makes
  // any downstream work easier
  std::vector<unsigned char> fHitTpc[3];
  std::vector<float> fHitWire[3], fHitTime[3], fHitCharge[3];

  // SpacePoint 3D positions, charge (uncalibrated!) and EM-like score
  std::vector<float> fPointX, fPointY, fPointZ;
  std::vector<float> fPointCharge;
  std::vector<float> fPointEmScore;

  // ******* fcl parameters *******
  art::InputTag fHitsModuleLabel;
  art::InputTag fSpacePointModuleLabel;
  art::InputTag fCnnModuleLabel;
};

ReadSpacePointAndCnn::ReadSpacePointAndCnn(Parameters const& config) : art::EDAnalyzer(config),
    fHitsModuleLabel(config().HitsModuleLabel()),
    fSpacePointModuleLabel(config().SpacePointModuleLabel()),
    fCnnModuleLabel(config().CnnModuleLabel())
{
}

void ReadSpacePointAndCnn::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs; // TTree's are created in the memory managed by ROOT (you don't delete them)

    fEventTree = tfs->make<TTree>("EventTree", "event by event info");
    fEventTree->Branch("event", &fEvNumber, "fEvNumber/I");
    fEventTree->Branch("nhits", &fNHits, "fNHits[3]/I");
    fEventTree->Branch("npoints", &fNPoints, "fNPoints/I");

    fEventTree->Branch("hit0tpc", &fHitTpc[0]);
    fEventTree->Branch("hit0w", &fHitWire[0]);
    fEventTree->Branch("hit0t", &fHitTime[0]);
    fEventTree->Branch("hit0q", &fHitCharge[0]);

    fEventTree->Branch("hit1tpc", &fHitTpc[1]);
    fEventTree->Branch("hit1w", &fHitWire[1]);
    fEventTree->Branch("hit1t", &fHitTime[1]);
    fEventTree->Branch("hit1q", &fHitCharge[1]);

    fEventTree->Branch("hit2tpc", &fHitTpc[2]);
    fEventTree->Branch("hit2w", &fHitWire[2]);
    fEventTree->Branch("hit2t", &fHitTime[2]);
    fEventTree->Branch("hit2q", &fHitCharge[2]);

    fEventTree->Branch("pointx", &fPointX);
    fEventTree->Branch("pointy", &fPointY);
    fEventTree->Branch("pointz", &fPointZ);
    fEventTree->Branch("pointq", &fPointCharge);
    fEventTree->Branch("emscore", &fPointEmScore);
}

void ReadSpacePointAndCnn::clear()
{
    for (size_t p = 0; p < 3; ++p)
    {
        fHitTpc[p].clear();
        fHitWire[p].clear();
        fHitTime[p].clear();
        fHitCharge[p].clear();
    }
}

void ReadSpacePointAndCnn::analyze(art::Event const & evt)
{
    clear();

    fEvNumber = evt.id().event();
    mf::LogVerbatim("ReadSpacePointAndCnn") << "ReadSpacePointAndCnn module on event #" << fEvNumber;

    // store 2D hits info, sorted by plane (but not by TPC, you need to select
    // on that later in your analysis)
    auto hitsHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitsModuleLabel);
    fNHits[0] = 0; fNHits[1] = 0; fNHits[2] = 0;

    for (const auto & h : *hitsHandle)
    {
        size_t p = h.WireID().Plane;
        fHitTpc[p].push_back(h.WireID().TPC);
        fHitWire[p].push_back(h.WireID().Wire);
        fHitTime[p].push_back(h.PeakTime());
        fHitCharge[p].push_back(h.Integral());
        fNHits[p]++;
    }

    // store SpacePoints basic info
    auto spHandle = evt.getValidHandle< std::vector<recob::SpacePoint> >(fSpacePointModuleLabel);
    auto qHandle = evt.getValidHandle< std::vector<recob::PointCharge> >(fSpacePointModuleLabel);
    if (spHandle->size() != qHandle->size())
    {
        throw cet::exception("tutorial::ReadSpacePointAndCnn")
            << "size of point and charge containers must be equal" << std::endl;
    }

    fNPoints = spHandle->size();
    fPointX.resize(fNPoints); fPointY.resize(fNPoints); fPointZ.resize(fNPoints);
    fPointCharge.resize(fNPoints);
    fPointEmScore.resize(fNPoints);

    for (size_t i = 0; i < spHandle->size(); ++i)
    {
        fPointX[i] = (*spHandle)[i].XYZ()[0];
        fPointY[i] = (*spHandle)[i].XYZ()[1];
        fPointZ[i] = (*spHandle)[i].XYZ()[2];
        fPointCharge[i] = (*qHandle)[i].charge();
        fPointEmScore[i] = 0.5; // neutral P(EM-like) value
    }

    // and tag SpacePoints with EM-like scores calculated by CNN on the cluster level
    // connect these scores to SpacePoints throuh hits, and choose the score of the
    // largest cluster containing the hit associated to SpacePoint;
    // other option could be e.g.: use CNN scores on the hit level, multiplying scores
    // from hits in all planes
    auto cluResults = anab::MVAReader<recob::Cluster, MVA_LENGTH>::create(evt, fCnnModuleLabel);
    if (cluResults)
    {
        size_t emLikeIdx = cluResults->getIndex("em"); // at which index EM-like is stored in CNN output vector

        const art::FindManyP<recob::Hit> hitsFromClusters(cluResults->dataHandle(), evt, cluResults->dataTag());
        const art::FindManyP<recob::SpacePoint> spFromHits(hitsHandle, evt, fSpacePointModuleLabel);

        std::vector<size_t> sizeScore(fNPoints, 0); // keep track of the max size of a cluster containing hit associated to spacepoint

        for (size_t c = 0; c < cluResults->size(); ++c)
        {
            //const recob::Cluster & clu = cluResults->item(c); // one can get the cluster object in this way, here not really used

     	    const std::vector< art::Ptr<recob::Hit> > & hits = hitsFromClusters.at(c);
    	    std::array<float, MVA_LENGTH> cnn_out = cluResults->getOutput(c);

            for (const auto & hptr : hits)
            {
                const std::vector< art::Ptr<recob::SpacePoint> > & sp = spFromHits.at(hptr.key());
                for (const auto & spptr : sp) // with SpacePointSolver should be just 1 hit, but be prepared for any other algorithm
                {
                    if (hits.size() > sizeScore[spptr.key()])
                    {
                        sizeScore[spptr.key()] = hits.size();
                        fPointEmScore[spptr.key()] = cnn_out[emLikeIdx];
                    }
                }
             }
    	}
    }

    fEventTree->Fill();

    clear();
}

void ReadSpacePointAndCnn::endJob()
{
    mf::LogVerbatim("ReadSpacePointAndCnn") << "ReadSpacePointAndCnn finished job";
}

} // tutorial namespace

DEFINE_ART_MODULE(tutorial::ReadSpacePointAndCnn)
