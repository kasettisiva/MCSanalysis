////////////////////////////////////////////////////////////////////////
// Class:       BeamHitFinder
// Plugin Type: producer (art v2_11_02)
// File:        BeamHitFinder_module.cc
//
// Generated at Tue Jun 12 19:56:05 2018 by Tingjun Yang using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ROOT includes
#include "TVector3.h"

#include <memory>

namespace pdune {
  class BeamHitFinder;
}


class pdune::BeamHitFinder : public art::EDProducer {
public:
  explicit BeamHitFinder(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BeamHitFinder(BeamHitFinder const &) = delete;
  BeamHitFinder(BeamHitFinder &&) = delete;
  BeamHitFinder & operator = (BeamHitFinder const &) = delete;
  BeamHitFinder & operator = (BeamHitFinder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  geo::GeometryCore const* fGeom;
  detinfo::DetectorProperties const* fDetProp;

  const art::InputTag fHitModuleLabel;

  TVector3 fBeamPos;
  TVector3 fBeamDir;
  double   fTolerance;
  double   fToleranceWire;

  double tickToDist;

  bool findClosestSP(std::vector< art::Ptr<recob::Hit> > hits, 
                     art::Ptr<recob::Hit> hit, 
                     art::FindManyP< recob::SpacePoint > &spFromHit,
                     TVector3 &point);

};


pdune::BeamHitFinder::BeamHitFinder(fhicl::ParameterSet const & p)
  : EDProducer(p)
  , fHitModuleLabel(p.get< art::InputTag>("HitModuleLabel"))
  , fTolerance(p.get< double >("Tolerance"))
  , fToleranceWire(p.get< double >("ToleranceWire"))
{
  std::vector<double> fpos = p.get< std::vector< double > >("BeamPos");
  if (fpos.size()!=3){
    throw cet::exception("BeamHitFinder")<< "BeamPos needs to be a vector of size 3";
  }
  fBeamPos[0] = fpos[0];
  fBeamPos[1] = fpos[1];
  fBeamPos[2] = fpos[2];

  std::vector<double> fdir = p.get< std::vector< double > >("BeamDir");
  if (fdir.size()!=3){
    throw cet::exception("BeamHitFinder")<< "BeamDir needs to be a vector of size 3";
  }
  fBeamDir[0] = fdir[0];
  fBeamDir[1] = fdir[1];
  fBeamDir[2] = fdir[2];
  fBeamDir = fBeamDir.Unit();

  // let HitCollectionCreator declare that we are going to produce
  // hits and associations with wires and raw digits
  // (with no particular product label)
  recob::HitCollectionCreator::declare_products(*this);
  
  // will also copy associations of SpacePoints to original hits
  produces<art::Assns<recob::Hit, recob::SpacePoint>>();

  // will also copy SpacePoints themselves
  produces< std::vector<recob::SpacePoint>>();

  fGeom = &*(art::ServiceHandle<geo::Geometry>());
  fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  tickToDist = fDetProp->DriftVelocity(fDetProp->Efield(),fDetProp->Temperature());
  tickToDist *= 1.e-3 * fDetProp->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns                                

}

void pdune::BeamHitFinder::produce(art::Event & evt)
{
  auto hitsHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
  
  // also get the associated wires and raw digits;
  // we assume they have been created by the same module as the hits
  art::FindOneP<raw::RawDigit> channelHitRawDigits(hitsHandle, evt, fHitModuleLabel);
  art::FindOneP<recob::Wire>   channelHitWires    (hitsHandle, evt, fHitModuleLabel);

  // this object contains the hit collection
  // and its associations to wires and raw digits:
  recob::HitCollectionCreator hcol(*this, evt,
                                   channelHitWires.isValid(), // doWireAssns
                                   channelHitRawDigits.isValid() // doRawDigitAssns
                                   );

  // here is the copy of associations to hits, based on original hit assns
  auto assns = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>(); 

  // here is the copy of space points
  std::unique_ptr<std::vector<recob::SpacePoint> > spcol(new std::vector<recob::SpacePoint>);

  // all hits in the collection
  std::vector< art::Ptr<recob::Hit> > hits;
  art::fill_ptr_vector(hits, hitsHandle);

  art::FindManyP< recob::SpacePoint > spFromHit(hitsHandle, evt, fHitModuleLabel);

  auto const hitPtrMaker = art::PtrMaker<recob::Hit>(evt);

  std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit>>> hitmap;
  std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit>>> beamhitmap;

  //Check if the associated space point is in the cone
  for (size_t i = 0; i < hits.size(); ++i){

    hitmap[geo::PlaneID(hits[i]->WireID())].push_back(hits[i]);

    TVector3 point; //3D point corresponding to the hit

    auto &sps = spFromHit.at(i);
    if (sps.size()){//associated space point
      point[0] = sps[0]->XYZ()[0];
      point[1] = sps[0]->XYZ()[1];
      point[2] = sps[0]->XYZ()[2];

      double dis = ((point-fBeamPos).Cross(fBeamDir)).Mag();
      if (dis<fTolerance){
        beamhitmap[geo::PlaneID(hits[i]->WireID())].push_back(hits[i]);
        //save hit and associations
        recob::HitCreator new_hit(*(hits[i]));
        hcol.emplace_back(new_hit.move(), channelHitWires.at(i), channelHitRawDigits.at(i));
        if (sps.size()){
          auto hitPtr = hitPtrMaker(hcol.size() - 1);
          for (auto const & spPtr : sps){
            assns->addSingle(hitPtr, spPtr);
            spcol->push_back(*spPtr);
          }
        }
      }//dis<fTolerance
    }//find associated space point
  }//loop over all hits

  //now check reminding hits, see if they are close to beam hits
  for (size_t i = 0; i < hits.size(); ++i){
    auto &sps = spFromHit.at(i);
    if (!sps.size()){//no associated space point
      double wirePitch = fGeom->WirePitch(hits[i]->WireID());
      double UnitsPerTick = tickToDist / wirePitch;
      double x0 = hits[i]->WireID().Wire;
      double y0 = hits[i]->PeakTime() * UnitsPerTick;
      double mindis = DBL_MAX;
      for (auto &hit : beamhitmap[geo::PlaneID(hits[i]->WireID())]){
        double x1 = hit->WireID().Wire;
        double y1 = hit->PeakTime() * UnitsPerTick;
        double dis = sqrt(pow(x1-x0,2)+pow(y1-y0,2));
        if (dis<mindis){
          mindis = dis;
        }
      }
      if (mindis < fToleranceWire){

        //now make sure the hit is not close to a spacepoint out of cone
        TVector3 point; //3D point corresponding to the hit
        if (!findClosestSP(hitmap[geo::PlaneID(hits[i]->WireID())], hits[i], spFromHit, point)) continue;
        double dis = ((point-fBeamPos).Cross(fBeamDir)).Mag();
        if (dis<fTolerance){
          //save hit and associations
          recob::HitCreator new_hit(*(hits[i]));
          hcol.emplace_back(new_hit.move(), channelHitWires.at(i), channelHitRawDigits.at(i));
        }
      }
    }
  }

  // put the hit collection and associations into the event
  hcol.put_into(evt);
  evt.put(std::move(assns));
  evt.put(std::move(spcol));
}

void pdune::BeamHitFinder::beginJob()
{
}

bool pdune::BeamHitFinder::findClosestSP(std::vector< art::Ptr<recob::Hit> > hits, art::Ptr<recob::Hit> hit, art::FindManyP< recob::SpacePoint > &spFromHit, TVector3 &point){
  
  double wirePitch = fGeom->WirePitch(hit->WireID());
  double UnitsPerTick = tickToDist / wirePitch;

  double x0 = hit->WireID().Wire;
  double y0 = hit->PeakTime() * UnitsPerTick;

  double mindis = DBL_MAX;
  for (size_t j = 0; j < hits.size(); ++j){
    if (hits[j].key()==hit.key()) continue;
    if (geo::PlaneID(hit->WireID()) != geo::PlaneID(hits[j]->WireID())) continue;
    //find associated space points
    auto &sps = spFromHit.at(hits[j].key());
    if (!sps.size()) continue;
    double x1 = hits[j]->WireID().Wire;
    double y1 = hits[j]->PeakTime() * UnitsPerTick;
    double dis = sqrt(pow(x1-x0,2)+pow(y1-y0,2));
    if (dis<mindis){
      mindis = dis;
      point[0] = sps[0]->XYZ()[0];
      point[1] = sps[0]->XYZ()[1];
      point[2] = sps[0]->XYZ()[2];
    }
  }//loop over all hits
  
  if (mindis<DBL_MAX)
    return true;
  else
    return false;
}

DEFINE_ART_MODULE(pdune::BeamHitFinder)
