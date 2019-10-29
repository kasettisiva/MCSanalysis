////////////////////////////////////////////////////////////////////////
// Class:       SaveSpacePoints
// Plugin Type: analyzer (art v2_11_02)
// File:        SaveSpacePoints_module.cc
//
// Generated at Wed Jul 11 09:18:07 2018 by Tingjun Yang using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardataobj/RawData/RDTimeStamp.h"

#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"

namespace proto {
  class SaveSpacePoints;
}


class proto::SaveSpacePoints : public art::EDAnalyzer {
public:
  explicit SaveSpacePoints(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SaveSpacePoints(SaveSpacePoints const &) = delete;
  SaveSpacePoints(SaveSpacePoints &&) = delete;
  SaveSpacePoints & operator = (SaveSpacePoints const &) = delete;
  SaveSpacePoints & operator = (SaveSpacePoints &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  const art::InputTag fSpacePointModuleLabel;
  const art::InputTag fBeamModuleLabel;
  const art::InputTag fTrackModuleLabel;
  const art::InputTag fTimeDecoderModuleLabel;

  TTree *fTree;
  // Run information
  int run;
  int subrun;
  int event;
  int trigger;
  double evttime;

  // space point information
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> vz;
  std::vector<double> vcharge;
  std::vector<int> vtrackid;

  // beam information
  std::vector<double> beamPosx;
  std::vector<double> beamPosy;
  std::vector<double> beamPosz;
  
  std::vector<double> beamDirx;
  std::vector<double> beamDiry;
  std::vector<double> beamDirz;

  std::vector<double> beamMomentum;

  double tof;
  short ckov0status;
  short ckov1status;

};


proto::SaveSpacePoints::SaveSpacePoints(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fSpacePointModuleLabel(p.get< art::InputTag >("SpacePointModuleLabel")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),
  fTimeDecoderModuleLabel(p.get< art::InputTag >("TimeDecoderModuleLabel"))
{}

void proto::SaveSpacePoints::analyze(art::Event const & evt)
{
  // Implementation of required member function here.
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  art::Timestamp ts = evt.time();
  //std::cout<<ts.timeHigh()<<" "<<ts.timeLow()<<std::endl;
  if (ts.timeHigh() == 0){
    TTimeStamp tts(ts.timeLow());
    evttime = tts.AsDouble();
  }
  else{
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime = tts.AsDouble();
  }
  vx.clear();
  vy.clear();
  vz.clear();
  vcharge.clear();
  vtrackid.clear();
  beamPosx.clear();
  beamPosy.clear();
  beamPosz.clear();
  beamDirx.clear();
  beamDiry.clear();
  beamDirz.clear();
  beamMomentum.clear();

  // Access the trigger information
  trigger = -1;
  art::ValidHandle<std::vector<raw::RDTimeStamp>> timeStamps = evt.getValidHandle<std::vector<raw::RDTimeStamp>>(fTimeDecoderModuleLabel);

  // Check that we have good information
  if(timeStamps.isValid() && timeStamps->size() == 1){
    // Access the trigger information. Beam trigger flag = 0xc
    const raw::RDTimeStamp& timeStamp = timeStamps->at(0);
    trigger = timeStamp.GetFlags();
  }

  art::Handle< std::vector<recob::SpacePoint> > spsHandle;
  std::vector< art::Ptr<recob::SpacePoint> > sps;
  if (evt.getByLabel(fSpacePointModuleLabel, spsHandle))
    art::fill_ptr_vector(sps, spsHandle);

  art::Handle< std::vector<recob::PointCharge> > pcsHandle;
  std::vector< art::Ptr<recob::PointCharge> > pcs;
  if (evt.getByLabel(fSpacePointModuleLabel, pcsHandle))
    art::fill_ptr_vector(pcs, pcsHandle);

  for (size_t i = 0; i<sps.size(); ++i){
    vx.push_back(sps[i]->XYZ()[0]);
    vy.push_back(sps[i]->XYZ()[1]);
    vz.push_back(sps[i]->XYZ()[2]);
    vcharge.push_back(pcs[i]->charge());
    vtrackid.push_back(-1);
  }

  art::Handle< std::vector<recob::Track> > trkHandle;
  std::vector< art::Ptr<recob::Track> > trks;
  if (evt.getByLabel(fTrackModuleLabel, trkHandle))
    art::fill_ptr_vector(trks, trkHandle);

  for (size_t i = 0; i<trks.size(); ++i){
    auto & trk = trks[i];
    for (size_t j = 0; j<trk->NPoints(); ++j){
      if (trk->HasValidPoint(j)){
        vx.push_back(trk->TrajectoryPoint(j).position.X());
        vy.push_back(trk->TrajectoryPoint(j).position.Y());
        vz.push_back(trk->TrajectoryPoint(j).position.Z());
        vcharge.push_back(0);
        vtrackid.push_back(trk->ID());
      }
    }
  }

  art::Handle< std::vector<beam::ProtoDUNEBeamEvent> > pdbeamHandle;
  std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;
  if (evt.getByLabel(fBeamModuleLabel, pdbeamHandle))
    art::fill_ptr_vector(beaminfo, pdbeamHandle);
  else{
    std::cout<<"No beam information from "<<fBeamModuleLabel<<std::endl;
  }

  tof = -1;
  ckov0status = -1;
  ckov1status = -1;
  if (beaminfo.size()){
    if (beaminfo[0]->GetTimingTrigger() == 12){
      if (beaminfo[0]->CheckIsMatched()){
        //Get TOF info
        if (beaminfo[0]->GetTOFChan() != -1){//if TOFChan == -1, then there was not a successful match, if it's 0, 1, 2, or 3, then there was a good match corresponding to the different pair-wise combinations of the upstream and downstream channels
          tof = beaminfo[0]->GetTOF();
        }
        //Get beam particle trajectory info
        auto & tracks = beaminfo[0]->GetBeamTracks();
        for (size_t i = 0; i<tracks.size(); ++i){
          beamPosx.push_back(tracks[i].End().X());
          beamPosy.push_back(tracks[i].End().Y());
          beamPosz.push_back(tracks[i].End().Z());
          beamDirx.push_back(tracks[i].StartDirection().X());
          beamDiry.push_back(tracks[i].StartDirection().Y());
          beamDirz.push_back(tracks[i].StartDirection().Z());
        }
        //Get reconstructed beam momentum info
        auto & beammom = beaminfo[0]->GetRecoBeamMomenta();
        for (size_t i = 0; i<beammom.size(); ++i){
          beamMomentum.push_back(beammom[i]);
        }
      }
    }
    if (beaminfo[0]->GetBITrigger() == 1){
      //Get CKov status
      ckov0status = beaminfo[0]->GetCKov0Status();
      ckov1status = beaminfo[0]->GetCKov1Status();
    }
  }

  fTree->Fill();
}

void proto::SaveSpacePoints::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("spt","space point tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("trigger",&trigger,"trigger/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("vx",&vx);
  fTree->Branch("vy",&vy);
  fTree->Branch("vz",&vz);
  fTree->Branch("vcharge",&vcharge);
  fTree->Branch("vtrackid",&vtrackid);
  fTree->Branch("beamPosx",&beamPosx);
  fTree->Branch("beamPosy",&beamPosy);
  fTree->Branch("beamPosz",&beamPosz);
  fTree->Branch("beamDirx",&beamDirx);
  fTree->Branch("beamDiry",&beamDiry);
  fTree->Branch("beamDirz",&beamDirz);
  fTree->Branch("beamMomentum",&beamMomentum);
  fTree->Branch("tof", &tof, "tof/D");
  fTree->Branch("ckov0status", &ckov0status, "ckov0status/S");
  fTree->Branch("ckov1status", &ckov1status, "ckov1status/S");
}

DEFINE_ART_MODULE(proto::SaveSpacePoints)
