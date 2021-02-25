////////////////////////////////////////////////////////////////////////
// Class:       EMCNNCheck
// Plugin Type: analyzer (art v3_03_01)
// File:        EMCNNCheck_module.cc
//
// Generated at Wed Jan  8 21:50:23 2020 by Tingjun Yang using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamCuts.h"

#include "TTree.h"

#include <iostream>

namespace pdsp {
  class EMCNNCheck;
}


class pdsp::EMCNNCheck : public art::EDAnalyzer {
public:
  explicit EMCNNCheck(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EMCNNCheck(EMCNNCheck const&) = delete;
  EMCNNCheck(EMCNNCheck&&) = delete;
  EMCNNCheck& operator=(EMCNNCheck const&) = delete;
  EMCNNCheck& operator=(EMCNNCheck&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginJob() override;

private:

  //int fselectpdg;
  std::string fGeneratorTag;
  std::string fCNNTag;
  //fhicl::ParameterSet BeamCuts;
  protoana::ProtoDUNEBeamCuts beam_cuts;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;

  TTree *ftree;
  int run;
  int subrun;
  int event;
  int beampdg;
  double average_score_em;
  double average_score_trk;
  double average_score_mic;
  std::vector<short> channel;
  std::vector<short> tpc;
  std::vector<short> plane;
  std::vector<short> wire;
  std::vector<double> charge;
  std::vector<double> peakt;
  std::vector<double> score_em;
  std::vector<double> score_trk;
  std::vector<double> score_mic;
  std::vector<int> pdg;
  std::vector<int> origin;

};


pdsp::EMCNNCheck::EMCNNCheck(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
//fselectpdg(p.get<int>("selectpdg")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fCNNTag(p.get<std::string>("CNNTag")),
  beam_cuts(p.get<fhicl::ParameterSet>("BeamCuts")),
  fBeamlineUtils(p.get<fhicl::ParameterSet>("BeamlineUtils"))
  //BeamCuts(p.get<fhicl::ParameterSet>("BeamCuts"))
  {
    //beam_cuts = protoana::ProtoDUNEBeamCuts(BeamCuts);
  }

void pdsp::EMCNNCheck::analyze(art::Event const& e)
{

  run = e.run();
  subrun = e.subRun();
  event = e.id().event();
  beampdg = 0;
  average_score_em  = 0.;
  average_score_trk = 0.;
  average_score_mic = 0.;
  channel.clear();
  tpc.clear();
  plane.clear();
  wire.clear();
  charge.clear();
  peakt.clear();
  score_em.clear();
  score_trk.clear();
  score_mic.clear();
  pdg.clear();
  origin.clear();

  //Services
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
 
  art::Handle < std::vector < recob::Slice > > sliceListHandle;
  std::vector < art::Ptr < recob::Slice > > sliceList;
  if (e.getByLabel("pandora", sliceListHandle)) {
    art::fill_ptr_vector(sliceList, sliceListHandle);
  }
  else return;

  // Get all pfparticles
  art::Handle < std::vector < recob::PFParticle > > pfpListHandle;
  std::vector < art::Ptr < recob::PFParticle > > pfpList;
  if (e.getByLabel("pandora", pfpListHandle)) {
    art::fill_ptr_vector(pfpList, pfpListHandle);
  }

  // Get all clusters
  art::Handle < std::vector < recob::Cluster > > cluListHandle;
  std::vector < art::Ptr < recob::Cluster > > cluList;
  if (e.getByLabel("pandora", cluListHandle)) {
    art::fill_ptr_vector(cluList, cluListHandle);
  }

  // Get cluster-PFParticle association
  art::FindManyP<recob::Cluster> fmcpfp(pfpListHandle, e, "pandora");

  // Get hit-cluster association
  art::FindManyP<recob::Hit> fmhc(cluListHandle, e, "pandora");

  art::FindManyP <recob::Hit> hitsFromSlice(sliceListHandle, e, "pandora");

  anab::MVAReader<recob::Hit,4> hitResults(e, fCNNTag);

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  if (!e.isRealData()){
    // Get the truth utility to help us out
    protoana::ProtoDUNETruthUtils truthUtil;
    // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
    // simulation has fGeneratorTag = "generator"
    auto mcTruths = e.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
    // of these, so we pass the first element into the function to get the good particle
    const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],e);
    if(geantGoodParticle != 0x0){
      //std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() << std::endl;
//      if (fselectpdg==211){
//        if (geantGoodParticle->PdgCode()!=211 &&
//            geantGoodParticle->PdgCode()!=13){
//          return;
//        }
//      }
//      else{
//        if (geantGoodParticle->PdgCode()!=fselectpdg) return;
//      }
      beampdg = geantGoodParticle->PdgCode();
    }
  }
  else{
    //Access the Beam Event
    auto beamHandle = e.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("beamevent");
  
    std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
    if( beamHandle.isValid()){
      art::fill_ptr_vector(beamVec, beamHandle);
    }
    
    if (beamVec.empty()) return;

    const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); //Should just have one
    /////////////////////////////////////////////////////////////
  
  
    //Check the quality of the event
    std::cout << "Timing Trigger: " << beamEvent.GetTimingTrigger() << std::endl; 
    std::cout << "Is Matched: "     << beamEvent.CheckIsMatched() << std::endl << std::endl;
    
    if( !fBeamlineUtils.IsGoodBeamlineTrigger( e ) ){
      std::cout << "Failed quality check" << std::endl;
      return;
    }
    //Access PID
    std::vector< int > pids = fBeamlineUtils.GetPID( beamEvent, 1. );
//    bool foundparticle = false;
//    for( size_t i = 0; i < pids.size(); ++i ){
//      //std::cout << pids[i] << std::endl;
//      if (pids[i] == fselectpdg) foundparticle = true;
//    }
//    if (!foundparticle) return;
    if (pids.empty()) return;
    beampdg = pids[0];
  }

  if (!beampdg) return;

  //std::cout<<"Found pion"<<std::endl;

  // Get the PFParticle utility
  protoana::ProtoDUNEPFParticleUtils pfpUtil;

  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(e,"pandora");

  if(beamParticles.size() == 0){
    std::cerr << "We found no beam particles for this event... moving on" << std::endl;
    return;
  }

  // We can now look at these particles
  for(const recob::PFParticle* particle : beamParticles){

    const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,e,"pandora","pandoraTrack");
    const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,e,"pandora","pandoraShower");
    if (thisTrack){
      if (!beam_cuts.IsBeamlike(*thisTrack, e, "1")) return;
    }
    if (thisShower){
      if (!beam_cuts.IsBeamlike(*thisShower, e, "1")) return;
    }
  }

  //int sliceid = pfpUtil.GetBeamSlice(e, "pandora");
  
  //if (sliceid!=9999){
  if (fmcpfp.isValid()){
    // Get clusters associated with pfparticle
    auto const& clusters = fmcpfp.at(beamParticles[0]->Self());
    for (auto const & cluster : clusters){
      if (fmhc.isValid()){
        // Get hits associated with cluster
        auto const& hits = fmhc.at(cluster.key());
        //auto const& hits = hitsFromSlice.at(sliceid);
        //std::cout<<hits.size()<<std::endl;
        for (auto & hit : hits){
          std::array<float,4> cnn_out = hitResults.getOutput(hit);
          if (hit->WireID().Plane == 2){
            channel.push_back(hit->Channel());
            tpc.push_back(hit->WireID().TPC);
            plane.push_back(hit->WireID().Plane);
            wire.push_back(hit->WireID().Wire);
            charge.push_back(hit->Integral());
            peakt.push_back(hit->PeakTime());     
            score_em.push_back(cnn_out[hitResults.getIndex("em")]);
            score_trk.push_back(cnn_out[hitResults.getIndex("track")]);
            score_mic.push_back(cnn_out[hitResults.getIndex("michel")]);
            int this_pdg = 0;
            int this_origin = -1;
            if (!e.isRealData()){
              int TrackID = 0;
              std::map<int,double> trkide;
              std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToEveTrackIDEs(clockData, hit);
              for(size_t e = 0; e < TrackIDs.size(); ++e){
                trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
              }	    
              
              // Work out which IDE despoited the most charge in the hit if there was more than one.
              double maxe = -1;
              double tote = 0;
              for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
                tote += ii->second;
                if ((ii->second)>maxe){
                  maxe = ii->second;
                  TrackID = ii->first;
                }
              }
              // Now have trackID, so get PdG code and T0 etc.
              const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackID);
              if (particle){
                this_pdg = particle->PdgCode();
                this_origin = pi_serv->ParticleToMCTruth_P(particle)->Origin();
              }
            }
            pdg.push_back(this_pdg);
            origin.push_back(this_origin);
          }
        }
      }
    }
  }

  // Get the average of the collection plane scores
  unsigned int nCollectionHits = 0;
  for(unsigned int h = 0; h < plane.size(); ++h){
    if(plane.at(h) == 2){
      ++nCollectionHits;
      average_score_em += score_em.at(h);
      average_score_trk += score_trk.at(h);
      average_score_mic += score_mic.at(h);
    }
  }
  if(nCollectionHits > 0){
    average_score_em /= static_cast<double>(nCollectionHits);
    average_score_trk /= static_cast<double>(nCollectionHits);
    average_score_mic /= static_cast<double>(nCollectionHits);
  }

  if (!channel.empty()) ftree->Fill();
}

void pdsp::EMCNNCheck::beginJob(){

  art::ServiceHandle<art::TFileService> fileServiceHandle;
  ftree = fileServiceHandle->make<TTree>("ftree", "hit info");
  ftree->Branch("run", &run, "run/I");
  ftree->Branch("event", &event, "event/I");
  ftree->Branch("beampdg", &beampdg, "beampdg/I");
  ftree->Branch("average_score_em" , &average_score_em , "average_score_em/D");
  ftree->Branch("average_score_trk", &average_score_trk, "average_score_trk/D");
  ftree->Branch("average_score_mic", &average_score_mic, "average_score_mic/D");
  ftree->Branch("channel", &channel);
  ftree->Branch("tpc", &tpc);
  ftree->Branch("plane", &plane);
  ftree->Branch("wire", &wire);
  ftree->Branch("charge", &charge);
  ftree->Branch("peakt", &peakt);
  ftree->Branch("score_em", &score_em);
  ftree->Branch("score_trk", &score_trk);
  ftree->Branch("score_mic", &score_mic);
  ftree->Branch("pdg", &pdg);
  ftree->Branch("origin", &origin);

}


DEFINE_ART_MODULE(pdsp::EMCNNCheck)
