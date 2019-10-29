#ifndef PROTODUNE_TRUTH_UTILS_H
#define PROTODUNE_TRUTH_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNETruthUtils
//  - Class to help analysers access useful truth information
//
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////


#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "art/Framework/Principal/Event.h"

namespace protoana {

  struct MCParticleSharedHits {
    const simb::MCParticle * particle;
    size_t nSharedHits;
    size_t nSharedDeltaRayHits;
  };

  class ProtoDUNETruthUtils {

  public:

    ProtoDUNETruthUtils();
    ~ProtoDUNETruthUtils();

    // Get shared hits between reco objects and MCParticle.
    
    std::vector< const recob::Hit * > FillSharedHits( const simb::MCParticle & mcpart, 
      const std::vector< const recob::Hit * > & hitsVec, bool delta_ray ) const;

    std::vector<const recob::Hit*> GetSharedHits(const simb::MCParticle &mcpart,
      const recob::PFParticle &pfpart, const art::Event &evt, std::string pfparticleModule, bool delta_ray = false) const;
    std::vector<const recob::Hit*> GetSharedHits(const simb::MCParticle &mcpart,
      const recob::Track &track, const art::Event &evt, std::string trackModule, bool delta_ray = false) const;
    std::vector<const recob::Hit*> GetSharedHits(const simb::MCParticle &mcpart,
      const recob::Shower &shower, const art::Event &evt, std::string showerModule, bool delta_ray = false) const;

    // Get hits associated with an MCParticle.
    std::vector<const recob::Hit*> GetMCParticleHits(const simb::MCParticle &mcpart,
      const art::Event &evt, std::string hitModule) const;
    
    // Get completeness and purity of reconstructed objects.
    template <typename T>
    double GetCompleteness(const T &recobj, const art::Event &evt,
      std::string recoModule, std::string hitModule) const;

    double GetPurity(const recob::PFParticle &pfpart, const art::Event &evt,
      std::string pfparticleModule) const;
    double GetPurity(const recob::Track &track, const art::Event &evt,
      std::string trackModule) const;
    double GetPurity(const recob::Shower &shower, const art::Event &evt,
      std::string showerModule) const;

    // Get MCParticle list from a hit vector depending on whether the hits came from a shower or track.
    std::vector<std::pair<const simb::MCParticle*, double>> GetMCParticleListFromTrackHits
      (const std::vector<const recob::Hit*>& hitVec) const;
    std::vector<std::pair<const simb::MCParticle*, double>> GetMCParticleListFromShowerHits
      (const std::vector<const recob::Hit*>& hitVec) const;

    // Contributions of MCParticles to PFParticles and vice versa
    std::vector<std::pair<const simb::MCParticle*, double>> GetMCParticleListFromPFParticle
      (const recob::PFParticle &pfpart, art::Event const &evt, std::string pfparticleModule) const;
    std::vector<std::pair<const recob::PFParticle*, double>> GetPFParticleListFromMCParticle
      (const simb::MCParticle &part, art::Event const &evt, std::string pfparticleModule) const;

    // Contributions of MCParticles to tracks and vice versa
    std::vector<std::pair<const simb::MCParticle*, double>> GetMCParticleListFromRecoTrack
      (const recob::Track &track, art::Event const &evt, std::string trackModule) const;
    std::vector<std::pair<const recob::Track*, double>> GetRecoTrackListFromMCParticle
      (const simb::MCParticle &part, art::Event const &evt, std::string trackModule) const;

    // Contributions of MCParticles to showers and vice versa
    std::vector<std::pair<const simb::MCParticle*, double>> GetMCParticleListFromRecoShower
      (const recob::Shower &shower, art::Event const &evt, std::string showerModule) const;
    std::vector<std::pair<const recob::Shower*, double>> GetRecoShowerListFromMCParticle
      (const simb::MCParticle &part, art::Event const &evt, std::string showerModule) const;

    // Best match getters between MC and reconstruction
    const simb::MCParticle* GetMCParticleFromPFParticle
      (const recob::PFParticle &pfpart, art::Event const &evt, std::string pfparticleModule) const;
    const recob::PFParticle* GetPFParticleFromMCParticle
      (const simb::MCParticle &part, art::Event const &evt, std::string pfparticleModule) const;
    const simb::MCParticle* GetMCParticleFromRecoTrack
      (const recob::Track &track, art::Event const &evt, std::string trackModule) const;
    const recob::Track* GetRecoTrackFromMCParticle
      (const simb::MCParticle &part, art::Event const &evt, std::string trackModule) const;
    const simb::MCParticle* GetMCParticleFromRecoShower
      (const recob::Shower &shower, art::Event const &evt, std::string showerModule) const;
    const recob::Shower* GetRecoShowerFromMCParticle
      (const simb::MCParticle &part, art::Event const &evt, std::string showerModule) const;

    // General overloaded MCParticle match getters
    const simb::MCParticle* GetMCParticleFromReco
      (const recob::PFParticle &pfpart, art::Event const &evt, std::string pfparticleModule) const;
    const simb::MCParticle* GetMCParticleFromReco
      (const recob::Track &track, art::Event const &evt, std::string trackModule) const;
    const simb::MCParticle* GetMCParticleFromReco
      (const recob::Shower &shower, art::Event const &evt, std::string showerModule) const;

    std::vector<std::pair<const simb::MCParticle*, double>> GetMCParticleListFromReco
      (const recob::PFParticle &pfpart, art::Event const &evt, std::string pfparticleModule) const;
    std::vector<std::pair<const simb::MCParticle*, double>> GetMCParticleListFromReco
      (const recob::Track &track, art::Event const &evt, std::string trackModule) const;
    std::vector<std::pair<const simb::MCParticle*, double>> GetMCParticleListFromReco
      (const recob::Shower &shower, art::Event const &evt, std::string showerModule) const;

    // Match Reco to True by Number of Hits contributed as opposed to Energy contributed
    template <typename T>
    //std::vector< std::pair< const simb::MCParticle*, size_t > > GetMCParticleListByHits
    std::vector< MCParticleSharedHits > GetMCParticleListByHits
      (const T &recobj, const art::Event &evt, std::string recoModule, std::string hitModule) const;

    template <typename T>
    //const simb::MCParticle * GetMCParticleByHits
    const MCParticleSharedHits GetMCParticleByHits
      (const T &recobj, const art::Event &evt, std::string recoModule, std::string hitModule) const;


    const simb::MCParticle* MatchPduneMCtoG4( const simb::MCParticle & pDunePart, const art::Event & evt );
    const simb::MCParticle* GetGeantGoodParticle(const simb::MCTruth &genTruth, const art::Event &evt) const;

    // Converting times in LArSoft can be a bit of a minefield. These functions convert true times in ns
    // to pandora times in ns
    const float ConvertTrueTimeToPandoraTimeNano(const simb::MCParticle &part) const;
    const float ConvertTrueTimeToPandoraTimeNano(const float trueTime) const;
    // Microsecond versions
    const float ConvertTrueTimeToPandoraTimeMicro(const simb::MCParticle &part) const;
    const float ConvertTrueTimeToPandoraTimeMicro(const float trueTime) const;

    // Get interaction process key.
    int GetProcessKey(std::string process);

    // Get the MC truth deposited energy
    double GetDepEnergyMC(const art::Event &evt, geo::GeometryCore const * fGeom, int trackid, int whichview) const;

    // Get first trajectory point in TPC active volume
    int GetFirstTrajectoryPointInTPCActiveVolume(const simb::MCParticle& mcpart, double tpcactiveXlow, double tpcactiveXhigh, double tpcactiveYlow, double tpcactiveYhigh, double tpcactiveZlow, double tpcactiveZhigh);

    // Get MC Particle length in TPC active volume
    double GetMCParticleLengthInTPCActiveVolume(const simb::MCParticle& mcpart, double tpcactiveXlow, double tpcactiveXhigh, double tpcactiveYlow, double tpcactiveYhigh, double tpcactiveZlow, double tpcactiveZhigh);

    // Estimate last point energy loss
    // By default the kinetic energy at the last trajectory point is zero. Estimate the energy loss of the last trajectory by using the average energy loss of the track in the TPC active volume
    double GetDepEnergyAtLastTrajPoint(const simb::MCParticle& mcpart, double tpcactiveXlow, double tpcactiveXhigh, double tpcactiveYlow, double tpcactiveYhigh, double tpcactiveZlow, double tpcactiveZhigh);

    // Get particle's kinetic energy at the interaction vertex
    // It takes the particle's energy at the next to the final trajectory point and subtracts the estimated kinetic energy of the last trajectory point (by default is zero)
    double GetKinEnergyAtVertex(const simb::MCParticle& mcpart, double kinene_lastpoint=0.0);


    // Get the sim::IDEs from the MCParticle, organized by the trajectory points
    std::map< size_t, std::vector< const sim::IDE * > > GetSimIDEs( const simb::MCParticle & mcpart );

  private:


  };

  //Define this it the source file, but outside of class
  bool sort_IDEs( const sim::IDE * i1, const sim::IDE * i2)/*{
    return( i1->z < i2->z ); 
  }*/;

}

#endif
