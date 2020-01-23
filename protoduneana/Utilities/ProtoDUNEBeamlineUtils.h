#ifndef PROTODUNE_BEAMLINE_UTILS_H
#define PROTODUNE_BEAMLINE_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNEBeamlineUtils
//  - Class to create tracks and momentum from beamline 
//    information. Will be used when updated monitor positions
//    are available.
//   
// - Jake Calcutt (calcuttj@msu.edu) 
///////////////////////////////////////////////////////////////

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
#include "lardataobj/RecoBase/Track.h"

#include <set>
#include <map>
#include <string>
#include <ostream>

#include "TVector3.h"

namespace protoana {

  enum beamPDG{
    kElectron = 11,
    kMuon = 13,
    kPion = 211,
    kKaon = 321,
    kProton = 2212,
    kDeuteron = 1000010020
  };

  struct PossibleParticleCands {
    public:
      bool electron = false;
      bool muon = false;
      bool pion = false;
      bool kaon = false;
      bool proton = false;
      bool deuteron = false;

    inline PossibleParticleCands operator&&(const PossibleParticleCands & b) const {
        return {electron && b.electron, muon && b.muon, pion && b.pion, kaon && b.kaon, proton && b.proton, deuteron && b.deuteron};
    }
    inline PossibleParticleCands operator||(const PossibleParticleCands & b) const {
        return {electron || b.electron, muon || b.muon, pion || b.pion, kaon || b.kaon, proton || b.proton, deuteron || b.deuteron};
    }
    inline operator std::string () const { // overload cast to string
        std::string result = "PossibleParticleCands: [ ";
        if (electron) result += "e ";
        if (muon) result += "mu ";
        if (pion) result += "pi ";
        if (kaon) result += "k ";
        if (proton) result += "p ";
        if (deuteron) result += "e ";
        result += "]";
        return result;
    }
    inline std::vector<int> getPDGCodes() const {
        std::vector<int> result;
        if (electron) result.push_back(kElectron);
        if (muon) result.push_back(kMuon);
        if (pion) result.push_back(kPion);
        if (kaon) result.push_back(kKaon);
        if (proton) result.push_back(kProton);
        if (deuteron) result.push_back(kDeuteron);
        return result;
    }
  };

  class ProtoDUNEBeamlineUtils {

  public:

    ProtoDUNEBeamlineUtils(fhicl::ParameterSet const &pset);
    ~ProtoDUNEBeamlineUtils();

    void reconfigure(fhicl::ParameterSet const &pset);


    //GetBeamEvent(art::Event const & evt, std::string beamEventTag) 
    const beam::ProtoDUNEBeamEvent GetBeamEvent(art::Event const & evt) const;

    //GetFibers(art::Event const & evt, std::string beamEventTag, std::string monitorName)
    //Replace with return double
    void GetFibers( art::Event const & evt); 
    //GetCurrent(art::Event const & evt, std::string beamEventTag)
    //Replace with return double
    void GetCurrent( art::Event const & evt);

    //GetTracks()
    //GetMomenta()

    //Remove
    std::vector< recob::Track> MakeTracks( art::Event const & evt);

    //Remove
    double   GetPosition( short );
    //Remove
    TVector3 ConvertMonitorCoordinates( double, double, double, double );
    //Remove
    void     BeamMonitorBasisVectors();
    //Remove
    void     RotateMonitorVector(TVector3&);
    //Remove
    TVector3 ProjectToTPC(TVector3, TVector3);

    //Remove
    std::vector< double > MomentumSpec( art::Event const & evt);
    //Remove
    double MomentumCosTheta( double, double, double );

    /**
     * Get the particle ID from beamline instrumentation info (returns PDG IDs)
     *
     * Make sure to set fcl param UseCERNCalibSelection=true if the beam reco is done after ~v08_07_00
     * otherwise if using an old reco'd BeamEvent (e.g. you are using the first central reco) set it to false
     */
    std::vector< int > GetPID( beam::ProtoDUNEBeamEvent const & beamevt, double nominal_momentum );
    std::vector< int > GetPID( art::Event const & evt, double nominal_momentum );

    /**
     * Get the particle ID from beamline instrumentation info
     *
     * Make sure to set fcl param UseCERNCalibSelection=true if the beam reco is done after ~v08_07_00
     * otherwise if using an old reco'd BeamEvent (e.g. you are using the first central reco) set it to false
     */
    PossibleParticleCands GetPIDCandidates( beam::ProtoDUNEBeamEvent const & beamevt, double nominal_momentum );
    PossibleParticleCands GetPIDCandidates( art::Event const & beamevt, double nominal_momentum );

    /**
     * Compute what the beamline momentum SHOULD BE given a 
     * certain pdg and beamline TOF
     */
    double ComputeMomentum( int pdg, double tof );
    /**
     * Compute what the TOF SHOULD BE given a certain pdg 
     * and beamline momentum
     */
    double ComputeTOF     ( int pdg, double momentum );

    /**
     * Returns true if the beamline instrumentation has a good trigger
     * that matches the ProtoDUNE trigger.
     */
    bool IsGoodBeamlineTrigger(art::Event const & evt) const;
    bool IsGoodBeamlineTrigger(beam::ProtoDUNEBeamEvent const & beamEvent) const;

    /**
     * Returns true if only 1 fiber is activated in each of the 
     * momentum spectrometer monitors
     */
    bool HasPerfectBeamMomentum(art::Event const & evt) const;
    bool HasPerfectBeamMomentum(beam::ProtoDUNEBeamEvent const & beamEvent) const;

    /**
     * Uses the beamline momentum and time of flight measurements
     * to estimate the beamline particle mass in GeV/c^2
     *
     * Returns a list of all possible cominations (there can be 
     * multiple momenta).
     */
    std::vector<double> GetBeamlineMass(art::Event const & evt) const;

    /**
     * Uses the beamline momentum and time of flight measurements
     * to estimate the beamline particle mass^2 in (GeV/c^2)^2
     *
     * This can be more useful than the mass since measurement
     * fluctuations can make it negative.
     *
     * Returns a list of all possible cominations (there can be 
     * multiple momenta).
     */
    std::vector<double> GetBeamlineMassSquared(art::Event const & evt) const;

    /**
     *  Get reconstructed beamline momentum (in GeV/c), tof (in ns), and
     *  flags for if the ckov's fired. Will be < 0 if invalid.
     *
     *  C++ structured binding, so call like:
     *  const auto [momentum, tof, ckov0, ckov1] = dataUtils.GetBeamlineInformation(e);
     *
     *  then you have the normal float momentum, int ckov0 variables, etc. in the current scope.
     */
    const std::tuple<double,double,int,int> GetBeamlineVars(art::Event const & evt) const;

    /**
     *  Get reconstructed beamline momentum (in GeV/c), tof (in ns),
     *  flags for the tofChannel and ckov's. Also the timing trigger (12 means beam),
     *  the beam instrumentation trigger (1 means beam trigger), and whether
     *  the timing and beam instrumentation triggers are matched.
     *
     *  All values will be < 0 if invalid.
     *
     *  C++ structured binding, so call like:
     *  const auto [momentum, tof, tofChannel,ckov0,ckov1,ckov0Pressure,ckov1Pressure,timingTrigger,BITrigger,areBIAndTimingMatched] = dataUtils.GetBeamlineInformation(e);
     *
     *  then you have the normal float momentum, int ckov0 variables, etc. in the current scope. You can use (void) variable; lines
     *  to get rid of unused var warnings.
     */
    const std::tuple<double,double,int,int,int,double,double,int,int,bool> GetBeamlineVarsAndStatus(art::Event const & evt) const;

  private:

    /**
     * Get the particle ID from beamline instrumentation info (this uses the official CERN cuts and assumes
     * CERN calibrations)
     */
    std::vector< int > GetPID_CERNCalib( beam::ProtoDUNEBeamEvent const & beamevt, double nominal_momentum );
    /**
     * Get the particle ID from beamline instrumentation info (this uses the official CERN cuts and assumes
     * CERN calibrations) if fUseCERNCalibSelection is set, otherwise uses old values.
     */
    PossibleParticleCands GetPIDCandidates_CERNCalib( beam::ProtoDUNEBeamEvent const & beamevt, double nominal_momentum );

    art::InputTag fBeamEventTag;

    bool fUseCERNCalibSelection;

    //Remove
    std::map< std::string, std::vector< short > > ActiveFibers;

    //Remove
    double Current;
    
    //std::vector< recob::Track > Tracks;

    //Just to make it easier to go through all devices
    std::vector< std::string > AllDevices = { 
      "XBPF022707",
      "XBPF022708",
      "XBPF022716",
      "XBPF022717",
      "XBPF022697",
      "XBPF022701",
      "XBPF022702"    
    };

    //These will be static
    //
    //Tracking monitors
    std::string HorizUpstream   =  "XBPF022707";
    std::string VertUpstream    =  "XBPF022708";
    std::string HorizDownstream =  "XBPF022716";
    std::string VertDownstream  =  "XBPF022717";

    //Momentum monitors
    std::string BProf1 =  "XBPF022697";
    std::string BProf2 =  "XBPF022701";
    std::string BProf3 =  "XBPF022702";


    //Basis Vectors for transforming coords
    TVector3 MonitorBasisX;
    TVector3 MonitorBasisY;
    TVector3 MonitorBasisZ;

    bool rotated = false;

    double fRotateMonitorXZ;
    double fRotateMonitorYZ;
    double fRotateMonitorYX;

    double fFirstTrackingProfZ;
    double fSecondTrackingProfZ;
    double fNP04FrontZ;  
    double fBeamX, fBeamY, fBeamZ;


    double fBeamBend;
    double L1, L2, L3;
    double fBProf1Shift, fBProf2Shift, fBProf3Shift;

    //Hardware Parameters for magnetic field stuff
    double mag_P1 = 5.82044830e-3;
//    double mag_P2 = 0.;
    double mag_P3 = -4.68880000e-6;
    double mag_P4 = 324.573967;

    std::map< int, double > particle_mass = {
       {kElectron, .0005109989},
       {kMuon,     .1056583745},
       {kPion,     .13957018  },
       {kKaon,     .493677    },
       {kProton,   .9382720813},
       {kDeuteron, 2.013553213}
    };

    double c = 299792458.; // m/s
    const double fTOFDist = 28.575;  // m

    //////////////////////////////////////////////////////////////////
    // for Justin Hugon's old beamline selection from BeamlineUtils //
    //////////////////////////////////////////////////////////////////
    float fMomentumScaleFactor;
    float fMomentumOffset; // GeV/c
    float fTOFScaleFactor;
    float fTOFOffset; // ns
    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
  };

}

#endif

