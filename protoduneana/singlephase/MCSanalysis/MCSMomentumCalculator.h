// -*- mode: c++ -*-
#ifndef MCSMOMENTUMCALCULATOR_H
#define MCSMOMENTUMCALCULATOR_H

// Framework includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Trajectory.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "MCSAngleCalculator.h"

// Root includes
#include "TVector3.h"
#include "TMath.h"
#include "TRandom3.h"

// C++ includes
#include <iostream>
#include <string>
#include <assert.h>

namespace trkf {
  /**
     @author Hunter C. Meyer (LSU, hmeyer5@lsu.edu)
     @file protoduneana/singlephase/MCSanalysis/MCSMomentumCalculator.h
     @class trkf::MCSMomentumCalculator MCSMomentumCalculator.h "protoduneana/singlephase/MCSanalysis/MCSMomentumCalculator.h"
     @brief A calculator that will estimate the momentum of a particle's track.
     
     For the fhicl configuration, see `trkf::MCSMomentumCalculator::Config`.  This momentum calculator should be not used concurrently, see `fSigmaRES`
     
     @todo Add a `GetResult` method that will accept a vector for the initial direction.  Details in `MCSAngleCalculator` documentation.
     @todo Add support for more types of particles than just muons by making the mass be an easy parameter.
     @todo Add support for additional types of energy loss, including dEdx measurements or other formulas besides Bethe-Bloch (Landau, for example)
     @todo Add a direction-dependent sigmaRES values.
     @todo Document type defintions
     
     @see trkf::MCSAngleCalculator
     @see trkf::MCSMomentumCalculator
   */
  class MCSMomentumCalculator {
  public:
    // Declare local public types
    struct Result;

    // Type Definitions
    using MCSAngleResult = MCSAngleCalculator::MCSAngleResult;
    using MCSSegmentResult = MCSAngleCalculator::MCSSegmentResult;
    using AnglePair_t = std::pair<Float_t, Float_t>; // (theta_XZ', theta_YZ')
    using Angles_t = MCSAngleCalculator::Angles_t; // (theta_3D, theta_XZ', theta_YZ')
    using Point_t = MCSAngleCalculator::Point_t;

    /// fhicl Configuration used for initialization of `trkf::MCSMomentumCalculator`
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<MCSAngleCalculator::Config> angleCalculator { Name("angleCalculator") };
      fhicl::Atom<double> sigmaRES { Name("sigmaRES"), Comment("Detector Angular Resolution (rad)"), 0.003 };
      fhicl::Atom<double> momentumMinimum { Name("momentumMinimum"), Comment("Minimum Momentum Estiamte (GeV/c)"), 0.1 };
      fhicl::Atom<double> momentumMaximum { Name("momentumMaximum"), Comment("Maximum Momentum Estimate (GeV/c)"), 7.5 };
      fhicl::Atom<double> momentumStep { Name("momentumStep"), Comment("Momentum Estimate Increment (GeV/c)"), 0.01 };
      fhicl::Atom<double> particleMass { Name("particleMass"), Comment("Predicted Mass of particle (MeV/c^2)"), 105.66 }; // TODO: Change to GeV
      fhicl::Atom<double> kappa_a { Name("kappa_a"), Comment("Fit parameter 'a' used in modified highland formula (MeV)"), 0.105 }; // TODO: Add reference for default value
      fhicl::Atom<double> kappa_c { Name("kappa_c"), Comment("Fit parameter 'c' used in modified highland formula (MeV)"), 11.004 }; // TODO: Add reference for default value
    };

    // Public Constructors
    /// Default Constructor, using an empty `fhicl::ParameterSet` for the fhicl configuration, leading to default values being used in `trkf::MCSMomentumCalculator::Config`
    explicit MCSMomentumCalculator() : 
      MCSMomentumCalculator(fhicl::Table<Config>(fhicl::ParameterSet())) 
    {}
    /// Constructor where you pass a `fhicl::Table<Config>`
    explicit MCSMomentumCalculator(const fhicl::Table<Config> t) :
      MCSMomentumCalculator(
			    t().angleCalculator,
			    t().sigmaRES(),
			    t().momentumMinimum(),
			    t().momentumMaximum(),
			    t().momentumStep(),
			    t().particleMass(),
			    t().kappa_a(),
			    t().kappa_c()
			    ) 
    {}
    /**
       Constructor where you can pass each parameter in `trkf::MCSMomentumCalculator::Config` explicitly.
       
       @details MLE = Maximum Likelihood Estimation
       
       @param angleCalculatorConfig The configuration used for the `MCSAngleCalculator`.  `fAngleCalculator` will only be used in certain `GetResult` methods.
       @param sigamRES The detector angular resolution used in the MLE.
       @param momentumMinimum The minimum momentum to use in the MLE.
       @param momentumMaximum The maximum momentum to use in the MLE.
       @param momentumStep The momentum stepsize to use in the MLE.
       @param particleMass The assumed mass of the particle.
       @param kappa_a One of the fit-parameters used in the modified-highland formula.
       @param kappa_c One of the fit-parameters used in the modified-highland formula.
     */
    explicit MCSMomentumCalculator(fhicl::Table<MCSAngleCalculator::Config> angleCalculatorConfig, double sigmaRES, double momentumMinimum, double momentumMaximum, double momentumStep, double particleMass, double kappa_a, double kappa_c) : 
      fAngleCalculator(MCSAngleCalculator(angleCalculatorConfig)), 
      fSigmaRES(sigmaRES), 
      fMomentumMinimum(momentumMinimum),
      fMomentumMaximum(momentumMaximum),
      fMomentumStep(momentumStep),
      fParticleMass(particleMass),
      fKappa_a(kappa_a),
      fKappa_c(kappa_c) 
    {}

    // Public Methods
    /**
       
     */
    Result GetResult(const MCSAngleResult angleResult, Int_t angleMethod) const;

    /**
       
     */
    Result GetResult(const MCSSegmentResult segmentResult, const int angleMethod = gDefaultMethod) const;

    /**
       
     */
    Result GetResult(const std::vector<Point_t> points, const int angleMethod = gDefaultMethod) const;
    // Start Offet is the number of trajectory points to ignore at the start (i.e. they're outside of the detector)
    // End Offset is which trajectory point to stop at (if it leaves the detector for instance), default to -1, which will tell it to just use all of them. We could have an issue if the number of points is less than the endOffset, will have to add an assertion

    /**
       
     */
    Result GetResult(const simb::MCTrajectory trajectory, const int angleMethod = gDefaultMethod, size_t startOffset = 0, size_t endOffset = -1);

    /**
       
     */
    inline Result GetResult(const recob::Trajectory trajectory, const int angleMethod = gDefaultMethod) const {
      std::vector<Point_t> points;
      for(recob::tracking::Point_t position: trajectory.Positions())
	points.push_back(Point_t(position));
      return GetResult(points, angleMethod);
    }

    /**
       
     */
    inline Result GetResult(const recob::TrackTrajectory trackTrajectory, const int angleMethod = gDefaultMethod) const {
      return GetResult(trackTrajectory.Trajectory(), angleMethod);
    }

    /**
       
     */
    inline Result GetResult(const recob::Track track, const int angleMethod = gDefaultMethod) const {
      return GetResult(track.Trajectory(), angleMethod);
    }

    /**
       
     */
    inline Result GetResult(const simb::MCParticle particle, const int angleMethod = gDefaultMethod, size_t startOffset = 0, size_t endOffset = -1) {
      return GetResult(particle.Trajectory(), angleMethod, startOffset, endOffset);
    }

    // Define Result
    /**
       @struct trkf::MCSMomentumCalculator::Result MCSMomentumCalculator.h "protoduneana/singlephase/MCSanalysis/MCSMomentumCalculator.h"
       
       @brief A container for momentum information, including the angle result used to create this result.  Depending on the angle method that is given in the appropriate `GetResult` method, the appropriate segment lengths and angles will be used to create the vector of log-likelihoods, which is used to estimate the MCS Momentum.
       
       @todo Consider moving to lardataobj instead of larreco, as was done in the `TrajectoryMCSFitter`.
       @todo Document the relative lengths of vectors.
       @todo Momentum Uncertainty
       @todo Calculate the MCS momentum for every trajectory point.
       @todo Change the log-likelihood vector to instead be a vector of std::pairs
       @todo Rename to `MCSMomentumResult`
     */
    struct Result {
    public:
      /// Explicitly create a `Result`
      explicit Result(int angleMethod, MCSAngleResult mcsAngleResult, std::vector<Float_t> segmentLength_vec, std::vector<Float_t> logLikelihood_vec, std::vector<Float_t> MCSMomentum_vec) :
	fAngleMethod(angleMethod),
	fAngleResult(mcsAngleResult),
	fSegmentLength_vec(segmentLength_vec),
	fLogLikelihood_vec(logLikelihood_vec),
	fMCSMomentum_vec(MCSMomentum_vec)
      {}
      /// Get the angle method that was used to create this momentum result.
      inline int GetAngleMethod() const { return fAngleMethod; }
      /// Get the momentum result's angle result.
      inline MCSAngleResult GetAngleResult() const { return fAngleResult; }
      /// Get the result's vector of segment lengths used in the momentum estimation.
      inline std::vector<Float_t> GetSegmentLength_vec() const { return fSegmentLength_vec; }
      /// Get the result's vector of log-likelihoods as a function of momentum that were used in the momentum estimation.
      inline std::vector<Float_t> GetLogLikelihood_vec() const { return fLogLikelihood_vec; }
      /// Get the vector of estimation MCS Momentum values for each segment.
      inline Float_t GetMCSMomentum(int i = 0) const { return fMCSMomentum_vec.at(i); }
    private:
      const int fAngleMethod; /// The angle method used to create this result.
      const MCSAngleResult fAngleResult; /// The angle result used to create this momentum result.
      const std::vector<Float_t> fSegmentLength_vec; /// The result's vector of segment lengths.
      const std::vector<Float_t> fLogLikelihood_vec; /// The result's vector of log-likelihoods
      const std::vector<Float_t> fMCSMomentum_vec; /// The result's vector of MCS Momentum values
    }; // struct Result

    // Public tests
    /// @todo Create MCS tests.
    void RunTest();

  private:
    // Private properties
    const MCSAngleCalculator fAngleCalculator; /// The angle calculator used in some `GetResult` methods.  If `GetResult(const MCSAngleResult angleREsult, const int segmentMethod)` is called, this will nobe used. Otherwise, it will be.

    double fSigmaRES; /// The sigmaRES value to be used in the momentum estimation. This cannot be const because it is temporarily set to 0 if the true MCS Momentum estimation is used.  Due to this mutable state, this momentum calculator should not be used concurrently.
    const double fMomentumMinimum; /// The minimum momentum value to be used in the MLE.
    const double fMomentumMaximum; /// The maximum momentum value to be used in the MLE.
    const double fMomentumStep; /// The momentum step size to use in the MLE.
    const double fParticleMass; /// The assumed particle mass.
    const double fKappa_a; /// One of the modified-highland formula fit-parameters.
    const double fKappa_c; /// One of the modified-highland formula fit-parameters.

    /// The default angle method to use.
    static const int gDefaultMethod = 0;
    
    /// Calculate the vector of log-likelihoods as a function of momentum (by index) using the projected angle values.
    std::vector<Float_t> SetLogLikelihood_vec(std::vector<AnglePair_t> anglePair_vec, std::vector<Float_t> segmentLength_vec) const;
    /// Calculate the vector of log-likelihoods as a function of momentum (by index) using the 3D angle values.
    /// @todo Rename rawAngle_vec to theta3D_vec
    std::vector<Float_t> SetLogLikelihood_vec(std::vector<Float_t> rawAngle_vec, std::vector<Float_t> segmentLength_vec) const;
    /// Calculate the vector of estimated MCS Momentum for each segment.
    std::vector<Float_t> SetMCSMomentum_vec(std::vector<Float_t> logLikelihood_vec, std::vector<Float_t> segmentLength_vec) const;
    /// Calculate the angle result for the trajPoint_vec, parsing the angle method to use the correct segment method.
    MCSAngleResult SetAngleResultForAngleMethod(const std::vector<Point_t> trajPoint_vec, const Int_t angleMethod) const;
    /// Calculate the angle result for the segmentResult, parsing the angle method to use the ocrrect segment method.
    MCSAngleResult SetAngleResultForAngleMethod(const MCSSegmentResult segmentResult, const Int_t angleMethod) const;

    // Functions used by other functions
    template<typename Point> Bool_t inDetector(Point point) const;

    // TODO: Remove these two once the tests have been moved, since that functionality is in MCSAC
    TVector3 transform(TVector3 currentSegment, TVector3 previousSegment);
    AnglePair_t calculateAnglePair(TVector3 currentSegment, TVector3 previousSegment);
    std::vector<AnglePair_t> ExtractProjectedAnglePair_vec(std::vector<Angles_t> angles_vec) const;
    std::vector<Float_t> ExtractTheta3D_vec(std::vector<Angles_t> angles_vec) const;
    double highland(double p, double l, double particleMass) const;
    double deltaP(double p, double l, double particleMass) const;

    // Private Tests
    // TODO: Move to separate test file
    Bool_t test_transform();
    Bool_t test_calculateAnglePair();
  }; // class MCSMomentumCalculator
} // namespace trkf

#endif
