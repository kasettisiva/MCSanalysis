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

// using namespace recob::tracking;

namespace trkf {

  class MCSMomentumCalculator {
  public:
    // Declare local public types
    struct Config;
    struct Result;

    // Type Definitions
    // using Parameters = fhicl::Table<Config>;
    using MCSAngleResult = MCSAngleCalculator::MCSAngleResult;
    using MCSSegmentResult = MCSAngleCalculator::MCSSegmentResult;
    using AnglePair_t = std::pair<Float_t, Float_t>; // (theta_XZ', theta_YZ')
    using Angles_t = MCSAngleCalculator::Angles_t; // (theta_3D, theta_XZ', theta_YZ')
    using Point_t = MCSAngleCalculator::Point_t;

    // Define fhicl Configuration
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
    explicit MCSMomentumCalculator() : 
      MCSMomentumCalculator(fhicl::Table<Config>(fhicl::ParameterSet())) 
    {}
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
    // Primary GetResult methods
    // TODO: Clean up this documentation
    // All we technically need is a positions vector
    // We can have our code parse either a recob::Track, recob::TrackTrajectory, recob::Trajectory, simb::MCParticle, or simb::MCTrajectory to just form a std::vector<Point_t> object that is passed to get everything we need.
    Result GetResult(const MCSAngleResult angleResult, Int_t angleMethod) const;
    Result GetResult(const MCSSegmentResult segmentResult, const int angleMethod = gDefaultMethod) const;
    Result GetResult(const std::vector<Point_t> points, const int angleMethod = gDefaultMethod) const;
    // Start Offet is the number of trajectory points to ignore at the start (i.e. they're outside of the detector)
    // End Offset is which trajectory point to stop at (if it leaves the detector for instance), default to -1, which will tell it to just use all of them. We could have an issue if the number of points is less than the endOffset, will have to add an assertion
    Result GetResult(const simb::MCTrajectory trajectory, const int angleMethod = gDefaultMethod, size_t startOffset = 0, size_t endOffset = -1);

    // Secondary GetResult methods that call other GetResult methods
    inline Result GetResult(const recob::Trajectory trajectory, const int angleMethod = gDefaultMethod) const {
      std::vector<Point_t> points;
      for(recob::tracking::Point_t position: trajectory.Positions())
	points.push_back(Point_t(position));
      return GetResult(points, angleMethod);
    }
    inline Result GetResult(const recob::TrackTrajectory trackTrajectory, const int angleMethod = gDefaultMethod) const {
      return GetResult(trackTrajectory.Trajectory(), angleMethod);
    }
    inline Result GetResult(const recob::Track track, const int angleMethod = gDefaultMethod) const {
      return GetResult(track.Trajectory(), angleMethod);
    }
    inline Result GetResult(const simb::MCParticle particle, const int angleMethod = gDefaultMethod, size_t startOffset = 0, size_t endOffset = -1) {
      return GetResult(particle.Trajectory(), angleMethod, startOffset, endOffset);
    }

    // Define Result
    struct Result {
    public:
      explicit Result(int angleMethod, MCSAngleResult mcsAngleResult, std::vector<Float_t> segmentLength_vec, std::vector<Float_t> logLikelihood_vec, std::vector<Float_t> MCSMomentum_vec) :
	fAngleMethod(angleMethod),
	fAngleResult(mcsAngleResult),
	fSegmentLength_vec(segmentLength_vec),
	fLogLikelihood_vec(logLikelihood_vec),
	fMCSMomentum_vec(MCSMomentum_vec)
      {}
      // Will correspond to either the track segment lengths or polygonal lengths, depending on the momentum estimation method.
      inline int GetAngleMethod() const { return fAngleMethod; }
      inline MCSAngleResult GetAngleResult() const { return fAngleResult; }
      inline std::vector<Float_t> GetSegmentLength_vec() const { return fSegmentLength_vec; }
      inline std::vector<Float_t> GetLogLikelihood_vec() const { return fLogLikelihood_vec; }
      inline Float_t GetMCSMomentum(int i = 0) const { return fMCSMomentum_vec.at(i); }
    private:
      const int fAngleMethod;
      const MCSAngleResult fAngleResult;
      const std::vector<Float_t> fSegmentLength_vec;
      const std::vector<Float_t> fLogLikelihood_vec;
      const std::vector<Float_t> fMCSMomentum_vec;
    };

    // Public tests
    // TODO: Move this to a test file.
    void RunTest();

  private:
    // Private properties
    const MCSAngleCalculator fAngleCalculator;
    double fSigmaRES; // This cannot be const so that it can be reset during true MCS calculations
    const double fMomentumMinimum;
    const double fMomentumMaximum;
    const double fMomentumStep;
    const double fParticleMass;
    const double fKappa_a;
    const double fKappa_c;

    // Static private properties (used as defaults for function parameters)
    static const int gDefaultMethod = 0;
    
    // Private Set Functions
    // Momentum Estimation
    std::vector<Float_t> SetLogLikelihood_vec(std::vector<AnglePair_t> anglePair_vec, std::vector<Float_t> segmentLength_vec) const;
    std::vector<Float_t> SetLogLikelihood_vec(std::vector<Float_t> rawAngle_vec, std::vector<Float_t> segmentLength_vec) const;
    std::vector<Float_t> SetMCSMomentum_vec(std::vector<Float_t> logLikelihood_vec, std::vector<Float_t> segmentLength_vec) const;
    MCSAngleResult SetAngleResultForAngleMethod(const std::vector<Point_t> trajPoint_vec, const Int_t angleMethod) const;
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
