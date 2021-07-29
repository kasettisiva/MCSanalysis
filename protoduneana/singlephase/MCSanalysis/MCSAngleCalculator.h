// -*- mode: c++ -*-
#ifndef MCSANGLECALCULATOR_H
#define MCSANGLECALCULATOR_H

// Framework Includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Trajectory.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "MCSSegmentCalculator.h"

// Root includes
#include "TVector3.h"
#include "TMath.h"
#include "TRandom3.h"

// C++ includes
#include <iostream>

namespace trkf {
  class MCSAngleCalculator {
  public:
    // Declare local public types
    struct Config;
    struct MCSAngleResult;

    // Type Defintions
    using MCSSegmentResult = MCSSegmentCalculator::MCSSegmentResult;
    using Barycenter_t = recob::tracking::Point_t;
    using Angles_t = std::tuple<Float_t, Float_t, Float_t>; // (theta_3D, theta_XZ', theta_YZ')
    //using Point_t = recob::tracking::Point_t;
    using Point_t = MCSSegmentCalculator::Point_t;
    using Vector_t = recob::tracking::Vector_t;

    // Define fhicl Configuration
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<MCSSegmentCalculator::Config> segmentCalculator { Name("segmentCalculator") };
    };

    // Public Constructors
    explicit MCSAngleCalculator() : 
      MCSAngleCalculator(fhicl::Table<Config>(fhicl::ParameterSet())) 
    {}
    explicit MCSAngleCalculator(const fhicl::Table<Config> t) : 
      MCSAngleCalculator(
			 t().segmentCalculator
			 ) 
    {}
    explicit MCSAngleCalculator(const fhicl::Table<MCSSegmentCalculator::Config> segmentCalculatorConfig) : 
      fSegmentCalculator(MCSSegmentCalculator(segmentCalculatorConfig))
    {}

    // Public Methods
    // Primary GetResult methods
    // TODO: Add a GetResult method that will also calculate the angles of the firstSegment relative to a vector that is passed.  
    //     For example, if you know the initial direction of your particle from beam-info, then there's no need to only use the first-segment as the initial direction.
    //     We should add a flag if this occured to the MCSAngleResult so that the MCSMomentumCalculator doesn't add back that first segment momentum.
    MCSAngleResult GetResult(const MCSSegmentResult segmentResult, const int angleMethod) const;
    MCSAngleResult GetResult(const std::vector<Point_t> points, const int angleMethod) const;
    MCSAngleResult GetResult(const simb::MCTrajectory trajectory, const int angleMethod) const;

    // Secondary GetResult methods that call other GetResult Methods
    inline MCSAngleResult GetResult(const simb::MCParticle particle, const int angleMethod) const {
      return GetResult(particle.Trajectory(), angleMethod);
    }
    inline MCSAngleResult GetResult(const recob::Trajectory trajectory, const int angleMethod) const {
      std::vector<Point_t> points;
      for(recob::tracking::Point_t position: trajectory.Positions())
	points.push_back(Point_t(position));
      return GetResult(points, angleMethod);
    }
    inline MCSAngleResult GetResult(const recob::TrackTrajectory trackTrajectory, const int angleMethod) const {
      return GetResult(trackTrajectory.Trajectory(), angleMethod);
    }
    inline MCSAngleResult GetResult(const recob::Track track, const int angleMethod) const {
      return GetResult(track.Trajectory(), angleMethod);
    }

    // TODO: Consider moving to lardataobj, similar to the TrajectoryMCSFitter's recob::MCSFitResult
    // MCSAngleResult
    struct MCSAngleResult {
    public:
      // Public MCSAngleResult Constructor
      explicit MCSAngleResult(
			      int segmentMethod,
			      MCSSegmentResult segmentResult,
			      std::vector<Float_t> segmentLength_vec,
			      std::vector<Vector_t> segmentVector_vec,
			      std::vector<Angles_t> angles_vec
			      ) :
	fSegmentMethod(segmentMethod),
	fSegmentResult(segmentResult),
	fSegmentLength_vec(segmentLength_vec),
	fSegmentVector_vec(segmentVector_vec),
	fAngles_vec(angles_vec)
      {}
      // Public MCSAngleResult Methods
      inline int GetSegmentMethod() const { return fSegmentMethod; }
      inline MCSSegmentResult GetSegmentResult() const { return fSegmentResult; }
      inline std::vector<Float_t> GetSegmentLength_vec() const { return fSegmentLength_vec; }
      inline std::vector<Vector_t> GetSegmentVector_vec() const { return fSegmentVector_vec; }
      inline std::vector<Angles_t> GetAngles_vec() const { return fAngles_vec; }

    private:
      // Private MCSAngleResult 
      const int fSegmentMethod;
      const MCSSegmentResult fSegmentResult;
      const std::vector<Float_t> fSegmentLength_vec;
      const std::vector<Vector_t> fSegmentVector_vec;
      const std::vector<Angles_t> fAngles_vec;
    }; // struct MCSAngleResult

    /* Public Helper Functions */

    // Transform and return the currentSegment into a coordinate system defined by the previouSegment, where the z-axis will lie on the previousSegment.  The phi rotation is ambiguous.
    TVector3 transform(TVector3 currentSegment, TVector3 previousSegment) const;

    // Calculate and return the 3D angle and projected 2D angles of the currentSegment with respect to the previousSegment.
    Angles_t calculateAngles(TVector3 currentSegment, TVector3 previousSegment) const;

    // Determine if the provided generic Point is located in the detector.
    // Implementation Comment: MCSSegmentCalculator::inDetector is called.
    // template<typename Point> Bool_t inDetector(Point point) const;
    template<typename Point> inline Bool_t inDetector(Point point) const {
      return fSegmentCalculator.inDetector(point);
      // return 
      // 	point.Z() <= 700 && point.Z() >= 0 &&
      // 	point.X() <= 300  && point.X() >= -300 &&
      // 	point.Y() <= 600 && point.Y() >= 0;
    }

    // Convert the vector from the `From` type to the `To` type.
    // Implementation Comment: MCSSegmentCalculator::convert is called.
    // template<typename From, typename To> To convert(From vector) const;
    template<typename From, typename To> inline To convert(From vector) const {
      return fSegmentCalculator.convert<From, To>(vector);
      // return To(vector.X(), vector.Y(), vector.Z());
    }
    
  private:
    // Private Properties
    const MCSSegmentCalculator fSegmentCalculator;

    // Calculate and return the vector of segment angles (3D angle & 2D projected angles) between the provided segment vectors.
    std::vector<Angles_t> SetAngles_vec(std::vector<Vector_t> segmentVector_vec) const;
  }; // class MCSAngleCalculator
} // namespace trkf

#endif

/* // TODO: This may be the ideal future: */
/* struct TrackSegment { */
/* private: */
/*   int fStartPoint; // Either this or change Point_t to some custom Type */
/*   std::vector<recob::Tracking::Point_t> fTrajectoryPoints; */
/*   Point_t fBarycenter; */
/*   double fSegmentLength; */
/* }; */
