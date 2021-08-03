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
  /**
     @author Hunter C. Meyer (LSU, hmeyer5@lsu.edu)
     @file protoduneana/singlephase/MCSanalysis/MCSAngleCaluclator.h
     @class trkf::MCSAngleCalculator MCSAngleCalculator.h "protoduneana/singlephase/MCSanalysis/MCSAngleCalculator.h"
     @brief A calculator that will measure angles between segments of a true or reconstructed trajectory.
     
     For the fhicl configuration, see `trkf::MCSAngleCalculator::Config`
     
     @todo Add a `GetResult` method that will accept a vector for the initial direction. This will be useful if you know the initial direction of your particle from beam-fino, then the first segment is not used as the initial direction.  A flag should be added to `MCSAngleResult` so that the `MCSMomentumCalculator` doesn't add back that first segment momentum.
     
     @see trkf::MCSSegmentCalculator
     @see trkf::MCSMomentumCalculator
     @see trkf::MCSAngleCalculator::Config
   */
  class MCSAngleCalculator {
  public:
    // Declare local public types
    struct Config;
    struct MCSAngleResult;

    // Type Defintions
    using MCSSegmentResult = MCSSegmentCalculator::MCSSegmentResult;
    using Barycenter_t = MCSSegmentCalculator::Barycenter_t;
    using Point_t = MCSSegmentCalculator::Point_t;
    using Vector_t = MCSSegmentCalculator::Vector_t;
    using Angles_t = std::tuple<Float_t, Float_t, Float_t>; // (theta_3D, theta_XZ', theta_YZ')

    /// fhicl Configuration used for initialization of `trkf::MCSAngleCalculator`
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<MCSSegmentCalculator::Config> segmentCalculator { Name("segmentCalculator") };
    }; // struct Config

    // Public Constructors
    /// Default Constructor, using an empty `fhicl::ParameterSet` for the fhicl configuration, leading to default values being used in `trkf::MCSAngleCalculator::Config`
    explicit MCSAngleCalculator() : 
      MCSAngleCalculator(fhicl::Table<Config>(fhicl::ParameterSet())) 
    {}
    /// Constructor where you pass a `fhicl::Table<Config>`
    explicit MCSAngleCalculator(const fhicl::Table<Config> t) : 
      MCSAngleCalculator(
			 t().segmentCalculator
			 ) 
    {}
    /// Constructor where you can pass each parameter in `trkf::MCSAngleCalculator::Config` explicity.
    /// @param segmentCalculatorConfig The configuration used for the `MCSSegmentCalculator`.  `fSegmentCalculator` will only be used for certain `GetResult` methods.
    explicit MCSAngleCalculator(const fhicl::Table<MCSSegmentCalculator::Config> segmentCalculatorConfig) : 
      fSegmentCalculator(MCSSegmentCalculator(segmentCalculatorConfig))
    {}

    // Public Methods

    /**
       Retrieve an `MCSAngleResult` using an `MCSSegmentResult`. This is the base method that all other `GetResult` methods call.
       
       @param segmentResult The segment result used to create the `MCSAngleResult`
       @param segmentMethod Determines which segments to use to calculate the angles.
           0: Angles between linear fit segments will be calculated.
	   1: Angles between polygonal segments will be calculated.
       @todo Add more segment methods, as they are added to the `MCSSegmentCalculator`
     */
    MCSAngleResult GetResult(const MCSSegmentResult segmentResult, const int segmentMethod) const;

    /**
       Retrieve an `MCSAngleResult` using a vector of points. 
       
       @details This will create an `MCSSegmentResult` and call `GetResult(const MCSSegmentResult segmentResult, const int segmentMethod)`.
       
       @param points The points used to create the `MCSSegmentResult`.
       @param segmentMethod See `GetResult(const MCSSegmentResult segmentResult, const int segmentMethod)` for segment method documentation.
       @see GetResult(const MCSSegmentResult segmentResult, const int segmentMethod)
    */
    MCSAngleResult GetResult(const std::vector<Point_t> points, const int segmentMethod) const;

    /**
       Retrieve an `MCSAngleResult` using a `simb::MCTrajectory`.
       
       @param points The points used to create the `MCSAngleResult`.
       @param segmentMethod See `GetResult(const MCSSegmentResult segmentResult, const int segmentMethod)` for segment method documentation.
       @see GetResult(const MCSSegmentResult segmentResult, const int segmentMethod)
    */
    MCSAngleResult GetResult(const simb::MCTrajectory trajectory, const int segmentMethod) const;

    /**
       Retrieve an `MCSAngleResult` using a `simb::MCParticle`.
       
       @param particle The particle used to create `MCSAngleResult`.
       @param segmentMethod See `GetResult(const MCSSegmentResult segmentResult, const int segmentMethod)` for segment method documentation.
       @see GetResult(const simb::MCTrajectory, const int segmentMethod)
     */
    inline MCSAngleResult GetResult(const simb::MCParticle particle, const int segmentMethod) const {
      return GetResult(particle.Trajectory(), segmentMethod);
    }

    /**
       Retrieve an `MCSAngleResult` using a `recob::Trajectory`.
       
       @details A vector of Points will be created by looping through `trajectory.Positions()`, each with an invalid momentum.
       
       @param trajectory The trajectory used to create the `MCSAngleResult`.
       @param segmentMethod See `GetResult(const MCSSegmentResult segmentResult, const int segmentMethod)` for segment method documentation.
       @see GetResult(const std::vector<Point_t> points, const int segmentMethod)
     */
    inline MCSAngleResult GetResult(const recob::Trajectory trajectory, const int segmentMethod) const {
      std::vector<Point_t> points;
      for(recob::tracking::Point_t position: trajectory.Positions())
	points.push_back(Point_t(position));
      return GetResult(points, segmentMethod);
    }

    /**
       Retrieve an `MCSAngleResult` using a `recob::TrackTrajectory`.
       
       @param trackTrajectory The trackTrajectory used to create the `MCSAngleResult`.
       @param segmentMethod See `GetResult(const MCSSegmentResult segmentResult, const int segmentMethod)` for segment method documentation.
       @see GetResult(const recob::Trajectory trajectory, const int segmentMethod)
       @todo Consider using `trackTrajectory.NextValidPoint` methodology, which would quire removing `GetResult(const recob::Trajectory trajectory, const int segmentMethod)`
     */
    inline MCSAngleResult GetResult(const recob::TrackTrajectory trackTrajectory, const int segmentMethod) const {
      return GetResult(trackTrajectory.Trajectory(), segmentMethod);
    }

    /**
       Retrieve an `MCSAngleResult` using a `recob::Track`.
       
       @param track The track used to create the `MCSAngleResult`.
       @param segmentMethod See `GetResult(const MCSSegmentResult segmentResult, const int segmentMethod)` for segment method documentation.
       @see GetResult(const recob::TrackTrajectory trackTrajectory, const int segmentMethod)
     */
    inline MCSAngleResult GetResult(const recob::Track track, const int segmentMethod) const {
      return GetResult(track.Trajectory(), segmentMethod);
    }

    /**
       @struct trkf::MCSAngleCalculator::MCSAngleResult MCSAngleCalculator.h "protoduneana/singlephase/MCSanalysis/MCSAngleCalculator.h"
       
       @brief A container for angle information, including the parsed vectors and segment lengths.  Depending on the segmentMethod that was given in the appropriate `GetResult` method, the appropriate segment vectors/segment lengths will be retrieved from an `MCSSegmentResult`.
       
       @todo 
       @todo Consider moving to lardataobj instead of larreco, as was done in the `TrajectoryMCSFitter`
       @todo Document the relative lengths of vectors.
       @todo Should the lengths just be stored in the length of vectors?
     */
    struct MCSAngleResult {
    public:
      /// Explicitly create an `MCSAngleResult`
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
      // Public Get Methods
      /// Get the segment method that was used to create this angle result.
      inline int GetSegmentMethod() const { return fSegmentMethod; }
      /// Get the angle result's segment result.
      inline MCSSegmentResult GetSegmentResult() const { return fSegmentResult; }
      /// Get the result's vector of segment length.
      inline std::vector<Float_t> GetSegmentLength_vec() const { return fSegmentLength_vec; }
      /// Get the result's vector of segment vectors.
      inline std::vector<Vector_t> GetSegmentVector_vec() const { return fSegmentVector_vec; }
      /// Get the result's vector of angles between segment vectors.
      inline std::vector<Angles_t> GetAngles_vec() const { return fAngles_vec; }

    private:
      // Private MCSAngleResult 
      const int fSegmentMethod; /// The segment method used to create this result.
      const MCSSegmentResult fSegmentResult; /// The segment result used to create this angle result.
      const std::vector<Float_t> fSegmentLength_vec; /// This result's vector of segment lengths.
      const std::vector<Vector_t> fSegmentVector_vec; /// This result's vector of segment vectors.
      const std::vector<Angles_t> fAngles_vec; /// This result's vector of angles.
    }; // struct MCSAngleResult

    /// Transform and return the currentSegment into a coordinate system defined by the previouSegment, where the z-axis will lie on the previousSegment.  The phi rotation is ambiguous.
    TVector3 transform(TVector3 currentSegment, TVector3 previousSegment) const;

    /// Calculate and return the 3D angle and projected 2D angles of the currentSegment with respect to the previousSegment.
    Angles_t calculateAngles(TVector3 currentSegment, TVector3 previousSegment) const;

    /// Determine if the provided generic Point is located in the detector.
    /// @see `trkf::MCSSegmentCalculator::inDetector`
    template<typename Point> inline Bool_t inDetector(Point point) const {
      return fSegmentCalculator.inDetector(point);
    }

    /// A helper method to convert vectors from type `From` to type `To`
    /// @see `trkf::MCSSegmentCalculator::convert`
    template<typename From, typename To> inline To convert(From vector) const {
      return fSegmentCalculator.convert<From, To>(vector);
    }
    
  private:
    const MCSSegmentCalculator fSegmentCalculator; /// The segment calculator used in some `GetResult` methods.  If `GetResult(const MCSSegmentResult segmentResult, const int segmentMethod)` is called, this will not be used. Otherwise, it will be.

    /// Calculate the vector of segment angles (3D angle & 2D projected angles) between the provided vector of segment vectors.
    std::vector<Angles_t> SetAngles_vec(std::vector<Vector_t> segmentVector_vec) const;
  }; // class MCSAngleCalculator
} // namespace trkf

#endif

/* TODO Potential Future:
  struct TrackSegment {
  private:
     std::vector<Point_t> fPoints;
     Point_t fBarycenter;
     Double_t fRawSegmentLength;
     Vector_t fLinearFit;
  }
*/
