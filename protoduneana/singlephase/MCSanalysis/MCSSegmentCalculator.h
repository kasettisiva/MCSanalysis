// -*- mode: c++ -*-
#ifndef MCSSEGMENTCALCULATOR_H
#define MCSSEGMENTCALCULATOR_H

// Framework Includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

// LArSoft Includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Trajectory.h"
#include "lardataobj/RecoBase/TrackingTypes.h"

// LArSoft Includes for SCE offset
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larcore/Geometry/Geometry.h"

// C++ Includes
#include <iostream>

namespace trkf {
  /**
     @author Hunter C. Meyer (LSU, hmeyer5@lsu.edu)
     @file protoduneana/singlephase/MCSanalysis/MCSSegmentCalculator.h
     @class trkf::MCSSegmentCalculator MCSSegmentCalculator.h "protoduneana/singlephase/MCSanalysis/MCSSegmentCalculator.h"
     @brief A calculator that will generate segments of a true or reconstructed trajectory.
     
     SCE corrections will be made for reconstructed tracks if `services.SpaceCharge.EnableCalSpatialSCE` is set to true in the fhicl configuration.
     
     For the fhicl configuration, see `trkf::MCSSegmentCalculator::Config`.
     
     @todo Move to `larreco/RecoAlg/MCSSegmentCalculator.h` and update documentation.
     @todo For helper functions, consider making static functions, or moving to an MCSUtilities class/struct?  Or possibly in `protoduneana/Utilities/`?
     @todo Check that documentation has been generated correctly. Add to documentation groups
     @todo Consider making SetX_vec functions into template functions.
     @todo Document the various segment defintions. Add more segment defintions.
     @todo Consider if virtual points should be declared in the constructor (via a private property fShouldCreateVirtualPoints and used in the `Config`) or in the `GetResult` methods.
     
     @see trkf::MCSAngleCalculator
     @see trkf::MCSMomentumCalculator
     @see trkf::MCSSegmentCalculator::Config
   */
  class MCSSegmentCalculator {
  public:

    struct MCSSegmentResult;
    struct Point_tt;

    // Public Type Defintions
    using Point_t = Point_tt;
    using Vector_t = recob::tracking::Vector_t;
    using Barycenter_t = recob::tracking::Point_t;
    using Segment_t = std::vector<Point_t>;

    /**
       @brief A custom Point type, that contains position information and momentum.
       
       @details This is required due to virtual points being added.  If virtual points are removed from the code, this could be removed
       and the typealias \code using Point_t = Point_tt \endcode could be changed to \code using Point_t = recob::tracking::Point_t \endcode
       
       @see trkf::MCSSegmentCalculator::Point_t
     */
    struct Point_tt {
    
      /// Explicitly create point.
      explicit Point_tt(Double_t x, Double_t y, Double_t z, Double_t momentum) :
	fX(x),
	fY(y),
	fZ(z),
	fMomentum(momentum)
      {}

      /// Create a point with only an X, Y, and Z coordinate.
      /// An invalid momentum of -999 will be set.
      explicit Point_tt(Double_t x, Double_t y, Double_t z) :
	Point_tt(x, y, z, -999)
      {}

      /// Create a point using a tracking::Point_t, this can be done implicitly.
      /// An invalid momentum of -999 will be set.
      Point_tt(recob::tracking::Point_t point) :
	Point_tt(point.X(), point.Y(), point.Z())
      {}
      
      /// Update this Point's Momentum in GeV/c.
      inline void SetP(Double_t momentum) { fMomentum = momentum; }
      /// Update this Point's X coordinate in cm.
      inline void SetX(Double_t x) { fX = x; }
      /// Update this Point's Y coordinate in cm.
      inline void SetY(Double_t y) { fY = y; }
      /// Update this Point's Z coordinate in cm.
      inline void SetZ(Double_t z) { fZ = z; }
      
      /// Get this Point's Momentum in GeV/c, which is only valid in a "true" context.  In the case of reconstruction, -999 is returned.
      inline Double_t P() const { return fMomentum; }
      /// Get this Point's X coordinate in cm.
      inline Double_t X() const { return fX; }
      /// Get this Point's Y coordinate in cm.
      inline Double_t Y() const { return fY; }
      /// Get this Point's Z coordinate in cm.
      inline Double_t Z() const { return fZ; }
      /// Access this Point's X, Y, and Z coordinates as a `recob::tracking::Point_t`.
      inline recob::tracking::Point_t underlying() const { return recob::tracking::Point_t(fX, fY, fZ); }
      
    private:

      Double_t fX; /// This Point's X coordinate (cm).
      Double_t fY; /// This Point's Y coordinate (cm).
      Double_t fZ; /// This Point's Z coordinate (cm).
      Double_t fMomentum; /// This Point's Momentum (GeV/c).

    }; // Point_tt

    /// fhicl Configuration used for initialization of `trkf::MCSSegmentCalculator`
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<Float_t> segmentLength { Name("segmentLength"), Comment("Length of raw segments that the trajectory will be broken into. (cm)"), 14.0 };
    }; // struct Config
    
    // Public Constructors
    /// Default Constructor, using an empty `fhicl::ParameterSet` for the fhicl configuration, leading to default values being used in `trkf::MCSSegmentCalculator::Config`
    explicit MCSSegmentCalculator() :
      MCSSegmentCalculator(fhicl::Table<Config>(fhicl::ParameterSet()))
    {}
    /// Constructor where you can pass a `fhicl::Table<Config>`
    explicit MCSSegmentCalculator(const fhicl::Table<Config> t) :
      MCSSegmentCalculator(
			   t().segmentLength()
			   )
    {}
    /// Constructor where you can pass each parameter in `trkf::MCSSegmentCalculator::Config` explicitly.
    /// @param segmentLength The length of each raw-segment. See `trkf::MCSSegmentCalculator::CreateSegments`.
    explicit MCSSegmentCalculator(const Float_t segmentLength) :
      fSegmentLength(segmentLength),
      fSCE(lar::providerFrom<spacecharge::SpaceChargeService>()),
      fGeom(art::ServiceHandle<geo::Geometry const>())
    {}

    // Public Methods

    /// Retrieve an `MCSSegmentResult` using a vector of points.  This is the base method that all other `GetResult` methods call.
    /// @param points The points to use to create the `MCSSegmentResult`.
    /// @param shouldCreateVirtualPoints A boolean tag that specifies whether or not virtual points should be created with the segments, forcing the segments to be 14-cm.
    MCSSegmentResult GetResult(const std::vector<Point_t> points, const Bool_t shouldCreateVirtualPoints) const;

    /// Retrieve an `MCSSegmentResult` using a `simb::MCTrajectory`.
    /// 
    /// @details A vector of Points will be created from the MCTrajectory and used to create a `MCSSegmentResult`
    /// 
    /// @param trajectory The trajectory to use to create the `MCSSegmentResult`.
    /// @param shouldCreateVirtualPoints A boolean tag that specifies whether or not virtual points should be created with the segments, forcing the segments to be 14-cm.
    /// @see GetResult(const std::vector<Point_t> points, const Bool_t shouldCreateVirtualPoints) const
    MCSSegmentResult GetResult(const simb::MCTrajectory trajectory, const Bool_t shouldCreateVirtualPoints) const;

    /// Retrieve an `MCSSegmentResult` using a `simb::MCParticle`.
    /// @param MCParticle The particle used to create the `MCSSegmentResult`.
    /// @param shouldCreateVirtualPoints A boolean tag that specifies whether or not virtual points should be created with the segments, forcing the segments to be 14-cm.
    /// @see trkf::MCSSegmentCalculator::GetResult(const simb::MCTrajectory trajectory, const Bool_t shouldCreateVirtualPoints) const
    inline MCSSegmentResult GetResult(const simb::MCParticle particle, const Bool_t shouldCreateVirtualPoints) const {
      return GetResult(particle.Trajectory(), shouldCreateVirtualPoints);
    }

    /// Retrieve an `MCSSegmentResult` using a `recob::Trajectory`.
    /// 
    /// @details A vector of Points will be created by looping through `trajectory.Positions()`, each with an invalid momentum.
    /// SCE corrections are made if `services.SpaceCharge.EnableCalSpatialSCE` is set to true in the fhicl configuration.
    /// 
    /// @todo Delete commented out code.  Currently kept for if it is needed.
    /// 
    /// @param trajectory The trajectory used to create the `MCSSegmentResult`.
    /// @param shouldCreateVirtualPoints A boolean tag that specifies whether or not virtual points should be created with the segments, forcing the segments to be 14-cm.
    /// @see GetResult(const std::vector<Point_t> points, const Bool_t shouldCreateVirtualPoints)
    inline MCSSegmentResult GetResult(const recob::Trajectory trajectory, const Bool_t shouldCreateVirtualPoints) const {

      std::vector<Point_t> points;

      for(recob::tracking::Point_t position: trajectory.Positions()) {

	// This should not be required, as if this was false, then an offset of 0 is used.  I've included this conditional it for clarity.
	if(fSCE->EnableCalSpatialSCE()) {

	  // private members fGeom & fSCE set in constructor.
	  // spacecharge::SpaceCharge const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
	  // art::ServiceHandle<geo::Geometry const> geom;

	  geo::Point_t point = convert<recob::tracking::Point_t, geo::Point_t>(position);
	
	  auto tpc = fGeom->FindTPCAtPosition(point).TPC;
	  geo::Vector_t offset = fSCE->GetCalPosOffsets(point, tpc);

	  Double_t xOffset = -offset.X(); // not sure why this is negative.
	  Double_t yOffset = offset.Y();
	  Double_t zOffset = offset.Z();

	  /*
	    std::cout << std::endl;
	    std::cout << "Point = (" << point.X() << ", " << point.Y() << ", " << point.Z() << ")" << std::endl;
	    std::cout << "TPC = " << tpc << std::endl;
	    std::cout << "SCE Offset = (" << xOffset << ", " << yOffset << ", " << zOffset << ")" << std::endl;
	  // */

	  position.SetX(position.X() + xOffset);
	  position.SetY(position.Y() + yOffset);
	  position.SetZ(position.Z() + zOffset);

	} // if SCE corrections should be applied

	points.push_back(Point_t(position));
      }

      return GetResult(points, shouldCreateVirtualPoints);
    }

    /**
       Retrieve an `MCSSegmentResult` using a `recob::TrackTrajectory`.
       @param trackTrajectory The track-trajectory to use to create the `MCSSegmentResult`
       @param shouldCreateVirtualPoints A boolean tag that specifies whether or not virtual points should be created with the segments, forcing the segments to be 14-cm.
       @see GetResult(const recob::Trajectory trajectory, const Bool_t shouldCreateVirtualPoints)
       @todo Consider using `trackTrajectory.NextValidPoint` methodology, which would require removing `GetResult(const recob::Trajectory trajectory, const Bool_t shouldCreateVirtualPoints)`. 
     */
    inline MCSSegmentResult GetResult(const recob::TrackTrajectory trackTrajectory, const Bool_t shouldCreateVirtualPoints) const {
      return GetResult(trackTrajectory.Trajectory(), shouldCreateVirtualPoints);
    }

    /**
       Retrieve an `MCSSegmentResult` using a `recob::Track`.
       @param track The track to use to create the `MCSSegmentResult`
       @param shouldCreateVirtualPoints A boolean tag that specifies whether or not virtual points should be created with the segments, forcing the segments to be 14-cm.
       @see GetResult(const recob::TrackTrajectory trackTrajectory, const Bool_t shouldCreateVirtualPoints) const 
     */
    inline MCSSegmentResult GetResult(const recob::Track track, const Bool_t shouldCreateVirtualPoints) const {
      return GetResult(track.Trajectory(), shouldCreateVirtualPoints);
    }

    /**
       A helper method to convert vectors from type `From` to type `To`
       
       Example usage:
       \code
       TVector3 vector(1, 2, 3);
       recob::tracking::Point_t point = convert<TVector3, recob::tracking::Point_t>(vector);
       \endcode
    */
    template<typename From, typename To> inline To convert(From vector) const {
      return To(vector.X(), vector.Y(), vector.Z());
    }

    /**
       Determine if the provided generic Point is located in the detector.
       @todo Change implmentation to use a Geometry Service? Or take a provided geometry service in the constructor.
    */
    template<typename Point> inline Bool_t inDetector(Point point) const {
      return 
	point.Z() <= 700 && point.Z() >= 0 &&
	point.X() <= 300  && point.X() >= -300 &&
	point.Y() <= 600 && point.Y() >= 0;
    }

    /// Calculate a barycenter (geometric average) using the provided list of points in a segment.
    template<typename Point> Barycenter_t CalculateBarycenter(const std::vector<Point> segment) const;

    /// Calculate a linear fit vector using the provided list of points in a segment.
    template<typename Point> Vector_t CalculateLinearFit(const std::vector<Point> segment) const;

    /// Calculate the raw segment length for a segment.  The raw segment length is used in the creation of a segment, but this is an independent way to verify that the segments were created correctly.
    /// @todo In `CreateSegment_vec`, check that the length is correct using this method before adding the vector to the list of segments.
    /// @todo Consider making this a template function \code template<typename Point> Float_t CalculateRawSegmentLength(const std::vector<Point> segment) const; \endcode
    Float_t CalculateRawSegmentLength(const Segment_t segment) const;

    /**
       @struct trkf::MCSSegmentCalculator::MCSSegmentResult MCSSegmentCalculator.h "protoduneana/singlephase/MCSanalysis/MCSSegmentCalculator.h"
       
       @brief A container for segment information, including the vectors that are formed for various segment defintions and the segment lengths for those segment defintions.
       
       @todo Should the details of the various segment defintions be outlined here?
       @todo Document the relative lengths of vectors.
       @todo Consider moving to lardataobj instead of larreco, as was done in the `TrajectoryMCSFitter`
       @todo Specify whether or not SCE corrections were applied.
     */
    struct MCSSegmentResult {
    public:
      /// Explicitly create an `MCSSegmentResult`.
      explicit MCSSegmentResult(
				// Basic Segment Info
				std::vector<Segment_t> segment_vec,
				std::vector<Barycenter_t> barycenter_vec,

				// Various Segment Vectors
				std::vector<Vector_t> linearFit_vec,
				std::vector<Vector_t> polygonalSegmentVector_vec,

				// Various Segment Length Defintions
				std::vector<Float_t> rawSegmentLength_vec,
				std::vector<Float_t> polygonalSegmentLength_vec
				) : 
	fSegment_vec(segment_vec),
	fBarycenter_vec(barycenter_vec),
	fLinearFit_vec(linearFit_vec),
	fPolygonalSegmentVector_vec(polygonalSegmentVector_vec),
	fRawSegmentLength_vec(rawSegmentLength_vec),
	fPolygonalSegmentLength_vec(polygonalSegmentLength_vec)
      {} // MCSSegmentResult constructor

      // Public Get Functions
      /// Get the result's vector of segments.
      inline std::vector<Segment_t> GetSegment_vec() const { return fSegment_vec; }
      /// Get the result's vector of barycenters.
      inline std::vector<Barycenter_t> GetBarycenter_vec() const { return fBarycenter_vec; }
      /// Get the result's vector of linear-fit vectors.
      inline std::vector<Vector_t> GetLinearFit_vec() const { return fLinearFit_vec; }
      /// Get the result's vector of polygonal-segment vectors.
      inline std::vector<Vector_t> GetPolygonalSegmentVector_vec() const { return fPolygonalSegmentVector_vec; }
      /// Get the result's vector of raw -segment lengths.
      inline std::vector<Float_t> GetRawSegmentLength_vec() const { return fRawSegmentLength_vec; }
      /// Get the result's polygonal-segment lengths.
      inline std::vector<Float_t> GetPolygonalSegmentLength_vec() const { return fPolygonalSegmentLength_vec; }
      /// Get the start momentum of segment i.  This should only be used in the context of a true segment result.
      inline Double_t GetMomentumForSegment(size_t i) const { return fSegment_vec.at(i).at(0).P(); }

    private:
      const std::vector<Segment_t> fSegment_vec; /// This result's vector of segments.
      const std::vector<Barycenter_t> fBarycenter_vec; /// This result's vector of barycenters.
      
      // List of Segment Vectors
      const std::vector<Vector_t> fLinearFit_vec; /// This result's vector of linear-fit vectors.
      const std::vector<Vector_t> fPolygonalSegmentVector_vec; /// This result's vector of polygonal-segment vectors.
      // const std::vector<Vector_t> fLastPointVector_vec;
      // const std::vector<Vector_t> fCrossover_vec;

      // List of Segment Lengths
      const std::vector<Float_t> fRawSegmentLength_vec; /// This result's vector of raw-segment lengths.
      const std::vector<Float_t> fPolygonalSegmentLength_vec; /// This result's vector of polygonal-segment lengths.
    }; // struct MCSSegmentResult

  private:
    // Set in constructor.
    const Float_t fSegmentLength; /// The length to use for segments, set during construction.
    const spacecharge::SpaceCharge* fSCE;
    const art::ServiceHandle<geo::Geometry const> fGeom;

    // Private Set Functions

    /// Calculate the vector of segments, optionally creating virtual points.
    std::vector<Segment_t> CreateSegment_vec(const std::vector<Point_t> points, Bool_t shouldCreateVirtualPoints) const;

    /// Calculate the vector of raw segment lengths using the provided vector of segments.
    std::vector<Float_t> SetRawSegmentLength_vec(const std::vector<Segment_t> segment_vec) const;

    /// Calculate the vector of barycenters for the provided vector of segments.
    std::vector<Barycenter_t> SetBarycenter_vec(const std::vector<Segment_t> segment_vec) const;
    
    /// Calculate the vector of polygonal segment lengths for the provided vector of barycenters, where the polygonal segment length is the distance between barycenters.  
    /// The first entry of the returned vector is the length between the provided startPoint (which should be the first point of the trajectory) and the first barycenter.
    std::vector<Float_t> SetPolygonalSegmentLength_vec(const std::vector<Barycenter_t> barycenter_vec, const Point_t startPoint) const;

    /// Calculate the vector linear fit vectors for the provided vector of segments.
    std::vector<Vector_t> SetLinearSegmentVector_vec(const std::vector<Segment_t> segment_vec) const;

    /// Calculate the vector of polygonal-segment vectors for the provided vector of points.
    /// The first entry of the returned vector is the vector connecting the startPoint to the first barycenter.
    std::vector<Vector_t> SetPolygonalSegmentVector_vec(const std::vector<Barycenter_t> barycenter_vec, const Point_t startPoint) const;

  }; // class MCSSegmentCalculator

} // namespace trkf

#endif
