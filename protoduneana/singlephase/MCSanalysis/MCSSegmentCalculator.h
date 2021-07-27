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

// C++ Includes
#include <iostream>

namespace trkf {
  class MCSSegmentCalculator {
  public:

    struct MCSSegmentResult;
    struct Point_tt;

    // Public Type Defintions
    // using Point_t = recob::tracking::Point_;
    // TODO: change all instances of Point_t to just Point_tt and use Point_t as recob::tracking::Point_t, Point_tt can just be an internal nomenclature and not overlap.
    using Point_t = Point_tt;
    using Vector_t = recob::tracking::Vector_t;
    using Barycenter_t = recob::tracking::Point_t;
    using Segment_t = std::vector<Point_t>;

    /**
       A custom Point type, that contains position information and momentum.

       This is required due to virtual points being added.  If virtual points are removed from the code, this could be removed
       and the typealias 'using Point_t = Point_tt' could be changed to 'using Point_t = recob::tracking::Point_t'
     */
    struct Point_tt {
    
      // Explicitly create point
      explicit Point_tt(Double_t x, Double_t y, Double_t z, Double_t momentum) :
	fX(x),
	fY(y),
	fZ(z),
	fMomentum(momentum)
      {}

      // Create a point with an invalid momentum
      explicit Point_tt(Double_t x, Double_t y, Double_t z) :
	Point_tt(x, y, z, -999)
      {}

      // Create a point using a tracking::Point_t, this can be done implicitly
      Point_tt(recob::tracking::Point_t point) :
	Point_tt(point.X(), point.Y(), point.Z())
      {}
      
      inline void SetP(Double_t momentum) { fMomentum = momentum; }
      inline void SetX(Double_t x) { fX = x; }
      inline void SetY(Double_t y) { fY = y; }
      inline void SetZ(Double_t z) { fZ = z; }
      
      inline Double_t P() const { return fMomentum; }
      inline Double_t X() const { return fX; }
      inline Double_t Y() const { return fY; }
      inline Double_t Z() const { return fZ; }
      inline recob::tracking::Point_t underlying() const { return recob::tracking::Point_t(fX, fY, fZ); }
      
    private:

      Double_t fX;
      Double_t fY;
      Double_t fZ;
      Double_t fMomentum;

    }; // Point_tt

    // fhicl Configuration
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<Float_t> segmentLength { Name("segmentLength"), Comment("Length of raw segments that the trajectory will be broken into. (cm)"), 14.0 };
    }; // struct Config

    // Public Constructors
    explicit MCSSegmentCalculator() :
      MCSSegmentCalculator(fhicl::Table<Config>(fhicl::ParameterSet()))
    {}
    explicit MCSSegmentCalculator(const fhicl::Table<Config> t) :
      MCSSegmentCalculator(
			   t().segmentLength()
			   )
    {}
    explicit MCSSegmentCalculator(const Float_t segmentLength) :
      fSegmentLength(segmentLength)
    {}

    // Public Methods
    MCSSegmentResult GetResult(const std::vector<Point_t> points, const Bool_t shouldCreateVirtualPoints) const;
    MCSSegmentResult GetResult(const simb::MCTrajectory trajectory, const Bool_t shouldCreateVirtualPoints) const;

    inline MCSSegmentResult GetResult(const simb::MCParticle particle, const Bool_t shouldCreateVirtualPoints) const {
      return GetResult(particle.Trajectory(), shouldCreateVirtualPoints);
    }
    inline MCSSegmentResult GetResult(const recob::Trajectory trajectory, const Bool_t shouldCreateVirtualPoints) const {
      std::vector<Point_t> points;
      for(recob::tracking::Point_t position: trajectory.Positions())
	points.push_back(Point_t(position));
      return GetResult(points, shouldCreateVirtualPoints);
    }
    inline MCSSegmentResult GetResult(const recob::TrackTrajectory trackTrajectory, const Bool_t shouldCreateVirtualPoints) const {
      return GetResult(trackTrajectory.Trajectory(), shouldCreateVirtualPoints);
    }
    inline MCSSegmentResult GetResult(const recob::Track track, const Bool_t shouldCreateVirtualPoints) const {
      return GetResult(track.Trajectory(), shouldCreateVirtualPoints);
    }
    //   template<typename From, typename To> To convert(From vector) const;
    template<typename From, typename To> inline To convert(From vector) const {
      return To(vector.X(), vector.Y(), vector.Z());
    }
    // Determine if the provided generic Point is located in the detector.
    // TODO: Change implmentation to use a Geometry Service? Or take a provided geometry service in the constructor.
    // TODO: Change to a static function, or move some of these more "helper" functions to a MCSUtils class/struct?
    // template<typename Point> Bool_t inDetector(Point point) const;
    template<typename Point> inline Bool_t inDetector(Point point) const {
      return 
	point.Z() <= 700 && point.Z() >= 0 &&
	point.X() <= 300  && point.X() >= -300 &&
	point.Y() <= 600 && point.Y() >= 0;
    }

    // Calculate a barycenter using the provided list of points in a segment.
    template<typename Point> Barycenter_t CalculateBarycenter(const std::vector<Point> segment) const;

    // Calculate a linear fit vector using the provided list of points in a segment.
    template<typename Point> Vector_t CalculateLinearFit(const std::vector<Point> segment) const;

    // Calculate the raw segment length for a segment.  The raw segment length is used in the creation of a segment, but this is an independent way to verify that the segments were created correctly.
    // TODO: Add this check to a test case?
    Float_t CalculateRawSegmentLength(const Segment_t segment) const;

    struct MCSSegmentResult {
    public:
      // Public MCSSegmentResult Constructor
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
      inline std::vector<Segment_t> GetSegment_vec() const { return fSegment_vec; }
      inline std::vector<Barycenter_t> GetBarycenter_vec() const { return fBarycenter_vec; }
      inline std::vector<Vector_t> GetLinearFit_vec() const { return fLinearFit_vec; }
      inline std::vector<Vector_t> GetPolygonalSegmentVector_vec() const { return fPolygonalSegmentVector_vec; }
      inline std::vector<Float_t> GetRawSegmentLength_vec() const { return fRawSegmentLength_vec; }
      inline std::vector<Float_t> GetPolygonalSegmentLength_vec() const { return fPolygonalSegmentLength_vec; }
      // This should only be used in the context of a true segment result
      inline Double_t GetMomentumForSegment(size_t i) const { return fSegment_vec.at(i).at(0).P(); }

    private:
      const std::vector<Segment_t> fSegment_vec;
      const std::vector<Barycenter_t> fBarycenter_vec;
      
      // List of Segment Vectors
      const std::vector<Vector_t> fLinearFit_vec;
      const std::vector<Vector_t> fPolygonalSegmentVector_vec;
      // const std::vector<Vector_t> fLastPointVector_vec;
      // const std::vector<Vector_t> fCrossover_vec;

      // List of Segment Lengths
      const std::vector<Float_t> fRawSegmentLength_vec;
      const std::vector<Float_t> fPolygonalSegmentLength_vec;
    }; // struct MCSSegmentResult

  private:
    Float_t fSegmentLength;
    // TODO: Segmentation Mode set in constructor? For example: Add lengths between point or just use diff between first and last point?
    // TODO: Indicate whether or not to use Virtual Points?

    // Private Set Functions

    // Calculate the vector of segments, optionally creating virtual points.
    std::vector<Segment_t> CreateSegment_vec(const std::vector<Point_t> points, Bool_t shouldCreateVirtualPoints) const;

    // Calculate the vector of raw segment lengths using the provided vector of segments.
    std::vector<Float_t> SetRawSegmentLength_vec(const std::vector<Segment_t> segment_vec) const;

    // Calculate the vector segment start points for the provided vector of trajectory points, also setting the provided vector of segment lengths.
    // std::vector<size_t> SetSegmentStartPoint_vec(const std::vector<Point_t> trajPoint_vec, std::vector<Float_t>& rawSegmentLength_vec) const;

    // Calculate the vector of barycenters for the provided vector of trajectory points, using the provided vector of segment start points as segment breakpoints.
    // TODO: Consider SetBarycenter_vec(std::vector<std::vector<Point_t> > where it's a vector of a segment where the segment is a vector of points
    // std::vector<Barycenter_t> SetBarycenter_vec(const std::vector<Point_t> trajPoint_vec, const std::vector<size_t> segmentStartPoint_vec) const;
    std::vector<Barycenter_t> SetBarycenter_vec(const std::vector<Segment_t> segment_vec) const;
    
    // Calculate the vector of polygonal segment lengths for the provided vector of barycenters, where the polygonal segment length is the distance between barycenters.  The first entry of the returned vector is the length between the provided startPoint (which should be the first point in the vector) and the first barycenter.
    std::vector<Float_t> SetPolygonalSegmentLength_vec(const std::vector<Barycenter_t> barycenter_vec, const Point_t startPoint) const;

    // Calculate the vector linear fit vectors for the provided vector of trajectory points, using the provided vector of segment start points as segment breakpoints.
    // std::vector<Vector_t> SetLinearFit_vec(const std::vector<Point_t> trajPoint_vec, const std::vector<size_t> segmentStartPoint_vec) const;
    std::vector<Vector_t> SetLinearSegmentVector_vec(const std::vector<Segment_t> segment_vec) const;

    // TODO: Documenatation
    std::vector<Vector_t> SetPolygonalSegmentVector_vec(const std::vector<Barycenter_t> barycenter_vec, const Point_t startPoint) const;

    // TODO: Implement these functions later:
    /*
    //    std::vector<Vector_t> SetLastPointVector_vec( idk what goes here

    // This will take a linear fit just the points on either side of a segment start point.
    std::vector<Vector_t> SetCrossover_vec(const std::vector<Point_t> trajPoint_vec, const std::vector<size_t> segmentStartPoint_vec) const;
    */

  }; // class MCSSegmentCalculator

} // namespace trkf

#endif
