#include "MCSSegmentCalculator.h"

using namespace trkf;

using Point_t = MCSSegmentCalculator::Point_t;
using Vector_t = MCSSegmentCalculator::Vector_t;
using Barycenter_t = MCSSegmentCalculator::Barycenter_t;
using MCSSegmentResult = MCSSegmentCalculator::MCSSegmentResult;
using Segment_t = MCSSegmentCalculator::Segment_t;

/// @todo Add assertions that the relative sizes of vectors that are created are correct.
MCSSegmentResult MCSSegmentCalculator::GetResult(const std::vector<Point_t> points, const Bool_t shouldCreateVirtualPoints) const {
  // Create the vector of segments
  std::vector<Segment_t> segment_vec = CreateSegment_vec(points, shouldCreateVirtualPoints);
  // Calculate the vector of the "raw length" of each segment.
  std::vector<Float_t> rawSegmentLength_vec = SetRawSegmentLength_vec(segment_vec);
  // Calculate the vector of barycenters
  std::vector<Barycenter_t> barycenter_vec = SetBarycenter_vec(segment_vec);
  // Get the first trajectory point of the trajectory, used by polygonal segments as the first-point.
  Point_t startPoint = segment_vec.at(0).at(0);
  // Calculate the vector of the "polygonal segment lengths" (distance between barycenters)
  std::vector<Float_t> polygonalSegmentLength_vec = SetPolygonalSegmentLength_vec(barycenter_vec, startPoint);
  // Create the vector of linear segment vectors. (Linear fit of each segment)
  std::vector<Vector_t> linearSegmentVector_vec = SetLinearSegmentVector_vec(segment_vec);
  // Calculate the vector of polygonal segment vectors. (Vector connecting barycenters of each segment)
  std::vector<Vector_t> polygonalSegmentVector_vec = SetPolygonalSegmentVector_vec(barycenter_vec, startPoint);

  // Create the MCSSegmentResult
  return MCSSegmentResult(
			  segment_vec,
			  barycenter_vec,
			  linearSegmentVector_vec,
			  polygonalSegmentVector_vec,
			  rawSegmentLength_vec,
			  polygonalSegmentLength_vec
			  );
} // MCSSegmentCalculator::GetResult(const std::vector<Point_t>, const Bool_t)

MCSSegmentResult MCSSegmentCalculator::GetResult(const simb::MCTrajectory trajectory, const Bool_t shouldCreateVirtualPoints) const {
  // The vector of trajectory points that's about to be filled.
  std::vector<Point_t> trajPoint_vec;

  // Loop through the trajectory points
  for(size_t i = 0; i < trajectory.size(); ++i) {
    // Get the position of the trajectory point.
    TVector3 position = trajectory.Position(i).Vect();
    // Conver the point to a Point_t
    Point_t point = convert<TVector3, Point_t>(position);
    // Set the momentum of that point
    point.SetP(trajectory.Momentum(i).P());
    // If the point is in the detector, add it to the vector of trajectory points
    if(inDetector(point))
      trajPoint_vec.push_back(point);
  }

  // Get the result for the vector of trajectory points, and return it.
  MCSSegmentResult result = GetResult(trajPoint_vec, shouldCreateVirtualPoints);
  return result;
} // MCSSegmentCalculator::GetResult(const simb::MCTrajectory, const Bool_t)

/// @todo Add assertions that the previous length is never greater than fSegmentLength.  Organize algorithm such that if that happens, it's handled properly.
/// @todo Are there other assertions that can be made DURING the algorithm, not just afterwards?  
/// Assertion and/or Test Possibilities: 
/// 1. Assert that every segment has more than 1 point.  
/// 2. Another one: The start of one segment is *always* the end of the previous segment.
/// 3. In the middle of the algoirhtm, if the number of points in a segment is 1, the length should be 0.
/// 4. No segments have duplicate points.
/// 
/// @todo Can this function be made into smaller, more testable functions?
/// @todo A maximum of 2 virtual points will ever be made per segment.  This should be necessarily the case if the gaps between points is larger than 3 segment lengths.
/// @todo Consider rewrite of segmentation algorithm.  Assumptions might be made that no longer apply.  Made be able to reorder the if statements to prevent bugs.
/// @todo This should be an easily unit-tested function.
/// @todo For virtual points, should the momentum be the closest point's momentum, the average momentum of the points on either side, or the previous point's momentum?  Currently, it's the previous point's momentum.
/// @todo Revisit diff >= 2*fSegmentLength.  Should this be diff >= 2*fSegmentLength or currLength >= 2*fSegmentLength, or perhaps an alternative? Leaning towards currLength but it's currently diff.
/// @todo Revisit diff >= 2*fSegmentLength where shouldCreateVirtualPoints is false.
/// @todo Consider what should happen with the last, non 14-cm segment.
/// @todo: Segmentation Mode set in constructor? For example: Add lengths between point or just use diff between first and last point?
std::vector<Segment_t> MCSSegmentCalculator::CreateSegment_vec(const std::vector<Point_t> points, Bool_t shouldCreateVirtualPoints) const {
  std::vector<Segment_t> segments;
  Float_t currLength = 0; // cm
  Float_t prevLength = 0; // cm

  Segment_t currSegment;
  
  // Loop through the trajectory points
  for(size_t i = 1; i < points.size(); ++i) {
    // Define the current and previous trajectory point.
    Point_t prevPoint = points.at(i-1);
    Point_t currPoint = points.at(i);

    // If this point is invalid.
    if(currPoint.X() == -999)
      break;

    // Define the difference between the current and previous trajectory point.
    Float_t diff = sqrt((currPoint.underlying() - prevPoint.underlying()).Mag2());
    // Update the length.
    currLength += diff;

    if(diff >= 2*fSegmentLength) {
      // If the current length is greater than 2*fSegmentLength, then a special case has occured (kinda rare).
      if(shouldCreateVirtualPoints) {
	// If shouldCreateVirtualPoints, multiple virtual points should be created.

	// Create the first virtual point
	Float_t multiplier1 = (fSegmentLength - prevLength) / diff;
	Float_t x1 = multiplier1*(currPoint.X() - prevPoint.X()) + prevPoint.X();
	Float_t y1 = multiplier1*(currPoint.Y() - prevPoint.Y()) + prevPoint.Y();
	Float_t z1 = multiplier1*(currPoint.Z() - prevPoint.Z()) + prevPoint.Z();
	Point_t virtualPoint1(x1, y1, z1);
	// Set the first virtual point's momentum.
	virtualPoint1.SetP(prevPoint.P());

	// Add the first virtual point to the current segment (to end the segment at 14-cm).
	currSegment.push_back(virtualPoint1);

	// Add the current segment to the segments vector and clear the current segment.
	segments.push_back(currSegment);
	currSegment.clear();

	// Create the second virtual point
	Float_t multiplier2 = (2*fSegmentLength - prevLength) / diff;
	Float_t x2 = multiplier2*(currPoint.X() - prevPoint.X()) + prevPoint.X();
	Float_t y2 = multiplier2*(currPoint.Y() - prevPoint.Y()) + prevPoint.Y();
	Float_t z2 = multiplier2*(currPoint.Z() - prevPoint.Z()) + prevPoint.Z();
	Point_t virtualPoint2(x2, y2, z2);
	virtualPoint2.SetP(prevPoint.P());

	// Add both virtual points to the current segment (legngth between them should be 14-cm).
	Segment_t virtualSegment;
	virtualSegment.push_back(virtualPoint1);
	virtualSegment.push_back(virtualPoint2);
	
	// Add the virtual segment to the segments vector (no need to clear the virtual segment, it will be lost in scope).
	segments.push_back(virtualSegment);

	// Add the first couple of points to the next segment, update the length.
	currSegment.push_back(virtualPoint2);
	currSegment.push_back(currPoint);
	currLength -= 2*fSegmentLength;
      } else {
	// If you should not create virtual points, the closest point will be chosen.
	// An exception will be made if the current segment only has one point, to prevent a segment of zero length from occuring.
	currSegment.push_back(currPoint);
	segments.push_back(currSegment);
	currSegment.clear();
	currSegment.push_back(currPoint);
	currLength = 0;
	/* Temporary removal for testing purposes.
	// With this commented out, it's the same as the previous segmentation methods, though that was actually wrong.
	if(currSegment.size() == 1) {
	  // Add the current point to the current segment.
	  currSegment.push_back(currPoint);

	  // Add the current segment to the vector of segments.  Clear the current segment.
	  segments.push_back(currSegment);
	  currSegment.clear();

	  // Add the current point to the next segment.
	  currSegment.push_back(currPoint);
	  currLength = 0;
	} else {
	  // End the current segment.
	  // Add the current segment to the vector of segments.  Clear the current segment.
	  segments.push_back(currSegment);
	  currSegment.clear();

	  // Create a segment using just these previous point and current point. (Remember we're in the context where diff > 2* segmentLength)
	  currSegment.push_back(prevPoint);
	  currSegment.push_back(currPoint);

	  // Add this two-point segment to the vector of segments.  Clear the current segment.
	  segments.push_back(currSegment);
	  currSegment.clear();

	  // Add the current point as the first point of the next segment, set the length to 0.
	  currSegment.push_back(currPoint);
	  currLength = 0;
	} // else
	*/
      } // else
    } else if((currLength > fSegmentLength && prevLength <= fSegmentLength) || prevLength > fSegmentLength) {
      if(shouldCreateVirtualPoints) {
	// Create a virutal point
	Float_t multiplier = (fSegmentLength - prevLength) / diff;
	Float_t x = multiplier*(currPoint.X() - prevPoint.X()) + prevPoint.X();
	Float_t y = multiplier*(currPoint.Y() - prevPoint.Y()) + prevPoint.Y();
	Float_t z = multiplier*(currPoint.Z() - prevPoint.Z()) + prevPoint.Z();
	Point_t virtualPoint(x, y, z);
	virtualPoint.SetP(prevPoint.P());

	// Add the virtual point as the last point of the current segment.
	currSegment.push_back(virtualPoint);

	// Add the current segment to the vector of segments.  Clear the current segment.
	segments.push_back(currSegment);
	currSegment.clear();

	// Add the first two points of the next segment.  Update the length.
	currSegment.push_back(virtualPoint);
	currSegment.push_back(currPoint);
	currLength -= fSegmentLength;
      } else {
	// If the current point would form a segment length closer to fSegmentLength than the previous point.
	if(fSegmentLength - prevLength >= currLength - fSegmentLength) {
	  // Add the current point to the current segment.
	  currSegment.push_back(currPoint);

	  // End the current segment with the last point being the previous point.
	  // Add the current segment to the vector of segments.  Clear the current segment.
	  segments.push_back(currSegment);
	  currSegment.clear();

	  // Add the first point to the next segment.  Update the length.
	  currSegment.push_back(currPoint);
	  currLength = 0;
	} else {
	  // End the current segment.
	  // Add the current segment to the vector of segments.  Clear the current segment.
	  segments.push_back(currSegment);
	  currSegment.clear();

	  // Add the first two points to the next segment.  Update the length.
	  currSegment.push_back(prevPoint);
	  currSegment.push_back(currPoint);
	  currLength = diff;
	} // else
      } // else
    } else {
      currSegment.push_back(currPoint);
    }

    // Update the oldLength.
    prevLength = currLength;
  } // for i

  // Add the last segment, which never reached 14-cm.
  segments.push_back(currSegment);

  return segments;
} // MCSSegmentCalculator::CreateSegment_vec

Float_t MCSSegmentCalculator::CalculateRawSegmentLength(const Segment_t segment) const {
  Float_t rawSegmentLength = 0;

  // Loop through the points in the segment.
  for(size_t i = 1; i < segment.size(); ++i) {
    // Get the current and previous points.
    Point_t prevPoint = segment.at(i-1);
    Point_t currPoint = segment.at(i);

    // Calculate the distance between the points and add it to the rawSegmentLength
    Float_t diff = sqrt((currPoint.underlying() - prevPoint.underlying()).Mag2());
    rawSegmentLength += diff;
  } // for i

  return rawSegmentLength;
} // MCSSegmentCalculator::CalculateRawSegmentLength

/// @details If n is the size of segment_vec, the returned vector should be of size n-1.
/// This is because the last segment was not fully formed.
std::vector<Float_t> MCSSegmentCalculator::SetRawSegmentLength_vec(const std::vector<Segment_t> segment_vec) const {
  // The vector of raw-segment lengths
  std::vector<Float_t> rawSegmentLength_vec;
  
  // Loop through the vector of segments, with the exception of the last segment, since it was not fully formed.
  for(size_t i = 0; i < segment_vec.size() - 1; ++i) {
    // Get this segment.
    Segment_t segment = segment_vec.at(i);
    // Calculate the raw segment length for this segment.
    Float_t rawSegmentLength = CalculateRawSegmentLength(segment);
    // Append this raw segment length to the vector of raw-segment lengths.
    rawSegmentLength_vec.push_back(rawSegmentLength);
  }
  
  // Return the vector of raw-segment lengths.
  return rawSegmentLength_vec;
} // MCSSegmentCalculator:SetRawSegmentLength_vec

/// @details If n is the size of segment_vec, the returned vector should be of size n-1.
/// This is because the last segment was not fully formed.
std::vector<Barycenter_t> MCSSegmentCalculator::SetBarycenter_vec(const std::vector<Segment_t> segment_vec) const {
  // The vector of barycenters.
  std::vector<Barycenter_t> barycenter_vec;

  // Loop through the vector of segments, with the exception of the last segment, since it was not fully formed.
  for(size_t i = 0; i < segment_vec.size() - 1; ++i) {
    // Get this segment.
    Segment_t segment = segment_vec.at(i);
    // Calculate the barycenter for this segment.
    Barycenter_t barycenter = CalculateBarycenter(segment);
    // Append this barycenter to the vector of barycenters.
    barycenter_vec.push_back(barycenter);
  }

  // Return the vector of barycenters
  return barycenter_vec;
} // MCSSegmentCalculator::SetBarycenter_vec

/// @details If n-1 is the size of barycenter_vec, the returned vector should be of size n-1.
/// @todo Add the distance between two points function to `MCSSegmentCalculator` or some Utilities class, and use it here.
std::vector<Float_t> MCSSegmentCalculator::SetPolygonalSegmentLength_vec(const std::vector<Barycenter_t> barycenter_vec, Point_t startPoint) const {
  // The vector of polygonal-segment lengths.
  std::vector<Float_t> polygonalSegmentLength_vec;
  
  // Calcluate the distance between the firstBarycenter and the startPoint.
  Barycenter_t previousBarycenter = barycenter_vec.at(0);
  Float_t x = previousBarycenter.X() - startPoint.X();
  Float_t y = previousBarycenter.Y() - startPoint.Y();
  Float_t z = previousBarycenter.Z() - startPoint.Z();
  Float_t diff = sqrt(x*x + y*y + z*z);

  // Add this first distance to the vector of polygonal-segment lengths.
  polygonalSegmentLength_vec.push_back(diff);

  // Loop through the barycenters, starting with i = 1 (the second barycenter).
  for(size_t i = 1; i < barycenter_vec.size(); i++) {
    // The barycenter for this segment.
    Barycenter_t currentBarycenter = barycenter_vec.at(i);

    // Calculate the distance between the current barycenter and the previous barycenter.
    x = currentBarycenter.X() - previousBarycenter.X();
    y = currentBarycenter.Y() - previousBarycenter.Y();
    z = currentBarycenter.Z() - previousBarycenter.Z();
    Float_t diff = sqrt(x*x + y*y + z*z);

    // Push back length
    polygonalSegmentLength_vec.push_back(diff);

    // Update the previousBarycenter
    previousBarycenter = currentBarycenter;
  }

  // Return the vector of polygonal-segment lengths
  return polygonalSegmentLength_vec;
} // MCSSegmentCalculator::SetPolygonalSegmentLength_vec

/// @details If n is the size of segment_vec, the returned vector should be of size n-1.
/// This is the because the last segment is not fully formed.
std::vector<Vector_t> MCSSegmentCalculator::SetLinearSegmentVector_vec(const std::vector<Segment_t> segment_vec) const {
  // The vector of linear-fit vectors.
  std::vector<Vector_t> linearFit_vec;

  // Loop through the vector of segments.
  for(size_t i = 0; i < segment_vec.size() - 1; ++i) {
    // The current segment.
    Segment_t segment = segment_vec.at(i);
    // The linear-fit vector for the current segment.
    Vector_t linearFit = CalculateLinearFit(segment);
    // Append the linear-fit vector to the vector of linear-fit vectors.
    linearFit_vec.push_back(linearFit);
  }

  // Return the vector of linear-fit vectors
  return linearFit_vec;
} // MCSSegmentCalculator::SetLinearSegmentVector_vec

/// @details If n-1 is the size of barycenter_vec, the returned vector should be of size n-1.
/// @todo Change all of these vectors to be of unit length, like the linear-fit vectors.
/// Alternatively, make the vectors be the length that they should be and add functions to segmentResult that will map the legnths for backwards compatibility.  This makes things a bit more straight forward.
std::vector<Vector_t> MCSSegmentCalculator::SetPolygonalSegmentVector_vec(const std::vector<Barycenter_t> barycenter_vec, const Point_t startPoint) const {
  std::vector<Vector_t> polygonalSegmentVector_vec;

  Barycenter_t previousBarycenter = barycenter_vec.at(0);

  Vector_t firstSegment = previousBarycenter - startPoint.underlying();
  polygonalSegmentVector_vec.push_back(firstSegment);

  for(size_t i = 1; i < barycenter_vec.size(); ++i) {
    Barycenter_t currentBarycenter = barycenter_vec.at(i);
    Vector_t currentSegment = currentBarycenter - previousBarycenter;
    polygonalSegmentVector_vec.push_back(currentSegment);

    previousBarycenter = currentBarycenter;
  } // for i

  return polygonalSegmentVector_vec;

} // MCSSegmentCalculator::SetPolygonalSegmentVector_vec

template<typename Point>
Barycenter_t MCSSegmentCalculator::CalculateBarycenter(const std::vector<Point> segment) const {
  // Declare initial sums of the componenets (for averaging)
  Double_t xSum = 0, ySum = 0, zSum = 0;

  // Loop through the points in this segment
  for(size_t i = 0; i < segment.size(); ++i) {

    // Define the Point_t of this point
    Point_t trajPoint = segment.at(i);
    
    xSum += trajPoint.X();
    ySum += trajPoint.Y();
    zSum += trajPoint.Z();

  } // for i

  // Calculate the average of each component
  Double_t n = segment.size();  
  Double_t xAvg = xSum / n;
  Double_t yAvg = ySum / n;
  Double_t zAvg = zSum / n;

  // Define the Barycenter for this segment.
  Barycenter_t barycenter = Barycenter_t(xAvg, yAvg, zAvg);

  // Return the barycenter.
  return barycenter;

} // MCSSegmentCalculator::CalculateBarycenter

template<typename Point>
Vector_t MCSSegmentCalculator::CalculateLinearFit(const std::vector<Point> segment) const {

  // Calculate the barycenter of this segment.
  Barycenter_t barycenter = CalculateBarycenter(segment);

  // Declare initial sums as 0. (used in linear-regression)
  Float_t SSxz = 0, SSxx = 0, SSzy = 0, SSzz = 0;

  // Loop through the points in this segment.
  for(size_t k = 0; k < segment.size(); ++k) {
    // Calculate values used in sum calculations.
    Point_t trajPoint = segment.at(k);
    Double_t xDiff = trajPoint.X() - barycenter.X();
    Double_t yDiff = trajPoint.Y() - barycenter.Y();
    Double_t zDiff = trajPoint.Z() - barycenter.Z();

    SSxz += xDiff * zDiff;
    SSxx += xDiff * xDiff;
    SSzy += yDiff * zDiff;
    SSzz += zDiff * zDiff;
  } // for k

  // Slopes of the linear fit, projected onto two planes
  Double_t A = SSxz / SSxx; // Source: z = Ax + B
  Double_t C = SSzy / SSzz; // Source: y = Cz + D

  // Calculate the components of the linear fit.
  Double_t x = 1;
  Double_t y = A*C;
  Double_t z = A;

  // Calculate linear fit vector, transform to unit length.
  Vector_t linearFit = Vector_t(x, y, z).Unit();

  // Return the linear fit vector.
  return linearFit;
} // MCSSegmentCalculator::CalculateLinearFit

