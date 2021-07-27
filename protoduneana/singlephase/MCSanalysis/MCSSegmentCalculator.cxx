#include "MCSSegmentCalculator.h"

using namespace trkf;

using Point_t = MCSSegmentCalculator::Point_t;
using Vector_t = MCSSegmentCalculator::Vector_t;
using Barycenter_t = MCSSegmentCalculator::Barycenter_t;
using MCSSegmentResult = MCSSegmentCalculator::MCSSegmentResult;
using Segment_t = MCSSegmentCalculator::Segment_t;

/**
Create an MCSSegmentResult using a vector of points, optionally create virtual points to force the segment length to be 14-cm.
 */
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

  // TODO: Allow for other segment defintions.

  // TODO: Add assertions/exceptions to make sure that the relative sizes of these vectors are correct.

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

// TODO: Documentation
MCSSegmentResult MCSSegmentCalculator::GetResult(const simb::MCTrajectory trajectory, const Bool_t shouldCreateVirtualPoints) const {

  std::vector<Point_t> trajPoint_vec;
  for(size_t i = 0; i < trajectory.size(); ++i) {
    TVector3 position = trajectory.Position(i).Vect();
    Point_t point = convert<TVector3, Point_t>(position);
    point.SetP(trajectory.Momentum(i).P());
    if(inDetector(point))
      trajPoint_vec.push_back(point);
  }

  MCSSegmentResult result = GetResult(trajPoint_vec, shouldCreateVirtualPoints);

  return result;
} // MCSSegmentCalculator::GetResult(const simb::MCTrajectory, const Bool_t)

// TODO: Documentation
// This will create segments, optionally with virtual points.
std::vector<Segment_t> MCSSegmentCalculator::CreateSegment_vec(const std::vector<Point_t> points, Bool_t shouldCreateVirtualPoints) const {
  std::vector<Segment_t> segments;
  Float_t currLength = 0; // cm
  Float_t prevLength = 0; // cm

  Segment_t currSegment;
  
  // Loop through the trajectory points
  // TODO: Consider change to NextValidPoint methodology.
  for(size_t i = 1; i < points.size(); ++i) {
    // Define the current and previous trajectory point.
    Point_t prevPoint = points.at(i-1);
    Point_t currPoint = points.at(i);

    // TODO: Add an assertion that prevLength is never greater than fSegmentLength.
    // TODO: More assertions of what should never happen?
    // TODO: How much of this function can be broken into smaller, more testable functions? Maybe none.

    // If this point is invalid.
    if(currPoint.X() == -999)
      break;

    // Define the difference between the current and previous trajectory point.
    Float_t diff = sqrt((currPoint.underlying() - prevPoint.underlying()).Mag2());
    // Update the length.
    currLength += diff;

    // TODO: This algorithm could be generally rewritten.  There may be some assumptions being made that no longer apply or we could possibly change the order of some of these cases to reduce the amount of code needed.  
    // TODO: Make a list of the possible scenarios and play through what SHOULD happen in each of them, then go through an see what actually would happen.

    // TODO: Should this be diff >= or currLength >=?
    if(diff >= 2*fSegmentLength) {
      // If the current length is greater than 2*fSegmentLength, then a special case has occured (kinda rare).
      if(shouldCreateVirtualPoints) {
	// If shouldCreateVirtualPoints, multiple virtual points should be created.
	// TODO: Currently, only two are made.  But technically, more could be made if the gap between the trajectory points was VERY large.

	// Create the first virtual point
	Float_t multiplier1 = (fSegmentLength - prevLength) / diff;
	Float_t x1 = multiplier1*(currPoint.X() - prevPoint.X()) + prevPoint.X();
	Float_t y1 = multiplier1*(currPoint.Y() - prevPoint.Y()) + prevPoint.Y();
	Float_t z1 = multiplier1*(currPoint.Z() - prevPoint.Z()) + prevPoint.Z();
	Point_t virtualPoint1(x1, y1, z1);
	// Set the first virtual point's momentum.  TODO: Should this be set as the closest point's momentum, the average momentum of the points on either side, or the previous point's momentum?
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

    // TODO: Possible assertions here:
    // 1. The number of points is not zero.
    // 2. The currLength is 0 if the number of points is 1 and greater than 0 if the number of points is greater than 1.

    // Update the oldLength.
    prevLength = currLength;
  } // for i

  // Add the last segment, which never reached 14-cm.
  segments.push_back(currSegment);

  return segments;
} // MCSSegmentCalculator::CreateSegment_vec

// TODO: Documentation
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

// TODO: Documentation
std::vector<Float_t> MCSSegmentCalculator::SetRawSegmentLength_vec(const std::vector<Segment_t> segment_vec) const {
  std::vector<Float_t> rawSegmentLength_vec;
  
  for(size_t i = 0; i < segment_vec.size() - 1; ++i) {
    Segment_t segment = segment_vec.at(i);
    Float_t rawSegmentLength = CalculateRawSegmentLength(segment);
    rawSegmentLength_vec.push_back(rawSegmentLength);
  }
  
  return rawSegmentLength_vec;
} // MCSSegmentCalculator:SetRawSegmentLength_vec

// TODO: Documentation
// std::vector<size_t> MCSSegmentCalculator::SetSegmentStartPoint_vec(const std::vector<Point_t> trajPoint_vec, std::vector<Float_t>& rawSegmentLength_vec) const {
//   std::vector<size_t> segmentStartPoint_vec;
//   Float_t length = 0; //cm
//   Float_t oldLength = 0; //cm
//   Float_t diff;
//   segmentStartPoint_vec.push_back(0); // First point in detector
//   // Loop through trajectory points to form segments as close to segmentLength as possible.
//   for(size_t i = 1; i < trajPoint_vec.size(); i++) {
//     // Once you get to the invalid trajectory points at the end, break.
//     if(trajPoint_vec.at(i).X() == -999)
//       break;

//     // Calculate the distance between two points, add that to the current raw segment length.
//     diff = pow((trajPoint_vec.at(i).underlying()-trajPoint_vec.at(i-1).underlying()).Mag2(),0.5);
//     length += diff;

//     if(diff >= 2*fSegmentLength) {
//       segmentStartPoint_vec.push_back(i);
//       rawSegmentLength_vec.push_back(oldLength);
//       length = 0;
//     } else if((oldLength < fSegmentLength && length >= fSegmentLength) || oldLength > fSegmentLength) {
//       if(fabs(oldLength - fSegmentLength) <= fabs(length - fSegmentLength)) {
// 	segmentStartPoint_vec.push_back(i-1);
// 	rawSegmentLength_vec.push_back(oldLength);
// 	length-=oldLength;
//       } else {
// 	segmentStartPoint_vec.push_back(i);
// 	rawSegmentLength_vec.push_back(length);
// 	length = 0;
//       }
//     }
//     // Update the oldLength
//     oldLength = length;
//   } // for i if less than # of traj pts

//   return segmentStartPoint_vec;
//   /* Structure of vector<int> segmentStartPoint_vec

//      N = segmentStartPoint_vec.size() = The number of full-length raw segments + 1.

//      index of vector, i, goes from 0 to N-1 where N is the number of segments.

//      index 0 is of the first point in the detector, the start point of the first raw segment.
//      index N-1 is the start point of the partial segment at the end. (not 14-cm)
//   */
// }

// // TODO: Documentation
// std::vector<Barycenter_t> MCSSegmentCalculator::SetBarycenter_vec(const std::vector<Point_t> trajPoint_vec, std::vector<size_t> segmentStartPoint_vec) const {

//   // Declare an initial vector of barycenters
//   std::vector<Barycenter_t> barycenter_vec;

//   // Loop through segment start points
//   for(size_t i = 0; i < segmentStartPoint_vec.size()-1; ++i) {

//     // Declare the list of points in this segment.
//     std::vector<Point_t> segment;
//     for(size_t j = segmentStartPoint_vec.at(i); j <= segmentStartPoint_vec.at(i+1); ++j)
//       segment.push_back(trajPoint_vec.at(j));

//     // Calculate the barycenter for this segment.
//     Barycenter_t barycenter = CalculateBarycenter(segment);
    
//     // Fill the barycenter_vec
//     barycenter_vec.push_back(barycenter);
//   } // for i

//   // Return the barycenter_vec.
//   return barycenter_vec;
//   // TODO: Update this for new Barycenter_t type.
//   /* Structure of barycenters vector<std::array<double,3> >

//      N = segmentStartPoint_vec.size() = The number of full-length raw segments + 1.

//      x of i'th barycenter = barycenters.at(i)[0]
//      y of i'th barycenter = barycenters.at(i)[1]
//      z of i'th barycenter = barycenters.at(i)[2]
//      index of vector, i, goes from 0 to N-1 where N is the number of segments.

//      The first barycenter is the barycenter of the first segment, i = 0 and is the average xyz of all trajectory points from segmentStartPoints[0] to [1]

//   */
// } // MCSSegmentCalculator::SetBarycenter_vec

// TODO: Documentation
std::vector<Barycenter_t> MCSSegmentCalculator::SetBarycenter_vec(const std::vector<Segment_t> segment_vec) const {
  std::vector<Barycenter_t> barycenter_vec;

  for(size_t i = 0; i < segment_vec.size() - 1; ++i) {
    Segment_t segment = segment_vec.at(i);
    Barycenter_t barycenter = CalculateBarycenter(segment);
    barycenter_vec.push_back(barycenter);
  }

  return barycenter_vec;
} // SetBarycenter_vec

// TODO: Documentation
// TODO: Update this to use distanceBetween (and move that method to this class)
std::vector<Float_t> MCSSegmentCalculator::SetPolygonalSegmentLength_vec(const std::vector<Barycenter_t> barycenter_vec, Point_t startPoint) const {
  std::vector<Float_t> polygonalSegmentLength_vec;
  Barycenter_t previousBarycenter = barycenter_vec.at(0);
  Float_t x = previousBarycenter.X() - startPoint.X();
  Float_t y = previousBarycenter.Y() - startPoint.Y();
  Float_t z = previousBarycenter.Z() - startPoint.Z();
  polygonalSegmentLength_vec.push_back(sqrt(x*x + y*y + z*z));
  for(size_t i = 1; i < barycenter_vec.size(); i++) {
    Barycenter_t currentBarycenter = barycenter_vec.at(i);

    x = currentBarycenter.X() - previousBarycenter.X();
    y = currentBarycenter.Y() - previousBarycenter.Y();
    z = currentBarycenter.Z() - previousBarycenter.Z();

    // Push back length
    polygonalSegmentLength_vec.push_back(sqrt(x*x + y*y + z*z));

    previousBarycenter = currentBarycenter;
  }

  return polygonalSegmentLength_vec;
}

// TODO: Documentation
// std::vector<Vector_t> MCSSegmentCalculator::SetLinearFit_vec(const std::vector<Point_t> trajPoint_vec, std::vector<size_t> segmentStartPoint_vec) const {

//   // Create an initial linearFit_vec.
//   std::vector<Vector_t> linearFit_vec;

//   // Loop through the segment start points.
//   for(size_t i = 0; i < segmentStartPoint_vec.size() - 1; i++) {

//     // Create the vector of Point_t for this segment.
//     std::vector<Point_t> segment;
//     for(size_t j = segmentStartPoint_vec.at(i); j <= segmentStartPoint_vec.at(i+1); ++j)
//       segment.push_back(trajPoint_vec.at(j));

//     // Calculate the linear fit for this segment
//     Vector_t linearFit = CalculateLinearFit(segment);
//     linearFit_vec.push_back(linearFit);

//   } // for i

//   // Return the vector of linear fits.
//   return linearFit_vec;
// }

// TODO: Documentation
std::vector<Vector_t> MCSSegmentCalculator::SetLinearSegmentVector_vec(const std::vector<Segment_t> segment_vec) const {
  std::vector<Vector_t> linearFit_vec;

  for(size_t i = 0; i < segment_vec.size() - 1; ++i) {
    Segment_t segment = segment_vec.at(i);
    Vector_t linearFit = CalculateLinearFit(segment);
    linearFit_vec.push_back(linearFit);
  }

  return linearFit_vec;
} // MCSSegmentCalculator::SetLinearSegmentVector_vec

// TODO: Documenation
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

// TODO: Documentation
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

// TODO: Documentation
// This will take the linear fit of the provided trajectory points.
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

