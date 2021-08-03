#include "MCSAngleCalculator.h"

using namespace trkf;

using Point_t = MCSAngleCalculator::Point_t;
using Vector_t = MCSAngleCalculator::Vector_t;
using Angles_t = MCSAngleCalculator::Angles_t;
using MCSAngleResult = MCSAngleCalculator::MCSAngleResult;
using MCSSegmentResult = MCSAngleCalculator::MCSSegmentResult;

/// @todo Add assertions that the relative sizes of vectors that are created are correct.
/// @todo Add more segmentation methods to the segmentMethod, once they are added to the `MCSSegmentCalculator`
MCSAngleResult MCSAngleCalculator::GetResult(const MCSSegmentResult segmentResult, const int segmentMethod) const {

  // These will be set in the following switch block.
  std::vector<Float_t> segmentLength_vec;
  std::vector<Vector_t> segmentVector_vec;
  std::vector<Angles_t> angles_vec;

  // Depending on the segmentMethod, retrieve the correct data from the segment result.
  // In the default case, the linear segment method will be used and a warning will be printed to the console.
  switch(segmentMethod) {
  default:
    std::cout << "WARNING: MCSMomentumCalculator::GetResult's method parameter is incorrectly defined.  Please specify a correct method.  The linear method will be used as a default." << std::endl;
      // Purposely no break; Will continue to linear method (case 0).
  case 0:
    {
      // Linear Segment Method
      segmentLength_vec = segmentResult.GetRawSegmentLength_vec();
      segmentVector_vec = segmentResult.GetLinearFit_vec();
      angles_vec = SetAngles_vec(segmentVector_vec);
    } // case 0
    break;
  case 1: 
    {
      // Polygonal Segment Method
      segmentLength_vec = segmentResult.GetPolygonalSegmentLength_vec();
      segmentVector_vec = segmentResult.GetPolygonalSegmentVector_vec();
      angles_vec = SetAngles_vec(segmentVector_vec);
    } // case 1
    break;
  } // switch method

  // Create the angle result object
  MCSAngleResult angleResult = MCSAngleResult(segmentMethod, segmentResult, segmentLength_vec, segmentVector_vec, angles_vec);

  // Return the MCSAngleResult
  return angleResult;
} // MCSAngleCalculator::GetResult(const MCSSegmentResult, const int)

MCSAngleResult MCSAngleCalculator::GetResult(const std::vector<Point_t> trajPoint_vec, const int segmentMethod) const {
  // Get the segment result using the provided vector of trajectory points.
  MCSSegmentResult segmentResult = fSegmentCalculator.GetResult(trajPoint_vec, false);

  // Get the angle result using the calculated segment result.
  MCSAngleResult angleResult = GetResult(segmentResult, segmentMethod);

  // Return the angle result
  return angleResult;
} // MCSAngleCalculator::GetResult(const std::vector<Point_t>, const int)

MCSAngleResult MCSAngleCalculator::GetResult(const simb::MCTrajectory trajectory, const int segmentMethod) const {

  // Define an empty vector of trajectory points.
  std::vector<Point_t> trajPoint_vec;

  // Loop through trajectory points.
  for(size_t i = 0; i < trajectory.size(); ++i) {
    // Extract position for each trajectory point
    TVector3 position = trajectory.Position(i).Vect();
    // Convert the point from a TVector3 to a Point_t
    Point_t point = convert<TVector3, Point_t>(position);
    point.SetP(trajectory.Momentum(i).P());

    // If the point is in the detector, add it to the vector of trajectory points
    if(inDetector(point))
      trajPoint_vec.push_back(point);
  } // for i

  // Get Result
  MCSAngleResult result = GetResult(trajPoint_vec, segmentMethod);

  // Return Result
  return result;
} // MCSAngleCalculator::GetResult(const simb::MCTrajectory, const int)

// Transform the currentSegment into the coordinate system where the z-axis lies on the previousSegment
TVector3 MCSAngleCalculator::transform(TVector3 currentSegment, TVector3 previousSegment) const {

  // Define the rotation that will occur as an angle and a vector to rotate around.
  Double_t polarAngle = previousSegment.Theta(); // The angle the currentSegment will be rotating to form vPrime.
  TVector3 rotateAround = TVector3(0,0,1).Cross(previousSegment).Unit(); // The unit vector that the currentSegment will rotate around to form vPrime, Perpendicular to the previousSegment and the z-axis.
  
  // Initialize vPrime and rotate.
  TVector3 vPrime = currentSegment;
  vPrime.Rotate(-polarAngle, rotateAround);

  return vPrime;
} // transform

// Calculate the angles of the currentSegment projected onto the previousSegment.
Angles_t MCSAngleCalculator::calculateAngles(TVector3 currentSegment, TVector3 previousSegment) const {
  // vPrime is now the currentSegment defined in the coordinate system where the previousSegment lies on the z-axis.
  Float_t theta3D = previousSegment.Angle(currentSegment);

  TVector3 vPrime = transform(currentSegment, previousSegment);

  // Calculate projection angles.
  Double_t thetaXZprime = TMath::ATan2(vPrime.X(), vPrime.Z());
  Double_t thetaYZprime = TMath::ATan2(vPrime.Y(), vPrime.Z());
  
  // Return std::pair of projection angles.
  return std::make_tuple(theta3D, thetaXZprime, thetaYZprime);
}

/// @details If n-1 is the size of the segmentVector_vec, the returned vector should be of size n-2.
/// This is because 2 segment vectors are required to calculate the angles between them.
/// @todo Add support for an optional start-vector, so that the first segment is measured with respect to that.
std::vector<Angles_t> MCSAngleCalculator::SetAngles_vec(std::vector<Vector_t> segmentVector_vec) const {

  // Define the angles_vec
  std::vector<Angles_t> angles_vec;

  // Initial previousSegment
  TVector3 previousSegment = convert<Vector_t, TVector3>(segmentVector_vec.at(0));

  // Loop through the segments
  for(size_t i = 1; i < segmentVector_vec.size(); ++i) {
    // Vector of the current segment
    TVector3 currentSegment = convert<Vector_t, TVector3>(segmentVector_vec.at(i));

    // Calculate the angles of the currentSegment projected into the reference frame of the previousSegment.
    Angles_t angles = calculateAngles(currentSegment, previousSegment);

    // Fill the angles_vec
    angles_vec.push_back(angles);

    // Update the previouSegment
    previousSegment = currentSegment;
  } // for i

  // Return the angles_vec
  return angles_vec;

} // MCSAngleCalculator::SetAngles_vec
