#include "MCSAngleCalculator.h"

using namespace trkf;

// TODO: Documentation
using Point_t = MCSAngleCalculator::Point_t;
using Vector_t = MCSAngleCalculator::Vector_t;
using Angles_t = MCSAngleCalculator::Angles_t;
using MCSAngleResult = MCSAngleCalculator::MCSAngleResult;
using MCSSegmentResult = MCSAngleCalculator::MCSSegmentResult;

// Get an MCSAngleResult using the provided MCSSegmentResult.
MCSAngleResult MCSAngleCalculator::GetResult(const MCSSegmentResult segmentResult, const int segmentMethod) const {

  // TODO: Run MCSSegmentResult assertions with proper error handling.

  // These will be set in the following switch block.
  std::vector<Float_t> segmentLength_vec;
  std::vector<Vector_t> segmentVector_vec;
  std::vector<Angles_t> angles_vec;

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
    // TODO: Add more cases for more segmentation methods (last point, crossover linear fit, straight line between first and last point, etc.).
  } // switch method

  // Create the angle result object
  MCSAngleResult angleResult = MCSAngleResult(segmentMethod, segmentResult, segmentLength_vec, segmentVector_vec, angles_vec);

  // TODO: Run MCSAngleResult assertions with proper error handling.

  // Return the MCSAngleResult
  return angleResult;
} // MCSAngleCalculator::GetResult(const MCSSegmentResult, const int)

// Get an MCSAngleResult using the provided vector of trajectory points.
// segmentMethod:
//     0: Linear Fit Segment Method
//     1: Polygonal Segment Method
//     TODO: Add more methods
MCSAngleResult MCSAngleCalculator::GetResult(const std::vector<Point_t> trajPoint_vec, const int segmentMethod) const {
  // Get the segment result using the provided vector of trajectory points.
  MCSSegmentResult segmentResult = fSegmentCalculator.GetResult(trajPoint_vec, false);

  // Get the angle result using the calculated segment result.
  MCSAngleResult angleResult = GetResult(segmentResult, segmentMethod);

  // Return the angle result
  return angleResult;
} // MCSAngleCalculator::GetResult(const std::vector<Point_t>, const int)

/*
// TODO: Documentation
MCSAngleResult MCSAngleCalculator::GetResult(const std::vector<Point_t> trajPoint_vec, const int segmentMethod) const {

  // TODO: Change to use segmentResult as a function parameter, add another function parameter that can be told to just use a method.
  //       This will make it so that not all of the methods need to be run over and over again.
  //       Eventually, it will make more sense as well to only perform the segmentation for the provided segmentation method.  Currently, all segments are formed but only the specified angles are calculated and only the momentum result is calculated for the specified method.
  MCSSegmentResult segmentResult = fSegmentCalculator.GetResult(trajPoint_vec, false); // TODO: The false specifies that the result should not use virtual points.

  // TODO: Begin by adding checks, else return a bad result or throw an error.

  // These will be set in the following switch block.
  std::vector<Float_t> segmentLength_vec;
  std::vector<Vector_t> segmentVector_vec;
  std::vector<Angles_t> angles_vec;

  switch(segmentMethod) {
  default:
    std::cout << "WARNING: MCSMomentumCalculator::GetResult's method parameter is incorrectly defined.  Please specify a correct method.  The linear method will be used as a default." << std::endl;
      // Purposely no break; Will use linear method.
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
    // TODO: Add more cases for more segmentation (last point, crossover linear fit, straight line between first and last point, etc.) methods.
  } // swtich method

  // std::cout << "segmentMethod = " << segmentMethod << std::endl;
  // std::cout << "segmentVector_vec.size() = " << segmentVector_vec.size() << std::endl;
  // std::cout << "segmentLength-vec.size() = " << segmentLength_vec.size() << std::endl;
  // std::cout << "angles_vec.size() = " << angles_vec.size() << std::endl;


  // TODO: Add assertions to avoid issues with size or things not filling properly.  Serves as an inline test.

  MCSAngleResult result = MCSAngleResult(segmentMethod, segmentResult, segmentLength_vec, segmentVector_vec, angles_vec);

  // TODO: Assertions and/or tests of the result?

  // Return the MCSAngleResult
  return result;
}
*/

// Get an MCSAngleResult for the provided MCTrajectory
// TODO: Add segmentMethod documentation
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

// template<typename Point>
// Bool_t MCSAngleCalculator::inDetector(Point point) const {
//   return fSegmentCalculator.inDetector(point);
// }

// template<typename From, typename To>
// To MCSAngleCalculator::convert(From vector) const {
//   return fSegmentCalculator.convert<From, To>(vector);
// }

/*
template<typename From, typename To>
To MCSAngleCalculator::convert(From vector) const {
  return To(vector.X(), vector.Y(), vector.Z());
}

// Determine if the point provided is in the detector coordinates.
// TODO: Surely there's a better way to do this using the detector geometry?
template<typename Point>
Bool_t MCSAngleCalculator::inDetector(Point point) const {
  return 
    point.Z() <= 700 && point.Z() >= 0 &&
    point.X() <= 300  && point.X() >= -300 &&
    point.Y() <= 600 && point.Y() >= 0;
}
*/

// Transform the currentSegment into the coordinate system where the z-axis lies on the previousSegment
TVector3 MCSAngleCalculator::transform(TVector3 currentSegment, TVector3 previousSegment) const {

  // Define the rotation that will occur as an angle and a vector to rotate around.
  Double_t polarAngle = previousSegment.Theta(); // The angle the currentSegment will be rotating to form vPrime.
  TVector3 rotateAround = TVector3(0,0,1).Cross(previousSegment).Unit(); // The unit vector that the currentSegment will rotate around to form vPrime, Perpendicular to the previousSegment and the z-axis.
  
  // Initialize vPrime and rotate.
  TVector3 vPrime = currentSegment;
  vPrime.Rotate(-polarAngle, rotateAround);

  // A random rotation around zPrime, from 0 -> Pi/2
  /*
  TRandom3 rng;
  rng.SetSeed();
  Double_t angle = rng.Rndm()*TMath::PiOver2();
  vPrime.RotateZ(angle);
  */
  
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

// TODO: Documentation
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

// // TODO: Documentation
// std::vector<size_t> MCSAngleCalculator::SetSegmentStartPoint_vec(std::vector<Point_t> trajPoint_vec, Float_t segmentLength, std::vector<Float_t>& rawSegmentLengths_vec) const {
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
//     diff = pow((trajPoint_vec.at(i)-trajPoint_vec.at(i-1)).Mag2(),0.5);
//     length += diff;
//     if(diff >= 2*segmentLength) {
//       segmentStartPoint_vec.push_back(i);
//       rawSegmentLengths_vec.push_back(oldLength);
//       length = 0;
//     } else if((oldLength < segmentLength && length >= segmentLength) || oldLength > segmentLength) {
//       if(fabs(oldLength - segmentLength) <= fabs(length-segmentLength)) {
// 	segmentStartPoint_vec.push_back(i-1);
// 	rawSegmentLengths_vec.push_back(oldLength);
// 	length-=oldLength;
//       } else {
// 	segmentStartPoint_vec.push_back(i);
// 	rawSegmentLengths_vec.push_back(length);
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
// std::vector<Barycenter_t> MCSAngleCalculator::SetBarycenter_vec(std::vector<Point_t> trajPoint_vec, std::vector<size_t> segmentStartPoint_vec) const {
//   std::vector<Barycenter_t> barycenter_vec;
//   double xSum, ySum, zSum;
//   // Loop through segments
//   for(size_t i = 0; i < segmentStartPoint_vec.size()-1; i++) {
//     Barycenter_t xyz;
//     xSum = 0;
//     ySum = 0;
//     zSum = 0;
//     // Loop through trajectory points in this segment
//     for(size_t j = segmentStartPoint_vec.at(i); j <= segmentStartPoint_vec.at(i+1); j++) {
//       xSum += trajPoint_vec.at(j).X();
//       ySum += trajPoint_vec.at(j).Y();
//       zSum += trajPoint_vec.at(j).Z();
//     } // looping through each point in each segment
//     Double_t xAvg = xSum / (segmentStartPoint_vec.at(i+1) - segmentStartPoint_vec.at(i)+1); // average x of this segment
//     Double_t yAvg = ySum / (segmentStartPoint_vec.at(i+1) - segmentStartPoint_vec.at(i)+1); // average y of this segment
//     Double_t zAvg = zSum / (segmentStartPoint_vec.at(i+1) - segmentStartPoint_vec.at(i)+1); // average z of this segment
//     xyz.SetXYZ(xAvg, yAvg, zAvg);

//     barycenter_vec.push_back(xyz);  // The x,y,and z array matched to the i'th segment's barycenter
//   } // looping through each segment

//   return barycenter_vec;
//   /* Structure of barycenters vector<std::array<double,3> >

//      N = segmentStartPoint_vec.size() = The number of full-length raw segments + 1.

//      x of i'th barycenter = barycenters.at(i)[0]
//      y of i'th barycenter = barycenters.at(i)[1]
//      z of i'th barycenter = barycenters.at(i)[2]
//      index of vector, i, goes from 0 to N-1 where N is the number of segments.

//      The first barycenter is the barycenter of the first segment, i = 0 and is the average xyz of all trajectory points from segmentStartPoints[0] to [1]
	  
//   */
// }

// // TODO: Documentation
// std::vector<Float_t> MCSAngleCalculator::SetPolygonalSegmentLength_vec(std::vector<Barycenter_t> barycenter_vec, Point_t startPoint) const {
//   std::vector<Float_t> polygonalSegmentLength_vec;
//   Barycenter_t previousBarycenter = barycenter_vec.at(0);
//   Float_t x = previousBarycenter.X() - startPoint.X();
//   Float_t y = previousBarycenter.Y() - startPoint.Y();
//   Float_t z = previousBarycenter.Z() - startPoint.Z();
//   polygonalSegmentLength_vec.push_back(sqrt(x*x + y*y + z*z));
//   for(size_t i = 1; i < barycenter_vec.size(); i++) {
//     Barycenter_t currentBarycenter = barycenter_vec.at(i);

//     x = currentBarycenter.X() - previousBarycenter.X();
//     y = currentBarycenter.Y() - previousBarycenter.Y();
//     z = currentBarycenter.Z() - previousBarycenter.Z();

//     // Push back length
//     polygonalSegmentLength_vec.push_back(sqrt(x*x + y*y + z*z));

//     previousBarycenter = currentBarycenter;
//   }

//   return polygonalSegmentLength_vec;
// }

// // TODO: Documentation
// std::vector<Vector_t> MCSAngleCalculator::SetLinearFit_vec(std::vector<Point_t> trajPoint_vec, std::vector<size_t> segmentStartPoint_vec) const {
//   std::vector<Vector_t> linearFit_vec;
//   for(size_t i = 0; i < segmentStartPoint_vec.size() - 1; i++) {
//     // segmentStartPoint_vec.at(segmentStartPoint_vec.size()-1) = the last segment start point which ISN'T a full segment and shouldn't be counted in linear fit angles and what not.

//     Float_t xAvg=0, yAvg=0, zAvg=0, SSxz = 0, SSxx = 0, SSzy = 0, SSzz = 0/*, SSyy = 0, rSq_xz, rSq_zy*/;
//     // Calculate linear fit using linear regression for each segment.
//     for(size_t j = segmentStartPoint_vec.at(i); j <= segmentStartPoint_vec.at(i+1); j++) {
//       xAvg += trajPoint_vec.at(j).X();
//       yAvg += trajPoint_vec.at(j).Y();
//       zAvg += trajPoint_vec.at(j).Z();
//     }
//     size_t n = segmentStartPoint_vec.at(i+1)-segmentStartPoint_vec.at(i)+1;
//     xAvg /= n;
//     yAvg /= n;
//     zAvg /= n;
//     for(size_t k = segmentStartPoint_vec.at(i); k <= segmentStartPoint_vec.at(i+1); k++) {
//       SSxz += ((trajPoint_vec.at(k).X()-xAvg)*(trajPoint_vec.at(k).Z()-zAvg));
//       SSxx += pow((trajPoint_vec.at(k).X()-xAvg),2);
//       SSzy += ((trajPoint_vec.at(k).Y()-yAvg)*(trajPoint_vec.at(k).Z()-zAvg));
//       SSzz += pow((trajPoint_vec.at(k).Z()-zAvg),2);
//     }

//     // Slopes of the projected linear fit
//     Double_t A = SSxz / SSxx; // Slope of z = Ax + B
//     Double_t C = SSzy / SSzz; // Slope of y = Cz + D

//     // Calculate 
//     Double_t x = 1;
//     Double_t y = A*C;
//     Double_t z = A;
    
//     // Push_back unit vector that is the linear fit
//     linearFit_vec.push_back(Vector_t(x, y, z).Unit());
//   }

//   return linearFit_vec;
// }

// NOTE TO SELF: DO NOT UNCOMMENT, THIS IS AN OLD IMPLEMENTATION
// TODO: Documentation
// std::vector<AnglePair_t> MCSAngleCalculator::SetLinearAngles_vec(std::vector<Vector_t> linearFit_vec) const {
//   // Define linearAngles_vec
//   std::vector<AnglePair_t> linearAngles_vec;

//   // Initial previousSegment (a linear fit of the first segment)
//   TVector3 previousSegment = convert(linearFit_vec.at(0));
//   // Loop through the segments
//   for(size_t i = 1; i < linearFit_vec.size(); ++i) {
//     // A linear fit of the current segment
//     TVector3 currentSegment = convert(linearFit_vec.at(i));
//     // Calculate the angles of the currentSegment projected into the reference frame of the previousSegment.
//     AnglePair_t anglePair = calculateAnglePair(currentSegment, previousSegment);
    
//     // Fill linearAngles_vec
//     linearAngles_vec.push_back(anglePair);

//     // Update the previousSegment
//     previousSegment = currentSegment;
//   }

//   return linearAngles_vec;
// }

// std::vector<Angles_t> MCSAngleCalculator::SetLinearAngles_vec(std::vector<Vector_t> linearFit_vec) const {
//   // Define linearAngles_vec
//   std::vector<Angles_t> linearAngles_vec;

//   // Initial previousSegment (a linear fit of the first segment)
//   TVector3 previousSegment = convert(linearFit_vec.at(0));
//   // Loop through the segments
//   for(size_t i = 1; i < linearFit_vec.size(); ++i) {
//     // A linear fit of the current segment
//     TVector3 currentSegment = convert(linearFit_vec.at(i));
//     // Calculate the angles of the currentSegment projected into the reference frame of the previousSegment.
//     Angles_t linearAngles = calculateAngles(currentSegment, previousSegment);
    
//     // Fill linearAngles_vec
//     linearAngles_vec.push_back(linearAngles);

//     // Update the previousSegment
//     previousSegment = currentSegment;
//   }

//   return linearAngles_vec;
// }

// TODO: Documentation
// std::vector<Float_t> MCSAngleCalculator::SetRawAngle_vec(std::vector<Vector_t> linearFit_vec) const {
//   std::vector<Float_t> rawAngle_vec;

//   // Initial Segment Vector
//   TVector3 previousSegment = convert(linearFit_vec.at(0));
//   // Loop through segments
//   for(size_t i = 1; i < linearFit_vec.size(); ++i) {
//     TVector3 currentSegment = convert(linearFit_vec.at(i));

//     // Calculate the raw angle
//     Float_t rawAngle = previousSegment.Angle(currentSegment);

//     // Fill rawAngle_vec
//     rawAngle_vec.push_back(rawAngle);

//     // Update the previousSegment
//     previousSegment = currentSegment;
//   }

//   return rawAngle_vec;
// }

// TODO: Documentation
// std::vector<AnglePair_t> MCSAngleCalculator::SetPolygonalAngles_vec(std::vector<Barycenter_t> barycenter_vec, Point_t startPoint) const {
//   // Define polygonalAngles_vec
//   std::vector<AnglePair_t> polygonalAngles_vec;

//   Barycenter_t previousBarycenter = barycenter_vec.at(0);
//   TVector3 previousSegment(previousBarycenter.X()-startPoint.X(), previousBarycenter.Y()-startPoint.Y(), previousBarycenter.Z()-startPoint.Z());
//   // Loop through the segments
//   for(size_t i = 1; i < barycenter_vec.size(); ++i) {
//     // The Current barycenter
//     Barycenter_t currentBarycenter = barycenter_vec.at(i);
//     // Vector connecting the barycenters
//     TVector3 currentSegment(currentBarycenter.X()-previousBarycenter.X(), currentBarycenter.Y()-previousBarycenter.Y(), currentBarycenter.Z()-previousBarycenter.Z());
//     // Calculate the angles of the currentSegment projected into the reference frame of the previousSegment
//     AnglePair_t anglePair = calculateAnglePair(currentSegment, previousSegment);

//     // Fill polygonalAngles_vec
//     polygonalAngles_vec.push_back(anglePair);

//     // Update the previousBarycenter and previousSegment
//     previousBarycenter = currentBarycenter;
//     previousSegment = currentSegment;
//   }

//   return polygonalAngles_vec;
// }

// std::vector<Angles_t> MCSAngleCalculator::SetPolygonalAngles_vec(std::vector<Barycenter_t> barycenter_vec, Point_t startPoint) const {
//   // Define polygonalAngles_vec
//   std::vector<Angles_t> polygonalAngles_vec;

//   Barycenter_t previousBarycenter = barycenter_vec.at(0);
//   TVector3 previousSegment(previousBarycenter.X()-startPoint.X(), previousBarycenter.Y()-startPoint.Y(), previousBarycenter.Z()-startPoint.Z());
//   // Loop through the segments
//   for(size_t i = 1; i < barycenter_vec.size(); ++i) {
//     // The Current barycenter
//     Barycenter_t currentBarycenter = barycenter_vec.at(i);
//     // Vector connecting the barycenters
//     TVector3 currentSegment(currentBarycenter.X()-previousBarycenter.X(), currentBarycenter.Y()-previousBarycenter.Y(), currentBarycenter.Z()-previousBarycenter.Z());
//     // Calculate the angles of the currentSegment projected into the reference frame of the previousSegment
//     Angles_t polygonalAngles = calculateAngles(currentSegment, previousSegment);

//     // Fill polygonalAngles_vec
//     polygonalAngles_vec.push_back(polygonalAngles);

//     // Update the previousBarycenter and previousSegment
//     previousBarycenter = currentBarycenter;
//     previousSegment = currentSegment;
//   }

//   return polygonalAngles_vec;
// }

// // Convert a Vector_t to a TVector3
// TVector3 MCSAngleCalculator::convert(Vector_t vector) const {
//   return TVector3(vector.X(), vector.Y(), vector.Z());
// }

// // Determine if the point provided is in the detector coordinates.
// // TODO: Surely there's a better way to do this using the detector geometry?
// template<typename Point>
// Bool_t MCSAngleCalculator::inDetector(Point point) const {
//   return 
//     point.Z() <= 700 && point.Z() >= 0 &&
//     point.X() <= 300  && point.X() >= -300 &&
//     point.Y() <= 600 && point.Y() >= 0;
// }
