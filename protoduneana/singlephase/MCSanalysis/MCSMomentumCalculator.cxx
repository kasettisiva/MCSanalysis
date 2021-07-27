#define _unused [[maybe_unused]]// Pre-processor directive in a C++ macro.

#include "MCSMomentumCalculator.h"

using namespace recob::tracking;
using namespace trkf;

using MCSAngleResult = MCSMomentumCalculator::MCSAngleResult;
using MCSSegmentResult = MCSMomentumCalculator::MCSSegmentResult;
using MCSMomentumResult = MCSMomentumCalculator::Result;
using AnglePair_t = MCSMomentumCalculator::AnglePair_t;
using Angles_t = MCSMomentumCalculator::Angles_t;

MCSMomentumResult MCSMomentumCalculator::GetResult(const MCSAngleResult angleResult, Int_t angleMethod) const {
  std::vector<Angles_t> angles_vec = angleResult.GetAngles_vec();
  std::vector<Float_t> segmentLength_vec = angleResult.GetSegmentLength_vec();

  // Will set in switch statement:
  std::vector<Float_t> logLikelihood_vec;

  switch(angleMethod) {
  default:
    std::cout << "WARNING: MCSMomentumCalculator::GetResult method is incorrectly defined.  Please you 0 for the linear method and 1 for the polygonal method. Going to default to the linear method.";
    // Purposely no break; Will use projected angle method.
    // TODO: Consider adding angleMethod AND segmentMethod to the GetResult function call.  Add as an enum?  Then we don't have to even have the option of the default case.
  case 0: case 1:
    {
      // Projected Angles Method
      std::vector<AnglePair_t> projectedAngles_vec = ExtractProjectedAnglePair_vec(angles_vec);
      logLikelihood_vec = SetLogLikelihood_vec(projectedAngles_vec, segmentLength_vec);
    } // case 0 or case 1
    break;
  case 2: case 3:
    {
      // 3D Angle Method
      std::vector<Float_t> theta3D_vec = ExtractTheta3D_vec(angles_vec);
      logLikelihood_vec = SetLogLikelihood_vec(theta3D_vec, segmentLength_vec);
    } // case 2 or case 3
    break;
  } // switch angleMethod

  // Calculate the MCSMomentum_vec
  std::vector<Float_t> MCSMomentum_vec = SetMCSMomentum_vec(logLikelihood_vec, segmentLength_vec);

  // Get the result
  MCSMomentumResult result = MCSMomentumResult(angleMethod, angleResult, segmentLength_vec, logLikelihood_vec, MCSMomentum_vec);

  // Return the result
  return result;
}

MCSMomentumResult MCSMomentumCalculator::GetResult(const MCSSegmentResult segmentResult, Int_t angleMethod) const {
  MCSAngleResult angleResult = SetAngleResultForAngleMethod(segmentResult, angleMethod);  
  MCSMomentumResult momentumResult = GetResult(angleResult, angleMethod);
  return momentumResult;
}

// Get an MCSMomentumCalculator::Result given a set of trajectory points and a method to use.
// Methods:
// linear projected angles: 0
// polygonal projected angles: 1
// linear 3D angles: 2
// polygonal 3D angles: 3
MCSMomentumResult MCSMomentumCalculator::GetResult(const std::vector<Point_t> trajPoint_vec, int angleMethod) const {

  // Calculate the MCS Angle Result.
  // TODO: This needs to change to parse the angleMethod and provide a segmentMethod to the MCSAngleCalculator.GetResult function.  Note, angleMethod and segmentMethod are DIFFERENT.
  // TODO: Allow the API user to optionally provide an MCSSegmentResult to this method?

  MCSAngleResult angleResult = SetAngleResultForAngleMethod(trajPoint_vec, angleMethod);
  MCSMomentumResult momentumResult = GetResult(angleResult, angleMethod);
  return momentumResult;

} // MCSMomentumCalculator::GetResult

// TODO: Documentation
MCSAngleResult MCSMomentumCalculator::SetAngleResultForAngleMethod(const MCSSegmentResult segmentResult, const Int_t angleMethod) const {
  switch(angleMethod) {
  default:
  case 0: case 2:
    return fAngleCalculator.GetResult(segmentResult, 0);
  case 1: case 3:
    return fAngleCalculator.GetResult(segmentResult, 1);
  } // switch angleMethod
} // MCSMomentumCalculator::SetAngleResultForAngleMethod(const MCSSegmentResult, const Int_t)

MCSAngleResult MCSMomentumCalculator::SetAngleResultForAngleMethod(const std::vector<Point_t> trajPoint_vec, const Int_t angleMethod) const {
  switch(angleMethod) {
  default:
    std::cout << "WARNING: MCSMomentumCalculator::SetAngleResultForAngleMethod: angleMethod is incorrectly defined.  Will default to the linear method." << std::endl;
  case 0: case 2:
    return fAngleCalculator.GetResult(trajPoint_vec, 0);
  case 1: case 3:
    return fAngleCalculator.GetResult(trajPoint_vec, 1);
  // case 2:
  //   return fAngleCalculator.GetResult(trajPoint_vec, 0);
  // case 3:
  //   return fAngleCalculator.GetResult(trajPoint_vec, 1);
  } // switch angleMethod
} // MCSMomentumCalculator::SetAngleResultForAngleMethod

// TODO: Add a note in documentation about async calls because of sigmaRES being reset.
// TODO: It would be better if this was a const function.
MCSMomentumCalculator::Result MCSMomentumCalculator::GetResult(const simb::MCTrajectory trajectory, const int angleMethod, size_t startOffset, size_t endOffset) {

  // Will reset after the result has been acquired.
  // This is because sigmaRES should be 0 for Monte Carlo Tracks,
  // though this could be changed for testing purposes.
  double currentSigmaRES = fSigmaRES;
  fSigmaRES = 0;

  // Alternatively, check if it's in the detector instead?
  // A couple of checks
  if(startOffset >= endOffset) {
    if(endOffset >= 0)
      throw "ERROR: startOffset >= endOffset in MCSMomentumCalculator::GetResult(trajectory, ...).";
    else
      endOffset = trajectory.size();
  }

  if(endOffset > trajectory.size()-1)
    std::cout << "WARNING: endOffset is larger than trajectory.size()-1.  Will resort to only looping to trajectory.size()." << std::endl;

  // Set trajPoint_vec
  std::vector<Point_t> trajPoint_vec;
  for(size_t i = startOffset; i < trajectory.size(); ++i) {
    if(i > endOffset)
      break;
    TVector3 position = trajectory.Position(i).Vect();
    Point_t point(position.X(), position.Y(), position.Z());
    point.SetP(trajectory.Momentum(i).P());
    if(inDetector(point))
      trajPoint_vec.push_back(point);
  }
  // Get Result
  MCSMomentumResult result = GetResult(trajPoint_vec, angleMethod);
  // Reset sigmaRES
  fSigmaRES = currentSigmaRES;
  // Return Result
  return result;
}

// Set Functions:

// Transform the currentSegment into the coordinate system where the z-axis lies on the previousSegment
TVector3 MCSMomentumCalculator::transform(TVector3 currentSegment, TVector3 previousSegment) {
  // Define the rotation that will occur as an angle and a vector to rotate around.
  Double_t polarAngle = previousSegment.Theta(); // The angle the currentSegment will be rotating to form vPrime.
  TVector3 rotateAround = TVector3(0,0,1).Cross(previousSegment).Unit(); // The unit vector that the currentSegment will rotate around to form vPrime, Perpendicular to the previousSegment and the z-axis.
  
  // Initialize vPrime and rotate.
  TVector3 vPrime = currentSegment;
  vPrime.Rotate(-polarAngle, rotateAround);
  
  return vPrime;
}

// Calculate the angles of the currentSegment projected onto the previousSegment.
AnglePair_t MCSMomentumCalculator::calculateAnglePair(TVector3 currentSegment, TVector3 previousSegment) {
  // vPrime is now the currentSegment defined in the coordinate system where the previousSegment lies on the z-axis.
  TVector3 vPrime = transform(currentSegment, previousSegment);

  // Calculate projection angles.
  Float_t thetaXZprime = TMath::ATan2(vPrime.X(), vPrime.Z());
  Float_t thetaYZprime = TMath::ATan2(vPrime.Y(), vPrime.Z());
  
  // Return std::pair of projection angles.
  return std::make_pair(thetaXZprime, thetaYZprime);
}

// TODO: Documentation
std::vector<Float_t> MCSMomentumCalculator::SetLogLikelihood_vec(std::vector<AnglePair_t> anglePair_vec, std::vector<Float_t> segmentLength_vec) const {
  std::vector<Float_t> likelihood_vec;
  Float_t llh; //temporary likelihood to be push_back at the end
  for(double p = fMomentumMinimum; p <= fMomentumMaximum; p += fMomentumStep) {
    double mom = p*1000;
    llh = /*0.5**/(anglePair_vec.size())*log(2*TMath::Pi());
    for(size_t i = 0; i < anglePair_vec.size(); i++) {
      if(mom > 0) {
	Float_t thetaXZprime = anglePair_vec.at(i).first, thetaYZprime = anglePair_vec.at(i).second, l = segmentLength_vec.at(i+1);

	double sigma = highland(mom/1000, l, fParticleMass);
	sigma = pow(1*(sigma*sigma + fSigmaRES*fSigmaRES), 0.5);

	// TODO: if(theta <= 3*sigma)
	if(thetaXZprime <= 3*sigma && thetaYZprime <= 3*sigma)
	  llh += 2*log(sigma) + 0.5*pow(thetaXZprime/sigma, 2) + 0.5*pow(thetaYZprime/sigma, 2);
	mom -= deltaP(mom, l, fParticleMass);
      } else {
	break;
      }
    }
    likelihood_vec.push_back(llh);
  } // p loop to get -ln(L) for each p
  return likelihood_vec;
}

// TODO: Documentation
std::vector<Float_t> MCSMomentumCalculator::SetLogLikelihood_vec(std::vector<Float_t> rawAngle_vec, std::vector<Float_t> segmentLength_vec) const {
  std::vector<Float_t> likelihood_vec;
  Float_t llh; //temporary likelihood to be push_back at the end
  for(double p = fMomentumMinimum; p <= fMomentumMaximum; p += fMomentumStep) {
    double mom = p*1000;
    llh = 0.5*(rawAngle_vec.size())*log(2*TMath::Pi());
    for(size_t i = 0; i < rawAngle_vec.size(); i++) {
      if(mom > 0) {
	Float_t rawAngle = rawAngle_vec.at(i), l = segmentLength_vec.at(i+1);

	double sigma = highland(mom/1000,l, fParticleMass);
	sigma = pow(1.72*(sigma*sigma + fSigmaRES*fSigmaRES), 0.5);

	// TODO: if(theta <= 3*sigma)
	if(rawAngle <= 3*sigma)
	  llh += log(sigma) + 0.5*pow(rawAngle/sigma,2);
	mom -= deltaP(mom,l, fParticleMass);
      } else {
	break;
      }
    }
    likelihood_vec.push_back(llh);
  } // p loop to get -ln(L) for each p
  return likelihood_vec;
}

std::vector<Float_t> MCSMomentumCalculator::SetMCSMomentum_vec(std::vector<Float_t> logLikelihood_vec, std::vector<Float_t> segmentLength_vec) const {
    int index = std::min_element(logLikelihood_vec.begin(), logLikelihood_vec.end()) - logLikelihood_vec.begin();

  // Start Momentum
  Float_t momentum = fMomentumMinimum + fMomentumStep*index; // TODO: Move these to function inputs, rather than retrieving from fVariables
  // Start Momentum Correction
  momentum += (deltaP(momentum*1000, segmentLength_vec.at(0), fParticleMass))/1000;

  // Vector of momentum values for each segment.
  std::vector<Float_t> MCSMomentum_vec;
  MCSMomentum_vec.push_back(momentum);

  for(int i = 0; i <= (int) segmentLength_vec.size() - 2; i++) {
    momentum -= (deltaP(momentum*1000, segmentLength_vec.at(i), fParticleMass))/1000;
    MCSMomentum_vec.push_back(momentum);
  }

  // Momentum at the start of each segment.
  return MCSMomentum_vec;
}

// TODO: Documentation
std::vector<AnglePair_t> MCSMomentumCalculator::ExtractProjectedAnglePair_vec(std::vector<Angles_t> angles_vec ) const {
  std::vector<AnglePair_t> anglePair_vec;

  for(size_t i = 0; i < angles_vec.size(); ++i) {
    Angles_t angles = angles_vec.at(i);
    AnglePair_t anglePair = std::make_pair(std::get<1>(angles), std::get<2>(angles));
    anglePair_vec.push_back(anglePair);
  }

  return anglePair_vec;
}

std::vector<Float_t> MCSMomentumCalculator::ExtractTheta3D_vec(std::vector<Angles_t> angles_vec) const {
  std::vector<Float_t> theta3D_vec;
  
  for(size_t i = 0; i < angles_vec.size(); ++i) {
    Angles_t angles = angles_vec.at(i);
    Float_t theta3D = std::get<0>(angles);
    theta3D_vec.push_back(theta3D);
  }

  return theta3D_vec;
}

// TODO: Documentation

// Highland formula. Provide p in GeV/c, returns 
double MCSMomentumCalculator::highland(double p, double l, double particleMass) const {
  double M = particleMass/1000;
  double kappa = (fKappa_a/(p*p)) + fKappa_c;
  double beta = p / sqrt(p*p + M*M);
  double radLength = 14.0;
  double eps = 0.038;

  return 0.001*(kappa/(p*beta))*std::sqrt(l/radLength)*(1 + eps*std::log(l/radLength));
} // p must be given in GeV

double MCSMomentumCalculator::deltaP(double p, double l, double particleMass) const {
  double M = particleMass; // MeV/c^2, mass of incident particle
  double initialE = pow(M*M+p*p,0.5); // MeV
  double Z = 18, A = 39.948; //For Ar, A = [g/mol]
  double K = 0.307075; // MeV*cm^2/mol
  double I = 188*pow(10,-6); //MeV
  double m = 0.511; // MeV/c^2, electron mass
  double bg = p/M, bb = (p*p)/(M*M+p*p), gg = (M*M+p*p)/(M*M); //Beta*Gamma, Beta^2, Gamma^2
  double rho = 1.396; // Liquid Argon density, g/cm^3
  double W = (2*m*bb*gg)/(1+(2*pow(gg,0.5)*m/M)+pow(m/M,2));
  double del, a = 0.19559, x0 = 0.2, x1 = 3, C = 5.2146, k = 3, x = log10(bg);
  if(x >= x1)
    del = 2*log(10)*x-C;
  else if(x >= x0 && x < x1)
    del = 2*log(10)*x-C+a*pow(x1-x,k);
  else
    del = 0;

  //First calculate dEdx (MeV)/cm
  double dEdx = rho*K*(Z/A)*(1/bb)*(0.5*log(2*m*bb*gg*W/(I*I))-bb-del/2);
  ////std::cout << "dEdx = " << dEdx << std::endl;

  return p-pow(pow(initialE - dEdx*l, 2) - M*M,0.5);
} // p must be given in MeV

// Some very crude unit tests, that could be integrated into a better system later.
// Currently called at beginJob() of my analyzer module.
Bool_t MCSMomentumCalculator::test_transform() {
  // Setup 1: No transform
  TVector3 previousSegment1(0,0,1);
  TVector3 currentSegment1(1,2,3); // TODO: This doesn't need to be unit does it?
  TVector3 transformedCurrentSegment1 = transform(currentSegment1, previousSegment1);
  if(currentSegment1 != transformedCurrentSegment1) {
    std::cerr << "test_transform() failed: Current segment is not equal to transform current segment when previous segment is the z-axis." << std::endl;
    return false;
  }
  
  // Setup 2: A 90 degree transform
  TVector3 previousSegment2(0,1,0); // The y-axis
  TVector3 currentSegment2(1,0,0); // The x-axis
  TVector3 transformedCurrentSegment2 = transform(currentSegment2, previousSegment2);
  if(currentSegment2 != transformedCurrentSegment2) {
    std::cerr << "test_transform() failed:: Current segment is not equal to transformed current segment when previousSegment and current segment are both 90 degrees from the z-axis." << std::endl;
    return false;
  }
  
  // Setup 3: A 30 degree transform in the x,y plane
  TVector3 previousSegment3(0,0,1);
  previousSegment3.SetMagThetaPhi(1, 0, TMath::Pi()/6);
  TVector3 currentSegment3(0,0,1);
  currentSegment3.SetMagThetaPhi(1, 0, TMath::Pi()/3);
  TVector3 transformedCurrentSegment3 = transform(currentSegment3, previousSegment3);
  if(transformedCurrentSegment3 != previousSegment3) {
    std::cerr << "test_transform() failed: Transformed Current Segment is not equal to (1,0,Pi/6) (sph. coords) after 30 degree transform." << std::endl;
    return false;
  }

  // Setup 4: A random transform, does the rawAngle between the two segments equal the polar angle of the transformed segment?
  TVector3 previousSegment4(0, 0, 1);
  previousSegment4.SetTheta(fabs(TRandom3().Gaus(0,0.25)*TMath::Pi()));
  previousSegment4.SetPhi(TRandom3().Rndm()*2*TMath::Pi());
  TVector3 currentSegment4(0, 0, 1);
  currentSegment4.SetTheta(fabs(TRandom3().Gaus(0,0.25)*TMath::Pi()));
  currentSegment4.SetPhi(TRandom3().Rndm()*2*TMath::Pi());
  Double_t rawAngle4 = currentSegment4.Angle(previousSegment4);
  TVector3 transformedCurrentSegment4 = transform(currentSegment4, previousSegment4);
  Double_t polarAngle4 = transformedCurrentSegment4.Theta();
  if(rawAngle4 != polarAngle4) {
    std::cerr << "test_transform() failed: Polar angle of transformed segment does not equal raw angle between randomly chosen segments." << std::endl;
    return false;
  }

  // TODO: What other setups can be made?  More random ones would be good.

  return true;
}

Bool_t MCSMomentumCalculator::test_calculateAnglePair() {
  // Setup 1
  TVector3 previousSegment1(0,1,0);
  TVector3 currentSegment1(1,0,0);
  AnglePair_t anglePair1;
  anglePair1 = calculateAnglePair(currentSegment1, previousSegment1);
  if((anglePair1.first - TMath::Pi()/2) > 10e-6) {
    std::cerr << "test_calculateAnglePair() failed in setup 1: thetaXZprime is not equal to pi/2 for given setup 1" << std::endl;
    return false;
  }
  if((anglePair1.second - 0) > 10e-6) {
    std::cerr << "test_calculateAnglePair() failed in setup 1: thetaYZprime is not equal to 0 for given setup 1." << std::endl;
    return false;
  }

  // TODO: What other possible setups could be added?

  return true;
}

void MCSMomentumCalculator::RunTest() {
  Bool_t test1 = test_transform();
  Bool_t test2 = test_calculateAnglePair();

  if(test1 && test2)
    std::cout << "MCSMomentumCalculator: All tests have run successfully." << std::endl;
  else {
    std::cout << "MCSMomentumCalculator: Error in tests." << std::endl;
    throw "MCSMomentumCalculator::RunTest() error";
  }
}

template<typename Point>
Bool_t MCSMomentumCalculator::inDetector(Point point) const {
  return 
    point.Z() <= 700 && point.Z() >= 0 &&
    point.X() <= 300  && point.X() >= -300 &&
    point.Y() <= 600 && point.Y() >= 0;
}
