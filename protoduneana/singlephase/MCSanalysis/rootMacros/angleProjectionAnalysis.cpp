// This file contains a root macro to compare angle projections as a result of transformations.

// Root Includes
#include "TVector3.h"

// Type declarations
using AnglePair = std::pair<Double_t, Double_t>;

// Helper function declarations
TVector3 transform(TVector3 currentSegment, TVector3 previousSegment, Bool_t reverseTransform = false);
AnglePair calculateAnglePair(TVector3 currentSegment, TVector3 previousSegment);
TVector3 GetRandomSegment(Double_t theta_stdev);

// Root Macro that is being run
void angleProjectionAnalysis() {
  // Constant Parameters used in analysis, the same for all setups.
  const Double_t PI = TMath::Pi();
  const Double_t TWOPI = 2*PI;
  const Double_t THETA_STDEV = PI/100; // Roughly 0.0314 rads = 31.4 mrads

  // Declare histograms for each setup
  TH1F* rawAngle1_HIST = new TH1F("rawAngles1_HIST", "#theta_1", 100, -0.15, 0.15);
  TH1F* thetaXZprime1_HIST = new TH1F("thetaXZprime1_HIST","#theta_{XZ1}'", 100, -0.15, 0.15);
  TH1F* thetaYZprime1_HIST = new TH1F("thetaYZprime1_HIST","#theta_{YZ1}'", 100, -0.15, 0.15);

  TH1F* rawAngle4_HIST = new TH1F("rawAngles4_HIST", "#theta_4", 100, -0.15, 0.15);
  TH1F* thetaXZprime4_HIST = new TH1F("thetaXZprime4_HIST","#theta_{XZ4}'", 100, -0.15, 0.15);
  TH1F* thetaYZprime4_HIST = new TH1F("thetaYZprime4_HIST","#theta_{YZ4}'", 100, -0.15, 0.15);

  // Modify histogram attributes for each setup
  // Setup 1
  thetaXZprime1_HIST->SetLineColor(kRed);
  thetaYZprime1_HIST->SetLineColor(kGreen+1);

  thetaXZprime1_HIST->SetFillStyle(0);
  thetaYZprime1_HIST->SetFillStyle(0);

  // Setup 4
  thetaXZprime4_HIST->SetLineColor(kRed);
  thetaYZprime4_HIST->SetLineColor(kGreen+1);

  thetaXZprime4_HIST->SetFillStyle(0);
  thetaYZprime4_HIST->SetFillStyle(0);
  
  gRandom->SetSeed(); // Set a random new seed.
  // Fill histograms for various setups
  for(int i = 0; i < 10e6; ++i) {

    // Calculate Data for each setup

    // Setup 1: PreviousSegment = z-axis, currentSegment = random.
    // This should follow the expected distribution.
    TVector3 previousSegment1(0,0,1);
    TVector3 currentSegment1 = GetRandomSegment(THETA_STDEV);
    Double_t rawAngle1 = currentSegment1.Angle(previousSegment1); // Should be the same currentSegment1's polar angle, since previousSegment = z-axis
    AnglePair anglePair1 = calculateAnglePair(currentSegment1, previousSegment1); // No transform should take place.
    Double_t thetaXZprime1 = anglePair1.first;
    Double_t thetaYZprime1 = anglePair1.second;
    
    // Setup 4: previousSegment = random, currentSegment = random but with respect to the first.  Current segment is made by making a random segment, then doing a reverse transform into the previousSegment's reference frame.
    // This should follow the expected distribution.
    TVector3 previousSegment4 = GetRandomSegment(THETA_STDEV);
    TVector3 randomSegment = GetRandomSegment(THETA_STDEV);
    TVector3 currentSegment4 = transform(randomSegment, previousSegment4, true); // Reverse Transformation
    Double_t rawAngle4 = currentSegment4.Angle(previousSegment4);
    AnglePair anglePair4 = calculateAnglePair(currentSegment4, previousSegment4);
    Double_t thetaXZprime4 = anglePair4.first;
    Double_t thetaYZprime4 = anglePair4.second;

    // TODO: Assertion that these angles are the same as the angles of randomSegment with respect to z-axis?

    // Fill Histograms for each setup
    // Setup 1
    rawAngle1_HIST->Fill(rawAngle1);
    thetaXZprime1_HIST->Fill(thetaXZprime1);
    thetaYZprime1_HIST->Fill(thetaYZprime1);

    // Setup 4
    rawAngle4_HIST->Fill(rawAngle4);
    thetaXZprime4_HIST->Fill(thetaXZprime4);
    thetaYZprime4_HIST->Fill(thetaYZprime4);
  }

  // Calculate Standard deviation & error in standard deviation of plots for setups 1 & 4
  Double_t rawAngle1_stdev = rawAngle1_HIST->GetStdDev();
  Double_t thetaXZprime1_stdev = thetaXZprime1_HIST->GetStdDev();
  Double_t thetaYZprime1_stdev = thetaYZprime1_HIST->GetStdDev();
  Double_t rawAngle1_stdeverr = rawAngle1_HIST->GetStdDevError();
  Double_t thetaXZprime1_stdeverr = thetaXZprime1_HIST->GetStdDevError();
  Double_t thetaYZprime1_stdeverr = thetaYZprime1_HIST->GetStdDevError();

  Double_t rawAngle4_stdev = rawAngle4_HIST->GetStdDev();
  Double_t thetaXZprime4_stdev = thetaXZprime4_HIST->GetStdDev();
  Double_t thetaYZprime4_stdev = thetaYZprime4_HIST->GetStdDev();
  Double_t rawAngle4_stdeverr = rawAngle4_HIST->GetStdDevError();
  Double_t thetaXZprime4_stdeverr = thetaXZprime4_HIST->GetStdDevError();
  Double_t thetaYZprime4_stdeverr = thetaYZprime4_HIST->GetStdDevError();

  // Print Expectation
  std::cout << std::endl;
  std::cout << "theta_stdev used in segment generation = " << THETA_STDEV << std::endl;
  std::cout << "ratio = measured_stdev / theta_stdev_generation" << std::endl;
  // Print Setup 1
  std::cout << "Setups 1 & 4: (should be correct results)" << std::endl;
  std::cout << "rawAngle1_stdev = " << rawAngle1_stdev << " +/- " << rawAngle1_stdeverr << std::endl;
  std::cout << "thetaXZprime1_stdev = " << thetaXZprime1_stdev << " +/- " << thetaXZprime1_stdeverr << std::endl;
  std::cout << "thetaYZprime1_stdev = " << thetaYZprime1_stdev << " +/- " << thetaYZprime1_stdeverr << std::endl;
  std::cout << "rawAngle1_ratio = " << rawAngle1_stdev / THETA_STDEV << std::endl;
  std::cout << "thetaXZprime1_ratio = " << thetaXZprime1_stdev / THETA_STDEV << std::endl;
  std::cout << "thetaYZprime1_ratio = " << thetaYZprime1_stdev / THETA_STDEV << std::endl;
  std::cout << "setup 1: rawSigma / XZsigma = " << rawAngle1_stdev / thetaXZprime1_stdev << std::endl;
  std::cout << "setup 1: rawSigma / YZsigma = " << rawAngle1_stdev / thetaYZprime1_stdev<< std::endl;
  std::cout << std::endl;
  // Print Setup 4
  std::cout << "rawAngle4_stdev = " << rawAngle4_stdev << " +/- " << rawAngle4_stdeverr << std::endl;
  std::cout << "thetaXZprime4_stdev = " << thetaXZprime4_stdev << " +/- " << thetaXZprime4_stdeverr << std::endl;
  std::cout << "thetaYZprime4_stdev = " << thetaYZprime4_stdev << " +/- " << thetaYZprime4_stdeverr << std::endl;
  std::cout << "rawAngle4_ratio = " << rawAngle4_stdev / THETA_STDEV << std::endl;
  std::cout << "thetaXZprime4_ratio = " << thetaXZprime4_stdev / THETA_STDEV << std::endl;
  std::cout << "thetaYZprime4_ratio = " << thetaYZprime4_stdev / THETA_STDEV << std::endl;
  std::cout << "setup 4: rawSigma / XZsigma = " << rawAngle4_stdev / thetaXZprime4_stdev  << std::endl;
  std::cout << "setup 4: rawSigma / YZsigma = " << rawAngle4_stdev / thetaYZprime4_stdev << std::endl;
  std::cout << std::endl;

  // Construct histogram stack
  THStack* histogramStack14 = new THStack("hs14","Distribution of angles, Setups 1 & 4; #theta (mrad); Counts/bin");
  histogramStack14->Add(rawAngle1_HIST);
  histogramStack14->Add(rawAngle4_HIST);
  histogramStack14->Add(thetaXZprime1_HIST);
  histogramStack14->Add(thetaXZprime4_HIST);
  histogramStack14->Add(thetaYZprime1_HIST);
  histogramStack14->Add(thetaYZprime4_HIST);

  // Draw histogram stack
  histogramStack14->Draw("nostack");
  gPad->BuildLegend(0.2, 0.6, 0.4, 0.8);
}

// Helper function definitions, taken directly from larreco/RecoAlg/MCSMomentumCalculator.cxx
// Transform the currentSegment into the coordinate system where the z-axis lies on the previousSegment
// larreco/RecoAlg/MCSMomentumCalculator::transform doesn't include the reverseTransform parameter.
TVector3 transform(TVector3 currentSegment, TVector3 previousSegment, Bool_t reverseTransform) {
  // Define the rotation that will occur as an angle and a vector to rotate around.
  Double_t polarAngle = previousSegment.Theta(); // The angle the currentSegment will be rotating to form vPrime.
  TVector3 rotateAround = TVector3(0,0,1).Cross(previousSegment).Unit(); // The unit vector that the currentSegment will rotate around to form vPrime, Perpendicular to the previousSegment and the z-axis.
  
  // Initialize vPrime and rotate.
  TVector3 vPrime = currentSegment;
  vPrime.Rotate(reverseTransform ? polarAngle : -polarAngle, rotateAround); // It used to just be -polarAngle
  
  return vPrime;
}

// Calculate the angles of the currentSegment projected onto the previousSegment.
AnglePair calculateAnglePair(TVector3 currentSegment, TVector3 previousSegment) {
  // vPrime is now the currentSegment defined in the coordinate system where the previousSegment lies on the z-axis.
  TVector3 vPrime = transform(currentSegment, previousSegment);

  // Calculate projection angles.
  Double_t thetaXZprime = TMath::ATan2(vPrime.X(), vPrime.Z());
  Double_t thetaYZprime = TMath::ATan2(vPrime.Y(), vPrime.Z());
  
  // Return std::pair of projection angles.
  return std::make_pair(thetaXZprime, thetaYZprime);
}

// GetRandomSegment, not in larreco/RecoAlg/MCSMomentumCalculator.cxx
TVector3 GetRandomSegment(Double_t theta_stdev) {
  const Double_t PI = TMath::Pi();
  const Double_t TWOPI = 2*PI;

  TVector3 segment(0,0,1);
  Double_t segmentTheta = gRandom->Gaus(0,theta_stdev);
  Double_t segmentPhi = gRandom->Rndm()*TWOPI;
  segment.SetTheta(segmentTheta);
  segment.SetPhi(segmentPhi);

  return segment;
}
