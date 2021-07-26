////////////////////////////////////////////////////////////////////////
// Class:       MCSAngleAnalysis
// Module Type: analyzer
// File:        MCSAngleAnalysis_module.cc
// Author: Hunter Meyer | hmeyer5@lsu.edu
// Advisor: Thomas Kutter
// Description: Analyze MCS 3D and projected angles
//
// Measure modified-highland formula fit parameters.
// Measure detector angular resolution.
// Plot 3D angles, projected angles, and compare distribution
// parameters with nominal values.
////////////////////////////////////////////////////////////////////////

#define _unused [[maybe_unused]]

// C++ Includes
#include <iostream>

// Root Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TVector3.h"

// Framework Includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

// LArSoft Inclues
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "MCSAngleCalculator.h"
#include "MCSSegmentCalculator.h"

// LSU Includes
#include "BackTrackerAlg.h"

class MCSAngleAnalysis;

class MCSAngleAnalysis : public art::EDAnalyzer {
public:

  // fhicl configuration
  struct Config {
    fhicl::Atom<art::InputTag> trackModuleLabel { fhicl::Name("TrackModuleLabel") };
    fhicl::Atom<Bool_t> shouldCreateVirtualPoints { fhicl::Name("ShouldCreateVirtualPoints"), fhicl::Comment("Determines if virtual points will be created for MCS segments"), false };

    fhicl::Table<trkf::MCSAngleCalculator::Config> angleCalculator { fhicl::Name("angleCalculator") };
  };
  // using Point_t = recob::tracking::Point_t;
  using Point_t = trkf::MCSSegmentCalculator::Point_t;
  using Vector_t = recob::tracking::Vector_t;
  using Segment_t = trkf::MCSSegmentCalculator::Segment_t;
  using Angles_t = trkf::MCSAngleCalculator::Angles_t;

  // Constructor
  explicit MCSAngleAnalysis(art::EDAnalyzer::Table<Config> const & t);

  // Plugins should not be copied or assigned.
  MCSAngleAnalysis(MCSAngleAnalysis const &) = delete;
  MCSAngleAnalysis(MCSAngleAnalysis &&) = delete;
  MCSAngleAnalysis & operator = (MCSAngleAnalysis const &) = delete;
  MCSAngleAnalysis & operator = (MCSAngleAnalysis &&) = delete;

  // Art Functions.
  void beginJob();
  void analyze(art::Event const & e) override;
  void endJob();

  // Local Functions
  Double_t highland(Double_t momentum, Double_t l) const;
  double deltaP(double p, double l) const;
  template<typename Point> Bool_t inDetector(Point point) const;
  TVector3 startingDirectionFor(simb::MCParticle particle) const;
  TVector3 startingPositionFor(simb::MCParticle particle) const;
  template<typename Vector1, typename Vector2> Double_t distanceBetween(Vector1 vector1, Vector2 vector2) const;
  template<typename Vector> size_t TrajectoryPointWithClosestZ(simb::MCParticle particle, Vector recoPoint) const;
  //std::vector<double> getThetaOverSigma(double p, std::vector<double> thetaXZprimevec, std::vector<double> thetaYZprimevec, std::vector<double> lengthsvec);

private:
  
  // Retrieved from fhicl parameters
  const art::InputTag fTrackModuleLabel;
  const Bool_t fShouldCreateVirtualPoints;
  //const trkf::MCSAngleCalculator::Parameters fMCSAngleCalculatorParameters;
  const fhicl::Table<trkf::MCSAngleCalculator::Config> fAngleCalculatorParameters;

  // TODO: Add reco angles counts
  // Number of angles within each sigma range.
  Float_t trueLinearXZ_oneSigmaCount = 0;
  Float_t trueLinearXZ_twoSigmaCount = 0;
  Float_t trueLinearXZ_threeSigmaCount = 0;
  Float_t trueLinearXZ_outsideThreeSigmaCount = 0;

  Float_t trueLinearYZ_oneSigmaCount = 0;
  Float_t trueLinearYZ_twoSigmaCount = 0;
  Float_t trueLinearYZ_threeSigmaCount = 0;
  Float_t trueLinearYZ_outsideThreeSigmaCount = 0;

  Float_t trueLinear3D_oneSigmaCount = 0;
  Float_t trueLinear3D_twoSigmaCount = 0;
  Float_t trueLinear3D_threeSigmaCount = 0;
  Float_t trueLinear3D_outsideThreeSigmaCount = 0;

  Float_t totalAngleCount = 0;

  int counter = 0;

  // True Data
  // =====================================================
  // True Linear 3D Angles
  TH2F* trueLinear_theta3DVsSegmentMomentum_HIST;
  // True Polygonal 3D Angles
  TH2F* truePolygonal_theta3DVsSegmentMomentum_HIST;
  // Linear Angles
  TH2F* trueLinear_thetaXZprimeVsSegmentMomentum_HIST;
  TH2F* trueLinear_thetaYZprimeVsSegmentMomentum_HIST;
  // Polygonal Angles
  TH2F* truePolygonal_thetaXZprimeVsSegmentMomentum_HIST;
  TH2F* truePolygonal_thetaYZprimeVsSegmentMomentum_HIST;
  // Compare 3D vs. 2D
  TH2F* trueLinear_theta3DVsTheta2DInQuadrature_HIST;
  TH2F* truePolygonal_theta3DVsTheta2DInQuadrature_HIST;
  // =====================================================

  // Reco Data
  // =====================================================
  // Reco Linear 3D Angles
  TH2F* recoLinear_theta3DVsSegmentMomentum_HIST;
  // Reco Polygonal 3D Angles
  TH2F* recoPolygonal_theta3DVsSegmentMomentum_HIST;
  // Linear Angles
  TH2F* recoLinear_thetaXZprimeVsSegmentMomentum_HIST;
  TH2F* recoLinear_thetaYZprimeVsSegmentMomentum_HIST;
  // Polygonal Angles
  TH2F* recoPolygonal_thetaXZprimeVsSegmentMomentum_HIST;
  TH2F* recoPolygonal_thetaYZprimeVsSegmentMomentum_HIST;
  // Compare 3D vs. 2D
  TH2F* recoLinear_theta3DVsTheta2DInQuadrature_HIST;
  TH2F* recoPolygonal_theta3DVsTheta2DInQuadrature_HIST;
  // =====================================================

  // Compare True and Reco
  // =====================================================
  TH2F* recoLinearThetaXZprimeVsTrueLinearThetaXZprime_HIST;
  TH2F* recoLinearThetaYZprimeVsTrueLinearThetaYZprime_HIST;
  // Polygonal
  // =====================================================

  // First Segment Analysis
  // =====================================================
  // First Reco Segment
  TH2F* first_recoDots_theta3DVsSegmentMomentum_HIST;
  TH2F* first_recoDots_thetaXZprimeVsSegmentMomentum_HIST;
  TH2F* first_recoDots_thetaYZprimeVsSegmentMomentum_HIST;

  TH2F* first_recoLinear_theta3DVsSegmentMomentum_HIST;
  TH2F* first_recoLinear_thetaXZprimeVsSegmentMomentum_HIST;
  TH2F* first_recoLinear_thetaYZprimeVsSegmentMomentum_HIST;
  
  TH2F* first_recoLinear_theta3DVsSegmentBBMomentum_HIST;
  TH2F* first_recoLinear_thetaXZprimeVsSegmentBBMomentum_HIST;
  TH2F* first_recoLinear_thetaYZprimeVsSegmentBBMomentum_HIST;

  // First True Segment
  TH2F* first_trueLastPoint_theta3DVsSegmentMomentum_HIST;
  TH2F* first_trueLastPoint_thetaXZprimeVsSegmentMomentum_HIST;
  TH2F* first_trueLastPoint_thetaYZprimeVsSegmentMomentum_HIST;

  TH2F* first_trueDots_theta3DVsSegmentMomentum_HIST;
  TH2F* first_trueDots_thetaXZprimeVsSegmentMomentum_HIST;
  TH2F* first_trueDots_thetaYZprimeVsSegmentMomentum_HIST;

  TH2F* first_trueLinear_theta3DVsSegmentMomentum_HIST;
  TH2F* first_trueLinear_thetaXZprimeVsSegmentMomentum_HIST;
  TH2F* first_trueLinear_thetaYZprimeVsSegmentMomentum_HIST;

  TH2F* first_trueLinear_theta3DVsSegmentBBMomentum_HIST;
  TH2F* first_trueLinear_thetaXZprimeVsSegmentBBMomentum_HIST;
  TH2F* first_trueLinear_thetaYZprimeVsSegmentBBMomentum_HIST;

  // TODO: Finish implementing alternative segment definitions in MCSSegmentCalculator and put MCSSegmentCalculator into MCSAngleCalculator. (refactor)

  // Distance between starting position of reco & true track.
  TH1F* reco_vertexDifference_HIST;

  // TODO: First polygonal
  // TODO: First theta vs. momentum AND vs. length (TH3F...)
  // =====================================================

  // =====================================================
  // Bethe-Bloch Verification
  TH1F* distanceBetweenSegmentStartPoints_HIST;
  TH1F* distanceBetweenRecoPointWithTruePointWithClosestZ_HIST;
  TH1F* distanceBetweenTruePointWithTruePointWithClosestZ_HIST;

  TH2F* distanceBetweenSegmentStartPointsVsSegmentNumber_HIST;
  TH2F* distanceBetweenRecoPointWithTruePointWithClosestZVsSegmentNumber_HIST;
  TH2F* distanceBetweenTruePointWithTruePointWithClosestZVsSegmentNumber_HIST;

  TH2F* true_rawSegmentLengthVsStraightLineDistance_HIST;
  TH2F* reco_rawSegmentLengthVsStraightLineDistance_HIST;

  TH2F* true_rawSegmentLengthVsSegmentMomentum_HIST;
  TH2F* reco_rawSegmentLengthVsSegmentMomentumClosestZ_HIST;

  TH2F* recoBBMomentumVsTrueMomentumClosestZ_HIST; // Using closestZ method for true and BB for reco
  TH2F* trueBBMomentumVsTrueMomentum_HIST; // Using BB with true trajectory points method
  // =====================================================
};

MCSAngleAnalysis::MCSAngleAnalysis(art::EDAnalyzer::Table<MCSAngleAnalysis::Config> const & t) : 
  art::EDAnalyzer(t),				     
  fTrackModuleLabel(t().trackModuleLabel()),
  fShouldCreateVirtualPoints(t().shouldCreateVirtualPoints()),
  fAngleCalculatorParameters(t().angleCalculator)
{}

// TODO: Verify that curvy tracks were no longer occuring, then show that they are now fixed...
// TODO: Plot the theta/sigma plots as well. (derived in the trkf::MCSMomentumCalculator::Result?)

void MCSAngleAnalysis::analyze(art::Event const & e) {
  std::cout << "Event number: " << e.event() << std::endl;

  // MCS SegmentCalculator
  trkf::MCSSegmentCalculator segmentCalculator = trkf::MCSSegmentCalculator(fAngleCalculatorParameters().segmentCalculator);
  // MCS Angle Calculator
  trkf::MCSAngleCalculator angleCalculator = trkf::MCSAngleCalculator(fAngleCalculatorParameters);
  // LSU Backtracker
  lsu::BackTrackerAlg backtracker;

  // TODO: if(!e.isRealData())
  auto particleListHandle = e.getValidHandle<std::vector<simb::MCParticle> >("largeant");
  std::vector<art::Ptr<simb::MCParticle> > particleList;
  art::fill_ptr_vector(particleList, particleListHandle);

  // Loop through MCParticles
  size_t nParticles = particleList.size();
  for(size_t iParticle = 0; iParticle < nParticles; ++iParticle) {
    art::Ptr<simb::MCParticle> particle = particleList.at(iParticle);
    if(particle->Trajectory().TotalLength() >= 100 && particle->PdgCode() == 13) {
      // Calculate MCS Angles
      // TODO: Update to just calculate the MCSSegmentResult here and pass to GetResult methods once those are added to prevent having to run the segmentation algorithms multiple times.
      // Get an MCSSegmentResult
      trkf::MCSSegmentCalculator::MCSSegmentResult trueSegmentResult = segmentCalculator.GetResult(*particle, fShouldCreateVirtualPoints);

      // Extract data from the true segment result.
      // std::vector<size_t> trueSegmentStartPoint_vec = trueSegmentResult.GetSegmentStartPoint_vec();
      std::vector<Float_t> trueRawSegmentLength_vec = trueSegmentResult.GetRawSegmentLength_vec();

      // Get an MCSAngleResult for various methods using the true segment result.
      trkf::MCSAngleCalculator::MCSAngleResult trueLinearMCS = angleCalculator.GetResult(trueSegmentResult, 0);
      trkf::MCSAngleCalculator::MCSAngleResult truePolygonalMCS = angleCalculator.GetResult(trueSegmentResult, 1);

      // Extract the data from the various angle results.
      std::vector<Angles_t> trueLinearAngles_vec = trueLinearMCS.GetAngles_vec();
      std::vector<Angles_t> truePolygonalAngles_vec = truePolygonalMCS.GetAngles_vec();

      // Loop through angle measurements for each segment.
      for(size_t i = 0; i < trueLinearAngles_vec.size() /* && i < trueSegmentStartPoint_vec.size() - 1 */; ++i) {
	std::cout << "i = " << i << std::endl;
	// The i+1 is here since our first angle is for the second segment.
	// _unused size_t segmentStartPoint = trueSegmentStartPoint_vec.at(i+1);
	// Momentum for segment.
	// Double_t segmentMomentum = particle->P(segmentStartPoint);
	Double_t segmentMomentum = trueSegmentResult.GetMomentumForSegment(i+1);
	Double_t segmentLength = trueRawSegmentLength_vec.at(i+1);

	// TODO: This momentum is only actual valid for the start of linear segments, but it's pretty close for the polygonal segments.
	//     Consider alternative: Bethe-Bloch energy loss for polyognal segments?
	// TODO: Use average segment momentum?

	// Angles for segment.
	Angles_t linearAngles = trueLinearAngles_vec.at(i);
	Angles_t polygonalAngles = truePolygonalAngles_vec.at(i);

	Float_t linearTheta3D = std::get<0>(linearAngles)*1000.0;
	Float_t linearThetaXZprime = std::get<1>(linearAngles)*1000.0;
	Float_t linearThetaYZprime = std::get<2>(linearAngles)*1000.0;
	Float_t polygonalTheta3D = std::get<0>(polygonalAngles)*1000.0;
	Float_t polygonalThetaXZprime = std::get<1>(polygonalAngles)*1000.0;
	Float_t polygonalThetaYZprime = std::get<2>(polygonalAngles)*1000.0;

	trueLinear_theta3DVsTheta2DInQuadrature_HIST->Fill(pow(linearThetaXZprime*linearThetaXZprime + linearThetaYZprime*linearThetaYZprime,0.5), linearTheta3D);
	truePolygonal_theta3DVsTheta2DInQuadrature_HIST->Fill(pow(polygonalThetaXZprime*polygonalThetaXZprime + polygonalThetaYZprime*polygonalThetaYZprime,0.5), polygonalTheta3D);

	// Fill histograms with projected angles vs. segment momentum.
	trueLinear_theta3DVsSegmentMomentum_HIST->Fill(segmentMomentum, linearTheta3D);
	trueLinear_thetaXZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, linearThetaXZprime);
	trueLinear_thetaYZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, linearThetaYZprime);
	truePolygonal_theta3DVsSegmentMomentum_HIST->Fill(segmentMomentum, polygonalTheta3D);
	truePolygonal_thetaXZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, polygonalThetaXZprime);
	truePolygonal_thetaYZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, polygonalThetaYZprime);

	// Update angles within sigma counts
	Double_t sigmaHL = 1000.0*highland(segmentMomentum, segmentLength);
	// LinearThetaXZprime
	if(fabs(linearThetaXZprime) <= sigmaHL)
	  ++trueLinearXZ_oneSigmaCount;
	if(fabs(linearThetaXZprime) <= 2*sigmaHL)
	  ++trueLinearXZ_twoSigmaCount;
	if(fabs(linearThetaXZprime) <= 3*sigmaHL)
	  ++trueLinearXZ_threeSigmaCount;
	else
	  ++trueLinearXZ_outsideThreeSigmaCount;

	// LinearThetaYZprime
	if(fabs(linearThetaYZprime) <= sigmaHL)
	  ++trueLinearYZ_oneSigmaCount;
	if(fabs(linearThetaYZprime) <= 2*sigmaHL)
	  ++trueLinearYZ_twoSigmaCount;
	if(fabs(linearThetaYZprime) <= 3*sigmaHL)
	  ++trueLinearYZ_threeSigmaCount;
	else
	  ++trueLinearYZ_outsideThreeSigmaCount;

	// Linear 3D angle
	if(fabs(linearTheta3D) <= sigmaHL)
	  ++trueLinear3D_oneSigmaCount;
	if(fabs(linearTheta3D) <= 2*sigmaHL)
	  ++trueLinear3D_twoSigmaCount;
	if(fabs(linearTheta3D) <= 3*sigmaHL)
	  ++trueLinear3D_threeSigmaCount;
	else
	  ++trueLinear3D_outsideThreeSigmaCount;

	// Total Angle Count
	++totalAngleCount;

      } // for i (segment angle)
    } // if TotalLength >= 100 cm
  } // for iParticle

  // recob::Track analysis
  auto trackListHandle = e.getValidHandle<std::vector<recob::Track> >(fTrackModuleLabel);
  std::vector<art::Ptr<recob::Track> > trackList;
  art::fill_ptr_vector(trackList, trackListHandle);
  std::cout << "Number of Tracks: " << trackList.size() << std::endl;

  // Loop through the list of reconstructed tracks
  size_t nTracks = trackList.size();
  for(size_t iTrack = 0; iTrack < nTracks; iTrack++) {

    // Get this specific track
    art::Ptr<recob::Track> track = trackList.at(iTrack);

    // (attempt to) Backtrack this track to an MCParticle
    simb::MCParticle particle = backtracker.getMCParticle(track, e, fTrackModuleLabel.label());

    // Only run this analysis if the track length is greater than 100 cm
    if(track->Length() >= 100) {
      // Calculate Reco MCSSegmentResult
      trkf::MCSSegmentCalculator::MCSSegmentResult recoSegmentResult = segmentCalculator.GetResult(*track, fShouldCreateVirtualPoints);

      // Extract data from the reco segment result
      std::vector<Float_t> rawSegmentLength_vec = recoSegmentResult.GetRawSegmentLength_vec();
      // std::vector<size_t> reco_segmentStartPoint_vec = recoSegmentResult.GetSegmentStartPoint_vec();
      std::vector<Segment_t> reco_segment_vec = recoSegmentResult.GetSegment_vec();
      std::vector<Vector_t> reco_linearFit_vec = recoSegmentResult.GetLinearFit_vec();

      // Calculate Reco MCSAngleResult for various methods
      trkf::MCSAngleCalculator::MCSAngleResult recoLinearMCS = angleCalculator.GetResult(recoSegmentResult, 0);
      trkf::MCSAngleCalculator::MCSAngleResult recoPolygonalMCS = angleCalculator.GetResult(recoSegmentResult, 1);

      // Extract data from the various reco angle results
      std::vector<Angles_t> recoLinearAngles_vec = recoLinearMCS.GetAngles_vec();
      std::vector<Angles_t> recoPolygonalAngles_vec = recoPolygonalMCS.GetAngles_vec();
      
      // Calculate the True MCSSegmentResult
      trkf::MCSSegmentCalculator::MCSSegmentResult trueSegmentResult = segmentCalculator.GetResult(particle, fShouldCreateVirtualPoints);

      // Extract data from the true segment result
      // std::vector<size_t> true_segmentStartPoint_vec = trueSegmentResult.GetSegmentStartPoint_vec();
      std::vector<Segment_t> true_segment_vec = trueSegmentResult.GetSegment_vec();
      std::vector<Float_t> true_rawSegmentLength_vec = trueSegmentResult.GetRawSegmentLength_vec();
      std::vector<Vector_t> true_linearFit_vec = trueSegmentResult.GetLinearFit_vec();

      // Calculate an True MCSAngleResult for various methods
      trkf::MCSAngleCalculator::MCSAngleResult trueLinearMCS = angleCalculator.GetResult(particle, 0);
      // TODO: truePolygonalMCS

      // Extract data from the various true angle results
      std::vector<Angles_t> trueLinearAngles_vec = trueLinearMCS.GetAngles_vec();

      // Loop through angle measurements for each segment.
      for(size_t i = 0; i < recoLinearAngles_vec.size(); ++i) {

	// Get the segment momentum of this reco segment using the "closest Z" method.
	// size_t recoSegmentStartPoint = reco_segmentStartPoint_vec.at(i);
	Segment_t recoSegment = reco_segment_vec.at(i+1);
	Point_t recoPoint = recoSegment.at(0);
	// TVector3 recoPoint = angleCalculator.convert<Point_t, TVector3>(track->LocationAtPoint(recoSegmentStartPoint));
	size_t truePointIndex = TrajectoryPointWithClosestZ(particle, recoPoint);
	Double_t segmentMomentum = particle.P(truePointIndex);

	// Get the collection of angles for this segment.
	Angles_t recoLinearAngles = recoLinearAngles_vec.at(i);
	Angles_t recoPolygonalAngles = recoPolygonalAngles_vec.at(i);

	// Extract the individual angles for this segment.
	Float_t linearTheta3D = std::get<0>(recoLinearAngles)*1000.0;
	Float_t linearThetaXZprime = std::get<1>(recoLinearAngles)*1000.0;
	Float_t linearThetaYZprime = std::get<2>(recoLinearAngles)*1000.0;
	Float_t polygonalTheta3D = std::get<0>(recoPolygonalAngles)*1000.0;
	Float_t polygonalThetaXZprime = std::get<1>(recoPolygonalAngles)*1000.0;
	Float_t polygonalThetaYZprime = std::get<2>(recoPolygonalAngles)*1000.0;

	// Fill hyistograms comparing 3D and 2D added in quadrature. (1-to-1) graph.
	recoLinear_theta3DVsTheta2DInQuadrature_HIST->Fill(pow(linearThetaXZprime*linearThetaXZprime + linearThetaYZprime*linearThetaYZprime,0.5), linearTheta3D);
	recoPolygonal_theta3DVsTheta2DInQuadrature_HIST->Fill(pow(polygonalThetaXZprime*polygonalThetaXZprime + polygonalThetaYZprime*polygonalThetaYZprime,0.5), polygonalTheta3D);

	// Fill histograms with projected angles vs. segment momentum.
	recoLinear_theta3DVsSegmentMomentum_HIST->Fill(segmentMomentum, linearTheta3D);
	recoLinear_thetaXZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, linearThetaXZprime);
	recoLinear_thetaYZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, linearThetaYZprime);
	recoPolygonal_theta3DVsSegmentMomentum_HIST->Fill(segmentMomentum, polygonalTheta3D);
	recoPolygonal_thetaXZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, polygonalThetaXZprime);
	recoPolygonal_thetaYZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, polygonalThetaYZprime);

      } // for i (segment angle)

      // Calculate and plot the distance between vertices the true and reconstructed vertex.
      TVector3 initialPosition = startingPositionFor(particle);
      recob::tracking::Point_t reco_vertex = track->Vertex();
      Double_t vertexDifference = distanceBetween(initialPosition, reco_vertex);
      reco_vertexDifference_HIST->Fill(vertexDifference);

      // Direct Comparison Analysis and First Segment Analysis
      // Only run for events with the same number of segments.
      if(vertexDifference < 1.0) {
	// ========================================================================
	// Direct Comparison Analysis across many segments
	// Only possible if they're the same size (same # of segments)
	// ========================================================================
	if(trueLinearAngles_vec.size() == recoLinearAngles_vec.size()) {
	  // Loop through segments
	  for(size_t i = 0; i < recoLinearAngles_vec.size(); ++i) {

	    Angles_t trueLinearAngles = trueLinearAngles_vec.at(i);
	    Angles_t recoLinearAngles = recoLinearAngles_vec.at(i);
	    // TODO: Polygonal and/or other segment definitions direct comparison.
	    
	    // True/Reco angles for this segment
	    Double_t trueLinearThetaXZprime = std::get<1>(trueLinearAngles)*1000.0;
	    Double_t trueLinearThetaYZprime = std::get<2>(trueLinearAngles)*1000.0;
	    Double_t recoLinearThetaXZprime = std::get<1>(recoLinearAngles)*1000.0;
	    Double_t recoLinearThetaYZprime = std::get<2>(recoLinearAngles)*1000.0;

	    // Fill true/reco comparison histograms.
	    recoLinearThetaXZprimeVsTrueLinearThetaXZprime_HIST->Fill(trueLinearThetaXZprime, recoLinearThetaXZprime);
	    recoLinearThetaYZprimeVsTrueLinearThetaYZprime_HIST->Fill(trueLinearThetaYZprime, recoLinearThetaYZprime);
	  } // for i in recoLinearAngle_vec.size()
	} // If the number of true segments == the number of reco segments

	// ========================================================================
	// First Segment Analysis
	// ========================================================================

	// Initial Info (prior to first segment being formed)
	Double_t initialMomentum = particle.P();
	TVector3 initialDirection = startingDirectionFor(particle);
	
	// Get initial vector of trueMCS / recoMCS first segment.

	// TODO: This is a test for the MCSSegemntCalculator.convert methods.
	// trkf::MCSSegmentCalculator segmentCalculator(14.0);

	// First Segment Vectors
	// TVector3 trueLastPoint_firstSegment = particle.Momentum(true_segmentStartPoint_vec.at(1)).Vect().Unit();
	// TODO: Alternative segment defintions
	TVector3 trueLinear_firstSegment = angleCalculator.convert<Vector_t, TVector3>(true_linearFit_vec.at(0));
	TVector3 recoLinear_firstSegment = angleCalculator.convert<Vector_t, TVector3>(reco_linearFit_vec.at(0));

	// Calculate first segment angles
	// Angles_t trueLastPoint_firstSegmentAngles = angleCalculator.calculateAngles(trueLastPoint_firstSegment, initialDirection);
	Angles_t trueLinear_firstSegmentAngles = angleCalculator.calculateAngles(angleCalculator.convert<Vector_t, TVector3>(true_linearFit_vec.at(1)), trueLinear_firstSegment);
	Angles_t recoLinear_firstSegmentAngles = angleCalculator.calculateAngles(recoLinear_firstSegment, initialDirection);
	
	// First True Segment Angles, individually
	// _unused Float_t trueLastPoint_theta3D = std::get<0>(trueLastPoint_firstSegmentAngles)*1000.0;
	// _unused Float_t trueLastPoint_thetaXZprime = std::get<1>(trueLastPoint_firstSegmentAngles)*1000.0;
	// _unused Float_t trueLastPoint_thetaYZprime = std::get<2>(trueLastPoint_firstSegmentAngles)*1000.0;

	_unused Float_t trueLinear_theta3D = std::get<0>(trueLinear_firstSegmentAngles)*1000.0;
	_unused Float_t trueLinear_thetaXZprime = std::get<1>(trueLinear_firstSegmentAngles)*1000.0;
	_unused Float_t trueLinear_thetaYZprime = std::get<2>(trueLinear_firstSegmentAngles)*1000.0;
	
	// TODO: true dots segments

	// First Reco Segment Angles, individually
	_unused Float_t recoLinear_theta3D = std::get<0>(recoLinear_firstSegmentAngles)*1000.0;
	_unused Float_t recoLinear_thetaXZprime = std::get<1>(recoLinear_firstSegmentAngles)*1000.0;
	_unused Float_t recoLinear_thetaYZprime = std::get<2>(recoLinear_firstSegmentAngles)*1000.0;

	// TODO: reco dots segments
	
	// Fill True Histograms
	// first_trueLastPoint_theta3DVsSegmentMomentum_HIST->Fill(initialMomentum, trueLastPoint_theta3D);
	// first_trueLastPoint_thetaXZprimeVsSegmentMomentum_HIST->Fill(initialMomentum, trueLastPoint_thetaXZprime);
	// first_trueLastPoint_thetaYZprimeVsSegmentMomentum_HIST->Fill(initialMomentum, trueLastPoint_thetaYZprime);

	// TODO: Alternative Segment Definitions
	// Fill True Histograms
	// True Momentum
	first_trueLinear_theta3DVsSegmentMomentum_HIST->Fill(initialMomentum, trueLinear_theta3D);
	first_trueLinear_thetaXZprimeVsSegmentMomentum_HIST->Fill(initialMomentum, trueLinear_thetaXZprime);
	first_trueLinear_thetaYZprimeVsSegmentMomentum_HIST->Fill(initialMomentum, trueLinear_thetaYZprime);

	// BB Momentum
	first_trueLinear_theta3DVsSegmentBBMomentum_HIST->Fill(initialMomentum, trueLinear_theta3D);
	first_trueLinear_thetaXZprimeVsSegmentBBMomentum_HIST->Fill(initialMomentum, trueLinear_thetaXZprime);
	first_trueLinear_thetaYZprimeVsSegmentBBMomentum_HIST->Fill(initialMomentum, trueLinear_thetaYZprime);

	// Fill Reco Histograms
	// True Momentum (closest Z method techcnially for later segments)
	first_recoLinear_theta3DVsSegmentMomentum_HIST->Fill(initialMomentum, recoLinear_theta3D);
	first_recoLinear_thetaXZprimeVsSegmentMomentum_HIST->Fill(initialMomentum, recoLinear_thetaXZprime);
	first_recoLinear_thetaYZprimeVsSegmentMomentum_HIST->Fill(initialMomentum, recoLinear_thetaYZprime);

	// BB Momentum
	first_recoLinear_theta3DVsSegmentBBMomentum_HIST->Fill(initialMomentum, recoLinear_theta3D);
	first_recoLinear_thetaXZprimeVsSegmentBBMomentum_HIST->Fill(initialMomentum, recoLinear_thetaXZprime);
	first_recoLinear_thetaYZprimeVsSegmentBBMomentum_HIST->Fill(initialMomentum, recoLinear_thetaYZprime);

	// Optionally, add to the first-segment analysis plots using the first few segments.
	// The momentum is not as well-defined, which is the point of the Bethe-Bloch verification that is also being done (code is after this section).
	// I'll update this documentation if we can conclude that the momentum IS well-defined and we're using the updated/correct method.
	// Ideally, we can use a reco-formula rather than the MC info, but we might need it instead (if we use the closest Z method, for instance).

	// Loop through the next few segments and add to the first segment plots.
	
	// We'll begin by defining the number of segments to be all of the segments.
	size_t numberOfTrueSegments = trueLinearAngles_vec.size();
	size_t numberOfRecoSegments = recoLinearAngles_vec.size();

	// Optionally reset the number of segments to some value to just look at the first couple of segments.
	// numberOfTrueSegments = 4;
	// numberOfRecoSegments = 4;

	// Calculate the inital true momentum
	// _unused Double_t initialTrueSegmentStartPoint = true_segmentStartPoint_vec.at(0);
	Double_t initialTrueMomentum = trueSegmentResult.GetMomentumForSegment(0);
	// Double_t initialTrueMomentum = particle.P(initialTrueSegmentStartPoint);

	// The segment momentum at the start of the first segment.
	Double_t trueSegmentBBMomentum = initialTrueMomentum;

	// Plot the first few true segments and perform BB verification.
	for(size_t i = 0; i < true_rawSegmentLength_vec.size() /* && i < true_segmentStartPoint_vec.size() - 1 */; ++i) { // TODO: Tempororary i < true_segmentStartPoint_vec.size()

	  std::cout << "i true = " << i << std::endl;

	  // Update the true segment BB momentum for this segment.
	  Double_t trueRawSegmentLength = true_rawSegmentLength_vec.at(i);
	  trueSegmentBBMomentum -= deltaP(trueSegmentBBMomentum*1000.0, trueRawSegmentLength)/1000.0;

	  // Calculate the true segment momentum for the start of this segment.
	  // _unused size_t trueSegmentStartPoint = true_segmentStartPoint_vec.at(i+1);
	  Double_t trueSegmentMomentum = trueSegmentResult.GetMomentumForSegment(i+1);
	  // Double_t trueSegmentMomentum = particle.P(trueSegmentStartPoint);

	  // Fill BB momentum verification histograms
	  trueBBMomentumVsTrueMomentum_HIST->Fill(trueSegmentMomentum, trueSegmentBBMomentum);

	  // Compare the straight-line distance between adjacent segment start points with this segment's raw segment length:
	  // size_t trueSegmentStartPoint_beginning = true_segmentStartPoint_vec.at(i);
	  Segment_t trueSegment = true_segment_vec.at(i);
	  Point_t truePoint_beginning = trueSegment.front();
	  Point_t truePoint_end = trueSegment.back();
	  // TVector3 truePoint_beginning = particle.Position(trueSegmentStartPoint_beginning).Vect();
	  // TVector3 truePoint_end = particle.Position(trueSegmentStartPoint).Vect();
	  Double_t straightLineDistance = distanceBetween(truePoint_beginning, truePoint_end);

	  // Fill segment length histogram
	  true_rawSegmentLengthVsStraightLineDistance_HIST->Fill(straightLineDistance, trueRawSegmentLength);
	  true_rawSegmentLengthVsSegmentMomentum_HIST->Fill(trueSegmentMomentum, trueRawSegmentLength);

	  // If we should add the first few segments to the first segment angles plots.
	  // numberOfTrueSegments will be 0 if not.
	  // TODO: Re-add the Bool?
	  if(i > 0 && i < numberOfTrueSegments) {

	    // Get the True Angles
	    // The i-1 is because of the relative sizes of the vectors and because we have the i > 0 condition.
	    Angles_t trueLinearAngles = trueLinearAngles_vec.at(i-1);

	    // Extract the True Angles
	    Float_t trueLinear_theta3D = std::get<0>(trueLinearAngles)*1000.0;
	    Float_t trueLinear_thetaXZprime = std::get<1>(trueLinearAngles)*1000.0;
	    Float_t trueLinear_thetaYZprime = std::get<2>(trueLinearAngles)*1000.0;

	    // Fill True Angle Histograms for/using trueSegmentMomentum
	    first_trueLinear_theta3DVsSegmentMomentum_HIST->Fill(trueSegmentMomentum, trueLinear_theta3D);
	    first_trueLinear_thetaXZprimeVsSegmentMomentum_HIST->Fill(trueSegmentMomentum, trueLinear_thetaXZprime);
	    first_trueLinear_thetaYZprimeVsSegmentMomentum_HIST->Fill(trueSegmentMomentum, trueLinear_thetaYZprime);

	    // Fill the True Angle Histograms for/using trueSegmentBBMomentum
	    first_trueLinear_theta3DVsSegmentBBMomentum_HIST->Fill(trueSegmentBBMomentum, trueLinear_theta3D);
	    first_trueLinear_thetaXZprimeVsSegmentBBMomentum_HIST->Fill(trueSegmentBBMomentum, trueLinear_thetaXZprime);
	    first_trueLinear_thetaYZprimeVsSegmentBBMomentum_HIST->Fill(trueSegmentBBMomentum, trueLinear_thetaYZprime);
	  } // if i > 0 && i < numberOfTrueSegments
	} // for i in true_rawSegmentLength_vec.size()

	// Define the initial recoSegmentBBMomentum
	Double_t recoSegmentBBMomentum = initialTrueMomentum;

	// Plot the first few reco segments.
	for(size_t i = 0; i < rawSegmentLength_vec.size() /* && i < reco_segmentStartPoint_vec.size() - 1 */; ++i) { // TODO: Tempororary i < reco_segmentStartPoint_vec.size()

	  std::cout << "i reco = " << i << std::endl;

	  // Calculate the reco segment momentum.
	  // This is currently using the TrajectoryPointWithClosestZ method, TODO: we should quantify the changes this has compared with using BB.
	  // TODO: Overlay sigmaHL/sigmaRMS vs. segment momentum using BB and closest Z momentum measurement.  We should see right-shifted sigma values, right?   Why do we still see sigmaHL > sigmaRMS at low p right now?  What happens if we get sigmaHL using BB momentum rather than particle.P(trueTrajectoryPoint)?  Solve on Friday.
	  // TODO: Combine Bethe-Bloch Verification / Analysis with this for loop, use those values to additionally plot: first_recoLinear_theta...SegmentBBMomentum_HIST->Fill(recoSegmentBBMomentum, recoLinear_theta...);

	  // Calculate the true segment momentum for the true point that is closest in Z to the reco segment start point
	  // size_t recoSegmentStartPoint = reco_segmentStartPoint_vec.at(i+1);
	  // Point_t recoPoint_t = track->LocationAtPoint(recoSegmentStartPoint);
	  // TVector3 recoPoint = TVector3(recoPoint_t.X(), recoPoint_t.Y(), recoPoint_t.Z());
	  Segment_t recoSegment = reco_segment_vec.at(i+1);
	  Point_t recoPoint = recoSegment.at(0);
	  _unused size_t truePointClosestZ_index = TrajectoryPointWithClosestZ(particle, recoPoint);
	  Double_t trueSegmentMomentumClosestZ = particle.P(truePointClosestZ_index);

	  // Update the reco segment BB momentum for the next segment.
	  Double_t recoRawSegmentLength = rawSegmentLength_vec.at(i);
	  recoSegmentBBMomentum -= deltaP(recoSegmentBBMomentum*1000.0, recoRawSegmentLength)/1000.0;

	  // Fill BB histograms
	  recoBBMomentumVsTrueMomentumClosestZ_HIST->Fill(trueSegmentMomentumClosestZ, recoSegmentBBMomentum);	  

	  // Compare the straight-line distance between adjacent segment start points with this segment's raw segment length:
	  Point_t recoPoint_beginning = recoSegment.front();
	  Point_t recoPoint_end = recoSegment.back();

	  // size_t recoSegmentStartPoint_beginning = reco_segmentStartPoint_vec.at(i);
	  // Point_t recoPoint_beginning = track->LocationAtPoint(recoSegmentStartPoint_beginning);
	  // Point_t recoPoint_end = track->LocationAtPoint(recoSegmentStartPoint);
	  Double_t straightLineDistance = distanceBetween(recoPoint_beginning, recoPoint_end);
	  // TVector3 recoPoint_beginning = angleCalculator.convert<Point_t, TVector3>(track->LocationAtPoint(recoSegmentStartPoint_beginning));
	  // TVector3 recoPoint_end = angleCalculator.convert<Point_t, TVector3>(track->LocationAtPoint(recoSegmentStartPoint));
	  // Point_t recoPoint_beginning_t = track->LocationAtPoint(recoSegmentStartPoint_beginning);
	  // Point_t recoPoint_end_t = track->LocationAtPoint(recoSegmentStartPoint);
	  // TVector3 recoPoint_beginning = TVector3(recoPoint_beginning_t.X(), recoPoint_beginning_t.Y(), recoPoint_beginning_t.Z());
	  // TVector3 recoPoint_end = TVector3(recoPoint_end_t.X(), recoPoint_end_t.Y(), recoPoint_end_t.Z());
	  // Double_t straightLineDistance = distanceBetween(recoPoint_beginning, recoPoint_end);

	  // Fill segment length histogram
	  reco_rawSegmentLengthVsStraightLineDistance_HIST->Fill(straightLineDistance, recoRawSegmentLength);
	  reco_rawSegmentLengthVsSegmentMomentumClosestZ_HIST->Fill(trueSegmentMomentumClosestZ, recoRawSegmentLength);

	  if(i > 0 && i < numberOfRecoSegments) {
	    // Get the Reco Angles
	    Angles_t recoLinearAngles = recoLinearAngles_vec.at(i-1);

	    // Extract the Reco Angles
	    Float_t recoLinear_theta3D = std::get<0>(recoLinearAngles)*1000.0;
	    Float_t recoLinear_thetaXZprime = std::get<1>(recoLinearAngles)*1000.0;
	    Float_t recoLinear_thetaYZprime = std::get<2>(recoLinearAngles)*1000.0;

	    // Fill the Reco Angle Histograms for/using trueSegmentMomentumClosestZ
	    first_recoLinear_theta3DVsSegmentMomentum_HIST->Fill(trueSegmentMomentumClosestZ, recoLinear_theta3D);
	    first_recoLinear_thetaXZprimeVsSegmentMomentum_HIST->Fill(trueSegmentMomentumClosestZ, recoLinear_thetaXZprime);
	    first_recoLinear_thetaYZprimeVsSegmentMomentum_HIST->Fill(trueSegmentMomentumClosestZ, recoLinear_thetaYZprime);

	    // Fill the Reco Angle Histograms for/using recoSegmentBBMomentum
	    first_recoLinear_theta3DVsSegmentBBMomentum_HIST->Fill(recoSegmentBBMomentum, recoLinear_theta3D);
	    first_recoLinear_thetaXZprimeVsSegmentBBMomentum_HIST->Fill(recoSegmentBBMomentum, recoLinear_thetaXZprime);
	    first_recoLinear_thetaYZprimeVsSegmentBBMomentum_HIST->Fill(recoSegmentBBMomentum, recoLinear_thetaYZprime);
	  } // if i > 0  && i < numberOfRecoSegments
	} // for i in numberOfRecoSegments

	// ========================================================================
	// Bethe-Bloch Verification
	// Some of these plots are a little more abstract (think segment analysis, not angle analysis or momentum analysis).
	// If we do the MCSSegmentCalcualtor, some of this can actually move to an analysis module there).
	// MCSSegmentAnalysis_module.cc, as opposed to MCSAngleAnalysis_module.cc or CompareMCSAnalyzer_module.cc (which should be renamed or moved to MCSMomentumAnalysis_module.cc).
	// ========================================================================
	
	// Segment Analysis:
	// TODO: Plot length of segment calculated using straight-line distance between first and last point, compare with raw segment length vec
	// Not technically Bethe-Bloch Verification, but it will tell us something about the segments that we can use to better interpret the BB verification.
	// Calculate the distance between segment start points for true and reco segment start points for tracks that have the same number of segments.
	// for(size_t i = 0; i < true_segmentStartPoint_vec.size() && i < reco_segmentStartPoint_vec.size(); ++i) {
	for(size_t i = 0; i < true_segment_vec.size() && i < reco_segment_vec.size(); ++i) {
	  std::cout << "i = " << i << std::endl;
	  // Define this segment.
	  Segment_t trueSegment = true_segment_vec.at(i);
	  std::cout << "in between" << std::endl;
	  Segment_t recoSegment = reco_segment_vec.at(i);

	  std::cout << "HERE" << std::endl;

	  // Define various points.
	  Point_t recoPoint = recoSegment.front();
	  // size_t recoSegmentStartPoint = reco_segmentStartPoint_vec.at(i);
	  // Point_t recoPoint = track->LocationAtPoint(recoSegmentStartPoint);
	  // Point_t recoPoint_t = track->LocationAtPoint(recoSegmentStartPoint);
	  // TVector3 recoPoint = TVector3(recoPoint_t.X(), recoPoint_t.Y(), recoPoint_t.Z());
	  
	  // size_t trueSegmentStartPoint = true_segmentStartPoint_vec.at(i);
	  // TVector3 truePoint = particle.Position(trueSegmentStartPoint).Vect();
	  Point_t truePoint = trueSegment.front();

	  size_t truePointWithClosestZToRecoPoint_index = TrajectoryPointWithClosestZ(particle, recoPoint);
	  std::cout << "Here?" << std::endl;
	  TVector3 truePointWithClosestZToRecoPoint = particle.Position(truePointWithClosestZToRecoPoint_index).Vect();

	  // Plot distances between various points
	  Double_t distanceBetweenSegmentStartPoints = distanceBetween(recoPoint, truePoint);
	  Double_t distanceBetweenRecoPointWithTruePointWithClosestZ = distanceBetween(recoPoint, truePointWithClosestZToRecoPoint);
	  Double_t distanceBetweenTruePointWithTruePointWithClosestZ = distanceBetween(truePoint, truePointWithClosestZToRecoPoint);

	  // Fill Histograms
	  distanceBetweenSegmentStartPoints_HIST->Fill(distanceBetweenSegmentStartPoints);
	  distanceBetweenRecoPointWithTruePointWithClosestZ_HIST->Fill(distanceBetweenRecoPointWithTruePointWithClosestZ);
	  distanceBetweenTruePointWithTruePointWithClosestZ_HIST->Fill(distanceBetweenTruePointWithTruePointWithClosestZ);

	  distanceBetweenSegmentStartPointsVsSegmentNumber_HIST->Fill(i, distanceBetweenSegmentStartPoints);
	  distanceBetweenRecoPointWithTruePointWithClosestZVsSegmentNumber_HIST->Fill(i, distanceBetweenRecoPointWithTruePointWithClosestZ);
	  distanceBetweenTruePointWithTruePointWithClosestZVsSegmentNumber_HIST->Fill(i, distanceBetweenTruePointWithTruePointWithClosestZ);

	  std::cout << "I'm confused" << std::endl;

	} // for i in trueSegmentStartPoints.size() && recoSegmentStartPoints.size()

	// TODO: Plot reco BB momentum vs. true momentum of true point with closest Z to the reco point
	// TODO: Plot true BB momentum vs. true momentum of true point

      } // if nTracks == 1
      else {
	counter++;
      }
    } // if length > 100
  } // for iTrack
} // analyze function

void MCSAngleAnalysis::beginJob() {
  art::ServiceHandle<art::TFileService> tfs; // Consider moving to private class scope, rather than function scope, if we do endJob analysis that needs to be saved (such as sigmaRES analysis).

  // TODO: You *must* expand the angle range today.

  // True Data
  // ========================================================================
  // True Linear 3D Angles
  trueLinear_theta3DVsSegmentMomentum_HIST = tfs->make<TH2F>("trueLinear_theta3DVsSegmentMomentum_HIST","True Linear #theta_{3D}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta'",500, -0.05, 5.05, 211, -10.5, 200.5);
  // True Polygonal 3D Angles
  truePolygonal_theta3DVsSegmentMomentum_HIST = tfs->make<TH2F>("truePolygonal_theta3DVsSegmentMomentum_HIST","True Polygonal #theta_{3D}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta'", 500, -0.05, 5.05, 211, -10.5, 200.5);
  // True Linear Angles
  trueLinear_thetaXZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("trueLinear_thetaXZprimeVsSegmentMomentum_HIST","True Linear #theta_{XZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{XZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  trueLinear_thetaYZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("trueLinear_thetaYZprimeVsSegmentMomentum_HIST","True Linear #theta_{YZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{YZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  // True Polygonal Angles
  truePolygonal_thetaXZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("truePolygonal_thetaXZprimeVsSegmentMomentum_HIST","True Polygonal #theta_{XZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{XZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  truePolygonal_thetaYZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("truePolygonal_thetaYZprimeVsSegmentMomentum_HIST","True Polygonal #theta_{YZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{YZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  // Compare 3D and 2D angles
  trueLinear_theta3DVsTheta2DInQuadrature_HIST = tfs->make<TH2F>("trueLinear_theta3DVsTheta2DInQuadrature_HIST", "True Linear #theta_{3D}' Vs. 2D Angles in Quad; 2D Angles Added in Quadrature; #theta_{3D}'", 500, 0, 500, 500, 0, 500);
  truePolygonal_theta3DVsTheta2DInQuadrature_HIST = tfs->make<TH2F>("truePolygonal_theta3DVsTheta2DInQuadrature_HIST", "True Polyognal #theta_{3D}' Vs. 2D Angles in Quad; 2D Angles Added in Quadrature; #theta_{3D}'", 500, 0, 500, 500, 0, 500);
  // ========================================================================

  // Reco Data
  // ========================================================================
  // Reco Linear 3D Angles
  recoLinear_theta3DVsSegmentMomentum_HIST = tfs->make<TH2F>("recoLinear_theta3DVsSegmentMomentum_HIST","Reco Linear #theta_{3D}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta'",500, -0.05, 5.05, 211, -10.5, 200.5);
  // Reco Polygonal 3D Angles
  recoPolygonal_theta3DVsSegmentMomentum_HIST = tfs->make<TH2F>("recoPolygonal_theta3DVsSegmentMomentum_HIST","Reco Polygonal #theta_{3D}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta'", 500, -0.05, 5.05, 211, -10.5, 200.5);
  // Reco Linear Angles
  recoLinear_thetaXZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("recoLinear_thetaXZprimeVsSegmentMomentum_HIST","Reco Linear #theta_{XZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{XZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  recoLinear_thetaYZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("recoLinear_thetaYZprimeVsSegmentMomentum_HIST","Reco Linear #theta_{YZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{YZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  // Reco Polygonal Angles
  recoPolygonal_thetaXZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("recoPolygonal_thetaXZprimeVsSegmentMomentum_HIST","Reco Polygonal #theta_{XZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{XZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  recoPolygonal_thetaYZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("recoPolygonal_thetaYZprimeVsSegmentMomentum_HIST","Reco Polygonal #theta_{YZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{YZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  // Compare 3D and 2D angles
  recoLinear_theta3DVsTheta2DInQuadrature_HIST = tfs->make<TH2F>("recoLinear_theta3DVsTheta2DInQuadrature_HIST", "Reco Linear #theta_{3D}' Vs. 2D Angles in Quad; 2D Angles Added in Quadrature; #theta_{3D}'", 500, 0, 500, 500, 0, 500);
  recoPolygonal_theta3DVsTheta2DInQuadrature_HIST = tfs->make<TH2F>("recoPolygonal_theta3DVsTheta2DInQuadrature_HIST", "Reco Polyognal #theta_{3D}' Vs. 2D Angles in Quad; 2D Angles Added in Quadrature; #theta_{3D}'", 500, 0, 500, 500, 0, 500);
  // ========================================================================

  // Compare True and Reco
  // ========================================================================
  recoLinearThetaXZprimeVsTrueLinearThetaXZprime_HIST = tfs->make<TH2F>("recoLinearThetaXZprimeVsTrueLinearThetaXZprime_HIST", "Reco Linear #theta_{XZ}' Vs. True Linear #theta_{XZ}'; True #theta_{XZ}'; Reco #theta_{XZ}'", 401, -200.5, 200.5, 401, -200.5, 200.5);
  recoLinearThetaYZprimeVsTrueLinearThetaYZprime_HIST = tfs->make<TH2F>("recoLinearThetaYZprimeVsTrueLinearThetaYZprime_HIST", "Reco Linear #theta_{YZ}' Vs. True Linear #theta_{YZ}'; True #theta_{YZ}'; Reco #theta_{YZ}'", 401, -200.5, 200.5, 401, -200.5, 200.5);

  // ========================================================================

  // First Segment Analysis
  // ========================================================================
  // TODO: After we've tested the relationship between recoDots & recoLinear angles, we can check the implementation of MCSMomentumCalculator using the dots segments
  // TODO: Rename dots segments to something a little more descriptive but also short.  border segments?
  first_recoDots_theta3DVsSegmentMomentum_HIST = tfs->make<TH2F>("first_recoDots_theta3DVsSegmentMomentum_HIST","First Segment Reco Dots #theta_{3D}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{3D}'",500, -0.05, 5.05, 211, -10.5, 200.5);
  first_recoDots_thetaXZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("first_recoDots_thetaXZprimeVsSegmentMomentum_HIST", "First Segment Reco Dots #theta_{XZ}' Vs. Segment Momentum; Segment Momentum (GeV/c); #theta_{XZ}'", 500, -0.05, 5.05, 401, -200.5, 200.5);
  first_recoDots_thetaYZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("first_recoDots_thetaYZprimeVsSegmentMomentum_HIST", "First Segment Reco Dots #theta_{YZ}' Vs. Segment Momentum; Segment Momentum (GeV/c); #theta_{YZ}'", 500, -0.05, 5.05, 401, -200.5, 200.5);

  // First Reco Linear Theta vs. Segment Momentum
  first_recoLinear_theta3DVsSegmentMomentum_HIST = tfs->make<TH2F>("first_recoLinear_theta3DVsSegmentMomentum_HIST","First Segment Reco Linear #theta_{3D}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{3D}'",500, -0.05, 5.05, 211, -10.5, 200.5);
  first_recoLinear_thetaXZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("first_recoLinear_thetaXZprimeVsSegmentMomentum_HIST","First Segment Reco Linear #theta_{XZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{XZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  first_recoLinear_thetaYZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("first_recoLinear_thetaYZprimeVsSegmentMomentum_HIST","First Segment Reco Linear #theta_{YZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{YZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);

  // First Reco Linear Theta vs. Segment BB Momentum
  first_recoLinear_theta3DVsSegmentBBMomentum_HIST = tfs->make<TH2F>("first_recoLinear_theta3DVsSegmentBBMomentum_HIST","First Segment Reco Linear #theta_{3D}' Vs Segment BB Momentum; Segment BB Momentum (GeV/c); #theta_{3D}'",500, -0.05, 5.05, 211, -10.5, 200.5);
  first_recoLinear_thetaXZprimeVsSegmentBBMomentum_HIST = tfs->make<TH2F>("first_recoLinear_thetaXZprimeVsSegmentBBMomentum_HIST","First Segment Reco Linear #theta_{XZ}' Vs Segment BB Momentum; Segment BB Momentum (GeV/c); #theta_{XZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  first_recoLinear_thetaYZprimeVsSegmentBBMomentum_HIST = tfs->make<TH2F>("first_recoLinear_thetaYZprimeVsSegmentBBMomentum_HIST","First Segment Reco Linear #theta_{YZ}' Vs Segment BB Momentum; Segment BB Momentum (GeV/c); #theta_{YZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);

  // True First Segments
  first_trueLastPoint_theta3DVsSegmentMomentum_HIST = tfs->make<TH2F>("first_trueLastPoint_theta3DVsSegmentMomentum_HIST","First Last Point #theta_{3D}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{3D}' (mrad)",500, -0.05, 5.05, 211, -10.5, 200.5);
  first_trueLastPoint_thetaXZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("first_trueLastPoint_thetaXZprimeVsSegmentMomentum_HIST","First Segment True Last Point #theta_{XZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{XZ}' (mrad)",500, -0.05, 5.05, 401, -200.5, 200.5);
  first_trueLastPoint_thetaYZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("first_trueLastPoint_thetaYZprimeVsSegmentMomentum_HIST","First Segment True Last Point #theta_{YZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{YZ}' (mrad)",500, -0.05, 5.05, 401, -200.5, 200.5);

  first_trueDots_theta3DVsSegmentMomentum_HIST = tfs->make<TH2F>("first_trueDots_theta3DVsSegmentMomentum_HIST","First Segment True Dots #theta_{3D}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{3D}' (mrad)",500, -0.05, 5.05, 211, -10.5, 200.5);
  first_trueDots_thetaXZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("first_trueDots_thetaXZprimeVsSegmentMomentum_HIST", "First Segment True Dots #theta_{XZ}' Vs. Segment Momentum; Segment Momentum (GeV/c); #theta_{XZ}'", 500, -0.05, 5.05, 401, -200.5, 200.5);
  first_trueDots_thetaYZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("first_trueDots_thetaYZprimeVsSegmentMomentum_HIST", "First Segment True Dots #theta_{YZ}' Vs. Segment Momentum; Segment Momentum (GeV/c); #theta_{YZ}'", 500, -0.05, 5.05, 401, -200.5, 200.5);

  // First True Linear Theta vs. Segment Momentum
  first_trueLinear_theta3DVsSegmentMomentum_HIST = tfs->make<TH2F>("first_trueLinear_theta3DVsSegmentMomentum_HIST","First Segment True Linear #theta_{3D}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{3D}' (mrad)",500, -0.05, 5.05, 211, -10.5, 200.5);
  first_trueLinear_thetaXZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("first_trueLinear_thetaXZprimeVsSegmentMomentum_HIST","First Segment True Linear #theta_{XZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{XZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  first_trueLinear_thetaYZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("first_trueLinear_thetaYZprimeVsSegmentMomentum_HIST","First Segment True Linear #theta_{YZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{YZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);

  // First True Linear Theta vs. Segment BB Momentum
  first_trueLinear_theta3DVsSegmentBBMomentum_HIST = tfs->make<TH2F>("first_trueLinear_theta3DVsSegmentBBMomentum_HIST","First Segment True Linear #theta_{3D}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{3D}' (mrad)",500, -0.05, 5.05, 211, -10.5, 200.5);
  first_trueLinear_thetaXZprimeVsSegmentBBMomentum_HIST = tfs->make<TH2F>("first_trueLinear_thetaXZprimeVsSegmentBBMomentum_HIST","First Segment True Linear #theta_{XZ}' Vs Segment BB Momentum; Segment BB Momentum (GeV/c); #theta_{XZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  first_trueLinear_thetaYZprimeVsSegmentBBMomentum_HIST = tfs->make<TH2F>("first_trueLinear_thetaYZprimeVsSegmentBBMomentum_HIST","First Segment True Linear #theta_{YZ}' Vs Segment BB Momentum; Segment BB Momentum (GeV/c); #theta_{YZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);

  reco_vertexDifference_HIST = tfs->make<TH1F>("reco_vertexDifference_HIST", "Reco Vertex Difference; Distance between vertices (cm); Counts / bin", 500, 0, 50);
  // ========================================================================

  // ========================================================================
  // Bethe-Bloch Verification / Segment Analysis
  distanceBetweenSegmentStartPoints_HIST = tfs->make<TH1F>("distanceBetweenSegmentStartPoints_HIST", "Distance Between Segment Start Points; Distance (cm); Counts / bin", 100, 0, 14);
  distanceBetweenRecoPointWithTruePointWithClosestZ_HIST = tfs->make<TH1F>("distanceBetweenRecoPointWithTruePointWithClosestZ_HIST", "Distance Between Reco Point with True Point with Closest Z; Distance (cm); Counts / bin", 100, 0, 14);
  distanceBetweenTruePointWithTruePointWithClosestZ_HIST = tfs->make<TH1F>("distanceBetweenTruePointWithTruePointWithClosestZ_HIST", "Distance Between True Point with True Point with Closest Z; Distance (cm); Counts / bin", 100, 0, 14);

  distanceBetweenSegmentStartPointsVsSegmentNumber_HIST = tfs->make<TH2F>("distanceBetweenSegmentStartPointsVsSegmentNumber_HIST", "Distance Between Segment Start Points vs. Segment Number; Segment Number; Distance (cm)", 51, -0.5, 50.5, 100, 0, 20);
  distanceBetweenRecoPointWithTruePointWithClosestZVsSegmentNumber_HIST = tfs->make<TH2F>("distanceBetweenRecoPointWithTruePointWithClosestZVsSegmentNumber_HIST", "Distance Between Reco Point with True Point with Closest Z vs. Segment Number; Segment Number; Distance (cm)", 51, -0.5, 50.5, 100, 0, 20);
  distanceBetweenTruePointWithTruePointWithClosestZVsSegmentNumber_HIST = tfs->make<TH2F>("distanceBetweenTruePointWithTruePointWithClosestZVsSegmentNumber_HIST", "Distance Between True Point with True Point with Closest Z vs. Segment Number; Segment Number; Distance (cm)", 51, -0.5, 50.5, 100, 0, 20);

  true_rawSegmentLengthVsStraightLineDistance_HIST = tfs->make<TH2F>("true_rawSegmentLengthVsStraightLineDistance_HIST", "True Raw Segment Length vs. Straight-Line Distance; Straight-Line Distance (cm); Raw Segment Length", 281, -0.5, 28.5, 281, -0.5, 28.5);
  reco_rawSegmentLengthVsStraightLineDistance_HIST = tfs->make<TH2F>("reco_rawSegmentLengthVsStraightLineDistance_HIST", "Reco Raw Segment Length vs. Straight-Line Distance; Straight-Line Distance (cm); Raw Segment Length", 281, -0.5, 28.5, 281, -0.5, 28.5);

  true_rawSegmentLengthVsSegmentMomentum_HIST = tfs->make<TH2F>("true_rawSegmentLengthVsSegmentMomentum_HIST", "True Raw Segment Length vs. Segment Momentum; Segment Momentum (GeV/c); Raw Segment Length (cm)", 451, -0.05, 4.55, 281, -0.5, 28.5);
  reco_rawSegmentLengthVsSegmentMomentumClosestZ_HIST = tfs->make<TH2F>("reco_rawSegmentLengthVsSegmentMomentumClosestZ_HIST", "Reco Raw Segment Length vs. Segment Momentum (closest z method); Segment Momentum (GeV/c); Raw Segment Length (cm)", 451, -0.05, 4.55, 281, -0.5, 28.5);

  recoBBMomentumVsTrueMomentumClosestZ_HIST = tfs->make<TH2F>("recoBBMomentumVsTrueMomentumClosestZ_HIST", "Reco Seg. Momentum vs. True Seg. Momentum (closest z method); True Momentum using Closest Z (GeV/c); Reco Momentum using BB (GeV/c)", 451, -0.05, 4.55, 451, -0.05, 4.55);
  trueBBMomentumVsTrueMomentum_HIST = tfs->make<TH2F>("trueBBMomentumVsTrueMomentum_HIST", "True BB Seg. Momentum vs. True Seg. Momentum; True Seg. Momentum (GeV/c); True BB Seg. Momentum (GeV/c)", 451, -0.05, 4.55, 451, -0.05, 4.55);
  // ========================================================================
}

void MCSAngleAnalysis::endJob() {

  std::cout << "counter = " << counter << std::endl;

  std::cout << std::endl;
  std::cout << "totalAngleCount: " << totalAngleCount << std::endl
	    << std::endl;

  std::cout << std::fixed << std::setprecision(3)
	    << "trueLinearXZ counts:" << std::endl
	    << "within one sigma = " << trueLinearXZ_oneSigmaCount << "\t" << 100*(trueLinearXZ_oneSigmaCount / totalAngleCount) << "%" << std::endl
	    << "within two sigma = " << trueLinearXZ_twoSigmaCount << "\t" << 100*(trueLinearXZ_twoSigmaCount / totalAngleCount) << "%"<< std::endl
	    << "within three sigma = " << trueLinearXZ_threeSigmaCount << "\t" << 100*(trueLinearXZ_threeSigmaCount / totalAngleCount) << "%"<< std::endl
	    << "outside three sigma = " << trueLinearXZ_outsideThreeSigmaCount << "\t" << 100*(trueLinearXZ_outsideThreeSigmaCount / totalAngleCount) << "%"<< std::endl
	    << std::endl;

  std::cout << std::fixed << std::setprecision(3)
	    << "trueLinearYZ counts:" << std::endl
	    << "within one sigma = " << trueLinearYZ_oneSigmaCount << "\t" << 100*(trueLinearYZ_oneSigmaCount / totalAngleCount) << "%" << std::endl
	    << "within two sigma = " << trueLinearYZ_twoSigmaCount << "\t" << 100*(trueLinearYZ_twoSigmaCount / totalAngleCount) << "%"<< std::endl
	    << "within three sigma = " << trueLinearYZ_threeSigmaCount << "\t" << 100*(trueLinearYZ_threeSigmaCount / totalAngleCount) << "%"<< std::endl
	    << "outside three sigma = " << trueLinearYZ_outsideThreeSigmaCount << "\t" << 100*(trueLinearYZ_outsideThreeSigmaCount / totalAngleCount) << "%"<< std::endl
	    << std::endl;

  std::cout << std::fixed << std::setprecision(3)
	    << "trueLinear3D counts:" << std::endl
	    << "within one sigma = " << trueLinear3D_oneSigmaCount << "\t" << 100*(trueLinear3D_oneSigmaCount / totalAngleCount) << "%" << std::endl
	    << "within two sigma = " << trueLinear3D_twoSigmaCount << "\t" << 100*(trueLinear3D_twoSigmaCount / totalAngleCount) << "%"<< std::endl
	    << "within three sigma = " << trueLinear3D_threeSigmaCount << "\t" << 100*(trueLinear3D_threeSigmaCount / totalAngleCount) << "%"<< std::endl
	    << "outside three sigma = " << trueLinear3D_outsideThreeSigmaCount << "\t" << 100*(trueLinear3D_outsideThreeSigmaCount / totalAngleCount) << "%"<< std::endl
	    << std::endl;
}

template<typename Point>
Bool_t MCSAngleAnalysis::inDetector(Point point) const {
  return 
    point.Z() <= 700 && point.Z() >= 0 &&
    point.X() <= 300  && point.X() >= -300 &&
    point.Y() <= 600 && point.Y() >= 0;
}

// Return the distance between two vectors
template<typename Vector1, typename Vector2>
Double_t MCSAngleAnalysis::distanceBetween(Vector1 vector1, Vector2 vector2) const {

  // Calculate the difference between each vector's components
  Double_t deltaX = vector1.X() - vector2.X();
  Double_t deltaY = vector1.Y() - vector2.Y();
  Double_t deltaZ = vector1.Z() - vector2.Z();

  // Calculate and return the difference between the two vectors by adding the difference in the components in quadrature.
  Double_t diff = sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
  return diff;
} // MCSAngleAnalysis::distanceBetween(Vector1, Vector2)

// Returns the first point in the detector for the specified MCParticle
TVector3 MCSAngleAnalysis::startingPositionFor(simb::MCParticle particle) const {
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); ++i) {
    if(inDetector(particle.Position(i)))
      return particle.Position(i).Vect();
  }

  return TVector3(-999,-999,-999);
}

// Will return a TVector3 with the starting direction of the provided MCParticle for the first trajectory point in the detector.
TVector3 MCSAngleAnalysis::startingDirectionFor(simb::MCParticle particle) const {
  // Loop through trajectory points
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); ++i) {
    // If the trajectory point is in the detector.
    if(inDetector(particle.Position(i))) {
      // Return the first trajectory point in the detector's direction unit vector.
      return particle.Momentum(i).Vect().Unit();
    }
  }

  return TVector3(0,0,0);
}

Double_t MCSAngleAnalysis::highland(Double_t p, Double_t l) const {
  Double_t kappa_a = 0.095; // MeV
  Double_t kappa_c = 11.968; // MeV
  Double_t kappa = kappa_c + kappa_a/(p*p); // MeV (p is unitless but given in GeV/c)
  
  Double_t mass = 0.106; // MeV/c
  Double_t pBetaC = (p*p)/sqrt(mass*mass + p*p);

  return 0.001*((1+0.038*log(l/14))*sqrt(l/14)*kappa)/pBetaC; // 0.001 is for converting mrads to rads
}

// p is given in MeV
// TODO: Change to provide p in GeV/c
Double_t MCSAngleAnalysis::deltaP(Double_t p, Double_t l) const {
  Double_t M = 105.66; // MeV/c^2, mass of muon
  Double_t initialE = pow(M*M+p*p,0.5);
  Double_t Z = 18, A = 39.948; //For Ar, do they use a specific Isotope?
  Double_t K = 0.307075; // MeV*cm^2/mol
  Double_t I = 188*pow(10,-6); //MeV
  Double_t m = 0.511; // MeV/c^2, electron mass
  Double_t bg = p/M, bb = (p*p)/(m*m+p*p), gg = (m*m+p*p)/(m*m); //Beta*Gamma, Beta^2, Gamma^2
  Double_t W = (2*m*bb*gg)/(1+(2*m*pow(gg,0.5)/M)+pow(m/M,2));
  Double_t del, a = 0.19559, x0 = 0.2, x1 = 3, C = 5.2146, k = 3, x = log10(bg);
  if(x >= x1)
    del = 0.5*(2*log(10)*x-C);
  else if(x >= x0 && x <= x1)
    del = 0.5*(2*log(10)*x-C+a*pow(x1+x,k));
  else
    del = 0;
  
  //First calculate dEdx (MeV)/cm
  Double_t dEdx = K*(Z/A)*(1/bb)*(0.5*log(2*m*bb*gg*W/(I*I))-bb-del);
  return p-pow(pow(initialE - dEdx*l,2)-M*M,0.5); //Convert MeV to GeV
}

// TODO: Documentation
template<typename Vector>
size_t MCSAngleAnalysis::TrajectoryPointWithClosestZ(simb::MCParticle particle, Vector recoPoint) const {

  // Define initial distance and index for the first point.
  size_t closestIndex = 0;
  TVector3 firstPoint = particle.Position(0).Vect();
  Double_t closestDistance = distanceBetween(firstPoint, recoPoint);

  // Loop through the rest of the points.
  for(size_t i = 1; i < particle.NumberTrajectoryPoints(); ++i) {
    // Define distance between this point and the recoPoint
    TVector3 truePoint = particle.Position(i).Vect();
    Double_t distance = distanceBetween(truePoint, recoPoint);

    // If the distance is less than the previuos closestDistance, update the closestDistance and closestIndex
    if(distance < closestDistance) {
      closestDistance = distance;
      closestIndex = i;
    } // if
  } // for i

  return closestIndex;
} // TrajectoryPointWithClosestZ

DEFINE_ART_MODULE(MCSAngleAnalysis)
