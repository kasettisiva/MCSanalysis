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
#include <fstream>

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
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
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
  //template<typename Point> Vector_t Calculate2DLinearFit(const std::vector<Point> segment) const;
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

  //Plot the number of points in segments
  TH1F* trueNumberOfSegmentPoints_HIST;
  TH1F* recoNumberOfSegmentPoints_HIST;

  TH2F* recoSegmentPointsXZ_HIST;
  TH2F* recoSegmentPointsYZ_HIST;

  TH2F* trueSegmentPointsXZ_HIST;
  TH2F* trueSegmentPointsYZ_HIST;

  TH1F* recoNumberOfTracks_HIST;

  TH2F* true_numberOfSegmentPointsVsSegmentMomentum_HIST;
  TH2F* reco_numberOfSegmentPointsVsSegmentMomentum_HIST;
  TH2F* alttrue_numberOfSegmentPointsVsSegmentMomentum_HIST;

  TH1F* trueSimIDEdistance_HIST;
  TH1F* simIDEtrackID_HIST;
  TH1F* simIDEedep_HIST;
  TH2F* simIDEedepFromStart_HIST;

  TH1F* deltaXZ_HIST;
  TH1F* deltaYZ_HIST;

  TH2F* trueLinear_deltaXZprimeVsSegmentMomentum_HIST;
  TH2F* trueLinear_deltaYZprimeVsSegmentMomentum_HIST;

  TH1F* trueLinear_transverseDistX_HIST;
  TH1F* trueLinear_transverseDistY_HIST;
  TH1F* recoLinear_transverseDistX_HIST;
  TH1F* recoLinear_transverseDistY_HIST;

  TH2F* trueLinear_transverseDistXVsSegmentMomentum_HIST;
  TH2F* trueLinear_transverseDistYVsSegmentMomentum_HIST;

  TH2F* recoLinear_transverseDistXVsSegmentMomentum_HIST;
  TH2F* recoLinear_transverseDistYVsSegmentMomentum_HIST;

  TH2F* altTrueLinear_ThetaXZprimeVsSegmentMomentum_HIST;
  TH2F* altTrueLinear_ThetaYZprimeVsSegmentMomentum_HIST;

  TH1F* altTrueNumberOfSegmentPoints_HIST;
  TH1F* truePtsRecoMatched_distance_HIST;
  TH1F* recoPts_distance_HIST;

  TH1F* trueMatchedReco_deltaX_HIST;
  TH1F* trueMatchedReco_deltaY_HIST;
  TH1F* trueMatchedReco_deltaZ_HIST;

  TH2F* alttrueSegmentPointsXZ_HIST;
  TH2F* alttrueSegmentPointsYZ_HIST;

  TH2F* reco_thetaXZVsNpts_HIST;
  TH2F* reco_thetaYZVsNpts_HIST;

  TH2F* alttrue_thetaXZVsNpts_HIST;
  TH2F* alttrue_thetaYZVsNpts_HIST;

  TH2F* alttrueSegmentPointsXZ_lt20PtsperSeg_HIST;
  TH2F* alttrueSegmentPointsYZ_lt20PtsperSeg_HIST;

  TH2F* reco_vs_trueX_HIST;
  TH2F* reco_vs_trueY_HIST;
  TH2F* reco_vs_trueZ_HIST;

  TH2F* distFromStartVsDeltaX_HIST;
  TH2F* distFromStartVsDeltaY_HIST;
  TH2F* distFromStartVsDeltaZ_HIST;

  TH2F* distFromStartIfDistAbove3VsDeltaX_HIST;
  TH2F* distFromStartIfDistAbove3VsDeltaY_HIST;
  TH2F* distFromStartIfDistAbove3VsDeltaZ_HIST;

  //TH2F* reco_distFromStartVsDeltaX_HIST;
  //TH2F* reco_distFromStartVsDeltaY_HIST;
  //TH2F* reco_distFromStartVsDeltaZ_HIST;

  TH2F* alttrue_distFromStartVsdistBtwPts_HIST;
  TH2F* reco_distFromStartVsdistBtwPts_HIST;

  TH2F* alttrue_pointIndexVsdistBtwPts_HIST;
  TH2F* reco_pointIndexVsdistBtwPts_HIST;

  TH1F* nSegwithLessthan10pts_HIST;
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

  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  // MCS SegmentCalculator
  trkf::MCSSegmentCalculator segmentCalculator = trkf::MCSSegmentCalculator(fAngleCalculatorParameters().segmentCalculator);
  // MCS Angle Calculator
  trkf::MCSAngleCalculator angleCalculator = trkf::MCSAngleCalculator(fAngleCalculatorParameters);
  // LSU Backtracker
  lsu::BackTrackerAlg backtracker;
  art::ServiceHandle<cheat::BackTrackerService> bt;
  art::ServiceHandle<geo::Geometry> geom;
  //geo::GeometryCore const * fGeometry = &*(art::ServiceHandle<geo::Geometry>());

  // TODO: if(!e.isRealData())
  auto particleListHandle = e.getValidHandle<std::vector<simb::MCParticle> >("largeant");
  std::vector<art::Ptr<simb::MCParticle> > particleList;
  art::fill_ptr_vector(particleList, particleListHandle);

  //std::vector<const sim::SimChannel*> simChannels;
  const std::vector<art::Ptr<sim::SimChannel> > simChannels=bt->SimChannels();
  //try{
    //e.getView("largeant", simChannels);
  //}
  //catch (art::Exception const&evt){
  //}
  //std::cout << " size of simchannels " << simChannels.size() << std::endl;

  Double_t inMom=-999;
  for(size_t iParticle = 0; iParticle < particleList.size(); ++iParticle) {
    art::Ptr<simb::MCParticle> particle = particleList.at(iParticle);
    if(particle->Trajectory().TotalLength() >= 100 && abs(particle->PdgCode()) == 13) {
      inMom = particle->P(0);
      std::cout << " initial mom " << inMom << std::endl;
    }
  }

  std::vector<Point_t> points; 
  std::vector<Point_t> tmp_points; //to make distance vs edep
  double mom=0;
  double eloss=0;
  double inE=TMath::Sqrt(inMom*inMom + 0.106*0.106);

/*
  for(auto const& simchannel :simChannels) {
     auto const& alltimeslices = simchannel->TDCIDEMap(); 
     //std::cout << " size of TDCIDEmap  " << alltimeslices.size() << " channel " << simchannel->Channel() << std::endl;
     //if(simchannel->Channel() != 2) continue;
     std::vector<geo::WireID> chan2wires = fGeometry->ChannelToWire(simchannel->Channel()); 
     //std::cout << " size of wires in each channel " << chan2wires.size() << std::endl;
     for(auto const& wireid : chan2wires) { 
      //std::cout << " Plane ID " << wireid.Plane << std::endl;
      if(wireid.Plane != 2) continue;
      for(auto const& tslice : alltimeslices){
        auto const& simide = tslice.second;
        //std::cout << " size of simIDE " << simide.size() << std::endl;
        // Loop over energy deposits
        for(auto const& eDep : simide){
	  //double mom=-999;
	  //if(eDep.trackID == 1 || eDep.trackID == -1) std::cout << " TrackID of simIDE " << eDep.trackID << std::endl;
	  //if(eDep.trackID == 1 || eDep.trackID == -1) continue;
          //std::cout << " TrackID of simIDE " << eDep.trackID << std::endl;
	  simIDEtrackID_HIST->Fill(eDep.trackID);
          if(eDep.trackID != 1) continue;
          simIDEedep_HIST->Fill(eDep.energy);
	  //std::cout << " energy deposited by muon " << eDep.energy << std::endl;
	  //mom+=TMath::Sqrt(eDep.energy*eDep.energy*1e-6 - 0.106*0.106);
	  eloss += eDep.energy*1e-3; //GeV
	  mom = TMath::Sqrt((inE-eloss)*(inE-eloss) - 0.106*0.106);	  	  
          //std::cout << " inMom:inE " << inMom << " : " << inE  << " mom " << mom << " energy loss " << eloss << " (inE-eloss) " << (inE-eloss) << std::endl;
          Point_t pt(eDep.x, eDep.y, eDep.z, mom);
	  points.push_back(pt); 
	  Point_t pt1(eDep.x, eDep.y, eDep.z, eDep.energy);
          tmp_points.push_back(pt1);
        }
     }
   }
  }
  //std::cout << " number of points using simIDE " << points.size() << std::endl;
  //trkf::MCSSegmentCalculator::MCSSegmentResult trueSegmentResult_fromSimIDE = segmentCalculator.GetResult(points, fShouldCreateVirtualPoints);

  //sort the vector of points created from simIDE's
  int temp1=0; int temp2=0;
  std::for_each(points.begin(), points.end(), [&temp1](Point_t point) {
        //std::cout << "Before sorting: " << temp1 << " x " << point.X() << " y " << point.Y() << " z " << point.Z() << "\n";
        ++temp1;
   });

  std::sort(points.begin(), points.end(),
            [](Point_t i1, Point_t i2) {return (i1.Z() < i2.Z());});

  std::sort(tmp_points.begin(), tmp_points.end(),
            [](Point_t i1, Point_t i2) {return (i1.Z() < i2.Z());});

  std::for_each(points.begin(), points.end(), [&temp2](Point_t point) {
	//std::cout << "After sorting: " << temp2 << " x " << point.X() << " y " << point.Y() << " z " << point.Z() << "\n";
        ++temp2;
   });
*/

  // Get SimEnergyDeposit
  std::vector<art::Ptr<sim::SimEnergyDeposit>> sedlist;
  auto sedListHandle = e.getHandle< std::vector<sim::SimEnergyDeposit> >("largeant:LArG4DetectorServicevolTPCActive");
  if (sedListHandle){
    art::fill_ptr_vector(sedlist, sedListHandle);
  }
  else{
    std::cout<<"sedListHandle invalid"<<std::endl;
    return;
  }

  for (size_t i = 0; i<sedlist.size(); ++i){
    double x = sedlist[i]->X();
    double y = sedlist[i]->Y();
    double z = sedlist[i]->Z();
    double edep = sedlist[i]->E();
    int trkid = sedlist[i]->TrackID();
    //int pdg = sedlist[i]->PdgCode();
    //cout<<x<<" "<<y<<" "<<z<<" "<<trkid<<" "<<pdg<<endl;
    if (trkid!=1) continue;
    eloss += edep*1e-3; //GeV
    mom = TMath::Sqrt((inE-eloss)*(inE-eloss) - 0.106*0.106);               
    //std::cout << " inMom:inE " << inMom << " : " << inE  << " mom " << mom << " energy loss " << eloss << " (inE-eloss) " << (inE-eloss) << std::endl;
    Point_t pt(x,y,z,mom);
    points.push_back(pt); 
  }

  //tmp_points.push_back(points[0])
  // Loop through the trajectory points
  for(size_t i = 1; i < points.size(); ++i) {
    // Define the current and previous trajectory point.
    Point_t prevPoint = points.at(i-1);
    Point_t currPoint = points.at(i);
    //std::cout << " i " << i << std::endl;
    // If this point is invalid.
    if(currPoint.X() == -999)
      break;
    // Define the difference between the current and previous trajectory point.
    Float_t diff = sqrt((currPoint.underlying() - prevPoint.underlying()).Mag2());
    trueSimIDEdistance_HIST->Fill(diff);
 
    //if(i%6==0) tmp_points.push_back(currPoint);

  }
/*
  for(size_t i = 1; i < tmp_points.size(); ++i) {
    // Define the current and previous trajectory point.
    Point_t prevPoint = tmp_points.at(0);
    Point_t currPoint = tmp_points.at(i);
    //std::cout << " i " << i << std::endl;
    // If this point is invalid.
    if(currPoint.X() == -999)
      break;
    // Define the difference between the current and previous trajectory point.
    Float_t diff = sqrt((currPoint.underlying() - prevPoint.underlying()).Mag2());
    simIDEedepFromStart_HIST->Fill(diff,currPoint.P());
  }
*/

  std::vector<art::Ptr<recob::SpacePoint>> splist;
  auto spListHandle = e.getHandle< std::vector<recob::SpacePoint> >("pandora");
  if (spListHandle){
    art::fill_ptr_vector(splist, spListHandle);
  }
  else{
    std::cout<<"SpacePointListHandle invalid"<<std::endl;
    return;
  }

  std::vector<Point_t> sp_points;
  for(size_t i = 0; i < splist.size(); ++i) {
     double x=splist[i]->XYZ()[0];
     double y=splist[i]->XYZ()[1];
     double z=splist[i]->XYZ()[2];
     Point_t pt(x,y,z);
     sp_points.push_back(pt);
  }

  // Loop through MCParticles
  size_t nParticles = particleList.size();
  for(size_t iParticle = 0; iParticle < nParticles; ++iParticle) {
    art::Ptr<simb::MCParticle> particle = particleList.at(iParticle);
    //std::cout << " true traj length " << particle->Trajectory().TotalLength() << " particle->PdgCode() "  << particle->PdgCode() << " momentum " << particle->P() << std::endl;
    if(particle->Trajectory().TotalLength() >= 100 && abs(particle->PdgCode()) == 13) {
      // Get an MCSSegmentResult
      //trkf::MCSSegmentCalculator::MCSSegmentResult trueSegmentResult = segmentCalculator.GetResult(*particle, fShouldCreateVirtualPoints);
      trkf::MCSSegmentCalculator::MCSSegmentResult trueSegmentResult = segmentCalculator.GetResult(points, fShouldCreateVirtualPoints);
      //trkf::MCSSegmentCalculator::MCSSegmentResult trueSegmentResult = segmentCalculator.GetResult(tmp_points, fShouldCreateVirtualPoints);

      std::vector<Segment_t> trueSegment_vec = trueSegmentResult.GetSegment_vec();
      std::cout << " trueSegment_vec size " << trueSegment_vec.size() << std::endl;
      for(size_t i=0;i<trueSegment_vec.size();++i) {
        Double_t trueSegmentMomentum = trueSegmentResult.GetMomentumForSegment(i);
        std::cout << " i " << i << " trueSegmentMomentum " << trueSegmentMomentum << std::endl;
        trueNumberOfSegmentPoints_HIST->Fill(trueSegment_vec.at(i).size());
        for(size_t j=0;j<trueSegment_vec.at(i).size();++j) {
          trueSegmentPointsXZ_HIST->Fill(trueSegment_vec.at(i).at(j).Z(), trueSegment_vec.at(i).at(j).X());
          trueSegmentPointsYZ_HIST->Fill(trueSegment_vec.at(i).at(j).Z(), trueSegment_vec.at(i).at(j).Y());
        }
      }

      Segment_t first_trueSegment = trueSegment_vec.at(0); 
      for(size_t i=1;i<trueSegment_vec.size();++i) {
          Double_t segmentMomentum = trueSegmentResult.GetMomentumForSegment(i);
	  Segment_t trueSegment = trueSegment_vec.at(i);
          Double_t atZ = trueSegment.at(trueSegment.size()-1).Z(); 	  
          //Double_t xatZ = trueSegment.at(trueSegment.size()-1).X();
          //Double_t yatZ = trueSegment.at(trueSegment.size()-1).Y();

          std::vector<Vector_t> fitVec_segment = segmentCalculator.LinearFit_2D(trueSegment);
          TVector3 fitVecX_segment = segmentCalculator.convert<Vector_t, TVector3>(fitVec_segment.at(0));
          TVector3 fitVecY_segment = segmentCalculator.convert<Vector_t, TVector3>(fitVec_segment.at(1));
          Double_t xatZ = fitVecX_segment.X()*atZ + fitVecX_segment.Z();
          Double_t yatZ = fitVecY_segment.Y()*atZ + fitVecY_segment.Z();

          //Vector_t fitVec = segmentCalculator.LinearFit_2D(first_trueSegment);
	  //TVector3 fitVecX = segmentCalculator.convert<Vector_t, TVector3>(segmentCalculator.LinearFit_2D(first_trueSegment));
          //TVector3 fitVecY = segmentCalculator.convert<Vector_t, TVector3>(segmentCalculator.LinearFit_2D(first_trueSegment,1));
          std::vector<Vector_t> fitVec = segmentCalculator.LinearFit_2D(first_trueSegment);
	  TVector3 fitVecX = segmentCalculator.convert<Vector_t, TVector3>(fitVec.at(0));
          TVector3 fitVecY = segmentCalculator.convert<Vector_t, TVector3>(fitVec.at(1));

	  Double_t x1atZ = fitVecX.X()*atZ + fitVecX.Z();
          Double_t y1atZ = fitVecY.Y()*atZ + fitVecY.Z();

	  Double_t transverseDistX = xatZ - x1atZ; 
          Double_t transverseDistY = yatZ - y1atZ;

	  trueLinear_transverseDistX_HIST->Fill(transverseDistX);
          trueLinear_transverseDistY_HIST->Fill(transverseDistY);

	  trueLinear_transverseDistXVsSegmentMomentum_HIST->Fill(segmentMomentum, transverseDistX);
          trueLinear_transverseDistYVsSegmentMomentum_HIST->Fill(segmentMomentum, transverseDistY);

          //std::cout << " deltaX " << transverseDistX << std::endl;

   	  first_trueSegment = trueSegment_vec.at(i);
      }

      // Extract data from the true segment result.
      std::vector<Float_t> trueRawSegmentLength_vec = trueSegmentResult.GetRawSegmentLength_vec();

      // Get an MCSAngleResult for various methods using the true segment result.
      trkf::MCSAngleCalculator::MCSAngleResult trueLinearMCS = angleCalculator.GetResult(trueSegmentResult, 0);
      trkf::MCSAngleCalculator::MCSAngleResult truePolygonalMCS = angleCalculator.GetResult(trueSegmentResult, 1);

      // Extract the data from the various angle results.
      std::vector<Angles_t> trueLinearAngles_vec = trueLinearMCS.GetAngles_vec();
      std::vector<Angles_t> truePolygonalAngles_vec = truePolygonalMCS.GetAngles_vec();

      // Loop through angle measurements for each segment.
      for(size_t i = 0; i < trueLinearAngles_vec.size(); ++i) {
	// TODO: Assumption: Relative sizes of vectors
	// Momentum for segment. The i+1 is here since our first angle is for the second segment.
	Double_t segmentMomentum = trueSegmentResult.GetMomentumForSegment(i+1);
	Double_t segmentLength = trueRawSegmentLength_vec.at(i+1);
	Segment_t trueSegment = trueSegment_vec.at(i+1);
	// TODO: Invalid momentum for polygonal segments.

	// Angles for segment.
	Angles_t linearAngles = trueLinearAngles_vec.at(i);
	Angles_t polygonalAngles = truePolygonalAngles_vec.at(i);

	// Get individual angles for this segment.
	Float_t linearTheta3D = std::get<0>(linearAngles)*1000.0;
	Float_t linearThetaXZprime = std::get<1>(linearAngles)*1000.0;
	Float_t linearThetaYZprime = std::get<2>(linearAngles)*1000.0;
	Float_t polygonalTheta3D = std::get<0>(polygonalAngles)*1000.0;
	Float_t polygonalThetaXZprime = std::get<1>(polygonalAngles)*1000.0;
	Float_t polygonalThetaYZprime = std::get<2>(polygonalAngles)*1000.0;

	Float_t deltaXZ = 14/TMath::Sqrt(1+std::pow(1/TMath::Tan(linearThetaXZprime),2)*std::pow(1/TMath::Cos(linearThetaYZprime),2));
        Float_t deltaYZ = 14/TMath::Sqrt(1+std::pow(1/TMath::Tan(linearThetaYZprime),2)*std::pow(1/TMath::Cos(linearThetaXZprime),2));

        deltaXZ_HIST->Fill((linearThetaXZprime>0)?deltaXZ:-deltaXZ);
        deltaYZ_HIST->Fill((linearThetaYZprime>0)?deltaYZ:-deltaYZ);

        trueLinear_deltaXZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, (linearThetaXZprime>0)?deltaXZ:-deltaXZ);
        trueLinear_deltaYZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, (linearThetaYZprime>0)?deltaYZ:-deltaYZ);

        true_numberOfSegmentPointsVsSegmentMomentum_HIST->Fill(segmentMomentum, trueSegment.size());
	// Directly Compare 3D angles with projected angles added in quadrature. Should be one-to-one.
	// TODO: Assumption: Small angle assumption
	trueLinear_theta3DVsTheta2DInQuadrature_HIST->Fill(pow(linearThetaXZprime*linearThetaXZprime + linearThetaYZprime*linearThetaYZprime,0.5), linearTheta3D);
	truePolygonal_theta3DVsTheta2DInQuadrature_HIST->Fill(pow(polygonalThetaXZprime*polygonalThetaXZprime + polygonalThetaYZprime*polygonalThetaYZprime,0.5), polygonalTheta3D);

	// Fill histograms with angles vs. segment momentum.
	trueLinear_theta3DVsSegmentMomentum_HIST->Fill(segmentMomentum, linearTheta3D);
	trueLinear_thetaXZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, linearThetaXZprime);
	trueLinear_thetaYZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, linearThetaYZprime);
	truePolygonal_theta3DVsSegmentMomentum_HIST->Fill(segmentMomentum, polygonalTheta3D);
	truePolygonal_thetaXZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, polygonalThetaXZprime);
	truePolygonal_thetaYZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, polygonalThetaYZprime);

	// Update counters that count the number of angles within n*sigma.
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
  //Get associations between tracks and hits
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(
    trackListHandle,
    e, //event
    fTrackModuleLabel);

  std::vector<Point_t> points_vec;
  std::vector<Point_t> points_vec_rec;
  std::ofstream myfile;
  myfile.open ("example.txt");

  std::cout << "Number of Tracks: " << trackList.size() << std::endl;
  // Loop through the list of reconstructed tracks
  size_t nTracks = trackList.size();
  recoNumberOfTracks_HIST->Fill(nTracks);

  for(size_t iTrack = 0; iTrack < nTracks; iTrack++) {
    points_vec.clear();
    points_vec_rec.clear();
    Int_t nSegwithLessthan10pts=0;
    // Get this specific track
    art::Ptr<recob::Track> track = trackList.at(iTrack);

    // (attempt to) Backtrack this track to an MCParticle
    // TODO: Backtracking error handling.
    simb::MCParticle particle = backtracker.getMCParticle(track, e, fTrackModuleLabel.label());

    std::cout << " track length " << track->Length() << std::endl;
    // Only run this analysis if the track length is greater than 100 cm
    if(track->Length() >= 100) {

    if (fmthm.isValid()) {
        auto vhit = fmthm.at(iTrack);
        auto vmeta = fmthm.data(iTrack);
        //std::cout << " size of recob::Track trajectory " << track->Trajectory().Trajectory().Positions().size() << " size of vhit " << vhit.size() << " " << vmeta.size() << std::endl; 
	std::cout << " size of recob::Track trajectory " << track->Trajectory().Trajectory().Positions().size() << " track->NPoints() " << track->NPoints() << std::endl;
        //for(recob::tracking::Point_t position: track->Trajectory().Trajectory().Positions()) {
        for(unsigned int position=0; position < track->Trajectory().Trajectory().Positions().size();++position) {
            //std::cout << " position " << position << std::endl;
	    Double_t recox = track->Trajectory().Trajectory().Positions().at(position).X();
            Double_t recoy = track->Trajectory().Trajectory().Positions().at(position).Y();
            Double_t recoz = track->Trajectory().Trajectory().Positions().at(position).Z();

            //std::cout << " reco position index " << position << " " << recox << " " << recoy << " " << recoz << std::endl;

	    for(size_t ii = 1; ii < vhit.size(); ++ii) { 
		//std::cout << " ii " << ii << " index is " << vmeta[ii]->Index() << std::endl;
                std::vector<const sim::IDE *> hittosimide = bt->HitToSimIDEs_Ps(clock_data,vhit.at(ii));
		//std::vector<int> hittotrkid = bt->HitToTrackIds(clock_data,*vhit.at(ii));

		if(hittosimide.size()==0) continue; //not sure why there are no simIDEs

                //Double_t _truex = bt->HitToXYZ(clock_data,vhit.at(ii))[0];
                //Double_t _truey = bt->HitToXYZ(clock_data,vhit.at(ii))[1];
                //Double_t _truez = bt->HitToXYZ(clock_data,vhit.at(ii))[2];
                //std::cout << " true position index " << ii << " matching index " << vmeta[ii]->Index() << " x "  << _truex << " " << _truey << " " << _truez << std::endl;

		if(vmeta[ii]->Index() == position) {
		    //std::cout << " found the correct hit " << std::endl;
		    //std::vector<const sim::IDE *> hittosimide = bt->HitToSimIDEs_Ps(clock_data,vhit.at(ii));
                    //std::cout << " simide size " << hittosimide.size() << std::endl;
		    //for(long unsigned int isim=0;isim<hittosimide.size();++isim) {
			//std::cout << " isim " << isim << " nelectrons " << hittosimide[isim]->numElectrons << std::endl;
		    //}

		    //for(long unsigned int id=0;id<hittotrkid.size();++id) {
 			//std::cout << " trackid " << hittotrkid[id] << " for reco point " << position << std::endl; 
   		    //}

		    Double_t truex = bt->HitToXYZ(clock_data,vhit.at(ii))[0];
                    Double_t truey = bt->HitToXYZ(clock_data,vhit.at(ii))[1];
                    Double_t truez = bt->HitToXYZ(clock_data,vhit.at(ii))[2];

		    Point_t tmp_pt(truex,truey,truez);
                    Point_t tmp_pt1(recox,recoy,recoz);
		    //Point_t tmp_pt(bt->HitToXYZ(clock_data,vhit[ii])[0],bt->HitToXYZ(clock_data,vhit[ii])[1],bt->HitToXYZ(clock_data,vhit[ii])[2]);
            	    //points_vec.push_back(Point_t(track->Trajectory().Trajectory().Positions().at(position)));

		    //std::cout << " matched " << recox << " " << recoy << " " << recoz << " " << truex << " " << truey << " " << truez << " wireID: " << std::string(vhit[ii]->WireID()) << std::endl;

                    if(abs(truex - recox)>5 && abs(truey - recoy)>5 && abs(truez - recoz)>5) std::cout << "Event number with any differnce > 5: " << e.event() << " trueX:trueY:trueZ " << truex << " : " << truey << " : " << truez << " recox:recoy:recoz " << recox << " : " << recoy << " : " << recoz << std::endl;

                    trueMatchedReco_deltaX_HIST->Fill(truex - recox); 
                    trueMatchedReco_deltaY_HIST->Fill(truey - recoy);                    
		    trueMatchedReco_deltaZ_HIST->Fill(truez - recoz);

	            alttrueSegmentPointsXZ_HIST->Fill(truez,truex);
          	    alttrueSegmentPointsYZ_HIST->Fill(truez,truey);

                    recoSegmentPointsXZ_HIST->Fill(recoz,recox);
                    recoSegmentPointsYZ_HIST->Fill(recoz,recoy);

                    reco_vs_trueX_HIST->Fill(truex,recox);
                    reco_vs_trueY_HIST->Fill(truey,recoy);
                    reco_vs_trueZ_HIST->Fill(truez,recoz);

		    points_vec.push_back(tmp_pt);
                    points_vec_rec.push_back(tmp_pt1);
		}
	    }
        }
    }
    std::cout << " points vec size " << points_vec.size() << std::endl;

    Double_t distFromStart=0;
    Double_t distFromStart_rec=0;
    for(size_t i=1;i<points_vec.size();++i) {
          //auto &loc=points_vec.at(i);
          //alttrueSegmentPointsXZ_HIST->Fill(loc.Z(), loc.X());
          //alttrueSegmentPointsYZ_HIST->Fill(loc.Z(), loc.Y());

          Point_t prevPoint = points_vec.at(i-1);
          Point_t currPoint = points_vec.at(i);
          Float_t distBtwPts = sqrt((currPoint.underlying() - prevPoint.underlying()).Mag2());
	  distFromStart += distBtwPts;

	  alttrue_distFromStartVsdistBtwPts_HIST->Fill(distFromStart,distBtwPts);
          alttrue_pointIndexVsdistBtwPts_HIST->Fill(i,distBtwPts);
	  recoPts_distance_HIST->Fill(distBtwPts);

	  Double_t truex = currPoint.X();
          Double_t truey = currPoint.Y();
          Double_t truez = currPoint.Z();

	  if((distBtwPts>12) && (distBtwPts<16)) std::cout << "Event number with distance btw pts around 14: " << e.event() << " distance is " << distBtwPts << std::endl;

          Point_t prevPoint_rec = points_vec_rec.at(i-1);
          Point_t currPoint_rec = points_vec_rec.at(i);
          Float_t distBtwPts_rec = sqrt((currPoint_rec.underlying() - prevPoint_rec.underlying()).Mag2());
          distFromStart_rec += distBtwPts_rec;

          reco_distFromStartVsdistBtwPts_HIST->Fill(distFromStart_rec,distBtwPts_rec);
          reco_pointIndexVsdistBtwPts_HIST->Fill(i,distBtwPts_rec);
          truePtsRecoMatched_distance_HIST->Fill(distBtwPts_rec);

          Double_t recox = currPoint_rec.X();
          Double_t recoy = currPoint_rec.Y();
          Double_t recoz = currPoint_rec.Z();

          distFromStartVsDeltaX_HIST->Fill(i, currPoint.X() - currPoint_rec.X());
          distFromStartVsDeltaY_HIST->Fill(i, currPoint.Y() - currPoint_rec.Y());
          distFromStartVsDeltaZ_HIST->Fill(i, currPoint.Z() - currPoint_rec.Z());

          if(i==1) std::cout << " point 0 "
                    << " true:reco X " << prevPoint.X() << " : " << prevPoint_rec.X() << " : " << prevPoint.X() - prevPoint_rec.X() << " "
                    << " true:reco Y " << prevPoint.Y() << " : " << prevPoint_rec.Y() << " : " << prevPoint.Y() - prevPoint_rec.Y() << " " 
                    << " true:reco Z " << prevPoint.Z() << " : " << prevPoint_rec.Z() << " : " << prevPoint.Z() - prevPoint_rec.Z()
                    << std::endl;

          std::cout << " true distance " << distBtwPts << " reco " << distBtwPts_rec
                    << " true:reco X " << truex << " : " << recox << " " << truex - recox
                    << " true:reco Y " << truey << " : " << recoy << " " << truey - recoy
                    << " true:reco Z " << truez << " : " << recoz << " " << truez - recoz
                    << std::endl;

	  std::cout << i << " " << truex << " " << truey << " " << truez << " " << recox << " " << recoy << " " << recoz << std::endl; 
          myfile << i << " " << truex << " " << truey << " " << truez << " " << recox << " " << recoy << " " << recoz << " " << (truex-recox) << " " << (truey-recoy) << " " << (truez-recoz) << "\n";

	  if(distBtwPts > (2*distBtwPts_rec)) { 
	        //std::cout << " == the distance btw true points > reco points: true distance " << distBtwPts << " reco " << distBtwPts_rec
		//          << " true:reco X " << truex << " : " << recox << " " << truex - recox
		//          << " true:reco Y " << truey << " : " << recoy << " " << truey - recoy
		//          << " true:reco Z " << truez << " : " << recoz << " " << truez - recoz
		//          << std::endl; 
                distFromStartIfDistAbove3VsDeltaX_HIST->Fill(i, currPoint.X() - currPoint_rec.X());
                distFromStartIfDistAbove3VsDeltaY_HIST->Fill(i, currPoint.Y() - currPoint_rec.Y());
                distFromStartIfDistAbove3VsDeltaZ_HIST->Fill(i, currPoint.Z() - currPoint_rec.Z());
   	  }

	  //std::cout << i << " distBtwPts " << distBtwPts << " distFromStart " << distFromStart << " distBtwPts_rec " << distBtwPts_rec << " distFromStart_rec " << distFromStart_rec << std::endl;

    }

      // Calculate Reco MCSSegmentResult
      trkf::MCSSegmentCalculator::MCSSegmentResult AltTrueSegmentResult = segmentCalculator.GetResult(points_vec, fShouldCreateVirtualPoints);
      // Calculate Reco MCSAngleResult for various methods
      trkf::MCSAngleCalculator::MCSAngleResult altTrueLinearAngleResult = angleCalculator.GetResult(AltTrueSegmentResult, 0);
      // Extract data from the reco segment result
      std::vector<Segment_t> alttrue_segment_vec = AltTrueSegmentResult.GetSegment_vec();
      std::vector<Vector_t> alttrue_linearFit_vec = AltTrueSegmentResult.GetLinearFit_vec();

      std::cout << " alttrue_segment_vec size " << alttrue_segment_vec.size() << std::endl;
      for(size_t i=0;i<alttrue_segment_vec.size()-1;++i) {
        altTrueNumberOfSegmentPoints_HIST->Fill(alttrue_segment_vec.at(i).size());

        std::cout << " Event: " << e.event()
                  << " segment number " << i
                  << " no.of alt.true points per segment " << alttrue_segment_vec.at(i).size()
                  << std::endl;
	if(alttrue_segment_vec.at(i).size()<10) {
		std::cout << " Event: " << e.event() 
			  << " segment number " << i
			  << " no.of alt.true points per segment less than 10 " << alttrue_segment_vec.at(i).size() 
			  << std::endl;
		++nSegwithLessthan10pts;
	}

        for(size_t j=1;j<alttrue_segment_vec.at(i).size();++j) {

 	  //alttrueSegmentPointsXZ_HIST->Fill(alttrue_segment_vec.at(i).at(j).Z(), alttrue_segment_vec.at(i).at(j).X());
          //alttrueSegmentPointsYZ_HIST->Fill(alttrue_segment_vec.at(i).at(j).Z(), alttrue_segment_vec.at(i).at(j).Y());
	  if(alttrue_segment_vec.at(i).size()<10) {
	     alttrueSegmentPointsXZ_lt20PtsperSeg_HIST->Fill(alttrue_segment_vec.at(i).at(j).Z(), alttrue_segment_vec.at(i).at(j).X());
	     alttrueSegmentPointsYZ_lt20PtsperSeg_HIST->Fill(alttrue_segment_vec.at(i).at(j).Z(), alttrue_segment_vec.at(i).at(j).Y());
	  }
/*
          Point_t prevPoint = alttrue_segment_vec.at(i).at(j-1);
          Point_t currPoint = alttrue_segment_vec.at(i).at(j);
          //std::cout << " i " << i << std::endl;
          // If this point is invalid.
          if(currPoint.X() == -999)
          break;
          // Define the difference between the current and previous trajectory point.
          Float_t diff = sqrt((currPoint.underlying() - prevPoint.underlying()).Mag2());
          truePtsRecoMatched_distance_HIST->Fill(diff);
	  */
        }
      } 

      // Extract data from the various reco angle results
      std::vector<Angles_t> altTrueLinearAngles_vec = altTrueLinearAngleResult.GetAngles_vec();
      // Loop through angle measurements for each segment.
      for(size_t i = 0; i < altTrueLinearAngles_vec.size(); ++i) {

        // Get the segment momentum of this reco segment using the "closest Z" method, described here: TODO: Move to general documentation.
        // Get the first reco point in this segment.  
        // Find the MCParticle trajectory point index that has the closest position with that reco point.
        // Get the momentum of that MCParticle traj. pt. index
        // TODO: Assumption: Relative Vector sizes.
        Segment_t altSegment = alttrue_segment_vec.at(i+1);
        Point_t altFirstPoint = altSegment.at(0);
        size_t truePointIndex = TrajectoryPointWithClosestZ(particle, altFirstPoint);
        Double_t segmentMomentum = particle.P(truePointIndex);
        //std::cout << " segmentMomentum " << segmentMomentum << std::endl;
        // Get the angles for this segment.
        Angles_t altTrueLinearAngles = altTrueLinearAngles_vec.at(i);

        // Extract the individual angles for this segment.
        Float_t linearThetaXZprime = std::get<1>(altTrueLinearAngles)*1000.0;
        Float_t linearThetaYZprime = std::get<2>(altTrueLinearAngles)*1000.0;

        // Fill histograms with angles vs. segment momentum.
        altTrueLinear_ThetaXZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, linearThetaXZprime);
        altTrueLinear_ThetaYZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, linearThetaYZprime);
        alttrue_numberOfSegmentPointsVsSegmentMomentum_HIST->Fill(segmentMomentum, altSegment.size());

	alttrue_thetaXZVsNpts_HIST->Fill(altSegment.size(), linearThetaXZprime);
        alttrue_thetaYZVsNpts_HIST->Fill(altSegment.size(), linearThetaYZprime);
     } 
     /// End of alternative true segment from reco trajectory points
     //============================


      // Calculate Reco MCSSegmentResult
      //trkf::MCSSegmentCalculator::MCSSegmentResult recoSegmentResult = segmentCalculator.GetResult(*track, fShouldCreateVirtualPoints);
      //trkf::MCSSegmentCalculator::MCSSegmentResult recoSegmentResult = segmentCalculator.GetResult(points_vec_rec, fShouldCreateVirtualPoints);
      trkf::MCSSegmentCalculator::MCSSegmentResult recoSegmentResult = segmentCalculator.GetResult(sp_points, fShouldCreateVirtualPoints);
      // Extract data from the reco segment result
      std::vector<Float_t> rawSegmentLength_vec = recoSegmentResult.GetRawSegmentLength_vec();
      std::vector<Segment_t> reco_segment_vec = recoSegmentResult.GetSegment_vec();
      std::vector<Vector_t> reco_linearFit_vec = recoSegmentResult.GetLinearFit_vec();

      std::cout << " reco_segment_vec size " << reco_segment_vec.size() << " track->NPoints() " << track->NPoints() << std::endl;
/*
      for(size_t i=0;i<track->NPoints()-1;++i) {
	  auto &loc = track->LocationAtPoint(i);
          recoSegmentPointsXZ_HIST->Fill(loc.Z(), loc.X());
          recoSegmentPointsYZ_HIST->Fill(loc.Z(), loc.Y());
      }
*/

      for(size_t i=0;i<reco_segment_vec.size()-1;++i) {
        recoNumberOfSegmentPoints_HIST->Fill(reco_segment_vec.at(i).size());
        std::cout << " Event:  " << e.event() 
                  << " segment number " << i
		  <<  " no.of points per reco segment " << reco_segment_vec.at(i).size() 
		  << std::endl;
/*        for(size_t j=1;j<reco_segment_vec.at(i).size();++j) {
          //recoSegmentPointsXZ_HIST->Fill(reco_segment_vec.at(i).at(j).Z(), reco_segment_vec.at(i).at(j).X());
          //recoSegmentPointsYZ_HIST->Fill(reco_segment_vec.at(i).at(j).Z(), reco_segment_vec.at(i).at(j).Y());

          Point_t prevPoint = reco_segment_vec.at(i).at(j-1);
          Point_t currPoint = reco_segment_vec.at(i).at(j);
          //std::cout << " i " << i << std::endl;
          // If this point is invalid.
          if(currPoint.X() == -999)
          break;
          // Define the difference between the current and previous trajectory point.
          Float_t diff = sqrt((currPoint.underlying() - prevPoint.underlying()).Mag2());
          recoPts_distance_HIST->Fill(diff);

        }
*/      }

      Segment_t first_recoSegment = reco_segment_vec.at(0); 
      for(size_t i=1;i<reco_segment_vec.size();++i) {
          Segment_t recoSegment = reco_segment_vec.at(i);
          Point_t recoFirstPoint = recoSegment.at(0);
          size_t truePointIndex = TrajectoryPointWithClosestZ(particle, recoFirstPoint);
          Double_t segmentMomentum = particle.P(truePointIndex);

          Double_t atZ = recoSegment.at(recoSegment.size()-1).Z();        
          //Double_t xatZ = recoSegment.at(recoSegment.size()-1).X();
          //Double_t yatZ = recoSegment.at(recoSegment.size()-1).Y();

          std::vector<Vector_t> fitVec_segment = segmentCalculator.LinearFit_2D(recoSegment);
          TVector3 fitVecX_segment = segmentCalculator.convert<Vector_t, TVector3>(fitVec_segment.at(0));
          TVector3 fitVecY_segment = segmentCalculator.convert<Vector_t, TVector3>(fitVec_segment.at(1));
  	  Double_t xatZ = fitVecX_segment.X()*atZ + fitVecX_segment.Z();
          Double_t yatZ = fitVecY_segment.Y()*atZ + fitVecY_segment.Z();

          //Vector_t fitVec = segmentCalculator.LinearFit_2D(first_recoSegment);
	  std::vector<Vector_t> fitVec = segmentCalculator.LinearFit_2D(first_recoSegment);
          TVector3 fitVecX = segmentCalculator.convert<Vector_t, TVector3>(fitVec.at(0));
          TVector3 fitVecY = segmentCalculator.convert<Vector_t, TVector3>(fitVec.at(1));
          //TVector3 fitVecY = segmentCalculator.convert<Vector_t, TVector3>(segmentCalculator.LinearFit_2D(first_recoSegment,0));

          Double_t x1atZ = fitVecX.X()*atZ + fitVecX.Z();
          Double_t y1atZ = fitVecY.Y()*atZ + fitVecY.Z();
	
	  //std::cout << " fitVecY.Y() " << fitVecY.Y() << " fitVecY.Z() " << fitVecY.Z() << " fitVecY.X() " << fitVecY.X() << " yatZ " << yatZ << " y1atZ " << y1atZ << " fitVecX.X() " << fitVecX.X() << " fitVecX.Y() " << fitVecX.Y() << " fitVecX.Z() " << fitVecX.Z() << " xatZ " << xatZ << " x1atZ " << x1atZ << std::endl;
		
          Double_t transverseDistX = xatZ - x1atZ;
          Double_t transverseDistY = yatZ - y1atZ;

          recoLinear_transverseDistX_HIST->Fill(transverseDistX);
          recoLinear_transverseDistY_HIST->Fill(transverseDistY);

          recoLinear_transverseDistXVsSegmentMomentum_HIST->Fill(segmentMomentum, transverseDistX);
          recoLinear_transverseDistYVsSegmentMomentum_HIST->Fill(segmentMomentum, transverseDistY);

          //std::cout << " deltaX " << transverseDistX << std::endl;
          
          first_recoSegment = reco_segment_vec.at(i);
      } 

      // Calculate Reco MCSAngleResult for various methods
      trkf::MCSAngleCalculator::MCSAngleResult recoLinearAngleResult = angleCalculator.GetResult(recoSegmentResult, 0);
      trkf::MCSAngleCalculator::MCSAngleResult recoPolygonalAngleResult = angleCalculator.GetResult(recoSegmentResult, 1);
      // TODO: Alternative segment defintions

      // Extract data from the various reco angle results
      std::vector<Angles_t> recoLinearAngles_vec = recoLinearAngleResult.GetAngles_vec();
      std::vector<Angles_t> recoPolygonalAngles_vec = recoPolygonalAngleResult.GetAngles_vec();
     
      //std::cout << " reco_segment_vec size " << reco_segment_vec.size() << " recoLinearAngles_vec size " << recoLinearAngles_vec.size() << std::endl;
 
      // Calculate the True MCSSegmentResult for the backtracked MCParticle.
      trkf::MCSSegmentCalculator::MCSSegmentResult trueSegmentResult = segmentCalculator.GetResult(particle, fShouldCreateVirtualPoints);

      // Extract data from the true segment result
      std::vector<Segment_t> true_segment_vec = trueSegmentResult.GetSegment_vec();
      std::vector<Float_t> true_rawSegmentLength_vec = trueSegmentResult.GetRawSegmentLength_vec();
      std::vector<Vector_t> true_linearFit_vec = trueSegmentResult.GetLinearFit_vec();

      // Calculate an True MCSAngleResult for various methods
      trkf::MCSAngleCalculator::MCSAngleResult trueLinearAngleResult = angleCalculator.GetResult(particle, 0);
      // TODO: truePolygonalAngleResult

      // Extract data from the various true angle results
      std::vector<Angles_t> trueLinearAngles_vec = trueLinearAngleResult.GetAngles_vec();

      // Loop through angle measurements for each segment.
      for(size_t i = 0; i < recoLinearAngles_vec.size(); ++i) {

	// Get the segment momentum of this reco segment using the "closest Z" method, described here: TODO: Move to general documentation.
	// Get the first reco point in this segment.  
	// Find the MCParticle trajectory point index that has the closest position with that reco point.
	// Get the momentum of that MCParticle traj. pt. index
	// TODO: Assumption: Relative Vector sizes.
	Segment_t recoSegment = reco_segment_vec.at(i+1);
	Point_t recoFirstPoint = recoSegment.at(0);
	size_t truePointIndex = TrajectoryPointWithClosestZ(particle, recoFirstPoint);
	Double_t segmentMomentum = particle.P(truePointIndex);
        //std::cout << " reco segmentMomentum " << segmentMomentum << std::endl;
	// Get the angles for this segment.
	Angles_t recoLinearAngles = recoLinearAngles_vec.at(i);
	Angles_t recoPolygonalAngles = recoPolygonalAngles_vec.at(i);

	// Extract the individual angles for this segment.
	Float_t linearTheta3D = std::get<0>(recoLinearAngles)*1000.0;
	Float_t linearThetaXZprime = std::get<1>(recoLinearAngles)*1000.0;
	Float_t linearThetaYZprime = std::get<2>(recoLinearAngles)*1000.0;
	Float_t polygonalTheta3D = std::get<0>(recoPolygonalAngles)*1000.0;
	Float_t polygonalThetaXZprime = std::get<1>(recoPolygonalAngles)*1000.0;
	Float_t polygonalThetaYZprime = std::get<2>(recoPolygonalAngles)*1000.0;

        reco_thetaXZVsNpts_HIST->Fill(recoSegment.size(),linearThetaXZprime);
        reco_thetaYZVsNpts_HIST->Fill(recoSegment.size(),linearThetaYZprime);

	// Directly Compare 3D angles with projected angles added in quadrature. Should be one-to-one.
	// TODO: Assumption: small angle assumption
	recoLinear_theta3DVsTheta2DInQuadrature_HIST->Fill(pow(linearThetaXZprime*linearThetaXZprime + linearThetaYZprime*linearThetaYZprime,0.5), linearTheta3D);
	recoPolygonal_theta3DVsTheta2DInQuadrature_HIST->Fill(pow(polygonalThetaXZprime*polygonalThetaXZprime + polygonalThetaYZprime*polygonalThetaYZprime,0.5), polygonalTheta3D);

	// Fill histograms with angles vs. segment momentum.
	recoLinear_theta3DVsSegmentMomentum_HIST->Fill(segmentMomentum, linearTheta3D);
	recoLinear_thetaXZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, linearThetaXZprime);
	recoLinear_thetaYZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, linearThetaYZprime);
	recoPolygonal_theta3DVsSegmentMomentum_HIST->Fill(segmentMomentum, polygonalTheta3D);
	recoPolygonal_thetaXZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, polygonalThetaXZprime);
	recoPolygonal_thetaYZprimeVsSegmentMomentum_HIST->Fill(segmentMomentum, polygonalThetaYZprime);

	//Fill segment momentum vs number of segment points
        reco_numberOfSegmentPointsVsSegmentMomentum_HIST->Fill(segmentMomentum, recoSegment.size());
      } // for i (segment angle)
      //std::cout << " Here 1 " << std::endl;
      // ========================================================================
      // Direct Angle Comparison and First Segment Analysis
      // ========================================================================

      // Calculate and plot the distance between vertices the true and reconstructed vertex.
      TVector3 initialPosition = startingPositionFor(particle); // The first point of the MCParticle that's within the detector.
      recob::tracking::Point_t reco_vertex = track->Vertex();
      Double_t vertexDifference = distanceBetween(initialPosition, reco_vertex);
      reco_vertexDifference_HIST->Fill(vertexDifference);

      // Only run for events with the same number of segments.
      // if(vertexDifference < 1.0) {
	// ========================================================================
	// Direct Comparison Analysis across many segments
	// Only possible if they're the same size (same # of segments)
	// ========================================================================
	if(trueLinearAngles_vec.size() == recoLinearAngles_vec.size()) {
	  // Loop through segments
	  for(size_t i = 0; i < recoLinearAngles_vec.size(); ++i) {

	    Angles_t trueLinearAngles = trueLinearAngles_vec.at(i);
	    Angles_t recoLinearAngles = recoLinearAngles_vec.at(i);
	    // TODO: Alternative Segment Defintions.
	    
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
        //std::cout << " Here 2 " << std::endl;
	// ========================================================================
	// First Segment Analysis
	// The first segment's linear-fit typically serves as the intial direction for the second segment.
	// The first angles in the (true/reco)LinearAngles_vecare calculated between the first and second segments.
	// We're about to calculate the angles of this segment relative to the initial direction.
	// We'll then optionally add the next few segment's angles, up to how many are stated in the "numberOf(True/Reco)Segments" variable.
	//
	// A note about the SegmentBBMomentum plots:
	// If numberOf(True/Reco)Segments == 0, these should be exactly the same as the regular SegmentMomentum plots.
	//
	// TODO: Alternative Segment Defintions.
	// ========================================================================

	// Initial Info (prior to first segment being formed)
	Double_t initialMomentum = particle.P();
	TVector3 initialDirection = startingDirectionFor(particle);
	
	// Get initial linear-fit vectors of first segment.

	// Calculate the TVector3 of the linear-fit of the first segment.
	TVector3 trueLinear_firstSegment = angleCalculator.convert<Vector_t, TVector3>(true_linearFit_vec.at(0));
	TVector3 recoLinear_firstSegment = angleCalculator.convert<Vector_t, TVector3>(reco_linearFit_vec.at(0));

	// Calculate angles of the first segment relative to the initial direction.
	Angles_t trueLinear_firstSegmentAngles = angleCalculator.calculateAngles(angleCalculator.convert<Vector_t, TVector3>(true_linearFit_vec.at(1)), trueLinear_firstSegment);
	Angles_t recoLinear_firstSegmentAngles = angleCalculator.calculateAngles(recoLinear_firstSegment, initialDirection);

	// Extract the individual true angles
	Float_t trueLinear_theta3D = std::get<0>(trueLinear_firstSegmentAngles)*1000.0;
	Float_t trueLinear_thetaXZprime = std::get<1>(trueLinear_firstSegmentAngles)*1000.0;
        Float_t trueLinear_thetaYZprime = std::get<2>(trueLinear_firstSegmentAngles)*1000.0;

	// Extract the individual reco angles
        Float_t recoLinear_theta3D = std::get<0>(recoLinear_firstSegmentAngles)*1000.0;
	Float_t recoLinear_thetaXZprime = std::get<1>(recoLinear_firstSegmentAngles)*1000.0;
        Float_t recoLinear_thetaYZprime = std::get<2>(recoLinear_firstSegmentAngles)*1000.0;

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

	// Loop through the next few segments and add to the first segment plots.
	
	// We'll begin by defining the number of segments to be all of the segments.
	size_t numberOfTrueSegmentsToIncludeInFirst = trueLinearAngles_vec.size();
	size_t numberOfRecoSegmentsToIncludeInFirst = recoLinearAngles_vec.size();

	// Optionally reset the number of segments to some value to just look at the first couple of segments.
	// Set to 0 to not add any additional segments to the first segment angles histograms
	// TODO: Consider making this a fhicl parameter.
	numberOfTrueSegmentsToIncludeInFirst = 0;
	numberOfRecoSegmentsToIncludeInFirst = 0;

	// Calculate the inital true momentum
	Double_t initialTrueMomentum = trueSegmentResult.GetMomentumForSegment(0);

	// The segment momentum at the start of the first segment.
	Double_t trueSegmentBBMomentum = initialTrueMomentum;

	// Perform BB verification and optionally add the first-few segment angles to the first segment angle graphs.
	for(size_t i = 0; i < true_rawSegmentLength_vec.size(); ++i) {
	  // TODO: Assumption: Relative Sizes of vectors. Should add about what we're looping through.

	  // Update the true segment BB momentum for this segment.
	  Double_t trueRawSegmentLength = true_rawSegmentLength_vec.at(i);
	  trueSegmentBBMomentum -= deltaP(trueSegmentBBMomentum*1000.0, trueRawSegmentLength)/1000.0;

	  // Calculate the true segment momentum for the start of this segment.
	  Double_t trueSegmentMomentum = trueSegmentResult.GetMomentumForSegment(i+1);

	  // Fill BB momentum verification histograms
	  trueBBMomentumVsTrueMomentum_HIST->Fill(trueSegmentMomentum, trueSegmentBBMomentum);

	  // Compare the straight-line distance between segment end points with this segment's raw segment length:
	  // Calculate the straight-line distance between the segment end points.
	  Segment_t trueSegment = true_segment_vec.at(i);
	  Point_t truePoint_beginning = trueSegment.front();
	  Point_t truePoint_end = trueSegment.back();
	  Double_t straightLineDistance = distanceBetween(truePoint_beginning, truePoint_end);

	  // Fill segment length comparison histograms
	  true_rawSegmentLengthVsStraightLineDistance_HIST->Fill(straightLineDistance, trueRawSegmentLength);
	  true_rawSegmentLengthVsSegmentMomentum_HIST->Fill(trueSegmentMomentum, trueRawSegmentLength);
	  //std::cout << " Here 3 " << std::endl;
	  // If we should add the first few segments to the first segment angles plots.
	  // numberOfTrueSegments will be 0 if not.
	  if(i > 0 && i < numberOfTrueSegmentsToIncludeInFirst+1) { // TODO: Assumption i-1 and numberOfTrueSegmentToIncludeInFirst+1 is to make sure the segment momentum is correct.

	    // Get the True Angles
	    // The i-1 is because of the relative sizes of the vectors and because we have the i > 0 condition.
	    // This was needed to the segment momentum is correct.
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
	for(size_t i = 0; i < rawSegmentLength_vec.size(); ++i) {

	  // TODO: Assumption: Relative vector sizes

	  // Calculate the true segment momentum for the true point that is closest in Z to the reco segment start point
	  // TODO: Assumption that this is the best way to have the "true momentum" for the start of the reco segment.  Can this be done using sim::IDE?
	  Segment_t recoSegment = reco_segment_vec.at(i+1);
	  Point_t recoPoint = recoSegment.at(0);
	  size_t truePointClosestZ_index = TrajectoryPointWithClosestZ(particle, recoPoint);
	  Double_t trueSegmentMomentumClosestZ = particle.P(truePointClosestZ_index);

	  // Update the reco segment BB momentum for the next segment.
	  Double_t recoRawSegmentLength = rawSegmentLength_vec.at(i);
	  recoSegmentBBMomentum -= deltaP(recoSegmentBBMomentum*1000.0, recoRawSegmentLength)/1000.0;

	  // Fill BB histograms
	  recoBBMomentumVsTrueMomentumClosestZ_HIST->Fill(trueSegmentMomentumClosestZ, recoSegmentBBMomentum);	  

	  // Compare the straight-line distance between adjacent segment start points with this segment's raw segment length:
	  // TODO: Assumption that the end of one segment is the start of the next segment.
	  Point_t recoPoint_beginning = recoSegment.front();
	  Point_t recoPoint_end = recoSegment.back();
	  Double_t straightLineDistance = distanceBetween(recoPoint_beginning, recoPoint_end);

	  // Fill segment length histogram
	  reco_rawSegmentLengthVsStraightLineDistance_HIST->Fill(straightLineDistance, recoRawSegmentLength);
	  reco_rawSegmentLengthVsSegmentMomentumClosestZ_HIST->Fill(trueSegmentMomentumClosestZ, recoRawSegmentLength);

	  // If we should add the first few segments to the first segment angles plots.
	  // numberOfRecoSegments will be 0 if not.
	  if(i > 0 && i < numberOfRecoSegmentsToIncludeInFirst+1) { // TODO: Assumption, same as true version
	    // Get the Reco Angles
	    Angles_t recoLinearAngles = recoLinearAngles_vec.at(i-1); // TODO: Assumption

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
	//std::cout << " Here 4 " << true_segment_vec.size() << " " << reco_segment_vec.size() << std::endl;	

	// Segment Analysis:
	// Tells us about our segment length measurements and the closest Z momentum ethod, which can be used to interpret the validity of the momentum estimation

	// Calculate the distance between segment start points for true and reco segment start points for tracks that have the same number of segments.
	for(size_t i = 0; i < true_segment_vec.size() && i < reco_segment_vec.size(); ++i) {
	  // Define this segment.
	  Segment_t trueSegment = true_segment_vec.at(i);
	  Segment_t recoSegment = reco_segment_vec.at(i);
        //std::cout << " Here 5 " << i  << std::endl;
	  // Define various points.
	  Point_t recoPoint = recoSegment.front();
	  Point_t truePoint = trueSegment.front();

	  size_t truePointWithClosestZToRecoPoint_index = TrajectoryPointWithClosestZ(particle, recoPoint);
        //std::cout << " Here 6 " << i  << std::endl;
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

	} // for i in trueSegmentStartPoints.size() && recoSegmentStartPoints.size()

	// TODO: Plot reco BB momentum vs. true momentum of true point with closest Z to the reco point
	// TODO: Plot true BB momentum vs. true momentum of true point

	// } // if vertexDifference < 1
      // else {
      // 	counter++;
      // }
    } // if length > 100
  nSegwithLessthan10pts_HIST->Fill(nSegwithLessthan10pts); //this is filled for every track
  } // for iTrack
  myfile.close();
} // analyze function

void MCSAngleAnalysis::beginJob() {
  art::ServiceHandle<art::TFileService> tfs; // Consider moving to private class scope, rather than function scope, if we do endJob analysis that needs to be saved (such as sigmaRES analysis).

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

  // reco tracklength
  recoNumberOfTracks_HIST=tfs->make<TH1F>("recoNumberOfTracks_HIST","Reco. Number of Tracks; Number of Tracks;Events",100,0,100);

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

  trueNumberOfSegmentPoints_HIST  = tfs->make<TH1F>("trueNumberOfSegmentPoints_HIST", "Number of points in a true segment;#Points in true segment;#segment", 500, 0, 1000);
  recoNumberOfSegmentPoints_HIST  = tfs->make<TH1F>("recoNumberOfSegmentPoints_HIST", "Number of points in a reco segment;#Points in reco segment;#segment", 100, 0, 200);

  true_numberOfSegmentPointsVsSegmentMomentum_HIST = tfs->make<TH2F>("true_numberOfSegmentPointsVsSegmentMomentum_HIST", "#trajectory points vs. Segment Momentum; Segment Momentum (GeV/c); #Trajectory points per segment (14cm)", 451, -0.05, 4.55, 100, 0, 200);
  reco_numberOfSegmentPointsVsSegmentMomentum_HIST = tfs->make<TH2F>("reco_numberOfSegmentPointsVsSegmentMomentum_HIST", "#trajectory points vs. Segment Momentum; Segment Momentum (GeV/c); #Trajectory points per segment (14cm)", 451, -0.05, 4.55, 100, 0, 200);
  alttrue_numberOfSegmentPointsVsSegmentMomentum_HIST = tfs->make<TH2F>("alttrue_numberOfSegmentPointsVsSegmentMomentum_HIST", "#trajectory points vs. Segment Momentum; Segment Momentum (GeV/c); #Trajectory points per segment (14cm)", 451, -0.05, 4.55, 100, 0, 200);

  recoSegmentPointsXZ_HIST = tfs->make<TH2F>("recoSegmentPointsXZ_HIST", "X vs Z; Z; X", 610,-10,600,600,-300,300);
  recoSegmentPointsYZ_HIST = tfs->make<TH2F>("recoSegmentPointsYZ_HIST", "Y vs Z; Z; Y", 610,-10,600,600,100,500);

  trueSegmentPointsXZ_HIST = tfs->make<TH2F>("trueSegmentPointsXZ_HIST", "X vs Z; Z; X", 610,-10,600,600,-300,300);
  trueSegmentPointsYZ_HIST = tfs->make<TH2F>("trueSegmentPointsYZ_HIST", "Y vs Z; Z; Y", 610,-10,600,600,100,500);

  alttrueSegmentPointsXZ_HIST = tfs->make<TH2F>("alttrueSegmentPointsXZ_HIST", "X vs Z; Z; X", 610,-10,600,600,-300,300);
  alttrueSegmentPointsYZ_HIST = tfs->make<TH2F>("alttrueSegmentPointsYZ_HIST", "Y vs Z; Z; Y", 610,-10,600,600,100,500);

  alttrueSegmentPointsXZ_lt20PtsperSeg_HIST = tfs->make<TH2F>("alttrueSegmentPointsXZ_lt20PtsperSeg_HIST", "X vs Z; Z; X", 1000,-10,600,1000,-300,300);
  alttrueSegmentPointsYZ_lt20PtsperSeg_HIST = tfs->make<TH2F>("alttrueSegmentPointsYZ_lt20PtsperSeg_HIST", "Y vs Z; Z; Y", 1000,-10,600,1000,100,500);

  trueSimIDEdistance_HIST = tfs->make<TH1F>("trueSimIDEdistance_HIST", "Distance between simIDE;Distance between simIDE (cm);", 100, 0, 100);
  simIDEtrackID_HIST = tfs->make<TH1F>("simIDEtrackID_HIST", "Track ID by simIDE;TrackID;", 10000, -5000, 5000);
  simIDEedep_HIST = tfs->make<TH1F>("simIDEedep_HIST", "Edep by simIDE;Energy deposited [MeV];", 1000, 0, 10);
  simIDEedepFromStart_HIST = tfs->make<TH2F>("simIDEedepFromStart_HIST", "Edep by simIDE;Distance from start [cm];Energy deposited [MeV]", 250,0,600,1000, 0, 10);
  deltaXZ_HIST = tfs->make<TH1F>("deltaXZ_HIST", ";#Delta_{XZ};", 200, -100, 100);
  deltaYZ_HIST = tfs->make<TH1F>("deltaYZ_HIST", ";#Delta_{YZ};", 200, -100, 100);

  trueLinear_deltaXZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("trueLinear_deltaXZprimeVsSegmentMomentum_HIST","True #Delta_{XZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #delta_{XZ}' (cm)", 500, -0.05, 5.05, 100, -50, 50);

  trueLinear_deltaYZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("trueLinear_deltaYZprimeVsSegmentMomentum_HIST","True #Delta_{YZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #delta_{YZ}' (cm)", 500, -0.05, 5.05, 100, -50, 50);

  trueLinear_transverseDistX_HIST = tfs->make<TH1F>("trueLinear_transverseDistX_HIST", ";Transverse Distance X;", 200, -5, 5);
  recoLinear_transverseDistX_HIST = tfs->make<TH1F>("recoLinear_transverseDistX_HIST", ";Transverse Distance X;", 200, -5, 5);
  trueLinear_transverseDistY_HIST = tfs->make<TH1F>("trueLinear_transverseDistY_HIST", ";Transverse Distance Y;", 200, -5, 5);
  recoLinear_transverseDistY_HIST = tfs->make<TH1F>("recoLinear_transverseDistY_HIST", ";Transverse Distance Y;", 200, -5, 5);

  trueLinear_transverseDistXVsSegmentMomentum_HIST = tfs->make<TH2F>("trueLinear_transverseDistXVsSegmentMomentum_HIST","True Transverse Distance in X Vs Segment Momentum; Segment Momentum (GeV/c); Transverse Distance in X (cm)", 500, -0.05, 5.05, 200, -5, 5);
  trueLinear_transverseDistYVsSegmentMomentum_HIST = tfs->make<TH2F>("trueLinear_transverseDistYVsSegmentMomentum_HIST","True Transverse Distance in Y Vs Segment Momentum; Segment Momentum (GeV/c); Transverse Distance in Y (cm)", 500, -0.05, 5.05, 200, -5, 5);

  recoLinear_transverseDistXVsSegmentMomentum_HIST = tfs->make<TH2F>("recoLinear_transverseDistXVsSegmentMomentum_HIST","Reco Transverse Distance in X Vs Segment Momentum; Segment Momentum (GeV/c); Transverse Distance in X (cm)", 500, -0.05, 5.05, 200, -5, 5);
  recoLinear_transverseDistYVsSegmentMomentum_HIST = tfs->make<TH2F>("recoLinear_transverseDistYVsSegmentMomentum_HIST","Reco Transverse Distance in Y Vs Segment Momentum; Segment Momentum (GeV/c); Transverse Distance in Y (cm)", 500, -0.05, 5.05, 200, -5, 5);

  altTrueLinear_ThetaXZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("altTrueLinear_ThetaXZprimeVsSegmentMomentum_HIST", "Alt. True Linear #theta_{XZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{XZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);
  altTrueLinear_ThetaYZprimeVsSegmentMomentum_HIST = tfs->make<TH2F>("altTrueLinear_ThetaYZprimeVsSegmentMomentum_HIST", "Alt. True Linear #theta_{YZ}' Vs Segment Momentum; Segment Momentum (GeV/c); #theta_{YZ}' (mrad)", 500, -0.05, 5.05, 401, -200.5, 200.5);

  altTrueNumberOfSegmentPoints_HIST  = tfs->make<TH1F>("altTrueNumberOfSegmentPoints_HIST", "Number of points in a true matched segment;#Points in true matched segment;#segment", 100, 0, 200);
  truePtsRecoMatched_distance_HIST = tfs->make<TH1F>("truePtsRecoMatched_distance_HIST", "Distance between True Pts Matched;Distance between True Pts Matched (cm);", 100, 0, 10);
  recoPts_distance_HIST = tfs->make<TH1F>("recoPts_distance_HIST", "Distance between Reco Pts;Distance between Reco Pts (cm);", 100, 0, 10);

  trueMatchedReco_deltaX_HIST = tfs->make<TH1F>("trueMatchedReco_deltaX_HIST", "True matched - Reco X;#DeltaX;", 1000, -50, 50);
  trueMatchedReco_deltaY_HIST = tfs->make<TH1F>("trueMatchedReco_deltaY_HIST", "True matched - Reco Y;#DeltaY;", 1000, -50, 50);
  trueMatchedReco_deltaZ_HIST = tfs->make<TH1F>("trueMatchedReco_deltaZ_HIST", "True matched - Reco Z;#DeltaZ;", 1000, -50, 50);

  reco_thetaXZVsNpts_HIST = tfs->make<TH2F>("reco_thetaXZVsNpts_HIST","Reco Linear #theta_{XZ}' Vs #Points per segment; Number of Points per segment; #theta_{XZ}' (mrad)", 100, 0, 100, 401, -200.5, 200.5);
  reco_thetaYZVsNpts_HIST = tfs->make<TH2F>("reco_thetaYZVsNpts_HIST","Reco Linear #theta_{YZ}' Vs #Points per segment; Number of Points per segment; #theta_{YZ}' (mrad)", 100, 0, 100, 401, -200.5, 200.5);

  alttrue_thetaXZVsNpts_HIST = tfs->make<TH2F>("alttrue_thetaXZVsNpts_HIST","True Linear #theta_{XZ}' Vs #Points per segment; Number of Points per segment; #theta_{XZ}' (mrad)", 100, 0, 100, 401, -200.5, 200.5);
  alttrue_thetaYZVsNpts_HIST = tfs->make<TH2F>("alttrue_thetaYZVsNpts_HIST","True Linear #theta_{YZ}' Vs #Points per segment; Number of Points per segment; #theta_{YZ}' (mrad)", 100, 0, 100, 401, -200.5, 200.5);

  reco_vs_trueX_HIST = tfs->make<TH2F>("reco_vs_trueX_HIST", "True X vs Reco X; True X; Reco X", 600,-300,300,600,-300,300);
  reco_vs_trueY_HIST = tfs->make<TH2F>("reco_vs_trueY_HIST", "True Y vs Reco Z; True Y; Reco Y", 600,100,500,600,100,500);
  reco_vs_trueZ_HIST = tfs->make<TH2F>("reco_vs_trueZ_HIST", "True Z vs Reco Z; True Z; Reco Z", 610,-10,600,610,-10,600);

  alttrue_distFromStartVsdistBtwPts_HIST = tfs->make<TH2F>("alttrue_distFromStartVsdistBtwPts_HIST", "Dist from start vs dist btw pts; Distance from start [cm]; Distance between Pts w.r.t. prev [cm]", 6500,0,650,200, 0, 20);
  reco_distFromStartVsdistBtwPts_HIST = tfs->make<TH2F>("reco_distFromStartVsdistBtwPts_HIST", "Dist from start vs dist btw pts; Distance from start [cm]; Distance between Pts w.r.t. prev [cm]", 6500,0,650,200, 0, 20);

  alttrue_pointIndexVsdistBtwPts_HIST = tfs->make<TH2F>("alttrue_pointIndexVsdistBtwPts_HIST", "Dist from start vs dist btw pts; Distance from start [cm]; Distance between Pts w.r.t. prev [cm]", 5000,0,5000,200, 0, 20);
  reco_pointIndexVsdistBtwPts_HIST = tfs->make<TH2F>("reco_pointIndexVsdistBtwPts_HIST", "Dist from start vs dist btw pts; Distance from start [cm]; Distance between Pts w.r.t. prev [cm]", 5000,0,5000,200, 0, 20);

  distFromStartVsDeltaX_HIST = tfs->make<TH2F>("distFromStartVsDeltaX_HIST", "Point Index from Start vs (true-reco) X; Point Index; #Delta X", 5000,0,5000, 1000, -50, 50);
  distFromStartVsDeltaY_HIST = tfs->make<TH2F>("distFromStartVsDeltaY_HIST", "Point Index from Start vs (true-reco) Y; Point Index; #Delta Y", 5000,0,5000, 1000, -50, 50);
  distFromStartVsDeltaZ_HIST = tfs->make<TH2F>("distFromStartVsDeltaZ_HIST", "Point Index from Start vs (true-reco) Z; Point Index; #Delta Z", 5000,0,5000, 1000, -50, 50);

  distFromStartIfDistAbove3VsDeltaX_HIST = tfs->make<TH2F>("distFromStartIfDistAbove3VsDeltaX_HIST", "Point Index from Start vs (true-reco) X; Point Index; #Delta X", 5000,0,5000, 1000, -50, 50);
  distFromStartIfDistAbove3VsDeltaY_HIST = tfs->make<TH2F>("distFromStartIfDistAbove3VsDeltaY_HIST", "Point Index from Start vs (true-reco) Y; Point Index; #Delta Y", 5000,0,5000, 1000, -50, 50);
  distFromStartIfDistAbove3VsDeltaZ_HIST = tfs->make<TH2F>("distFromStartIfDistAbove3VsDeltaZ_HIST", "Point Index from Start vs (true-reco) Z; Point Index; #Delta Z", 5000,0,5000, 1000, -50, 50);

  nSegwithLessthan10pts_HIST = tfs->make<TH1F>("nSegwithLessthan10pts_HIST", "#seg per track with <10 points;#seg per track;#tracks", 100, 0, 50);
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
/*
template<typename Point>
Vector_t MCSAngleAnalysis::Calculate2DLinearFit(const std::vector<Point> segment) const {

  // Calculate the barycenter of this segment.
  Barycenter_t barycenter = CalculateBarycenter(segment);

  // Declare initial sums as 0. (used in linear-regression)
  Float_t SSxz = 0, SSzz = 0;

  // Loop through the points in this segment.
  for(size_t k = 0; k < segment.size(); ++k) {
    // Calculate values used in sum calculations.
    Point_t trajPoint = segment.at(k);
    Double_t xDiff = trajPoint.X() - barycenter.X();
    Double_t zDiff = trajPoint.Z() - barycenter.Z();

    SSxz += xDiff * zDiff;
    SSzz += zDiff * zDiff;
  } // for k

  // Slopes of the linear fit, projected onto two planes
  Double_t A = SSxz / SSzz; // Source: x = Az + B
  Double_t B = barycenter.X() - A * barycenter.Z();

  // Calculate the components of the linear fit.
  Double_t x = A;
  Double_t y = 0;
  Double_t z = B;

  // Calculate linear fit vector, transform to unit length.
  //Vector_t linearFit = Vector_t(x, y, z);

  // Return the linear fit vector.
  return Vector_t(x, y, z);
} // MCSSegmentCalculator::CalculateLinearFit
*/
DEFINE_ART_MODULE(MCSAngleAnalysis)
