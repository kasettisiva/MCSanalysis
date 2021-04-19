#include "ProtoDUNECalibration.h"
#include "cetlib/search_path.h"
#include "cetlib/filesystem.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


protoana::ProtoDUNECalibration::ProtoDUNECalibration(const fhicl::ParameterSet & pset) 
    : //planeID(pset.get<unsigned int>("PlaneID")),
      betap(pset.get<double>("betap")),
      Rho(pset.get<double>("Rho")),
      Wion(pset.get<double>("Wion")),
      alpha(pset.get<double>("alpha")),
      //norm_factor(pset.get<double>("norm_factor")),
      //calib_factor(pset.get<double>("calib_factor")),
      X_correction_name(pset.get<std::string>("X_correction")),
      YZ_correction_name(pset.get<std::string>("YZ_correction")),
      E_field_correction_name(pset.get<std::string>("E_field_correction")) {

  
  /*
  X_correction_file = new TFile(X_correction_name.c_str(), "OPEN");
  YZ_correction_file = new TFile(YZ_correction_name.c_str(), "OPEN");
  E_field_file = new TFile(E_field_correction_name.c_str(), "OPEN");
  */
  X_correction_file = OpenFile(X_correction_name);
  YZ_correction_file = OpenFile(YZ_correction_name);
  E_field_file = OpenFile(E_field_correction_name);

  std::vector<fhicl::ParameterSet> PlaneParameters =
      pset.get<std::vector<fhicl::ParameterSet>>("PlaneParameters");

  for (size_t i = 0; i < PlaneParameters.size(); ++i) {
    size_t planeID = PlaneParameters[i].get<size_t>("PlaneID");
    norm_factors[planeID] = PlaneParameters[i].get<double>("norm_factor");
    calib_factors[planeID] = PlaneParameters[i].get<double>("calib_factor");

    std::string hist_name = "dqdx_X_correction_hist_" + std::to_string(planeID);
    X_correction_hists[planeID] = (TH1F*)X_correction_file->Get( hist_name.c_str() );

    std::string YZ_neg_name = "correction_dqdx_ZvsY_negativeX_hist_" +
                              std::to_string(planeID);
    YZ_neg_hists[planeID] = (TH2F*)YZ_correction_file->Get(YZ_neg_name.c_str());

    std::string YZ_pos_name = "correction_dqdx_ZvsY_positiveX_hist_" +
                              std::to_string(planeID);
    YZ_pos_hists[planeID] = (TH2F*)YZ_correction_file->Get(YZ_pos_name.c_str());
  }

  ex_neg = (TH3F*)E_field_file->Get("Reco_ElecField_X_Neg");
  ey_neg = (TH3F*)E_field_file->Get("Reco_ElecField_Y_Neg");
  ez_neg = (TH3F*)E_field_file->Get("Reco_ElecField_Z_Neg");
  ex_pos = (TH3F*)E_field_file->Get("Reco_ElecField_X_Pos");
  ey_pos = (TH3F*)E_field_file->Get("Reco_ElecField_Y_Pos");
  ez_pos = (TH3F*)E_field_file->Get("Reco_ElecField_Z_Pos");
 

  /*
  std::cout << "Calibration" << std::endl;
  std::cout << planeID << std::endl;
  std::cout << betap << std::endl;
  std::cout << Rho << std::endl;
  std::cout << Wion << std::endl;
  std::cout << norm_factor << std::endl;
  std::cout << calib_factor << std::endl;
  */
}

std::vector<float> protoana::ProtoDUNECalibration::GetCalibratedCalorimetry(
    const recob::Track &track, art::Event const &evt,
    const std::string trackModule, const std::string caloModule,
    size_t planeID, double negativeZFix) {

  std::vector<float> calibrated_dEdx;
  
  //If we can't find the plane ID in the configuration, return empty vector
  if (norm_factors.find(planeID) == norm_factors.end())
    return calibrated_dEdx;

  //Get the Calorimetry vector from the track
  std::vector< anab::Calorimetry > caloVector = trackUtil.GetRecoTrackCalorimetry( track, evt, trackModule, caloModule ); 
  
  size_t calo_position;
  bool found_plane = false;
  for( size_t i = 0; i < caloVector.size(); ++i ){
     unsigned int thePlane = caloVector.at(i).PlaneID().Plane;
     if( thePlane == planeID ){
       calo_position = i;
       found_plane = true;
       break;
     }
  }

  if( !found_plane ){
    std::cout << "Could not find the correct plane in the calorimetry vector" << std::endl;
    return calibrated_dEdx;
  }

  std::vector< float > dQdX = caloVector.at( calo_position).dQdx();
  auto theXYZPoints = caloVector.at( calo_position).XYZ();
  std::vector< float > resRange = caloVector.at( calo_position ).ResidualRange();

  //Get the hits from the track from a specific plane
  const std::vector< const recob::Hit* > hits = trackUtil.GetRecoTrackHitsFromPlane( track, evt, trackModule, planeID ); 
  if( hits.size() == 0 ){
    std::cout << "Got empty hits vector" << std::endl;
    return calibrated_dEdx;
  }

  if (negativeZFix > 0.) {
    return calibrated_dEdx;
  }

  double z_check = negativeZFix;

  //Do Ajib's correction 
  for( size_t i = 0; i < dQdX.size(); ++i ){ 
    float hit_x = theXYZPoints[i].X();
    float hit_y = theXYZPoints[i].Y();
    float hit_z = theXYZPoints[i].Z();

    //if( hit_y < 0. || hit_y > 600. ) continue;
    //if( hit_z < z_check || hit_z > 695. ) continue;
    if (hit_y < 0. || hit_y > 600. || hit_z < z_check || hit_z > 695.) {
      calibrated_dEdx.push_back(-999.);
      continue; 
    }

    //Set the z position to 0. for small (configurable) negative positions
    if (negativeZFix < hit_z && hit_z < 0.) {
      //std::cout << "Fixing: " << hit_z << " " << negativeZFix << std::endl;
      hit_z = 0.;
    }


    int X_bin = X_correction_hists[planeID]->FindBin(hit_x);
    float X_correction = X_correction_hists[planeID]->GetBinContent(X_bin);

    TH2F * YZ_hist = (hit_x < 0 ? 
                      YZ_neg_hists[planeID] :
                      YZ_pos_hists[planeID]);
    int YZ_bin = YZ_hist->FindBin(hit_z, hit_y);
    double YZ_correction = YZ_hist->GetBinContent(YZ_bin);


    float corrected_dq_dx = dQdX[i] * X_correction *
                            YZ_correction * norm_factors[planeID];
    float scaled_corrected_dq_dx = corrected_dq_dx / calib_factors[planeID];

    double Efield = tot_Ef( hit_x, hit_y, hit_z );


    float cal_de_dx = calc_dEdX( scaled_corrected_dq_dx,  betap,  Rho,  Efield,  Wion,  alpha );
 
    calibrated_dEdx.push_back( cal_de_dx );
  }


  return calibrated_dEdx;
}

std::vector<double> protoana::ProtoDUNECalibration::CalibratedQdX(
    const recob::Track &track, art::Event const &evt,
    const std::string trackModule, const std::string caloModule,
    size_t planeID, double negativeZFix) {
  std::vector<double> calibrated_dQdX;
  
  //If we can't find the plane ID in the configuration, return empty vector
  if (norm_factors.find(planeID) == norm_factors.end())
    return calibrated_dQdX;

  //Get the Calorimetry vector from the track
  std::vector< anab::Calorimetry > caloVector = trackUtil.GetRecoTrackCalorimetry( track, evt, trackModule, caloModule ); 
  
  size_t calo_position;
  bool found_plane = false;
  for( size_t i = 0; i < caloVector.size(); ++i ){
     unsigned int thePlane = caloVector.at(i).PlaneID().Plane;
     if( thePlane == planeID ){
       calo_position = i;
       found_plane = true;
       break;
     }
  }

  if( !found_plane ){
    std::cout << "Could not find the correct plane in the calorimetry vector" << std::endl;
    return calibrated_dQdX;
  }

  std::vector< float > dQdX = caloVector.at( calo_position).dQdx();
  auto theXYZPoints = caloVector.at( calo_position).XYZ();

  //Get the hits from the track from a specific plane
  const std::vector< const recob::Hit* > hits = trackUtil.GetRecoTrackHitsFromPlane( track, evt, trackModule, planeID ); 
  if( hits.size() == 0 ){
    std::cout << "Got empty hits vector" << std::endl;
    return calibrated_dQdX;
  }

  if (negativeZFix > 0.) {
    return calibrated_dQdX;
  }

  double z_check = negativeZFix;

  //Do Ajib's correction 
  for( size_t i = 0; i < dQdX.size(); ++i ){ 
    double hit_x = theXYZPoints[i].X();
    double hit_y = theXYZPoints[i].Y();
    double hit_z = theXYZPoints[i].Z();

    //if( hit_y < 0. || hit_y > 600. ) continue;
    //if( hit_z < z_check || hit_z > 695. ) continue;
    if (hit_y < 0. || hit_y > 600. || hit_z < z_check || hit_z > 695.) {
      calibrated_dQdX.push_back(-999.);
      continue; 
    }

    //Set the z position to 0. for small (configurable) negative positions
    if (negativeZFix < hit_z && hit_z < 0.) {
      hit_z = 0.;
    }


    int X_bin = X_correction_hists[planeID]->FindBin(hit_x);
    double X_correction = X_correction_hists[planeID]->GetBinContent(X_bin);

    TH2F * YZ_hist = (hit_x < 0 ? 
                      YZ_neg_hists[planeID] :
                      YZ_pos_hists[planeID]);
    int YZ_bin = YZ_hist->FindBin(hit_z, hit_y);
    double YZ_correction = YZ_hist->GetBinContent(YZ_bin);


    double corrected_dq_dx = dQdX[i] * X_correction *
                            YZ_correction * norm_factors[planeID];
    double scaled_corrected_dq_dx = corrected_dq_dx / calib_factors[planeID];

    calibrated_dQdX.push_back(scaled_corrected_dq_dx);
  }
  return calibrated_dQdX;
}

float protoana::ProtoDUNECalibration::calc_dEdX(double dqdx, double betap, double Rho, double Efield, double Wion, double alpha){
  return ( exp( dqdx * ( betap / ( Rho * Efield ) * Wion ) ) -alpha ) / ( betap / ( Rho*Efield ) );  
}

double protoana::ProtoDUNECalibration::tot_Ef( double x, double y, double z ){

  if( x >= 0 ){
    double ex = 0.5 + 0.5 * ex_pos->GetBinContent( ex_pos->FindBin( x, y, z ) );
    double ey = 0.5 * ey_pos->GetBinContent( ey_pos->FindBin( x, y, z ) );
    double ez = 0.5 * ez_pos->GetBinContent( ez_pos->FindBin( x, y, z ) );
    return sqrt( (ex*ex) + (ey*ey) + (ez*ez) );
  }
  else if( x < 0 ){
    double ex= 0.5 + 0.5 * ex_neg->GetBinContent( ex_neg->FindBin( x, y, z ) );
    double ey= 0.5 * ey_neg->GetBinContent( ey_neg->FindBin( x, y, z ) );
    double ez= 0.5 * ez_neg->GetBinContent( ez_neg->FindBin( x, y, z ) );
    return sqrt( (ex*ex) + (ey*ey) + (ez*ez) );
  }
  else return 0.5;
}

std::vector<double> protoana::ProtoDUNECalibration::GetEFieldVector(
    const recob::Track &track, art::Event const &evt,
    const std::string trackModule, const std::string caloModule,
    size_t planeID, double negativeZFix) {
  std::vector<double> results;

  //Get the Calorimetry vector from the track
  std::vector< anab::Calorimetry > caloVector = trackUtil.GetRecoTrackCalorimetry( track, evt, trackModule, caloModule ); 
  
  size_t calo_position;
  bool found_plane = false;
  for( size_t i = 0; i < caloVector.size(); ++i ){
     unsigned int thePlane = caloVector.at(i).PlaneID().Plane;
     if( thePlane == planeID ){
       calo_position = i;
       found_plane = true;
       break;
     }
  }

  if( !found_plane ){
    std::cout << "Could not find the correct plane in the calorimetry vector" << std::endl;
    return results;
  }

  auto theXYZPoints = caloVector.at( calo_position).XYZ();

  //Get the hits from the track from a specific plane
  const std::vector< const recob::Hit* > hits = trackUtil.GetRecoTrackHitsFromPlane( track, evt, trackModule, planeID ); 
  if( hits.size() == 0 ){
    std::cout << "Got empty hits vector" << std::endl;
    return results;
  }

  if (negativeZFix > 0.) {
    return results;
  }

  double z_check = negativeZFix;

  for( size_t i = 0; i < theXYZPoints.size(); ++i ){ 
    float hit_x = theXYZPoints[i].X();
    float hit_y = theXYZPoints[i].Y();
    float hit_z = theXYZPoints[i].Z();

    if( hit_y < 0. || hit_y > 600. || hit_z < z_check || hit_z > 695.) {
      results.push_back(-999.);
      continue;
    }

    //Set the z position to 0. for small (configurable) negative positions
    if (negativeZFix < hit_z && hit_z < 0.) {
      //std::cout << "Fixing: " << hit_z << " " << negativeZFix << std::endl;
      hit_z = 0.;
    }
    results.push_back(tot_Ef( hit_x, hit_y, hit_z ));
  }

  return results;
}

double protoana::ProtoDUNECalibration::HitToEnergy(
    const art::Ptr<recob::Hit> hit, double X, double Y, double Z,
    double recomb_factor) {

  //Only do collection plane
  //if( hit->View() != 2 ) return 0.;
  size_t planeID = hit->View();
  if (norm_factors.find(planeID) == norm_factors.end())
    return 0.;

  int X_bin = X_correction_hists[planeID]->FindBin(X);
  double X_factor = X_correction_hists[planeID]->GetBinContent(X_bin);
  double YZ_factor = 1.;
  if(X < 0.){
    int YZ_bin = YZ_neg_hists[planeID]->FindBin(Z, Y);
    YZ_factor = YZ_neg_hists[planeID]->GetBinContent(YZ_bin);
  }
  else{
    int YZ_bin = YZ_pos_hists[planeID]->FindBin(Z, Y);
    YZ_factor = YZ_pos_hists[planeID]->GetBinContent(YZ_bin);
  }

  double energy = hit->Integral();
  energy *= norm_factors[planeID];
  energy *= Wion/*23.6e-6*/;
  energy /= calib_factors[planeID];
  energy *= X_factor;
  energy *= YZ_factor;
  energy /= recomb_factor;
  
  return energy;

}

TFile * protoana::ProtoDUNECalibration::OpenFile(const std::string filename) {
  TFile * theFile = 0x0;
  mf::LogInfo("protoana::ProtoDUNECalibration::OpenFile") << "Searching for " << filename;
  if (cet::file_exists(filename)) {
    mf::LogInfo("protoana::ProtoDUNECalibration::OpenFile") << "File exists. Opening " << filename;
    theFile = new TFile(filename.c_str());
    if (!theFile ||theFile->IsZombie() || !theFile->IsOpen()) {
      delete theFile;
      theFile = 0x0;
      throw cet::exception("ProtoDUNECalibration.cxx") << "Could not open " << filename;
    }
  }
  else {
    mf::LogInfo("protoana::ProtoDUNECalibration::OpenFile") << "File does not exist here. Searching FW_SEARCH_PATH";
    cet::search_path sp{"FW_SEARCH_PATH"};
    std::string found_filename;
    auto found = sp.find_file(filename, found_filename);
    if (!found) {
      throw cet::exception("ProtoDUNECalibration.cxx") << "Could not find " << filename;
    }

    mf::LogInfo("protoana::ProtoDUNECalibration::OpenFile") << "Found file " << found_filename;
    theFile = new TFile(found_filename.c_str());
    if (!theFile ||theFile->IsZombie() || !theFile->IsOpen()) {
      delete theFile;
      theFile = 0x0;
      throw cet::exception("ProtoDUNECalibration.cxx") << "Could not open " << found_filename;
    }
  }
  return theFile;
}
