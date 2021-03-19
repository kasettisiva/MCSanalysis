#ifndef NewProtoDUNECalibration_h
#define NewProtoDUNECalibration_h

#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Track.h"
#include "art/Framework/Principal/Event.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"


namespace protoana{


  class ProtoDUNECalibration{

    public:
      ProtoDUNECalibration(){};
      ProtoDUNECalibration(const fhicl::ParameterSet & pset);
      std::vector<float> GetCalibratedCalorimetry(
          const recob::Track &track, art::Event const &evt,
          const std::string trackModule, const std::string caloModule,
          size_t planeID, double negativeZFix = 0.);
      double HitToEnergy(
          const art::Ptr<recob::Hit> hit, double X, double Y, double Z,
          double recomb_factor=.6417);
      std::vector<double> GetEFieldVector(
          const recob::Track &track, art::Event const &evt,
          const std::string trackModule, const std::string caloModule,
          size_t planeID, double negativeZFix = 0.);
      std::vector<double> CalibratedQdX(
          const recob::Track &track, art::Event const &evt,
          const std::string trackModule, const std::string caloModule,
          size_t planeID, double negativeZFix);
      float calc_dEdX(double dqdx, double betap, double Rho, double Efield, double Wion, double alpha);

    private:

      double tot_Ef( double, double, double );

      double betap;
      double Rho;
      //double Efield;
      double Wion;
      double alpha;

      //size_t planeID;
      std::map<size_t, double> norm_factors;
      std::map<size_t, double> calib_factors;
      std::map<size_t, TH1F *> X_correction_hists;
      std::map<size_t, TH2F *> YZ_neg_hists;
      std::map<size_t, TH2F *> YZ_pos_hists;

      std::string X_correction_name;
      TFile * X_correction_file;

      std::string YZ_correction_name;
      TFile * YZ_correction_file;

      std::string E_field_correction_name;
      TFile * E_field_file;


      TH3F * ex_neg;
      TH3F * ey_neg;
      TH3F * ez_neg;

      TH3F * ex_pos;
      TH3F * ey_pos;
      TH3F * ez_pos;

      ProtoDUNETrackUtils trackUtil;

      TFile * OpenFile(const std::string filename);

  };

}


#endif
