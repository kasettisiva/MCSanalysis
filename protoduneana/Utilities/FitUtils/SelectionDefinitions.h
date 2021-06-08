class new_interaction_topology {
 private:
  double fEndZLow, fEndZHigh;
  double fThreshold;
 public: 
  new_interaction_topology(double endz_low, double endz_high,
                           double threshold)
    : fEndZLow(endz_low), fEndZHigh(endz_high), fThreshold(threshold) {}

  int operator()(int pdg, double endZ,
                 std::string process, int nPi0,
                 std::vector<int> & true_daughter_pdg,
                 std::vector<double> & true_daughter_startP/*,
                 std::vector<double> & incidentEnergies*/) {

    int topology = -1;
    if (pdg == 211) {
      if (endZ < fEndZLow/*-.49375*/) {
        topology = 4;
      }
      //After FV
      else if (endZ > fEndZHigh/*222.10561*/) {
        topology = 6;
      }
      else if (process == "pi+Inelastic") {
        //daughters with & without thresholds
        bool has_pion_above_threshold = false;
        for (size_t i = 0; i < true_daughter_startP.size(); ++i) {
          if (abs(true_daughter_pdg[i]) == 211 &&
              true_daughter_startP[i] > fThreshold/*.150*/) {
            has_pion_above_threshold = true;
            break;
          }
        }

        if (has_pion_above_threshold) {
          topology = 3;
        }
        else if (nPi0 == 0) {
          topology = 1;
        }
        else if (nPi0 > 0) {
          topology = 2;
        }
      }
      else {
        topology = 7;
      }
    }
    else if (pdg == -13) {
      topology = 5;
    }
    else {
      topology = 7;
    }

    return topology;

  }
};

auto selection_ID = [](bool beam_is_track, bool ends_in_APA3,
                       bool no_pion_daughter,
                       bool beam_cuts, bool has_shower) {
  if (!beam_is_track) {
    return 6;
  }

  if (!beam_cuts) {
    return 5;
  }
  
  if (!ends_in_APA3) {
    return 4;
  }

  if (!no_pion_daughter) {
    return 3;
  }

  if (has_shower) {
    return 2;
  }
  else {
    return 1;
  }
};

auto daughter_PDG_types(const std::vector<int> bt_PDGs) {
  std::vector<int> results;
  for (const int PDG : bt_PDGs) {
    if (abs(PDG) == 211) {
      results.push_back(1);
    }
    else if (abs(PDG) == 13) {
      results.push_back(2);
    }
    else if (abs(PDG) == 2212) {
      results.push_back(3);
    }
    else if (abs(PDG) == 22) {
      results.push_back(4);
    }
    else if (abs(PDG) > 2212) {
      results.push_back(5);
    }
    else if (abs(PDG) == 11) {
      results.push_back(6);
    }
    else {
      results.push_back(7);
    }
  }
  return results;
}

auto categorize_daughters = [](
    const int beam_ID,
    const std::vector<int> bt_origins, const std::vector<int> bt_IDs,
    const std::vector<int> bt_PDGs, const std::vector<int> bt_ParIDs,
    const std::vector<int> bt_ParPDGs,
    const std::vector<int> true_daughters,
    const std::vector<int> true_grand_daughters) {
  std::vector<int> results;
  for (size_t i = 0; i < bt_origins.size(); ++i) {
    if (bt_IDs[i] == beam_ID) {
      results.push_back(1);
    }
    else if (bt_origins[i] == 2) {
      results.push_back(2);
    }
    else if (abs(bt_PDGs[i]) == 11 && (abs(bt_ParPDGs[i]) == 13)) {
      results.push_back(10); 
    }
    else if (std::find(true_daughters.begin(), true_daughters.end(), bt_IDs[i]) !=
             true_daughters.end()) {
      if (abs(bt_PDGs[i]) == 211) {
        results.push_back(3);
      }
      else if (abs(bt_PDGs[i]) == 13) {
        results.push_back(4);
      }
      else if (bt_PDGs[i] == 2212) {
        results.push_back(5);
      }
      else if (bt_PDGs[i] == 22) {
        results.push_back(6);
      }
      else if (bt_PDGs[i] > 2212) {
        results.push_back(7);
      }
      else {
        //std::cout << bt_PDGs[i] << " " << bt_ParPDGs[i] << " " << (abs(bt_PDGs[i]) == 11 && (abs(bt_ParPDGs[i]) == 13)) << std::endl;
        results.push_back(12);
      }
    }
    else if (std::find(true_grand_daughters.begin(), 
                       true_grand_daughters.end(), bt_IDs[i]) !=
             true_grand_daughters.end()) {
      if ((bt_PDGs[i] == 22 || abs(bt_PDGs[i]) == 11) && bt_ParPDGs[i] == 111) {
        results.push_back(11);
      }
      else {
        results.push_back(8);
      }
    }
    else if (std::find(true_grand_daughters.begin(), 
                       true_grand_daughters.end(), bt_ParIDs[i]) !=
             true_grand_daughters.end()) {
        results.push_back(9);
    }
    else {
        //std::cout << bt_PDGs[i] << " " << bt_ParPDGs[i] << " " << (abs(bt_PDGs[i]) == 11 && (abs(bt_ParPDGs[i]) == 13)) << std::endl;
        results.push_back(12);
    }
  }
  return results;
};

auto backtrack_beam = [](const std::string process,
                         const bool matched,
                         const int origin, const int PDG) {
  if (process == "primary" && matched && origin == 4 && PDG == 211) {
    return 1;
  }
  else if (process == "primary" && matched && origin == 4 && PDG == -13) {
    return 2;
  }
  else if (origin == 2) {
    return 3;
  }
  else if (process == "primaryBackground") {
    return 4;
  }
  else if (process.find("Inelastic") != std::string::npos) {
    return 5;
  }
  else if (process == "Decay") {
    return 6;
  }
  else {
    return 7;
  }
};


class truncatedMean_pos {
 private:
  double fLimit;
 public:
  truncatedMean_pos(double limit)
    : fLimit(limit) {}

  std::vector<double> operator()(
      std::vector<std::vector<double>> &vecs_dEdX) {
    //size_t size = 0;
    std::vector<double> trunc_mean;
    //std::vector<double> help_vec;
    double truncate_high = 1 - fLimit;
    int i_low = 0;
    int i_high = 0;

    //sort the dEdX vecotrs in matrix
    for(auto &&vec : vecs_dEdX){
       //size = vec.size();
       std::vector<double> help_vec;


       //check dEdX vector isn't empty!
       if(vec.empty()){
          trunc_mean.push_back(-9999.);
          continue;
       }

       else{
          std::vector<double> temp_vec;
          for (double d : vec) {
            if (d < 0.) {
              continue;
            }
            temp_vec.push_back(d);
          }
          for (double d : temp_vec) {
            if (d < 0.) {
              std::cout << d << std::endl;
            }
          }
          if (temp_vec.empty()) {
            trunc_mean.push_back(-9999.);
            continue;
          }
          //Sort Vector
          sort(temp_vec.begin(), temp_vec.end());
        
          //Discard upper and lower part of signal
          //rint rounds to integer
          i_low = rint ( temp_vec.size()*fLimit);
          i_high = rint( temp_vec.size()*truncate_high);
          
          
          //if (i_high >= temp_vec.size()) std::cout << "Warning: too high" << std::endl;
          for(int i = i_low; i </*=*/ i_high; i++){
            if (temp_vec[i] < 0) {
              std::cout << "added " << temp_vec[i] << " " << i << std::endl;
            }
            help_vec.push_back(temp_vec[i]);
          };

          //Mean of help vector

          trunc_mean.push_back(accumulate(help_vec.begin(), help_vec.end(), 0.0) / help_vec.size());
          if (trunc_mean.back() < 0 && trunc_mean.back() != -9999. && trunc_mean.back() != -999.) {
            std::cout << accumulate(help_vec.begin(), help_vec.end(), 0.0) << " " << help_vec.size() << std::endl;
            std::cout << temp_vec.size() << " " << i_low << " " << i_high << std::endl;
            for (size_t i = 0; i < help_vec.size(); ++i) {
              std::cout << "\t" << help_vec[i] << std::endl;
            }
          }
       }


    }
    return trunc_mean;
  }

};


class shower_dists {
 private:
  double fTrackScoreCut;
 public:
  shower_dists(double cut)
    : fTrackScoreCut(cut) {}
  std::vector<double> operator ()(const std::vector<double> &track_score,
                                  const std::vector<double> &shower_x,
                                  const std::vector<double> &shower_y,
                                  const std::vector<double> &shower_z,
                                  double & x, double & y, double & z) {
    std::vector<double> results;
    for(size_t i = 0; i < track_score.size(); ++i){
      if ((track_score[i] < fTrackScoreCut) &&
          (track_score[i] > 0.)) {
        double dist = sqrt(std::pow((shower_x[i] - x), 2) +
                           std::pow((shower_y[i] - y), 2) +
                           std::pow((shower_z[i] - z), 2));
        results.push_back(dist);
      }
      else {
        results.push_back(-999.);
      }
    }

    return results;
  }
};


class has_shower_dist_energy {
 private:
  double fTrackScoreCut;
 public:
  has_shower_dist_energy(double cut)
    : fTrackScoreCut(cut) {}
  bool operator()(const std::vector<double> &track_score,
                  const std::vector<double> &shower_x,
                  const std::vector<double> &shower_y,
                  const std::vector<double> &shower_z,
                  const std::vector<double> &energy,
                  double & x, double & y, double & z) {
    for(size_t i = 0; i < track_score.size(); ++i){
       double dist = sqrt(std::pow((shower_x[i] - x), 2) +
                          std::pow((shower_y[i] - y), 2) +
                          std::pow((shower_z[i] - z), 2));
       if ((track_score[i] < fTrackScoreCut) &&
           (track_score[i] > 0.) &&
           (dist > 5. && dist < 1000.) &&
           (energy[i] > 80. && energy[i] < 1000.)) {
         return true;
       }
    }

    return false;
  }
};

class isBeamType {
 private:
   bool fCheckCalo;
 public:
   isBeamType(bool check_calo)
     : fCheckCalo(check_calo) {}
   bool operator()(int reco_beam_type, std::vector<double> energies) {
     if (fCheckCalo) {
       return (energies.size() > 0 && reco_beam_type == 13);
     }

     return (reco_beam_type == 13);
   }
};

class endAPA3 {
 private:
  double fEndZCut;
 public:
  endAPA3(double cut)
    : fEndZCut(cut) {}
  
  bool operator()(double reco_beam_endZ) {
    return (reco_beam_endZ < fEndZCut);
  }
};

class secondary_noPion {
 private:
  double fTrackScoreCut;
  double fChi2Cut;
  double fdEdXLow, fdEdXMed, fdEdXHigh;
 public:
  secondary_noPion(double track_score_cut, double chi2_cut,
                   double dEdX_low, double dEdX_med, double dEdX_high)
    : fTrackScoreCut(track_score_cut),
      fChi2Cut(chi2_cut),
      fdEdXLow(dEdX_low),
      fdEdXMed(dEdX_med),
      fdEdXHigh(dEdX_high) {}
  bool operator()(
      const std::vector<double> &track_score, 
      const std::vector<int> &trackID,
      const std::vector<double> &dEdX,
      const std::vector<double> &chi2,
      const std::vector<int> &ndof) {
    for( size_t i = 0; i < track_score.size(); ++i ) {
      if ((trackID[i] != -999) && (track_score[i] > fTrackScoreCut)) {
        if (dEdX[i] < fdEdXMed/*2.8*/ && dEdX[i] > fdEdXLow/*0.5*/) {
          return false;
        }
        //else if (dEdX[i] > 2.8 && dEdX[i] < 3.4) {
        else if (dEdX[i] < fdEdXHigh) {
          if (ndof[i] > 0 && chi2[i]/ndof[i] > fChi2Cut) {
            return false;
          }
        }
      }
    }

    return true;
  }
};

auto data_beam_PID = [](const std::vector<int>& pidCandidates){
  auto pid_search = std::find(pidCandidates.begin(), pidCandidates.end(), 211);
  return (pid_search != pidCandidates.end());
};

class data_BI_quality {
 private:
   bool fDoTracks;
 public:
   data_BI_quality(bool do_tracks)
     : fDoTracks(do_tracks) {}
    bool operator()(int data_BI_nMomenta, int data_BI_nTracks) {

      if (fDoTracks) {
        return (data_BI_nMomenta == 1 && data_BI_nTracks == 1);
      }
      else {
        return (data_BI_nMomenta == 1);
      }
    }
};

class beam_cut_BI {
 private:
  double fXLow, fXHigh,
         fYLow, fYHigh,
         fZLow, fZHigh,
         fCosLow;
 public:
  beam_cut_BI(double x_low, double x_high,
              double y_low, double y_high,
              double z_low, double z_high,
              double cos_low)
    : fXLow(x_low), fXHigh(x_high),
      fYLow(y_low), fYHigh(y_high),
      fZLow(z_low), fZHigh(z_high),
      fCosLow(cos_low) {}
  bool operator()(double startX,
                  double startY,   double startZ,
                  double dirX,     double dirY,
                  double dirZ,     double BI_X,
                  double BI_Y,     double BI_dirX,
                  double BI_dirY,  double BI_dirZ) {

    double deltaX = startX - BI_X;
    double deltaY = startY - BI_Y;
    double cos = BI_dirX*dirX + BI_dirY*dirY +
                 BI_dirZ*dirZ;
    if( (deltaX < fXLow) || (deltaX > fXHigh) )
      return false;

    if ( (deltaY < fYLow) || (deltaY > fYHigh) )
      return false;

    if ( (startZ < fZLow) || (startZ > fZHigh) )
      return false;

    if (cos < fCosLow)
      return false;

    return true;
  }
};

class beam_cut_TPC {
 private:
  bool fDoAngle;
  //double fXCut, fYCut, fZCut, fCosLow;
  double fXYZCut, fCosCut;
  double fMeanX, fMeanY, fMeanZ;
  double fSigmaX, fSigmaY, fSigmaZ;
  double fMeanThetaX, fMeanThetaY, fMeanThetaZ;

 public:
  beam_cut_TPC(bool do_angle, double xyz_cut, double cos_cut,
               double mean_x, double mean_y, double mean_z,
               double sigma_x, double sigma_y, double sigma_z,
               double mean_theta_x, double mean_theta_y, double mean_theta_z)
   : fDoAngle(do_angle),
     fXYZCut(xyz_cut), fCosCut(cos_cut),
     fMeanX(mean_x), fMeanY(mean_y), fMeanZ(mean_z),
     fSigmaX(sigma_x), fSigmaY(sigma_y), fSigmaZ(sigma_z),
     fMeanThetaX(mean_theta_x), fMeanThetaY(mean_theta_y),
     fMeanThetaZ(mean_theta_z) {}

  bool operator()(double calo_beam_startX, double calo_beam_startY,
                  double calo_beam_startZ, double calo_beam_endX,
                  double calo_beam_endY, double calo_beam_endZ)   {

   double diffX = calo_beam_endX - calo_beam_startX;
   double diffY = calo_beam_endY - calo_beam_startY;
   double diffZ = calo_beam_endZ - calo_beam_startZ;
   double r = sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);

   double cosTrk_thetaX = diffX / r;
   double cosTrk_thetaY = diffY / r;
   double cosTrk_thetaZ = diffZ / r;

   double cosTheta = cos(fMeanThetaX*TMath::Pi()/180.)*cosTrk_thetaX +
                     cos(fMeanThetaY*TMath::Pi()/180.)*cosTrk_thetaY +
                     cos(fMeanThetaZ*TMath::Pi()/180.)*cosTrk_thetaZ;

   if (abs((calo_beam_startX - fMeanX)/fSigmaX) > fXYZCut)
      return false;

   if (abs((calo_beam_startY - fMeanY)/fSigmaY) > fXYZCut)
      return false;

   if (abs((calo_beam_startZ - fMeanZ)/fSigmaZ) > fXYZCut)
      return false;

   if (fDoAngle && (cosTheta < fCosCut))
      return false;

   return true;
   }
};
