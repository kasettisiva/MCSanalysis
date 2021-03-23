#include "AbsCexDriver.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

#include <iomanip>
#include <iostream> 

#include "ThinSliceDriverFactory.h"
DECLARE_THINSLICEDRIVER_FACTORY_NS(protoana::AbsCexDriver, protoana, AbsCexDriver)

protoana::AbsCexDriver::~AbsCexDriver() {}

protoana::AbsCexDriver::AbsCexDriver(
    const fhicl::ParameterSet & extra_options)
    : ThinSliceDriver(extra_options),
      fEnergyFix(extra_options.get<double>("EnergyFix")),
      fDoEnergyFix(extra_options.get<bool>("DoEnergyFix")),
      fPitch(extra_options.get<double>("WirePitch")),
      fZ0(extra_options.get<double>("Z0")),
      fEndZCut(extra_options.get<double>("EndZCut")),
      fSliceMethod(extra_options.get<std::string>("SliceMethod")) {
  fIn = new TFile("end_slices.root", "OPEN");
  fEndSlices = (TH2D*)fIn->Get("h2D")->Clone();
  std::cout << "Get end slices" << fEndSlices << " " <<
                fEndSlices->GetXaxis() << " " <<
                fEndSlices->GetXaxis()->GetNbins() << std::endl;
  for (int i = 1; i <= 735; ++i) {
    fMeans[i] = fEndSlices->ProjectionY("", i, i)->GetMean();
  }
  std::cout << "Get slices" << fEndSlices << " " <<
                fEndSlices->GetXaxis() << " " <<
                fEndSlices->GetXaxis()->GetNbins() << std::endl;
  if (fSliceMethod == "Traj") {
    fSliceCut = extra_options.get<int>("SliceCut");
  }
  else {
    fSliceCut = std::floor((fEndZCut - (fZ0 - fPitch/2.))/fPitch);
  }
}

void protoana::AbsCexDriver::BuildMCSamples(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::map<int, double> & nominal_fluxes,
    std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
    std::vector<double> & beam_energy_bins) {

  for (size_t i = 0; i < events.size(); ++i) {
    //std::cout << "Building from event " << i << std::endl;
    const ThinSliceEvent & event = events.at(i);

    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID();

    double true_beam_interactingEnergy = event.GetTrueInteractingEnergy();
    double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();
    double true_beam_endP = event.GetTrueEndP();
    double reco_beam_endZ = event.GetRecoEndZ();
    double true_beam_startP = event.GetTrueStartP();

    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();
    const std::vector<double> & true_beam_traj_Z
        = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE
        = event.GetTrueTrajKE();
    const std::vector<int> & true_beam_slices
        = event.GetTrueSlices();

    //double end_energy = (fExtraOptions.get<bool>("TrajSlices") ?
    //                     (*true_beam_traj_KE)[true_beam_traj_KE->size()-2] :
    //                     true_beam_interactingEnergy);
    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      //std::cout << "Chose alt " << fEndSlices << std::endl;
      //std::cout << "axis " << fEndSlices->GetXaxis() << std::endl;
      //std::cout << fEndSlices->GetXaxis()->GetNbins() << std::endl;
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z./*->*/back());
      //std::cout << "bin " << bin << std::endl;
      if (bin > 0)
        end_energy = fMeans.at(bin);
    }

    if (samples.find(sample_ID) == samples.end()) {
      std::cout << "Warning: skipping sample " << sample_ID << std::endl;
      continue;
    }

    //Build the true incident energy vector based on the slice cut
    std::vector<double> good_true_incEnergies;
    if (fSliceMethod == "Traj") {
      double next_slice_z = fExtraOptions.get<double>("TrajZStart");
      //slice_cut = fExtraOptions.get<int>("SliceCut");
      int next_slice_num = 0;
      //bool found = false;
      for (size_t j = 1; j < true_beam_traj_Z./*->*/size() - 1; ++j) {
        double z = (/***/true_beam_traj_Z)[j];
        double ke = (/***/true_beam_traj_KE)[j];

        if (z < fExtraOptions.get<double>("TrajZStart")) {
          //std::cout << "Skipping " << z << std::endl;
          continue;
        }

        if (z >= next_slice_z) {
          double temp_z = (/***/true_beam_traj_Z)[j-1];
          double temp_e = (/***/true_beam_traj_KE)[j-1];
          
          //if (/*end_energy > 800. && !found &&*/ j == true_beam_traj_Z->size() - 2) {
          //  //std::cout << "e, z i-1: " << temp_e << " " << temp_z <<
          //  //             " e, z i: " << (*true_beam_traj_KE)[j] << " " << (*true_beam_traj_Z)[j] <<
          //  //             " e, z i+1: " << (*true_beam_traj_KE)[j+1] << " " << (*true_beam_traj_Z)[j+1] << std::endl;
          //  //std::cout << "end e: " << end_energy << " slices: " << ((*true_beam_traj_Z)[j+1] - (*true_beam_traj_Z)[j])/pitch << std::endl;
          //  if (!(sample_ID == 4 || sample_ID == 5)) 
          //    h.Fill(end_energy, ((*true_beam_traj_Z)[j+1] - (*true_beam_traj_Z)[j])/pitch);
          //}
          //found = true;
          while (next_slice_z < z && next_slice_num < fSliceCut) {
            double sub_z = next_slice_z - temp_z;
            double delta_e = (/***/true_beam_traj_KE)[j-1] - ke;
            double delta_z = z - (/***/true_beam_traj_Z)[j-1]/* - z*/;
            temp_e -= (sub_z/delta_z)*delta_e;
          /*  if (end_energy > 800.) {
              std::cout << "delta_z: " << delta_z << std::endl;
              std::cout << (sub_z/delta_z)*delta_e << std::endl;
              std::cout << "\tNext: " << next_slice_z << " temp e: " << temp_e <<
                           " Z: " << z << " E: " << ke << " slice: " <<
                           next_slice_num << " " <<
                           (next_slice_z - fExtraOptions.get<double>("TrajZStart"))/next_slice_num <<std::endl;
              std::cout << "\t\tsub z: " << sub_z << " delta_z: " << delta_z << std::endl;
            }*/
            good_true_incEnergies.push_back(temp_e);
            temp_z = next_slice_z;
            next_slice_z += fPitch;
            ++next_slice_num;
          }
          //if (end_energy > 800.)
          //std::cout << "Next: " << next_slice_z << " Z: " << z << " slice: " << next_slice_num << std::endl;
        }
      }
    }
    else if (fSliceMethod == "Default") {
      for (size_t j = 0; j < true_beam_incidentEnergies./*->*/size(); ++j) {
        int slice = (/***/true_beam_slices)[j]; 
        if (slice > fSliceCut) continue;
        good_true_incEnergies.push_back((/***/true_beam_incidentEnergies)[j]);
      }
    }
    else if (fSliceMethod == "Alt") {
      //int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z->back());
      //std::cout << "alt 2 " << fEndSlices << std::endl;
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z.back());
      //std::cout << bin << std::endl;
      for (int i = 1; i <= bin; ++i) {
        good_true_incEnergies.push_back(fMeans.at(i));
      }
    }

    //Look for the coinciding energy bin
    int bin = -1;
    for (size_t j = 1; j < beam_energy_bins.size(); ++j) {
      if ((beam_energy_bins[j-1] <= 1.e3*true_beam_startP) &&
          (1.e3*true_beam_startP < beam_energy_bins[j])) {
        bin = j - 1;
        break;
      }
    }
    if (bin == -1) {
      std::string message = "Could not find beam energy bin for " +
                            std::to_string(true_beam_startP);
      throw std::runtime_error(message);
    }


    std::vector<ThinSliceSample> & samples_vec = samples.at(sample_ID)[bin];
    bool is_signal = signal_sample_checks.at(sample_ID);

    ThinSliceSample * this_sample = 0x0;
    if (!is_signal) {
      this_sample = &samples_vec.at(0);

      fluxes_by_sample[sample_ID][bin][0] += 1.;
    }
    else {
      //Iterate through the true bins and find the correct one
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {
          this_sample = &sample;

          fluxes_by_sample[sample_ID][bin][j] += 1.;
          found = true;
          break;
        }
      }
      if (!found) {
        //over/underflow here
        if (end_energy <= samples_vec[1].RangeLowEnd()) {
          this_sample = &samples_vec[0];

          fluxes_by_sample[sample_ID][bin][0] += 1.;
        }
        else if (end_energy >
                 samples_vec[samples_vec.size()-2].RangeHighEnd()) {
          this_sample = &samples_vec.back();
          fluxes_by_sample[sample_ID][bin].back() += 1.;
        }
        else {
          std::cout << "Warning: could not find true bin " <<
                       end_energy << std::endl;
        }
      }
    }

    int flux_type = this_sample->GetFluxType();
    nominal_fluxes[flux_type] += 1.;
    this_sample->AddFlux();
    double val[1] = {0};
    if (selection_ID == 4) {
      TH1D * selected_hist
          = (TH1D*)this_sample->GetSelectionHists().at(selection_ID);
      if (selected_hist->FindBin(reco_beam_endZ) == 0) {
        val[0] = selected_hist->GetBinCenter(1);
      }
      else if (selected_hist->FindBin(reco_beam_endZ) >
               selected_hist->GetNbinsX()) {
        val[0] = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
      }
      else {
        val[0] = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val[0] = .5;
    }
    else if (reco_beam_incidentEnergies./*->*/size()) {
      //this_sample->FillIncidentHist(*reco_beam_incidentEnergies);
      double energy[1] = {reco_beam_interactingEnergy};
      if (fExtraOptions.get<bool>("DoEnergyFix")) {
        for (size_t k = 1; k < reco_beam_incidentEnergies./*->*/size(); ++k) {
          double deltaE = ((/***/reco_beam_incidentEnergies)[k-1] -
                           (/***/reco_beam_incidentEnergies)[k]);
          if (deltaE > fExtraOptions.get<double>("EnergyFix")) {
            energy[0] += deltaE; 
          }
        }
      }

      TH1D * selected_hist
          = (TH1D*)this_sample->GetSelectionHists().at(selection_ID);
      if (selected_hist->FindBin(energy[0]) == 0) {
        val[0] = selected_hist->GetBinCenter(1);
      }
      else if (selected_hist->FindBin(energy[0]) >
               selected_hist->GetNbinsX()) {
        val[0] = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
      }
      else {
        val[0] = energy[0];
      }
    }
    else {
      TH1D * selected_hist
          = (TH1D*)this_sample->GetSelectionHists().at(selection_ID);
      val[0] = selected_hist->GetBinCenter(1);
    }
    this_sample->FillSelectionHist(selection_ID, val);

    //Fill the total incident hist with truth info
    this_sample->FillTrueIncidentHist(good_true_incEnergies);
    this_sample->AddIncidentEnergies(good_true_incEnergies);
    if (true_beam_incidentEnergies./*->*/size() > 0) {
    this_sample->AddESliceEnergies(
        {(/***/true_beam_incidentEnergies)[0],
         end_energy});
    }

  }

}

void protoana::AbsCexDriver::RefillMCSamples(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins,
    const std::map<int, std::vector<double>> & signal_pars,
    const std::map<int, double> & flux_pars,
    const std::map<std::string, ThinSliceSystematic> & syst_pars,
    bool fill_incident) {

  //Reset all samples
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        it->second[i][j].Reset();
      }
    }
  }

/*
  TH2D * h2D = 0x0;
  std::map<int, double> * means = 0x0;
  if (slice_method == "Alt") {
    tree->Draw("true_beam_traj_KE:true_beam_traj_Z>>h2D(734, -.49375, 351.63541, 200, 0, 2000)",
               "true_beam_PDG == 211 && new_interaction_topology != 4 && new_interaction_topology != 5 && true_beam_traj_Z > -.43975 && true_beam_traj_KE > 10.");
    h2D = (TH2D*)gDirectory->Get("h2D");
    means = new std::map<int, double>();
    for (int i = 1; i <= 735; ++i) {
      (*means)[i] = h2D->ProjectionY("", i, i)->GetMean();
    }
  }*/

  for (size_t i = 0; i < events.size(); ++i) {
    const ThinSliceEvent & event = events.at(i);
    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID();

    double true_beam_interactingEnergy = event.GetTrueInteractingEnergy();
    double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();
    double true_beam_endP = event.GetTrueEndP();
    double reco_beam_endZ = event.GetRecoEndZ();
    double true_beam_startP = event.GetTrueStartP();

    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();
    const std::vector<double> & true_beam_traj_Z
        = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE
        = event.GetTrueTrajKE();
    const std::vector<int> & true_beam_slices
        = event.GetTrueSlices();
    double beam_inst_P = event.GetBeamInstP();
    const std::vector<double> calibrated_dQdX
        = event.GetdQdXCalibrated();
    const std::vector<double> beam_EField
        = event.GetEField();
    const std::vector<double> track_pitch
        = event.GetTrackPitch();


    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z.back());
      if (bin > 0)
        end_energy = fMeans.at(bin);
    }

    if (samples.find(sample_ID) == samples.end()) {
      std::cout << "Warning: skipping sample " << sample_ID << std::endl;
      continue;
    }

    //Build the true incident energy vector based on the slice cut
    std::vector<double> good_true_incEnergies;
    if (fill_incident) {
    if (fSliceMethod == "Traj") {
      double next_slice_z = fExtraOptions.get<double>("TrajZStart");
      //slice_cut = fExtraOptions.get<int>("SliceCut");
      int next_slice_num = 0;
      for (size_t j = 1; j < true_beam_traj_Z.size() - 1; ++j) {
        double z = (true_beam_traj_Z)[j];
        double ke = (true_beam_traj_KE)[j];

        if (z < fExtraOptions.get<double>("TrajZStart")) {
          continue;
        }

        if (z >= next_slice_z) {
          double temp_z = (true_beam_traj_Z)[j-1];
          double temp_e = (true_beam_traj_KE)[j-1];
          
          while (next_slice_z < z && next_slice_num < fSliceCut) {
            double sub_z = next_slice_z - temp_z;
            double delta_e = (true_beam_traj_KE)[j-1] - ke;
            double delta_z = z - (true_beam_traj_Z)[j-1]/* - z*/;
            temp_e -= (sub_z/delta_z)*delta_e;
            good_true_incEnergies.push_back(temp_e);
            temp_z = next_slice_z;
            next_slice_z += fPitch;
            ++next_slice_num;
          }
        }
      }
    }
    else if (fSliceMethod == "Default") {
      for (size_t j = 0; j < true_beam_incidentEnergies.size(); ++j) {
        int slice = (true_beam_slices)[j]; 
        if (slice > fSliceCut) continue;
        good_true_incEnergies.push_back((true_beam_incidentEnergies)[j]);
      }
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z.back());
      for (int i = 1; i <= bin; ++i) {
        good_true_incEnergies.push_back(fMeans.at(i));
      }
    }
    }

    //Look for the coinciding energy bin
    int bin = -1;
    for (size_t j = 1; j < beam_energy_bins.size(); ++j) {
      if ((beam_energy_bins[j-1] <= 1.e3*true_beam_startP) &&
          (1.e3*true_beam_startP < beam_energy_bins[j])) {
        bin = j - 1;
        break;
      }
    }
    if (bin == -1) {
      std::string message = "Could not find beam energy bin for " +
                            std::to_string(true_beam_startP);
      throw std::runtime_error(message);
    }

    //Weight for the event
    //Possibly affected by signal, flux, and syst parameters
    double weight = 1.;

    std::vector<ThinSliceSample> & samples_vec = samples.at(sample_ID)[bin];
    bool is_signal = signal_sample_checks.at(sample_ID);

    ThinSliceSample * this_sample = 0x0;
    if (!is_signal) {
      this_sample = &samples_vec.at(0);
    }
    else {
      //Iterate through the true bins and find the correct one
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {
          this_sample = &sample;

          found = true;

          weight *= signal_pars.at(sample_ID)[j-1];

          break;
        }
      }
      if (!found) {
        //over/underflow here
        if (end_energy <= samples_vec[1].RangeLowEnd()) {
          this_sample = &samples_vec[0];

        }
        else if (end_energy >
                 samples_vec[samples_vec.size()-2].RangeHighEnd()) {
          this_sample = &samples_vec.back();
        }
        else {
          std::cout << "Warning: could not find true bin " <<
                       end_energy << std::endl;
        }
      }
    }

    int flux_type = this_sample->GetFluxType();
    if (flux_pars.find(flux_type) != flux_pars.end()) {
      weight *= flux_pars.at(flux_type);
    }

    double val[1] = {0};
    TH1D * selected_hist
        = (TH1D*)this_sample->GetSelectionHists().at(selection_ID);
    if (selection_ID == 4) {
      if (selected_hist->FindBin(reco_beam_endZ) == 0) {
        val[0] = selected_hist->GetBinCenter(1);
      }
      else if (selected_hist->FindBin(reco_beam_endZ) >
               selected_hist->GetNbinsX()) {
        val[0] = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
      }
      else {
        val[0] = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val[0] = .5;
    }
    else if (reco_beam_incidentEnergies.size()) {

      double energy[1] = {0.};
      if (syst_pars.find("dEdX_Cal") != syst_pars.end()) {
        //std::cout << "Doing dEdX_Cal" << std::endl;
        energy[0] = sqrt(beam_inst_P*beam_inst_P*1.e6 + 139.57*139.57) -
                        139.57;
        //limits?
        for (size_t k = 0; k < calibrated_dQdX.size()-1; ++k) {
          if ((calibrated_dQdX)[k] < 0.) continue;

          double dedx = (1./(syst_pars.at("dEdX_Cal").GetValue()/**fNominalCCal*/));
          dedx *= (calibrated_dQdX)[k];
          dedx *= (fBetaP / ( fRho * (beam_EField)[k] ) * fWion);
          dedx = exp(dedx);
          dedx -= fAlpha;
          dedx *= ((fRho*(beam_EField)[k])/fBetaP);

          if (dedx*(track_pitch)[k] > fEnergyFix)
            continue;
          energy[0] -= dedx*(track_pitch)[k];
          //std::cout << "Energy: " << energy[0] << " dedx: " << dedx <<
          //             std::endl;
        }
        //double int_energy = {reco_beam_interactingEnergy};
        //if (fDoEnergyFix) {
        //  for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
        //    double deltaE = ((reco_beam_incidentEnergies)[k-1] -
        //                     (reco_beam_incidentEnergies)[k]);
        //    if (deltaE > fEnergyFix) {
        //      int_energy += deltaE; 
        //    }
        //  }
        //}
        //std::cout << syst_pars.at("dEdX_Cal").GetValue() << " " << energy[0] << " " << int_energy << std::endl;
      }

      else {
        energy[0] = {reco_beam_interactingEnergy};
        if (fDoEnergyFix) {
          for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
            double deltaE = ((reco_beam_incidentEnergies)[k-1] -
                             (reco_beam_incidentEnergies)[k]);
            if (deltaE > fEnergyFix) {
              energy[0] += deltaE; 
            }
          }
        }
      }

      //TH1D * selected_hist
      //    = (TH1D*)this_sample->GetSelectionHists().at(selection_ID);
      if (selected_hist->FindBin(energy[0]) == 0) {
        val[0] = selected_hist->GetBinCenter(1);
      }
      else if (selected_hist->FindBin(energy[0]) >
               selected_hist->GetNbinsX()) {
        val[0] = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
      }
      else {
        val[0] = energy[0];
      }
    }
    else {
      //TH1D * selected_hist
      //    = (TH1D*)this_sample->GetSelectionHists().at(selection_ID);
      val[0] = selected_hist->GetBinCenter(1);
    }

    if (syst_pars.find("dEdX_Cal_Spline") != syst_pars.end()) {
      //std::cout << "Doing dEdX_Cal_spline" << std::endl;
      int bin = selected_hist->FindBin(val[0]);
      TSpline3 * spline
          = fFullSelectionSplines["dEdX_Cal_Spline"][selection_ID][bin-1];
      weight *= spline->Eval(syst_pars.at("dEdX_Cal_Spline").GetValue());
    }

    this_sample->FillSelectionHist(selection_ID, val, weight);

    //Fill the total incident hist with truth info
    if (fill_incident) {
      this_sample->FillTrueIncidentHist(good_true_incEnergies, weight);
      this_sample->AddIncidentEnergies(good_true_incEnergies, weight);
      if (true_beam_incidentEnergies.size() > 0) {
      this_sample->AddESliceEnergies(
          {(true_beam_incidentEnergies)[0],
           end_energy}, weight);
      }
    }

    this_sample->AddVariedFlux(weight);
  }

}

void protoana::AbsCexDriver::SetupSysts(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins,
    const std::map<std::string, ThinSliceSystematic> & pars,
    TFile & output_file) {

  if (pars.find("dEdX_Cal") != pars.end()) {
    //std::cout << "Found par dEdX_Cal" << std::endl;
    fhicl::ParameterSet cal_set
        = pars.at("dEdX_Cal").GetOption<fhicl::ParameterSet>("Cal_set");
    fBetaP = cal_set.get<double>("betap");
    fRho = cal_set.get<double>("Rho");
    fWion = cal_set.get<double>("Wion");
    fAlpha = cal_set.get<double>("alpha");

    std::vector<fhicl::ParameterSet> PlanePars
        = cal_set.get<std::vector<fhicl::ParameterSet>>("PlaneParameters");
    bool found_collection = false;
    for (auto & p : PlanePars) {
      if (p.get<int>("PlaneID") == 2) {
        fNominalCCal = p.get<double>("calib_factor");
        found_collection = true;
        break;
      }
    }
  
    if (!found_collection) {
      std::string message = "Could not find collection plane calibration factor";
      throw std::runtime_error(message);
    }
  }
  else if (pars.find("dEdX_Cal_Spline") != pars.end()) {
    fhicl::ParameterSet cal_set
        = pars.at("dEdX_Cal_Spline").GetOption<fhicl::ParameterSet>("Cal_set");
    fBetaP = cal_set.get<double>("betap");
    fRho = cal_set.get<double>("Rho");
    fWion = cal_set.get<double>("Wion");
    fAlpha = cal_set.get<double>("alpha");

    std::vector<fhicl::ParameterSet> PlanePars
        = cal_set.get<std::vector<fhicl::ParameterSet>>("PlaneParameters");
    bool found_collection = false;
    for (auto & p : PlanePars) {
      if (p.get<int>("PlaneID") == 2) {
        fNominalCCal = p.get<double>("calib_factor");
        found_collection = true;
        break;
      }
    }
  
    if (!found_collection) {
      std::string message = "Could not find collection plane calibration factor";
      throw std::runtime_error(message);
    }
    SystRoutine_dEdX_Cal(events, samples, pars, output_file);
  }
}

void protoana::AbsCexDriver::BuildSystSamples(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins
    /*,
    std::map<int, double> & nominal_fluxes,
    std::map<int, std::vector<double>> & fluxes_by_sample*/) {
  std::vector<fhicl::ParameterSet> routines
      = fExtraOptions.get<std::vector<fhicl::ParameterSet>>("Systematics");
  for (size_t i = 0; i < routines.size(); ++i) {
    fhicl::ParameterSet routine = routines[i];
    std::string name = routine.get<std::string>("Name");
    if (name == "G4RW") {
      SystRoutine_G4RW(tree, samples, signal_sample_checks, routine);
    }
    //if (name == "dEdX_Cal") {
    //  SystRoutine_dEdX_Cal(tree, samples, signal_sample_checks, routine,
    //                       beam_energy_bins);
    //}
    else {
      std::string message = "Could not find systematics routine " + name;
      throw std::runtime_error(message);
    }
  }
}

void protoana::AbsCexDriver::SystRoutine_G4RW(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    const fhicl::ParameterSet & routine) {
  return;
}

void protoana::AbsCexDriver::SystRoutine_dEdX_Cal(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<std::string, ThinSliceSystematic> & pars,
    TFile & output_file) {

  //Get the systematic variations to the calibration constant
  //then build systematic shift hists
  std::vector<double> C_cal_vars
      = pars.at("dEdX_Cal_Spline").GetOption<std::vector<double>>("C_cal_vars");

  //make sure the number of shifts are even
  //such that you have same number +/-
  if (C_cal_vars.size() % 2) {
    std::string message = "SystRoutine_dEdX_Cal Error: ";
    message += "odd number of variations to cal constant";
    throw std::runtime_error(message);
  }

  //Get the first sample and get the selection hists
  //also make full hist
  std::map<int, TH1D*> full_hists;

  ThinSliceSample & temp_sample = samples.begin()->second[0][0];
  const std::map<int, TH1*> & sel_hists = temp_sample.GetSelectionHists();
  for (auto it = sel_hists.begin(); it != sel_hists.end(); ++it) {
    std::string sel_hist_name = it->second->GetName();
    sel_hist_name += "Syst_dEdX_Cal_Spline";
    int shift_number = -2;

    fFullSelectionVars["dEdX_Cal_Spline"][it->first] = std::vector<TH1D*>();
    for (size_t k = 0; k < C_cal_vars.size(); ++k) {
      std::string shift_name = sel_hist_name;
      shift_name += std::to_string(shift_number);
      fFullSelectionVars["dEdX_Cal_Spline"][it->first].push_back((TH1D*)it->second->Clone(shift_name.c_str()));
      fFullSelectionVars["dEdX_Cal_Spline"][it->first].back()->Reset();

      ++shift_number;
      if (shift_number == 0)
        ++shift_number;
    }

    sel_hist_name += "_FullVar";
    full_hists[it->first] = (TH1D*)it->second->Clone(sel_hist_name.c_str());
    full_hists[it->first]->Reset();
  }

  for (size_t i = 0; i < events.size(); ++i) {
    const ThinSliceEvent & event = events.at(i);
    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID();

    double reco_beam_endZ = event.GetRecoEndZ();

    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    double beam_inst_P = event.GetBeamInstP();
    const std::vector<double> calibrated_dQdX
        = event.GetdQdXCalibrated();
    const std::vector<double> beam_EField
        = event.GetEField();
    const std::vector<double> track_pitch
        = event.GetTrackPitch();


    if (samples.find(sample_ID) == samples.end()) {
      std::cout << "Warning: skipping sample " << sample_ID << std::endl;
      continue;
    }

    std::vector<double> vals(C_cal_vars.size(), 0.);
    if (selection_ID == 4) {
      TH1D * selected_hist
          = fFullSelectionVars["dEdX_Cal_Spline"][selection_ID][0];
      if (selected_hist->FindBin(reco_beam_endZ) == 0) {
        for (double & v : vals) v = selected_hist->GetBinCenter(1);
      }
      else if (selected_hist->FindBin(reco_beam_endZ) >
               selected_hist->GetNbinsX()) {
        for (double & v : vals) v = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
      }
      else {
        for (double & v : vals) v = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      for (double & v : vals) v = .5;
    }
    else if (reco_beam_incidentEnergies.size()) {
      for (size_t j = 0; j < C_cal_vars.size(); ++j) {
        double energy = sqrt(beam_inst_P*beam_inst_P*1.e6 + 139.57*139.57) -
                        139.57;
        //limits?
        for (size_t k = 0; k < calibrated_dQdX.size()-1; ++k) {
          if ((calibrated_dQdX)[k] < 0.) continue;

          double dedx = (pars.at("dEdX_Cal_Spline").GetCentral()/C_cal_vars[j]);
          dedx *= (calibrated_dQdX)[k];
          dedx *= (fBetaP / ( fRho * (beam_EField)[k] ) * fWion);
          dedx = exp(dedx);
          dedx -= fAlpha;
          dedx *= ((fRho*(beam_EField)[k])/fBetaP);

          if (dedx*(track_pitch)[k] > fEnergyFix)
            continue;
          energy -= dedx*(track_pitch)[k];

          TH1D * selected_hist
              = fFullSelectionVars["dEdX_Cal_Spline"][selection_ID][0];
          if (selected_hist->FindBin(energy) == 0) {
            vals[j] = selected_hist->GetBinCenter(1);
          }
          else if (selected_hist->FindBin(energy) >
                   selected_hist->GetNbinsX()) {
            vals[j] = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
          }
          else {
            vals[j] = energy;
          }
        }
      }
    }
    else {
      TH1D * selected_hist
          = fFullSelectionVars["dEdX_Cal_Spline"][selection_ID][0];
      for (double & v : vals) v = selected_hist->GetBinCenter(1);
    }
    for (size_t j = 0; j < vals.size(); ++j) {
      fFullSelectionVars["dEdX_Cal_Spline"][selection_ID][j]->Fill(vals[j]);
    }
  }

  C_cal_vars.insert(C_cal_vars.begin() + C_cal_vars.size()/2,
                    pars.at("dEdX_Cal_Spline").GetCentral());

  TDirectory * dir = output_file.mkdir("dEdX_Cal_Spline_Syst");
  dir->cd();

  //Take the vars and make into ratios, then turn into splines
  for (auto it = fFullSelectionVars["dEdX_Cal_Spline"].begin(); 
       it != fFullSelectionVars["dEdX_Cal_Spline"].end(); ++it) {
    int selection_ID = it->first;

    //Build the full hist
    for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
      for (size_t i = 0; i < it2->second.size(); ++i) {
        for (size_t j = 0; j < it2->second[i].size(); ++j) {
          full_hists[selection_ID]->Add(it2->second[i][j].GetSelectionHist(selection_ID));
        }
      }
    }

    for (size_t i = 0; i < it->second.size(); ++i) {
      it->second[i]->Write();
      it->second[i]->Divide(full_hists[selection_ID]);
    }
    
    fFullSelectionSplines["dEdX_Cal_Spline"][selection_ID] = std::vector<TSpline3*>();
    for (int i = 1; i <= full_hists[selection_ID]->GetNbinsX(); ++i) {
      std::vector<double> vals;
      for (size_t j = 0; j < it->second.size(); ++j) {
        vals.push_back(it->second[j]->GetBinContent(i));
      }
      vals.insert(vals.begin() + vals.size()/2, 1.);

      std::string spline_name = full_hists[selection_ID]->GetName();
      spline_name += "_Spline" + std::to_string(i);

      fFullSelectionSplines["dEdX_Cal_Spline"][selection_ID].push_back(
        new TSpline3(spline_name.c_str(), &C_cal_vars[0], &vals[0], vals.size()));
      TCanvas c(spline_name.c_str(), "");
      fFullSelectionSplines["dEdX_Cal_Spline"][selection_ID].back()->Draw();
      c.Write();
    }
  }
}

void protoana::AbsCexDriver::BuildDataHists(
    TTree * tree, ThinSliceDataSet & data_set, double & flux) {
  int selection_ID; 
  double reco_beam_interactingEnergy, reco_beam_endZ;
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                        &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                        &reco_beam_incidentEnergies);

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  flux = tree->GetEntries();

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
          double energy = reco_beam_interactingEnergy;
          if (fExtraOptions.get<bool>("DoEnergyFix")) {
            for (size_t k = 1; k < reco_beam_incidentEnergies->size(); ++k) {
              double deltaE = ((*reco_beam_incidentEnergies)[k-1] -
                               (*reco_beam_incidentEnergies)[k]);
              if (deltaE > fExtraOptions.get<double>("EnergyFix")) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val);
  }

}

void protoana::AbsCexDriver::BuildFakeData(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales) {
  std::string routine = fExtraOptions.get<std::string>("FakeDataRoutine");
  if (routine == "SampleScales") {
    FakeDataSampleScales(tree, samples, signal_sample_checks, data_set, flux,
                         sample_scales);
  }
  else if (routine == "G4RW") {
    FakeDataG4RW(tree, samples, signal_sample_checks, data_set, flux,
                 sample_scales);
  }
  else if (routine == "BinnedScales") {
    FakeDataBinnedScales(tree, samples, signal_sample_checks, data_set, flux,
                         sample_scales);
  }
}

void protoana::AbsCexDriver::FakeDataSampleScales(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales) {


  //Build the map for fake data scales
  std::vector<std::pair<int, double>> temp_vec
      = fExtraOptions.get<std::vector<std::pair<int, double>>>(
          "FakeDataScales");
  std::map<int, double> fake_data_scales(temp_vec.begin(), temp_vec.end()); 

  int sample_ID, selection_ID; 
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  double true_beam_endP;
  std::vector<double> * reco_beam_incidentEnergies = 0x0,
                      * true_beam_traj_Z = 0x0;
  double reco_beam_endZ;
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  tree->SetBranchAddress("new_interaction_topology", &sample_ID);
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("true_beam_interactingEnergy",
                         &true_beam_interactingEnergy);
  tree->SetBranchAddress("true_beam_endP", &true_beam_endP);
  tree->SetBranchAddress("true_beam_traj_Z", &true_beam_traj_Z);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                         &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                         &reco_beam_incidentEnergies);

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  double new_flux = 0.;
  //flux = 0.;
  flux = tree->GetEntries(); 

  tree->Draw("true_beam_traj_KE:true_beam_traj_Z>>h2D(734, -.49375, 351.63541, 200, 0, 2000)", "true_beam_PDG == 211 && new_interaction_topology != 4 && new_interaction_topology != 5 && true_beam_traj_Z > -.43975 && true_beam_traj_KE > 10.");
  TH2D * h2D = (TH2D*)gDirectory->Get("h2D");
  std::map<int, double> means;
  for (int i = 1; i <= 735; ++i) {
    means[i] = h2D->ProjectionY("", i, i)->GetMean();
  }

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    //double end_energy = (fExtraOptions.get<bool>("TrajSlices") ?
    //                     sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57 :
    //                     true_beam_interactingEnergy);
    std::string slice_method = fExtraOptions.get<std::string>("SliceMethod");
    double end_energy = true_beam_interactingEnergy;
    if (slice_method == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (slice_method == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (slice_method == "Alt") {
      int bin = h2D->GetXaxis()->FindBin(true_beam_traj_Z->back());
      if (bin > 0)
        end_energy = means.at(bin);
    }

    //Add under/overflow check here
    if (samples.find(sample_ID) == samples.end())
      continue;

    double scale = (fake_data_scales.find(sample_ID) != fake_data_scales.end() ?
                    fake_data_scales.at(sample_ID) : 1.);

    //If it's signal
    //Determine if the over/underflow bin 
    bool is_signal = signal_sample_checks.at(sample_ID);
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][0];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {     
          found = true;
          sample_scales[sample_ID][j] += scale;
          nominal_samples[sample_ID][j] += 1.;
          break;
        }
      }
      if (!found) {//If in the under/overflow, just set to 1. 
        scale = 1.;
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          sample_scales[sample_ID][0] += scale;
          nominal_samples[sample_ID][0] += 1.;
        }
        else {
          sample_scales[sample_ID].back() += scale;
          nominal_samples[sample_ID].back() += 1.;
        }
      }
    }
    else {
      sample_scales[sample_ID][0] += scale;
      nominal_samples[sample_ID][0] += 1.;
    }

    new_flux += scale; //1 or scaled

    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
          double energy = reco_beam_interactingEnergy;
          if (fExtraOptions.get<bool>("DoEnergyFix")) {
            for (size_t k = 1; k < reco_beam_incidentEnergies->size(); ++k) {
              double deltaE = ((*reco_beam_incidentEnergies)[k-1] -
                               (*reco_beam_incidentEnergies)[k]);
              if (deltaE > fExtraOptions.get<double>("EnergyFix")) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val, scale);
  }

  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (it->second[i] > 0.) {
        it->second[i] /= nominal_samples[it->first][i];
      }
      else {
        it->second[i] = 1.;
      }
      it->second[i] *= (flux/new_flux);
    }
  }

  incident_hist.Scale(flux/new_flux);
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Scale(flux/new_flux);
  }
}

void protoana::AbsCexDriver::FakeDataBinnedScales(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales) {

  //Build the map for fake data scales
  std::vector<std::pair<int, std::vector<double>>> temp_vec
      = fExtraOptions.get<std::vector<std::pair<int, std::vector<double>>>>(
          "FakeDataBinnedScales");
  std::map<int, std::vector<double>> fake_data_scales(temp_vec.begin(), temp_vec.end()); 

  int sample_ID, selection_ID; 
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  double true_beam_endP;
  std::vector<double> * reco_beam_incidentEnergies = 0x0,
                      * true_beam_traj_Z = 0x0;
  double reco_beam_endZ;
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  tree->SetBranchAddress("new_interaction_topology", &sample_ID);
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("true_beam_interactingEnergy",
                         &true_beam_interactingEnergy);
  tree->SetBranchAddress("true_beam_endP", &true_beam_endP);
  tree->SetBranchAddress("true_beam_traj_Z", &true_beam_traj_Z);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                         &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                         &reco_beam_incidentEnergies);

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();


  //Add check here 
  tree->Draw("true_beam_traj_KE:true_beam_traj_Z>>h2D(734, -.49375, 351.63541, 200, 0, 2000)", "true_beam_PDG == 211 && new_interaction_topology != 4 && new_interaction_topology != 5 && true_beam_traj_Z > -.43975 && true_beam_traj_KE > 10.");
  TH2D * h2D = (TH2D*)gDirectory->Get("h2D");
  std::map<int, double> means;
  for (int i = 1; i <= 735; ++i) {
    means[i] = h2D->ProjectionY("", i, i)->GetMean();
  }

  double new_flux = 0.;
  //flux = 0.;
  flux = tree->GetEntries(); 

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    std::string slice_method = fExtraOptions.get<std::string>("SliceMethod");
    double end_energy = true_beam_interactingEnergy;
    if (slice_method == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (slice_method == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (slice_method == "Alt") {
      int bin = h2D->GetXaxis()->FindBin(true_beam_traj_Z->back());
      if (bin > 0) {
        end_energy = means.at(bin);
      }
    }

    //Add under/overflow check here
    if (samples.find(sample_ID) == samples.end())
      continue;
    double scale = 1.;

    bool is_scaled = (fake_data_scales.find(sample_ID) !=
                      fake_data_scales.end());

    //If it's signal
    //Determine if the over/underflow bin 
    bool is_signal = signal_sample_checks.at(sample_ID);
    if (is_signal && is_scaled) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][0];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {     
          found = true;
          scale = fake_data_scales.at(sample_ID)[j];
          sample_scales[sample_ID][j] += scale;
          nominal_samples[sample_ID][j] += 1.;
          break;
        }
      }
      if (!found) {
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          scale = fake_data_scales.at(sample_ID)[0];
          sample_scales[sample_ID][0] += scale;
          nominal_samples[sample_ID][0] += 1.;
        }
        else {
          scale = fake_data_scales.at(sample_ID).back();
          sample_scales[sample_ID].back() += scale;
          nominal_samples[sample_ID].back() += 1.;
        }
      }
    }

    if (!is_signal && is_scaled) {
      scale = fake_data_scales.at(sample_ID)[0];
      sample_scales[sample_ID][0] += scale;
      nominal_samples[sample_ID][0] += 1.;
    }


    new_flux += scale; //1 or scaled

/*
    if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j], scale);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4) {
          double energy = reco_beam_interactingEnergy;
          if (fExtraOptions.get<bool>("DoEnergyFix")) {
            for (size_t k = 1; k < reco_beam_incidentEnergies->size(); ++k) {
              double deltaE = ((*reco_beam_incidentEnergies)[k-1] -
                               (*reco_beam_incidentEnergies)[k]);
              if (deltaE > fExtraOptions.get<double>("EnergyFix")) {
                energy += deltaE; 
              }
            }
          }
          selected_hists[selection_ID]->Fill(energy, scale);
        }
        else {
          selected_hists[selection_ID]->Fill(reco_beam_endZ, scale);
        }
      }
    }
    */
    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
          double energy = reco_beam_interactingEnergy;
          if (fExtraOptions.get<bool>("DoEnergyFix")) {
            for (size_t k = 1; k < reco_beam_incidentEnergies->size(); ++k) {
              double deltaE = ((*reco_beam_incidentEnergies)[k-1] -
                               (*reco_beam_incidentEnergies)[k]);
              if (deltaE > fExtraOptions.get<double>("EnergyFix")) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val, scale);
  }

  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (it->second[i] > 0.) {
        it->second[i] /= nominal_samples[it->first][i];
      }
      else {
        it->second[i] = 1.;
      }
      it->second[i] *= (flux/new_flux);
    }
  }

  std::cout << "Flux: " << flux << " new_flux: " << new_flux << std::endl;

  incident_hist.Scale(flux/new_flux);
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Scale(flux/new_flux);
  }
}

void protoana::AbsCexDriver::FakeDataG4RW(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales) {

  //Build the map for fake data scales
  fhicl::ParameterSet g4rw_options 
      = fExtraOptions.get<fhicl::ParameterSet>("FakeDataG4RW");

  int sample_ID, selection_ID; 
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  double true_beam_endP;
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  double reco_beam_endZ;
  std::vector<double> * g4rw_weight = 0x0;
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  tree->SetBranchAddress("new_interaction_topology", &sample_ID);
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("true_beam_interactingEnergy",
                         &true_beam_interactingEnergy);
  tree->SetBranchAddress("true_beam_endP", &true_beam_endP);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                         &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                         &reco_beam_incidentEnergies);
  std::vector<double> * true_beam_traj_Z = 0x0;
  tree->SetBranchAddress("true_beam_traj_Z", &true_beam_traj_Z);

  if (g4rw_options.get<int>("Shift") > 0) {
    if (g4rw_options.get<bool>("Full")) {
      std::cout << "Full weight" << std::endl;
      tree->SetBranchAddress("g4rw_full_primary_plus_sigma_weight", &g4rw_weight);
    }
    else {
      tree->SetBranchAddress("g4rw_alt_primary_plus_sigma_weight", &g4rw_weight);
    }
  }
  else {
    if (g4rw_options.get<bool>("Full")) {
      std::cout << "Full weight" << std::endl;
      tree->SetBranchAddress("g4rw_full_primary_minus_sigma_weight", &g4rw_weight);
    }
    else {
      tree->SetBranchAddress("g4rw_alt_primary_minus_sigma_weight", &g4rw_weight);
    }
  }

  size_t g4rw_pos = g4rw_options.get<size_t>("Position");

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  double new_flux = 0.;
  flux = tree->GetEntries();
  

  tree->Draw("true_beam_traj_KE:true_beam_traj_Z>>h2D(734, -.49375, 351.63541, 200, 0, 2000)", "true_beam_PDG == 211 && new_interaction_topology != 4 && new_interaction_topology != 5 && true_beam_traj_Z > -.43975 && true_beam_traj_KE > 10.");
  TH2D * h2D = (TH2D*)gDirectory->Get("h2D");
  std::map<int, double> means;
  for (int i = 1; i <= 735; ++i) {
    means[i] = h2D->ProjectionY("", i, i)->GetMean();
  }

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    if (samples.find(sample_ID) == samples.end())
      continue;

    std::string slice_method = fExtraOptions.get<std::string>("SliceMethod");
    double end_energy = true_beam_interactingEnergy;
    if (slice_method == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (slice_method == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (slice_method == "Alt") {
      int bin = h2D->GetXaxis()->FindBin(true_beam_traj_Z->back());
      if (bin > 0) {
        end_energy = means.at(bin);
      }
    }

    double scale = 1.;
    if (g4rw_weight->size() > 0) {
      scale = g4rw_weight->at(g4rw_pos);
    }

    bool is_signal = signal_sample_checks.at(sample_ID);
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][0];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {     
          found = true;
          sample_scales[sample_ID][j] += scale;
          nominal_samples[sample_ID][j] += 1.;
          break;
        }
      }
      if (!found) {
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          sample_scales[sample_ID][0] += scale;
          nominal_samples[sample_ID][0] += 1.;
        }
        else {
          sample_scales[sample_ID].back() += scale;
          nominal_samples[sample_ID].back() += 1.;
        }
      }
    }
    else {
      sample_scales[sample_ID][0] += scale;
      nominal_samples[sample_ID][0] += 1.;
    }

    new_flux += scale; //1 or scaled
    /*
    if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j], scale);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4) {
          double energy = reco_beam_interactingEnergy;
          if (fExtraOptions.get<bool>("DoEnergyFix")) {
            for (size_t k = 1; k < reco_beam_incidentEnergies->size(); ++k) {
              double deltaE = ((*reco_beam_incidentEnergies)[k-1] -
                               (*reco_beam_incidentEnergies)[k]);
              if (deltaE > fExtraOptions.get<double>("EnergyFix")) {
                energy += deltaE; 
              }
            }
          }
          //selected_hists[selection_ID]->Fill(reco_beam_interactingEnergy, scale);
          selected_hists[selection_ID]->Fill(energy, scale);
        }
        else {
          selected_hists[selection_ID]->Fill(reco_beam_endZ, scale);
        }
      }
    }
    */
    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
          double energy = reco_beam_interactingEnergy;
          if (fExtraOptions.get<bool>("DoEnergyFix")) {
            for (size_t k = 1; k < reco_beam_incidentEnergies->size(); ++k) {
              double deltaE = ((*reco_beam_incidentEnergies)[k-1] -
                               (*reco_beam_incidentEnergies)[k]);
              if (deltaE > fExtraOptions.get<double>("EnergyFix")) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val, scale);
  }

  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (it->second[i] > 0.) {
        it->second[i] /= nominal_samples[it->first][i];
      }
      else {
        it->second[i] = 1.;
      }
      it->second[i] *= (flux/new_flux);
    }
  }

  incident_hist.Scale(flux/new_flux);
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Scale(flux/new_flux);
  }

  //std::cout << "Fluxes: " << flux << " " << new_flux << std::endl;
}

std::pair<double, size_t> protoana::AbsCexDriver::CalculateChi2(
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    ThinSliceDataSet & data_set) {

  double chi2 = 0.;
  size_t nPoints = 0;

  std::map<int, TH1 *> & selected_data_hists = data_set.GetSelectionHists();
  double total_data = 0., total_mc = 0.;
  for (auto it = selected_data_hists.begin();
       it != selected_data_hists.end(); ++it) {
    TH1D * data_hist = (TH1D*)it->second;
    int selection_ID = it->first;
    if (data_hist->GetBinContent(0) > 0.) {
      std::cout << "Warning: underflow bin of " << selection_ID <<
                   " has " << data_hist->GetBinContent(0) << " events" <<
                   std::endl;
    }
    else if (data_hist->GetBinContent(data_hist->GetNbinsX()+1) > 0.) {
      std::cout << "Warning: overflow bin of " << selection_ID <<
                   " has " << data_hist->GetBinContent(data_hist->GetNbinsX()+1) <<
                   " events" <<
                   std::endl;
    }
    for (int i = 1; i <= data_hist->GetNbinsX(); ++i) {
      double data_val = data_hist->GetBinContent(i);
      //double data_err = data_hist->GetBinError(i);

      double mc_val = 0.;
      //Go through all the samples and get the values from mc
      for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
        std::vector<std::vector<ThinSliceSample>> & samples_vec_2D = it2->second;
        for (size_t j = 0; j < samples_vec_2D.size(); ++j) {
          std::vector<ThinSliceSample> & samples_vec = samples_vec_2D[j];
          for (size_t k = 0; k < samples_vec.size(); ++k) {
            ThinSliceSample & sample = samples_vec[k];
            mc_val += sample.GetSelectionHist(selection_ID)->GetBinContent(i);

            TH1D * mc_hist = (TH1D*)sample.GetSelectionHist(selection_ID);
            if (mc_hist->GetBinContent(0) > 0.) {
              std::cout << "Warning: underflow bin of " << selection_ID <<
                           " has " << mc_hist->GetBinContent(0) << " events" <<
                           std::endl;
            }
            else if (mc_hist->GetBinContent(mc_hist->GetNbinsX()+1) > 0.) {
              std::cout << "Warning: overflow bin of " << selection_ID <<
                           " has " << mc_hist->GetBinContent(mc_hist->GetNbinsX()+1) <<
                           " events" <<
                           std::endl;
            }
          }
        }
      }
      //chi2 += (std::pow((data_val - mc_val), 2) / std::pow(data_err, 2));
      chi2 += 2*data_val*std::log(data_val/mc_val);
      //std::cout << "data val: " << data_val << " mc val: " << mc_val << " current chi2: " << chi2 << std::endl;
      ++nPoints;
      total_mc += mc_val;
      total_data += data_val;
    }
  }
  /*
  std::cout << "totals (data, mc): " << total_data << " " << total_mc << std::endl;
  if (total_data != total_mc) {
    std::cout << "Warning: data != mc " << total_data << " " << total_mc <<
                 std::endl;
  }*/

  return {chi2, nPoints};
}

void protoana::AbsCexDriver::CompareSelections(
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    ThinSliceDataSet & data_set, TFile & output_file,
    std::vector<std::pair<int, int>> plot_style,
    bool plot_rebinned,
    bool post_fit, int nPars) {

  output_file.cd();
  std::map<int, TH1*> data_hists
      = (plot_rebinned ?
         data_set.GetRebinnedSelectionHists() :
         data_set.GetSelectionHists());
  for (auto it = data_hists.begin(); it != data_hists.end(); ++it) {
    int selection_ID = it->first;
    TH1D * data_hist = (TH1D*)it->second;
    data_hist->SetLineColor(kBlack);
    data_hist->SetMarkerColor(kBlack);
    data_hist->SetMarkerStyle(20);

    std::string canvas_name = "c";
    canvas_name += (post_fit ? "PostFit" : "Nominal") +
                   data_set.GetSelectionName(selection_ID);
    TCanvas cSelection(canvas_name.c_str(), "");
    cSelection.SetTicks();

    std::map<int, std::vector<TH1D *>> temp_hists;
    for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
      temp_hists[it2->first] = std::vector<TH1D *>();
      std::vector<ThinSliceSample> & vec = it2->second[0];
      for (size_t i = 0; i < vec.size(); ++i) {
        temp_hists[it2->first].push_back((TH1D*)(
            plot_rebinned ?
            vec[i].GetRebinnedSelectionHist(selection_ID) :
            vec[i].GetSelectionHist(selection_ID))->Clone());
      }
      for (size_t i = 1; i < it2->second.size(); ++i) {
        for (size_t j = 0; j < it2->second[i].size(); ++j) {
          temp_hists[it2->first][j]->Add((TH1D*)(
              plot_rebinned ?
              it2->second[i][j].GetRebinnedSelectionHist(selection_ID) :
              it2->second[i][j].GetSelectionHist(selection_ID)));
          }
      }
    }
    //Build the mc stack
    std::string stack_name = (post_fit ? "PostFit" : "Nominal") +
                             data_set.GetSelectionName(selection_ID) +
                             "Stack";
    THStack mc_stack(stack_name.c_str(), "");
    TLegend leg;
    std::vector<TH1D*> temp_vec;
    size_t iColor = 0;
    //need to add second loop with temp hists
    for (auto it2 = temp_hists.begin(); it2 != temp_hists.end(); ++it2) {
      for (size_t i = 0; i < it2->second.size(); ++i) {
        TH1D * sel_hist = it2->second.at(i);
        std::pair<int, int> color_fill = GetColorAndStyle(iColor, plot_style);
        sel_hist->SetFillColor(color_fill.first);
        sel_hist->SetFillStyle(color_fill.second);
        sel_hist->SetLineColor(kBlack);
        mc_stack.Add(sel_hist);
        temp_vec.push_back(sel_hist);
        ++iColor;
      }
    }
    mc_stack.Write();

    for (auto it = temp_vec.rbegin(); it != temp_vec.rend(); ++it) {
      leg.AddEntry(*it);
    }
    leg.AddEntry(data_hist, "Data");
    
    std::pair<double, size_t> chi2 = CalculateChi2(samples, data_set);
    std::string chi2_str = "#chi^{2}/ndof = " +
                           std::to_string(chi2.first) + "/" +
                           std::to_string(chi2.second - nPars);
    leg.AddEntry((TObject*)0x0, chi2_str.c_str(), "");

    mc_stack.Draw("hist");
    double max_mc = mc_stack.GetHistogram()->GetMaximum();
    int max_data_bin = data_hist->GetMaximumBin();
    double max_data = data_hist->GetBinContent(max_data_bin) +
                      data_hist->GetBinError(max_data_bin);
    std::string title = data_set.GetSelectionName(selection_ID);
    title += ";";
    title += data_hist->GetXaxis()->GetTitle();
    //if (selection_ID == 4) {
    //  title += ";Reconstructed End Z (cm)";
    //}
    //else if (selection_ID < 4) {
    //  title += ";Reconstructed KE (MeV)";
    //}
    mc_stack.SetTitle(title.c_str());

    mc_stack.GetHistogram()->SetTitleSize(.04, "X");
    mc_stack.SetMaximum((max_data > max_mc ? max_data : max_mc));
    mc_stack.Draw("hist");
    data_hist->Draw("e1 same");
    leg.Draw("same");

    cSelection.Write();

    //Get the full incident hist from stack
    TList * l = (TList*)mc_stack.GetHists();
    TH1D * hMC = (TH1D*)l->At(0)->Clone();
    for (int i = 1; i < l->GetSize(); ++i) {
      hMC->Add((TH1D*)l->At(i));
    }

    std::string ratio_name = data_set.GetSelectionName(selection_ID) + "Ratio" +
                             (post_fit ? "PostFit" : "Nominal");
    TH1D * hRatio
        = (TH1D*)data_hist->Clone(ratio_name.c_str());
    hRatio->Divide(hMC);
    hRatio->Write(); 

    canvas_name += "Ratio";
    TCanvas cRatio(canvas_name.c_str(), "");
    cRatio.SetTicks();
    TPad p1((canvas_name + "pad1").c_str(), "", 0, 0.3, 1., 1.);
    p1.SetBottomMargin(0.1);
    p1.Draw();
    p1.cd();
    mc_stack.Draw("hist");
    mc_stack.GetHistogram()->SetTitle("Abs;;");
    for (int i = 1; i < mc_stack.GetHistogram()->GetNbinsX(); ++i) {
      hRatio->GetXaxis()->SetBinLabel(i, mc_stack.GetHistogram()->GetXaxis()->GetBinLabel(i));
      mc_stack.GetHistogram()->GetXaxis()->SetBinLabel(i, "");
    }
    mc_stack.Draw("hist");
    data_hist->Draw("e1 same");

    cRatio.cd();
    TPad p2((canvas_name + "pad2").c_str(), "", 0, 0, 1., 0.3);
    p2.SetTopMargin(0.1);
    p2.SetBottomMargin(.2);
    p2.Draw();
    p2.cd();
    hRatio->Sumw2();
    std::string r_title = "";/*(selection_ID != 4 ?
                           ";Reconstructed KE (MeV)" :
                           ";Reconstructed End Z (cm)");*/
    //if (selection_ID == 4) {
    //  r_title += ";Reconstructed End Z (cm)";
    //}
    //else if (selection_ID < 4) {
    //  r_title += ";Reconstructed KE (MeV)";
    //}
    //r_title += ";Data/MC";
    //hRatio->SetTitle(r_title.c_str());
    hRatio->GetYaxis()->SetTitle("Data/MC");
    hRatio->SetTitleSize(.04, "XY");
    hRatio->Draw("ep");

    cRatio.Write();

  }
}

void protoana::AbsCexDriver::GetCurrentHists(
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    ThinSliceDataSet & data_set,
    std::map<int, std::vector<TH1*>> & throw_hists,
    bool plot_rebinned) {
  
  //Use data set as a template
  std::map<int, TH1*> data_hists
      = (plot_rebinned ?
         data_set.GetRebinnedSelectionHists() :
         data_set.GetSelectionHists());

  for (auto it = data_hists.begin(); it != data_hists.end(); ++it) {
    std::vector<double> bins;
    std::string name = data_set.GetSelectionName(it->first) + "Throw" +
                       std::to_string(throw_hists[it->first].size());

    TH1D * temp_hist = (TH1D*)it->second->Clone(name.c_str());
    temp_hist->Reset();
    throw_hists[it->first].push_back(temp_hist);
  }

  //Iterate over samples
  for (auto it = throw_hists.begin(); it != throw_hists.end(); ++it) {
    int selection_ID = it->first;
    for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
      std::vector<std::vector<ThinSliceSample>> & samples_vec_2D = it2->second;
      for (size_t j = 0; j < samples_vec_2D.size(); ++j) {
        std::vector<ThinSliceSample> & samples_vec = samples_vec_2D[j];
        for (size_t k = 0; k < samples_vec.size(); ++k) {
          ThinSliceSample & sample = samples_vec[k];
          for (int i = 1; i <= it->second.back()->GetNbinsX(); ++i) {
            it->second.back()->AddBinContent(i,
                sample.GetSelectionHist(selection_ID)->GetBinContent(i));
          }
        }
      }
    }
  }
}

void protoana::AbsCexDriver::GetCurrentTruthHists(
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    std::map<int, std::vector<TH1*>> & throw_hists,
    std::map<int, std::vector<TH1*>> & throw_inc_hists,
    std::map<int, std::vector<TH1*>> & throw_xsec_hists,
    const std::vector<int> & incident_samples,
    const std::map<int, std::vector<double>> & signal_bins) {
  //Loop over the samples
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    //Get the number of bins from the first entry of the beam energy bins
    std::vector<std::vector<ThinSliceSample>> & samples_vec_2D = it->second;
    size_t nBins = samples_vec_2D[0].size();
    std::string name = it->second[0][0].GetName() + "Throw" +
                       std::to_string(throw_hists[it->first].size());
    TH1D * temp_hist = new TH1D(name.c_str(), "", nBins, 0, nBins); 
    for (size_t i = 0; i < samples_vec_2D.size(); ++i) {
      std::vector<ThinSliceSample> & samples_vec = samples_vec_2D[i];
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        temp_hist->AddBinContent(j+1, samples_vec[j].GetVariedFlux());
      }
    }
    throw_hists[it->first].push_back(temp_hist);
  }

  for (auto it = throw_inc_hists.begin(); it != throw_inc_hists.end(); ++it) {
    int s = it->first;
    auto & samples_vec_2D = samples[s];
    const std::vector<double> & bins = signal_bins.at(s);
    std::string name = samples_vec_2D[0][0].GetName();
    name += "IncidentThrow" +
             std::to_string(throw_inc_hists[it->first].size());
    TH1D * temp_inc_hist = new TH1D(name.c_str(), "", bins.size() - 1, &bins[0]); 
    

    name = samples_vec_2D[0][0].GetName();
    name += "XSecThrow" +
             std::to_string(throw_inc_hists[it->first].size());
    TH1D * temp_xsec_hist = new TH1D(name.c_str(), "", bins.size() - 1,
                                     &bins[0]);
    for (auto i_s : incident_samples) {
      auto & incident_vec_2D = samples[i_s];
      for (size_t i = 0; i < incident_vec_2D.size(); ++i) {
        for (size_t j = 0; j < incident_vec_2D[i].size(); ++j) {
          if (fExtraOptions.get<std::string>("SliceMethod") == "E") {
            incident_vec_2D[i][j].FillESliceHist(*temp_inc_hist);
          }
          else {
            incident_vec_2D[i][j].FillHistFromIncidentEnergies(*temp_inc_hist);
          }
        }
      }
    }
    throw_inc_hists[s].push_back(temp_inc_hist);

    for (int i = 1; i <= temp_xsec_hist->GetNbinsX(); ++i) {
      temp_xsec_hist->SetBinContent(
          i, throw_hists[s].back()->GetBinContent(i+1));
    }
    temp_xsec_hist->Divide(temp_inc_hist);
    throw_xsec_hists[s].push_back(temp_xsec_hist);
  }
}

void protoana::AbsCexDriver::PlotThrows(
    ThinSliceDataSet & data_set, std::map<int, std::vector<TH1*>> & throw_hists,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    size_t nThrows,
    std::map<int, std::vector<TH1*>> & truth_throw_hists,
    std::map<int, std::vector<TH1*>> & truth_inc_hists,
    std::map<int, std::vector<TH1*>> & truth_xsec_hists,
    std::map<int, TH1*> & best_fit_incs,
    std::map<int, TH1*> & best_fit_xsecs,
    std::map<int, TH1*> & nominal_incs,
    std::map<int, TH1*> & nominal_xsecs,
    TFile & output_file, bool plot_rebinned,
    std::map<int, std::vector<double>> * sample_scales) {
  std::map<int, TH1*> data_hists
      = (plot_rebinned ?
         data_set.GetRebinnedSelectionHists() :
         data_set.GetSelectionHists());

  //Build best fit hists and get bins for covariance 
  std::map<int, TH1D*> best_fit_selection_hists;
  int nBins = 0;
  for (auto it = data_hists.begin(); it != data_hists.end(); ++it ) {
    TH1D * best_fit_hist = (TH1D*)it->second->Clone();
    best_fit_hist->Reset();
    for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
      for (size_t i = 0; i < it2->second.size(); ++i) {
        for (size_t j = 0; j < it2->second[i].size(); ++j) {
          it2->second[i][j].SetFactorToBestFit();
          best_fit_hist->Add(
              (TH1D*)(plot_rebinned ?
                      it2->second[i][j].GetRebinnedSelectionHist(it->first) :
                      it2->second[i][j].GetSelectionHist(it->first)));
        }
      }
    }
    best_fit_selection_hists[it->first] = best_fit_hist;
    nBins += best_fit_hist->GetNbinsX();
  }

  TH2D selection_cov("SelectionCov", "", nBins, 0, nBins, nBins, 0, nBins);

  nBins = 0;
  std::map<int, size_t> sample_bins;
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    nBins += it->second[0].size();
    sample_bins[it->first] = it->second[0].size();
  }

  std::map<int, std::vector<double>> best_fit_truth;
  std::map<int, std::vector<double>> best_fit_errs;

  for (auto it = samples.begin(); it != samples.end(); ++it) {
    best_fit_truth[it->first]
        = std::vector<double>(sample_bins[it->first], 0.);
    best_fit_errs[it->first]
        = std::vector<double>(sample_bins[it->first], 0.);
   
    for (size_t i = 0; i < sample_bins[it->first]; ++i) {
      double best_fit_val_i = 0.;
      for (size_t j = 0; j < it->second.size(); ++j) {
        best_fit_val_i += it->second[j][i].GetVariedFlux();
      }

      best_fit_truth[it->first][i] = best_fit_val_i;
    }
  }

  TH2D interaction_cov("interaction_cov", "", nBins, 0, nBins, nBins, 0, nBins);
  std::map<int, std::vector<double>> best_fit_inc_truth;
  std::map<int, std::vector<double>> best_fit_xsec_truth;
  std::map<int, std::vector<double>> best_fit_inc_errs;
  std::map<int, std::vector<double>> best_fit_xsec_errs;

  nBins = 0;
  std::map<int, size_t> xsec_bins;
  for (auto it = best_fit_incs.begin(); it != best_fit_incs.end(); ++it) {
    int s = it->first;
    nBins += it->second->GetNbinsX();
    xsec_bins[s] = it->second->GetNbinsX();

    best_fit_inc_truth[s] = std::vector<double>(xsec_bins[s], 0.);
    best_fit_xsec_truth[s] = std::vector<double>(xsec_bins[s], 0.);
    best_fit_inc_errs[s] = std::vector<double>(xsec_bins[s], 0.);
    best_fit_xsec_errs[s] = std::vector<double>(xsec_bins[s], 0.);
    
    for (size_t i = 0; i < xsec_bins[s]; ++i) {
      best_fit_inc_truth[s][i] = it->second->GetBinContent(i+1);
      best_fit_xsec_truth[s][i] = best_fit_xsecs[s]->GetBinContent(i+1);
    }
  }

  //TH2D incident_cov("incident_cov", "", nBins, 0, nBins, nBins, 0, nBins);
  TH2D xsec_cov("xsec_cov", "", nBins, 0, nBins, nBins, 0, nBins);

  for (size_t z = 0; z < nThrows; ++z) {
    int bin_i = 1;
    for (auto it = best_fit_selection_hists.begin();
         it != best_fit_selection_hists.end(); ++it) {
      TH1D * best_fit = it->second;
      int selection_ID = it->first;
      std::vector<TH1*> & temp_throws = throw_hists[selection_ID];
      for (int i = 1; i <= best_fit->GetNbinsX(); ++i) {
        double best_fit_val_i = best_fit->GetBinContent(i);
        int bin_j = 1;
        for (auto it2 = best_fit_selection_hists.begin();
             it2 != best_fit_selection_hists.end(); ++it2) {

          TH1D * best_fit_2 = it2->second;
          int selection_ID_2 = it2->first;
          std::vector<TH1*> & temp_throws_2 = throw_hists[selection_ID_2];
          for (int j = 1; j <= best_fit_2->GetNbinsX(); ++j) {
            double best_fit_val_j = best_fit_2->GetBinContent(j);
            double val = (best_fit_val_i - temp_throws[z]->GetBinContent(i))*
                         (best_fit_val_j - temp_throws_2[z]->GetBinContent(j));
            selection_cov.SetBinContent(
                bin_i, bin_j, (val/temp_throws.size() +
                               selection_cov.GetBinContent(bin_i, bin_j)));
            ++bin_j;
          }
        }
        ++bin_i;
      }
    }

    bin_i = 1;
    for (auto it = samples.begin(); it != samples.end(); ++it) {
      std::vector<TH1 *> throw_hists_i = truth_throw_hists[it->first];
     
      for (size_t i = 0; i < sample_bins[it->first]; ++i) {
        double best_fit_val_i = best_fit_truth[it->first][i];

        int bin_j = 1;
        for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
          std::vector<TH1 *> throw_hists_j = truth_throw_hists[it2->first];
          for (size_t j = 0; j < sample_bins[it2->first]; ++j) {
            double best_fit_val_j = best_fit_truth[it2->first][j];

            double val
                = (throw_hists_i[z]->GetBinContent(i+1) - best_fit_val_i)*
                  (throw_hists_j[z]->GetBinContent(j+1) - best_fit_val_j);
            interaction_cov.SetBinContent(
                bin_i, bin_j,
                (interaction_cov.GetBinContent(bin_i, bin_j) +
                 val/throw_hists_i.size()));
            if (bin_i == bin_j && (z == nThrows - 1)) {
              best_fit_errs[it->first][i]
                  = sqrt(interaction_cov.GetBinContent(bin_i, bin_j));
            }
            ++bin_j;
          }
        }

        ++bin_i;
      }
    }

    bin_i = 1;
    for (auto it = truth_inc_hists.begin(); it != truth_inc_hists.end(); ++it) {
      //std::vector<TH1 *> inc_hists_i = it->second;
      std::vector<TH1 *> xsec_hists_i = truth_xsec_hists[it->first];

      for (size_t i = 0; i < xsec_bins[it->first]; ++i) {
        //double best_fit_inc_i = best_fit_inc_truth[it->first][i];
        double best_fit_xsec_i = best_fit_xsec_truth[it->first][i];

        int bin_j = 1;
        for (auto it2 = truth_inc_hists.begin(); it2 != truth_inc_hists.end();
             ++it2) {
          std::vector<TH1 *> xsec_hists_j = truth_xsec_hists[it2->first];
          for (size_t j = 0; j < xsec_bins[it2->first]; ++j) {
            double best_fit_xsec_j = best_fit_xsec_truth[it2->first][j];

            double val
                = (xsec_hists_i[z]->GetBinContent(i+1) - best_fit_xsec_i)*
                  (xsec_hists_j[z]->GetBinContent(j+1) - best_fit_xsec_j);
            xsec_cov.SetBinContent(
                bin_i, bin_j,
                (xsec_cov.GetBinContent(bin_i, bin_j) +
                 val/nThrows));
            if (bin_i == bin_j && (z == nThrows - 1)) {
              best_fit_xsec_errs[it->first][i]
                  = sqrt(xsec_cov.GetBinContent(bin_i, bin_j));
            }
          }
        }
      }
    }
  }


  output_file.cd("Throws");
  selection_cov.Write();
  interaction_cov.Write();
  xsec_cov.Write();

  int bin_count = 0;
  for (auto it = data_hists.begin(); it != data_hists.end(); ++it) {
    int selection_ID = it->first;
    std::vector<TH1*> hists = throw_hists.at(selection_ID);

    std::string canvas_name = "cThrow" +
                              data_set.GetSelectionName(selection_ID);
    TCanvas cThrow(canvas_name.c_str(), "");
    cThrow.SetTicks();

    std::string name = "Throw" + data_set.GetSelectionName(selection_ID);
    auto data_hist = it->second;
    std::vector<double> xs, xs_width;
    std::vector<double> ys, errs;
    for (int i = 1;
         i <= best_fit_selection_hists[it->first]->GetNbinsX(); ++i) {
      ys.push_back(
          best_fit_selection_hists[it->first]->GetBinContent(i));
      errs.push_back(
          sqrt(selection_cov.GetBinContent(bin_count+i, bin_count+i)));
      xs.push_back(data_hist->GetBinCenter(i));
      xs_width.push_back(data_hist->GetBinWidth(i)/2.);
    } 

    TGraphAsymmErrors throw_gr(data_hist->GetNbinsX(),
                               &xs[0], &ys[0], 
                               &xs_width[0], &xs_width[0], &errs[0], &errs[0]);

    throw_gr.SetFillStyle(3144);
    throw_gr.SetFillColor(kRed);
    throw_gr.Draw("a2");
    data_hist->Draw("same e1");
    output_file.cd("Throws");
    cThrow.Write();

    bin_count += data_hist->GetNbinsX();
  }

  bin_count = 0;
  for (auto it = truth_throw_hists.begin(); it != truth_throw_hists.end(); ++it) {
    int sample_ID = it->first;

    std::vector<double> xs, xs_width;
    for (size_t i = 0; i < sample_bins[it->first]; ++i) {
      xs.push_back(i + 0.5);
      xs_width.push_back(.5);
    }

    std::string name = "hNominal" + samples[sample_ID][0][0].GetName();
    TH1D temp_nominal(name.c_str(), "", xs.size(), 0, xs.size());
    std::vector<std::vector<ThinSliceSample>> & samples_vec_2D
        = samples[sample_ID];
    for (size_t i = 0; i < samples_vec_2D.size(); ++i) {
      for (size_t j = 0; j < samples_vec_2D[i].size(); ++j) {
        temp_nominal.AddBinContent(j+1, samples_vec_2D[i][j].GetNominalFlux());
      }
    }

    double max = -999.;
    for (size_t i = 0; i < sample_bins[it->first]; ++i) {
      if ((best_fit_truth[sample_ID][i] + best_fit_errs[sample_ID][i]) > max)
        max = (best_fit_truth[sample_ID][i] + best_fit_errs[sample_ID][i]);

      if (temp_nominal.GetBinContent(i+1) > max)
        max = temp_nominal.GetBinContent(i+1);
    }

    output_file.cd("Throws");
    std::string canvas_name = "cTruthThrow" + samples[sample_ID][0][0].GetName();
    TCanvas cThrow(canvas_name.c_str(), "");
    cThrow.SetTicks();
    TGraphAsymmErrors throw_gr(xs.size(),
                                &xs[0], &best_fit_truth[it->first][0], 
                                &xs_width[0], &xs_width[0],
                                &best_fit_errs[it->first][0],
                                &best_fit_errs[it->first][0]);
    throw_gr.SetFillStyle(3144);
    throw_gr.SetFillColor(kRed);
    throw_gr.SetMinimum(0.);
    throw_gr.SetMaximum(1.5*max);
    throw_gr.Draw("a2");
    throw_gr.Draw("p");

    temp_nominal.SetMarkerColor(kBlue);
    temp_nominal.SetMarkerStyle(20);
    temp_nominal.Draw("same p");

    TLegend leg;
    leg.AddEntry(&throw_gr, "Throws", "lpf");
    leg.AddEntry(&temp_nominal, "Nominal", "p");

    if (sample_scales) {
      name = "hVaried" + samples[sample_ID][0][0].GetName();
      TH1D * temp_varied = (TH1D*)temp_nominal.Clone(name.c_str());
      for (size_t i = 0; i < xs.size(); ++i) {
        temp_varied->SetBinContent(
            i+1, temp_varied->GetBinContent(i+1)*(*sample_scales)[sample_ID][i]);
      }
      temp_varied->SetMarkerColor(kBlack);
      temp_varied->SetMarkerStyle(20);
      temp_varied->Draw("same p");
      leg.AddEntry(temp_varied, "Fake Data", "p");
    }

    leg.Draw();
    cThrow.Write();

    bin_count += xs.size();
  }

  for (auto it = best_fit_xsec_truth.begin(); it != best_fit_xsec_truth.end();
       ++it) {
    int sample_ID = it->first;

    std::vector<double> xs, xs_width;
    for (size_t i = 0; i < xsec_bins[sample_ID]; ++i) {
      xs.push_back(i + 0.5);
      xs_width.push_back(.5);
    }

    output_file.cd("Throws");
    std::string canvas_name = "cXSecThrow" + samples[sample_ID][0][0].GetName();
    TCanvas cThrow(canvas_name.c_str(), "");
    cThrow.SetTicks();
    TGraphAsymmErrors throw_gr(xs.size(),
                               &xs[0], &best_fit_xsec_truth[it->first][0], 
                               &xs_width[0], &xs_width[0],
                               &best_fit_xsec_errs[it->first][0],
                               &best_fit_xsec_errs[it->first][0]);
    throw_gr.SetFillStyle(3144);
    throw_gr.SetFillColor(kRed);
    throw_gr.SetMinimum(0.);
    throw_gr.Draw("a2");
    throw_gr.Draw("p");

    cThrow.Write();
  }
}
