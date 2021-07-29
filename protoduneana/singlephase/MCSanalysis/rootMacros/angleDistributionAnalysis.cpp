#include "HelperFunctions.cpp"
#include "PlottingFunctions.cpp"
#include "AngleFittingFunctions.cpp"

// Root MACRO Function
void angleDistributionAnalysis() {

  // TODO: Set this to kTRUE
  gROOT->SetBatch(kFALSE);

  /////////////////////////////////////////////////////////////////////////////////////
  // Get plots from root file
  /////////////////////////////////////////////////////////////////////////////////////

  gSystem->Exec("rm -r ADA_plots");
  gSystem->Exec("mkdir ADA_plots");

  // Get Root File
  //const char* path = "/dune/app/users/hmeyer5/curvyTrackFix2/srcs/dunetpc/dune/LSU/500_mu_0.5-4.5_GeV_start_beam-entry_dir_beam-dir_standardReco.root";
  const char* path = "../MCSAngleAnalysis_hist.root";
  TFile* file = new TFile(path);

  // TODO: Add polygonal 3D angles
  // TODO: Verify that all plots are made in the exact same way.
  // TODO: Fit the 3D angles using a different equation than the 2D angles, use factor that SHOULD fit to sqrt(2)*sqrt(1-2/pi) or something like that...
  // TODO: Formally compare 3D and 2D angles over a variety of new metrics.
  // TODO: Verify angle distributions.  Are they Gaussian/Half-Gaussian? How to quantify this?
  // TODO: Extract similar code to be done in separate functions, similar to what is done in momentumAnalysis.cpp?

  // Get each histogram from the root file
  // True Linear 3D Plots
  TH2F* trueLinear_theta3DVsSegmentMomentum_HIST = Get2DHist(file, "trueLinear_theta3DVsSegmentMomentum_HIST");
  TH2F* trueLinear_thetaXZprimeVsSegmentMomentum_HIST = Get2DHist(file, "trueLinear_thetaXZprimeVsSegmentMomentum_HIST");
  TH2F* trueLinear_thetaYZprimeVsSegmentMomentum_HIST = Get2DHist(file, "trueLinear_thetaYZprimeVsSegmentMomentum_HIST");
  TH2F* truePolygonal_theta3DVsSegmentMomentum_HIST = Get2DHist(file, "truePolygonal_theta3DVsSegmentMomentum_HIST");
  TH2F* truePolygonal_thetaXZprimeVsSegmentMomentum_HIST = Get2DHist(file, "truePolygonal_thetaXZprimeVsSegmentMomentum_HIST");
  TH2F* truePolygonal_thetaYZprimeVsSegmentMomentum_HIST = Get2DHist(file, "truePolygonal_thetaYZprimeVsSegmentMomentum_HIST");

  TH2F* recoLinear_theta3DVsSegmentMomentum_HIST = Get2DHist(file, "recoLinear_theta3DVsSegmentMomentum_HIST");
  TH2F* recoLinear_thetaXZprimeVsSegmentMomentum_HIST = Get2DHist(file, "recoLinear_thetaXZprimeVsSegmentMomentum_HIST");
  TH2F* recoLinear_thetaYZprimeVsSegmentMomentum_HIST = Get2DHist(file, "recoLinear_thetaYZprimeVsSegmentMomentum_HIST");
  TH2F* recoPolygonal_theta3DVsSegmentMomentum_HIST = Get2DHist(file, "recoPolygonal_theta3DVsSegmentMomentum_HIST");
  TH2F* recoPolygonal_thetaXZprimeVsSegmentMomentum_HIST = Get2DHist(file, "recoPolygonal_thetaXZprimeVsSegmentMomentum_HIST");
  TH2F* recoPolygonal_thetaYZprimeVsSegmentMomentum_HIST = Get2DHist(file, "recoPolygonal_thetaYZprimeVsSegmentMomentum_HIST");

  /////////////////////////////////////////////////////////////////////////////////////
  // SECTION 1: Create/modify plots
  /////////////////////////////////////////////////////////////////////////////////////

  //
  // Subsection 1.1: Rebin each Histogram using this rebinFactor
  //
  const Int_t rebinFactor = 25;
  trueLinear_theta3DVsSegmentMomentum_HIST->RebinX(rebinFactor);
  trueLinear_thetaXZprimeVsSegmentMomentum_HIST->RebinX(rebinFactor);
  trueLinear_thetaYZprimeVsSegmentMomentum_HIST->RebinX(rebinFactor);
  truePolygonal_theta3DVsSegmentMomentum_HIST->RebinX(rebinFactor);
  truePolygonal_thetaXZprimeVsSegmentMomentum_HIST->RebinX(rebinFactor);
  truePolygonal_thetaYZprimeVsSegmentMomentum_HIST->RebinX(rebinFactor);

  recoLinear_theta3DVsSegmentMomentum_HIST->RebinX(rebinFactor);
  recoLinear_thetaXZprimeVsSegmentMomentum_HIST->RebinX(rebinFactor);
  recoLinear_thetaYZprimeVsSegmentMomentum_HIST->RebinX(rebinFactor);
  recoPolygonal_theta3DVsSegmentMomentum_HIST->RebinX(rebinFactor);
  recoPolygonal_thetaXZprimeVsSegmentMomentum_HIST->RebinX(rebinFactor);
  recoPolygonal_thetaYZprimeVsSegmentMomentum_HIST->RebinX(rebinFactor);
  
  //
  // Subsection 1.2: Combine XZ & YZ histograms
  //
  TH2F* trueLinear_thetaPrimeVsSegmentMomentum_HIST = CombineHistograms(trueLinear_thetaXZprimeVsSegmentMomentum_HIST, trueLinear_thetaYZprimeVsSegmentMomentum_HIST, "trueLinear_thetaPrimeVsSegmentMomentum_HIST", "True Linear #theta' Vs. Segment Momentum; Segment Momentum (GeV/c); #theta' (mrad)");
  TH2F* truePolygonal_thetaPrimeVsSegmentMomentum_HIST = CombineHistograms(truePolygonal_thetaXZprimeVsSegmentMomentum_HIST, truePolygonal_thetaYZprimeVsSegmentMomentum_HIST, "truePolygonal_thetaPrimeVsSegmentMomentum_HIST", "True Polygonal #theta' Vs. Segment Momentum; Segment Momentum (GeV/c); #theta' (mrad)");

  TH2F* recoLinear_thetaPrimeVsSegmentMomentum_HIST = CombineHistograms(recoLinear_thetaXZprimeVsSegmentMomentum_HIST, recoLinear_thetaYZprimeVsSegmentMomentum_HIST, "recoLinear_thetaPrimeVsSegmentMomentum_HIST", "Reco Linear #theta' Vs. Segment Momentum; Segment Momentum (GeV/c); #theta' (mrad)");
  TH2F* recoPolygonal_thetaPrimeVsSegmentMomentum_HIST = CombineHistograms(recoPolygonal_thetaXZprimeVsSegmentMomentum_HIST, recoPolygonal_thetaYZprimeVsSegmentMomentum_HIST, "recoPolygonal_thetaPrimeVsSegmentMomentum_HIST", "Reco Polygonal #theta' Vs. Segment Momentum; Segment Momentum (GeV/c); #theta' (mrad)");

  //
  // Subsection 1.3: Get Graph of sigmaHL/sigmaRMS vs. segment momentum for 2D angles. (with y-uncertainties)
  //

  // Subsubsection 1.3.1: Ideal future, if there is no direction discrepancy. (for reco, this should already be the case for true)
  TGraphErrors* trueLinear_sigmaHLVsSegmentMomentum = GetSigmaGraph(trueLinear_thetaPrimeVsSegmentMomentum_HIST, "trueLinear_sigmaHLVsSegmentMomentum", "True Linear #sigma_{HL} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{HL} (mrad)");
  TGraphErrors* truePolygonal_sigmaHLVsSegmentMomentum = GetSigmaGraph(truePolygonal_thetaPrimeVsSegmentMomentum_HIST, "truePolygonal_sigmaHLVsSegmentMomentum", "True Polygonal #sigma_{HL} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{HL} (mrad)");

  TGraphErrors* recoLinear_sigmaRMSVsSegmentMomentum = GetSigmaGraph(recoLinear_thetaPrimeVsSegmentMomentum_HIST, "recoLinear_sigmaRMSVsSegmentMomentum", "Reco Linear #sigma_{RMS} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{RMS} (mrad)");
  TGraphErrors* recoPolygonal_sigmaRMSVsSegmentMomentum = GetSigmaGraph(recoPolygonal_thetaPrimeVsSegmentMomentum_HIST, "recoPolygonal_sigmaRMSVsSegmentMomentum", "Reco Polygonal #sigma_{RMS} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{RMS} (mrad)");

  // Subsubsection 1.3.2: Split into two projection directions, due to direction discrepancy.  Should only need to do reco, but we're also doing true to cover all bases.
  TGraphErrors* trueLinearXZ_sigmaHLVsSegmentMomentum = GetSigmaGraph(trueLinear_thetaXZprimeVsSegmentMomentum_HIST, "trueLinearXZ_sigmaHLVsSegmentMomentum", "True Linear XZ #sigma_{HL} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{HL} (mrad)");
  TGraphErrors* trueLinearYZ_sigmaHLVsSegmentMomentum = GetSigmaGraph(trueLinear_thetaYZprimeVsSegmentMomentum_HIST, "trueLinearYZ_sigmaHLVsSegmentMomentum", "True Linear YZ #sigma_{HL} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{HL} (mrad)");
  TGraphErrors* truePolygonalXZ_sigmaHLVsSegmentMomentum = GetSigmaGraph(truePolygonal_thetaXZprimeVsSegmentMomentum_HIST, "truePolygonalXZ_sigmaHLVsSegmentMomentum", "True Polygonal XZ #sigma_{HL} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{HL} (mrad)");
  TGraphErrors* truePolygonalYZ_sigmaHLVsSegmentMomentum = GetSigmaGraph(truePolygonal_thetaYZprimeVsSegmentMomentum_HIST, "truePolygonalYZ_sigmaHLVsSegmentMomentum", "True Polgyonal YZ #sigma_{HL} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{HL} (mrad)");

  TGraphErrors* recoLinearXZ_sigmaRMSVsSegmentMomentum = GetSigmaGraph(recoLinear_thetaXZprimeVsSegmentMomentum_HIST, "recoLinearXZ_sigmaRMSVsSegmentMomentum", "Reco Linear XZ #sigma_{RMS} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{RMS} (mrad)");
  TGraphErrors* recoLinearYZ_sigmaRMSVsSegmentMomentum = GetSigmaGraph(recoLinear_thetaYZprimeVsSegmentMomentum_HIST, "recoLinearXZ_sigmaRMSVsSegmentMomentum", "Reco Linear YZ #sigma_{RMS} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{RMS} (mrad)");
  TGraphErrors* recoPolygonalXZ_sigmaRMSVsSegmentMomentum = GetSigmaGraph(recoPolygonal_thetaXZprimeVsSegmentMomentum_HIST, "recoPolygonalXZ_sigmaRMSVsSegmentMomentum", "Reco Polygonal XZ #sigma_{RMS} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{RMS} (mrad)");
  TGraphErrors* recoPolygonalYZ_sigmaRMSVsSegmentMomentum = GetSigmaGraph(recoPolygonal_thetaYZprimeVsSegmentMomentum_HIST, "recoPolygonalXZ_sigmaRMSVsSegmentMomentum", "Reco Polygonal YZ #sigma_{RMS} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{RMS} (mrad)");

  // Subsubsection 1.3.3: Get Graph of sigma vs. segment momentum for 3D angles. (with y-uncertainties)
  // TODO: Changed from GetRawSigmaGaph to GetSigmaGraph, should we changed back?
  TGraphErrors* trueLinear3D_sigmaVsSegmentMomentum = GetSigmaGraph(trueLinear_theta3DVsSegmentMomentum_HIST, "trueLinear3D_sigmaVsSegmentMomentum", "True Linear 3D #sigma Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma (mrad)");
  TGraphErrors* truePolygonal3D_sigmaVsSegmentMomentum = GetSigmaGraph(truePolygonal_theta3DVsSegmentMomentum_HIST, "truePolygonal3D_sigmaVsSegmentMomentum", "True Polygonal 3D #sigma Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma (mrad)");

  TGraphErrors* recoLinear3D_sigmaVsSegmentMomentum = GetSigmaGraph(recoLinear_theta3DVsSegmentMomentum_HIST, "recoLinear3D_sigmaVsSegmentMomentum", "Reco Linear 3D #sigma Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma (mrad)");
  TGraphErrors* recoPolygonal3D_sigmaVsSegmentMomentum = GetSigmaGraph(recoPolygonal_theta3DVsSegmentMomentum_HIST, "recoPolygonal3D_sigmaVsSegmentMomentum", "Reco Polygonal 3D #sigma Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma (mrad)");

  //
  // Subsection 1.4: Get Ratio of 3D and 2D angle distributions.
  // TODO: Add graph titles, use prefixes like I did in momentumAnalysis.cpp
  //

  // Subsubsection 1.4.1: Ideal future, if there is no direction discrepancy. (for reeco, this hsould already be the case for true)
  TGraph* trueLinear_sigmaRatio = GetSigmaRatioGraph(trueLinear3D_sigmaVsSegmentMomentum, trueLinear_sigmaHLVsSegmentMomentum, "trueLinear_sigmaRatio", "");
  TGraph* truePolygonal_sigmaRatio = GetSigmaRatioGraph(trueLinear3D_sigmaVsSegmentMomentum, truePolygonal_sigmaHLVsSegmentMomentum, "truePolygonal_sigmaRatio", "");

  TGraph* recoLinear_sigmaRatio = GetSigmaRatioGraph(recoLinear3D_sigmaVsSegmentMomentum, recoLinear_sigmaRMSVsSegmentMomentum, "recoLinear_sigmaRatio", "");
  TGraph* recoPolygonal_sigmaRatio = GetSigmaRatioGraph(recoPolygonal3D_sigmaVsSegmentMomentum, recoPolygonal_sigmaRMSVsSegmentMomentum, "recoPolygonal_sigmaRatio", "");
 
  // Subsubsection 1.4.2: Split into two projection directions, due to direction discrepancy.  Should only need to do reco, but we're also diong true to cover all bases.
  TGraph* trueLinearXZ_sigmaRatio = GetSigmaRatioGraph(trueLinear3D_sigmaVsSegmentMomentum, trueLinearXZ_sigmaHLVsSegmentMomentum, "trueLinearXZ_sigmaRatio", "");
  TGraph* trueLinearYZ_sigmaRatio = GetSigmaRatioGraph(trueLinear3D_sigmaVsSegmentMomentum, trueLinearYZ_sigmaHLVsSegmentMomentum, "trueLinearYZ_sigmaRatio", "");
  TGraph* truePolygonalXZ_sigmaRatio = GetSigmaRatioGraph(truePolygonal3D_sigmaVsSegmentMomentum, truePolygonalXZ_sigmaHLVsSegmentMomentum, "truePolygonalXZ_sigmaRatio", "");
  TGraph* truePolygonalYZ_sigmaRatio = GetSigmaRatioGraph(truePolygonal3D_sigmaVsSegmentMomentum, truePolygonalYZ_sigmaHLVsSegmentMomentum, "truePolygonalYZ_sigmaRatio", "");

  TGraph* recoLinearXZ_sigmaRatio = GetSigmaRatioGraph(recoLinear3D_sigmaVsSegmentMomentum, recoLinearXZ_sigmaRMSVsSegmentMomentum, "recoLinearXZ_sigmaRatio", "");
  TGraph* recoLinearYZ_sigmaRatio = GetSigmaRatioGraph(recoLinear3D_sigmaVsSegmentMomentum, recoLinearYZ_sigmaRMSVsSegmentMomentum, "recoLinearYZ_sigmaRatio", "");
  TGraph* recoPolygonalXZ_sigmaRatio = GetSigmaRatioGraph(recoLinear3D_sigmaVsSegmentMomentum, recoPolygonalXZ_sigmaRMSVsSegmentMomentum, "recoPolygonalXZ_sigmaRatio", "");
  TGraph* recoPolygonalYZ_sigmaRatio = GetSigmaRatioGraph(recoLinear3D_sigmaVsSegmentMomentum, recoPolygonalYZ_sigmaRMSVsSegmentMomentum, "recoPolygonalYZ_sigmaRatio", "");

  // TODO: Delete this or make use of it in some other way:
  // // print out comparison between projected sigmas added in quadrature with raw sigma
  // for(Int_t point = 0; point < trueLinear3D_sigmaVsSegmentMomentum->GetN(); ++point) {
  //   std::cout << "mom = " << trueLinear3D_sigmaVsSegmentMomentum->GetPointX(point) << " GeV/c" << "\t"
  // 	      << "rawSigma = " << trueLinear3D_sigmaVsSegmentMomentum->GetPointY(point) << " mrad" << "\t"
  // 	      << "quad = " << pow(pow(trueLinearXZ_sigmaHLVsSegmentMomentum->GetPointY(point),2) + pow(trueLinearYZ_sigmaHLVsSegmentMomentum->GetPointY(point),2),0.5) << " mrad" << "\t"
  // 	      << "scaled proj = " << trueLinear_sigmaHLVsSegmentMomentum->GetPointY(point) / sqrt(2) << " mrad" << "\t"
  // 	      << std::endl;
  // }

  // // Scale raw angle sigma
  // for(Int_t point = 0; point < trueLinear3D_sigmaVsSegmentMomentum->GetN(); ++point) {
  //   Float_t segmentMomentum = trueLinear3D_sigmaVsSegmentMomentum->GetPointX(point);
  //   Float_t rawSigma = trueLinear3D_sigmaVsSegmentMomentum->GetPointY(point);
  //   Float_t rawSigma_HalfGaussScaled = rawSigma / sqrt(1-2/TMath::Pi()); // Scaled from half-normal distribution
  //   Float_t rawSigma_projScaled = rawSigma_HalfGaussScaled / sqrt(2); // At this point, the value should be the same as the projected angle distribution.
  //   Float_t projectedSigma = trueLinear_sigmaHLVsSegmentMomentum->GetPointY(point);
  //   trueLinear3D_sigmaVsSegmentMomentum->SetPoint(point, segmentMomentum, rawSigma_projScaled);

  //   std::cout << "mom = " << segmentMomentum << " GeV/c" << "\t"
  // 	      << "rawSigma = " << rawSigma << " mrad" << "\t"
  // 	      << "rawSigma_HGS = " << rawSigma_HalfGaussScaled << " mrad" << "\t"
  // 	      << "rawSigma_PS = " << rawSigma_projScaled << " mrad" << "\t"
  // 	      << "projSigma = " << projectedSigma << " mrad" << "\t"
  // 	      << std::endl;
  // }

  std::cout << "trueLinear_sigmaRatio_mean = " << trueLinear_sigmaRatio->GetMean(2) << std::endl;
  std::cout << "truePolygonal_sigmaRatio_mean = " << truePolygonal_sigmaRatio->GetMean(2) << std::endl;
  std::cout << "recoLinear_sigmaRatio_mean = " << recoLinear_sigmaRatio->GetMean(2) << std::endl;
  std::cout << "recoPolygonal_sigmaRatio_mean = " << recoPolygonal_sigmaRatio->GetMean(2) << std::endl;

  std::cout << "trueLinearXZ_sigmaRatio_mean = " << trueLinearXZ_sigmaRatio->GetMean(2) << std::endl;
  std::cout << "trueLinearYZ_sigmaRatio_mean = " << trueLinearYZ_sigmaRatio->GetMean(2) << std::endl;
  std::cout << "truePolygonalXZ_sigmaRatio_mean = " << truePolygonalXZ_sigmaRatio->GetMean(2) << std::endl;
  std::cout << "truePolygonalYZ_sigmaRatio_mean = " << truePolygonalYZ_sigmaRatio->GetMean(2) << std::endl;

  std::cout << "recoLinearXZ_sigmaRatio_mean = " << recoLinearXZ_sigmaRatio->GetMean(2) << std::endl;
  std::cout << "recoLinearYZ_sigmaRatio_mean = " << recoLinearYZ_sigmaRatio->GetMean(2) << std::endl;
  std::cout << "recoPolygonalXZ_sigmaRatio_mean = " << recoPolygonalXZ_sigmaRatio->GetMean(2) << std::endl;
  std::cout << "recoPolygonalYZ_sigmaRatio_mean = " << recoPolygonalYZ_sigmaRatio->GetMean(2) << std::endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // Section 2: Fit Plots to get sigmaHL, sigmaRES, w_0, parameters.
  /////////////////////////////////////////////////////////////////////////////////////

  // TODO: Automatic range calculation?
  const Double_t fitMin = 0.2;
  const Double_t fitMax = 4.1;

  //
  // Subsection 2.1: True sigmaHL fits.
  //

  // Subsubsection 2.1.1: Fit sigma of combined angle directions.  These fit parameters are used in other fits.
  TF1* trueLinear_sigmaHL_function = GetSigmaHL_function("trueLinear_sigmaHL_function", fitMin, fitMax);
  TFitResultPtr trueLinear_sigmaHL_fit = trueLinear_sigmaHLVsSegmentMomentum->Fit(trueLinear_sigmaHL_function, "QS0");

  TF1* truePolygonal_sigmaHL_function = GetSigmaHL_function("truePolygonal_sigmaHL_function", fitMin, fitMax);
  TFitResultPtr truePolygonal_sigmaHL_fit = truePolygonal_sigmaHLVsSegmentMomentum->Fit(truePolygonal_sigmaHL_function, "QS0");

  // Subsubsection 2.1.2: Fit sigma of individual angle directions.  These fit parameters are *not* used in other fits.
  TF1* trueLinearXZ_sigmaHL_function = GetSigmaHL_function("trueLinearXZ_sigmaHL_function", fitMin, fitMax);
  TFitResultPtr trueLinearXZ_sigmaHL_fit = trueLinearXZ_sigmaHLVsSegmentMomentum->Fit(trueLinearXZ_sigmaHL_function, "QS0");

  TF1* trueLinearYZ_sigmaHL_function = GetSigmaHL_function("trueLinearYZ_sigmaHL_function", fitMin, fitMax);
  TFitResultPtr trueLinearYZ_sigmaHL_fit = trueLinearYZ_sigmaHLVsSegmentMomentum->Fit(trueLinearYZ_sigmaHL_function, "QS0");

  TF1* truePolygonalXZ_sigmaHL_function = GetSigmaHL_function("truePolygonalXZ_sigmaHL_function", fitMin, fitMax);
  TFitResultPtr truePolygonalXZ_sigmaHL_fit = truePolygonalXZ_sigmaHLVsSegmentMomentum->Fit(truePolygonalXZ_sigmaHL_function, "QS0");

  TF1* truePolygonalYZ_sigmaHL_function = GetSigmaHL_function("truePolygonalYZ_sigmaHL_function", fitMin, fitMax);
  TFitResultPtr truePolygonalYZ_sigmaHL_fit = truePolygonalYZ_sigmaHLVsSegmentMomentum->Fit(truePolygonalYZ_sigmaHL_function, "QS0");

  // Subsubsection 2.1.2: Define Fit Result Parameters
  Double_t trueLinear_kappa_a = trueLinear_sigmaHL_fit->Value(0);
  Double_t trueLinear_kappa_c = trueLinear_sigmaHL_fit->Value(1);

  Double_t truePolygonal_kappa_a = truePolygonal_sigmaHL_fit->Value(0);
  Double_t truePolygonal_kappa_c = truePolygonal_sigmaHL_fit->Value(1);

  //
  // Subsection 2.2: True sigma3D fits.
  //

  // Subsubsection 2.2.1: Fit true sigma3D.
  TF1* trueLinear3D_sigmaHL_function = GetSigmaHL_3D_function("trueLinear3D_sigmaHL_function", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c);
  TFitResultPtr trueLinear3D_sigmaHL_fit = trueLinear3D_sigmaVsSegmentMomentum->Fit(trueLinear3D_sigmaHL_function, "QS0");

  TF1* truePolygonal3D_sigmaHL_function = GetSigmaHL_3D_function("truePolygonal3D_sigmaHL_function", fitMin, fitMax, truePolygonal_kappa_a, truePolygonal_kappa_c);
  TFitResultPtr truePolygonal3D_sigmaHL_fit = truePolygonal3D_sigmaVsSegmentMomentum->Fit(truePolygonal3D_sigmaHL_function, "QS0");

  // Subsubsection 2.2.2: Define Fit Result Parameters
  Double_t trueLinear_gamma3D = trueLinear3D_sigmaHL_fit->Value(3);
  Double_t truePolygonal_gamma3D = trueLinear3D_sigmaHL_fit->Value(3);

  //
  // Subsection 2.3: Reco sigmaRMS fits.
  //
  
  // Subsubsection 2.3.1: Fit sigma of combined angle directions.
  TF1* recoLinear_sigmaRMS_function = GetSigmaRMS_function("recoLinear_sigmaRMS_function", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c);
  TFitResultPtr recoLinear_sigmaRMS_fit = recoLinear_sigmaRMSVsSegmentMomentum->Fit(recoLinear_sigmaRMS_function, "QS0");

  TF1* recoLinear_sigmaRMSraw_function = GetSigmaRMSraw_function("recoLinear_sigmaRMSraw_function", fitMin, fitMax);
  TFitResultPtr recoLinear_sigmaRMSraw_fit = recoLinear_sigmaRMSVsSegmentMomentum->Fit(recoLinear_sigmaRMSraw_function, "QS0");

  TF1* recoPolygonal_sigmaRMS_function = GetSigmaRMS_function("recoPolygonal_sigmaRMS_function", fitMin, fitMax, truePolygonal_kappa_a, truePolygonal_kappa_c);
  TFitResultPtr recoPolygonal_sigmaRMS_fit = recoPolygonal_sigmaRMSVsSegmentMomentum->Fit(recoPolygonal_sigmaRMS_function, "QS0");
  
  // Subsubsection 2.3.2: Fit sigma of individual angle directions.
  TF1* recoLinearXZ_sigmaRMS_function = GetSigmaRMS_function("recoLinearXZ_sigmaRMS_function", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c);
  TFitResultPtr recoLinearXZ_sigmaRMS_fit = recoLinearXZ_sigmaRMSVsSegmentMomentum->Fit(recoLinearXZ_sigmaRMS_function, "QS0");

  TF1* recoLinearYZ_sigmaRMS_function = GetSigmaRMS_function("recoLinearYZ_sigmaRMS_function", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c);
  TFitResultPtr recoLinearYZ_sigmaRMS_fit = recoLinearYZ_sigmaRMSVsSegmentMomentum->Fit(recoLinearYZ_sigmaRMS_function, "QS0");

  TF1* recoPolygonalXZ_sigmaRMS_function = GetSigmaRMS_function("recoPolygonalXZ_sigmaRMS_function", fitMin, fitMax, truePolygonal_kappa_a, truePolygonal_kappa_c);
  TFitResultPtr recoPolygonalXZ_sigmaRMS_fit = recoPolygonalXZ_sigmaRMSVsSegmentMomentum->Fit(recoPolygonalXZ_sigmaRMS_function, "QS0");

  TF1* recoPolygonalYZ_sigmaRMS_function = GetSigmaRMS_function("recoPolygonalYZ_sigmaRMS_function", fitMin, fitMax, truePolygonal_kappa_a, truePolygonal_kappa_c);
  TFitResultPtr recoPolygonalYZ_sigmaRMS_fit = recoPolygonalYZ_sigmaRMSVsSegmentMomentum->Fit(recoPolygonalYZ_sigmaRMS_function, "QS0");

  // Subsubsection 2.3.3: Define Fit Result Parameters
  
  Double_t recoLinear_sigmaRES = recoLinear_sigmaRMS_fit->Value(2);
  Double_t recoPolygonal_sigmaRES = recoPolygonal_sigmaRMS_fit->Value(2);

  Double_t recoLinearXZ_sigmaRES = recoLinearXZ_sigmaRMS_fit->Value(2);
  Double_t recoLinearYZ_sigmaRES = recoLinearYZ_sigmaRMS_fit->Value(2);

  Double_t recoPolygonalXZ_sigmaRES = recoPolygonalXZ_sigmaRMS_fit->Value(2);
  Double_t recoPolygonalYZ_sigmaRES = recoPolygonalYZ_sigmaRMS_fit->Value(2);

  //
  // Subsection 2.4: Reco sigma3D fits.
  //

  // Subsubsection 2.4.1: Fit reco sigma3D
  TF1* recoLinear3D_sigmaRMS_function = GetSigmaRMS_3D_function("recoLinear3D_sigmaRMS_function", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c);
  TFitResultPtr recoLinear3D_sigmaRMS_fit = recoLinear3D_sigmaVsSegmentMomentum->Fit(recoLinear3D_sigmaRMS_function, "QS0");

  TF1* recoPolygonal3D_sigmaRMS_function = GetSigmaRMS_3D_function("recoPolygonal3D_sigmaRMS_function", fitMin, fitMax, truePolygonal_kappa_a, truePolygonal_kappa_c);
  TFitResultPtr recoPolygonal3D_sigmaRMS_fit = recoPolygonal3D_sigmaVsSegmentMomentum->Fit(recoPolygonal3D_sigmaRMS_function, "QS0");

  // Subsubsection 2.4.2: Fit reco sigma3D, fixing either sigmaRES, to see if gamma3D is consistent

  TF1* recoLinear3D_sigmaRMS_function_fixSigmaRES = GetSigmaRMS_3D_function_fixSigmaRES("recoLinear3D_sigmaRMS_function_fixSigmaRES", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c, recoLinear_sigmaRES);
  TFitResultPtr recoLinear3D_sigmaRMS_fit_fixSigmaRES = recoLinear3D_sigmaVsSegmentMomentum->Fit(recoLinear3D_sigmaRMS_function_fixSigmaRES, "QS0");

  TF1* recoPolygonal3D_sigmaRMS_function_fixSigmaRES = GetSigmaRMS_3D_function_fixSigmaRES("recoPolygonal3D_sigmaRMS_function_fixSigmaRES", fitMin, fitMax, truePolygonal_kappa_a, truePolygonal_kappa_c, recoPolygonal_sigmaRES);
  TFitResultPtr recoPolygonal3D_sigmaRMS_fit_fixSigmaRES = recoPolygonal3D_sigmaVsSegmentMomentum->Fit(recoPolygonal3D_sigmaRMS_function_fixSigmaRES, "QS0");

  // Subsubsection 2.4.3: Fit reco sigma3D, fixing gamma3D, to see if sigmaRES is consistent

  TF1* recoLinear3D_sigmaRMS_function_fixGamma3D = GetSigmaRMS_3D_function_fixGamma3D("recoLinear3D_sigmaRMS_function_fixGamma3D", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c, trueLinear_gamma3D);
  TFitResultPtr recoLinear3D_sigmaRMS_fit_fixGamma3D = recoLinear3D_sigmaVsSegmentMomentum->Fit(recoLinear3D_sigmaRMS_function_fixGamma3D, "QS0");

  TF1* recoPolygonal3D_sigmaRMS_function_fixGamma3D = GetSigmaRMS_3D_function_fixGamma3D("recoPolygonal3D_sigmaRMS_function_fixGamma3D", fitMin, fitMax, truePolygonal_kappa_a, truePolygonal_kappa_c, truePolygonal_gamma3D);
  TFitResultPtr recoPolygonal3D_sigmaRMS_fit_fixGamma3D = recoPolygonal3D_sigmaVsSegmentMomentum->Fit(recoPolygonal3D_sigmaRMS_function_fixGamma3D, "QS0");

  // Subsubsection 2.4.4: Define Fit Result Parameters
  Double_t recoLinear3D_sigmaRES = recoLinear3D_sigmaRMS_fit->Value(2);
  Double_t recoLinear_gamma3D = recoLinear3D_sigmaRMS_fit->Value(3);

  Double_t recoPolygonal3D_sigmaRES = recoPolygonal3D_sigmaRMS_fit->Value(2);
  Double_t recoPolygonal_gamma3D = recoPolygonal3D_sigmaRMS_fit->Value(3);

  // TODO: Define the final values from 2.4.2 & 2.4.3




  /*
    OLD:
   */
  /*
  // True LinearXZ sigmaHL Fit.
  TF1* linearXZ_sigmaHL_function = new TF1("linearXZ_sigmaHL_function", SigmaHL, fitMin, fitMax, 2);
  linearXZ_sigmaHL_function->SetLineColor(kGreen);
  TFitResultPtr linearXZ_sigmaHL_fit = trueLinearXZ_sigmaHLVsSegmentMomentum->Fit(linearXZ_sigmaHL_function, "QS0");
  Double_t linearXZ_kappa_a = linearXZ_sigmaHL_fit->Value(0);
  Double_t linearXZ_kappa_c = linearXZ_sigmaHL_fit->Value(1);

  // True LinearYZ sigmaHL Fit.
  TF1* linearYZ_sigmaHL_function = new TF1("linearYZ_sigmaHL_function", SigmaHL, fitMin, fitMax, 2);
  linearYZ_sigmaHL_function->SetLineColor(kGreen);
  TFitResultPtr linearYZ_sigmaHL_fit = trueLinearYZ_sigmaHLVsSegmentMomentum->Fit(linearYZ_sigmaHL_function, "QS0");
  Double_t linearYZ_kappa_a = linearYZ_sigmaHL_fit->Value(0);
  Double_t linearYZ_kappa_c = linearYZ_sigmaHL_fit->Value(1);

  // True Polygonal sigmaHL Fit.
  TF1* polygonal_sigmaHL_function = new TF1("polygonal_sigmaHL_function", SigmaHL, fitMin, fitMax, 2);
  polygonal_sigmaHL_function->SetLineColor(kGreen);
  TFitResultPtr polygonal_sigmaHL_fit = truePolygonal_sigmaHLVsSegmentMomentum->Fit(polygonal_sigmaHL_function, "QS0");
  Double_t polygonal_kappa_a = polygonal_sigmaHL_fit->Value(0);
  Double_t polygonal_kappa_c = polygonal_sigmaHL_fit->Value(1);

  // True PolygonalXZ sigmaHL Fit.
  TF1* polygonalXZ_sigmaHL_function = new TF1("polygonalXZ_sigmaHL_function", SigmaHL, fitMin, fitMax, 2);
  polygonalXZ_sigmaHL_function->SetLineColor(kGreen);
  TFitResultPtr polygonalXZ_sigmaHL_fit = truePolygonalXZ_sigmaHLVsSegmentMomentum->Fit(polygonalXZ_sigmaHL_function, "QS0");
  Double_t polygonalXZ_kappa_a = polygonalXZ_sigmaHL_fit->Value(0);
  Double_t polygonalXZ_kappa_c = polygonalXZ_sigmaHL_fit->Value(1);

  // True PolygonalYZ sigmaHL Fit.
  TF1* polygonalYZ_sigmaHL_function = new TF1("polygonalYZ_sigmaHL_function", SigmaHL, fitMin, fitMax, 2);
  polygonalYZ_sigmaHL_function->SetLineColor(kGreen);
  TFitResultPtr polygonalYZ_sigmaHL_fit = truePolygonalYZ_sigmaHLVsSegmentMomentum->Fit(polygonalYZ_sigmaHL_function, "QS0");
  Double_t polygonalYZ_kappa_a = polygonalYZ_sigmaHL_fit->Value(0);
  Double_t polygonalYZ_kappa_c = polygonalYZ_sigmaHL_fit->Value(1);




  // TODO: Documentation for 3D fits. Include this stuff:
  //     Two fits per method
  //     Method one: Fix kappa_a & kappa_c with regular sigmaHL parameters, get the scaling factor (SUPPOSED to be sqrt(2)*sqrt(1-2/pi) but it's not quite that, I think due to the distributions being non-Gaussian)
  //     Method two: Fit kappa_a & kappa_c & scaling factor.
  //    What should we call 3D scaling factor? gamma_3D

  // True Linear3D sigmaHL Fit. Just fit gamma_3D, fix the other two.
  TF1* trueLinear3D_sigmaHL_function_1 = new TF1("trueLinear3D_sigmaHL_function_1", SigmaHL_3D, fitMin, fitMax, 3);
  trueLinear3D_sigmaHL_function_1->SetLineColor(kGreen);
  trueLinear3D_sigmaHL_function_1->FixParameter(1, linear_kappa_a);
  trueLinear3D_sigmaHL_function_1->FixParameter(2, linear_kappa_c);
  TFitResultPtr trueLinear3D_sigmaHL_fit_1 = trueLinear3D_sigmaVsSegmentMomentum->Fit(trueLinear3D_sigmaHL_function_1, "QS0");
  Double_t trueLinear_gamma3D_1 = trueLinear3D_sigmaHL_fit_1->Value(0);
  Double_t trueLinear_gamma3D_error_1 = trueLinear3D_sigmaHL_fit_1->Error(0);
  std::cout << "trueLinear_gamma3D_1 = " << trueLinear_gamma3D_1 << " +/- " << trueLinear_gamma3D_error_1 << std::endl;

  // True Linear3D sigmaHL Fit. Fit all three parameters.
  TF1* trueLinear3D_sigmaHL_function_2 = new TF1("trueLinear3D_sigmaHL_function_2", SigmaHL_3D, fitMin, fitMax, 3);
  trueLinear3D_sigmaHL_function_2->SetLineColor(kGreen);
  trueLinear3D_sigmaHL_function_2->SetParameters(0.8525, linear_kappa_a, linear_kappa_c); // initial parameter guesses
  TFitResultPtr trueLinear3D_sigmaHL_fit_2 = trueLinear3D_sigmaVsSegmentMomentum->Fit(trueLinear3D_sigmaHL_function_2, "QS0");
  
  // True Polygonal3D sigmaHL Fit.  Just fit gamma_3D, fix the other two.
  TF1* truePolygonal3D_sigmaHL_function_1 = new TF1("truePolygonal3D_sigmaHL_function_1", SigmaHL_3D, fitMin, fitMax, 3);
  truePolygonal3D_sigmaHL_function_1->SetLineColor(kGreen);
  truePolygonal3D_sigmaHL_function_1->FixParameter(1, polygonal_kappa_a);
  truePolygonal3D_sigmaHL_function_1->FixParameter(2, polygonal_kappa_c);
  TFitResultPtr truePolygonal3D_sigmaHL_fit_1 = truePolygonal3D_sigmaVsSegmentMomentum->Fit(truePolygonal3D_sigmaHL_function_1, "QS0");
  
  // True Polygonal3D sigmaHL Fit.  Fit all three parameters.
  TF1* truePolygonal3D_sigmaHL_function_2 = new TF1("truePolygonal3D_sigmaHL_function_2", SigmaHL_3D, fitMin, fitMax, 3);
  truePolygonal3D_sigmaHL_function_2->SetLineColor(kGreen);
  truePolygonal3D_sigmaHL_function_2->SetParameters(0.8525, polygonal_kappa_a, polygonal_kappa_c); // Initial Parameter Guesses
  TFitResultPtr truePolygonal3D_sigmaHL_fit_2 = truePolygonal3D_sigmaVsSegmentMomentum->Fit(truePolygonal3D_sigmaHL_function_2, "QS0");

  // TODO: Fit 3D sigma to different functions that include scaling factor in fit.  Or modifiy sigmaHL function to include that scaling factor and fix the parameter for 2D fits.
  // TODO: Fit Scaled 2D functions as well? Or Scale the other way too?
  // Linear theta3D sigmaHL fitting function and fit
  // TF1* linear3D_sigma_function = new TF1("linear3D_sigma_function", SigmaHL, fitMin, fitMax, 2);
  // linear3D_sigma_function->SetLineColor(kGreen);
  // TFitResultPtr linear3D_sigma_fit = trueLinear3D_sigmaVsSegmentMomentum->Fit(linear3D_sigma_function, "QS0");
  // 

  // Reco Linear sigmaRMS Fit. Fix parameters using linear sigmaHL fit.
  TF1* linear_sigmaRMS_function = new TF1("linear_sigmaRMS_function", SigmaRMS, fitMin, fitMax, 3);
  // Fix SigmaRMS parameters using SigmaHL fit parameters.
  linear_sigmaRMS_function->FixParameter(1, linear_kappa_a);
  linear_sigmaRMS_function->FixParameter(2, linear_kappa_c);
  TFitResultPtr linear_sigmaRMS_fit = recoLinear_sigmaRMSVsSegmentMomentum->Fit(linear_sigmaRMS_function, "QS0");
  Double_t linear_sigmaRES = linear_sigmaRMS_fit->Value(0);

  // Reco LinearXZ sigmaRMS Fit. Fix parameters using linear sigmaHL fit.
  TF1* linearXZ_sigmaRMS_function = new TF1("linearXZ_sigmaRMS_function", SigmaRMS, fitMin, fitMax, 3);
  linearXZ_sigmaRMS_function->FixParameter(1, linear_kappa_a);
  linearXZ_sigmaRMS_function->FixParameter(2, linear_kappa_c);
  TFitResultPtr linearXZ_sigmaRMS_fit = recoLinearXZ_sigmaRMSVsSegmentMomentum->Fit(linearXZ_sigmaRMS_function, "QS0");
  Double_t linearXZ_sigmaRES = linearXZ_sigmaRMS_fit->Value(0);

  // Reco LinearYZ sigmaRMS Fit. Fix parameters using linear sigmaHL fit.
  TF1* linearYZ_sigmaRMS_function = new TF1("linearYZ_sigmaRMS_function", SigmaRMS, fitMin, fitMax, 3);
  linearYZ_sigmaRMS_function->FixParameter(1, linear_kappa_a);
  linearYZ_sigmaRMS_function->FixParameter(2, linear_kappa_c);
  TFitResultPtr linearYZ_sigmaRMS_fit = recoLinearYZ_sigmaRMSVsSegmentMomentum->Fit(linearYZ_sigmaRMS_function, "QS0");
  Double_t linearYZ_sigmaRES = linearYZ_sigmaRMS_fit->Value(0);

  // Reco Polygonal sigmaRMS Fit. Fix Parameters using polygonal sigmaHL fit.
  TF1* polygonal_sigmaRMS_function = new TF1("polygonal_sigmaRMS_function", SigmaRMS, fitMin, fitMax, 3);
  polygonal_sigmaRMS_function->FixParameter(1, polygonal_kappa_a);
  polygonal_sigmaRMS_function->FixParameter(2, polygonal_kappa_c);
  TFitResultPtr polygonal_sigmaRMS_fit = recoPolygonal_sigmaRMSVsSegmentMomentum->Fit(polygonal_sigmaRMS_function, "QS0");
  Double_t polygonal_sigmaRES = polygonal_sigmaRMS_fit->Value(0);

  // Reco PolygonalXZ sigmaRMS Fit. Fix parameters using polygonal sigmaHL fit.
  TF1* polygonalXZ_sigmaRMS_function = new TF1("polyXZ_sigmaRMS_function", SigmaRMS, fitMin, fitMax, 3);
  // Fix a & c parameters used in the polygonal sigmaRMS fitting function
  polygonalXZ_sigmaRMS_function->FixParameter(1, polygonal_kappa_a);
  polygonalXZ_sigmaRMS_function->FixParameter(2, polygonal_kappa_c);
  TFitResultPtr polygonalXZ_sigmaRMS_fit = recoPolygonalXZ_sigmaRMSVsSegmentMomentum->Fit(polygonalXZ_sigmaRMS_function, "QS0");
  Double_t polygonalXZ_sigmaRES = polygonalXZ_sigmaRMS_fit->Value(0);

  // Reco PolygonalYZ sigmaRMS Fit. Fix parameters using polygonal sigmaHL fit.
  TF1* polygonalYZ_sigmaRMS_function = new TF1("polyYZ_sigmaRMS_function", SigmaRMS, fitMin, fitMax, 3);
  polygonalYZ_sigmaRMS_function->FixParameter(1, polygonal_kappa_a);
  polygonalYZ_sigmaRMS_function->FixParameter(2, polygonal_kappa_c);
  TFitResultPtr polygonalYZ_sigmaRMS_fit = recoPolygonalYZ_sigmaRMSVsSegmentMomentum->Fit(polygonalYZ_sigmaRMS_function, "QS0");
  Double_t polygonalYZ_sigmaRES = polygonalYZ_sigmaRMS_fit->Value(0);

  */

  // TODO: Reco Linear & reco Polygonal 3D fits.

  /////////////////////////////////////////////////////////////////////////////////////
  // Plot graphs in individual canvases. Save to files.
  /////////////////////////////////////////////////////////////////////////////////////

  // TODO: Save to canvases
  //recoPolygonalYZ_sigmaRMSVsSegmentMomentum->Draw("AP");
  //recoPolygonalXZ_sigmaRMSVsSegmentMomentum->Draw("P same");

  // TODO: How much of this can be placed into functions?

  // Range parameters
  Double_t yMax = 85.0;
  Double_t yMin = 0.0;

  // Plot linear theta3D sigma vs. segment momentum
  // Plot fit over that
  /*
  TCanvas* linear3D_sigmaVsSegmentMomentum_CANVAS = new TCanvas("", "", 1500, 1500);
  linear3D_sigmaVsSegmentMomentum_CANVAS->cd();
  trueLinear3D_sigmaVsSegmentMomentum->GetYaxis()->SetRangeUser(yMin, yMax);
  trueLinear3D_sigmaVsSegmentMomentum->Draw("AP"); // TODO: Change back to true
  TLegend* linear3D_sigmaVsSegmentMomentum_LEGEND = new TLegend(0.6, 0.75, 0.9, 0.9);
  linear3D_sigmaVsSegmentMomentum_LEGEND->SetHeader("Fit Parameters", "C");
  std::stringstream linear3D_kappa_a_ss;
  linear3D_kappa_a_ss << std::fixed << std::setprecision(3) << "#kappa_{a} = " << linear3D_sigma_fit->Value(0) << " MeV";
  std::stringstream linear3D_kappa_c_ss;
  linear3D_kappa_c_ss << std::fixed << std::setprecision(3) << "#kappa_{c} = " << linear3D_sigma_fit->Value(1) << " MeV";
  std::stringstream linear3D_sigma_fitChi2_ss;
  linear3D_sigma_fitChi2_ss << std::fixed << std::setprecision(3) << "#Chi^{2} = " << linear3D_sigma_fit->Chi2();
  linear3D_sigmaVsSegmentMomentum_LEGEND->AddEntry((TObject*)0, linear3D_kappa_a_ss.str().c_str(), "");
  linear3D_sigmaVsSegmentMomentum_LEGEND->AddEntry((TObject*)0, linear3D_kappa_c_ss.str().c_str(), "");
  linear3D_sigmaVsSegmentMomentum_LEGEND->AddEntry((TObject*)0, linear3D_sigma_fitChi2_ss.str().c_str(), "");
  linear3D_sigmaVsSegmentMomentum_LEGEND->Draw();
  linear3D_sigma_function->Draw("same");
  */

  // TODO: Documentation
  plotAndSave_anglesVsSegmentMomentum(trueLinear_theta3DVsSegmentMomentum_HIST, "trueLinear3D", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(trueLinear_thetaPrimeVsSegmentMomentum_HIST, "trueLinear", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(trueLinear_thetaXZprimeVsSegmentMomentum_HIST, "trueLinearXZ", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(trueLinear_thetaYZprimeVsSegmentMomentum_HIST, "trueLinearYZ", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(truePolygonal_theta3DVsSegmentMomentum_HIST,"truePolygonal3D", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(truePolygonal_thetaPrimeVsSegmentMomentum_HIST, "truePolygonal", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(truePolygonal_thetaXZprimeVsSegmentMomentum_HIST, "truePolygonalXZ", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(truePolygonal_thetaYZprimeVsSegmentMomentum_HIST, "truePolygonalYZ", "ADA_plots/");

  plotAndSave_anglesVsSegmentMomentum(recoLinear_theta3DVsSegmentMomentum_HIST, "recoLinear3D", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoLinear_thetaPrimeVsSegmentMomentum_HIST, "recoLinear", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoLinear_thetaXZprimeVsSegmentMomentum_HIST, "recoLinearXZ", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoLinear_thetaYZprimeVsSegmentMomentum_HIST, "recoLinearYZ", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoPolygonal_theta3DVsSegmentMomentum_HIST,"recoPolygonal3D", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoPolygonal_thetaPrimeVsSegmentMomentum_HIST, "recoPolygonal", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoPolygonal_thetaXZprimeVsSegmentMomentum_HIST, "recoPolygonalXZ", "ADA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoPolygonal_thetaYZprimeVsSegmentMomentum_HIST, "recoPolygonalYZ", "ADA_plots/");

  // TODO: Documentation
  //TLegend* linear_sigmaHL_LEGEND = GetSigmaHL_LEGEND(linear_sigmaHL_fit);
  // TLegend* linearXZ_sigmaHL_LEGEND = GetSigmaHL_LEGEND(linearXZ_sigmaHL_fit);
  // TLegend* linearYZ_sigmaHL_LEGEND = GetSigmaHL_LEGEND(linearYZ_sigmaHL_fit);
  // TLegend* polygonal_sigmaHL_LEGEND = GetSigmaHL_LEGEND(polygonal_sigmaHL_fit);
  // TLegend* polygonalXZ_sigmaHL_LEGEND = GetSigmaHL_LEGEND(polygonalXZ_sigmaHL_fit);
  // TLegend* polygonalYZ_sigmaHL_LEGEND = GetSigmaHL_LEGEND(polygonalYZ_sigmaHL_fit);

  // TLegend* trueLinear3D_sigmaHL_LEGEND_1 = GetSigmaHL_3D_LEGEND(trueLinear3D_sigmaHL_fit_1);
  // TLegend* trueLinear3D_sigmaHL_LEGEND_2 = GetSigmaHL_3D_LEGEND(trueLinear3D_sigmaHL_fit_2);

  // TLegend* truePolygonal3D_sigmaHL_LEGEND_1 = GetSigmaHL_3D_LEGEND(truePolygonal3D_sigmaHL_fit_1);
  // TLegend* truePolygonal3D_sigmaHL_LEGEND_2 = GetSigmaHL_3D_LEGEND(truePolygonal3D_sigmaHL_fit_2);

  // new legends
  TLegend* trueLinear_sigmaHL_LEGEND = GetSigma_LEGEND(trueLinear_sigmaHL_fit, "sg");
  TLegend* trueLinearXZ_sigmaHL_LEGEND = GetSigma_LEGEND(trueLinearXZ_sigmaHL_fit, "sg");
  TLegend* trueLinearYZ_sigmaHL_LEGEND = GetSigma_LEGEND(trueLinearYZ_sigmaHL_fit, "sg");

  TLegend* truePolygonal_sigmaHL_LEGEND = GetSigma_LEGEND(truePolygonal_sigmaHL_fit, "sg");
  TLegend* truePolygonalXZ_sigmaHL_LEGEND = GetSigma_LEGEND(truePolygonalXZ_sigmaHL_fit, "sg");
  TLegend* truePolygonalYZ_sigmaHL_LEGEND = GetSigma_LEGEND(truePolygonalYZ_sigmaHL_fit, "sg");

  TLegend* trueLinear3D_sigmaHL_LEGEND = GetSigma_LEGEND(trueLinear3D_sigmaHL_fit, "acs");
  TLegend* truePolygonal3D_sigmaHL_LEGEND = GetSigma_LEGEND(truePolygonal3D_sigmaHL_fit, "acs");

  TLegend* recoLinear_sigmaRMS_LEGEND = GetSigma_LEGEND(recoLinear_sigmaRMS_fit, "acg");
  TLegend* recoLinear_sigmaRMSraw_LEGEND = GetSigma_LEGEND(recoLinear_sigmaRMSraw_fit, "g");
  TLegend* recoLinearXZ_sigmaRMS_LEGEND = GetSigma_LEGEND(recoLinearXZ_sigmaRMS_fit, "acg");
  TLegend* recoLinearYZ_sigmaRMS_LEGEND = GetSigma_LEGEND(recoLinearYZ_sigmaRMS_fit, "acg");

  TLegend* recoPolygonal_sigmaRMS_LEGEND = GetSigma_LEGEND(recoPolygonal_sigmaRMS_fit, "acg");
  TLegend* recoPolygonalXZ_sigmaRMS_LEGEND = GetSigma_LEGEND(recoPolygonalXZ_sigmaRMS_fit, "acg");
  TLegend* recoPolygonalYZ_sigmaRMS_LEGEND = GetSigma_LEGEND(recoPolygonalYZ_sigmaRMS_fit, "acg");

  TLegend* recoLinear3D_sigmaRMS_LEGEND = GetSigma_LEGEND(recoLinear3D_sigmaRMS_fit, "ac");
  TLegend* recoPolygonal3D_sigmaRMS_LEGEND = GetSigma_LEGEND(recoPolygonal3D_sigmaRMS_fit, "ac");

  TLegend* recoLinear3D_sigmaRMS_fixSigmaRES_LEGEND = GetSigma_LEGEND(recoLinear3D_sigmaRMS_fit_fixSigmaRES, "acs");
  TLegend* recoPolygonal3D_sigmaRMS_fixSigmaRES_LEGEND = GetSigma_LEGEND(recoPolygonal3D_sigmaRMS_fit_fixSigmaRES, "acs");

  TLegend* recoLinear3D_sigmaRMS_fixGamma3D_LEGEND = GetSigma_LEGEND(recoLinear3D_sigmaRMS_fit_fixGamma3D, "acg");
  TLegend* recoPolygonal3D_sigmaRMS_fixGamma3D_LEGEND = GetSigma_LEGEND(recoPolygonal3D_sigmaRMS_fit_fixGamma3D, "acg");

  // TOOD: Documentation
  // plotAndSave_sigmaFitVsSegmentMomentum(trueLinear_sigmaHLVsSegmentMomentum, linear_sigmaHL_function, linear_sigmaHL_LEGEND, yMin, yMax, "trueLinear_sigmaHL");
  // plotAndSave_sigmaFitVsSegmentMomentum(trueLinearXZ_sigmaHLVsSegmentMomentum, linearXZ_sigmaHL_function, linearXZ_sigmaHL_LEGEND, yMin, yMax, "trueLinearXZ_sigmaHL");
  // plotAndSave_sigmaFitVsSegmentMomentum(trueLinearYZ_sigmaHLVsSegmentMomentum, linearYZ_sigmaHL_function, linearYZ_sigmaHL_LEGEND, yMin, yMax, "trueLinearYZ_sigmaHL");
  // plotAndSave_sigmaFitVsSegmentMomentum(truePolygonal_sigmaHLVsSegmentMomentum, polygonal_sigmaHL_function, polygonal_sigmaHL_LEGEND, yMin, yMax, "truePolygonal_sigmaHL");
  // plotAndSave_sigmaFitVsSegmentMomentum(truePolygonalXZ_sigmaHLVsSegmentMomentum, polygonalXZ_sigmaHL_function, polygonalXZ_sigmaHL_LEGEND, yMin, yMax, "truePolygonalXZ_sigmaHL");
  // plotAndSave_sigmaFitVsSegmentMomentum(truePolygonalYZ_sigmaHLVsSegmentMomentum, polygonalYZ_sigmaHL_function, polygonalYZ_sigmaHL_LEGEND, yMin, yMax, "truePolygonalYZ_sigmaHL");

  // plotAndSave_sigmaFitVsSegmentMomentum(trueLinear3D_sigmaVsSegmentMomentum, trueLinear3D_sigmaHL_function_1, trueLinear3D_sigmaHL_LEGEND_1, yMin, yMax, "trueLinear3D_1_sigmaHL");
  // plotAndSave_sigmaFitVsSegmentMomentum(trueLinear3D_sigmaVsSegmentMomentum, trueLinear3D_sigmaHL_function_2, trueLinear3D_sigmaHL_LEGEND_2, yMin, yMax, "trueLinear3D_2_sigmaHL");

  // plotAndSave_sigmaFitVsSegmentMomentum(truePolygonal3D_sigmaVsSegmentMomentum, truePolygonal3D_sigmaHL_function_1, truePolygonal3D_sigmaHL_LEGEND_1, yMin, yMax, "truePolygonal3D_1_sigmaHL");
  // plotAndSave_sigmaFitVsSegmentMomentum(truePolygonal3D_sigmaVsSegmentMomentum, truePolygonal3D_sigmaHL_function_2, truePolygonal3D_sigmaHL_LEGEND_2, yMin, yMax, "truePolygonal3D_2_sigmaHL");

  // True Linear 2D
  plotAndSave_sigmaFitVsSegmentMomentum(trueLinear_sigmaHLVsSegmentMomentum, trueLinear_sigmaHL_function, trueLinear_sigmaHL_LEGEND, yMin, yMax, "trueLinear_sigmaHL", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(trueLinearXZ_sigmaHLVsSegmentMomentum, trueLinearXZ_sigmaHL_function, trueLinearXZ_sigmaHL_LEGEND, yMin, yMax, "trueLinearXZ_sigmaHL", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(trueLinearYZ_sigmaHLVsSegmentMomentum, trueLinearYZ_sigmaHL_function, trueLinearYZ_sigmaHL_LEGEND, yMin, yMax, "trueLinearYZ_sigmaHL", "ADA_plots/");

  // True Polygonal 2D
  plotAndSave_sigmaFitVsSegmentMomentum(truePolygonal_sigmaHLVsSegmentMomentum, truePolygonal_sigmaHL_function, truePolygonal_sigmaHL_LEGEND, yMin, yMax, "truePolygonal_sigmaHL", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(truePolygonalXZ_sigmaHLVsSegmentMomentum, truePolygonalXZ_sigmaHL_function, truePolygonalXZ_sigmaHL_LEGEND, yMin, yMax, "truePolygonalXZ_sigmaHL", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(truePolygonalYZ_sigmaHLVsSegmentMomentum, truePolygonalYZ_sigmaHL_function, truePolygonalYZ_sigmaHL_LEGEND, yMin, yMax, "truePolygonalYZ_sigmaHL", "ADA_plots/");

  // True Linear/Polygonal 3D
  plotAndSave_sigmaFitVsSegmentMomentum(trueLinear3D_sigmaVsSegmentMomentum, trueLinear3D_sigmaHL_function, trueLinear3D_sigmaHL_LEGEND, yMin, yMax, "trueLinear3D_sigmaHL", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(truePolygonal3D_sigmaVsSegmentMomentum, truePolygonal3D_sigmaHL_function, truePolygonal3D_sigmaHL_LEGEND, yMin, yMax, "truePolygonal3D_sigmaHL", "ADA_plots/");

  // Reco Linear 2D
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear_sigmaRMSVsSegmentMomentum, recoLinear_sigmaRMS_function, recoLinear_sigmaRMS_LEGEND, yMin, yMax, "recoLinear_sigmaRMS", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear_sigmaRMSVsSegmentMomentum, recoLinear_sigmaRMSraw_function, recoLinear_sigmaRMSraw_LEGEND, yMin, yMax, "recoLinear_sigmaRMSraw", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinearXZ_sigmaRMSVsSegmentMomentum, recoLinearXZ_sigmaRMS_function, recoLinearXZ_sigmaRMS_LEGEND, yMin, yMax, "recoLinearXZ_sigmaRMS", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinearYZ_sigmaRMSVsSegmentMomentum, recoLinearYZ_sigmaRMS_function, recoLinearYZ_sigmaRMS_LEGEND, yMin, yMax, "recoLinearYZ_sigmaRMS", "ADA_plots/");

  // Reco Polygonal 2D
  plotAndSave_sigmaFitVsSegmentMomentum(recoPolygonal_sigmaRMSVsSegmentMomentum, recoPolygonal_sigmaRMS_function, recoPolygonal_sigmaRMS_LEGEND, yMin, yMax, "recoPolygonal_sigmaRMS", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(recoPolygonalXZ_sigmaRMSVsSegmentMomentum, recoPolygonalXZ_sigmaRMS_function, recoPolygonalXZ_sigmaRMS_LEGEND, yMin, yMax, "recoPolygonalXZ_sigmaRMS", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(recoPolygonalYZ_sigmaRMSVsSegmentMomentum, recoPolygonalYZ_sigmaRMS_function, recoPolygonalYZ_sigmaRMS_LEGEND, yMin, yMax, "recoPolygonalYZ_sigmaRMS", "ADA_plots/");

  // Reco Linear/Polygonal 3D
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear3D_sigmaVsSegmentMomentum, recoLinear3D_sigmaRMS_function, recoLinear3D_sigmaRMS_LEGEND, yMin, yMax, "recoLinear3D_sigmaRMS", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(recoPolygonal3D_sigmaVsSegmentMomentum, recoPolygonal3D_sigmaRMS_function, recoPolygonal3D_sigmaRMS_LEGEND, yMin, yMax, "recoPolygonal3D_sigmaRMS", "ADA_plots/");

  // Reco Linear/Polygonal 3D - Fix sigmaRES
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear3D_sigmaVsSegmentMomentum, recoLinear3D_sigmaRMS_function_fixSigmaRES, recoLinear3D_sigmaRMS_fixSigmaRES_LEGEND, yMin, yMax, "recoLinear3D_fixSigmaRES_sigmaRMS", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(recoPolygonal3D_sigmaVsSegmentMomentum, recoPolygonal3D_sigmaRMS_function_fixSigmaRES, recoPolygonal3D_sigmaRMS_fixSigmaRES_LEGEND, yMin, yMax, "recoPolygonal3D_fixSigmaRES_sigmaRMS", "ADA_plots/");

  // Reco Linear/Polygonal 3D - Fix gamma3D
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear3D_sigmaVsSegmentMomentum, recoLinear3D_sigmaRMS_function_fixGamma3D, recoLinear3D_sigmaRMS_fixGamma3D_LEGEND, yMin, yMax, "recoLinear3D_fixGamma3D_sigmaRMS", "ADA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(recoPolygonal3D_sigmaVsSegmentMomentum, recoPolygonal3D_sigmaRMS_function_fixGamma3D, recoPolygonal3D_sigmaRMS_fixGamma3D_LEGEND, yMin, yMax, "recoPolygonal3D_fixGamma3D_sigmaRMS", "ADA_plots/");

  // Plot sigmaHL & sigmaRMS overlayed over each other.
  plotAndSave_sigmaHL_sigmaRMS_overlay(trueLinear_sigmaHLVsSegmentMomentum, trueLinear_sigmaHL_function, recoLinear_sigmaRMSVsSegmentMomentum, recoLinear_sigmaRMS_function, yMin, yMax, "linear", "ADA_plots/");
  plotAndSave_sigmaHL_sigmaRMS_overlay(trueLinearXZ_sigmaHLVsSegmentMomentum, trueLinearXZ_sigmaHL_function, recoLinearXZ_sigmaRMSVsSegmentMomentum, recoLinearXZ_sigmaRMS_function, yMin, yMax, "linearXZ", "ADA_plots/");
  plotAndSave_sigmaHL_sigmaRMS_overlay(trueLinearYZ_sigmaHLVsSegmentMomentum, trueLinearYZ_sigmaHL_function, recoLinearYZ_sigmaRMSVsSegmentMomentum, recoLinearYZ_sigmaRMS_function, yMin, yMax, "linearYZ", "ADA_plots/");

  // // TODO: Documentation
  // TLegend* linear_sigmaRMS_LEGEND = GetSigmaRMS_LEGEND(linear_sigmaRMS_fit);
  // TLegend* linearXZ_sigmaRMS_LEGEND = GetSigmaRMS_LEGEND(linearXZ_sigmaRMS_fit);
  // TLegend* linearYZ_sigmaRMS_LEGEND = GetSigmaRMS_LEGEND(linearYZ_sigmaRMS_fit);
  // TLegend* polygonal_sigmaRMS_LEGEND = GetSigmaRMS_LEGEND(polygonal_sigmaRMS_fit);
  // TLegend* polygonalXZ_sigmaRMS_LEGEND = GetSigmaRMS_LEGEND(polygonalXZ_sigmaRMS_fit);
  // TLegend* polygonalYZ_sigmaRMS_LEGEND = GetSigmaRMS_LEGEND(polygonalYZ_sigmaRMS_fit);

  // // TODO: Documentation
  // plotAndSave_sigmaFitVsSegmentMomentum(recoLinear_sigmaRMSVsSegmentMomentum, linear_sigmaRMS_function, linear_sigmaRMS_LEGEND, yMin, yMax, "recoLinear_sigmaRMS");
  // plotAndSave_sigmaFitVsSegmentMomentum(recoLinearXZ_sigmaRMSVsSegmentMomentum, linearXZ_sigmaRMS_function, linearXZ_sigmaRMS_LEGEND, yMin, yMax, "recoLinearXZ_sigmaRMS");
  // plotAndSave_sigmaFitVsSegmentMomentum(recoLinearYZ_sigmaRMSVsSegmentMomentum, linearYZ_sigmaRMS_function, linearYZ_sigmaRMS_LEGEND, yMin, yMax, "recoLinearYZ_sigmaRMS");
  // plotAndSave_sigmaFitVsSegmentMomentum(recoPolygonal_sigmaRMSVsSegmentMomentum, polygonal_sigmaRMS_function, polygonal_sigmaRMS_LEGEND, yMin, yMax, "recoPolygonal_sigmaRMS");
  // plotAndSave_sigmaFitVsSegmentMomentum(recoPolygonalXZ_sigmaRMSVsSegmentMomentum, polygonalXZ_sigmaRMS_function, polygonalXZ_sigmaRMS_LEGEND, yMin, yMax, "recoPolygonalXZ_sigmaRMS");
  // plotAndSave_sigmaFitVsSegmentMomentum(recoPolygonalYZ_sigmaRMSVsSegmentMomentum, polygonalYZ_sigmaRMS_function, polygonalYZ_sigmaRMS_LEGEND, yMin, yMax, "recoPolygonalYZ_sigmaRMS");

  /////////////////////////////////////////////////////////////////////////////////////
  // Lots of TODO items
  /////////////////////////////////////////////////////////////////////////////////////

  // TODO: Get sigmaRES discretely

  // TODO: Calculate sigma via Gaussian fit rather than using Stdev & StdevError TH1F methods?
  //     I can just use TH1F's Fit method w/ a Gaus function, should be quite easy.
  //     With this method, I can plot the Chi2 of each of the fits as well.

  // TODO: Plot raw sigma over sigmaHL/sigmaRMS (shouldn't matter, but do for both), vs. Segment Momentum
  //    sigmaRatioVsSegmentMomentum? Should be constant at Half-Gaussian vs. sqrt(2) calculation
  //    For half-normal: stdev measured = sigma_gen * sqrt(1-2/pi)
  //    For projection, stdev measured = sigma_gen / sqrt(2)
  //    Ratio of rawSigma / projectedSigma = sqrt(1-2/pi)*sqrt(2)

  // TODO: Draw graphs, save to TTree?  Then we could save the individual TH1F calculated in GetSigmaGraph function.
  // Directly compare RANSAC and standard results
  // Possible TODO at a later date: We could easily make a TGraph2DErrors that incorporates segmentLength, then we could easily do a multi dimensional fit. (provide TF2 instead of TF1)
}

// TODO: Documentation
TLegend* GetSigmaRMS_LEGEND(TFitResultPtr sigmaRMS_fit) {

  TLegend* legend = new TLegend(0.4, 0.70, 0.9, 0.9);
  legend->SetTextSize(0.035);

  std::stringstream sigmaRES_ss;
  sigmaRES_ss << std::fixed << std::setprecision(3) << "#sigma_{RES} = " << sigmaRMS_fit->Value(0) << " mrad";
  std::stringstream chi2_ss;
  chi2_ss << std::fixed << std::setprecision(3) << "#Chi^{2} = " << sigmaRMS_fit->Chi2();

  //legend->AddEntry(linear_sigmaHL_function, "Linear #sigma_{HL} fit"); // TODO: Consider adding sigmaHL fitted function for sigmaRMS overlay.
  //legend->AddEntry(sigmaRMS_function, "Linear #sigma_{RMS} fit, with:"); // TODO: If sigmaHL fit is added, then we need to have legend entries for both functions.
  legend->AddEntry((TObject*)0, sigmaRES_ss.str().c_str(), "");
  legend->AddEntry((TObject*)0, chi2_ss.str().c_str(), "");

  return legend;
}

// TOOD: Documentation
// TODO: This function could be wrapped into plotAndSave_sigmaFitVsSegmentMomentum if the sigmaHL & sigmaRMS fit functions were the same.  Not sure how this would also work with 3D fits.  Perhaps all in one and it just fix a bunch of parameters by extracting functions.
TLegend* GetSigmaHL_LEGEND(TFitResultPtr sigmaHL_fit) {

  TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.9);
  legend->SetHeader("Fit Parameters", "C");

  std::stringstream kappa_a_ss;
  kappa_a_ss << std::fixed << std::setprecision(3) << "#kappa_{a} = " << sigmaHL_fit->Value(0) << " MeV";
  std::stringstream kappa_c_ss;
  kappa_c_ss << std::fixed << std::setprecision(3) << "#kappa_{c} = " << sigmaHL_fit->Value(1) << " MeV";
  std::stringstream chi2_ss;
  chi2_ss << std::fixed << std::setprecision(3) << "#Chi^{2} = " << sigmaHL_fit->Chi2();

  legend->AddEntry((TObject*)0, kappa_a_ss.str().c_str(), "");
  legend->AddEntry((TObject*)0, kappa_c_ss.str().c_str(), "");
  legend->AddEntry((TObject*)0, chi2_ss.str().c_str(), "");

  return legend;
}

// TOOD: Documentation
TLegend* GetSigmaHL_3D_LEGEND(TFitResultPtr sigmaHL_3D_fit) {
  TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.9);
  legend->SetHeader("Fit Parameters", "C");

  std::stringstream gamma_3D_ss;
  gamma_3D_ss << std::fixed << std::setprecision(3) << "#gamma_{3D} = " << sigmaHL_3D_fit->Value(0);
  std::stringstream kappa_a_ss;
  kappa_a_ss << std::fixed << std::setprecision(3) << "#kappa_{a} = " << sigmaHL_3D_fit->Value(1) << " MeV";
  std::stringstream kappa_c_ss;
  kappa_c_ss << std::fixed << std::setprecision(3) << "#kappa_{c} = " << sigmaHL_3D_fit->Value(2) << " MeV";
  std::stringstream chi2_ss;
  chi2_ss << std::fixed << std::setprecision(3) << "Chi^{2} = " << sigmaHL_3D_fit->Chi2();

  legend->AddEntry((TObject*)0, gamma_3D_ss.str().c_str(), "");
  legend->AddEntry((TObject*)0, kappa_a_ss.str().c_str(), "");
  legend->AddEntry((TObject*)0, kappa_c_ss.str().c_str(), "");
  legend->AddEntry((TObject*)0, chi2_ss.str().c_str(), "");

  return legend;
}



// // Specifically to test fitting
// // Get a TGraphErrors of sigma vs. x for the provided histogram, where the y-errors are the standard deviation of each x-bin of the histogram.
// TGraphErrors* GetRawSigmaGraph(TH2F* hist2D, const char* name, const char* title) {
//   // Define graph
//   TGraphErrors* graph = new TGraphErrors();
//   graph->SetNameTitle(name, title);

//   // Define x and y axes
//   TAxis* xAxis = hist2D->GetXaxis();
//   TAxis* yAxis = hist2D->GetYaxis();

//   // Loop through x-bins, to fill graph
//   Int_t nXbins = hist2D->GetNbinsX();
//   for(Int_t xBin = 0; xBin < nXbins; ++xBin) {    // Define x for this bin
//     Double_t x = xAxis->GetBinCenter(xBin);

//     // Define y-parameters for this x-bin (same for all x-bins)
//     Int_t nYbins = hist2D->GetNbinsY();
//     Double_t yMax = yAxis->GetXmax();
//     Double_t yMin = yAxis->GetXmin();
//     // Define 1D histogram for this x-bin
//     TH1F* hist1D = new TH1F("", "", nYbins, yMin, yMax);
//     //hist1D->GetXaxis()->Set(nYbins, yMin, yMax);
//     // hist1D->GetXaxis()->SetRangeUser(yMin, yMax);
//     //TH1F* hist = new TH1F("hist name", "hist title", nYbins, yMax, yMin);
    
//     // Loop through ybins for this x-bin to fill the 1D histogram
//     for(Int_t yBin = 0; yBin < nYbins; ++yBin) {
//       Double_t y = yAxis->GetBinCenter(yBin);
//       Double_t weight = hist2D->GetBinContent(xBin, yBin);

//       // Fill histogram with values
//       hist1D->Fill(y, weight);
//     } // Loop through y-bins
    
//     // TODO: Add an option to get by fit.
//     Double_t sigma = hist1D->GetStdDev();
//     Double_t sigmaError = hist1D->GetStdDevError();

//     // if(sigma != 0) {
//     //   //std::cout << std::fixed << std::setprecision(3) << "sigma = " << sigma << "\t" << "sigmaMLE = " << sigmaMLE << std::endl;
//     //   std::cout << "fit half-gauss for bin: " << xBin << std::endl;
//     //   TFitResultPtr fit = hist1D->Fit("gaus","QS","", 0, 250);
//     //   std::cout << "chi2 = " << fit->Chi2() << std::endl;

//     //   std::cout << "by half-gauss fit: " << fit->Value(2) << std::endl;
//     //   std::cout << "by GetStdDev: " << sigma << std::endl;
//     // }

//     // Fill Graph
//     // TODO: Change to a number of points w/in the histogram cutoff, maybe like 50?
//     if(sigma != 0 /* && (x <= 2.0 && x >= 0.2) */) {

//       TFitResultPtr fit = hist1D->Fit("gaus", "QS0", "", 0, 250);
//       sigma = fit->Value(2);
//       sigmaError = fit->Error(2);

//       TF1* gaus = new TF1("myGaus", Gaus, 0, 250, 3);
//       TF1* halfGaus = new TF1("myHalfGaus", HalfGaus, 0, 250, 2);

//       TFitResultPtr gaus_fit = hist1D->Fit(gaus, "QS0");
//       TFitResultPtr halfGaus_fit = hist1D->Fit(halfGaus, "QS0");

//       std::cout << "sigma = " << sigma << "\t"
// 		<< "gaus_fit = " << gaus_fit->Value(2) << "\t"
// 		<< "halfGaus_fit = " << halfGaus_fit->Value(1) << "\t"
// 		<< std::endl;

//       Int_t pointNum = graph->GetN();
//       graph->SetPoint(pointNum, x, sigma);
//       graph->SetPointError(pointNum, 0, sigmaError);
//     }
//   }

//   // Stylize graph
//   graph->SetMarkerColor(kRed);
//   graph->SetMarkerStyle(kFullCross);
//   graph->SetLineColor(kBlue);

//   // Return graph
//   return graph;
// }

// TODO: Documentation on functions, explain parameters better.

// // Function used to fit sigmaRMS vs. Segment Momentum
// // paras[0] = sigmaRES, paras[1] = a, paras[2] = c (modified highland formula fit parameters)
// // In normal fit, paras 1 & 2 are held constant, provided by the fit from SigmaHL
// // In reco raw fit, paras 0 is held constant
// Double_t SigmaRMS(Double_t* x, Double_t* paras) {
//   Double_t* kappaParameters = (Double_t[]) { paras[1], paras[2] };
//   Double_t sigmaHL = SigmaHL(x, kappaParameters);

//   // Add in quadrature and return
//   return std::pow(sigmaHL*sigmaHL + paras[0]*paras[0],0.5);
// }

// // Function used to fit sigmaHL vs. Segment Momentum
// Double_t SigmaHL(Double_t* x, Double_t* paras) {
//   Double_t momentum = x[0]; // GeV/c
//   Double_t momentum2 = momentum*momentum; // Momentum Squared
//   Double_t kappa = paras[0] / momentum2 + paras[1];
//   Double_t pbetac = momentum2 / pow(0.106*0.106+momentum2,0.5); // p*beta*c

//   return kappa/pbetac; // mrad
// }

// // Function used to fit sigmaHL_3D vs. Segment Momentum.  kappa_a & kappa_c may or may not be fixed during the fit.
// Double_t SigmaHL_3D(Double_t* x, Double_t* paras) {

//   Double_t gamma_3D = paras[0];
//   Double_t kappa_a = paras[1];
//   Double_t kappa_c = paras[2];

//   Double_t* kappaParameters = (Double_t[]) { kappa_a, kappa_c };
//   Double_t sigmaHL = SigmaHL(x, kappaParameters);

//   return gamma_3D*sigmaHL;
// }

// Double_t Gaus(Double_t* x, Double_t* paras) {
//   Double_t x0 = x[0];

//   Double_t amp = paras[0];
//   Double_t mu = paras[1];
//   Double_t sigma = paras[2];

//   return amp*exp(-0.5*pow((x0-mu)/sigma,2));
// }

// Double_t HalfGaus(Double_t* x, Double_t* paras) {
//   Double_t x0 = x[0];
  
//   Double_t amp = paras[0];
//   Double_t sigma = paras[1];

//   return amp*exp(-0.5*pow((x0/sigma),2));
// }
