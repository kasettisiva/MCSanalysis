#include "HelperFunctions.cpp"
#include "PlottingFunctions.cpp"
#include "AngleFittingFunctions.cpp"

// Root MACRO Function
void firstSegmentAnalysis() {

  // TODO: Set this to kTRUE
  gROOT->SetBatch(kFALSE);

  /////////////////////////////////////////////////////////////////////////////////////
  // Get plots from root file
  /////////////////////////////////////////////////////////////////////////////////////

  // Reset the output directory.
  gSystem->Exec("rm -r FSA_plots");
  gSystem->Exec("mkdir FSA_plots");

  // Get Root File
  const char* path = "../MCSAngleAnalysis_hist.root";
  TFile* file = new TFile(path);

  // TODO: Add to GitHub Issues:
  // TODO: Functionalize the flow to allow for other segmentation methods easily.
  // TODO: Can we put individual plots that are made into individual function calls?
  // TODO: For functionalized fits, can we add a Properties struct that contains ALL of the possible changes that can be made into the fits? This allows for us to make only ONE function that just takes a BUNCH of parameters wrapped into one struct.
  // TODO: Are 3D angles being fit to all possible scenarios?
  // TODO: Are reco angles being fit to all possible scenarios?
  // TODO: Polygonal Fits
  // TODO: Plot the fit-statistics of each TH1F fit. Chi2, "Gaussian"-ness, etc.
  // TODO: Fit the 3D angles using a half-gaussian metric.
  // TODO: Combine firstSegmentAnalysis with angleDistributionAnalysis.

  // Get each histogram from the root file

  // True 'Linear' Plots
  TH2F* trueLinear_theta3DVsSegmentMomentum_HIST = Get2DHist(file, "first_trueLinear_theta3DVsSegmentMomentum_HIST");
  TH2F* trueLinear_thetaXZprimeVsSegmentMomentum_HIST = Get2DHist(file, "first_trueLinear_thetaXZprimeVsSegmentMomentum_HIST");
  TH2F* trueLinear_thetaYZprimeVsSegmentMomentum_HIST = Get2DHist(file, "first_trueLinear_thetaYZprimeVsSegmentMomentum_HIST");

  // Reco 'Linear' Plots
  TH2F* recoLinear_theta3DVsSegmentMomentum_HIST = Get2DHist(file, "first_recoLinear_theta3DVsSegmentMomentum_HIST");
  TH2F* recoLinear_thetaXZprimeVsSegmentMomentum_HIST = Get2DHist(file, "first_recoLinear_thetaXZprimeVsSegmentMomentum_HIST");
  TH2F* recoLinear_thetaYZprimeVsSegmentMomentum_HIST = Get2DHist(file, "first_recoLinear_thetaYZprimeVsSegmentMomentum_HIST");

  // True 'Linear' Plots (BB Momentum)
  TH2F* trueLinear_theta3DVsSegmentBBMomentum_HIST = Get2DHist(file, "first_trueLinear_theta3DVsSegmentBBMomentum_HIST");
  TH2F* trueLinear_thetaXZprimeVsSegmentBBMomentum_HIST = Get2DHist(file, "first_trueLinear_thetaXZprimeVsSegmentBBMomentum_HIST");
  TH2F* trueLinear_thetaYZprimeVsSegmentBBMomentum_HIST = Get2DHist(file, "first_trueLinear_thetaYZprimeVsSegmentBBMomentum_HIST");

  // Reco 'Linear' Plots (BB Momentum)
  TH2F* recoLinear_theta3DVsSegmentBBMomentum_HIST = Get2DHist(file, "first_recoLinear_theta3DVsSegmentBBMomentum_HIST");
  TH2F* recoLinear_thetaXZprimeVsSegmentBBMomentum_HIST = Get2DHist(file, "first_recoLinear_thetaXZprimeVsSegmentBBMomentum_HIST");
  TH2F* recoLinear_thetaYZprimeVsSegmentBBMomentum_HIST = Get2DHist(file, "first_recoLinear_thetaYZprimeVsSegmentBBMomentum_HIST");
  
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

  recoLinear_theta3DVsSegmentMomentum_HIST->RebinX(rebinFactor);
  recoLinear_thetaXZprimeVsSegmentMomentum_HIST->RebinX(rebinFactor);
  recoLinear_thetaYZprimeVsSegmentMomentum_HIST->RebinX(rebinFactor);

  trueLinear_theta3DVsSegmentBBMomentum_HIST->RebinX(rebinFactor);
  trueLinear_thetaXZprimeVsSegmentBBMomentum_HIST->RebinX(rebinFactor);
  trueLinear_thetaYZprimeVsSegmentBBMomentum_HIST->RebinX(rebinFactor);

  recoLinear_theta3DVsSegmentBBMomentum_HIST->RebinX(rebinFactor);
  recoLinear_thetaXZprimeVsSegmentBBMomentum_HIST->RebinX(rebinFactor);
  recoLinear_thetaYZprimeVsSegmentBBMomentum_HIST->RebinX(rebinFactor);
  
  //
  // Subsection 1.2: Combine XZ & YZ histograms
  //

  // True Momentum (closestZ for reco)
  TH2F* trueLinear_thetaPrimeVsSegmentMomentum_HIST = CombineHistograms(trueLinear_thetaXZprimeVsSegmentMomentum_HIST, trueLinear_thetaYZprimeVsSegmentMomentum_HIST, "trueLinear_thetaPrimeVsSegmentMomentum_HIST", "True Linear #theta' Vs. Segment Momentum; Segment Momentum (GeV/c); #theta' (mrad)");

  TH2F* recoLinear_thetaPrimeVsSegmentMomentum_HIST = CombineHistograms(recoLinear_thetaXZprimeVsSegmentMomentum_HIST, recoLinear_thetaYZprimeVsSegmentMomentum_HIST, "recoLinear_thetaPrimeVsSegmentMomentum_HIST", "Reco Linear #theta' Vs. Segment Momentum; Segment Momentum (GeV/c); #theta' (mrad)");

  // BB Momentum
  TH2F* trueLinear_thetaPrimeVsSegmentBBMomentum_HIST = CombineHistograms(trueLinear_thetaXZprimeVsSegmentBBMomentum_HIST, trueLinear_thetaYZprimeVsSegmentBBMomentum_HIST, "trueLinear_thetaPrimeVsSegmentBBMomentum_HIST", "True Linear #theta' Vs. Segment BB Momentum; Segment BB Momentum (GeV/c); #theta' (mrad)");

  TH2F* recoLinear_thetaPrimeVsSegmentBBMomentum_HIST = CombineHistograms(recoLinear_thetaXZprimeVsSegmentBBMomentum_HIST, recoLinear_thetaYZprimeVsSegmentBBMomentum_HIST, "recoLinear_thetaPrimeVsSegmentBBMomentum_HIST", "Reco Linear #theta' Vs. Segment BB Momentum; Segment BB Momentum (GeV/c); #theta' (mrad)");

  //
  // Subsection 1.3: Get Graph of sigmaHL/sigmaRMS vs. segment momentum for 2D angles. (with y-uncertainties)
  //

  // TODO: Finish incorporating lastPoint plots.

  // Subsubsection 1.3.1: Ideal future, if there is no direction discrepancy. (for reco, this should already be the case for true)
  // True Momentum (closestZ for reco)
  TGraphErrors* trueLinear_sigmaHLVsSegmentMomentum = GetSigmaGraph(trueLinear_thetaPrimeVsSegmentMomentum_HIST, "trueLinear_sigmaHLVsSegmentMomentum", "True Linear #sigma_{HL} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{HL} (mrad)");

  TGraphErrors* recoLinear_sigmaRMSVsSegmentMomentum = GetSigmaGraph(recoLinear_thetaPrimeVsSegmentMomentum_HIST, "recoLinear_sigmaRMSVsSegmentMomentum", "Reco Linear #sigma_{RMS} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{RMS} (mrad)");

  // BB Momentum
  TGraphErrors* trueLinear_sigmaHLVsSegmentBBMomentum = GetSigmaGraph(trueLinear_thetaPrimeVsSegmentBBMomentum_HIST, "trueLinear_sigmaHLVsSegmentBBMomentum", "True Linear #sigma_{HL} Vs. Segment BB Momentum; Segment BB Momentum (GeV/c); #sigma_{HL} (mrad)");

  TGraphErrors* recoLinear_sigmaRMSVsSegmentBBMomentum = GetSigmaGraph(recoLinear_thetaPrimeVsSegmentBBMomentum_HIST, "recoLinear_sigmaRMSVsSegmentBBMomentum", "Reco Linear #sigma_{RMS} Vs. Segment BB Momentum; Segment BB Momentum (GeV/c); #sigma_{RMS} (mrad)");

  // Subsubsection 1.3.2: Split into two projection directions, due to direction discrepancy.  Should only need to do reco, but we're also doing true to cover all bases.
  // True Momentum (closestZ for reco)
  TGraphErrors* trueLinearXZ_sigmaHLVsSegmentMomentum = GetSigmaGraph(trueLinear_thetaXZprimeVsSegmentMomentum_HIST, "trueLinearXZ_sigmaHLVsSegmentMomentum", "True Linear XZ #sigma_{HL} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{HL} (mrad)");
  TGraphErrors* trueLinearYZ_sigmaHLVsSegmentMomentum = GetSigmaGraph(trueLinear_thetaYZprimeVsSegmentMomentum_HIST, "trueLinearYZ_sigmaHLVsSegmentMomentum", "True Linear YZ #sigma_{HL} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{HL} (mrad)");

  TGraphErrors* recoLinearXZ_sigmaRMSVsSegmentMomentum = GetSigmaGraph(recoLinear_thetaXZprimeVsSegmentMomentum_HIST, "recoLinearXZ_sigmaRMSVsSegmentMomentum", "Reco Linear XZ #sigma_{RMS} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{RMS} (mrad)");
  TGraphErrors* recoLinearYZ_sigmaRMSVsSegmentMomentum = GetSigmaGraph(recoLinear_thetaYZprimeVsSegmentMomentum_HIST, "recoLinearXZ_sigmaRMSVsSegmentMomentum", "Reco Linear YZ #sigma_{RMS} Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma_{RMS} (mrad)");

  // BB Momentum
  TGraphErrors* trueLinearXZ_sigmaHLVsSegmentBBMomentum = GetSigmaGraph(trueLinear_thetaXZprimeVsSegmentBBMomentum_HIST, "trueLinearXZ_sigmaHLVsSegmentBBMomentum", "True Linear XZ #sigma_{HL} Vs. Segment BB Momentum; Segment BB Momentum (GeV/c); #sigma_{HL} (mrad)");
  TGraphErrors* trueLinearYZ_sigmaHLVsSegmentBBMomentum = GetSigmaGraph(trueLinear_thetaYZprimeVsSegmentBBMomentum_HIST, "trueLinearYZ_sigmaHLVsSegmentBBMomentum", "True Linear YZ #sigma_{HL} Vs. Segment BB Momentum; Segment BB Momentum (GeV/c); #sigma_{HL} (mrad)");

  TGraphErrors* recoLinearXZ_sigmaRMSVsSegmentBBMomentum = GetSigmaGraph(recoLinear_thetaXZprimeVsSegmentBBMomentum_HIST, "recoLinearXZ_sigmaRMSVsSegmentBBMomentum", "Reco Linear XZ #sigma_{RMS} Vs. Segment BB Momentum; Segment BB Momentum (GeV/c); #sigma_{RMS} (mrad)");
  TGraphErrors* recoLinearYZ_sigmaRMSVsSegmentBBMomentum = GetSigmaGraph(recoLinear_thetaYZprimeVsSegmentBBMomentum_HIST, "recoLinearXZ_sigmaRMSVsSegmentBBMomentum", "Reco Linear YZ #sigma_{RMS} Vs. Segment BB Momentum; Segment BB Momentum (GeV/c); #sigma_{RMS} (mrad)");

  // Subsubsection 1.3.3: Get Graph of sigma vs. segment momentum for 3D angles. (with y-uncertainties)
  // True Momentum (closestZ for reco)
  TGraphErrors* trueLinear3D_sigmaVsSegmentMomentum = GetSigmaGraph(trueLinear_theta3DVsSegmentMomentum_HIST, "trueLinear3D_sigmaVsSegmentMomentum", "True Linear 3D #sigma Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma (mrad)");
  TGraphErrors* recoLinear3D_sigmaVsSegmentMomentum = GetSigmaGraph(recoLinear_theta3DVsSegmentMomentum_HIST, "recoLinear3D_sigmaVsSegmentMomentum", "Reco Linear 3D #sigma Vs. Segment Momentum; Segment Momentum (GeV/c); #sigma (mrad)");

  // BB Momentum
  TGraphErrors* trueLinear3D_sigmaVsSegmentBBMomentum = GetSigmaGraph(trueLinear_theta3DVsSegmentBBMomentum_HIST, "trueLinear3D_sigmaVsSegmentBBMomentum", "True Linear 3D #sigma Vs. Segment BB Momentum; Segment BB Momentum (GeV/c); #sigma (mrad)");
  TGraphErrors* recoLinear3D_sigmaVsSegmentBBMomentum = GetSigmaGraph(recoLinear_theta3DVsSegmentBBMomentum_HIST, "recoLinear3D_sigmaVsSegmentBBMomentum", "Reco Linear 3D #sigma Vs. Segment BB Momentum; Segment BB Momentum (GeV/c); #sigma (mrad)");

  //
  // Subsection 1.4: Get Ratio of 3D and 2D angle distributions.
  // TODO: Add graph titles, use prefixes like I did in momentumAnalysis.cpp
  //

  // Subsubsection 1.4.1: Ideal future, if there is no direction discrepancy. (for reeco, this should already be the case for true)
  // True Momentum (closestZ for reco)
  TGraph* trueLinear_sigmaRatioVsSegmentMomentum = GetSigmaRatioGraph(trueLinear3D_sigmaVsSegmentMomentum, trueLinear_sigmaHLVsSegmentMomentum, "trueLinear_sigmaRatioVsSegmentMomentum", "");
  TGraph* recoLinear_sigmaRatioVsSegmentMomentum = GetSigmaRatioGraph(recoLinear3D_sigmaVsSegmentMomentum, recoLinear_sigmaRMSVsSegmentMomentum, "recoLinear_sigmaRatioVsSegmentMomentum", "");

  // BB Momentum
  // TODO: Note to self: Naming scheme had to change for BB compared with before this point, since the names of these objects don't say VsSegmentMomentum in them.  Perhaps they should for verbosity and consistency?
  TGraph* trueLinear_sigmaRatioVsSegmentBBMomentum = GetSigmaRatioGraph(trueLinear3D_sigmaVsSegmentBBMomentum, trueLinear_sigmaHLVsSegmentBBMomentum, "trueLinear_sigmaRatioVsSegmentBBMomentum", "");
  TGraph* recoLinear_sigmaRatioVsSegmentBBMomentum = GetSigmaRatioGraph(recoLinear3D_sigmaVsSegmentBBMomentum, recoLinear_sigmaRMSVsSegmentBBMomentum, "recoLinear_sigmaRatioVsSegmentBBMomentum", "");
 
  // Subsubsection 1.4.2: Split into two projection directions, due to direction discrepancy.  Should only need to do reco, but we're also diong true to cover all bases.
  // True Momentum (closestZ for reco)
  TGraph* trueLinearXZ_sigmaRatioVsSegmentMomentum = GetSigmaRatioGraph(trueLinear3D_sigmaVsSegmentMomentum, trueLinearXZ_sigmaHLVsSegmentMomentum, "trueLinearXZ_sigmaRatioVsSegmentMomentum", "");
  TGraph* trueLinearYZ_sigmaRatioVsSegmentMomentum = GetSigmaRatioGraph(trueLinear3D_sigmaVsSegmentMomentum, trueLinearYZ_sigmaHLVsSegmentMomentum, "trueLinearYZ_sigmaRatioVsSegmentMomentum", "");

  TGraph* recoLinearXZ_sigmaRatioVsSegmentMomentum = GetSigmaRatioGraph(recoLinear3D_sigmaVsSegmentMomentum, recoLinearXZ_sigmaRMSVsSegmentMomentum, "recoLinearXZ_sigmaRatioVsSegmentMomentum", "");
  TGraph* recoLinearYZ_sigmaRatioVsSegmentMomentum = GetSigmaRatioGraph(recoLinear3D_sigmaVsSegmentMomentum, recoLinearYZ_sigmaRMSVsSegmentMomentum, "recoLinearYZ_sigmaRatioVsSegmentMomentum", "");

  // BB Momentum
  TGraph* trueLinearXZ_sigmaRatioVsSegmentBBMomentum = GetSigmaRatioGraph(trueLinear3D_sigmaVsSegmentBBMomentum, trueLinearXZ_sigmaHLVsSegmentBBMomentum, "trueLinearXZ_sigmaRatioVsSegmentBBMomentum", "");
  TGraph* trueLinearYZ_sigmaRatioVsSegmentBBMomentum = GetSigmaRatioGraph(trueLinear3D_sigmaVsSegmentBBMomentum, trueLinearYZ_sigmaHLVsSegmentBBMomentum, "trueLinearYZ_sigmaRatioVsSegmentBBMomentum", "");

  TGraph* recoLinearXZ_sigmaRatioVsSegmentBBMomentum = GetSigmaRatioGraph(recoLinear3D_sigmaVsSegmentBBMomentum, recoLinearXZ_sigmaRMSVsSegmentBBMomentum, "recoLinearXZ_sigmaRatioVsSegmentBBMomentum", "");
  TGraph* recoLinearYZ_sigmaRatioVsSegmentBBMomentum = GetSigmaRatioGraph(recoLinear3D_sigmaVsSegmentBBMomentum, recoLinearYZ_sigmaRMSVsSegmentBBMomentum, "recoLinearYZ_sigmaRatioVsSegmentBBMomentum", "");

  // Print out sigma ratio mean values.  Ideally, they'd be the expected values, but they aren't quite.
  // True Momentum (closestZ for reco)
  std::cout << "trueLinear_sigmaRatio_mean = " << trueLinear_sigmaRatioVsSegmentMomentum->GetMean(2) << std::endl;
  std::cout << "recoLinear_sigmaRatio_mean = " << recoLinear_sigmaRatioVsSegmentMomentum->GetMean(2) << std::endl;

  std::cout << "trueLinearXZ_sigmaRatio_mean = " << trueLinearXZ_sigmaRatioVsSegmentMomentum->GetMean(2) << std::endl;
  std::cout << "trueLinearYZ_sigmaRatio_mean = " << trueLinearYZ_sigmaRatioVsSegmentMomentum->GetMean(2) << std::endl;

  std::cout << "recoLinearXZ_sigmaRatio_mean = " << recoLinearXZ_sigmaRatioVsSegmentMomentum->GetMean(2) << std::endl;
  std::cout << "recoLinearYZ_sigmaRatio_mean = " << recoLinearYZ_sigmaRatioVsSegmentMomentum->GetMean(2) << std::endl;

  // BB Momentum
  std::cout << "trueLinear_sigmaRatioBB_mean = " << trueLinear_sigmaRatioVsSegmentBBMomentum->GetMean(2) << std::endl;
  std::cout << "recoLinear_sigmaRatioBB_mean = " << recoLinear_sigmaRatioVsSegmentBBMomentum->GetMean(2) << std::endl;

  std::cout << "trueLinearXZ_sigmaRatioBB_mean = " << trueLinearXZ_sigmaRatioVsSegmentBBMomentum->GetMean(2) << std::endl;
  std::cout << "trueLinearYZ_sigmaRatioBB_mean = " << trueLinearYZ_sigmaRatioVsSegmentBBMomentum->GetMean(2) << std::endl;

  std::cout << "recoLinearXZ_sigmaRatioBB_mean = " << recoLinearXZ_sigmaRatioVsSegmentBBMomentum->GetMean(2) << std::endl;
  std::cout << "recoLinearYZ_sigmaRatioBB_mean = " << recoLinearYZ_sigmaRatioVsSegmentBBMomentum->GetMean(2) << std::endl;

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
  // True Momentum  
  TF1* trueLinear_sigmaHL_function = GetSigmaHL_function("trueLinear_sigmaHL_function", fitMin, fitMax);
  TFitResultPtr trueLinear_sigmaHL_fit = trueLinear_sigmaHLVsSegmentMomentum->Fit(trueLinear_sigmaHL_function, "QS0");

  // BB Momentum
  TF1* trueLinearBB_sigmaHL_function = GetSigmaHL_function("trueLinearBB_sigmaHL_function", fitMin, fitMax);
  TFitResultPtr trueLinearBB_sigmaHL_fit = trueLinear_sigmaHLVsSegmentBBMomentum->Fit(trueLinearBB_sigmaHL_function, "QS0"); 

  // Subsubsection 2.1.2: Fit sigma of individual angle directions.  These fit parameters are *not* used in other fits.
  // True Momentum
  TF1* trueLinearXZ_sigmaHL_function = GetSigmaHL_function("trueLinearXZ_sigmaHL_function", fitMin, fitMax);
  TFitResultPtr trueLinearXZ_sigmaHL_fit = trueLinearXZ_sigmaHLVsSegmentMomentum->Fit(trueLinearXZ_sigmaHL_function, "QS0");

  TF1* trueLinearYZ_sigmaHL_function = GetSigmaHL_function("trueLinearYZ_sigmaHL_function", fitMin, fitMax);
  TFitResultPtr trueLinearYZ_sigmaHL_fit = trueLinearYZ_sigmaHLVsSegmentMomentum->Fit(trueLinearYZ_sigmaHL_function, "QS0");

  // BB Momentum
  TF1* trueLinearXZBB_sigmaHL_function = GetSigmaHL_function("trueLinearXZBB_sigmaHL_function", fitMin, fitMax);
  TFitResultPtr trueLinearXZBB_sigmaHL_fit = trueLinearXZ_sigmaHLVsSegmentBBMomentum->Fit(trueLinearXZBB_sigmaHL_function, "QS0");

  TF1* trueLinearYZBB_sigmaHL_function = GetSigmaHL_function("trueLinearYZBB_sigmaHL_function", fitMin, fitMax);
  TFitResultPtr trueLinearYZBB_sigmaHL_fit = trueLinearYZ_sigmaHLVsSegmentBBMomentum->Fit(trueLinearYZBB_sigmaHL_function, "QS0");

  // Subsubsection 2.1.2: Define Fit Result Parameters
  // True Momentum
  Double_t trueLinear_kappa_a = trueLinear_sigmaHL_fit->Value(0);
  Double_t trueLinear_kappa_c = trueLinear_sigmaHL_fit->Value(1);

  Double_t trueLinearXZ_kappa_a = trueLinearXZ_sigmaHL_fit->Value(0);
  Double_t trueLinearXZ_kappa_c = trueLinearXZ_sigmaHL_fit->Value(1);

  // BB Momentum
  Double_t trueLinearBB_kappa_a = trueLinearBB_sigmaHL_fit->Value(0);
  Double_t trueLinearBB_kappa_c = trueLinearBB_sigmaHL_fit->Value(1);

  Double_t trueLinearXZBB_kappa_a = trueLinearXZBB_sigmaHL_fit->Value(0);
  Double_t trueLinearXZBB_kappa_c = trueLinearXZBB_sigmaHL_fit->Value(1);

  //
  // Subsection 2.2: True sigma3D fits.
  //

  // Subsubsection 2.2.1: Fit true sigma3D.
  // True Momentum
  TF1* trueLinear3D_sigmaHL_function = GetSigmaHL_3D_function("trueLinear3D_sigmaHL_function", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c);
  TFitResultPtr trueLinear3D_sigmaHL_fit = trueLinear3D_sigmaVsSegmentMomentum->Fit(trueLinear3D_sigmaHL_function, "QS0");

  // BB Momentum
  TF1* trueLinear3DBB_sigmaHL_function = GetSigmaHL_3D_function("trueLinear3DBB_sigmaHL_function", fitMin, fitMax, trueLinearBB_kappa_a, trueLinearBB_kappa_c);
  TFitResultPtr trueLinear3DBB_sigmaHL_fit = trueLinear3D_sigmaVsSegmentMomentum->Fit(trueLinear3DBB_sigmaHL_function, "QS0");

  // Subsubsection 2.2.2: Define Fit Result Parameters
  // True Momentum
  Double_t trueLinear_gamma3D = trueLinear3D_sigmaHL_fit->Value(3);

  // BB Momentum
  Double_t trueLinearBB_gamma3D = trueLinear3DBB_sigmaHL_fit->Value(3);

  //
  // Subsection 2.3: Reco sigmaRMS fits.
  //
  
  // Subsubsection 2.3.1: Fit sigma of combined angle directions.
  // True Momentum  
  TF1* recoLinear_sigmaRMS_function = GetSigmaRMS_function("recoLinear_sigmaRMS_function", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c);
  TFitResultPtr recoLinear_sigmaRMS_fit = recoLinear_sigmaRMSVsSegmentMomentum->Fit(recoLinear_sigmaRMS_function, "QS0");
  
  // BB Momentum
  TF1* recoLinearBB_sigmaRMS_function = GetSigmaRMS_function("recoLinearBB_sigmaRMS_function", fitMin, fitMax, trueLinearBB_kappa_a, trueLinearBB_kappa_c);
  TFitResultPtr recoLinearBB_sigmaRMS_fit = recoLinear_sigmaRMSVsSegmentBBMomentum->Fit(recoLinearBB_sigmaRMS_function, "QS0");

  // Subsubsection 2.3.2: Fit sigma of individual angle directions.
  // True Momentum
  TF1* recoLinearXZ_sigmaRMS_function = GetSigmaRMS_function("recoLinearXZ_sigmaRMS_function", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c);
  TFitResultPtr recoLinearXZ_sigmaRMS_fit = recoLinearXZ_sigmaRMSVsSegmentMomentum->Fit(recoLinearXZ_sigmaRMS_function, "QS0");

  TF1* recoLinearYZ_sigmaRMS_function = GetSigmaRMS_function("recoLinearYZ_sigmaRMS_function", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c);
  TFitResultPtr recoLinearYZ_sigmaRMS_fit = recoLinearYZ_sigmaRMSVsSegmentMomentum->Fit(recoLinearYZ_sigmaRMS_function, "QS0");

  // BB Momentum
  TF1* recoLinearXZBB_sigmaRMS_function = GetSigmaRMS_function("recoLinearXZBB_sigmaRMS_function", fitMin, fitMax, trueLinearBB_kappa_a, trueLinearBB_kappa_c);
  TFitResultPtr recoLinearXZBB_sigmaRMS_fit = recoLinearXZ_sigmaRMSVsSegmentBBMomentum->Fit(recoLinearXZBB_sigmaRMS_function, "QS0");

  TF1* recoLinearYZBB_sigmaRMS_function = GetSigmaRMS_function("recoLinearYZBB_sigmaRMS_function", fitMin, fitMax, trueLinearBB_kappa_a, trueLinearBB_kappa_c);
  TFitResultPtr recoLinearYZBB_sigmaRMS_fit = recoLinearYZ_sigmaRMSVsSegmentBBMomentum->Fit(recoLinearYZBB_sigmaRMS_function, "QS0");

  // Subsubsection 2.3.3: Define Fit Result Parameters
  
  // True Momentum
  Double_t recoLinear_sigmaRES = recoLinear_sigmaRMS_fit->Value(2);

  Double_t recoLinearXZ_sigmaRES = recoLinearXZ_sigmaRMS_fit->Value(2);
  Double_t recoLinearYZ_sigmaRES = recoLinearYZ_sigmaRMS_fit->Value(2);

  // BB Momentum
  Double_t recoLinearBB_sigmaRES = recoLinearBB_sigmaRMS_fit->Value(2);

  Double_t recoLinearXZBB_sigmaRES = recoLinearXZBB_sigmaRMS_fit->Value(2);
  Double_t recoLinearYZBB_sigmaRES = recoLinearYZBB_sigmaRMS_fit->Value(2);

  //
  // Subsection 2.4: Reco sigma3D fits.
  //

  // Subsubsection 2.4.1: Fit reco sigma3D
  // True Momentum
  TF1* recoLinear3D_sigmaRMS_function = GetSigmaRMS_3D_function("recoLinear3D_sigmaRMS_function", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c);
  TFitResultPtr recoLinear3D_sigmaRMS_fit = recoLinear3D_sigmaVsSegmentMomentum->Fit(recoLinear3D_sigmaRMS_function, "QS0");

  // BB Momentum
  TF1* recoLinear3DBB_sigmaRMS_function = GetSigmaRMS_3D_function("recoLinear3DBB_sigmaRMS_function", fitMin, fitMax, trueLinearBB_kappa_a, trueLinearBB_kappa_c);
  TFitResultPtr recoLinear3DBB_sigmaRMS_fit = recoLinear3D_sigmaVsSegmentBBMomentum->Fit(recoLinear3DBB_sigmaRMS_function, "QS0");

  // Subsubsection 2.4.2: Fit reco sigma3D, fixing either sigmaRES, to see if gamma3D is consistent
  // True Momentum
  TF1* recoLinear3D_sigmaRMS_function_fixSigmaRES = GetSigmaRMS_3D_function_fixSigmaRES("recoLinear3D_sigmaRMS_function_fixSigmaRES", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c, recoLinear_sigmaRES);
  TFitResultPtr recoLinear3D_sigmaRMS_fit_fixSigmaRES = recoLinear3D_sigmaVsSegmentMomentum->Fit(recoLinear3D_sigmaRMS_function_fixSigmaRES, "QS0");

  // BB Momentum
  TF1* recoLinear3DBB_sigmaRMS_function_fixSigmaRES = GetSigmaRMS_3D_function_fixSigmaRES("recoLinear3DBB_sigmaRMS_function_fixSigmaRES", fitMin, fitMax, trueLinearBB_kappa_a, trueLinearBB_kappa_c, recoLinearBB_sigmaRES);
  TFitResultPtr recoLinear3DBB_sigmaRMS_fit_fixSigmaRES = recoLinear3D_sigmaVsSegmentBBMomentum->Fit(recoLinear3DBB_sigmaRMS_function_fixSigmaRES, "QS0");

  // Subsubsection 2.4.3: Fit reco sigma3D, fixing gamma3D, to see if sigmaRES is consistent
  // True Momentum
  TF1* recoLinear3D_sigmaRMS_function_fixGamma3D = GetSigmaRMS_3D_function_fixGamma3D("recoLinear3D_sigmaRMS_function_fixGamma3D", fitMin, fitMax, trueLinear_kappa_a, trueLinear_kappa_c, trueLinear_gamma3D);
  TFitResultPtr recoLinear3D_sigmaRMS_fit_fixGamma3D = recoLinear3D_sigmaVsSegmentMomentum->Fit(recoLinear3D_sigmaRMS_function_fixGamma3D, "QS0");

  // BB Momentum
  TF1* recoLinear3DBB_sigmaRMS_function_fixGamma3D = GetSigmaRMS_3D_function_fixGamma3D("recoLinear3DBB_sigmaRMS_function_fixGamma3D", fitMin, fitMax, trueLinearBB_kappa_a, trueLinearBB_kappa_c, trueLinearBB_gamma3D);
  TFitResultPtr recoLinear3DBB_sigmaRMS_fit_fixGamma3D = recoLinear3D_sigmaVsSegmentBBMomentum->Fit(recoLinear3DBB_sigmaRMS_function_fixGamma3D, "QS0");

  // Subsubsection 2.4.4: Define Fit Result Parameters (results when parameters were not fixed)
  // True Momentum
  Double_t recoLinear3D_sigmaRES = recoLinear3D_sigmaRMS_fit->Value(2);
  Double_t recoLinear_gamma3D = recoLinear3D_sigmaRMS_fit->Value(3);

  // BB Momentum
  Double_t recoLinear3DBB_sigmaRES = recoLinear3DBB_sigmaRMS_fit->Value(2);
  Double_t recoLinearBB_gamma3D = recoLinear3DBB_sigmaRMS_fit->Value(3);

  // TODO: Define the final values from 2.4.2 & 2.4.3 (results when certain parameters were fixed)

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

  // TODO: Documentation
  // True Momentum
  plotAndSave_anglesVsSegmentMomentum(trueLinear_theta3DVsSegmentMomentum_HIST, "trueLinear3D", "FSA_plots/");
  plotAndSave_anglesVsSegmentMomentum(trueLinear_thetaPrimeVsSegmentMomentum_HIST, "trueLinear", "FSA_plots/");
  plotAndSave_anglesVsSegmentMomentum(trueLinear_thetaXZprimeVsSegmentMomentum_HIST, "trueLinearXZ", "FSA_plots/");
  plotAndSave_anglesVsSegmentMomentum(trueLinear_thetaYZprimeVsSegmentMomentum_HIST, "trueLinearYZ", "FSA_plots/");

  plotAndSave_anglesVsSegmentMomentum(recoLinear_theta3DVsSegmentMomentum_HIST, "recoLinear3D", "FSA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoLinear_thetaPrimeVsSegmentMomentum_HIST, "recoLinear", "FSA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoLinear_thetaXZprimeVsSegmentMomentum_HIST, "recoLinearXZ", "FSA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoLinear_thetaYZprimeVsSegmentMomentum_HIST, "recoLinearYZ", "FSA_plots/");

  // BB Momentum
  plotAndSave_anglesVsSegmentMomentum(trueLinear_theta3DVsSegmentBBMomentum_HIST, "trueLinear3DBB", "FSA_plots/");
  plotAndSave_anglesVsSegmentMomentum(trueLinear_thetaPrimeVsSegmentBBMomentum_HIST, "trueLinearBB", "FSA_plots/");
  plotAndSave_anglesVsSegmentMomentum(trueLinear_thetaXZprimeVsSegmentBBMomentum_HIST, "trueLinearXZBB", "FSA_plots/");
  plotAndSave_anglesVsSegmentMomentum(trueLinear_thetaYZprimeVsSegmentBBMomentum_HIST, "trueLinearYZBB", "FSA_plots/");

  plotAndSave_anglesVsSegmentMomentum(recoLinear_theta3DVsSegmentMomentum_HIST, "recoLinear3DBB", "FSA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoLinear_thetaPrimeVsSegmentMomentum_HIST, "recoLinearBB", "FSA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoLinear_thetaXZprimeVsSegmentMomentum_HIST, "recoLinearXZBB", "FSA_plots/");
  plotAndSave_anglesVsSegmentMomentum(recoLinear_thetaYZprimeVsSegmentMomentum_HIST, "recoLinearYZBB", "FSA_plots/");

  // new legends
  // True Momentum
  TLegend* trueLinear_sigmaHL_LEGEND = GetSigma_LEGEND(trueLinear_sigmaHL_fit, "sg");
  TLegend* trueLinearXZ_sigmaHL_LEGEND = GetSigma_LEGEND(trueLinearXZ_sigmaHL_fit, "sg");
  TLegend* trueLinearYZ_sigmaHL_LEGEND = GetSigma_LEGEND(trueLinearYZ_sigmaHL_fit, "sg");
  TLegend* trueLinear3D_sigmaHL_LEGEND = GetSigma_LEGEND(trueLinear3D_sigmaHL_fit, "acs")
;
  TLegend* recoLinear_sigmaRMS_LEGEND = GetSigma_LEGEND(recoLinear_sigmaRMS_fit, "acg");
  TLegend* recoLinearXZ_sigmaRMS_LEGEND = GetSigma_LEGEND(recoLinearXZ_sigmaRMS_fit, "acg");
  TLegend* recoLinearYZ_sigmaRMS_LEGEND = GetSigma_LEGEND(recoLinearYZ_sigmaRMS_fit, "acg");
  TLegend* recoLinear3D_sigmaRMS_LEGEND = GetSigma_LEGEND(recoLinear3D_sigmaRMS_fit, "ac");
  TLegend* recoLinear3D_sigmaRMS_fixSigmaRES_LEGEND = GetSigma_LEGEND(recoLinear3D_sigmaRMS_fit_fixSigmaRES, "acs");
  TLegend* recoLinear3D_sigmaRMS_fixGamma3D_LEGEND = GetSigma_LEGEND(recoLinear3D_sigmaRMS_fit_fixGamma3D, "acg");

  // BB Momentum
  TLegend* trueLinearBB_sigmaHL_LEGEND = GetSigma_LEGEND(trueLinearBB_sigmaHL_fit, "sg");
  TLegend* trueLinearXZBB_sigmaHL_LEGEND = GetSigma_LEGEND(trueLinearXZBB_sigmaHL_fit, "sg");
  TLegend* trueLinearYZBB_sigmaHL_LEGEND = GetSigma_LEGEND(trueLinearYZBB_sigmaHL_fit, "sg");
  TLegend* trueLinear3DBB_sigmaHL_LEGEND = GetSigma_LEGEND(trueLinear3DBB_sigmaHL_fit, "acs")
;
  TLegend* recoLinearBB_sigmaRMS_LEGEND = GetSigma_LEGEND(recoLinearBB_sigmaRMS_fit, "acg");
  TLegend* recoLinearXZBB_sigmaRMS_LEGEND = GetSigma_LEGEND(recoLinearXZBB_sigmaRMS_fit, "acg");
  TLegend* recoLinearYZBB_sigmaRMS_LEGEND = GetSigma_LEGEND(recoLinearYZBB_sigmaRMS_fit, "acg");
  TLegend* recoLinear3DBB_sigmaRMS_LEGEND = GetSigma_LEGEND(recoLinear3DBB_sigmaRMS_fit, "ac");
  TLegend* recoLinear3DBB_sigmaRMS_fixSigmaRES_LEGEND = GetSigma_LEGEND(recoLinear3DBB_sigmaRMS_fit_fixSigmaRES, "acs");
  TLegend* recoLinear3DBB_sigmaRMS_fixGamma3D_LEGEND = GetSigma_LEGEND(recoLinear3DBB_sigmaRMS_fit_fixGamma3D, "acg");

  // True Linear 2D
  // True Momentum
  plotAndSave_sigmaFitVsSegmentMomentum(trueLinear_sigmaHLVsSegmentMomentum, trueLinear_sigmaHL_function, trueLinear_sigmaHL_LEGEND, yMin, yMax, "trueLinear_sigmaHL", "FSA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(trueLinearXZ_sigmaHLVsSegmentMomentum, trueLinearXZ_sigmaHL_function, trueLinearXZ_sigmaHL_LEGEND, yMin, yMax, "trueLinearXZ_sigmaHL", "FSA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(trueLinearYZ_sigmaHLVsSegmentMomentum, trueLinearYZ_sigmaHL_function, trueLinearYZ_sigmaHL_LEGEND, yMin, yMax, "trueLinearYZ_sigmaHL", "FSA_plots/");

  // BB Momentum
  plotAndSave_sigmaFitVsSegmentMomentum(trueLinear_sigmaHLVsSegmentBBMomentum, trueLinearBB_sigmaHL_function, trueLinearBB_sigmaHL_LEGEND, yMin, yMax, "trueLinearBB_sigmaHL", "FSA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(trueLinearXZ_sigmaHLVsSegmentBBMomentum, trueLinearXZBB_sigmaHL_function, trueLinearXZBB_sigmaHL_LEGEND, yMin, yMax, "trueLinearXZBB_sigmaHL", "FSA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(trueLinearYZ_sigmaHLVsSegmentBBMomentum, trueLinearYZBB_sigmaHL_function, trueLinearYZBB_sigmaHL_LEGEND, yMin, yMax, "trueLinearYZBB_sigmaHL", "FSA_plots/");

  // True Linear/Polygonal 3D
  // True Momentum
  plotAndSave_sigmaFitVsSegmentMomentum(trueLinear3D_sigmaVsSegmentMomentum, trueLinear3D_sigmaHL_function, trueLinear3D_sigmaHL_LEGEND, yMin, yMax, "trueLinear3D_sigmaHL", "FSA_plots/");

  // BB Momentum
  plotAndSave_sigmaFitVsSegmentMomentum(trueLinear3D_sigmaVsSegmentBBMomentum, trueLinear3DBB_sigmaHL_function, trueLinear3DBB_sigmaHL_LEGEND, yMin, yMax, "trueLinear3DBB_sigmaHL", "FSA_plots/");

  // Reco Linear 2D
  // True Momentum
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear_sigmaRMSVsSegmentMomentum, recoLinear_sigmaRMS_function, recoLinear_sigmaRMS_LEGEND, yMin, yMax, "recoLinear_sigmaRMS", "FSA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinearXZ_sigmaRMSVsSegmentMomentum, recoLinearXZ_sigmaRMS_function, recoLinearXZ_sigmaRMS_LEGEND, yMin, yMax, "recoLinearXZ_sigmaRMS", "FSA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinearYZ_sigmaRMSVsSegmentMomentum, recoLinearYZ_sigmaRMS_function, recoLinearYZ_sigmaRMS_LEGEND, yMin, yMax, "recoLinearYZ_sigmaRMS", "FSA_plots/");

  // BB Momentum
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear_sigmaRMSVsSegmentBBMomentum, recoLinearBB_sigmaRMS_function, recoLinearBB_sigmaRMS_LEGEND, yMin, yMax, "recoLinearBB_sigmaRMS", "FSA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinearXZ_sigmaRMSVsSegmentBBMomentum, recoLinearXZBB_sigmaRMS_function, recoLinearXZBB_sigmaRMS_LEGEND, yMin, yMax, "recoLinearXZBB_sigmaRMS", "FSA_plots/");
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinearYZ_sigmaRMSVsSegmentBBMomentum, recoLinearYZBB_sigmaRMS_function, recoLinearYZBB_sigmaRMS_LEGEND, yMin, yMax, "recoLinearYZBB_sigmaRMS", "FSA_plots/");

  // Reco Linear/Polygonal 3D
  // True Momentum
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear3D_sigmaVsSegmentMomentum, recoLinear3D_sigmaRMS_function, recoLinear3D_sigmaRMS_LEGEND, yMin, yMax, "recoLinear3D_sigmaRMS", "FSA_plots/");

  // BB Momentum
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear3D_sigmaVsSegmentBBMomentum, recoLinear3DBB_sigmaRMS_function, recoLinear3DBB_sigmaRMS_LEGEND, yMin, yMax, "recoLinear3DBB_sigmaRMS", "FSA_plots/");

  // Reco Linear/Polygonal 3D - Fix sigmaRES
  // True Momentum
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear3D_sigmaVsSegmentMomentum, recoLinear3D_sigmaRMS_function_fixSigmaRES, recoLinear3D_sigmaRMS_fixSigmaRES_LEGEND, yMin, yMax, "recoLinear3D_fixSigmaRES_sigmaRMS", "FSA_plots/");

  // BB Momentum
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear3D_sigmaVsSegmentBBMomentum, recoLinear3DBB_sigmaRMS_function_fixSigmaRES, recoLinear3DBB_sigmaRMS_fixSigmaRES_LEGEND, yMin, yMax, "recoLinear3DBB_fixSigmaRES_sigmaRMS", "FSA_plots/");

  // Reco Linear/Polygonal 3D - Fix gamma3D
  // True Momentum
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear3D_sigmaVsSegmentMomentum, recoLinear3D_sigmaRMS_function_fixGamma3D, recoLinear3D_sigmaRMS_fixGamma3D_LEGEND, yMin, yMax, "recoLinear3D_fixGamma3D_sigmaRMS", "FSA_plots/");

  // BB Momentum
  plotAndSave_sigmaFitVsSegmentMomentum(recoLinear3D_sigmaVsSegmentBBMomentum, recoLinear3DBB_sigmaRMS_function_fixGamma3D, recoLinear3DBB_sigmaRMS_fixGamma3D_LEGEND, yMin, yMax, "recoLinear3DBB_fixGamma3D_sigmaRMS", "FSA_plots/");

  /////////////////////////////////////////////////////////////////////////////////////
  // Overlay sigmaHL and sigmaRMS plots. TODO: overlay theta3D plots too?
  /////////////////////////////////////////////////////////////////////////////////////
  // True Momentum
  plotAndSave_sigmaHL_sigmaRMS_overlay(trueLinear_sigmaHLVsSegmentMomentum, trueLinear_sigmaHL_function, recoLinear_sigmaRMSVsSegmentMomentum, recoLinear_sigmaRMS_function, yMin, yMax, "linear", "FSA_plots/");
  plotAndSave_sigmaHL_sigmaRMS_overlay(trueLinearXZ_sigmaHLVsSegmentMomentum, trueLinearXZ_sigmaHL_function, recoLinearXZ_sigmaRMSVsSegmentMomentum, recoLinearXZ_sigmaRMS_function, yMin, yMax, "linearXZ", "FSA_plots/");
  plotAndSave_sigmaHL_sigmaRMS_overlay(trueLinearYZ_sigmaHLVsSegmentMomentum, trueLinearYZ_sigmaHL_function, recoLinearYZ_sigmaRMSVsSegmentMomentum, recoLinearYZ_sigmaRMS_function, yMin, yMax, "linearYZ", "FSA_plots/");

  // BB Momentum
  plotAndSave_sigmaHL_sigmaRMS_overlay(trueLinear_sigmaHLVsSegmentBBMomentum, trueLinearBB_sigmaHL_function, recoLinear_sigmaRMSVsSegmentBBMomentum, recoLinearBB_sigmaRMS_function, yMin, yMax, "linearBB", "FSA_plots/");
  plotAndSave_sigmaHL_sigmaRMS_overlay(trueLinearXZ_sigmaHLVsSegmentBBMomentum, trueLinearXZBB_sigmaHL_function, recoLinearXZ_sigmaRMSVsSegmentBBMomentum, recoLinearXZBB_sigmaRMS_function, yMin, yMax, "linearXZBB", "FSA_plots/");
  plotAndSave_sigmaHL_sigmaRMS_overlay(trueLinearYZ_sigmaHLVsSegmentBBMomentum, trueLinearYZBB_sigmaHL_function, recoLinearYZ_sigmaRMSVsSegmentBBMomentum, recoLinearYZBB_sigmaRMS_function, yMin, yMax, "linearYZBB", "FSA_plots/");

  // True Momentum and BB Momentum overlay
  // sigmaHL
  plotAndSave_sigmaHL_trueMomentum_BBMomentum_overlay(trueLinear_sigmaHLVsSegmentMomentum, trueLinear_sigmaHL_function, trueLinear_sigmaHLVsSegmentBBMomentum, trueLinearBB_sigmaHL_function, yMin, yMax, "linear", "FSA_plots/");
  plotAndSave_sigmaHL_trueMomentum_BBMomentum_overlay(trueLinearXZ_sigmaHLVsSegmentMomentum, trueLinearXZ_sigmaHL_function, trueLinearXZ_sigmaHLVsSegmentBBMomentum, trueLinearXZBB_sigmaHL_function, yMin, yMax, "linearXZ", "FSA_plots/");
  plotAndSave_sigmaHL_trueMomentum_BBMomentum_overlay(trueLinearYZ_sigmaHLVsSegmentMomentum, trueLinearYZ_sigmaHL_function, trueLinearYZ_sigmaHLVsSegmentBBMomentum, trueLinearYZBB_sigmaHL_function, yMin, yMax, "linearYZ", "FSA_plots/");

  // sigmaRMS
  plotAndSave_sigmaRMS_trueMomentum_BBMomentum_overlay(recoLinear_sigmaRMSVsSegmentMomentum, recoLinear_sigmaRMS_function, recoLinear_sigmaRMSVsSegmentBBMomentum, recoLinearBB_sigmaRMS_function, yMin, yMax, "linear", "FSA_plots/");
  plotAndSave_sigmaRMS_trueMomentum_BBMomentum_overlay(recoLinearXZ_sigmaRMSVsSegmentMomentum, recoLinearXZ_sigmaRMS_function, recoLinearXZ_sigmaRMSVsSegmentBBMomentum, recoLinearXZBB_sigmaRMS_function, yMin, yMax, "linearXZ", "FSA_plots/");
  plotAndSave_sigmaRMS_trueMomentum_BBMomentum_overlay(recoLinearYZ_sigmaRMSVsSegmentMomentum, recoLinearYZ_sigmaRMS_function, recoLinearYZ_sigmaRMSVsSegmentBBMomentum, recoLinearYZBB_sigmaRMS_function, yMin, yMax, "linearYZ", "FSA_plots/");

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
