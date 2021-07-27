#include "HelperFunctions.cpp"
#include "PlottingFunctions.cpp"

TProfile* GetProfileX(TH2F* hist, const char* name, const char* title, const char* options);

void momentumAnalysis() {

  // TODO: save MCSMomentumVsTrueMomentum plots to canvases for each method.

  gSystem->Exec("rm -r MA_plots");
  gSystem->Exec("mkdir MA_plots");

  // Get Root File
  // const char* path = "/dune/app/users/hmeyer5/work_mcs_analysis/srcs/protoduneana/protoduneana/singlephase/MCSanalysis/CompareMCS_hist.root";
  const char* path = "../CompareMCS_hist.root";
  TFile* file = new TFile(path);

  file->ls();
  file->cd("comparemcsanalyzer");
  file->ls();

  // TODO: Contained vs. exiting frac bias & resolution plots.
  // TODO: TMC & TMF frac bias & resolution plots.

  // Get each MCS Momentum vs. True Momentum Histogram from root file.
  TH2F* trueLinear3D_MCSMomentumVsTrueMomentum_HIST = Get2DHist(file, "trueLinear3D_MCSMomentumVsTrueMomentum_HIST", "comparemcsanalyzer/");
  TH2F* truePolygonal3D_MCSMomentumVsTrueMomentum_HIST = Get2DHist(file, "truePolygonal3D_MCSMomentumVsTrueMomentum_HIST", "comparemcsanalyzer/");
  TH2F* trueLinear_MCSMomentumVsTrueMomentum_HIST = Get2DHist(file, "trueLinear_MCSMomentumVsTrueMomentum_HIST", "comparemcsanalyzer/");
  TH2F* truePolygonal_MCSMomentumVsTrueMomentum_HIST = Get2DHist(file, "truePolygonal_MCSMomentumVsTrueMomentum_HIST", "comparemcsanalyzer/");

  TH2F* recoLinear3D_MCSMomentumVsTrueMomentum_HIST = Get2DHist(file, "recoLinear3D_MCSMomentumVsTrueMomentum_HIST", "comparemcsanalyzer/");
  TH2F* recoPolygonal3D_MCSMomentumVsTrueMomentum_HIST = Get2DHist(file, "recoPolygonal3D_MCSMomentumVsTrueMomentum_HIST", "comparemcsanalyzer/");
  TH2F* recoLinear_MCSMomentumVsTrueMomentum_HIST = Get2DHist(file, "recoLinear_MCSMomentumVsTrueMomentum_HIST", "comparemcsanalyzer/");
  TH2F* recoPolygonal_MCSMomentumVsTrueMomentum_HIST = Get2DHist(file, "recoPolygonal_MCSMomentumVsTrueMomentum_HIST", "comparemcsanalyzer/");

  // Plot and save MCS Momentum vs. True Momentum histograms
  plotAndSave_MCSMomentumVsTrueMomentum(trueLinear3D_MCSMomentumVsTrueMomentum_HIST, "trueLinear3D", "MA_plots/");
  plotAndSave_MCSMomentumVsTrueMomentum(truePolygonal3D_MCSMomentumVsTrueMomentum_HIST, "truePolygonal3D", "MA_plots/");
  plotAndSave_MCSMomentumVsTrueMomentum(trueLinear_MCSMomentumVsTrueMomentum_HIST, "trueLinear", "MA_plots/");
  plotAndSave_MCSMomentumVsTrueMomentum(truePolygonal_MCSMomentumVsTrueMomentum_HIST, "truePolygonal", "MA_plots/");

  plotAndSave_MCSMomentumVsTrueMomentum(recoLinear3D_MCSMomentumVsTrueMomentum_HIST, "recoLinear3D", "MA_plots/");
  plotAndSave_MCSMomentumVsTrueMomentum(recoPolygonal3D_MCSMomentumVsTrueMomentum_HIST, "recoPolygonal3D", "MA_plots/");
  plotAndSave_MCSMomentumVsTrueMomentum(recoLinear_MCSMomentumVsTrueMomentum_HIST, "recoLinear", "MA_plots/");
  plotAndSave_MCSMomentumVsTrueMomentum(recoPolygonal_MCSMomentumVsTrueMomentum_HIST, "recoPolygonal", "MA_plots/");

  // Get each fractional difference histogram from root file
  TH2F* trueLinear3D_fracDiffVsTrueMomentum_HIST = Get2DHist(file, "trueLinear3D_fracDiffVsTrueMomentum_HIST", "comparemcsanalyzer/");
  TH2F* truePolygonal3D_fracDiffVsTrueMomentum_HIST = Get2DHist(file, "truePolygonal3D_fracDiffVsTrueMomentum_HIST", "comparemcsanalyzer/");
  TH2F* trueLinear_fracDiffVsTrueMomentum_HIST = Get2DHist(file, "trueLinear_fracDiffVsTrueMomentum_HIST", "comparemcsanalyzer/");
  TH2F* truePolygonal_fracDiffVsTrueMomentum_HIST = Get2DHist(file, "truePolygonal_fracDiffVsTrueMomentum_HIST", "comparemcsanalyzer/");

  TH2F* recoLinear3D_fracDiffVsTrueMomentum_HIST = Get2DHist(file, "recoLinear3D_fracDiffVsTrueMomentum_HIST", "comparemcsanalyzer/");
  TH2F* recoPolygonal3D_fracDiffVsTrueMomentum_HIST = Get2DHist(file, "recoPolygonal3D_fracDiffVsTrueMomentum_HIST", "comparemcsanalyzer/");
  TH2F* recoLinear_fracDiffVsTrueMomentum_HIST = Get2DHist(file, "recoLinear_fracDiffVsTrueMomentum_HIST", "comparemcsanalyzer/");
  TH2F* recoPolygonal_fracDiffVsTrueMomentum_HIST = Get2DHist(file, "recoPolygonal_fracDiffVsTrueMomentum_HIST", "comparemcsanalyzer/");

  // Rebin histograms
  trueLinear3D_fracDiffVsTrueMomentum_HIST->RebinX(25);
  truePolygonal3D_fracDiffVsTrueMomentum_HIST->RebinX(25);
  trueLinear_fracDiffVsTrueMomentum_HIST->RebinX(25);
  truePolygonal_fracDiffVsTrueMomentum_HIST->RebinX(25);
  recoLinear3D_fracDiffVsTrueMomentum_HIST->RebinX(25);
  recoPolygonal3D_fracDiffVsTrueMomentum_HIST->RebinX(25);
  recoLinear_fracDiffVsTrueMomentum_HIST->RebinX(25);
  recoPolygonal_fracDiffVsTrueMomentum_HIST->RebinX(25);

  // TODO: Consider allowing to set y-axis bounds to ensure that axes are consistent across plots.

  plotAndSave_fractionalBiasAndResolution(trueLinear3D_fracDiffVsTrueMomentum_HIST, "trueLinear3D", "True Linear 3D", "MA_plots/");
  plotAndSave_fractionalBiasAndResolution(truePolygonal3D_fracDiffVsTrueMomentum_HIST, "truePolygonal3D", "True Polygonal 3D", "MA_plots/");
  plotAndSave_fractionalBiasAndResolution(trueLinear_fracDiffVsTrueMomentum_HIST, "trueLinear", "True Linear", "MA_plots/");
  plotAndSave_fractionalBiasAndResolution(truePolygonal_fracDiffVsTrueMomentum_HIST, "truePolygonal", "True Polygonal", "MA_plots/");

  plotAndSave_fractionalBiasAndResolution(recoLinear3D_fracDiffVsTrueMomentum_HIST, "recoLinear3D", "Reco Linear 3D", "MA_plots/");
  plotAndSave_fractionalBiasAndResolution(recoPolygonal3D_fracDiffVsTrueMomentum_HIST, "recoPolygonal3D", "Reco Polygonal 3D", "MA_plots/");
  plotAndSave_fractionalBiasAndResolution(recoLinear_fracDiffVsTrueMomentum_HIST, "recoLinear", "Reco Linear", "MA_plots/");
  plotAndSave_fractionalBiasAndResolution(recoPolygonal_fracDiffVsTrueMomentum_HIST, "recoPolygonal", "Reco Polygonal", "MA_plots/");
}
