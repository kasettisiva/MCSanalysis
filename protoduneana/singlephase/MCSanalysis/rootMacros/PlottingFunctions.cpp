/*
  Functions used in plotting by other root macros.
 */

// Provide a Fit Result with 4 parameters: kappa_a, kappa_c, sigmaRES, gamma3D
// Options: Whichever options are included will correspond to which entry will fall under "fit parameters" or "fixed parameters"
//   a: kappa_a
//   c: kappa_c
//   s: sigmaRES
//   g: gamma3D
// Default Options are "".
TLegend* GetSigma_LEGEND(TFitResultPtr sigma_fit, const char* options) {
  // Define legend
  TLegend* legend = new TLegend(0.5, 0.55, 0.9, 0.9);
  legend->SetHeader("Fit Parameters", "C");
  legend->SetTextSize(0.035);

  // Define Legend Entry stringstreams
  std::stringstream kappa_a;
  kappa_a << std::fixed << std::setprecision(3) << "#kappa_{a} = " << sigma_fit->Value(0) << " +/- " << sigma_fit->Error(0) << " MeV";
  std::stringstream kappa_c;
  kappa_c << std::fixed << std::setprecision(3) << "#kappa_{c} = " << sigma_fit->Value(1) << " +/- " << sigma_fit->Error(1) << " MeV";
  std::stringstream sigmaRES;
  sigmaRES << std::fixed << std::setprecision(3) << "#sigma_{RES} = " << sigma_fit->Value(2) << " +/- " << sigma_fit->Error(2) << " mrad";
  std::stringstream gamma3D;
  gamma3D << std::fixed << std::setprecision(3) << "#gamma_{3D} = " << sigma_fit->Value(3) << " +/- " << sigma_fit->Error(3);

  std::string opts(options);

  // Add legend entries for fit parameters
  if(opts.find('a') == std::string::npos)
    legend->AddEntry((TObject*)0, kappa_a.str().c_str(), "");
  if(opts.find('c') == std::string::npos)
    legend->AddEntry((TObject*)0, kappa_c.str().c_str(), "");
  if(opts.find('s') == std::string::npos)
    legend->AddEntry((TObject*)0, sigmaRES.str().c_str(), "");
  if(opts.find('g') == std::string::npos)
    legend->AddEntry((TObject*)0, gamma3D.str().c_str(), "");
  
  /*
  // Add legend entries for fixed parameters
  if(options != "") // If there is at least one option, meaning at least one parameter is fixed.
    legend->AddEntry((TObject*)0, "Fixed Parameters:", "");
  if(opts.find('a') != std::string::npos)
    legend->AddEntry((TObject*)0, kappa_a.str().c_str(), "");
  if(opts.find('c') != std::string::npos)
    legend->AddEntry((TObject*)0, kappa_c.str().c_str(), "");
  if(opts.find('s') != std::string::npos)
   legend->AddEntry((TObject*)0, sigmaRES.str().c_str(), "");
  if(opts.find('g') != std::string::npos)
    legend->AddEntry((TObject*)0, gamma3D.str().c_str(), "");
  //*/
  std::stringstream chi2;
  chi2 << std::fixed << std::setprecision(3) << "#Chi^{2} = " << sigma_fit->Chi2();

  legend->AddEntry((TObject*)0, "Fit Info:", "");
  legend->AddEntry((TObject*)0, chi2.str().c_str(), "");

  // Return legend
  return legend;
}

// TODO: Documentation
void plotAndSave_sigmaFitVsSegmentMomentum(TGraphErrors* sigma_GRAPH, TF1* function, TLegend* legend, Double_t yMin, Double_t yMax, const char* namePrefix, const char* filePrefix) {
  TCanvas* canvas = new TCanvas("", "", 1500, 1500);

  sigma_GRAPH->GetYaxis()->SetRangeUser(yMin, yMax);
  sigma_GRAPH->Draw("AP");

  legend->Draw();

  function->Draw("same");

  std::stringstream fileName_ss;
  fileName_ss << filePrefix << namePrefix << "VsSegmentMomentum.png";

  canvas->SaveAs(fileName_ss.str().c_str());
}

// TODO: Documentation
void plotAndSave_anglesVsSegmentMomentum(TH2F* hist, const char* namePrefix, const char* filePrefix) {

  // Plot Angles vs. Segment Momentum Histogram
  // Overlay with Profile showing the stdev of each momentum bin.

  hist->SetStats(kFALSE);
  TCanvas* canvas = new TCanvas("", "", 1500, 1500);
  hist->Draw("colz");

  TProfile* profile = hist->ProfileX();
  profile->BuildOptions(1, -1, "s");
  profile->SetLineColor(kRed);
  profile->SetLineWidth(2);
  profile->Draw("same");

  std::stringstream fileName_ss;
  fileName_ss << filePrefix << namePrefix << "_anglesVsSegmentMomentum.png";

  canvas->SaveAs(fileName_ss.str().c_str());
}

// TODO: DOcumentation
void plotAndSave_MCSMomentumVsTrueMomentum(TH2F* hist, const char* namePrefix, const char* filePrefix) {
  TCanvas* canvas = new TCanvas("", "", 1500, 1500);
  hist->GetXaxis()->SetRangeUser(0, 5);
  hist->GetYaxis()->SetRangeUser(0,5);
  hist->SetStats(0);
  hist->Draw("box");
  TLine* line = new TLine(0, 0, 5, 5);
  line->SetLineColor(kRed);
  line->Draw("same");

  std::stringstream fileName_ss;
  fileName_ss << filePrefix << namePrefix << "_MCSMomentumVsTrueMomentum.png";

  canvas->SaveAs(fileName_ss.str().c_str());
}

// TODO: Documentation
void plotAndSave_fractionalBiasAndResolution(TH2F* hist, const char* namePrefix, const char* titlePrefix, const char* filePrefix) {

  std::stringstream fracBiasName_ss;
  fracBiasName_ss << namePrefix << "_fracBiasVsTrueMomentum_HIST";
  std::stringstream fracBiasTitle_ss;
  fracBiasTitle_ss << titlePrefix << " Fractional Bias vs. True Momentum; True Momentum (GeV/c); Fractional Bias";

  TGraph* line = new TGraph(2);
  line->SetLineStyle(kDashed);
  line->SetPoint(0, 0, 0);
  line->SetPoint(1, 5, 0);
  

  TProfile* fracBias_PROFILE = GetProfileX(hist, fracBiasName_ss.str().c_str(), fracBiasTitle_ss.str().c_str(), "");
  
  std::stringstream fracResolutionName_ss;
  fracResolutionName_ss << namePrefix << "_fracResolutionVsTrueMomentum_HIST";
  std::stringstream fracResolutionTitle_ss;
  fracResolutionTitle_ss << titlePrefix << " Fractional Resolution vs. True Momentum; True Momentum (GeV/c); Fractional Resolution";

  TGraphErrors* fracResolution_GRAPH = GetSigmaGraph(hist, fracResolutionName_ss.str().c_str(), fracResolutionTitle_ss.str().c_str());
  fracResolution_GRAPH->GetYaxis()->SetRangeUser(0.0, fracResolution_GRAPH->GetYaxis()->GetXmax());

  fracBias_PROFILE->GetXaxis()->SetRangeUser(fracResolution_GRAPH->GetXaxis()->GetXmin(), fracResolution_GRAPH->GetXaxis()->GetXmax());

  TCanvas* canvas = new TCanvas("", "", 1500, 1500);
  canvas->Divide(1, 2);
  canvas->cd(1);
  fracBias_PROFILE->Draw("E1X0");
  line->Draw("l");
  canvas->cd(2);
  fracResolution_GRAPH->Draw("AP");
  
  std::stringstream fileName_ss;
  fileName_ss << filePrefix << namePrefix << "_fractionalBiasAndResolutionVsTrueMomentum.png";

  canvas->SaveAs(fileName_ss.str().c_str());
}

void plotAndSave_sigmaHL_sigmaRMS_overlay(TGraphErrors* sigmaHL_GRAPH, TF1* sigmaHL_FUNC, TGraphErrors* sigmaRMS_GRAPH, TF1* sigmaRMS_FUNC, Double_t yMin, Double_t yMax, const char* namePrefix, const char* filePrefix) {
  TCanvas* canvas = new TCanvas("", "", 1500, 1500);

  sigmaHL_GRAPH->SetMarkerStyle(kFullCircle);
  sigmaHL_GRAPH->SetMarkerColor(kGreen);
  sigmaHL_GRAPH->SetLineColor(kOrange);

  TMultiGraph* multiGraph = new TMultiGraph();

  // Multigraph Title
  std::stringstream title;
  title << namePrefix << " #sigma_{HL} / #sigma_{RMS} vs. segment momentum overlay; Segment Momentum (GeV/c); #sigma";
  multiGraph->SetTitle(title.str().c_str());

  // Add graphs to multi graph
  multiGraph->Add(sigmaHL_GRAPH);
  multiGraph->Add(sigmaRMS_GRAPH);

  // For some reason... TMultiGraph sets the axis limits in different ways depending on the axis, doesn't work the same as normal TGraph.
  multiGraph->GetXaxis()->SetLimits(sigmaHL_GRAPH->GetXaxis()->GetXmin(), sigmaHL_GRAPH->GetXaxis()->GetXmax());
  multiGraph->SetMinimum(yMin);
  multiGraph->SetMaximum(yMax);

  // Draw graph and fit-functions
  multiGraph->Draw("AP");
  sigmaHL_FUNC->Draw("same");
  sigmaRMS_FUNC->Draw("same");

  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  // legend->SetTextSize(0.025);
  legend->AddEntry(sigmaHL_FUNC, "#sigma_{HL}");
  legend->AddEntry(sigmaRMS_FUNC, "#sigma_{RMS}");
  legend->Draw();

  // Save to file
  std::stringstream fileName;
  fileName << filePrefix << namePrefix << "_sigmaHL_sigmaRMS_overlay.png";
  canvas->SaveAs(fileName.str().c_str());
}

// TODO: Documentation
void plotAndSave_sigmaHL_trueMomentum_BBMomentum_overlay(TGraphErrors* sigmaHL_GRAPH, TF1* sigmaHL_FUNC, TGraphErrors* sigmaHLBB_GRAPH, TF1* sigmaHLBB_FUNC, Double_t yMin, Double_t yMax, const char* namePrefix, const char* filePrefix) {
  TCanvas* canvas = new TCanvas("", "", 1500, 1500);

  sigmaHL_GRAPH->SetMarkerStyle(kFullCircle);
  sigmaHL_GRAPH->SetMarkerColor(kOrange);
  sigmaHL_GRAPH->SetLineColor(kGreen);
  sigmaHL_FUNC->SetLineColor(kOrange);

  TMultiGraph* multiGraph = new TMultiGraph();

  // Multigraph Title
  std::stringstream title;
  title << namePrefix << " #sigma_{HL} vs. true momentum / BB momentum overlay; Momentum (GeV/c); #sigma";
  multiGraph->SetTitle(title.str().c_str());

  // Add graphs to multi graph
  multiGraph->Add(sigmaHL_GRAPH);
  multiGraph->Add(sigmaHLBB_GRAPH);

  // For some reason... TMultiGraph sets the axis limits in different ways depending on the axis, doesn't work the same as normal TGraph.
  multiGraph->GetXaxis()->SetLimits(sigmaHL_GRAPH->GetXaxis()->GetXmin(), sigmaHL_GRAPH->GetXaxis()->GetXmax());
  multiGraph->SetMinimum(yMin);
  multiGraph->SetMaximum(yMax);

  // Draw graph and fit-functions
  multiGraph->Draw("AP");
  sigmaHL_FUNC->Draw("same");
  sigmaHLBB_FUNC->Draw("same");

  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  // legend->SetTextSize(0.025);
  legend->AddEntry(sigmaHL_FUNC, "#sigma_{HL}");
  legend->AddEntry(sigmaHLBB_FUNC, "#sigma_{HL,BB}");
  legend->Draw();

  // Save to file
  std::stringstream fileName;
  fileName << filePrefix << namePrefix << "_sigmaHL_trueMomentum_BBMomentum_overlay.png";
  canvas->SaveAs(fileName.str().c_str());
}

// TODO: Documentation
void plotAndSave_sigmaRMS_trueMomentum_BBMomentum_overlay(TGraphErrors* sigmaRMS_GRAPH, TF1* sigmaRMS_FUNC, TGraphErrors* sigmaRMSBB_GRAPH, TF1* sigmaRMSBB_FUNC, Double_t yMin, Double_t yMax, const char* namePrefix, const char* filePrefix) {
  TCanvas* canvas = new TCanvas("", "", 1500, 1500);

  sigmaRMS_GRAPH->SetMarkerStyle(kFullCircle);
  sigmaRMS_GRAPH->SetMarkerColor(kOrange);
  sigmaRMS_GRAPH->SetLineColor(kGreen);
  sigmaRMS_FUNC->SetLineColor(kOrange);

  TMultiGraph* multiGraph = new TMultiGraph();

  // Multigraph Title
  std::stringstream title;
  title << namePrefix << " #sigma_{RMS} vs. true momentum / BB momentum overlay; Momentum (GeV/c); #sigma";
  multiGraph->SetTitle(title.str().c_str());

  // Add graphs to multi graph
  multiGraph->Add(sigmaRMS_GRAPH);
  multiGraph->Add(sigmaRMSBB_GRAPH);

  // For some reason... TMultiGraph sets the axis limits in different ways depending on the axis, doesn't work the same as normal TGraph.
  multiGraph->GetXaxis()->SetLimits(sigmaRMS_GRAPH->GetXaxis()->GetXmin(), sigmaRMS_GRAPH->GetXaxis()->GetXmax());
  multiGraph->SetMinimum(yMin);
  multiGraph->SetMaximum(yMax);

  // Draw graph and fit-functions
  multiGraph->Draw("AP");
  sigmaRMS_FUNC->Draw("same");
  sigmaRMSBB_FUNC->Draw("same");

  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  // legend->SetTextSize(0.025);
  legend->AddEntry(sigmaRMS_FUNC, "#sigma_{RMS}");
  legend->AddEntry(sigmaRMSBB_FUNC, "#sigma_{RMS,BB}");
  legend->Draw();

  // Save to file
  std::stringstream fileName;
  fileName << filePrefix << namePrefix << "_sigmaRMS_trueMomentum_BBMomentum_overlay.png";
  canvas->SaveAs(fileName.str().c_str());
}
