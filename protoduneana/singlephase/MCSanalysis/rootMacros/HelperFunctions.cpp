/*
  Generic Functions used by other root macros.
  These functions are all related to direct histogram manipulation.
 */

TGraphErrors* GetSigmaGraph(TH2F* hist2D, const char* name, const char* title);
TGraph* GetSigmaRatioGraph(TGraph* rawSigmaGraph, TGraph* projectedSigmaGraph, const char* name, const char* title);
TH2F* CombineHistograms(TH2F* hist1, TH2F* hist2, const char* name, const char* title);
TH2F* Get2DHist(TFile* file, const char* name, const char* prefix = "mcsAngleAnalysis/");

// Combine two TH2F* into one TH2F*
// Each of the TH2F* must have the same shape (# of x/y bins), though this is not checked.
// You'll probably just run into unexpected issues.
TH2F* CombineHistograms(TH2F* hist1, TH2F* hist2, const char* name, const char* title) {
  // Copy hist1 into a new histogram
  TH2F* newHistogram = new TH2F(*hist1);
  // Update the new histogram's name and title
  newHistogram->SetNameTitle(name, title);

  // This just uses the ROOT histogram adder.
  newHistogram->Add(hist2);
  // TODO: Add checks to see if the bins are the same size/shape
  // If not, then you should throw an error.
  /*
  // Loop through x-bins
  for(Int_t xBin = 0; xBin < hist2->GetNbinsX(); ++xBin) {
    Double_t x = hist2->GetXaxis()->GetBinCenter(xBin);
    for(Int_t yBin = 0; yBin < hist2->GetNbinsY(); ++yBin) {
      Double_t y = hist2->GetYaxis()->GetBinCenter(yBin);
      Double_t weight = hist2->GetBinContent(xBin, yBin);
      newHistogram->Fill(x, y, weight);
    }
  }
  */

  // Return newHistogram
  return newHistogram;
}

// Get a 2D hist from the file, with the specified prefix (default given for backwards compatibility.
TH2F* Get2DHist(TFile* file, const char* name, const char* prefix) {
  //const char* prefix = "mcsAngleAnalysis/";
  char* s = new char[strlen(prefix) + strlen(name) + 1];
  strcpy(s, prefix);
  strcat(s, name);
  return (TH2F*) file->Get(s);
}

// Get a TGraphErrors of sigma vs. x for the provided histogram, where the y-errors are the standard deviation of each x-bin of the histogram.
TGraphErrors* GetSigmaGraph(TH2F* hist2D, const char* name, const char* title) {
  // Define graph
  TGraphErrors* graph = new TGraphErrors();
  graph->SetNameTitle(name, title);

  // Define x and y axes
  TAxis* xAxis = hist2D->GetXaxis();
  TAxis* yAxis = hist2D->GetYaxis();

  // Loop through x-bins, to fill graph
  Int_t nXbins = hist2D->GetNbinsX();
  for(Int_t xBin = 0; xBin < nXbins; ++xBin) {    // Define x for this bin
    Double_t x = xAxis->GetBinCenter(xBin);

    // Define y-parameters for this x-bin (same for all x-bins)
    Int_t nYbins = hist2D->GetNbinsY();
    Double_t yMax = yAxis->GetXmax();
    Double_t yMin = yAxis->GetXmin();

    // Define 1D histogram for this x-bin
    TH1F* hist1D = new TH1F("", "", nYbins, yMin, yMax);
    //hist1D->GetXaxis()->Set(nYbins, yMin, yMax);
    // hist1D->GetXaxis()->SetRangeUser(yMin, yMax);
    //TH1F* hist = new TH1F("hist name", "hist title", nYbins, yMax, yMin);
    
    // Count is the number of entries within this x-bin, summed over all y-bins
    Int_t count = 0;
    // Loop through ybins for this x-bin to fill the 1D histogram
    for(Int_t yBin = 0; yBin < nYbins; ++yBin) {
      // Double_t y = yAxis->GetBinCenter(yBin);
      Double_t weight = hist2D->GetBinContent(xBin, yBin);
      // Double_t binError = hist2D->GetBinError(xBin, yBin);

      // Fill histogram with values
      // hist1D->Fill(y, weight);
      hist1D->SetBinContent(yBin, weight);
      // hist1D->SetBinError(yBin, binError);
      count += weight;
    } // Loop through y-bins

    // TODO: Add an option to get by fit.
    Double_t sigma = hist1D->GetStdDev();
    Double_t sigmaError = hist1D->GetStdDevError();
    
    // Fill Graph
    if(sigma != 0 /* && (x <= 2.0 && x >= 0.2) */ && x >= 0.2 && count >= 20) { // TODO: Add xMin function paramete

      // TFitResultPtr fit = hist1D->Fit("gaus","QS0","");
      // sigma = fit->Value(2);
      // sigmaError = fit->Error(2);

      Int_t pointNum = graph->GetN();
      graph->SetPoint(pointNum, x, sigma);
      graph->SetPointError(pointNum, 0, sigmaError);
    }
  }

  // Stylize graph
  graph->SetMarkerColor(kRed);
  graph->SetMarkerStyle(kFullCross);
  graph->SetLineColor(kBlue);

  // Return graph
  return graph;
}

// TODO: Documentation
TGraph* GetSigmaRatioGraph(TGraph* sigma3D_GRAPH, TGraph* sigma2D_GRAPH, const char* name, const char* title) {
  // TODO: Assert that they're the correct relative length.

  // Declare Graph
  TGraph* graph = new TGraph();
  graph->SetNameTitle(name, title);

  // Loop through points on each graph, take the ratio, and fill graph.
  Int_t numPoints = sigma3D_GRAPH->GetN();
  for(Int_t i = 0; i < numPoints; ++i) {
    Double_t segmentMomentum = sigma3D_GRAPH->GetPointX(i);
    Double_t sigma3D = sigma3D_GRAPH->GetPointY(i);
    Double_t sigma2D = sigma2D_GRAPH->GetPointY(i);
    graph->SetPoint(graph->GetN(), segmentMomentum, sigma3D / sigma2D); // We expect this to be of order sqrt(1-2/pi)*sqrt(2) ~ 0.8525 for all points.
  }

  // Stylize Graph
  graph->SetMarkerColor(kRed);
  graph->SetMarkerStyle(kFullCross);

  // Return the graph
  return graph;
}

// TODO: Documentation
TProfile* GetProfileX(TH2F* hist, const char* name, const char* title, const char* options) {
  TProfile* profile = hist->ProfileX(name, 0 , -1, options);
  profile->SetTitle(title);
  profile->SetLineColor(kRed);
  profile->SetStats(0);

  profile->GetXaxis()->SetRangeUser(0, 5);

  profile->SetMarkerStyle(kFullCircle);
  profile->GetXaxis()->SetLabelSize(0.06);
  profile->GetYaxis()->SetLabelSize(0.06);
  profile->GetXaxis()->SetTitleSize(0.05);
  profile->GetYaxis()->SetTitleSize(0.05);
  profile->GetYaxis()->SetTitleOffset(0.85);
  profile->GetXaxis()->SetTitleOffset(1.1);

  return profile;
}
