#ifndef PROTODUNEFITUTILS_h
#define PROTODUNEFITUTILS_h

#include <iostream>
#include <vector>
#include <string>

#include <TString.h>
#include <TTree.h>
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>

#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooWorkspace.h>
#include <RooFitResult.h>
#include <RooPlot.h>

namespace protoana{
  namespace ProtoDUNEFitUtils{

    // Get vector with all systematic histograms from file
    std::vector<TH1*> GetSystHistograms(std::string name);

    // Get histogram with MC statistical uncertainty
    TH1* GetStatsSystHistogram(TH1* nominal);

    // Get systematic histograms
    TH1* GetSystematicHistoFromNominal(TH1* nominal, TH1* syst, TString name);

    // Check if histogram has one or multiple bins
    bool IsSingleBinHisto(TH1* histo);

    // Get fit covariance matrix
    TH2* GetFitCovariance(RooFitResult* result);

    // Get fit correlation matrix
    TH2* GetFitCorrelationMatrix(RooFitResult* result);

    // Save/access workspace
    void SaveSnapshot(RooWorkspace* ws, TString snapshotname);
    bool LoadSnapshot(RooWorkspace* ws, TString snapshotname);
    void SaveWorkspace(RooWorkspace* ws, TString outFileName);

    // Set interpolation code
    void SetInterpolationCode(RooWorkspace* ws, int code);

    // Remove empty bins from RooPlot
    void RemoveEmptyBins(RooPlot* frame);

    // Calculate chi2 comparing the data and MC distributions
    double GetDataMCChi2(RooWorkspace *work, TString channelname, RooAbsData* data=NULL);

    // Vector of plots with data and pdfs
    std::vector<TCanvas*> PlotDatasetsAndPdfs(
        RooWorkspace *work, TString name, TString error, TString plottodraw,
        std::vector<TString> binnames, std::vector<double> recobins,
        std::vector<TString> incidentBinNames,
        std::vector<TString> sidebandBinNames, TString measurement="PDFit",
        bool doNegativeReco=false,RooAbsData* data=NULL,
        RooFitResult* result=NULL);

    std::vector<TH1 *> PlotXSecs(
        RooWorkspace * work, std::string name, /*std::string error,*/
        std::vector<TString> binnames, std::vector<double> recobins,
        std::vector<TString> incidentBinNames, RooAbsData * data = 0x0,
        RooFitResult * result = 0x0);
    
    // Vector of plots for the NLL
    std::vector<TCanvas*> PlotNLL(RooWorkspace *work, TString name, RooFitResult* result, bool plotPLL=false);

    // Plot the nuisance parameters impact on POI
    std::vector<TCanvas*> PlotNuisanceParametersImpact(RooWorkspace *work, RooFitResult* result, TString snapshotload, std::vector<std::string> SystsToConsider);

    // Plot pull mean and sigma of all the POI
    TCanvas* PlotParametersPull(TTree* tree, RooWorkspace* ws);

    // Plot nuisance parameters pulls
    TCanvas* PlotNuisanceParametersPull(TTree* tree, RooWorkspace* ws);

    // Plot nuisance parameters
    TCanvas* PlotNuisanceParameters(TTree* tree, RooWorkspace* ws);

    // Plot the average number of events from all toys. If fit the data the fit results will be shown
    TCanvas* PlotAverageResultsFromToys(TTree* tree, RooWorkspace* ws, TString channelname, TString catname);

    // Take list of poi
    RooArgList GetPostfitPOIList(RooArgList paramsfit, bool print=false);

    // Make all nuisance parameters in workspace constant. One exception is the exceptPar
    void MakeNuisanceParamsConstant(RooWorkspace* ws, TString exceptPar);
  
    // Reset all parameters
    void ResetValues(RooWorkspace* ws, const RooArgList& parList);
    
    // Reset all parameters to their nominal
    void ResetValuesToNominal(RooWorkspace* ws, const RooArgSet& parSet);

    // Reset error for all parameters
    void ResetError(RooWorkspace* ws, const RooArgList& parList);

    // Reset all values and errors
    void ResetAllValuesAndErrors(RooWorkspace* ws);

    // Get argon number density
    double GetArgonNumberDensity(double argon_density=1.3973, double argon_molecularmass=39.948);
    //                                                @87K

  }
}

#endif
