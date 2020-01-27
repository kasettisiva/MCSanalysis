#ifndef MCToyGenerationAndFit_h
#define MCToyGenerationAndFit_h

#include <iostream>
#include <vector>
#include <fstream>
#include <string>

// ROOT
#include "TTree.h"
#include "TString.h"

#include "RooAbsData.h"
#include "RooWorkspace.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooMultiVarGaussian.h"
#include "RooArgList.h"
#include "RooArgSet.h"

#include "RooStats/HypoTestResult.h"

namespace protoana{

  class MCToyGenerationAndFit{

  public:

    MCToyGenerationAndFit();
    MCToyGenerationAndFit(std::string minimizer, int fitstrategy, bool minoserror, double cl = 0.95);
    ~MCToyGenerationAndFit();

    // Generate a toy MC. Choose to generate either the expected events from data or MC
    RooAbsData* GenerateToyMC(RooWorkspace* ws, bool datanorm=true);

    // Generate and fit a number of nexp toys. Results will be written in a ttree  
    TTree* GenerateAndFit(RooWorkspace* ws, int nexp);

    // Fit the data
    RooFitResult* FitData(RooWorkspace* ws, bool isWeighted=false);
    
    // Fit Asimov dataset
    RooFitResult* FitAsimovData(RooWorkspace* w);
    
    // Fit MC toy dataset
    RooFitResult* FitToyData(RooWorkspace* w, RooAbsData* obsdata);

    // Fit toys from a list of workspaces
    TTree* FitToyMCFromWorkspace(std::vector<RooWorkspace*> wsvec, bool fitdata);

    // Write Roofitresult to a tree
    TTree* RooFitResultToTTree(RooWorkspace* ws, RooFitResult* res);

    // Write the fit results stored in RooMCStudy to a ttree
    TTree* RooMCStudyToTTree(RooMCStudy* mc);

    // Write a RooDataSet into a ttree
    TTree* RooDataSetToTTree(RooAbsData* data, TString treename);

    // Get chi2 between the fit models and a dataset. 
    RooArgSet* GetChi2Set(RooWorkspace *work, RooAbsData* data);

    // Generate tree with chi2 fit results
    TTree* GenerateChi2Tree(RooWorkspace* ws, int nexp, int seed);

    //HypoTestResult* GetHypothesisTest(RooWorkspace* ws, Int_t Ntoys, Int_t calctype, Int_t stattype);

    // Option to change the minimiser - default is Minuit2
    void SetMinimiser(std::string min) {_minimizer = TString(min.c_str());}

    // Change algorithm - default is migrad
    void SetAlgorithm(std::string min) {_algorithm = TString(min.c_str());}

    // Option to change fit strategy - default is 1
    void SetFitStrategy(int fs) {_fitstrategy = fs;}

    // Option to enable Minos error analysis - default is off
    void EnableMinosError() {_minoserror = true;}

    // Option to set CL - default is 95%
    void SetConfLevel(double cl) {_conflevel = cl;}
    
  protected:
    
    // Minimiser
    TString _minimizer;

    // Algorithm
    TString _algorithm;
    
    // Fit strategy
    int _fitstrategy;
    
    // Flag to enable minos error analysis
    bool _minoserror;

    // Confidence level
    double _conflevel;
    
  };

}

#endif
