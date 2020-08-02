#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TProfile.h"
#include "TH1.h"
#include "TFile.h"
#include "TH3.h"
#include "TObject.h"
#include <TLine.h>
#include "TAxis.h"

void OpVis_peripheral()
{
  //TFile *f1= new TFile("side_effs.root");
  TFile *infile = new TFile("opvis_hist.root");
  TFile *outfile = new TFile("opvis_canvases.root","recreate");
  infile->cd("opvisana");

  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  //loop counter(s)
  int i=0;
  int j=0;
 
  TIter nextkey(gDirectory->GetListOfKeys());
  while (auto key = (TKey*)nextkey()) {
    ++i;
    cout << key->GetName() << " " << key->GetClassName() << endl;
    auto obj = key->ReadObj();
    const char* name = obj->GetName();
    TCanvas c1(name,name);
    obj->Draw("colz");

    TLine middleline(0,0,700,0);
    middleline.SetLineWidth(2);
    middleline.SetLineColor(kBlack);
    middleline.Draw();

    TLine apaline1(235,600,235,-600);
    apaline1.SetLineWidth(2);
    apaline1.SetLineColor(kBlack);
    apaline1.Draw();

    TLine apaline2(465,600,465,-600);
    apaline2.SetLineWidth(2);
    apaline2.SetLineColor(kBlack);
    apaline2.Draw();
   
    c1.Update();
    outfile->cd();
    c1.Write();
    infile->cd();
    c1.Close();
  }
  
  infile->Close();
  outfile->Close();
}
