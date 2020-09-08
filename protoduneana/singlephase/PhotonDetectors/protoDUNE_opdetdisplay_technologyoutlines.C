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

    // TLine *apaline3=new TLine(0,0,710,0);
    // apaline3->SetLineColor(kBlack);
    // apaline3->Draw();

    TBox ara1(10.5,360,229.5,420);
    ara1.SetLineColor(3);
    ara1.SetLineStyle(10);
    ara1.SetFillStyle(0);
    ara1.Draw();

    TBox ara2(240.5,-360,459.5,-300);
    ara2.SetLineColor(3);
    ara2.SetLineStyle(10);
    ara2.SetFillStyle(0);
    ara2.Draw();

    TBox** APA5 = new TBox*[10];

    for(int i=0; i<10;i++){
      if(i==4)continue;
      APA5[i]=new TBox(240.5,-599+i*60,459.5,-541+i*60);
      APA5[i]->SetLineColor(kOrange+((i+2)%2)*10);
      APA5[i]->SetLineStyle(4+((i+2)%2));
      APA5[i]->SetFillStyle(0);
      APA5[i]->Draw();
    }

    TBox** APA6 = new TBox*[10];

    for(int i=0; i<10;i++){
     
      APA6[i]=new TBox(470.5,-599+i*60,689.5,-541+i*60);
      APA6[i]->SetLineColor(6+((i+2)%2));
      APA6[i]->SetLineStyle(2+((i+2)%2));
      APA6[i]->SetFillStyle(0);
      APA6[i]->Draw();
    }

    TBox** APA1 = new TBox*[10];

    for(int i=0; i<10;i++){
      if(i==6) continue;
      APA1[i]=new TBox(10.5,1+i*60,229.5,59+i*60);
      APA1[i]->SetLineColor(6+((i+2)%2));
      APA1[i]->SetLineStyle(2+((i+2)%2));
      APA1[i]->SetFillStyle(0);
      APA1[i]->Draw();
    }

    TBox** APA2 = new TBox*[10];

    for(int i=0; i<10;i++){
      if(i==6) continue;
      APA2[i]=new TBox(240.5,1+i*60,459.5,59+i*60);
      APA2[i]->SetLineColor(6+((i+2)%2));
      APA2[i]->SetLineStyle(2+((i+2)%2));
      APA2[i]->SetFillStyle(0);
      APA2[i]->Draw();
    }

    TBox** APA3 = new TBox*[10];

    for(int i=0; i<10;i++){
      if(i==6) continue;
      APA3[i]=new TBox(470.5,1+i*60,689.5,59+i*60);
      APA3[i]->SetLineColor(6+((i+2)%2));
      APA3[i]->SetLineStyle(2+((i+2)%2));
      APA3[i]->SetFillStyle(0);
      APA3[i]->Draw();
    }

    
    TBox** APA4 = new TBox*[10];

    for(int i=0; i<10;i++){
      APA4[i]=new TBox(10.5,-601+i*60,229.5,-539+i*60);
      if (i<4){
	APA4[i]->SetLineColor(6+((i+2)%2));
	APA4[i]->SetLineStyle(2+((i+2)%2));
      }
      else{
	APA4[i]->SetLineColor(kOrange+((i+2)%2)*10);
	APA4[i]->SetLineStyle(4+((i+2)%2));
      }
      APA4[i]->SetFillStyle(0);
      APA4[i]->Draw();
    }
    
    APA4[1]->SetLineColor(kOrange+10);
    APA4[1]->SetLineStyle(5);
    APA4[5]->SetLineColor(7);
    APA4[5]->SetLineStyle(3);

    // TBox YAxisCover(0, -20, -75, -618);
    // YAxisCover.SetFillColor(0);
    // YAxisCover.Draw();

    // TGaxis LowAxis(0,-600,0,0,0,600,3,"");
    // LowAxis.SetLabelSize(0.03);
    // LowAxis.Draw();

    
    
    c1.Update();
    outfile->cd();
    c1.Write();
    infile->cd();
    c1.Close();
  }
  
  infile->Close();
  outfile->Close();
}
