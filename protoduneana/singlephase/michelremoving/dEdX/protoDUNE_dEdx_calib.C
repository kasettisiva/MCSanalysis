#define protoDUNE_dEdx_calib_cxx
#include "protoDUNE_dEdx_calib.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <string>
#include <TImage.h>
#include <iomanip>
#include <TSpline.h>
#include <TText.h>
#include <TFrame.h>
#include <TMinuit.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"
#include "cetlib/search_path.h"
#include "cetlib/filesystem.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

using namespace std;

int hitplane=1;
bool usemap=true;
//const double calib_factor =4.50e-3; // ******************* change the calibraion factor here with the normalisation ****************** //
double calib_factor =1.029e-3;
double pitchvalue=0.65;
double Rho = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia) 
double betap = 0.212;//(kV/cm)(g/cm^2)/MeV
double alpha = 0.93;//parameter from ArgoNeuT experiment at 0.481kV/cm 
double Wion = 23.6e-6;//parameter from ArgoNeuT experiment at 0.481kV/cm
double Efield = 0.50;//kV/cm protoDUNE electric filed
double normalisation_factor[3]={0.992,0.990,0.989};//for plane 0
const int Z=18; //Atomic number of Argon
const double A=39.948; // g/mol Atomic mass of Argon
const double I=188.0e-6; // ev
const double K=0.307; // Mev.cm^2 / mol
const double Mmu=105.658; // Mev for Mu
const double Me=0.51; // Mev for electron
const double rho=1.396;//g/cm^3
string outfile_name;
TString mn = "2";

double spline_KE[13];
double spline_Range[13];
TFile * OpenFile(const std::string filename);

////getting the variable Efield using data driven maps
TFile *ef/*=new TFile("SCE_DataDriven_180kV_v4.root")*/;
TH3F *xneg=(TH3F*)ef->Get("Reco_ElecField_X_Neg");
TH3F *yneg=(TH3F*)ef->Get("Reco_ElecField_Y_Neg");
TH3F *zneg=(TH3F*)ef->Get("Reco_ElecField_Z_Neg");
TH3F *xpos=(TH3F*)ef->Get("Reco_ElecField_X_Pos");
TH3F *ypos=(TH3F*)ef->Get("Reco_ElecField_Y_Pos");
TH3F *zpos=(TH3F*)ef->Get("Reco_ElecField_Z_Pos");
 float tot_Ef(float xval,float yval,float zval){
  float E0value=0.4867;
  if(!usemap) return E0value;
  if(xval>=0){
    float ex=E0value+E0value*xpos->Interpolate(xval,yval,zval);
    float ey=E0value*ypos->Interpolate(xval,yval,zval);
    float ez=E0value*zpos->Interpolate(xval,yval,zval);
    return sqrt(ex*ex+ey*ey+ez*ez);
  }
if(xval<0){
  float ex=E0value+E0value*xneg->Interpolate(xval,yval,zval);
    float ey=E0value*yneg->Interpolate(xval,yval,zval);
    float ez=E0value*zneg->Interpolate(xval,yval,zval);
    return sqrt(ex*ex+ey*ey+ez*ez);
  }
 else return E0value;
}

////********************************************///


double beta(double gamma){
  double value=TMath::Sqrt(1-(1.0/(gamma*gamma)));
  return value;
}

double gamma(double KE,double mass){
  double value=(double(KE)/mass)+1;
  return value;
}

double KE=266.;
double g=gamma(KE,Mmu);
double b=beta(g);

const double C=-5.2146;
const double X0=0.2;
const double X1=3.0;
const double a=0.19559;
const double m=3.0;
const double N=2*TMath::Log(10);

double density(double bg){//replaced x by x1
  double value;
  double x = TMath::Log10(bg);
  if(x<X0) return 0;
  if(x>X1) return N*x + C;
  value=a*(TMath::Power((X1-x),m));
  return N*x + C + value;
}

double Wmax(double KE,double mass){
  double g=gamma(KE,mass);
  double b=beta(g);
  double num=2*Me*(TMath::Power(b*g,2));
  double den=1+double(2*g*Me)/mass + TMath::Power((double(Me)/mass),2);
  double value=double(num)/den;
  return value;
}

double dEdx(double KE,double mass){
  double g=gamma(KE,mass);
  double b=beta(g);
  double f=K*(double(Z)/A)*(TMath::Power(1.0/b,2));
  double wmax=Wmax(KE,mass);
  double a0=0.5*(TMath::Log(double(2*Me*(TMath::Power(b*g,2))*wmax)/(I*I)));
  double dens=density(b*g);
  double value=f*rho*(a0-b*b-double(dens)/2);
  return value;
}

double dpdx(double KE,double x,double mass){
  double g=gamma(KE,mass);
  double b=beta(g);
  double epsilon=(double(K)/2)*(double(Z)/A)*(double(x*rho)/(b*b));
  double A0=double(2*Me*(TMath::Power((b*g),2)))/I;
  double A1=double(epsilon)/I;
  double value=(1.0/x)*epsilon*((TMath::Log(A0)) + TMath::Log(A1) + 0.2 - b*b - density(b*g));
  return value;
}

///////////////////// End of Landau-Vavilov function /////////////////////////////////////////////

///////////////////////////////////// Function definition ////////////////////////////////

Double_t langaufun(Double_t *x, Double_t *par) {
  Double_t invsq2pi = 0.398942280401;// Control constants
  //Double_t mpshift = -0.22278298; 
  Double_t np = 500.0;
  Double_t sc = 5.0;// convolution extends to +-sc Gaussian sigmas   
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;

  //mpc = par[1]- mpshift * par[0];
  mpc=par[1];
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow)/np;
 
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

//////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////// Function definition ////////////////////////////////

 TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, Int_t *Status)
// TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
  Int_t i;
  Char_t FunName[100];

  sprintf(FunName,"Fitfcn_%s",his->GetName());

  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MPV","Area","GSigma");
   
  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

//  his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
  TFitResultPtr fitres = his->Fit(FunName,"RBOS"); // fit within specified range, use ParLimits, do not plot /////////////////// Initial code use the mode "RBO" (commented by VARUNA) ///////////

  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
    fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }
 
  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf
  Status[0] = fitres->CovMatrixStatus();
 
  return (ffit);              // return fit function
}

/////////////////////////////// Function definition /////////////////////////////////////

Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {
  Double_t p,x,fy,fxr,fxl;
  Double_t step;
  Double_t l,lold;
  Int_t i = 0;
  Int_t MAXCALLS = 10000;

  // Search for maximum

  p = params[1] - 0.1 * params[0];
  step = 0.05 * params[0];
  lold = -2.0;
  l = -1.0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = langaufun(&x,params);
    if (l < lold)
      step = -step/10;
    p += step;
  }

  if (i == MAXCALLS)
    return (-1);

  maxx = x;
  fy = l/2;

  // Search for right x location of fy

  p = maxx + params[0];
  step = params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
      step = -step/10;
    p += step;
  }

  if (i == MAXCALLS)
    return (-2);

  fxr = x;

  // Search for left x location of fy

  p = maxx - 0.5 * params[0];
  step = -params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
      step = -step/10;
    p += step;
  }

  if (i == MAXCALLS)
    return (-3);

  fxl = x;
  FWHM = fxr - fxl;
  return (0);
}

////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////// Function definition ////////////////////////////////

Float_t Dedx(float dqdx, float Ef){
  return (exp(dqdx*(betap/(Rho*Ef)*Wion))-alpha)/(betap/(Rho*Ef));
}

//////////////////////////////////////////////////////////////////////////////////////


void protoDUNE_dEdx_calib::Loop()
{
  std::cout << "******************************* Calibration.C is running *******************************" << std::endl;
  TH2F *fhist_dedx= new TH2F("dedx_vs_residual_range",Form("plane %d calibrated dE/dx vs residual range;residual range in cm;dEdx in MeV/cm",hitplane),200,0,200,50,0,5);
  TH1F *fhist_trkpitch = new TH1F("trkpitch",Form("trkpitch plane %d;trkpitch in cm",hitplane),30,0,3);
  TH1F *dqdx_rat=new TH1F("dqdx_ratio",Form("plane %d",hitplane),100,0,10);
  TH1F *my_hist = new TH1F("my_hist",Form("plane %d;Kinetic Energy (MeV);MPV dEdx (MeV/cm)",hitplane),10,0,500);
  TH2F *fhist_dqdxcal= new TH2F("fhist_dqdxcal",Form("plane %d calibrated dqdx vs residual range;residual range in cm;dqdx in ADC/cm",hitplane),200,0,200,500,0,1000);
  TH2F *fhist_dqdxuncal= new TH2F("fhist_dqdxuncal",Form("plane %d uncalibrated dq/dx vs residual range;residual range in cm;dqdx in ADC/cm",hitplane),200,0,200,500,0,1000);
  TH1F *hdedx= new TH1F("hdedx",Form("plane %d calibrated dE/dx;dEdx in MeV/cm;entries",hitplane),200,0,10);
  hdedx->Sumw2();
  TH1F *dqdx_uncal=new TH1F("dqdx_uncal","Uncalibrated dQ/dx;dQ/dx in ADC/cm;no of entries",500,0,200e3);
  TH1F *dqdx_cal=new TH1F("dqdx_cal","calibrated dQ/dx;dQ/dx in ADC/cm;no of entries",500,0,200e3);                                                                    
  my_hist->SetLineColor(kCyan);
  my_hist->SetFillColor(kCyan);
  my_hist->SetBinContent(6,10);
  my_hist->SetBinContent(7,10);
  my_hist->SetBinContent(8,10);
  my_hist->SetBinContent(9,10);
  my_hist->GetYaxis()->SetRangeUser(0.0,6.0); // what does this line mean?
  ofstream outfile;
  outfile.open(Form("MCsceon_plane_%d.txt",hitplane),std::ios_base::app);
  ofstream myfile;
  myfile.open(Form("trkinfo_%d.txt",hitplane));
  int nbin=40;
  //int nbin=20;
  int binsize=5;
  //int  binsize=5; // cm
  TH1D *dedx[nbin];
  //TCanvas *d[nbin];
  for (int i=0; i<nbin; ++i){
    if(i==0) dedx[i] = new TH1D(Form("dedx_%d",i),Form("dedx_%d",i),300,0.0,15); 
    if(i!=0) dedx[i] = new TH1D(Form("dedx_%d",i),Form("dedx_%d",i),200,0.0,10);
    dedx[i]->SetLineColor(kBlack); 
    dedx[i]->Sumw2();
  }
 
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  // gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");
 
  ///////////////// Make any changes to the Y and Z bin sizes here ///////////////

  //int x_bin_size = 5; // 148 bins in x direction
  fChain->GetEntry(0); 

  /////////////////////Importing X fractional corrections//////////////////////

  TFile my_file0(Form("YZcalo_mich%s_r%d.root",mn.Data(), run));
  TFile my_file2(Form("Xcalo_mich%s_r%d.root",mn.Data(), run));
  TH1F *X_correction_hist = (TH1F*)my_file2.Get(Form("dqdx_X_correction_hist_%d",hitplane));
  TH2F *YZ_correction_neg_hist=(TH2F*)my_file0.Get(Form("correction_dqdx_ZvsY_negativeX_hist_%d",hitplane));
  TH2F *YZ_correction_pos_hist=(TH2F*)my_file0.Get(Form("correction_dqdx_ZvsY_positiveX_hist_%d",hitplane));
 
  TFile *file = new TFile(Form("dedx_plane_%d_r%d.root",hitplane,run),"recreate");
 
 
  ////////////////////////////////////////////////////////////////////////////////// 
  int total_trks=0;
  int used_trks=0;
 
 if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//   for (Long64_t jentry=0; jentry<10000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry%100==0) cout<<jentry<<"/"<<nentries<<endl;
    vector<float> res, dq, first5dq, last5dq;
    for(int i=0; i<cross_trks; ++i){
      total_trks++; //total tracks filtered by the module
      if(trkstartx[i]*trkendx[i]>0) continue;
      /***************finding correct endx*****************************************/
      //float endx = trkendx[i];
      //if(trkendy[i]<trkstarty[i]) endx=trkendx[i];
      //if(trkendy[i]>trkstarty[i]) endx=trkstartx[i];
    
      if(peakT_min[i]<100||peakT_max[i]>5900||trklen[i]<100||trklen[i]>700||(trkendz[i]>226 && trkendz[i]<236)||(trkstartz[i]>226 && trkstartz[i]<236)||(trkendz[i]>456 && trkendz[i]<472)||(trkstartz[i]>456 && trkstartz[i]<472)) continue;//filter for plane 2
      if(hitplane==2 && ((TMath::Abs(trackthetaxz[i])>1.13 && TMath::Abs(trackthetaxz[i])<2.0)||(TMath::Abs(trackthetayz[i])>1.22 && TMath::Abs(trackthetayz[i])<1.92))) continue;
      //if(!((TMath::Abs(trkstartx[i])>350||trkstarty[i]<50||trkstarty[i]>550||trkstartz[i]<50||trkstartz[i]>645)&&(TMath::Abs(trkendx[i])>350||trkendy[i]<50||trkendy[i]>550||trkendz[i]<50||trkendz[i]>645))) continue;

      if(adjacent_hits[i]!=0 || dist_min[i]>5) continue;
      if(lastwire[i]<=5 || lastwire[i]>=475) continue;  //for plane 2  
 

      for(int j=0;j<ntrkhits[i][hitplane];j++){
	res.push_back(trkresrange[i][hitplane][j]);
	dq.push_back(trkdqdx[i][hitplane][j]);
      }//ntrkhits
      /**********************end of buffer filling*****************************/

      if(res.size()==0) continue;
      /***********************removed empty tracks to avoid segmentation fault******************/

      int siz1=ntrkhits[i][hitplane];
      float max=*max_element(res.begin(),res.end());
      //float min=*min_element(res.begin(),res.end());
      //int siz=res.size();
     
      /************************flipping wrongly ordered residual range values****************************/
      bool test=true;
      if((trkhity[i][hitplane][siz1-1]<trkhity[i][hitplane][0] && trkresrange[i][hitplane][siz1-1]>trkresrange[i][hitplane][0])||(trkhity[i][hitplane][0]<trkhity[i][hitplane][siz1-1] && trkresrange[i][hitplane][0]>trkresrange[i][hitplane][siz1-1])){
	test=false;
	for(int i1=0;i1<ntrkhits[i][hitplane];i1++){
	  trkresrange[i][hitplane][i1]=res[siz1-i1-1];
	}
	cout<<"This is a flipped track"<<endl;
      }
     
      /***************calculating the ratio of dQdx for first 5cm and last 5 cm of a track********************/ 
      for(size_t k=0;k<res.size();k++){
	if(trkresrange[i][hitplane][k]<5) first5dq.push_back(dq[k]);
	if(trkresrange[i][hitplane][k]>max-5) last5dq.push_back(dq[k]);
      }
     
      if(first5dq.size()<5){
	res.erase(res.begin(),res.end());
	dq.erase(dq.begin(),dq.end());
	first5dq.clear();
	last5dq.clear();
	continue;
	}

      
	 float med1 = TMath::Median(first5dq.size(), &first5dq[0]);
	 float med2= TMath::Median(last5dq.size(), &last5dq[0]);
	 dqdx_rat->Fill(med1/med2);
	 res.erase(res.begin(),res.end());
	 first5dq.erase(first5dq.begin(),first5dq.end());
	 last5dq.erase(last5dq.begin(),last5dq.end());
	 dq.erase(dq.begin(),dq.end());
	   if(!((med1/med2)>1.4)) continue;
	   if(!test) continue;
	 used_trks++;


     
      for(int j=0; j<TMath::Min(ntrkhits[i][hitplane],3000); ++j){
	fhist_trkpitch->Fill(trkpitch[i][hitplane][j]);
	if(trkpitch[i][hitplane][j]>=0.5 && trkpitch[i][hitplane][j]<=0.8){
	  if(trkhity[i][hitplane][j]>0 && trkhity[i][hitplane][j]<600){
	    if(trkhitz[i][hitplane][j]>0 && trkhitz[i][hitplane][j]<695){
	      if(trkhitx[i][hitplane][j]>-360 && trkhitx[i][hitplane][j]<0){ //negative X direction
		if(hitplane == 1 && abs(180/3.14*trackthetaxz[i])>130 && !(abs(180/3.14*trackthetayz[i])>80 && abs(180/3.14*trackthetayz[i])<100)){ //plane 1
		int x_bin=X_correction_hist->FindBin(trkhitx[i][hitplane][j]);
                float Cx=X_correction_hist->GetBinContent(x_bin);
                float Cyzneg=YZ_correction_neg_hist->GetBinContent(YZ_correction_neg_hist->FindBin(trkhitz[i][hitplane][j],trkhity[i][hitplane][j]));
                float corrected_dq_dx=trkdqdx[i][hitplane][j]*Cx*normalisation_factor[hitplane]*Cyzneg;
                float scaled_corrected_dq_dx=float(corrected_dq_dx)/calib_factor;
		float cal_de_dx=Dedx(scaled_corrected_dq_dx,tot_Ef(trkhitx[i][hitplane][j],trkhity[i][hitplane][j],trkhitz[i][hitplane][j]));
//		float cal_de_dx=-1;
//     		if(usemap) float cal_de_dx=Dedx(scaled_corrected_dq_dx,tot_Ef(trkhitx[i][hitplane][j],trkhity[i][hitplane][j],trkhitz[i][hitplane][j]));
//                else float cal_de_dx=Dedx(scaled_corrected_dq_dx,0.4867);
                hdedx->Fill(cal_de_dx);
                int bin=int(trkresrange[i][hitplane][j])/binsize;
                fhist_dedx->Fill(trkresrange[i][hitplane][j],cal_de_dx);
                fhist_dqdxcal->Fill(trkresrange[i][hitplane][j],corrected_dq_dx);
                fhist_dqdxuncal->Fill(trkresrange[i][hitplane][j],trkdqdx[i][hitplane][j]);
                dqdx_cal->Fill(corrected_dq_dx/calib_factor);
                dqdx_uncal->Fill(trkdqdx[i][hitplane][j]/calib_factor);
		 if(bin<nbin){
                  dedx[bin]->Fill(cal_de_dx);
                 } // bin <40
		}
		else if(hitplane == 0 && abs(180/3.14*trackthetaxz[i])<40 && !(abs(180/3.14*trackthetayz[i])>80 && abs(180/3.14*trackthetayz[i])<100)){ //plane 0
		int x_bin=X_correction_hist->FindBin(trkhitx[i][hitplane][j]);
		float Cx=X_correction_hist->GetBinContent(x_bin);
		float Cyzneg=YZ_correction_neg_hist->GetBinContent(YZ_correction_neg_hist->FindBin(trkhitz[i][hitplane][j],trkhity[i][hitplane][j]));
		float corrected_dq_dx=trkdqdx[i][hitplane][j]*Cx*normalisation_factor[hitplane]*Cyzneg;
		float scaled_corrected_dq_dx=float(corrected_dq_dx)/calib_factor;
		float cal_de_dx=Dedx(scaled_corrected_dq_dx,tot_Ef(trkhitx[i][hitplane][j],trkhity[i][hitplane][j],trkhitz[i][hitplane][j]));
//		float cal_de_dx=-1;
//		if(usemap) float cal_de_dx=Dedx(scaled_corrected_dq_dx,tot_Ef(trkhitx[i][hitplane][j],trkhity[i][hitplane][j],trkhitz[i][hitplane][j]));
//		else float cal_de_dx=Dedx(scaled_corrected_dq_dx,0.4867);
		hdedx->Fill(cal_de_dx);
	       	int bin=int(trkresrange[i][hitplane][j])/binsize;
		fhist_dedx->Fill(trkresrange[i][hitplane][j],cal_de_dx);
		fhist_dqdxcal->Fill(trkresrange[i][hitplane][j],corrected_dq_dx);
		fhist_dqdxuncal->Fill(trkresrange[i][hitplane][j],trkdqdx[i][hitplane][j]);
		dqdx_cal->Fill(corrected_dq_dx/calib_factor);
		dqdx_uncal->Fill(trkdqdx[i][hitplane][j]/calib_factor);
		if(bin<nbin){
		  dedx[bin]->Fill(cal_de_dx);
		} // bin <40
	       }
		else if(hitplane == 2){ //plane 2
		int x_bin=X_correction_hist->FindBin(trkhitx[i][hitplane][j]);
                float Cx=X_correction_hist->GetBinContent(x_bin);
                float Cyzneg=YZ_correction_neg_hist->GetBinContent(YZ_correction_neg_hist->FindBin(trkhitz[i][hitplane][j],trkhity[i][hitplane][j]));
                float corrected_dq_dx=trkdqdx[i][hitplane][j]*Cx*normalisation_factor[hitplane]*Cyzneg;
                float scaled_corrected_dq_dx=float(corrected_dq_dx)/calib_factor;
		float cal_de_dx=Dedx(scaled_corrected_dq_dx,tot_Ef(trkhitx[i][hitplane][j],trkhity[i][hitplane][j],trkhitz[i][hitplane][j]));
//		float cal_de_dx=-1;
//                if(usemap) float cal_de_dx=Dedx(scaled_corrected_dq_dx,tot_Ef(trkhitx[i][hitplane][j],trkhity[i][hitplane][j],trkhitz[i][hitplane][j]));
//                else float cal_de_dx=Dedx(scaled_corrected_dq_dx,0.4867);
                hdedx->Fill(cal_de_dx);
                int bin=int(trkresrange[i][hitplane][j])/binsize;
                fhist_dedx->Fill(trkresrange[i][hitplane][j],cal_de_dx);
                fhist_dqdxcal->Fill(trkresrange[i][hitplane][j],corrected_dq_dx);
                fhist_dqdxuncal->Fill(trkresrange[i][hitplane][j],trkdqdx[i][hitplane][j]);
                dqdx_cal->Fill(corrected_dq_dx/calib_factor);
                dqdx_uncal->Fill(trkdqdx[i][hitplane][j]/calib_factor);
                if(bin<nbin){
                  dedx[bin]->Fill(cal_de_dx);
                } // bin <40
               }
	      }
	      if(trkhitx[i][hitplane][j]>0 && trkhitx[i][hitplane][j]<360){ //positive X direction
	        if(hitplane == 1 && abs(180/3.14*trackthetaxz[i])<40 && !(abs(180/3.14*trackthetayz[i])>80 && abs(180/3.14*trackthetayz[i])<100)){ //plane 1
		int x_bin=X_correction_hist->FindBin(trkhitx[i][hitplane][j]);
                float Cx=X_correction_hist->GetBinContent(x_bin);
                float Cyzpos=YZ_correction_pos_hist->GetBinContent(YZ_correction_pos_hist->FindBin(trkhitz[i][hitplane][j],trkhity[i][hitplane][j]));
                float corrected_dq_dx=trkdqdx[i][hitplane][j]*Cx*normalisation_factor[hitplane]*Cyzpos;
                float scaled_corrected_dq_dx=float(corrected_dq_dx)/calib_factor;
		float cal_de_dx=Dedx(scaled_corrected_dq_dx,tot_Ef(trkhitx[i][hitplane][j],trkhity[i][hitplane][j],trkhitz[i][hitplane][j]));
//		float cal_de_dx=-1;
//                if(usemap) float cal_de_dx=Dedx(scaled_corrected_dq_dx,tot_Ef(trkhitx[i][hitplane][j],trkhity[i][hitplane][j],trkhitz[i][hitplane][j]));
//                else float cal_de_dx=Dedx(scaled_corrected_dq_dx,0.4867);
                hdedx->Fill(cal_de_dx);
                int bin=int(trkresrange[i][hitplane][j])/binsize;
                fhist_dedx->Fill(trkresrange[i][hitplane][j],cal_de_dx);
                fhist_dqdxcal->Fill(trkresrange[i][hitplane][j],corrected_dq_dx);
                fhist_dqdxuncal->Fill(trkresrange[i][hitplane][j],trkdqdx[i][hitplane][j]);
                dqdx_cal->Fill(corrected_dq_dx/calib_factor);
                dqdx_uncal->Fill(trkdqdx[i][hitplane][j]/calib_factor);
                if(bin<nbin){
                  dedx[bin]->Fill(cal_de_dx);
                } // x containment....
                }
	
		else if(hitplane == 0 && abs(180/3.14*trackthetaxz[i])>130 && !(abs(180/3.14*trackthetayz[i])>80 && abs(180/3.14*trackthetayz[i])<100)){ //plane 0
		int x_bin=X_correction_hist->FindBin(trkhitx[i][hitplane][j]);
		float Cx=X_correction_hist->GetBinContent(x_bin);
		float Cyzpos=YZ_correction_pos_hist->GetBinContent(YZ_correction_pos_hist->FindBin(trkhitz[i][hitplane][j],trkhity[i][hitplane][j]));
		float corrected_dq_dx=trkdqdx[i][hitplane][j]*Cx*normalisation_factor[hitplane]*Cyzpos;
		float scaled_corrected_dq_dx=float(corrected_dq_dx)/calib_factor;
		float cal_de_dx=Dedx(scaled_corrected_dq_dx,tot_Ef(trkhitx[i][hitplane][j],trkhity[i][hitplane][j],trkhitz[i][hitplane][j]));
//		float cal_de_dx=-1;
//                if(usemap) float cal_de_dx=Dedx(scaled_corrected_dq_dx,tot_Ef(trkhitx[i][hitplane][j],trkhity[i][hitplane][j],trkhitz[i][hitplane][j]));
//                else float cal_de_dx=Dedx(scaled_corrected_dq_dx,0.4867);
		hdedx->Fill(cal_de_dx);
		int bin=int(trkresrange[i][hitplane][j])/binsize;
		fhist_dedx->Fill(trkresrange[i][hitplane][j],cal_de_dx);
		fhist_dqdxcal->Fill(trkresrange[i][hitplane][j],corrected_dq_dx);
		fhist_dqdxuncal->Fill(trkresrange[i][hitplane][j],trkdqdx[i][hitplane][j]);
		dqdx_cal->Fill(corrected_dq_dx/calib_factor);
		dqdx_uncal->Fill(trkdqdx[i][hitplane][j]/calib_factor);
		if(bin<nbin){
		  dedx[bin]->Fill(cal_de_dx);
		} // x containment....
		}
		else if(hitplane == 2){ //plane 2
		int x_bin=X_correction_hist->FindBin(trkhitx[i][hitplane][j]);
                float Cx=X_correction_hist->GetBinContent(x_bin);
                float Cyzpos=YZ_correction_pos_hist->GetBinContent(YZ_correction_pos_hist->FindBin(trkhitz[i][hitplane][j],trkhity[i][hitplane][j]));
                float corrected_dq_dx=trkdqdx[i][hitplane][j]*Cx*normalisation_factor[hitplane]*Cyzpos;
                float scaled_corrected_dq_dx=float(corrected_dq_dx)/calib_factor;
		float cal_de_dx=Dedx(scaled_corrected_dq_dx,tot_Ef(trkhitx[i][hitplane][j],trkhity[i][hitplane][j],trkhitz[i][hitplane][j]));
//		float cal_de_dx=-1;
//                if(usemap) float cal_de_dx=Dedx(scaled_corrected_dq_dx,tot_Ef(trkhitx[i][hitplane][j],trkhity[i][hitplane][j],trkhitz[i][hitplane][j]));
//                else float cal_de_dx=Dedx(scaled_corrected_dq_dx,0.4867);
                hdedx->Fill(cal_de_dx);
                int bin=int(trkresrange[i][hitplane][j])/binsize;
                fhist_dedx->Fill(trkresrange[i][hitplane][j],cal_de_dx);
                fhist_dqdxcal->Fill(trkresrange[i][hitplane][j],corrected_dq_dx);
                fhist_dqdxuncal->Fill(trkresrange[i][hitplane][j],trkdqdx[i][hitplane][j]);
                dqdx_cal->Fill(corrected_dq_dx/calib_factor);
                dqdx_uncal->Fill(trkdqdx[i][hitplane][j]/calib_factor);
                if(bin<nbin){
                  dedx[bin]->Fill(cal_de_dx);
                } // x containment....
                }
	      } // bin <40		 
	    } // track pitch cut
	  } // z containment.......
	} // y containment.....
      } // loop over hits....
    } // loop over crossing trks.......
  } // loop over jentries...........



  std::cout << "***** Number of muons passing the cut : " << used_trks << std::endl;

  std::cout << "************************** Fitting Landau + Gaussian function to the histograms *****************************" << std::endl;

  ////////////////////////////////////// Fitting Landau+Gaussian function to the histogram ////////////////////////////////

  //ifstream in_stream_1;
  //in_stream_1.open("output_ke_range.txt");
  ofstream myfile1;
  //double KE[13];
  //double Range[13];
  double range_new[nbin];
  double dedxtheory[nbin];
  //double keng;
  //double krange;
  //for(int i=0; i<13; i++){
  //  in_stream_1 >> keng >> krange;
  //  KE[i]=keng;
  //  Range[i]=krange;
  //}
  //in_stream_1.close();

  //TSpline3 *sp = new TSpline3("Cubic Spline",Range,KE,13,"b2e2",0,0);
  TSpline3 *sp = new TSpline3("Cubic Spline", &spline_Range[0], &spline_KE[0],13,"b2e2",0,0);

  vector<double> range;
  vector<double> erange;
  vector<double> energy;
  vector<double> eenergy;
  float range_measured[nbin];
  float energy_measured[nbin];

  vector<double> Thke;
  vector<double> Theng;
  vector<double> Theke;
  vector<double> Theeng;

  vector<double> Thaveng;

  vector<double> chi_denominator;
  vector<double> chi_numerator;
  int dof=0;

  myfile1.open(Form("muon_mpv_0.50_r%d.txt",run));
  for (int i=0; i<nbin; i++){
    std::cout << "Fitting ************** " << i << std::endl;
    Thke.push_back(sp->Eval(i*binsize+double(binsize)/2));
    Theng.push_back(dpdx(sp->Eval(i*binsize+double(binsize)/2),pitchvalue,Mmu));
    myfile1<<i*binsize+double(binsize)/2<<"  "<<dpdx(sp->Eval(i*binsize+double(binsize)/2),pitchvalue,Mmu)<<endl;
    range_new[i]=i*binsize+2.5;
    dedxtheory[i]=dpdx(sp->Eval(i*binsize+double(binsize)/2),pitchvalue,Mmu);
    Theke.push_back(double(sp->Eval((i+1)*binsize)-sp->Eval(i*binsize))/2);
    Theeng.push_back(0);
    Thaveng.push_back(dEdx(sp->Eval(i*binsize+double(binsize)/2),Mmu));
    //d[i]=new TCanvas(Form("d_%d",i),Form("d_%d",i));
    Double_t fr[2];
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    fr[0]=1.1;//this was originally 0.
    fr[1]=10.;
    if(i==0){
      fr[0]=2.2;
      fr[1]=15;
    }
    if (dedx[i]->GetMean()<10){
      //sv[0]=0.13*dedx[i]->GetRMS(); sv[1]=0.8*dedx[i]->GetMean(); sv[2]=dedx[i]->GetEntries()*0.1; sv[3]=.3;
      sv[0]=0.1; sv[1]=1.66; sv[2]=dedx[i]->GetEntries()*0.05; sv[3]=0.05;
      if(i==0){ sv[0]=0.2; sv[1]=4.7; sv[2]=20; sv[3]=.01;}
      if(i==1){ sv[0]=0.2; sv[1]=3.0; sv[2]=10; sv[3]=.01;}
      if(i==2){ sv[1]=2.5;}
      if(i==3){ sv[1]=2.0;}
      if(i==4){ sv[1]=2.0;}

    }
    else{
      sv[0]=0.16*dedx[i]->GetRMS(); sv[1]=0.9*dedx[i]->GetMean(); sv[2]=dedx[i]->GetEntries()*100; sv[3]=dedx[i]->GetRMS()/5.;
    }
    for(int j=0; j<4; ++j){
      pllo[j] = 0.01*sv[j];
      plhi[j] = 100*sv[j];
    }
    Double_t chisqr;
    Int_t    ndf;
    Int_t    status;
//    TF1 *fitsnr = langaufit(dedx[i],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
    TF1 *fitsnr = langaufit(dedx[i],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status);
//    cout <<"************ Fit status (gMinuit): " << gMinuit << ", "<< gMinuit->fCstatu.Data() <<" *********"<<endl;
    cout <<"************ Fit status (FitPtr): " << status << " *********"<<endl;
    fitsnr->SetLineColor(kRed);
    std::cout << "************** MPV : " << fitsnr->GetParameter(1) << " +/- " << fitsnr->GetParError(1) << std::endl;
    std::cout << "************** Chi^2/NDF : " << fitsnr->GetChisquare()/fitsnr->GetNDF() << std::endl;
    std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$  KE : " << sp->Eval((i*binsize+double(binsize)/2)) << "   MPV : " << fitsnr->GetParameter(1) << " $$$$$$$$$$$$" << std::endl;
//    if(status&&dedx[i]->GetEntries()>100){
    if( (dedx[i]->GetEntries()>100)){
//      TString test =  gMinuit->fCstatu.Data(); 
      if(fitsnr->GetNDF() != 0 && status > 2){
//      if(test.EqualTo("CONVERGED ") && (fitsnr->GetParError(1)<1000) && ((fitsnr->GetChisquare()/fitsnr->GetNDF()<10))/* && (fitsnr->GetParError(0)<0.1)*/){
      if((fitsnr->GetParError(1)<1000) && ((fitsnr->GetChisquare()/fitsnr->GetNDF()<10))/* && (fitsnr->GetParError(0)<0.1)*/){
	 range.push_back(sp->Eval((i*binsize+double(binsize)/2)));
	 erange.push_back(double(sp->Eval((i+1)*binsize)-sp->Eval(i*binsize))/2);
	 range_measured[i]=i*binsize+double(binsize)/2;
	 energy_measured[i]=fitsnr->GetParameter(1);
	 cout<<" i "<<i<<" res range "<<i*binsize+double(binsize)/2<<"  KE "<<sp->Eval(i*binsize+double(binsize)/2)<<endl;
	 energy.push_back(fitsnr->GetParameter(1));
	 eenergy.push_back(fitsnr->GetParError(1));
	  
	      
	  ////////////////////////////////////////////// Chi 2 calculation ////////////////////////////
	      
	  if(sp->Eval((i*binsize+double(binsize)/2))<=450 && sp->Eval((i*binsize+double(binsize)/2))>=250){
	    double mpv_err=TMath::Power(fitsnr->GetParError(1),2);
	    //double recom_err=TMath::Power(fitsnr->GetParameter(1)*0.015,2);
	    //double meth_err=TMath::Power(fitsnr->GetParameter(1)*0.01,2);
	    double tot_err=mpv_err;//+recom_err;//+meth_err;
	    chi_denominator.push_back(tot_err);
	    double num=dpdx(sp->Eval(i*binsize+double(binsize)/2),pitchvalue,Mmu)-fitsnr->GetParameter(1);
	    num=TMath::Power(num,2);
	    chi_numerator.push_back(num);
	    dof++;
	  }
	   
	  ///////////////////////////////////////////// End of chi 2 cal. //////////////////////////
	}
      }
    }
  }

  double sum=0;

  for(size_t j=0; j<chi_numerator.size(); j++){
    double ratio=double(chi_numerator[j])/chi_denominator[j];
    cout<<"chi 2 num and den "<<chi_numerator[j]<<"  "<<chi_denominator[j]<<endl;
    sum=sum+ratio;
  }

  std::cout << "$$$$$$$ Chi squared : " << sum << " $$$$$$$$$$$$" << std::endl;
  outfile<<calib_factor<<"\t"<<sum<<std::endl;
  std::cout << "$$$$$$$ Chi squared/NDF : " << double(sum)/dof << " $$$$$$$$$" << std::endl;

  TGraphErrors *th_mpv_graph = new TGraphErrors(Thke.size(),&Thke[0],&Theng[0],&Theke[0],&Theeng[0]);
  th_mpv_graph->SetMarkerStyle(21);
  th_mpv_graph->SetMarkerColor(kBlue);

  TGraphErrors *th_ave_graph = new TGraphErrors(Thke.size(),&Thke[0],&Thaveng[0],&Theke[0],&Theeng[0]);
  th_ave_graph->SetMarkerStyle(21);
  th_ave_graph->SetMarkerColor(kGreen);

  TGraphErrors *e_graph = new TGraphErrors(range.size(),&range[0],&energy[0],&erange[0],&eenergy[0]);
  e_graph->GetXaxis()->SetTitle("Kinetic energy (MeV/cm)");
  e_graph->GetYaxis()->SetTitle("MPV Energy (MeV/cm)");
  e_graph->SetTitle("Kinetic energy Vs MPV energy");
  e_graph->SetMarkerStyle(21);
  e_graph->SetMarkerColor(kRed);

  TText *lable_1 = new TText();
  lable_1-> SetNDC();
  lable_1-> SetTextFont(1);
  lable_1-> SetTextColor(kRed);
  lable_1-> SetTextSize(0.03);
  lable_1-> SetTextAlign(22);
  lable_1-> SetTextAngle(0);

  TCanvas *c1 = new TCanvas("c1","c1");
  TCanvas *cdedx=new TCanvas("cdedx","cdedx");
  TGraph *theoretical = new TGraph(nbin,range_new,dedxtheory);

  c1->SetGrid();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(12);
  c1->cd();
  my_hist->Draw();
  e_graph->Draw("samePE1");
  th_mpv_graph->Draw("samePE1");
  TLegend *legend0=new TLegend(0.2,0.7,0.5,0.8);
  legend0->SetTextFont(72);
  legend0->SetTextSize(0.04);
  legend0->AddEntry(th_mpv_graph,"Theory","lep");
  legend0->AddEntry(e_graph,"Measured","lep");     	 
  legend0->Draw("same");
  lable_1-> DrawText(0.66, 0.58, "Region for chi^2 minimization (250 MeV - 450 MeV)");
  //th_ave_graph->Draw("samePE1");
  c1->Draw();
  c1->Write();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  fhist_dedx->Write();
  dqdx_rat->Write();
  fhist_trkpitch->Write();
  myfile1.close();
  cdedx->cd();
  gStyle->SetPalette(1,0);
  fhist_dedx->Draw("COLZ");
  theoretical->SetMarkerColor(kBlack);
  theoretical->SetMarkerSize(2);
  theoretical->SetMarkerStyle(5);
  theoretical->Draw("sameCP");
  

  TLegend *legdedx=new TLegend(0.2,0.7,0.8,0.8);
  legdedx->SetTextFont(72);
  legdedx->SetTextSize(0.03);
  legdedx->AddEntry(theoretical,"Theoretical most probable value (Landau-Vavilov theory)","LP");
  float X[11]= {0.9833,1.786,3.321,6.598,10.58,30.84,42.50,67.32,106.3,172.5,238.5};
  float Y[11]= {5.687,4.461,3.502,2.731,2.340,1.771,1.670,1.570,1.519,1.510,1.526};
  for(int i=0;i<11;i++){

    X[i]/=1.396;
    Y[i]*=1.396;
  }
  TGraph *gr2 = new TGraph(11,X,Y);
  gr2->SetLineWidth(3);
  gr2->SetLineColor(12);
  gr2->SetMarkerStyle(32);
  legdedx->Draw("same");
  //gr2->Draw("sameCP");
  TGraph *gr3 = new TGraph(nbin,range_measured,energy_measured); 
  gr3->SetLineWidth(3);
  gr3->SetLineColor(kBlue);
  gr3->SetMarkerStyle(33);
  gr3->Draw("sameCP");
  legdedx->AddEntry(gr3,"Measured using Gauss-Landau fit","LP");     	 
  legdedx->Draw("same");
 
  cdedx->Draw();

  cdedx->Write();
  fhist_dqdxcal->Write();
  fhist_dqdxuncal->Write();
  fhist_trkpitch->Draw();
  hdedx->Write();
  dqdx_cal->Write();
  dqdx_uncal->Write();
  gr3->Write("measured_dedxrr");
  theoretical->Write("theoretical");
  file->Write();
  file->Close();
 

  cout<<"Total tracks used"<<used_trks<<endl;

  std::cout << "************************* Calibration.C has ended ***************************" << std::endl;

  gStyle->SetPalette(1,0);
  gStyle->SetNumberContours(64);
  fhist_dqdxcal->Draw("colz");
 

}

int main(int argc, char *argv[]) {

  bool found_input = false,
       found_fcl = false,
       found_out = false;
  //string infile = argv[1];
  string infile;
  std::string fcl_file;

  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     fcl_file = argv[++iArg];
     found_fcl = true;
    }
    if (!strcasecmp(argv[iArg],"-i")) {
     infile = argv[++iArg];
     found_input = true;
    }
    if (!strcasecmp(argv[iArg],"-o")) {
      outfile_name = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << "Usage: dEdX_calibration " <<
                   "-i <input root file OR input list>" <<
                   "-c fclfile.fcl " << 
                   "-o outputfile.txt " << std::endl;
      return 1;
    }
  }

  if (!found_input) {
    cout << "Error: No input file was provided! Please provide with '-i'" << endl;
    return 0;
  }
  if (!found_fcl) {
    cout << "Error: No fcl file was provided! Please provide with '-c'" << endl;
    return 0;
  }
  if (!found_out) {
    cout << "Error: No fcl file was provided! Please provide with '-o'" << endl;
    return 0;
  }

  ////Setting up fcl parameters
  fhicl::ParameterSet pset;
  char const* fhicl_env = getenv("FHICL_FILE_PATH");
  std::string search_path;

  if (!fhicl_env) {
    std::cerr << "Expected environment variable FHICL_FILE_PATH is missing " <<
                 "or empty: using \".\"\n";
    search_path = ".";
  }
  else {
    search_path = std::string{fhicl_env};
  }

  cet::filepath_first_absolute_or_lookup_with_dot lookupPolicy{search_path};
  fhicl::make_ParameterSet(fcl_file, lookupPolicy, pset);
  /////

  //here
  std::string field_file = pset.get<std::string>("FieldMap");
  //ef = new TFile(field_file.c_str(), "OPEN");
  ef = OpenFile(field_file);

  std::vector<std::pair<int, double>> KE_Range
      = pset.get<std::vector<std::pair<int, double>>>("KE_Range");
  //std::vector<int> spline_KE;
  //std::vector<double> spline_Range;
  for (size_t i = 0; i < KE_Range.size(); ++i) {
    spline_KE[i] = KE_Range[i].first;
    spline_Range[i] = KE_Range[i].second;
  }

  usemap=true;

  TChain* shtree = new TChain("Event");
  TChain* shtree1 = new TChain("Event");
  TChain* shtree2 = new TChain("Event");

  if (infile.substr(infile.find_last_of(".") + 1) == "root"){
    shtree->Add(Form("%s/michelremoving2/Event", infile.c_str()));
    shtree1->Add(Form("%s/michelremoving2/Event", infile.c_str()));
    shtree2->Add(Form("%s/michelremoving2/Event", infile.c_str()));
  }
	
  else /*if(infile.substr(infile.find_last_of(".") + 1) == "txt")*/{
    std::ifstream in;
    in.open(infile.c_str());
    char line[1024];

    while(1){
     in.getline(line,1024);
      if (!in.good()) break;
      shtree->Add(Form("%s/michelremoving2/Event", line));
      shtree1->Add(Form("%s/michelremoving2/Event", line));
      shtree2->Add(Form("%s/michelremoving2/Event", line));
    }
    in.close();
    in.clear();
  }
  
  hitplane = 0;
  protoDUNE_dEdx_calib *t=new protoDUNE_dEdx_calib(shtree);
  
  double plane_0_low = pset.get<double>("Plane0Start");
  double plane_0_high = pset.get<double>("Plane0End");
  double plane_0_diff = pset.get<double>("Plane0Diff");
  double plane_1_low = pset.get<double>("Plane1Start");
  double plane_1_high = pset.get<double>("Plane1End");
  double plane_1_diff = pset.get<double>("Plane1Diff");
  double plane_2_low = pset.get<double>("Plane2Start");
  double plane_2_high = pset.get<double>("Plane2End");
  double plane_2_diff = pset.get<double>("Plane2Diff");
  //for(calib_factor = 1.028e-3; calib_factor<1.029e-3; calib_factor+=.001e-3) t->Loop();
  for(calib_factor = plane_0_low; calib_factor<plane_0_high; calib_factor+=plane_0_diff) t->Loop();
  delete t;
  std::cout << "************************* Start of hitplane 1 ***************************" << std::endl;
  
  hitplane = 1;
  protoDUNE_dEdx_calib *t1=new protoDUNE_dEdx_calib(shtree1);
  //for(calib_factor = 1.025e-3;  calib_factor<1.026e-3; calib_factor+=.001e-3) t1->Loop();
  for(calib_factor = plane_1_low; calib_factor<plane_1_high; calib_factor+=plane_1_diff) t1->Loop();
  delete t1;
  std::cout << "************************* Start of hitplane 2 ***************************" << std::endl;

  hitplane = 2;
  protoDUNE_dEdx_calib *t2=new protoDUNE_dEdx_calib(shtree2);
  //for(calib_factor = 1.011e-3; calib_factor<1.012e-3; calib_factor+=.001e-3) t2->Loop();
  for(calib_factor = plane_2_low; calib_factor<plane_2_high; calib_factor+=plane_2_diff) t2->Loop();
  delete t2; 

} // main

TFile * OpenFile(const std::string filename) {
  TFile * theFile = 0x0;
  mf::LogInfo("OpenFile") << "Searching for " << filename;
  if (cet::file_exists(filename)) {
    mf::LogInfo("OpenFile") << "File exists. Opening " << filename;
    theFile = new TFile(filename.c_str());
    if (!theFile ||theFile->IsZombie() || !theFile->IsOpen()) {
      delete theFile;
      theFile = 0x0;
      throw cet::exception("ProtoDUNECalibration.cxx") << "Could not open " << filename;
    }
  }
  else {
    mf::LogInfo("OpenFile") << "File does not exist here. Searching FW_SEARCH_PATH";
    cet::search_path sp{"FW_SEARCH_PATH"};
    std::string found_filename;
    auto found = sp.find_file(filename, found_filename);
    if (!found) {
      throw cet::exception("ProtoDUNECalibration.cxx") << "Could not find " << filename;
    }

    mf::LogInfo("OpenFile") << "Found file " << found_filename;
    theFile = new TFile(found_filename.c_str());
    if (!theFile ||theFile->IsZombie() || !theFile->IsOpen()) {
      delete theFile;
      theFile = 0x0;
      throw cet::exception("ProtoDUNECalibration.cxx") << "Could not open " << found_filename;
    }
  }
  return theFile;
};
