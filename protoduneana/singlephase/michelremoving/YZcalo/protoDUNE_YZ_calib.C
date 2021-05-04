#define protoDUNE_YZ_calib_cxx
#include "protoDUNE_YZ_calib.h"
#include <TH2.h>
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
#include <fstream>

using namespace std;

void protoDUNE_YZ_calib::Loop(TString mn)
{
  ///////////////// Make any changes to the Y and Z bin sizes here ///////////////
  //int x_bin_size=5;
  int y_bin_size = 5; // nbiny bins in y direction
  int z_bin_size = 5; // nbinz bins in z direction

  ///////////////////////////1//////////////////////////////////////////////////
  int zmax=695;
  int ymax=600;
  //int nbinx=720/x_bin_size;
  int nbiny=ymax/y_bin_size;
  int nbinz=zmax/z_bin_size;
  int filtered_tracks=0;
 
  ////plane_2 histograms 
  TH2F *dqdx_ZvsY_negativeX_hist_2 = new TH2F("dqdx_ZvsY_negativeX_hist_2","plane_2_negativeX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
  TH2F *correction_dqdx_ZvsY_negativeX_hist_2 = new TH2F("correction_dqdx_ZvsY_negativeX_hist_2","plane_2_negativeX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
  
  TH2F *corrected_dqdx_negativeX_hist_2 = new TH2F("corrected_dqdx_negativeX_hist_2","plane_2_negativeX;Z Coordinate (cm);Y Coordinate (cm);",nbinz,0,zmax,nbiny,0,ymax);
  TH2F *dqdx_ZvsY_positiveX_hist_2 = new TH2F("dqdx_ZvsY_positiveX_hist_2","plane_2_positiveX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
 
  TH2F *correction_dqdx_ZvsY_positiveX_hist_2 = new TH2F("correction_dqdx_ZvsY_positiveX_hist_2","plane_2_positiveX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
  TH2F *corrected_dqdx_positiveX_hist_2 = new TH2F("corrected_dqdx_positiveX_hist_2","plane_2_positiveX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
 
  ///////////////////////////////////////////////////////////////////////////////
  //plane_1 histograms
  TH2F *dqdx_ZvsY_negativeX_hist_1 = new TH2F("dqdx_ZvsY_negativeX_hist_1","plane_1_negativeX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
  TH2F *correction_dqdx_ZvsY_negativeX_hist_1 = new TH2F("correction_dqdx_ZvsY_negativeX_hist_1","plane_1_negativeX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
  
  TH2F *corrected_dqdx_negativeX_hist_1 = new TH2F("corrected_dqdx_negativeX_hist_1","plane_1_negativeX;Z Coordinate (cm);Y Coordinate (cm);",nbinz,0,zmax,nbiny,0,ymax);
  TH2F *dqdx_ZvsY_positiveX_hist_1 = new TH2F("dqdx_ZvsY_positiveX_hist_1","plane_1_positiveX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
 
  TH2F *correction_dqdx_ZvsY_positiveX_hist_1 = new TH2F("correction_dqdx_ZvsY_positiveX_hist_1","plane_1_positiveX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
  TH2F *corrected_dqdx_positiveX_hist_1 = new TH2F("corrected_dqdx_positiveX_hist_1","plane_1_positiveX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
 
  ///plane_0 histograms
  TH2F *dqdx_ZvsY_negativeX_hist_0 = new TH2F("dqdx_ZvsY_negativeX_hist_0","plane_0_negativeX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
  TH2F *correction_dqdx_ZvsY_negativeX_hist_0 = new TH2F("correction_dqdx_ZvsY_negativeX_hist_0","plane_0_negativeX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
  
  TH2F *corrected_dqdx_negativeX_hist_0 = new TH2F("corrected_dqdx_negativeX_hist_0","plane_0_negativeX;Z Coordinate (cm);Y Coordinate (cm);",nbinz,0,zmax,nbiny,0,ymax);
  TH2F *dqdx_ZvsY_positiveX_hist_0 = new TH2F("dqdx_ZvsY_positiveX_hist_0","plane_0_positiveX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
 
  TH2F *correction_dqdx_ZvsY_positiveX_hist_0 = new TH2F("correction_dqdx_ZvsY_positiveX_hist_0","plane_0_positiveX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
  TH2F *corrected_dqdx_positiveX_hist_0 = new TH2F("corrected_dqdx_positiveX_hist_0","plane_0_positiveX;Z Coordinate (cm);Y Coordinate (cm)",nbinz,0,zmax,nbiny,0,ymax);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //plane_2
  vector<vector<vector<float>>> dqdx_value_negativeX_2;
  vector<float> all_dqdx_value_negativeX_2;
  vector<vector<vector<float>>> dqdx_frac_correction_negativeX_2;
  dqdx_value_negativeX_2.resize(nbinz);
  dqdx_frac_correction_negativeX_2.resize(nbinz);
  for(size_t i = 0; i < dqdx_value_negativeX_2.size(); i++){
    dqdx_value_negativeX_2[i].resize(nbiny);
    dqdx_frac_correction_negativeX_2[i].resize(nbiny);
  }

  vector<vector<vector<float>>> dqdx_value_positiveX_2;
  vector<float> all_dqdx_value_positiveX_2;
  vector<vector<vector<float>>> dqdx_frac_correction_positiveX_2;
  dqdx_value_positiveX_2.resize(nbinz);
  dqdx_frac_correction_positiveX_2.resize(nbinz);
  for(size_t i = 0; i < dqdx_value_positiveX_2.size(); i++){
    dqdx_value_positiveX_2[i].resize(nbiny);
    dqdx_frac_correction_positiveX_2[i].resize(nbiny);
  }
  ////////////////////////////////////////////////////
  //plane_1
  vector<vector<vector<float>>> dqdx_value_negativeX_1;
  vector<float> all_dqdx_value_negativeX_1;
  vector<vector<vector<float>>> dqdx_frac_correction_negativeX_1;
  dqdx_value_negativeX_1.resize(nbinz);
  dqdx_frac_correction_negativeX_1.resize(nbinz);
  for(size_t i = 0; i < dqdx_value_negativeX_1.size(); i++){
    dqdx_value_negativeX_1[i].resize(nbiny);
    dqdx_frac_correction_negativeX_1[i].resize(nbiny);
  }

  vector<vector<vector<float>>> dqdx_value_positiveX_1;
  vector<float> all_dqdx_value_positiveX_1;
  vector<vector<vector<float>>> dqdx_frac_correction_positiveX_1;
  dqdx_value_positiveX_1.resize(nbinz);
  dqdx_frac_correction_positiveX_1.resize(nbinz);
  for(size_t i = 0; i < dqdx_value_positiveX_1.size(); i++){
    dqdx_value_positiveX_1[i].resize(nbiny);
    dqdx_frac_correction_positiveX_1[i].resize(nbiny);
  }
  /////////////////////////////////////////////////////
  //plane_0
  vector<vector<vector<float>>> dqdx_value_negativeX_0;
  vector<float> all_dqdx_value_negativeX_0;
  vector<vector<vector<float>>> dqdx_frac_correction_negativeX_0;
  dqdx_value_negativeX_0.resize(nbinz);
  dqdx_frac_correction_negativeX_0.resize(nbinz);
  for(size_t i = 0; i < dqdx_value_negativeX_0.size(); i++){
    dqdx_value_negativeX_0[i].resize(nbiny);
    dqdx_frac_correction_negativeX_0[i].resize(nbiny);
  }

  vector<vector<vector<float>>> dqdx_value_positiveX_0;
  vector<float> all_dqdx_value_positiveX_0;
  vector<vector<vector<float>>> dqdx_frac_correction_positiveX_0;
  dqdx_value_positiveX_0.resize(nbinz);
  dqdx_frac_correction_positiveX_0.resize(nbinz);
  for(size_t i = 0; i < dqdx_value_positiveX_0.size(); i++){
    dqdx_value_positiveX_0[i].resize(nbiny);
    dqdx_frac_correction_positiveX_0[i].resize(nbiny);
  }
  //////////////////////////////////////////////////////////

 
  //int used_tracks_negativeX=0;
  //int used_tracks_positiveX=0;
  if (fChain == 0) return;

  fChain->GetEntry(0);
  std::cout<<"Process Run "<<run<<std::endl;
  TFile *file = new TFile(Form("YZcalo_mich%s_r%d.root",mn.Data(), run), "recreate");

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry%100==0) cout<<jentry<<"/"<<nentries<<endl;
    for(int i=0; i<cross_trks; i++){
      if(!((TMath::Abs(trkstartx[i])>350||trkstarty[i]<50||trkstarty[i]>550||trkstartz[i]<50||trkstartz[i]>645)&&(TMath::Abs(trkendx[i])>350||trkendy[i]<50||trkendy[i]>550||trkendz[i]<50||trkendz[i]>645))) continue;
      filtered_tracks++;
      ///filling histograms for plane_2
      if(!((abs(180/3.14*trackthetaxz[i])>60 && abs(180/3.14*trackthetaxz[i])<120)||abs(180/3.14*trackthetaxz[i])<10||(abs(180/3.14*trackthetayz[i])>80 && abs(180/3.14*trackthetayz[i])<100))){
	for(int j=0; j<TMath::Min(ntrkhits[i][2],3000); ++j){
	  if((trkhitx[i][2][j]<0)&&(trkhitx[i][2][j]>-360)){
	    if((trkhity[i][2][j]<ymax)&&(trkhity[i][2][j]>0)){
	      if((trkhitz[i][2][j]<zmax)&&(trkhitz[i][2][j]>0)){
		int z_bin = int(trkhitz[i][2][j])/z_bin_size; 
		int y_bin= int(trkhity[i][2][j])/y_bin_size;
		dqdx_value_negativeX_2[z_bin][y_bin].push_back(trkdqdx[i][2][j]);
	      } // Z containment
	    } // Y containment
	  } // X containiment
	} // loop over hits of the track in the given plane
      } // theta xz and yz angle cut
      
      if(!((abs(180/3.14*trackthetaxz[i])>60 && abs(180/3.14*trackthetaxz[i])<120)||abs(180/3.14*trackthetaxz[i])<10||(abs(180/3.14*trackthetayz[i])>80 && abs(180/3.14*trackthetayz[i])<100))){
	for(int j=0; j<TMath::Min(ntrkhits[i][2],3000); ++j){
	  if((trkhitx[i][2][j]>0)&&(trkhitx[i][2][j]<360)){
	    if((trkhity[i][2][j]<ymax)&&(trkhity[i][2][j]>0)){
	      if((trkhitz[i][2][j]<zmax)&&(trkhitz[i][2][j]>0)){
		int z_bin = int(trkhitz[i][2][j])/z_bin_size; 
		int y_bin= int(trkhity[i][2][j])/y_bin_size;
		dqdx_value_positiveX_2[z_bin][y_bin].push_back(trkdqdx[i][2][j]);
	      } // Z containment
	    } // Y containment
	  } // X containiment
	} // loop over hits of the track in the given plane
      } // theta xz and yz angle cut

      ////for plane_1
      if(abs(180/3.14*trackthetaxz[i])>130 && !(abs(180/3.14*trackthetayz[i])>80 && abs(180/3.14*trackthetayz[i])<100)){   
	for(int j=0; j<TMath::Min(ntrkhits[i][1],3000); ++j){
	  if((trkhitx[i][1][j]<0)&&(trkhitx[i][1][j]>-360)){
	    if((trkhity[i][1][j]<ymax)&&(trkhity[i][1][j]>0)){
	      if((trkhitz[i][1][j]<zmax)&&(trkhitz[i][1][j]>0)){
		int z_bin = int(trkhitz[i][1][j])/z_bin_size; 
		int y_bin= int(trkhity[i][1][j])/y_bin_size;
		dqdx_value_negativeX_1[z_bin][y_bin].push_back(trkdqdx[i][1][j]);
	      } // Z containment
	    } // Y containment
	  } // X containiment
	} // loop over hits of the track in the given plane
      } // theta xz and yz angle cut
      if(abs(180/3.14*trackthetaxz[i])<40 && !(abs(180/3.14*trackthetayz[i])>80 && abs(180/3.14*trackthetayz[i])<100)){ 
	for(int j=0; j<TMath::Min(ntrkhits[i][1],3000); ++j){
	  if((trkhitx[i][1][j]>0)&&(trkhitx[i][1][j]<360)){
	    if((trkhity[i][1][j]<ymax)&&(trkhity[i][1][j]>0)){
	      if((trkhitz[i][1][j]<zmax)&&(trkhitz[i][1][j]>0)){
		int z_bin = int(trkhitz[i][1][j])/z_bin_size; 
		int y_bin= int(trkhity[i][1][j])/y_bin_size;
		dqdx_value_positiveX_1[z_bin][y_bin].push_back(trkdqdx[i][1][j]);
	      } // Z containment
	    } // Y containment
	  } // X containiment
	} // loop over hits of the track in the given plane
      } // theta xz and yz angle cut
    
      ///filling histograms for plane_0
      if(abs(180/3.14*trackthetaxz[i])<40 && !(abs(180/3.14*trackthetayz[i])>80 && abs(180/3.14*trackthetayz[i])<100)){
	  for(int j=0; j<TMath::Min(ntrkhits[i][0],3000); ++j){
	    if((trkhitx[i][0][j]<0)&&(trkhitx[i][0][j]>-360)){
	      if((trkhity[i][0][j]<ymax)&&(trkhity[i][0][j]>0)){
		if((trkhitz[i][0][j]<zmax)&&(trkhitz[i][0][j]>0)){
		  int z_bin = int(trkhitz[i][0][j])/z_bin_size; 
		  int y_bin= int(trkhity[i][0][j])/y_bin_size;
		  dqdx_value_negativeX_0[z_bin][y_bin].push_back(trkdqdx[i][0][j]);
		} // Z containment
	      } // Y containment
	    } // X containiment
	  } // loop over hits of the track in the given plane
	} // theta xz and yz angle cut
      if(abs(180/3.14*trackthetaxz[i])>130 && !(abs(180/3.14*trackthetayz[i])>80 && abs(180/3.14*trackthetayz[i])<100)){
	  for(int j=0; j<TMath::Min(ntrkhits[i][0],3000); ++j){
	    if((trkhitx[i][0][j]>0)&&(trkhitx[i][0][j]<360)){
	      if((trkhity[i][0][j]<ymax)&&(trkhity[i][0][j]>0)){
		if((trkhitz[i][0][j]<zmax)&&(trkhitz[i][0][j]>0)){
		  int z_bin = int(trkhitz[i][0][j])/z_bin_size; 
		  int y_bin= int(trkhity[i][0][j])/y_bin_size;
		  dqdx_value_positiveX_0[z_bin][y_bin].push_back(trkdqdx[i][0][j]);
		} // Z containment
	      }// Y containment
	    } // X containiment
	  } // loop over hits of the track in the given plane
      }// theta xz and yz angle cut

    }// loop over crossing tracks in the event
    } // loop over jentries

    std::cout << "*************** Calculating the local median dQ/dx values for each Y-Z cell ******************" << std::endl;
 
    for(size_t i=0; i<dqdx_value_negativeX_2.size(); i++){
      for(size_t j=0; j<dqdx_value_negativeX_2[i].size(); j++){
	if(dqdx_value_negativeX_2[i][j].size()>5){
	  for(size_t k=0; k<dqdx_value_negativeX_2[i][j].size(); k++){
	    all_dqdx_value_negativeX_2.push_back(dqdx_value_negativeX_2[i][j][k]); 
	  }
       
	  float local_median_dqdx_negativeX_2=TMath::Median(dqdx_value_negativeX_2[i][j].size(),&dqdx_value_negativeX_2[i][j][0]);
	  dqdx_ZvsY_negativeX_hist_2->SetBinContent(i+1,j+1,local_median_dqdx_negativeX_2);
	}
      }
    }
    float global_median_dqdx_negativeX_2=TMath::Median(all_dqdx_value_negativeX_2.size(),&all_dqdx_value_negativeX_2[0]);
    for(size_t i=0; i<dqdx_value_positiveX_2.size(); i++){
      for(size_t j=0; j<dqdx_value_positiveX_2[i].size(); j++){
	if(dqdx_value_positiveX_2[i][j].size()>5){
	  for(size_t k=0; k<dqdx_value_positiveX_2[i][j].size(); k++){
	    all_dqdx_value_positiveX_2.push_back(dqdx_value_positiveX_2[i][j][k]); 
	  }
	  float local_median_dqdx_positiveX_2=TMath::Median(dqdx_value_positiveX_2[i][j].size(),&dqdx_value_positiveX_2[i][j][0]);
	  dqdx_ZvsY_positiveX_hist_2->SetBinContent(i+1,j+1,local_median_dqdx_positiveX_2);
	}
      }
    }
    float global_median_dqdx_positiveX_2=TMath::Median(all_dqdx_value_positiveX_2.size(),&all_dqdx_value_positiveX_2[0]);

    /////////////////////////////////////////////////////////////////////////////////////
    std::cout << "********************** Calculating fractional dQ/dx corrections for each Y-Z cell ********************" << std::endl;
    //////////////////// Calculating the fractional corrections in each YZ cell /////////////
    for(size_t i=0; i<dqdx_value_negativeX_2.size(); i++){
      for(size_t j=0; j<dqdx_value_negativeX_2[i].size(); j++){
	if(dqdx_value_negativeX_2[i][j].size()>5){
	  float local_median_dqdx_negativeX_2=TMath::Median(dqdx_value_negativeX_2[i][j].size(),&dqdx_value_negativeX_2[i][j][0]);
	  float fractional_dqdx_negativeX_2=float(global_median_dqdx_negativeX_2)/local_median_dqdx_negativeX_2; 
	  correction_dqdx_ZvsY_negativeX_hist_2->SetBinContent(i+1,j+1,fractional_dqdx_negativeX_2);
	  dqdx_frac_correction_negativeX_2[i][j].push_back(fractional_dqdx_negativeX_2);
	}
      }
    }
    for(size_t i=0; i<dqdx_value_positiveX_2.size(); i++){
      for(size_t j=0; j<dqdx_value_positiveX_2[i].size(); j++){
	if(dqdx_value_positiveX_2[i][j].size()>5){
	  float local_median_dqdx_positiveX_2=TMath::Median(dqdx_value_positiveX_2[i][j].size(),&dqdx_value_positiveX_2[i][j][0]);
	  float fractional_dqdx_positiveX_2=float(global_median_dqdx_positiveX_2)/local_median_dqdx_positiveX_2; 
	  correction_dqdx_ZvsY_positiveX_hist_2->SetBinContent(i+1,j+1,fractional_dqdx_positiveX_2);
	  dqdx_frac_correction_positiveX_2[i][j].push_back(fractional_dqdx_positiveX_2);
	}
      }
    }
    ///////////////////////////////////////////////////////////////////////////////////
    std::cout << "******************** Calculating corrected dQ/dx value for each Y-Z cell **********************" << std::endl;
    //////////////// How corrected YZ dqdx distribution looks like ///////////////////
    for(size_t i=0; i<dqdx_value_negativeX_2.size(); i++){
      for(size_t j=0; j<dqdx_value_negativeX_2[i].size(); j++){
	if(dqdx_value_negativeX_2[i][j].size()>5){
	  float local_median_dqdx_negativeX_2=TMath::Median(dqdx_value_negativeX_2[i][j].size(),&dqdx_value_negativeX_2[i][j][0]);
	  float corrected_dqdx_negativeX_2=local_median_dqdx_negativeX_2*dqdx_frac_correction_negativeX_2[i][j][0];
	  corrected_dqdx_negativeX_hist_2->SetBinContent(i+1,j+1,corrected_dqdx_negativeX_2);
	} 
      }
    }
    for(size_t i=0; i<dqdx_value_positiveX_2.size(); i++){
      for(size_t j=0; j<dqdx_value_positiveX_2[i].size(); j++){
	if(dqdx_value_positiveX_2[i][j].size()>5){
	  float local_median_dqdx_positiveX_2=TMath::Median(dqdx_value_positiveX_2[i][j].size(),&dqdx_value_positiveX_2[i][j][0]);
	  float corrected_dqdx_positiveX_2=local_median_dqdx_positiveX_2*dqdx_frac_correction_positiveX_2[i][j][0];
	  corrected_dqdx_positiveX_hist_2->SetBinContent(i+1,j+1,corrected_dqdx_positiveX_2);
	} 
      }
    }
    /////////////////////////////////////////////////////////////////////////////////


    ////similar calculation for plane_1
    for(size_t i=0; i<dqdx_value_negativeX_1.size(); i++){
      for(size_t j=0; j<dqdx_value_negativeX_1[i].size(); j++){
	if(dqdx_value_negativeX_1[i][j].size()>5){
	  for(size_t k=0; k<dqdx_value_negativeX_1[i][j].size(); k++){
	    all_dqdx_value_negativeX_1.push_back(dqdx_value_negativeX_1[i][j][k]); 
	  }
       
	  float local_median_dqdx_negativeX_1=TMath::Median(dqdx_value_negativeX_1[i][j].size(),&dqdx_value_negativeX_1[i][j][0]);
	  dqdx_ZvsY_negativeX_hist_1->SetBinContent(i+1,j+1,local_median_dqdx_negativeX_1);
	}
      }
    }
    float global_median_dqdx_negativeX_1=TMath::Median(all_dqdx_value_negativeX_1.size(),&all_dqdx_value_negativeX_1[0]);
    for(size_t i=0; i<dqdx_value_positiveX_1.size(); i++){
      for(size_t j=0; j<dqdx_value_positiveX_1[i].size(); j++){
	if(dqdx_value_positiveX_1[i][j].size()>5){
	  for(size_t k=0; k<dqdx_value_positiveX_1[i][j].size(); k++){
	    all_dqdx_value_positiveX_1.push_back(dqdx_value_positiveX_1[i][j][k]); 
	  }
	  float local_median_dqdx_positiveX_1=TMath::Median(dqdx_value_positiveX_1[i][j].size(),&dqdx_value_positiveX_1[i][j][0]);
	  dqdx_ZvsY_positiveX_hist_1->SetBinContent(i+1,j+1,local_median_dqdx_positiveX_1);
	}
      }
    }
    float global_median_dqdx_positiveX_1=TMath::Median(all_dqdx_value_positiveX_1.size(),&all_dqdx_value_positiveX_1[0]);

    /////////////////////////////////////////////////////////////////////////////////////
    std::cout << "********************** Calculating fractional dQ/dx corrections for each Y-Z cell ********************" << std::endl;
    //////////////////// Calculating the fractional corrections in each YZ cell /////////////
    for(size_t i=0; i<dqdx_value_negativeX_1.size(); i++){
      for(size_t j=0; j<dqdx_value_negativeX_1[i].size(); j++){
	if(dqdx_value_negativeX_1[i][j].size()>5){
	  float local_median_dqdx_negativeX_1=TMath::Median(dqdx_value_negativeX_1[i][j].size(),&dqdx_value_negativeX_1[i][j][0]);
	  float fractional_dqdx_negativeX_1=float(global_median_dqdx_negativeX_1)/local_median_dqdx_negativeX_1; 
	  correction_dqdx_ZvsY_negativeX_hist_1->SetBinContent(i+1,j+1,fractional_dqdx_negativeX_1);
	  dqdx_frac_correction_negativeX_1[i][j].push_back(fractional_dqdx_negativeX_1);
	}
      }
    }
    for(size_t i=0; i<dqdx_value_positiveX_1.size(); i++){
      for(size_t j=0; j<dqdx_value_positiveX_1[i].size(); j++){
	if(dqdx_value_positiveX_1[i][j].size()>5){
	  float local_median_dqdx_positiveX_1=TMath::Median(dqdx_value_positiveX_1[i][j].size(),&dqdx_value_positiveX_1[i][j][0]);
	  float fractional_dqdx_positiveX_1=float(global_median_dqdx_positiveX_1)/local_median_dqdx_positiveX_1; 
	  correction_dqdx_ZvsY_positiveX_hist_1->SetBinContent(i+1,j+1,fractional_dqdx_positiveX_1);
	  dqdx_frac_correction_positiveX_1[i][j].push_back(fractional_dqdx_positiveX_1);
	}
      }
    }
    ///////////////////////////////////////////////////////////////////////////////////
    std::cout << "******************** Calculating corrected dQ/dx value for each Y-Z cell **********************" << std::endl;
    //////////////// How corrected YZ dqdx distribution looks like ///////////////////
    for(size_t i=0; i<dqdx_value_negativeX_1.size(); i++){
      for(size_t j=0; j<dqdx_value_negativeX_1[i].size(); j++){
	if(dqdx_value_negativeX_1[i][j].size()>5){
	  float local_median_dqdx_negativeX_1=TMath::Median(dqdx_value_negativeX_1[i][j].size(),&dqdx_value_negativeX_1[i][j][0]);
	  float corrected_dqdx_negativeX_1=local_median_dqdx_negativeX_1*dqdx_frac_correction_negativeX_1[i][j][0];
	  corrected_dqdx_negativeX_hist_1->SetBinContent(i+1,j+1,corrected_dqdx_negativeX_1);
	} 
      }
    }
    for(size_t i=0; i<dqdx_value_positiveX_1.size(); i++){
      for(size_t j=0; j<dqdx_value_positiveX_1[i].size(); j++){
	if(dqdx_value_positiveX_1[i][j].size()>5){
	  float local_median_dqdx_positiveX_1=TMath::Median(dqdx_value_positiveX_1[i][j].size(),&dqdx_value_positiveX_1[i][j][0]);
	  float corrected_dqdx_positiveX_1=local_median_dqdx_positiveX_1*dqdx_frac_correction_positiveX_1[i][j][0];
	  corrected_dqdx_positiveX_hist_1->SetBinContent(i+1,j+1,corrected_dqdx_positiveX_1);
	} 
      }
    }
    ///////////////////////////////////////////////////////////////////////////////
    ////similar calculations for plane 0
    for(size_t i=0; i<dqdx_value_negativeX_0.size(); i++){
      for(size_t j=0; j<dqdx_value_negativeX_0[i].size(); j++){
	if(dqdx_value_negativeX_0[i][j].size()>5){
	  for(size_t k=0; k<dqdx_value_negativeX_0[i][j].size(); k++){
	    all_dqdx_value_negativeX_0.push_back(dqdx_value_negativeX_0[i][j][k]); 
	  }
       
	  float local_median_dqdx_negativeX_0=TMath::Median(dqdx_value_negativeX_0[i][j].size(),&dqdx_value_negativeX_0[i][j][0]);
	  dqdx_ZvsY_negativeX_hist_0->SetBinContent(i+1,j+1,local_median_dqdx_negativeX_0);
	}
      }
    }
    float global_median_dqdx_negativeX_0=TMath::Median(all_dqdx_value_negativeX_0.size(),&all_dqdx_value_negativeX_0[0]);
    for(size_t i=0; i<dqdx_value_positiveX_0.size(); i++){
      for(size_t j=0; j<dqdx_value_positiveX_0[i].size(); j++){
	if(dqdx_value_positiveX_0[i][j].size()>5){
	  for(size_t k=0; k<dqdx_value_positiveX_0[i][j].size(); k++){
	    all_dqdx_value_positiveX_0.push_back(dqdx_value_positiveX_0[i][j][k]); 
	  }
	  float local_median_dqdx_positiveX_0=TMath::Median(dqdx_value_positiveX_0[i][j].size(),&dqdx_value_positiveX_0[i][j][0]);
	  dqdx_ZvsY_positiveX_hist_0->SetBinContent(i+1,j+1,local_median_dqdx_positiveX_0);
	}
      }
    }
    float global_median_dqdx_positiveX_0=TMath::Median(all_dqdx_value_positiveX_0.size(),&all_dqdx_value_positiveX_0[0]);

    /////////////////////////////////////////////////////////////////////////////////////
    std::cout << "********************** Calculating fractional dQ/dx corrections for each Y-Z cell ********************" << std::endl;
    //////////////////// Calculating the fractional corrections in each YZ cell /////////////
    for(size_t i=0; i<dqdx_value_negativeX_0.size(); i++){
      for(size_t j=0; j<dqdx_value_negativeX_0[i].size(); j++){
	if(dqdx_value_negativeX_0[i][j].size()>5){
	  float local_median_dqdx_negativeX_0=TMath::Median(dqdx_value_negativeX_0[i][j].size(),&dqdx_value_negativeX_0[i][j][0]);
	  float fractional_dqdx_negativeX_0=float(global_median_dqdx_negativeX_0)/local_median_dqdx_negativeX_0; 
	  correction_dqdx_ZvsY_negativeX_hist_0->SetBinContent(i+1,j+1,fractional_dqdx_negativeX_0);
	  dqdx_frac_correction_negativeX_0[i][j].push_back(fractional_dqdx_negativeX_0);
	}
      }
    }
    for(size_t i=0; i<dqdx_value_positiveX_0.size(); i++){
      for(size_t j=0; j<dqdx_value_positiveX_0[i].size(); j++){
	if(dqdx_value_positiveX_0[i][j].size()>5){
	  float local_median_dqdx_positiveX_0=TMath::Median(dqdx_value_positiveX_0[i][j].size(),&dqdx_value_positiveX_0[i][j][0]);
	  float fractional_dqdx_positiveX_0=float(global_median_dqdx_positiveX_0)/local_median_dqdx_positiveX_0; 
	  correction_dqdx_ZvsY_positiveX_hist_0->SetBinContent(i+1,j+1,fractional_dqdx_positiveX_0);
	  dqdx_frac_correction_positiveX_0[i][j].push_back(fractional_dqdx_positiveX_0);
	}
      }
    }
    ///////////////////////////////////////////////////////////////////////////////////
    std::cout << "******************** Calculating corrected dQ/dx value for each Y-Z cell **********************" << std::endl;
    //////////////// How corrected YZ dqdx distribution looks like ///////////////////
    for(size_t i=0; i<dqdx_value_negativeX_0.size(); i++){
      for(size_t j=0; j<dqdx_value_negativeX_0[i].size(); j++){
	if(dqdx_value_negativeX_0[i][j].size()>05){
	  float local_median_dqdx_negativeX_0=TMath::Median(dqdx_value_negativeX_0[i][j].size(),&dqdx_value_negativeX_0[i][j][0]);
	  float corrected_dqdx_negativeX_0=local_median_dqdx_negativeX_0*dqdx_frac_correction_negativeX_0[i][j][0];
	  corrected_dqdx_negativeX_hist_0->SetBinContent(i+1,j+1,corrected_dqdx_negativeX_0);
	} 
      }
    }
    for(size_t i=0; i<dqdx_value_positiveX_0.size(); i++){
      for(size_t j=0; j<dqdx_value_positiveX_0[i].size(); j++){
	if(dqdx_value_positiveX_0[i][j].size()>5){
	  float local_median_dqdx_positiveX_0=TMath::Median(dqdx_value_positiveX_0[i][j].size(),&dqdx_value_positiveX_0[i][j][0]);
	  float corrected_dqdx_positiveX_0=local_median_dqdx_positiveX_0*dqdx_frac_correction_positiveX_0[i][j][0];
	  corrected_dqdx_positiveX_hist_0->SetBinContent(i+1,j+1,corrected_dqdx_positiveX_0);
	} 
      }
    }
    /////plane_0 done

    dqdx_ZvsY_negativeX_hist_2->Write();
    dqdx_ZvsY_positiveX_hist_2->Write();
    correction_dqdx_ZvsY_negativeX_hist_2->Write();
    correction_dqdx_ZvsY_positiveX_hist_2->Write();
    corrected_dqdx_negativeX_hist_2->Write();
    corrected_dqdx_positiveX_hist_2->Write();

    dqdx_ZvsY_negativeX_hist_1->Write();
    dqdx_ZvsY_positiveX_hist_1->Write();
    correction_dqdx_ZvsY_negativeX_hist_1->Write();
    correction_dqdx_ZvsY_positiveX_hist_1->Write();
    corrected_dqdx_negativeX_hist_1->Write();
    corrected_dqdx_positiveX_hist_1->Write();

    dqdx_ZvsY_negativeX_hist_0->Write();
    dqdx_ZvsY_positiveX_hist_0->Write();
    correction_dqdx_ZvsY_negativeX_hist_0->Write();
    correction_dqdx_ZvsY_positiveX_hist_0->Write();
    corrected_dqdx_negativeX_hist_0->Write();
    corrected_dqdx_positiveX_hist_0->Write();
    dqdx_ZvsY_negativeX_hist_2->Draw("colz");
    std::cout<<"crossing tracks "<<filtered_tracks<<std::endl;
    file->Close();
    std::cout << "*************** Y_Z_Correction_make_class.C macro has ended ******************" << std::endl;  
  }

int main(int argc, char *argv[]) {
  
  if (!argv[2]) {
    cout << "Error: No input file or michelremoving tree number was provided!" << endl;
    cout << "Usage: " << endl;
    cout << "make_yz_correction root_file_or_list [michelremoving_tree_number]" << endl; 
    return 0;
  }
  
  string infile = argv[1];
  string michelnumber = argv[2];
  cout <<michelnumber <<endl;

  if (!(michelnumber == "0"||michelnumber == "1"||michelnumber == "2")){
    cout << "Error: Michel tree number must be 0,1, or 2" << endl;
    return 0;
    }

  if (michelnumber=="0") michelnumber = "";
  cout << Form("michelremoving%s/Event", michelnumber.c_str()) << endl;

 // string runident = argv[3];
  TChain* shtree = new TChain("Event");

  if (infile.substr(infile.find_last_of(".") + 1) == "root"){
    shtree->Add(Form("%s/michelremoving%s/Event", infile.c_str(), michelnumber.c_str()));
  }
  else if(infile.substr(infile.find_last_of(".") + 1) == "list"){
    std::ifstream in;
    in.open(infile.c_str());
    char line[1024];

    while(1){
      in.getline(line,1024);
      if (!in.good()) break;
      shtree->Add(Form("%s/michelremoving%s/Event", line, michelnumber.c_str()));
    }
    in.close();
    in.clear();
  }

  protoDUNE_YZ_calib t(shtree);

 
  t.Loop(michelnumber.c_str());
} // main
