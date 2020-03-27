//README
//This macro is an example of how to loop over the 
//PandoraBeam TTree that was create by the LArsoft model ProtoDUNEelectronAnaTree_module.cc
//It would make a selection sample based on the tech note 
//to run do $make -f Makefile, then execute the binary ./doAna mydatasample.root myoutputfile.root 

#include <iostream>
#include "include.h"
#include "AnalysisTree.h"
using namespace std;

/***************modified box model parameters and function*****************/
double Rho = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia) 
double betap = 0.212;//(kV/cm)(g/cm^2)/MeV
double alpha = 0.93;//parameter from ArgoNeuT experiment at 0.481kV/cm 
double Wion = 23.6e-6;//parameter from ArgoNeuT experiment at 0.481kV/cm
double normalisation_factor=1.0;
double calib_factor = 5.6e-3;
double recombination =0.6417;

double dEdx(float dqdx, float Efield) {  
  return (exp(dqdx*(betap/(Rho*Efield)*Wion))-alpha)/(betap/(Rho*Efield)); 
}

   TFile *ef=new TFile("/dune/data/users/higuera/data/SCE_DataDriven_180kV_v3.root");
   TH3F *xneg=(TH3F*)ef->Get("Reco_ElecField_X_Neg");
   TH3F *yneg=(TH3F*)ef->Get("Reco_ElecField_Y_Neg");
   TH3F *zneg=(TH3F*)ef->Get("Reco_ElecField_Z_Neg");
   TH3F *xpos=(TH3F*)ef->Get("Reco_ElecField_X_Pos");
   TH3F *ypos=(TH3F*)ef->Get("Reco_ElecField_Y_Pos");
   TH3F *zpos=(TH3F*)ef->Get("Reco_ElecField_Z_Pos");

double tot_Ef(double xval,double yval,double zval){
     if(xval>=0){
       double ex=0.5+0.5*xpos->GetBinContent(xpos->FindBin(xval,yval,zval));
       double ey=0.5*ypos->GetBinContent(ypos->FindBin(xval,yval,zval));
       double ez=0.5*zpos->GetBinContent(zpos->FindBin(xval,yval,zval));
       return sqrt(ex*ex+ey*ey+ez*ez);
     }
     else{
       double ex=0.5+0.5*xneg->GetBinContent(xneg->FindBin(xval,yval,zval));
       double ey=0.5*yneg->GetBinContent(yneg->FindBin(xval,yval,zval));
       double ez=0.5*zneg->GetBinContent(zneg->FindBin(xval,yval,zval));
       return sqrt(ex*ex+ey*ey+ez*ez);
     }
}


int EventSelection( const char *filename, const char *outfile ){

   TFile *fyz_corr=new TFile("/dune/data/users/higuera/data/YZcalo_ProtoDUNE.root");
   TH2F *YZ_neg=(TH2F*)fyz_corr->Get("run_YZ_negativeX_2");
   TH2F *YZ_pos=(TH2F*)fyz_corr->Get("run_YZ_positiveX_2"); 

   //X Corr.
   TFile *fx_corr=new TFile("/dune/data/users/higuera/data/Xcalo_ProtoDUNE.root"); 
   TH2F *X_correction= (TH2F*)fx_corr->Get("X_correction_2");

   TFile *f_data = new TFile(filename,"READ");
   TFile *f =new TFile(outfile,"recreate");

   TH1D *h_costheta = new TH1D("h_costheta",";cos#theta",100,-1,1);
   TH1D *h_nhits = new TH1D("h_nhits",";",500,0,2000);
   TH1D *h_dEdx = new TH1D("h_dEdx",";",50,0,10);
   TH1D *h_sh_length = new TH1D("h_sh_length",";",100,0,400);
   TH1D *h_sh_length_afcut = new TH1D("h_sh_length_afcut",";",100,0,400);
   TH1D *h_P = new TH1D("h_P","; Beamline Momentum [GeV/c]",300,0,3);
   TH1D *h_ratio = new TH1D("h_ratio",";(Shower Energy - Beamline Momentum)/Beamline Momentum",100,-1,1);
   TH1D *h_Ecal = new TH1D("h_Ecal",";Shower Energy [Gev]",300,0,3);
 
   //========================================================
   TTree *tree =(TTree*)f_data->Get("PandoraBeam");
   AnalysisTree *signal = new AnalysisTree(tree);
   int n_entries = tree->GetEntries();


   for(long int jentry=0; jentry<n_entries;jentry++) {

      signal->GetEntry(jentry); 
      if( signal->beamtrackMomentum == -999.0 ) continue;
      signal->beamtrackMomentum *= 1000.0;
      //============= electrons =========
      if( signal->cerenkovStatus[1] == 1 ){
        //=======shower ===============================
        if( signal->primaryIsshower ==1 ){
          h_nhits->Fill(signal->primaryNHits);        
          h_sh_length->Fill(signal->primaryLength);
          //select complete showers
          if( signal->primaryNHits < 200 ) continue;

          TVector3 zaxis(signal->beamtrackDir[0],signal->beamtrackDir[1],signal->beamtrackDir[2]);
          TVector3 shdir(signal->primaryStartDirection[0],signal->primaryStartDirection[1],signal->primaryStartDirection[2]);
          double costheta = (shdir.Dot(zaxis))/(shdir.Mag()*zaxis.Mag());
          h_costheta->Fill(costheta);
          //remove potential cosmic-ray background 
          if( costheta < 0.9 ) continue;

          //calculate shower energy 
          double ave_Y=0;
          double ave_Z=0;
          int n_hits_y =0;
          int n_hits_z =0;
          //calculate average Y and average Z for cases we could not find a space point for a given hit
          for( size_t h=0; h<signal->primaryShower_nHits; h++){
             if( signal->primaryShower_hit_Y[h] != -999 ){ 
               ave_Y +=  signal->primaryShower_hit_Y[h];
               n_hits_y ++;
             }
             if(signal->primaryShower_hit_Z[h] != -999 ){
               ave_Z +=  signal->primaryShower_hit_Z[h];
               n_hits_z ++;
             }
          }
          ave_Y /= n_hits_y;
          ave_Z /= n_hits_z;
          double YZcorrection_factor =1.0;
          double Xcorrection_factor=1.0;
          double Ecal =0;
          //loop over hits in the shower
          for( size_t h=0; h<signal->primaryShower_nHits; h++){
             Xcorrection_factor=X_correction->GetBinContent(X_correction->FindBin(5834,signal->primaryShower_hit_X[h]));
             if( signal->primaryShower_hit_Y[h] != -999 && signal->primaryShower_hit_Z[h] != -999 ){
               if(signal->primaryShower_hit_X[h]<0) { YZcorrection_factor=YZ_neg->GetBinContent(YZ_neg->FindBin(5834,signal->primaryShower_hit_Z[h],signal->primaryShower_hit_Y[h])); }
               if(signal->primaryShower_hit_X[h]>=0) { YZcorrection_factor=YZ_pos->GetBinContent(YZ_pos->FindBin(5834,signal->primaryShower_hit_Z[h],signal->primaryShower_hit_Y[h])); }
             }
             else {
               if(signal->primaryShower_hit_X[h]<0) { YZcorrection_factor=YZ_neg->GetBinContent(5834,YZ_neg->FindBin(ave_Z,ave_Y)); }
               if(signal->primaryShower_hit_X[h]>=0) { YZcorrection_factor=YZ_pos->GetBinContent(5834,YZ_pos->FindBin(ave_Z,ave_Y)); }
             }
             double E = (Xcorrection_factor*YZcorrection_factor*signal->primaryShower_hit_q[h]*Wion)/(calib_factor*recombination); 
             E *= normalisation_factor;
             Ecal += E;
          }
          //this part calculates dE/dx at the beginning of the shower 
          int idx=0;
          double dEdx_tmp[signal->primarynCal];
          //tempral vectors to be used as vertex
          //SCE corrections do not modify the shower vertex, therefore we need to use the first (X,Y,Z) of the dQ/dx calorimetry point as vertex 
          vector<double> tmp_x; 
          vector<double> tmp_y;
          vector<double> tmp_z;
          for( int k=0; k<signal->primarynCal-2; ++k){  //skip first two hits due to SCE
             k +=2;
             //skip calorimetry points where calorimetry information is not available
             if( signal->primarydQdx[k] == 0 ) continue;
             if( signal->primary_calX[k] == 0. && signal->primary_calY[k] == 0. && signal->primary_calZ[k] ==0 ) continue; 
             if( signal->primary_calZ[k] < 0 ) continue; 
              
             tmp_x.push_back(signal->primary_calX[k]); 
             tmp_y.push_back(signal->primary_calY[k]); 
             tmp_z.push_back(signal->primary_calZ[k]); 
          }   
          //we found a "vertex" 
          for( int k=0; k<signal->primarynCal-2; ++k){
             k += 2;
             if( signal->primarydQdx[k] == 0 ) continue;
             if( signal->primary_calZ[k] == 0. && signal->primary_calY[k] == 0. && signal->primary_calZ[k] ==0 ) continue; 
             if( signal->primary_cal_pitch[k] > 0 ) 
             double YZcorrection_factor =1.0;
             double Xcorrection_factor=1.0;
             if(signal->primary_calX[k]<0) { YZcorrection_factor=YZ_neg->GetBinContent(YZ_neg->FindBin(5834,signal->primary_calZ[k],signal->primary_calY[k])); }
             if(signal->primary_calX[k]>=0) { YZcorrection_factor=YZ_pos->GetBinContent(YZ_pos->FindBin(5834,signal->primary_calZ[k],signal->primary_calY[k])); }
             if( signal->primary_calZ[k] < 0 ) continue; //no YZcorrection_factor;
             
             Xcorrection_factor=X_correction->GetBinContent(X_correction->FindBin(5834,signal->primary_calX[k]));
             double corrected_dQdx=signal->primarydQdx[k]*normalisation_factor*Xcorrection_factor*YZcorrection_factor;

             double scaled_dQdx=corrected_dQdx/calib_factor;
             double Efield=tot_Ef(signal->primary_calX[k],signal->primary_calY[k],signal->primary_calZ[k]);
             double cali_dEdx=dEdx(scaled_dQdx, Efield); //this is your final dE/dx value
             //calculate dE/dx at the first 4 cm of the shower 
             double newvtx_x = tmp_x.at(0);
             double newvtx_y = tmp_y.at(0);
             double newvtx_z = tmp_z.at(0);
             double phi = TMath::ATan2(signal->primaryStartDirection[1],signal->primaryStartDirection[0]);
             double theta = TMath::ACos(signal->primaryStartDirection[2]);
             double x = 4.0*sin(theta)*cos(phi)+newvtx_x;
             double y = 4.0*sin(theta)*sin(phi)+newvtx_y;
             double z = 4.0*cos(theta)+newvtx_z;
            
             TVector3 x1(newvtx_x,newvtx_y,newvtx_z);
             TVector3 x2(x,y,z);
             TVector3 x0(signal->primary_calX[k],signal->primary_calY[k],signal->primary_calZ[k]);
             TVector3 diff1 = x2-x1;
             TVector3 diff2 = x1-x0;
             TVector3 cross = diff1.Cross(diff2);
             double r= cross.Mag()/diff1.Mag();
             double d = sqrt( pow(signal->primary_calX[k]-newvtx_x,2)+ pow(signal->primary_calY[k]-newvtx_y,2) + pow(signal->primary_calZ[k]-newvtx_z,2));
             //save dE/dx for the first 4 cm of the shower
             if( r < 1.0 && d < 4.0){
                dEdx_tmp[idx]=cali_dEdx;
                idx ++;
             } 


          }
          //calculate the median
          double dEdx_m =0;
          dEdx_m = TMath::Median(idx,dEdx_tmp);
         
          if( dEdx_m == 0 ) continue; 
          h_dEdx->Fill(dEdx_m);	
          h_sh_length_afcut->Fill(signal->primaryLength);
          h_P->Fill(signal->beamtrackMomentum/1000.0); 
          h_Ecal->Fill(Ecal/1000.0);
          h_ratio->Fill((Ecal-signal->beamtrackMomentum)/signal->beamtrackMomentum); 

        }//showers
      }//electrons
   } 
 
  f->Write();
  f->Close();
  return 0;
}

int main( int argc, char *argv[] ){
  cout<<"**************************** "<<endl;
  cout<<"*    WELCOME TO JAMAICA    * "<<endl;
  cout<<"*      HAVE A NICE DAY     * "<<endl;
  cout<<"**************************** "<<endl;

  if( argc == 3 ) return EventSelection(argv[1],argv[2]);
  else {
    cout<<"please provide an input and output file "<<endl;
    return 0;
  }

}
