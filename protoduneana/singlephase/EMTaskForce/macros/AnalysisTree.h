//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 25 21:55:36 2020 by ROOT version 6.18/04
// from TTree PandoraBeam/Beam events reconstructed with Pandora
// found on file: /dune/data/users/higuera/data/r5834_ana.root
//////////////////////////////////////////////////////////

#ifndef AnalysisTree_h
#define AnalysisTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;
class AnalysisTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   Double_t        timestamp;
   Int_t           Nactivefembs[5];
   Int_t           beamtrigger;
   Int_t           beamCheckIsMatched;
   Double_t        tof;
   Int_t           cerenkovStatus[2];
   Double_t        cerenkovTime[2];
   Double_t        cerenkovPressure[2];
   Double_t        beamtrackMomentum;
   Double_t        beamtrackEnergy;
   Double_t        beamtrackPos[3];
   Double_t        beamtrackEndPos[3];
   Double_t        beamtrackDir[3];
   Double_t        beamtrackTime;
   Int_t           beamtrackPdg;
   Int_t           beamtrackID;
   Int_t           beam_ntrjPoints;
   Int_t           beamtrackNDaughters;
   Double_t        beamtrackPos_at[1][3];   //[beam_ntrjPoints]
   Double_t        beamtrackMom_at[1][4];   //[beam_ntrjPoints]
   vector<string>  *beamtrackEndProcess;
   Double_t        primaryVertex[3];
   Int_t           primaryIsBeamparticle;
   Int_t           primaryIstrack;
   Int_t           primaryIsshower;
   Double_t        primaryBDTScore;
   Int_t           primaryNHits;
   Double_t        primaryTheta;
   Double_t        primaryPhi;
   Double_t        primaryLength;
   Double_t        primaryMomentum;
   Double_t        primaryEndMomentum;
   Double_t        primaryEndPosition[3];
   Double_t        primaryStartPosition[3];
   Double_t        primaryEndDirection[3];
   Double_t        primaryStartDirection[3];
   Double_t        primaryOpeningAngle;
   Int_t           primaryID;
   Double_t        primaryTruth_E;
   Double_t        primaryTruth_vtx[3];
   Int_t           primaryTruth_pdg;
   Int_t           primaryTruth_trkID;
   Int_t           primaryShowerBestPlane;
   Double_t        primaryShowerCharge;
   Double_t        primaryShowerEnergy;
   Double_t        primaryShowerMIPEnergy;
   Int_t           primaryShower_nHits;
   Double_t        primaryShower_hit_q[4847];   //[primaryShower_nHits]
   Int_t           primaryShower_hit_w[4847];   //[primaryShower_nHits]
   Double_t        primaryShower_hit_t[4847];   //[primaryShower_nHits]
   Double_t        primaryShower_hit_X[4847];   //[primaryShower_nHits]
   Double_t        primaryShower_hit_Y[4847];   //[primaryShower_nHits]
   Double_t        primaryShower_hit_Z[4847];   //[primaryShower_nHits]
   Double_t        primaryShower_hit_pitch[4847];   //[primaryShower_nHits]
   Double_t        primaryShower_hit_cnn[4847];   //[primaryShower_nHits]
   Double_t        primaryMomentumByRangeProton;
   Double_t        primaryMomentumByRangeMuon;
   Double_t        primaryKineticEnergy[3];
   Double_t        primaryRange[3];
   Int_t           primarynCal;
   Double_t        primarydQdx[4847];   //[primarynCal]
   Double_t        primary_calX[4847];   //[primarynCal]
   Double_t        primary_calY[4847];   //[primarynCal]
   Double_t        primary_calZ[4847];   //[primarynCal]
   Double_t        primary_cal_pitch[4847];   //[primarynCal]
   Double_t        primarydEdx[4847];   //[primarynCal]
   Double_t        primaryResidualRange[4847];   //[primarynCal]
   Double_t        primaryT0;
   Int_t           NDAUGHTERS;
   Double_t        daughterVertex[3];
   Int_t           daughterIstrack[27];   //[NDAUGHTERS]
   Int_t           daughterIsshower[27];   //[NDAUGHTERS]
   Int_t           daughterNHits[27];   //[NDAUGHTERS]
   Double_t        daughterTheta[27];   //[NDAUGHTERS]
   Double_t        daughterPhi[27];   //[NDAUGHTERS]
   Double_t        daughterLength[27];   //[NDAUGHTERS]
   Double_t        daughterEndPosition[27][3];   //[NDAUGHTERS]
   Double_t        daughterStartPosition[27][3];   //[NDAUGHTERS]
   Double_t        daughterStartDirection[27][3];   //[NDAUGHTERS]
   Double_t        daughterEndDirection[27][3];   //[NDAUGHTERS]
   Double_t        daughterOpeningAngle[27];   //[NDAUGHTERS]
   Double_t        daughterShowerBestPlane[27];   //[NDAUGHTERS]
   Int_t           daughterID[27];   //[NDAUGHTERS]
   Double_t        daughterT0[27];   //[NDAUGHTERS]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_timestamp;   //!
   TBranch        *b_Nactivefembs;   //!
   TBranch        *b_beamtrigger;   //!
   TBranch        *b_beamCheckIsMatched;   //!
   TBranch        *b_tof;   //!
   TBranch        *b_cerenkovStatus;   //!
   TBranch        *b_cerenkovTime;   //!
   TBranch        *b_cerenkovPressure;   //!
   TBranch        *b_beamtrackMomentum;   //!
   TBranch        *b_beamtrackEnergy;   //!
   TBranch        *b_beamtrackPos;   //!
   TBranch        *b_beamtrackEndPos;   //!
   TBranch        *b_beamtrackDir;   //!
   TBranch        *b_beamtrackTime;   //!
   TBranch        *b_beamtrackPdg;   //!
   TBranch        *b_beamtrackID;   //!
   TBranch        *b_beam_ntrjPoints;   //!
   TBranch        *b_beamtrackNDaughters;   //!
   TBranch        *b_beamtrackPos_at;   //!
   TBranch        *b_beamtrackMom_at;   //!
   TBranch        *b_beamtrackEndProcess;   //!
   TBranch        *b_primaryVertex;   //!
   TBranch        *b_primaryIsBeamparticle;   //!
   TBranch        *b_primaryIstrack;   //!
   TBranch        *b_primaryIsshower;   //!
   TBranch        *b_primaryBDTScore;   //!
   TBranch        *b_primaryNHits;   //!
   TBranch        *b_primaryTheta;   //!
   TBranch        *b_primaryPhi;   //!
   TBranch        *b_primaryLength;   //!
   TBranch        *b_primaryMomentum;   //!
   TBranch        *b_primaryEndMomentum;   //!
   TBranch        *b_primaryEndPosition;   //!
   TBranch        *b_primaryStartPosition;   //!
   TBranch        *b_primaryEndDirection;   //!
   TBranch        *b_primaryStartDirection;   //!
   TBranch        *b_primaryOpeningAngle;   //!
   TBranch        *b_primaryID;   //!
   TBranch        *b_primaryTruth_E;   //!
   TBranch        *b_primaryTruth_vtx;   //!
   TBranch        *b_primaryTruth_pdg;   //!
   TBranch        *b_primaryTruth_trkID;   //!
   TBranch        *b_primaryShowerBestPlane;   //!
   TBranch        *b_primaryShowerCharge;   //!
   TBranch        *b_primaryShowerEnergy;   //!
   TBranch        *b_primaryShowerMIPEnergy;   //!
   TBranch        *b_primaryShower_nHits;   //!
   TBranch        *b_primaryShower_hit_q;   //!
   TBranch        *b_primaryShower_hit_w;   //!
   TBranch        *b_primaryShower_hit_t;   //!
   TBranch        *b_primaryShower_hit_X;   //!
   TBranch        *b_primaryShower_hit_Y;   //!
   TBranch        *b_primaryShower_hit_Z;   //!
   TBranch        *b_primaryShower_hit_pitch;   //!
   TBranch        *b_primaryShower_hit_cnn;   //!
   TBranch        *b_primaryMomentumByRangeProton;   //!
   TBranch        *b_primaryMomentumByRangeMuon;   //!
   TBranch        *b_primaryKineticEnergy;   //!
   TBranch        *b_primaryRange;   //!
   TBranch        *b_primarynCal;   //!
   TBranch        *b_primarydQdx;   //!
   TBranch        *b_primary_calX;   //!
   TBranch        *b_primary_calY;   //!
   TBranch        *b_primary_calZ;   //!
   TBranch        *b_primary_cal_pitch;   //!
   TBranch        *b_primarydEdx;   //!
   TBranch        *b_primaryResidualRange;   //!
   TBranch        *b_primaryT0;   //!
   TBranch        *b_NDAUGHTERS;   //!
   TBranch        *b_daughterVertex;   //!
   TBranch        *b_daughterIstrack;   //!
   TBranch        *b_daughterIsshower;   //!
   TBranch        *b_daughterNHits;   //!
   TBranch        *b_daughterTheta;   //!
   TBranch        *b_daughterPhi;   //!
   TBranch        *b_daughterLength;   //!
   TBranch        *b_daughterEndPosition;   //!
   TBranch        *b_daughterStartPosition;   //!
   TBranch        *b_daughterStartDirection;   //!
   TBranch        *b_daughterEndDirection;   //!
   TBranch        *b_daughterOpeningAngle;   //!
   TBranch        *b_daughterShowerBestPlane;   //!
   TBranch        *b_daughterID;   //!
   TBranch        *b_daughterT0;   //!

   AnalysisTree(TTree *tree=0);
   virtual ~AnalysisTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalysisTree_cxx
AnalysisTree::AnalysisTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/dune/data/users/higuera/data/r5834_ana.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/dune/data/users/higuera/data/r5834_ana.root");
      }
      f->GetObject("PandoraBeam",tree);

   }
   Init(tree);
}

AnalysisTree::~AnalysisTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalysisTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalysisTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnalysisTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   beamtrackEndProcess = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   fChain->SetBranchAddress("Nactivefembs", Nactivefembs, &b_Nactivefembs);
   fChain->SetBranchAddress("beamtrigger", &beamtrigger, &b_beamtrigger);
   fChain->SetBranchAddress("beamCheckIsMatched", &beamCheckIsMatched, &b_beamCheckIsMatched);
   fChain->SetBranchAddress("tof", &tof, &b_tof);
   fChain->SetBranchAddress("cerenkovStatus", cerenkovStatus, &b_cerenkovStatus);
   fChain->SetBranchAddress("cerenkovTime", cerenkovTime, &b_cerenkovTime);
   fChain->SetBranchAddress("cerenkovPressure", cerenkovPressure, &b_cerenkovPressure);
   fChain->SetBranchAddress("beamtrackMomentum", &beamtrackMomentum, &b_beamtrackMomentum);
   fChain->SetBranchAddress("beamtrackEnergy", &beamtrackEnergy, &b_beamtrackEnergy);
   fChain->SetBranchAddress("beamtrackPos", beamtrackPos, &b_beamtrackPos);
   fChain->SetBranchAddress("beamtrackEndPos", beamtrackEndPos, &b_beamtrackEndPos);
   fChain->SetBranchAddress("beamtrackDir", beamtrackDir, &b_beamtrackDir);
   fChain->SetBranchAddress("beamtrackTime", &beamtrackTime, &b_beamtrackTime);
   fChain->SetBranchAddress("beamtrackPdg", &beamtrackPdg, &b_beamtrackPdg);
   fChain->SetBranchAddress("beamtrackID", &beamtrackID, &b_beamtrackID);
   fChain->SetBranchAddress("beam_ntrjPoints", &beam_ntrjPoints, &b_beam_ntrjPoints);
   fChain->SetBranchAddress("beamtrackNDaughters", &beamtrackNDaughters, &b_beamtrackNDaughters);
   fChain->SetBranchAddress("beamtrackPos_at", &beamtrackPos_at, &b_beamtrackPos_at);
   fChain->SetBranchAddress("beamtrackMom_at", &beamtrackMom_at, &b_beamtrackMom_at);
   fChain->SetBranchAddress("beamtrackEndProcess", &beamtrackEndProcess, &b_beamtrackEndProcess);
   fChain->SetBranchAddress("primaryVertex", primaryVertex, &b_primaryVertex);
   fChain->SetBranchAddress("primaryIsBeamparticle", &primaryIsBeamparticle, &b_primaryIsBeamparticle);
   fChain->SetBranchAddress("primaryIstrack", &primaryIstrack, &b_primaryIstrack);
   fChain->SetBranchAddress("primaryIsshower", &primaryIsshower, &b_primaryIsshower);
   fChain->SetBranchAddress("primaryBDTScore", &primaryBDTScore, &b_primaryBDTScore);
   fChain->SetBranchAddress("primaryNHits", &primaryNHits, &b_primaryNHits);
   fChain->SetBranchAddress("primaryTheta", &primaryTheta, &b_primaryTheta);
   fChain->SetBranchAddress("primaryPhi", &primaryPhi, &b_primaryPhi);
   fChain->SetBranchAddress("primaryLength", &primaryLength, &b_primaryLength);
   fChain->SetBranchAddress("primaryMomentum", &primaryMomentum, &b_primaryMomentum);
   fChain->SetBranchAddress("primaryEndMomentum", &primaryEndMomentum, &b_primaryEndMomentum);
   fChain->SetBranchAddress("primaryEndPosition", primaryEndPosition, &b_primaryEndPosition);
   fChain->SetBranchAddress("primaryStartPosition", primaryStartPosition, &b_primaryStartPosition);
   fChain->SetBranchAddress("primaryEndDirection", primaryEndDirection, &b_primaryEndDirection);
   fChain->SetBranchAddress("primaryStartDirection", primaryStartDirection, &b_primaryStartDirection);
   fChain->SetBranchAddress("primaryOpeningAngle", &primaryOpeningAngle, &b_primaryOpeningAngle);
   fChain->SetBranchAddress("primaryID", &primaryID, &b_primaryID);
   fChain->SetBranchAddress("primaryTruth_E", &primaryTruth_E, &b_primaryTruth_E);
   fChain->SetBranchAddress("primaryTruth_vtx", primaryTruth_vtx, &b_primaryTruth_vtx);
   fChain->SetBranchAddress("primaryTruth_pdg", &primaryTruth_pdg, &b_primaryTruth_pdg);
   fChain->SetBranchAddress("primaryTruth_trkID", &primaryTruth_trkID, &b_primaryTruth_trkID);
   fChain->SetBranchAddress("primaryShowerBestPlane", &primaryShowerBestPlane, &b_primaryShowerBestPlane);
   fChain->SetBranchAddress("primaryShowerCharge", &primaryShowerCharge, &b_primaryShowerCharge);
   fChain->SetBranchAddress("primaryShowerEnergy", &primaryShowerEnergy, &b_primaryShowerEnergy);
   fChain->SetBranchAddress("primaryShowerMIPEnergy", &primaryShowerMIPEnergy, &b_primaryShowerMIPEnergy);
   fChain->SetBranchAddress("primaryShower_nHits", &primaryShower_nHits, &b_primaryShower_nHits);
   fChain->SetBranchAddress("primaryShower_hit_q", primaryShower_hit_q, &b_primaryShower_hit_q);
   fChain->SetBranchAddress("primaryShower_hit_w", primaryShower_hit_w, &b_primaryShower_hit_w);
   fChain->SetBranchAddress("primaryShower_hit_t", primaryShower_hit_t, &b_primaryShower_hit_t);
   fChain->SetBranchAddress("primaryShower_hit_X", primaryShower_hit_X, &b_primaryShower_hit_X);
   fChain->SetBranchAddress("primaryShower_hit_Y", primaryShower_hit_Y, &b_primaryShower_hit_Y);
   fChain->SetBranchAddress("primaryShower_hit_Z", primaryShower_hit_Z, &b_primaryShower_hit_Z);
   fChain->SetBranchAddress("primaryShower_hit_pitch", primaryShower_hit_pitch, &b_primaryShower_hit_pitch);
   fChain->SetBranchAddress("primaryShower_hit_cnn", primaryShower_hit_cnn, &b_primaryShower_hit_cnn);
   fChain->SetBranchAddress("primaryMomentumByRangeProton", &primaryMomentumByRangeProton, &b_primaryMomentumByRangeProton);
   fChain->SetBranchAddress("primaryMomentumByRangeMuon", &primaryMomentumByRangeMuon, &b_primaryMomentumByRangeMuon);
   fChain->SetBranchAddress("primaryKineticEnergy", primaryKineticEnergy, &b_primaryKineticEnergy);
   fChain->SetBranchAddress("primaryRange", primaryRange, &b_primaryRange);
   fChain->SetBranchAddress("primarynCal", &primarynCal, &b_primarynCal);
   fChain->SetBranchAddress("primarydQdx", primarydQdx, &b_primarydQdx);
   fChain->SetBranchAddress("primary_calX", primary_calX, &b_primary_calX);
   fChain->SetBranchAddress("primary_calY", primary_calY, &b_primary_calY);
   fChain->SetBranchAddress("primary_calZ", primary_calZ, &b_primary_calZ);
   fChain->SetBranchAddress("primary_cal_pitch", primary_cal_pitch, &b_primary_cal_pitch);
   fChain->SetBranchAddress("primarydEdx", primarydEdx, &b_primarydEdx);
   fChain->SetBranchAddress("primaryResidualRange", primaryResidualRange, &b_primaryResidualRange);
   fChain->SetBranchAddress("primaryT0", &primaryT0, &b_primaryT0);
   fChain->SetBranchAddress("NDAUGHTERS", &NDAUGHTERS, &b_NDAUGHTERS);
   fChain->SetBranchAddress("daughterVertex", daughterVertex, &b_daughterVertex);
   fChain->SetBranchAddress("daughterIstrack", daughterIstrack, &b_daughterIstrack);
   fChain->SetBranchAddress("daughterIsshower", daughterIsshower, &b_daughterIsshower);
   fChain->SetBranchAddress("daughterNHits", daughterNHits, &b_daughterNHits);
   fChain->SetBranchAddress("daughterTheta", daughterTheta, &b_daughterTheta);
   fChain->SetBranchAddress("daughterPhi", daughterPhi, &b_daughterPhi);
   fChain->SetBranchAddress("daughterLength", daughterLength, &b_daughterLength);
   fChain->SetBranchAddress("daughterEndPosition", daughterEndPosition, &b_daughterEndPosition);
   fChain->SetBranchAddress("daughterStartPosition", daughterStartPosition, &b_daughterStartPosition);
   fChain->SetBranchAddress("daughterStartDirection", daughterStartDirection, &b_daughterStartDirection);
   fChain->SetBranchAddress("daughterEndDirection", daughterEndDirection, &b_daughterEndDirection);
   fChain->SetBranchAddress("daughterOpeningAngle", daughterOpeningAngle, &b_daughterOpeningAngle);
   fChain->SetBranchAddress("daughterShowerBestPlane", daughterShowerBestPlane, &b_daughterShowerBestPlane);
   fChain->SetBranchAddress("daughterID", daughterID, &b_daughterID);
   fChain->SetBranchAddress("daughterT0", daughterT0, &b_daughterT0);
   Notify();
}

Bool_t AnalysisTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalysisTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalysisTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalysisTree_cxx
