//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul 17 15:13:47 2019 by ROOT version 6.08/06
// from TTree orm/analysis results from raw ORM test data
// found on file: ../Trees/ormN_101_ormS_54.root
//////////////////////////////////////////////////////////

#ifndef ormAna_h
#define ormAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ormAna {
public :
   TTree          *fTree;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           ormN;
   Int_t           ormS;
   Float_t         effN;
   Float_t         effS;
   Int_t           mac5;
   UChar_t         chConfig[32];
   Double_t        gain[32];
   Double_t        gainErr[32];
   Double_t        gainPed[32];
   Double_t        gainPedErr[32];
   Double_t        gainChisqr[32];
   Double_t        gainNDF[32];
   Double_t        pedMean[32];
   Double_t        pedMeanErr[32];
   Double_t        pedSigma[32];
   Double_t        pedSigmaErr[32];
   Double_t        pedConst[32];
   Double_t        pedConstErr[32];
   Int_t           lyNFit[32];
   Double_t        lyPeak[32];
   Double_t        lyPeakErr[32];
   Double_t        lyFWHM[32];
   Double_t        lyLandauWidth[32];
   Double_t        lyLandauWidthErr[32];
   Double_t        lyLandauMPV[32];
   Double_t        lyLandauMPVErr[32];
   Double_t        lyNorm[32];
   Double_t        lyNormErr[32];
   Double_t        lyGaussSigma[32];
   Double_t        lyGaussSigmaErr[32];

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_ormN;   //!
   TBranch        *b_ormS;   //!
   TBranch        *b_effN;   //!
   TBranch        *b_effS;   //!
   TBranch        *b_mac;   //!
   TBranch        *b_chConfig;   //!
   TBranch        *b_gain;   //!
   TBranch        *b_gainErr;   //!
   TBranch        *b_gainPed;   //!
   TBranch        *b_gainPedErr;   //!
   TBranch        *b_gainChisqr;   //!
   TBranch        *b_gainNDF;   //!
   TBranch        *b_pedMean;   //!
   TBranch        *b_pedMeanErr;   //!
   TBranch        *b_pedSigma;   //!
   TBranch        *b_pedSigmaErr;   //!
   TBranch        *b_pedConst;   //!
   TBranch        *b_pedConstErr;   //!
   TBranch        *b_lyNFit;   //!
   TBranch        *b_lyPeak;   //!
   TBranch        *b_lyPeakErr;   //!
   TBranch        *b_lyFWHM;   //!
   TBranch        *b_lyLandauWidth;   //!
   TBranch        *b_lyLandauWidthErr;   //!
   TBranch        *b_lyLandauMPV;   //!
   TBranch        *b_lyLandauMPVErr;   //!
   TBranch        *b_lyNorm;   //!
   TBranch        *b_lyNormErr;   //!
   TBranch        *b_lyGaussSigma;   //!
   TBranch        *b_lyGaussSigmaErr;   //!

   ormAna(TString fname);
   Int_t    GetEntry(Long64_t entry);
   void     Init();

   //user functions
   void lyVsGain(int mac);

};

#endif

#ifdef ormAna_cxx
ormAna::ormAna(TString fname) : fTree(0) 
{
   if (gSystem->AccessPathName(fname)) 
      cout << "file name given by user does not exist! no tree loaded!" << endl;

   else  {
      TFile* f = new TFile("../Trees/ormN_101_ormS_54.root");
      f->GetObject("orm",fTree);
      Init();
   }
}


Int_t ormAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fTree) return 0;
   return fTree->GetEntry(entry);
}

void ormAna::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.

   // Set branch addresses and branch pointers
   fTree->SetBranchAddress("run", &run, &b_run);
   fTree->SetBranchAddress("ormN", &ormN, &b_ormN);
   fTree->SetBranchAddress("ormS", &ormS, &b_ormS);
   fTree->SetBranchAddress("effN", &effN, &b_effN);
   fTree->SetBranchAddress("effS", &effS, &b_effS);
   fTree->SetBranchAddress("mac5", &mac5, &b_mac);
   fTree->SetBranchAddress("chConfig", chConfig, &b_chConfig);
   fTree->SetBranchAddress("gain", gain, &b_gain);
   fTree->SetBranchAddress("gainErr", gainErr, &b_gainErr);
   fTree->SetBranchAddress("gainPed", gainPed, &b_gainPed);
   fTree->SetBranchAddress("gainPedErr", gainPedErr, &b_gainPedErr);
   fTree->SetBranchAddress("gainChisqr", gainChisqr, &b_gainChisqr);
   fTree->SetBranchAddress("gainNDF", gainNDF, &b_gainNDF);
   fTree->SetBranchAddress("pedMean", pedMean, &b_pedMean);
   fTree->SetBranchAddress("pedMeanErr", pedMeanErr, &b_pedMeanErr);
   fTree->SetBranchAddress("pedSigma", pedSigma, &b_pedSigma);
   fTree->SetBranchAddress("pedSigmaErr", pedSigmaErr, &b_pedSigmaErr);
   fTree->SetBranchAddress("pedConst", pedConst, &b_pedConst);
   fTree->SetBranchAddress("pedConstErr", pedConstErr, &b_pedConstErr);
   fTree->SetBranchAddress("lyNFit", lyNFit, &b_lyNFit);
   fTree->SetBranchAddress("lyPeak", lyPeak, &b_lyPeak);
   fTree->SetBranchAddress("lyPeakErr", lyPeakErr, &b_lyPeakErr);
   fTree->SetBranchAddress("lyFWHM", lyFWHM, &b_lyFWHM);
   fTree->SetBranchAddress("lyLandauWidth", lyLandauWidth, &b_lyLandauWidth);
   fTree->SetBranchAddress("lyLandauWidthErr", lyLandauWidthErr, &b_lyLandauWidthErr);
   fTree->SetBranchAddress("lyLandauMPV", lyLandauMPV, &b_lyLandauMPV);
   fTree->SetBranchAddress("lyLandauMPVErr", lyLandauMPVErr, &b_lyLandauMPVErr);
   fTree->SetBranchAddress("lyNorm", lyNorm, &b_lyNorm);
   fTree->SetBranchAddress("lyNormErr", lyNormErr, &b_lyNormErr);
   fTree->SetBranchAddress("lyGaussSigma", lyGaussSigma, &b_lyGaussSigma);
   fTree->SetBranchAddress("lyGaussSigmaErr", lyGaussSigmaErr, &b_lyGaussSigmaErr);
}

#endif // #ifdef ormAna_cxx
