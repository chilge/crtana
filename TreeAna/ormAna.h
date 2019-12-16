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
#include <vector>

// Header file for the classes stored in the TTree if any.

class ormAna {
public :
   TTree          *fTree;   //!pointer to the analyzed TTree or TChain
   TList          *fList;   //!point to list of lists of histograms

   std::vector<int>    fMacs;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
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

   // List of branches
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

   ormAna(TString fname);
   Int_t    GetEntry(Long64_t entry);
   void     Init();
   Int_t    GetNEntries();
   Int_t    LoadMacs();
   Int_t    FillList(TFile* f);
};

#endif

#ifdef ormAna_cxx
ormAna::ormAna(TString fname) : fTree(0) 
{
   if (!gSystem->AccessPathName(fname))
      cout << "input file found!" << endl; 
      //cout << "file name given by user does not exist! no tree loaded!" << endl;

   else  {
      TFile* f = new TFile("../Trees/crt_north_thr240_11dec2019_.root");
      f->GetObject("orm",fTree);

      Init();
      LoadMacs();
      if(fMacs.size()==0) cout << "FEB mac5 list is empty!" << endl;
      FillList(f);
      if(fList->GetSize()==0) cout << "Hist list is empty!" << endl;
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
}

Int_t ormAna::GetNEntries(){
   return fTree->GetEntriesFast();
}

Int_t ormAna::LoadMacs(){
        Int_t nentries = this->GetNEntries();
        for(Int_t ientry=0; ientry<nentries; ientry++){
                this->GetEntry(ientry);
                fMacs.push_back(this->mac5);
        }
        return fMacs.size();
}

Int_t ormAna::FillList(TFile* f){
	if(fMacs.size()==0) 
		cout << "error in FillList: FEB Mac5 list is empty!" << endl;
	fList = new TList();
	for(int i=0; i<fMacs.size(); i++) {
		string name="spectra";
		if(fMacs[i]<10) name+="00";
		else if(fMacs[i]<100) name+="0";
		name+=to_string(fMacs[i]);
		TList* l = (TList*)f->FindObjectAny(name.c_str());
		if(l->IsEmpty()) {
			cout << name << " is empty?" << endl;
			continue;
		}
		l->SetName(name.c_str());
		fList->Add(l);
	}
	return fList->GetSize();
}

#endif // #ifdef ormAna_cxx
