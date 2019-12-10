#define ormAna_cxx
#include "ormAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <set>

//   In a ROOT session, you can do:
//      root> .L ormAna.C
//      root> ormAna a("path/to/file")
//      root> a.UserFunction1();
//      root> a.UserFunction2();  
//

using namespace std;

//we know what FEBs we have, but his could change...be aware
//there's a better way to do this, but I'm lazy and this works for now
set<int> febs = {85,86,87,213,214,215};

void ormAna::lyVsGain(int mac=-1) {

	if (fTree == 0) {
		cout << "no tree found!" << endl;
		return;
	}

	if(mac!=-1&&febs.find(mac)==febs.end()) {
		cout << "no data found for given mac5 address" << endl;
		return;
	}

	const int nentries = fTree->GetEntriesFast();

	//histogram title
	TString htitle = "LY vs. Gain, ";
	if(mac==-1)
		htitle+="All FEBs";
	else
		htitle+="mac5 "+to_string(mac);

	//histogram range, binning
	int gainLow = 40, gainHigh=90;
	int lyLow = 10, lyHigh = 30;
	int gainBinning = 1;
	int lyBinning = 1;
	int nbinsGain = (gainHigh-gainLow)/gainBinning;
	int nbinsLy = (lyHigh-lyLow)/lyBinning;

	TH2F* h = new TH2F("h",htitle,nbinsGain,gainLow,gainHigh,nbinsLy,lyLow,lyHigh);

	//histrogram styling
	h->GetXaxis()->SetTitle("gain [ADC/PE]");
	h->GetYaxis()->SetTitle("light yield [PE]");
	h->SetMarkerStyle(8);

	//loop over tree entries (6 per ORM pair, 1 per FEB)
	for (int ientry=0; ientry<nentries; ientry++) {

		GetEntry(ientry);

		//if not all FEBs selected and current entry is not FEB of interest, skip to next entry
		if(mac!=-1&&mac5!=mac)
			continue;

		//loop over FEB channels
		for(int ch=0; ch<32; ch++) {
			//if current channel not connected, skip to next
			if(!chConfig[ch])
				continue;
			h->Fill(gain[ch],lyPeak[ch]);
		}//for FEB channels
	}//for tree entries

	h->Draw();
	
}
