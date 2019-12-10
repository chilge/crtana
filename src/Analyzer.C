#ifndef Analyzer_cxx
#define Analyzer_cxx

#include "mppc.C"

#include <string>
#include <iostream>

using namespace std;

//checks if input char is a number, returns it as an int
int charToInt(char c) {
	int i = (int)(c-'0');
	if(i<0||i>9)
		return -1;
	return i;
}

//checks if input string is an integer
bool isInt(string s){
	for(size_t i=0; i<s.length(); i++)
		if(charToInt(s[i])<0)
			return false;
	return true;
}

// checks if input string is a float
bool isFloat(string s){
	bool isfloat=true;
	if(s.find(".")==string::npos) {
		if(!isInt(s))
			isfloat=false;
	}
	else {
		if(!isInt(s.substr(0,s.find("."))))
			isfloat=false;
		if(!isInt(s.substr(s.find(".")+1,s.length())))
			isfloat=false;
	}
	if(!isfloat)
		cout << "non-numeric entry given...expected float" << endl;
	return isfloat;

}

//main code body
int main() {

	string in;

	cout << "." << endl;
	cin >> in;

	gROOT->Reset();

	mppc m();
	if(m.fTreeHod==0||m.fTreeSelf==0) {
		cout << "ERROR: failed to initialize trees. Exiting..." << endl;
		return 0;
	}

	TString fname = "Trees/"+m.ormDirName+".root"; //output file path

	//check if file already exists, exit if it does
	if(!gSystem->AccessPathName(fname)) {
		cout << "ERROR: these data appear to already have been processed..." << '\n' 
			<< "    output tree found: " << fname << '\n'
			<< "    Exiting..." << endl;
		return 1;
	}

	m.ParseConfig();
	cout << "Found " << febchans.size() << " FEBs in config file" << endl;

	m.VerifyMacAddresses();

	Double_t gain[32], gainerr[32], gped[32], gpederr[32], gchisqr[32], gndf[32];
	Double_t pmean[32], pmeanerr[32], psig[32], psigerr[32], pconst[32], pconsterr[32];
	Int_t lynfit[32];
	Double_t lypeak[32], lypeakerr[32], lyfwhm[32], lynorm[32], lynormerr[32]; 
	Double_t lyw[32], lywerr[32], lympv[32], lympverr[32], lygsig[32], lygsigerr[32];
	Bool_t chconfig[32];
	int mac;

	//output file 
	cout << "Opened output file '" << fname << "'" << endl;
	TFile f(fname,"RECREATE");
	TTree* tree = new TTree("orm","analysis results from raw ORM test data");
	
	//run info
	//tree->Branch("run",              &run,     "run/I");
	tree->Branch("mac5",             &mac,     "mac/I");
	tree->Branch("chConfig",         chconfig, "chConfig[32]/b");

	//gain fit info
	tree->Branch("gain",             gain,     "gain[32]/D");
	tree->Branch("gainErr",          gainerr,  "gainErr[32]/D");
	tree->Branch("gainPed",          gped,     "gainPed[32]/D");
	tree->Branch("gainPedErr",       gpederr,  "gainPedErr[32]/D");
	tree->Branch("gainChisqr",       gchisqr,  "gainChisqr[32]/D");
	tree->Branch("gainNDF",          gndf,     "gainNDF[32]/D");

	//ped fit info
	tree->Branch("pedMean",          pmean,    "pedMean[32]/D");
	tree->Branch("pedMeanErr",       pmeanerr, "pedMeanErr[32]/D");
	tree->Branch("pedSigma",         psig,     "pedSigma[32]/D");
	tree->Branch("pedSigmaErr",      psigerr,  "pedSigmaErr[32]/D");
	tree->Branch("pedConst",         pconst,   "pedConst[32]/D");
	tree->Branch("pedConstErr",      pconsterr,"pedConstErr[32]/D");

	//LY fit info
	tree->Branch("lyNFit",           lynfit,   "lyNFit[32]/I");
	tree->Branch("lyPeak",           lypeak,   "lyPeak[32]/D");
	tree->Branch("lyPeakErr",        lypeakerr,"lyPeakErr[32]/D");
	tree->Branch("lyFWHM",           lyfwhm,   "lyFWHM[32]/D");
	tree->Branch("lyLandauWidth",    lyw,      "lyLandauWidth[32]/D");
	tree->Branch("lyLandauWidthErr", lywerr,   "lyLandauWidthErr[32]/D");
	tree->Branch("lyLandauMPV",      lympv,    "lyLandauMPV[32]/D");
	tree->Branch("lyLandauMPVErr",   lympverr, "lyLandauMPVErr[32]/D");
	tree->Branch("lyNorm",           lynorm,   "lyNorm[32]/D");
	tree->Branch("lyNormErr",        lynormerr,"lyNormErr[32]/D");
	tree->Branch("lyGaussSigma",     lygsig,   "lyGaussSigma[32]/D");
	tree->Branch("lyGaussSigmaErr",  lygsigerr,"lyGaussSigmaErr[32]/D");

	for(auto feb=febchans.begin(); feb!=febchans.end(); feb++) 
	{
		mac = (*feb).first;

		for(size_t i=0; i<32; i++)
			chconfig[i]=kFALSE;
		for(auto ich=(*feb).second.begin(); ich!=(*feb).second.end(); ich++)
			chconfig[*ich] = kTRUE;

		cout << "processing data for mac5 " << mac << endl;
		cout << " calibrating..." << endl;
		m.Cal(mac); //pedestal then gain calibration
		cout << " LY fit..." << endl;
		m.PlotSingleChanCutLoop('h',mac,2.5); //LY fits, single chan cut @ thr=2.5 PE

		//filling local variables from cal/ana data files
		cout << " extracting analysis data..." << endl;
		m.FillPedArr(mac,pmean,pmeanerr,psig,psigerr,pconst,pconsterr);
		m.FillGainArr(mac,gain,gainerr,gped,gpederr,gchisqr,gndf);
		m.FillLYArr(mac,lynfit,lypeak,lyfwhm,lyw,lywerr,lympv,lympverr,lynorm,lynormerr,lygsig,lygsigerr);
		for(size_t ch=0; ch<32; ch++) {
			if(lynfit[ch]==0)
				continue;
			lypeakerr[ch] = (lyfwhm[ch]/2.355)/sqrt(lynfit[ch]);
		}

		cout << " filling tree..." << endl;
		tree->Fill();
	}
	cout << "writing tree..." << endl;
	tree->Write();	
	f.Close();
	return 0;
}

#endif
