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

string formatMac(int mac){
	string out;
	if(mac<10) out="00";
	else if(mac<100) out="0";
	out+=to_string(mac);	

	return out;
}

//main code body
int main() {

	string in, file1, file2;

	cout << "Welcome to the analyzer module!" << '\n'
	     << "You will be asked to enter the full path" << '\n'
	     << " to each of a pair of root ntuples, one"  << '\n'
	     << " for each CRT layer..." <<  endl;
	cout << "Enter the first file path." << endl;
	cin >> file1;
	cout << "Enter the second file path." << endl;
	cin >> file2;
	cout << "Attempting to analyze..." << endl;

	gROOT->Reset();

	mppc m(file1,file2);
	string fname = file1;
	while(fname.find("/")!=UINT64_MAX)
        	fname=fname.substr(fname.find("/")+1,fname.size());
	fname=fname.substr(0,fname.find("inner"));
	fname=fname.substr(0,fname.find("outer"));
	fname="./Trees/"+fname+".root";

	//check if file already exists, exit if it does
	if(!gSystem->AccessPathName(fname.c_str())) {
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
	Bool_t chconfig[32];
	int mac;

	//output file 
	cout << "Opened output file '" << fname << "'" << endl;
	TFile f(fname.c_str(),"RECREATE");
	TTree* tree = new TTree("orm","analysis results from raw ORM test data");
	
	//run info
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

		//filling local variables from cal/ana data files
		cout << " extracting analysis data..." << endl;
		m.FillPedArr(mac,pmean,pmeanerr,psig,psigerr,pconst,pconsterr);
		m.FillGainArr(mac,gain,gainerr,gped,gpederr,gchisqr,gndf);

		cout << " generating charge spectra..." << endl;
		TList *list = new TList();
		for(int ch=0; ch<32; ch++){
			if(chconfig[ch])
				list->Add(m.PlotSpectrum(mac,ch,0,0));
		}
		string listname="spectra"+formatMac(mac);
		list->Write(listname.c_str(),TObject::kSingleKey);

		cout << " filling tree..." << endl;
		tree->Fill();
	}
	cout << "writing tree..." << endl;
	tree->Write();	
	f.Close();
	return 0;
}

#endif
