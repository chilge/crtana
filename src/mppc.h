//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////

#ifndef mppc_h
#define mppc_h

#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TCollection.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <sys/stat.h>
using namespace std;

class mppc 
{

public :
	const TString HOME = "/home/chris/Analysis/north_wall/anacode/";
	const TString DATA = "/home/chris/Data/north_wall/";

	const TString PLOTS = HOME+"Plots/";
	const TString CAL = HOME+"CalData/";
	const TString ANA = HOME+"AnaData/";

	//const TString FILEBASE = "data_ORMtest_North_";
	//directories to be determined in constructor
	TString plotDir;    //top level of plot directory
	TString calPlotDir; //subdir in plot dir for cal plots
	TString lyPlotDir;  //subdir in plot dir for ly plots
	TString calDir;     //subdir of Cal dir for ORM specific gain/ped data
	TString anaDir;     //subdir of AnaData dir form ORM specific ly fits
	TString ormDirName;

	TTree	  *fTreeInner = 0;   //pointer to hod-trig data tree
	TTree	  *fTreeOuter = 0;  //pointer to self-trig data tree

	// Declaration of leaf types
	UChar_t	 mac5;
	UShort_t chg[32];
	UInt_t	 ts0;
	UInt_t	 ts1;
	//UInt_t	 ts0_ref;
	//UInt_t	 ts1_ref;

	// List of branches
	TBranch	*b_mac5;
	TBranch	*b_chg;
	TBranch	*b_ts0;
	TBranch	*b_ts1;
  	//TBranch	*b_ts0_ref;
	//TBranch	*b_ts1_ref;

	//private function declarations
	mppc(TString file1, TString file2); //constructor
	virtual ~mppc(); //destructor
	virtual void     Init(TTree *tree, char type);
	virtual Bool_t   Notify();

	//public function declarations
	bool ParseConfig();
	bool VerifyMacAddresses();
	bool IsMacConfigured(int mac);
	bool IsChanConfigured(int mac, int chan);
	string FormatMacString(int mac);
	void PlotSpectrum(Int_t mac, Int_t chan, Int_t save_opt);
	void PlotSingleChanCutLoop(Int_t mac, Float_t cut);
	TH1F* PlotSingleChanCut(Int_t mac, Int_t chan, Float_t cut, Int_t save_opt);
	void ThresholdStats(char runtype, Float_t thresh);
	void ThresholdCorrelation(char runtype, Int_t ch, Int_t mac, Float_t thresh);
	Double_t* LangausFit(TH1F *h, Int_t save_opt);
	Double_t * PlotGainFit(Int_t mac, Int_t chan, Int_t save_opt, Double_t gain_seed);
	Double_t * FitGain(Int_t mac, TH1F *hs, Int_t chan, Int_t save_opt, Double_t gain_seed);
	Double_t * PlotPedFit(Int_t mac, Int_t chan , Int_t save_opt, bool useInit );
	Double_t * InitGain(Int_t mac);
	Double_t * InitPed(Int_t mac);
	Double_t * GetGain(Int_t mac);
	Double_t * GetPed(Int_t mac);
        Double_t * GetPedWidth(Int_t mac);
	void Cal(Int_t mac);
	void FillLYArr(Int_t mac, Int_t (&n)[32], Double_t (&peak)[32], Double_t (&fwhm)[32], Double_t(&w)[32],
		Double_t (&werr)[32], Double_t (&mpv)[32], Double_t (&mpverr)[32], Double_t (&norm)[32],
		Double_t (&normerr)[32], Double_t (&gsig)[32], Double_t (&gsigerr)[32] );
	void FillGainArr( Int_t mac, Double_t (&g)[32], Double_t (&gerr)[32], Double_t (&gped)[32],
		Double_t (&gpederr)[32], Double_t (&chi)[32], Double_t (&ndf)[32]);
	void FillPedArr(Int_t mac, Double_t (&mean)[32], Double_t (&meanerr)[32],
		Double_t (&sig)[32], Double_t (&sigerr)[32], Double_t (&cons)[32], Double_t (&conserr)[32]);
	vector<Int_t>* GetRawArr(char runtype, Int_t mac);
	vector<Double_t>* GetPEArr(char runtype, Int_t mac);
	Int_t GetNentries(char runtype, Int_t mac);
	Float_t GetMaxInRange(TH1F* htemp, Float_t low, Float_t high);
	Int_t GetMaxBinInRange(TH1F* htemp, Float_t low, Float_t high);
	void CountEmptyEvents(char runtype, int mac );
	void DumpEvents(char runtype, Int_t mac);
	void RunTime(char runtype, int input);
	void CorrelatedChannels(char runtype );
	void CheckEventPacket(char runtype );
	void TimeResolution(char runtype, Int_t mac);
	void MissedTimeStamps(char runtype, Int_t mac);
	void CountClockResets(char runtype, Int_t mac);
	void PlotRawSpectrum(Int_t mac, Int_t chan, Int_t save_opt);
	void CountOverflows(char runtype );
	void PlotADCTimeSlice(char runtype, Int_t mac, Int_t adclow, Int_t adchigh);
};

//constuctor
mppc::mppc(TString file1, TString file2)
{
	TString finner, fouter;
        if( (file1.Contains("inner")&&file2.Contains("inner")) ||
             (file1.Contains("outer")&&file2.Contains("outer")) ||
             (!file1.Contains("outer")&&!file1.Contains("inner")) ||
             (!file2.Contains("outer")&&!file2.Contains("inner")) ) {
		
		cout << "FATAL ERROR: expected one 'inner' and one 'outer'" 
                   << " file while both were the same or one type was missing...aborting" << endl;
		abort();
	}

	if(file1.Contains("inner")) {
		finner = file1;
		fouter = file2;
        }
	else {
		finner = file2;
		fouter = file1;
	}

	TFile *fin = TFile::Open(finner);
	TFile *fout = TFile::Open(fouter);
	string finstr = (string)finner;
	string foutstr = (string)fouter;
	while(finstr.find("/")!=UINT64_MAX)
		finstr=finstr.substr(finstr.find("/")+1,finstr.size());
	finstr=finstr.substr(0,finstr.find(".root"));
        while(foutstr.find("/")!=UINT64_MAX)
                foutstr=foutstr.substr(foutstr.find("/")+1,foutstr.size());
        foutstr=foutstr.substr(0,foutstr.find(".root"));

	/*TSystemDirectory dir(DATA,DATA);
	TList *files = dir.GetListOfFiles();
	TSystemFile *file;
	TString fname; 
	TIter nxt(files);

	TString prefix = FILEBASE + to_string(n) + "_South_" + to_string(s);
	TString pathHod=DATA, pathSelf=DATA;
	size_t ctrHod=0, ctrSelf=0;
	int base = 1e5;

	while ((file=(TSystemFile*)nxt())) { 
		fname = file->GetName(); 
		if (!file->IsDirectory() && fname.EndsWith(".root")&&fname.BeginsWith(prefix)) { 
			if(fname.Contains("HOD_trig")) {
				ctrHod++;
				cout << "found hod trig file: " << fname.Data() << endl; 
				pathHod+=fname;
			}
			if(fname.Contains("Self_trig")) {
				ctrSelf++;
				cout << "found self trig file: " << fname.Data() << endl; 
				pathSelf+=fname;
			}	
				
		}
	}

	//if user didn't shorten ORM ID from SN00#-00### to #00###
	if(ctrHod==0 && ctrSelf==0 && (n%base<n||s%base<s)) {
		nxt.Reset();
		prefix = FILEBASE;
		if(n%base<n){
			prefix+="SN00"+to_string(n/base)+"-00";
			if(n%base<100&&n%base>9)
				prefix+="0"+to_string(n%base);
			if(n%base<10)
				prefix+="00"+to_string(n%base);
		}
		else
			prefix+=to_string(n);
		prefix+="_South_";
		if(s%base<s){
			prefix+="SN00"+to_string(s/base)+"-00";
			if(s%base<100&&s%base>9)
				 prefix+="0"+to_string(s%base);
			if(s%base<10)
                                prefix+="00"+to_string(s%base);
		}
		else
			prefix+=to_string(s);

		pathHod=DATA; pathSelf=DATA;

		while ((file=(TSystemFile*)nxt())) { 
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(".root")&&fname.BeginsWith(prefix)) { 
				if(fname.Contains("HOD_trig")) {
					ctrHod++;
					cout << "found hod trig file: " << fname.Data() << endl; 
					pathHod+=fname;
				}
				if(fname.Contains("Self_trig")) {
					ctrSelf++;
					cout << "found self trig file: " << fname.Data() << endl; 
					pathSelf+=fname;
				}
			}
		}
	}

	if(ctrHod==0) cout << "WARNING: HOD TRIG RUN NOT FOUND using prefix " << prefix << " !" << endl;
	if(ctrSelf==0) cout << "WARNING: SELF TRIG RUN NOT FOUND using prefix " << prefix << " !" << endl;
	if(ctrHod>1) cout << "WARNING: MULTIPLE HOD TRIG RUNS FOUND!" << endl;
	if(ctrSelf>1) cout << "WARNING: MULTIPLE SELF TRIG RUNS FOUND!" << endl;

	TFile *fhod=0, *fself=0;
	if(ctrHod>0)  fhod = (TFile*)gROOT->GetListOfFiles()->FindObject(pathHod);
	if(ctrSelf>0) fself = (TFile*)gROOT->GetListOfFiles()->FindObject(pathSelf);

	if (!fhod || !fhod->IsOpen())
	{
 		if(ctrHod>0) fhod = new TFile(pathHod);
	}
	if (!fself || !fself->IsOpen())
	{
		if(ctrSelf>0) fself = new TFile(pathSelf);
	}

*/
	TTree* treeInner = 0;
	TTree* treeOuter = 0;
	/*if(ctrHod>0)*/  treeInner = (TTree*)fin->FindObjectAny("events");//,treeInner);
	/*if(ctrSelf>0)*/ treeOuter = (TTree*)fout->FindObjectAny("events");//,treeOuter);

	Init(treeInner,'i');
	Init(treeOuter,'o');

	struct stat sb;
	string com = "";
	TString dirName = (TString)finstr+"_AND_"+(TString)foutstr;
	TString plotDir = PLOTS + dirName;
	TString calPlotDir = plotDir + "/Cal";
	TString lyPlotDir = plotDir + "/LY";
	TString calDir = CAL + dirName;
	TString anaDir = ANA + dirName;

	if(1/*ctrHod>0&&ctrSelf>0i*/) {
		cout << "generating output directories with base " << dirName << "..." << endl;

		//check to see if plot directory exists, if not then create
		if (stat(string(plotDir).c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
			com = "mkdir "+(string)plotDir; std::system(com.c_str());
			cout << "   created directory '" << plotDir << "'" << endl;
		}
		else
			cout << "   plot directory already exists...skipping to next" << endl;

		//check to see if cal plot directory exists, if not then create
		if (stat(string(calPlotDir).c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
			com = "mkdir "+(string)calPlotDir; std::system(com.c_str());
			cout << "   created directory '" << calPlotDir << "'" << endl;
		}
		else 
			cout << "   calibration plot directory already exists...skipping to next" << endl;

		//check to see if LY plot directory exists, if not then create
		if (stat(string(lyPlotDir).c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
			com = "mkdir "+(string)lyPlotDir; std::system(com.c_str());
			cout << "   created directory '" << lyPlotDir << "'" << endl;
		}
		else
			cout << "   LY plot directory already exists...skipping to next" << endl;

		//check to see if calibration directory exists, if not then create
		if (stat(string(calDir).c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
			com = "mkdir "+(string)calDir; std::system(com.c_str());
			cout << "   created directory '" << calDir << "'" << endl;
		}
		else
			cout << "   calibration directory already exists...skipping to next" << endl;

		//check to see if analysis data directory exists, if not then create
		if (stat(string(anaDir).c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
			com = "mkdir "+(string)anaDir; std::system(com.c_str());
			cout << "   created directory '" << anaDir << "'" << endl;
		}
		else
			cout << "   analysis data directory already exists...skipping to next" << endl;
		cout << "directory generation complete" << endl;
	}
	else
		cout << "WARNING: missing one or both data files -> "
			<< "output directories not generated" << endl;

}

mppc::~mppc()
{
	if (fTreeInner)
		delete fTreeInner->GetCurrentFile();
	if (fTreeOuter)
		delete fTreeOuter->GetCurrentFile();
	return;
}

void mppc::Init(TTree *tree, char type)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

	// Set branch addresses and branch pointers
	if (!tree) return;

	if(type=='i') {
		fTreeInner = tree;

		fTreeInner->SetBranchAddress("mac5", &mac5, &b_mac5);
		fTreeInner->SetBranchAddress("adc", &chg, &b_chg);
		fTreeInner->SetBranchAddress("ts0", &ts0, &b_ts0);
		fTreeInner->SetBranchAddress("ts1", &ts1, &b_ts1);
		//fTreeInner->SetBranchAddress("ts0_ref", &ts0_ref, &b_ts0_ref);
		//fTreeInner->SetBranchAddress("ts1_ref", &ts1_ref, &b_ts1_ref);
	}
	if(type=='o') {
		fTreeOuter = tree;

		fTreeOuter->SetBranchAddress("mac5", &mac5, &b_mac5);
		fTreeOuter->SetBranchAddress("adc", &chg, &b_chg);
		fTreeOuter->SetBranchAddress("ts0", &ts0, &b_ts0);
		fTreeOuter->SetBranchAddress("ts1", &ts1, &b_ts1);
		//fTreeOuter->SetBranchAddress("ts0_ref", &ts0_ref, &b_ts0_ref);
		//fTreeOuter->SetBranchAddress("ts1_ref", &ts1_ref, &b_ts1_ref);
	}
   
	Notify();
}

Bool_t mppc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif
