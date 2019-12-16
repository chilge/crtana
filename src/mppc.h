//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////

#ifndef mppc_h
#define mppc_h

#include <TSystem.h>
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
	UShort_t   flags;
	//UInt_t	 ts0_ref;
	//UInt_t	 ts1_ref;
	ULong_t  run_start_time;
	ULong_t  this_poll_start;
	ULong_t  this_poll_end;
	ULong_t  last_poll_start;
	ULong_t  last_poll_end;
	ULong_t  fragment_timestamp;

	// List of branches
	TBranch	*b_mac5;
	TBranch	*b_chg;
	TBranch	*b_ts0;
	TBranch	*b_ts1;
	TBranch *b_flags;
  	//TBranch	*b_ts0_ref;
	//TBranch	*b_ts1_ref;
	TBranch *b_run_start_time;
	TBranch *b_this_poll_start;
	TBranch *b_this_poll_end;
	TBranch *b_last_poll_start;
	TBranch *b_last_poll_end;
	TBranch *b_fragment_timestamp;

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
	TH1F* PlotSpectrum(Int_t mac, Int_t chan, Int_t save_opt, Bool_t overlay);
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
	void MacSpectraOverlay(Int_t mac);
	void PlotTS0(Int_t mac);
	void EventRate(int mac);
};

//constuctor
mppc::mppc(TString file1, TString file2)
{
	TString finner, fouter;

	if(gSystem->AccessPathName(file1)) {
		cout << file1 << " not found!" << endl;
		throw std::exception();
	}
        if(gSystem->AccessPathName(file2)) {
                cout << file2 << " not found!" << endl;
                throw std::exception();
        }

        if( (file1.Contains("inner")&&file2.Contains("inner")) ||
             (file1.Contains("outer")&&file2.Contains("outer")) ||
             (!file1.Contains("outer")&&!file1.Contains("inner")) ||
             (!file2.Contains("outer")&&!file2.Contains("inner")) ) {
		
		cout << "expected one 'inner' and one 'outer'" <<
                   " file while both were the same or one type was missing...aborting" << endl;
		throw std::exception();
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

	TTree* treeInner = 0;
	TTree* treeOuter = 0;
	treeInner = (TTree*)fin->FindObjectAny("events");
	treeOuter = (TTree*)fout->FindObjectAny("events");

	Init(treeInner,'i');
	Init(treeOuter,'o');

	struct stat sb;
	string com = "";
	TString dirName = (TString)finstr+"_AND_"+(TString)foutstr;
	plotDir = PLOTS + dirName;
	calPlotDir = plotDir + "/Cal";
	lyPlotDir = plotDir + "/LY";
	calDir = CAL + dirName;
	anaDir = ANA + dirName;

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
		fTreeInner->SetBranchAddress("flags",&flags,&b_flags);
		//fTreeInner->SetBranchAddress("ts0_ref", &ts0_ref, &b_ts0_ref);
		//fTreeInner->SetBranchAddress("ts1_ref", &ts1_ref, &b_ts1_ref);
		fTreeInner->SetBranchAddress("run_start_time",     &run_start_time,     &b_run_start_time);
		fTreeInner->SetBranchAddress("this_poll_start",    &this_poll_start,    &b_this_poll_start);
		fTreeInner->SetBranchAddress("this_poll_end",      &this_poll_end,      &b_this_poll_end);
		fTreeInner->SetBranchAddress("last_poll_start",    &last_poll_start,    &b_last_poll_start);
		fTreeInner->SetBranchAddress("last_poll_end",      &last_poll_end,      &b_last_poll_end);
		fTreeInner->SetBranchAddress("fragment_timestamp", &fragment_timestamp, &b_&fragment_timestamp);

	}
	if(type=='o') {
		fTreeOuter = tree;

		fTreeOuter->SetBranchAddress("mac5", &mac5, &b_mac5);
		fTreeOuter->SetBranchAddress("adc", &chg, &b_chg);
		fTreeOuter->SetBranchAddress("ts0", &ts0, &b_ts0);
		fTreeOuter->SetBranchAddress("ts1", &ts1, &b_ts1);
		fTreeOuter->SetBranchAddress("flags",&flags,&b_flags);
		//fTreeOuter->SetBranchAddress("ts0_ref", &ts0_ref, &b_ts0_ref);
		//fTreeOuter->SetBranchAddress("ts1_ref", &ts1_ref, &b_ts1_ref);
		fTreeOuter->SetBranchAddress("run_start_time",     &run_start_time,     &b_run_start_time);
		fTreeOuter->SetBranchAddress("this_poll_start",    &this_poll_start,    &b_this_poll_start);
		fTreeOuter->SetBranchAddress("this_poll_end",      &this_poll_end,      &b_this_poll_end);
		fTreeOuter->SetBranchAddress("last_poll_start",    &last_poll_start,    &b_last_poll_start);
		fTreeOuter->SetBranchAddress("last_poll_end",      &last_poll_end,      &b_last_poll_end);
		fTreeOuter->SetBranchAddress("fragment_timestamp", &fragment_timestamp, &b_&fragment_timestamp);
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
