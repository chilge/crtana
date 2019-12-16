#define mppc_cxx

#include "mppc.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <iostream>
#include <fstream>
#include <TPaveText.h>
#include <TSystem.h>
#include <TF1.h>
#include <TMath.h>
#include <TImage.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <vector>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/types.h>

using namespace std;

//****************************************************************
//****************************************************************
//** Description: mppc class library for use with ORM tests     **
//** at Wideband test stand, FNAL, V1.0		                **
//** Last Edit:   10 July 2019	 				**
//** Authors:     Chris.Hilgenberg@colostate.edu                **
//** Intstitution: Colorado State University			**
//**							        **
//** Release notes: initial release			        **
//** 								**
//** To Use:					         	**
//** 1) cd <path_to_working_directory>			        **
//** 2) start ROOT session				        **
//** 3) enter the following in the command line		        **
//**      ROOT$  [1] .L mppc.C			                **
//**      ROOT$  [2] mppc m(northORM, southORM)	                **
//**      ROOT$  [3] m.<function_from_mppc_class_library>       **
//**      ROOT$  [4] m.<another_function...>		        **
//**      ROOT$  [#] .q					        **
//****************************************************************
//****************************************************************


//**************************************************************************************
//**************************************************************************************
//******		 configuration/validation methods	                  ******
//**************************************************************************************
//**************************************************************************************

//global variable storing mac5 addresses and active channels
//as specified in configuration file config.txt
std::map<int,std::set<int>> febchans;
std::map<int,char> mactolay;

//open configuration file and parse values
bool mppc::ParseConfig()
{ 
	febchans.clear();
	ifstream fin("./config.txt");
	if(!fin.good()) {
		cout << "configuration file not found!" << endl;
		return false;
	}
	string line;
	while( getline(fin,line)) {
		stringstream ss(line);
		string tmp;
		int mac=0, chan=0;
		ss >> tmp;

		if(tmp[tmp.size()-1]==':') {
			mac = atoi(tmp.c_str());
			if(mac<0 || mac>256) {
				cout << "in config file, mac address out of range: " << mac << endl;
				return false;
			}
		}
		else {
			cout << "unknown format in config file! expected '<feb>:'" << endl;
			return false;
		}

		if(ss.good()) while (ss.good()) {
			ss >> tmp;
			chan = atoi(tmp.c_str());
			if(chan<0||chan>31) {
				cout << "in config file, bad channel for mac5 " << mac << ": " << chan << endl;
				return false;
			}
			febchans[mac].insert(chan);
		}
		else {
			cout << "in config file, no channels specified for mac5 " << mac << endl;
			return false;
		}

	} //while line

	if(febchans.size()==0)
		return false;

	return true;
}

//VERIFY MUST BE CALLED PROIR TO CALLING ANY ANALYSIS FUNCTION!
//here is where the layer for each mac address is determined
bool mppc::VerifyMacAddresses(){

	if(febchans.size()==0) {
		cout << "FEB mac5 addresses not yet loaded from config file!"
			 << " Could not verify mac5s." << endl;
		return false;
	}

	bool verified = true;
	set<int> macs, inmacs, outmacs;
	int nentriesIn = fTreeInner->GetEntriesFast();
	int nentriesOut = fTreeOuter->GetEntriesFast();

	for(int i=0; i<0.5e6&&i<nentriesIn; i++){
		fTreeInner->GetEntry(i);
		macs.insert(mac5);
		inmacs.insert(mac5);
	}
        for(int i=0; i<0.5e6&&i<nentriesOut; i++){
                fTreeOuter->GetEntry(i);
                macs.insert(mac5);
		outmacs.insert(mac5);
        }

	for(auto it=inmacs.begin(); it!=inmacs.end(); it++)
		mactolay[*it] = 'i';
	for(auto it=outmacs.begin(); it!=outmacs.end(); it++)
		mactolay[*it] = 'o';

	auto itconfig = febchans.begin();
	auto itread = macs.begin();

	while(itconfig != febchans.end()){
		if(macs.find((*itconfig).first)==macs.end()) {
			cout << "mac5 " << (*itconfig).first << " not found in data file!" << endl;
			verified = false;
		}
		itconfig++;
	}
	while(itread != macs.end()) {
		if(febchans.find(*itread)==febchans.end()) {
			cout << "found mac5 not listed in config file: " << *itread << endl;
			verified = false;
		}
		itread++;
	} 

	if(verified) {
		itconfig = febchans.begin();
		cout << "verified FEBs: ";
		while(itconfig!=febchans.end()) {
			cout << (*itconfig).first;
			itconfig++;
			if(itconfig!=febchans.end())
				cout << ", ";
			else
				cout << endl;
		}
			
	}
	return verified;
}

bool mppc::IsMacConfigured(int mac) {
	if(febchans.size()==0)
		this->ParseConfig();
	return (febchans.find(mac) != febchans.end());
}

bool mppc::IsChanConfigured(int mac, int chan) {
	if(!IsMacConfigured(mac)) {
		cout << "ERROR: mac5 not configured!" << endl;
		return false;
	}
	return (febchans[mac].find(chan)!=febchans[mac].end());
}

string mppc::FormatMacString(int mac) {
	string ret;
	if(mac<0||mac>256) {
		cout << "ERROR in FormatMacString: given mac address out of range!" << endl;
		return "";
	}

	if(mac/10==0) ret = "00";
	else if (mac/100==0) ret = "0";
	ret+=to_string(mac);
	return ret;
}

//**************************************************************************************
//**************************************************************************************
//******                 BEGIN ANALYSIS/PLOT FUNCTION DEFINITIONS                ******
//**************************************************************************************
//**************************************************************************************

void mppc::PlotRawSpectrum(Int_t mac=0, Int_t chan=1, Int_t save_opt=0)
{
	if(!IsMacConfigured(mac)||mactolay.size()==0)
		cout << "WARNING in PlotRawSpectrum: mac5 address passed which is not configured" << endl;

	TTree* tree = 0;	
	if(mactolay[mac]=='i') tree = fTreeInner;
	else if(mactolay[mac]=='o') tree = fTreeOuter;
	else {
		cout << "ERROR in PlotRawSpectrum: unknown layer type found - " << mactolay[mac] << endl;
		return;
	}
	Int_t nentries = tree->GetEntriesFast();

	TString htitle = "mac5 "+to_string(mac)+", ch. "+to_string(chan);

	Float_t low = 1;
	Float_t high = 4071;
	Float_t binning = 5; //3/conv[1];
	Int_t nbins = (high-low)/binning;
	TString binstr="counts / "; binstr+=binning; binstr+=" ADC";

	TH1F *h = new TH1F("h",htitle,nbins,low,high);
	h->GetXaxis()->SetTitle("signal amplitude [ADC]");
	h->GetYaxis()->SetTitle(binstr);
	h->SetLineWidth(2);
	h->SetLineColor(kBlue);

	for ( Int_t ientry=0; ientry<nentries; ientry++ )
	{
		tree->GetEntry(ientry);

		if(mac5==mac) h->Fill(chg[chan]);
	}


	TCanvas *c = new TCanvas();
	c->SetGrid();
	c->SetLogy();
	h->Draw();

	/*TImage *img = TImage::Create();
	ofstream fout;
	TString fname = PLOTS; fname+="12feb/mac"; 
	fname+=mac; fname+="_ch"; fname+=chan; fname+="_raw.gif";
*/
	if(save_opt==1) {
		TImage *img = TImage::Create();
		ofstream fout;
		TString fname = PLOTS; fname+="12feb/mac";
		fname+=mac; fname+="_ch"; fname+=chan; fname+="_raw.gif";
		cout << "attempting to save plot as " << fname << endl;
		fout.open(fname);
		img->FromPad(c);
		img->WriteImage(fname);
		fout.close();
	}
}



TH1F* mppc::PlotSpectrum(Int_t mac=0, Int_t chan=1, Int_t save_opt=0, Bool_t overlay=0)
{
	if(!IsChanConfigured(mac,chan)||mactolay.size()==0) {
		cout << "ERROR in PlotSpectrum: channel passed which is not configured!" << endl;
		return nullptr;
	}

	TTree* tree;
	if(mactolay[mac]=='i') tree = fTreeInner;
	else if(mactolay[mac]=='o') tree = fTreeOuter;
	else {
		cout << "ERROR in PlotSpectrum: unknown runtype provided - " << mactolay[mac] << endl;
		return nullptr;
	}

	Double_t *conv = GetGain(mac);
	Double_t *ped = GetPed(mac);

	Int_t nentries = tree->GetEntriesFast();

	TString htitle = "ch. ";
	TString hobj = "hspectrum_"+FormatMacString(mac)+"_"; 
	htitle+=chan;
	if(chan<10) hobj+="0";
	hobj+=chan;

	Float_t low = -1.5;	
	Float_t high = 75.5;
	Float_t binning = 0.1; 
	Int_t nbins = (high-low)/binning;
	TString binningstr = binning; binningstr.Resize(4);
	TString binstr="counts / "+binningstr+" PE";

	TH1F *h = new TH1F(hobj,htitle,nbins,low,high);
	h->GetXaxis()->SetTitle("signal amplitude (PE)");
	h->GetYaxis()->SetTitle(binstr);
	//h->GetYaxis()->SetTitle("Hz / PE");
	h->SetLineWidth(2);
	//if(mac==85) h->SetLineColor(kRed);
	//if(mac==213) 
	h->SetLineColor(kBlue);

	for ( Int_t ientry=0; ientry<nentries; ientry++ )
	{
		tree->GetEntry(ientry);

		if(mac5==mac) h->Fill((chg[chan]-ped[chan])/conv[chan]);
	}

	//h->Sumw2();
	//h->Scale(1.0/(180-h->Integral()*2.2e-6));

	if(!overlay) {
		TCanvas *c = new TCanvas();
		c->SetGrid();
		c->SetLogy();
		h->Draw("hist");
	}
	return h;
}

void mppc::MacSpectraOverlay(Int_t mac=1) {

	Int_t colors[32];
	for(int i=0; i<32; i++) {
		if(i<23) colors[i] = 49-i;
		else colors[i] = i-22;
	}

	TCanvas* call = new TCanvas();
	bool first = true;
	for(int ch=0; ch<32; ch++) {
		if(!IsChanConfigured(mac,ch)) continue;
		TH1F* h = mppc::PlotSpectrum(mac,ch,0,1);
		call->cd();
		h->SetLineColor(colors[ch]);
		if(first) {
			h->Draw();
			first=false;
		}
		else
			h->Draw("same");
		//delete h;
	}
}

//******************************************************************************
//** Function plots converted signals for all 20 fibers in each FEB requiring
//** for each histogram, only events where all other channels are below given threshold
//*****************************************************************************
void mppc::PlotSingleChanCutLoop(Int_t mac=85, Float_t cut=1.0) 
{
	TTree* tree;
	if(mactolay[mac]=='i') tree = fTreeInner;
	else if(mactolay[mac]=='o') tree = fTreeOuter;
	else {
		cout << "error in PlotSingleChanCutLoop: unknown runtype provided - " << mactolay[mac] << endl;
		return;
	}

	if(!IsMacConfigured(mac)) {
		cout << "ERROR in PlotSingleChanCutLoop: mac not configured" << endl;
		return;
	}

	TString fname = anaDir+"/data_lyfit_thr"+to_string(cut).substr(0,3)+"_mac";
	fname+=FormatMacString(mac)+"_allchan.txt";
	ofstream fout;
	if (gSystem->AccessPathName(fname))
		fout.open(fname);

	else {
		cout << fname << " already exists -> new fit data not written" << endl;
		fout.close();
	}

	for(size_t ch=0; ch<32; ch++) {
		if(!IsChanConfigured(mac,ch)) {
			fout << 0 << '\n';
			continue;
		}
		this->PlotSingleChanCut(mac,ch,cut,1);
		if(fout.good()) {
			TString fana = anaDir+"/data_lyfit_thr"+to_string(cut).substr(0,3)+"_mac";
			fana+=FormatMacString(mac)+"_ch"+to_string(ch)+".txt";
			if (!gSystem->AccessPathName(fana)) {
				ifstream fin;
				fin.open(fana);
				while(fin.good()) {
					string tmp;
					fin >> tmp;
					fout << tmp << " ";
				}
				fout << '\n';
			}
			else 
				cout << "single chan ana file '" << fana << "not found...skipping" << endl;
			
		}
	}

	fout.close();
	return;

}//end def.

//******************************************************************************
//** Function plots converted signals for all 20 fibers in each FEB requiring
//** for each histogram, only events where all other 
//*****************************************************************************
TH1F* mppc::PlotSingleChanCut(Int_t mac=85, Int_t chan=1, Float_t cut=1.0, Int_t save_opt = 0) 
{
	TTree* tree;
	if(mactolay[mac]=='i') tree = fTreeInner;
	else if(mactolay[mac]=='o') tree = fTreeOuter;
	else {
		cout << "error in PlotSingleChanCut: unknown runtype provided - " << mactolay[mac] << endl;
		return 0;
	}

    if(!IsChanConfigured(mac,chan)) {
        cout << "ERROR in PlotSingleChanCut: channel not configured (if not mac)" << endl;
        return nullptr;
    }

	//readecalibration files each array[20]
	const Double_t* const conv = GetGain(mac);
	const Double_t* const ped = GetPed(mac);
	if(conv==nullptr) {
		cout << "ERROR in PlotSingleChanCut: could not retreive gain values!" << endl;
		return nullptr;
		}
	if(ped==nullptr) {
		cout << "ERROR in PlotSingleChanCut: could not retreive ped values!" << endl;
		return nullptr;
	}

	Int_t nentries = tree->GetEntriesFast(); //get number of entries in tree

	//Set histogram low/high values
	Float_t low = -1.5;
	Float_t high = 60.5;
	Float_t binning = 1.0; 

	//set binning: divide by desired # p.e.u. / bin
	Int_t bins = ( high - low ) / binning;

	TString ytitle = "counts / " + to_string(binning).substr(0,4) + " p.e.u.";
	TString htitle = "LY@4m: mac5 "+FormatMacString(mac)+", ch. "+to_string(chan);
		
	TH1F *h = new TH1F("h", htitle, bins, low, high);
	h->GetXaxis()->SetTitle("sigal amplitude [p.e.u.]");
	h->GetYaxis()->SetTitle(ytitle);
	h->GetYaxis()->SetTitleOffset(1.2);	
	h->SetLineWidth(3);

	//loop over events in tree (applying single chan cut)
	for (Int_t ientry=0; ientry < nentries; ientry++ )
	{
		tree->GetEntry(ientry); //get event number ientry
		if(mac5!=mac) //examin only events with given mac5 
			continue;
		bool pass = true;
		//check all but given channel are below specified threshold
		for ( auto ch = febchans[mac].begin(); ch!=febchans[mac].end(); ch++ ) 
			if(*ch!=chan && (chg[*ch]-ped[*ch])/conv[*ch]>cut)
				pass=false;
		if(pass) h->Fill((chg[chan]-ped[chan])/conv[chan]);
	}//end loop over events
		
	//text box to display fit max, FWHM
	TPaveText *stat2 = new TPaveText(0.79,0.19,0.89,0.34,"NDC");
	stat2->SetName("mystat");
	stat2->SetBorderSize(0);
	stat2->SetFillColor(kWhite);
	stat2->SetTextSize(0.035);
	stat2->SetTextAlign(31);

	TCanvas *cfit = new TCanvas();
	cfit->SetGrid();

	//Set y-range to include error bars (assuming Poisson statistics)
	int max_bin = GetMaxBinInRange(h,1.5,high);
	double max_val = GetMaxInRange(h,1.5,high);
	int max_err = h->GetBinError(max_bin);
	double range_max = 1.01*(max_val + max_err);
	h->GetYaxis()->SetRangeUser(0,range_max);

	gStyle->Reset("Modern");
	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(111);
	gStyle->SetStatW(0.25);
	gStyle->SetStatH(0.4);
	gStyle->SetStatX(0.95);
	gStyle->SetStatY(0.95);

	h->Draw("E0hist"); //draw error bars assuming Poisson statistics
	cfit->Update(); //this makes all the difference in avoiding seg. viol.!!!

	const Double_t* const lanarr = LangausFit(h,save_opt); //size 10
	if(lanarr[1]<1.0) cout << "bad fit!" << endl;

	TString line = "fit max [PE]: "+to_string(lanarr[0]).substr(0,4); 
	stat2->AddText(line);
	line = "fit FWHM [PE]: "+to_string(lanarr[1]).substr(0,4);
	stat2->AddText(line);
	stat2->Draw();
	
	gPad->Modified();
	gPad->Update();
	cfit->Update();//this makes all the difference in avoiding seg. viol.!!!

	if ( save_opt == 1 )
	{
		//setup to save histograms to file
		int ctr=1;
		TString suff=".png";
		TString fname = lyPlotDir + "/";
		fname += "ly_mac"+FormatMacString(mac)+"_ch"+to_string(chan);
		fname+= "_thr_"+to_string(cut).substr(0,3)+"pe_1"+suff;

		while (!gSystem->AccessPathName(fname))
		{
			fname.Remove(fname.Length()-suff.Length(),fname.Length());
			fname.Remove(fname.Length()-1,fname.Length());
			ctr++;
			fname+=to_string(ctr)+suff;
		}

		//write image to file
		TImage *img = TImage::Create();
		img->FromPad(cfit);
		img->WriteImage(fname);

		delete img;
		delete cfit;
		delete stat2;

		TString fana = anaDir+"/data_lyfit_thr"+to_string(cut).substr(0,3)+"_mac";
		fana+=FormatMacString(mac)+"_ch"+to_string(chan)+".txt";
		if (gSystem->AccessPathName(fana)){

			ofstream fout;
			fout.open(fana);
			fout << h->Integral(h->FindBin(5.5),h->FindBin(38.5));
			fout << " ";
			for (Int_t k=0; k<10; k++) fout << lanarr[k] << " "; fout << "\n";
			fout.close();
		}
		
		else
			cout << fana << " already exists -> new fit data not written" << endl;

		delete[] lanarr;
		delete h;
		return nullptr;
	}//end save opt

	return h;

}//end def.

//*****************************************************************************************
//
//****************************************************************************************
void mppc::ThresholdStats(char runtype='h', Float_t thresh = 0.5)
{
	TTree* tree;
	if(runtype=='h') tree = fTreeInner;
	else if(runtype=='s') tree = fTreeOuter;
	else {
		cout << "error in ThresholdStats: unknown runtype provided - " << runtype << endl;
		return;
	}
	//read calibration files
	Double_t * conv85 = GetGain(85);
	Double_t * ped85 = GetPed(85);
	Double_t * conv213 = GetGain(213);
	Double_t * ped213 = GetPed(213);

	Bool_t hit[21];

	Int_t nentries = tree->GetEntriesFast(); //get number of entries in tree

	Float_t val=0.0;
	//Float_t thresh = 0.0;
	Int_t nabove = 0;
	TH1F *h85 = new TH1F("h85","above threshold multiplicity",20,0,20);
	TH1F *h213 = new TH1F("h213","above threshold multiplicity",20,0,20);
	TH2F *h2d85 = new TH2F("h2d85","south FEB above thresh mult by strip",20,1,21,20,0,20);
	TH2F *h2d213 = new TH2F("h2d213","north FEB above thresh mult by strip",20,1,21,20,0,20);

	for ( Int_t i=0; i<nentries; i++ )
	{
		tree->GetEntry(i);
		nabove = 0;
		for (Int_t j=0; j<21; j++) hit[j] = kFALSE;
		for ( Int_t chan=1; chan<21; chan++ )
		{
			if (mac5 == 85) val = (chg[chan]-ped85[chan])/conv85[chan];
			if (mac5== 213) val = (chg[chan]-ped213[chan])/conv213[chan];
			if (val>thresh){nabove+=1; hit[chan]= kTRUE;}

		}//end loop over channels

		if (mac5==85) h85->Fill(nabove);
		if (mac5==213) h213->Fill(nabove);

		for (Int_t chan=1; chan<21; chan++)
		{
			if (hit[chan])
			{			
				if (mac5==85) h2d85->Fill(chan,nabove);
				if (mac5==213) h2d213->Fill(chan,nabove);
			}
			else if (nabove==0)
			{			
				if(mac5==85) h2d85->Fill(chan,0);
				if(mac5==213) h2d213->Fill(chan,0);
			}
		}

	}//end loop over events

	h85->SetLineColor(kBlue);
	h213->SetLineColor(kRed);
	h85->GetXaxis()->SetTitle("# channels");
	h213->GetXaxis()->SetTitle("# channels");
	h2d85->GetYaxis()->SetTitle("# channels");
	h2d213->GetYaxis()->SetTitle("# channels");
	h2d85->GetXaxis()->SetTitle("channel");
	h2d213->GetXaxis()->SetTitle("channel");

	//if (hod_pos==0||hod_pos==1) {h85->Draw(); h213->Draw("same");}
	//else {h213->Draw(); h85->Draw("same");}

	new TCanvas(); h2d85->Draw("COLZ");
	new TCanvas(); h2d213->Draw("COLZ");

}

void mppc::ThresholdCorrelation(char runtype='h', Int_t ch=31, Int_t mac=85, Float_t thresh=4.5)
{
    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
	cout << "error in ThresholdCorrelation: unknown runtype provided - " << runtype << endl;
	return; 
    } 

	//read calibration files
	Double_t * conv85 = GetGain(85);
	Double_t * ped85 = GetPed(85);
	Double_t * conv213 = GetGain(213);
	Double_t * ped213 = GetPed(213);

	Float_t val=0.0;//,val2;
	Int_t nabove=0;
	Bool_t hit[21];

	Int_t nentries = tree->GetEntriesFast(); //get number of entries in tree

	TString title = "events with ch. "; title+=ch; title+=" > "; title+=thresh; title+=" p.e.";
	if (mac==85) title+=" - south FEB";
	if (mac==213) title+=" - north FEB";

	TH2F *h = new TH2F("h",title,20,1,21,19,2,21);

	for (Int_t ientry=0; ientry<nentries; ientry++)
	{
		tree->GetEntry(ientry);

		for(Int_t j=0; j<21; j++) hit[j] = kFALSE;
		nabove=0;

		for(Int_t chan=1; chan<21; chan++)
		{
			if(mac5==mac&&mac==85 )  val = (chg[chan]-ped85[chan])/conv85[chan];
			if(mac5==mac&&mac==213) val = (chg[chan]-ped213[chan])/conv213[chan];
			if (val>thresh) {nabove++; hit[chan]=kTRUE;}
		}

		if(hit[ch]) {for(Int_t chan=1; chan<21; chan++) {if (hit[chan]) h->Fill(chan,nabove); else h->Fill(chan,0);}}
	}

	h->GetXaxis()->SetTitle("strip");
	h->GetYaxis()->SetTitle("number of strips above threshold");
	h->SetStats(0);

	TCanvas *c = new TCanvas();
	c->SetLogz();
	h->Draw("COLZ");

}
//**************************************************************************************
//**************************************************************************************
//******		   END ANALYSIS/PLOT FUNCTION DEFINITIONS		 ******
//**************************************************************************************
//**************************************************************************************

//**************************************************************************************
//**************************************************************************************
//******		     BEGIN Fit FUNCTION DEFINITIONS		       ******
//**************************************************************************************
//**************************************************************************************

//********************************************************************
//** Description: function to perform convolution of landau w/gauss **
//** the function is auxilliary to the langaus fit		  **
//********************************************************************
Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

	// Numeric constants
	Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
	Double_t mpshift  = -0.22278298;       // Landau maximum location

	// Control constants
	Double_t np = 200.0;      // number of convolution steps
	Double_t sc =   2.0;      // convolution extends to +-sc Gaussian sigmas

	// Variables
	Double_t xx;
	Double_t mpc;
	Double_t fland;
	Double_t sum = 0.0;
	Double_t xlow,xupp,step;
	Double_t sig;

	// MP shift correction
	mpc = par[1] - mpshift * par[0];
	//sig = TMath::Sqrt(x[0]);//par[3]*x[0]+par[4]);

	// Range of convolution integral
	xlow = x[0] - sc * par[3];
	xupp = x[0] + sc * par[3];

	step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
	for(Int_t i=1; i<=np/2; i++) 
	{
		//sig = par[3]*x[0]+par[4];

		xx = xlow + (i-.5) * step;
		fland = TMath::Landau(xx,mpc,par[0]);// / par[0];
		sum += fland * TMath::Gaus(x[0],xx,par[3]);
		//sum += fland * TMath::Gaus(x[0],xx,sig);
	
		xx = xupp - (i-.5) * step;
		fland = TMath::Landau(xx,mpc,par[0]);// / par[0];
		sum += fland * TMath::Gaus(x[0],xx,par[3]);
		//sum += fland * TMath::Gaus(x[0],xx,sig);
      }

      return (par[2] * sum); //step * sum * invsq2pi / par[3]);

}//end function definition

//***************************************************************************
//** Description: convoluted Gaussian-Landau fit function.		 **
//** Fit parameters:						       **
//**  par[0] = Width (scale) parameter of Landau density		   **
//**  par[1] = Most Probable (MP, location) parameter of Landau density    **
//**  par[2] = Total area (integral -inf to inf, normalization constant)   **
//**  par[3]= Width (sigma) of convoluted Gaussian function		**
//**								       **
//** Variables for langaufit call:					 **
//**  his	     histogram to fit				     **
//**  fitrange[2]     lo and hi boundaries of fit range		    **
//**  startvalues[4]  reasonable start values for the fit		  **
//**  parlimitslo[4]  lower parameter limits			       **
//**  parlimitshi[4]  upper parameter limits			       **
//**  fitparams[4]    returns the final fit parameters		     **
//**  fiterrors[4]    returns the final fit errors			 **
//**  ChiSqr	  returns the chi square			       **
//**  NDF	     returns ndf					  **
//***************************************************************************
TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");

   for (i=0; i<4; i++) {
   //for (i=0; i<5; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

  // ffit->FixParameter(3,0.01);

   his->Fit(FunName,"R0MQL");   // fit within specified range, use ParLimits, do not plot
   Float_t chi2 = ffit->GetChisquare();
   Int_t ndf = ffit->GetNDF();
   Float_t chi2_reduced = chi2/ndf;


   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
   //for (i=0; i<5; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();	   // obtain ndf

   return (ffit);	      // return fit function

}//end function definition

//*********************************************************************
//**  Seaches for the location (x value) at the maximum of the       **
//**  Landau-Gaussian convolution and its full width at half-maximum **
//*********************************************************************
Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 100000;


   // Search for maximum

   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = langaufun(&x,params);

      if (l < lold)
	 step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-1);

   maxx = x;

   fy = l/2;


   // Search for right x location of fy

   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
	 step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-2);

   fxr = x;


   // Search for left x location of fy

   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
	 step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-3);


   fxl = x;

   FWHM = fxr - fxl;
   return (0);
}//end function definition

//*************************************************************************
//** Description: function to fit convoluted Gaussian-Landau function to **
//** user-specified histogram (TH1F* hfit) and channel (Int_t chan)      **
//** The generated image is written to file in current working directory **
//**								     **
//** Fit parameters:						     **
//**  par[0] = Width (scale) parameter of Landau density		 **
//**  par[1] = Most Probable (MP, location) parameter of Landau density  **
//**  par[2] = Total area (integral -inf to inf, normalization constant) **
//**  par[3]= Width (sigma) of convoluted Gaussian function	      **
//**								     **
//**  SNRPeak = apparent fit peak position				 **
//**  SNRFWHM = apparent fit FWHM					**
//*************************************************************************
Double_t* mppc::LangausFit(TH1F *hfit, Int_t save_opt=0)
{
	Double_t *statarr = new Double_t[10];

	Float_t minpe = 5.5, maxpe = 39.5;

	//Float_t hist_max = hfit->GetBinLowEdge(hfit->GetMaximumBin())+0.5;
	Float_t hist_max = hfit->GetBinLowEdge(GetMaxBinInRange(hfit,minpe,maxpe))+0.5;
	if (hist_max<minpe||hist_max>maxpe) cout << "hist_max out of range! = " << hist_max << endl;
	Float_t hist_width = TMath::Sqrt(hist_max);

	//declare vars for paramater limits/initialization
	Double_t fr[2];
	Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];

	//lower, upper bound for fit
	//if(hist_max-1.75*hist_width > minpe) fr[0] = hist_max-1.75*hist_width;
	//else 
	fr[0] = minpe;
	//if(hist_max+2.25*hist_width < maxpe) fr[1]=hist_max+2.25*hist_width;
	//else 
	fr[1] = maxpe;

	pllo[0]=0.001; plhi[0]=8*hfit->GetRMS(); sv[0]=hfit->GetRMS(); //low,high,IV Landau width

	//low,high,IV Landau MPV
	//if(hist_max-hist_width > fr[0]) pllo[1]=hist_max-hist_width;
	//else 
	pllo[1]=fr[0];
	//if(hist_max+hist_width < fr[1]) plhi[1]=hist_max+hist_width;
	//else 
	plhi[1]=fr[1]; 
	sv[1]=hist_max; 

	pllo[2]=1.0; plhi[2]=5*hfit->Integral(); sv[2]=hfit->Integral(); //low,high,IV normalization constant
	pllo[3]=1.0; plhi[3]=2*TMath::Sqrt(sv[1]); sv[3]= TMath::Sqrt(sv[1]);//hist_width; //low,high,IV Gauss sigma

	Double_t chisqr, chisqr0;
	Int_t    ndf;
	TF1 *fitsnr = langaufit(hfit,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf); //`auxiallry fit function (defined above)

	Double_t SNRPeak, SNRFWHM, sv1i;
	langaupro(fp,SNRPeak,SNRFWHM);
	Int_t ctr=-1;
	Float_t dsv[12] = {-0.5,0.5,-1.0,1.0,-1.5,1.5,-2.0,2.0,-2.5,2.5,-3.0,3.0};

	if (SNRFWHM < 1.0 || chisqr/ndf > 5.0)
	{
		if (SNRFWHM > 1.0) chisqr0 = chisqr;
		else chisqr0 = 100.0;
		sv1i = sv[1];

		for(Int_t i=0; i<12; i++)
		{
			if ((i%2==0&&sv1i+dsv[i]>fr[0]+0.5)||(i%2!=0&&sv1i+dsv[i]<fr[1]-0.5))
			{
				sv[1] = sv1i+dsv[i];

				if(sv[1]-0.5*hist_width > fr[0]) pllo[1]=sv[1]-0.5*hist_width;
				else pllo[1]=fr[0];

				if(sv[1]+0.5*hist_width < fr[1]) plhi[1]=sv[1]+0.5*hist_width;
				else plhi[1]=fr[1];

				fitsnr = langaufit(hfit,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
				langaupro(fp,SNRPeak,SNRFWHM);
				if (chisqr<chisqr0&&SNRFWHM>1.0) {chisqr0=chisqr; ctr=i;}
			}
		}
		
		if (ctr==-1) sv[1] = sv1i;
		else sv[1] = sv1i+dsv[ctr];

		if(sv[1]-0.5*hist_width > fr[0]) pllo[1]=sv[1]-0.5*hist_width;
		else pllo[1]=fr[0];

		if(sv[1]+0.5*hist_width < fr[1]) plhi[1]=sv[1]+0.5*hist_width;
		else plhi[1]=fr[1];

		fitsnr = langaufit(hfit,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
		langaupro(fp,SNRPeak,SNRFWHM);
	}

	TString peak; peak += SNRPeak;
	TString fwhm; fwhm += SNRFWHM;
	peak.Remove(4); fwhm.Remove(4);

	hfit->Draw(); //draw original histogram
	hfit->Draw("E0,SAME"); //draw error bars
	fitsnr->Draw("l,same"); //draw fit curve over hist.
	
	statarr[0] = SNRPeak; statarr[1] = SNRFWHM;                                //fit peak, FWHM
	statarr[2] = fitsnr->GetParameter(0); statarr[3] = fitsnr->GetParError(0); //landau width
	statarr[4] = fitsnr->GetParameter(1); statarr[5] = fitsnr->GetParError(1); //landau MPV
	statarr[6] = fitsnr->GetParameter(2); statarr[7] = fitsnr->GetParError(2); //norm const
	statarr[8] = fitsnr->GetParameter(3); statarr[9] = fitsnr->GetParError(3); //gauss sigma

	return statarr;

}//end def.

//**************************************************************************************
//**************************************************************************************
//******		     END FIT FUNCTION DEFINITIONS			 ******
//**************************************************************************************
//**************************************************************************************

//**************************************************************************************
//**************************************************************************************
//******		BEGIN CALIBRATION FUNCTION DEFINITIONS		    ******
//**************************************************************************************
//**************************************************************************************
//**************************************************************************************
void mppc::Cal(Int_t mac)
{
	if(!IsMacConfigured(mac)){
		cout << "error in Cal: not configured to accept mac5 = " << mac << endl;
		return;
	}

	Double_t* garr = new Double_t[6];
	Double_t* parr = new Double_t[6];
	bool hasPed=false, hasGain=false;

	//setup files for writing
	ofstream fgout, fpout;
	TString fgname=calDir+"/gains", fpname= calDir+"/peds";

	fgname+=FormatMacString(mac)+".txt";
	fpname+=FormatMacString(mac)+".txt";

	//check if files already exist
	if(!gSystem->AccessPathName(fpname)) {
		cout << "skipping pedestal calibration for mac5 " << FormatMacString(mac) 
			<< "! Existing files found." << endl;
	}
	else
	{
		fpout.open(fpname);

		//loop over FEB channels, perform pedestal and gain fits, save to file
		for (int i = 0; i < 32; i++)
		{
			//perform pedestal calibration for all channels
			parr = PlotPedFit(mac,i,1,false);
			// mean, mean err, sigma, sigma err, const, const err
			fpout << parr[2] << " " << parr[3] << " " << parr[4] << " " 
				<< parr[5] << " " << parr[0] << " " << parr[1] <<'\n';
		}
		fpout.close();
		cout << fpname << " written" << endl;
	}

	if(!gSystem->AccessPathName(fgname))
		cout << "skipping calibration for mac5 " << FormatMacString(mac) 
			<< "! Existing files found." << endl;
	else {
		fgout.open(fgname);

		for (int i = 0; i < 32; i++)
		{
			//perform gain calibration only for configured channels
			if (!IsChanConfigured(mac,i)){
				fgout << 0 << '\n';
			}
			else {
				garr = PlotGainFit(mac,i,1,55.5);
				// gain, gain err, ped, ped err, chi-2, ndf
				fgout << garr[0] << " " << garr[1] << " " << garr[2] << " " 
					<< garr[3] << " " << garr[4] << " " << garr[5] << '\n';
			}
		}//loop over FEB channels

		fgout.close();
		cout << fgname << " written" << endl;

	}//end else cal gain

	delete[] garr;
	delete[] parr;

}//end def. Cal

//**********************************************************************************
//** Function to plot fitted calibration spectrum - used to extract gain values   **
//** Int_t chan -> channel to be plotted					  **
//** Int_t low, high -> lower / upper bounds of fit range			 **
//** Float_t binning -> # ADC / bin					       **
//** Float_t save_opt = 0 (no save) or = 1 (save to working directy)	      **
//**
//** edit selection cut applied during histo filling as needed
//**********************************************************************************
Double_t* mppc::PlotGainFit(Int_t mac=80, Int_t chan=30, Int_t save_opt = 1, Double_t gain_seed=-1.0)
{
	///make sure we expected given mac5, channel
	if(!IsChanConfigured(mac,chan)) {
		cout << "error in PlotGainFit: not configured for channel passed: mac5 " << mac << ", chan " << chan << endl;
		return nullptr;
	}

        TTree* tree;
        if(mactolay[mac]=='i') tree = fTreeInner;
        else if(mactolay[mac]=='o') tree = fTreeOuter;
	const int nentries = tree->GetEntriesFast();

        Double_t *peds =GetPed(mac);

	Int_t lowbin, highbin;
	if (gain_seed==-1.0) {
		lowbin=400;
		highbin=700;
	}
	else{
		lowbin=peds[chan]+4.0*gain_seed;
		highbin=peds[chan]+8.0*gain_seed;
	}

	//histogram setup
	const Int_t binning = 7;
	const Int_t low = lowbin;//400;
	const Int_t high = highbin;//700;	
	const Int_t bins = (high-low)/binning;

	TString htitle = "mac5 "+FormatMacString(mac)+", ch. "+to_string(chan);
	TString ytitle = "counts / "+to_string(binning)+" ADC";

	TH1F *h = new TH1F("h", htitle, bins, low, high);
	h->GetXaxis()->SetTitle("signal amplitude [ADC]");
	h->SetLineWidth(3);
	h->GetYaxis()->SetTitleOffset(1.2);

	Double_t p1;
	Double_t * statarr = new Double_t[6];

	//fill histogram with charge values for given mac5, channel
	for (int ientry=0; ientry < nentries; ientry++ )
	{
		tree->GetEntry(ientry);
		if ( mac5 == mac ) h->Fill(chg[chan]);
	}

	gStyle->Reset("Modern");
	TCanvas *c = new TCanvas();
	c->SetGrid();

	//stat box dimensions
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.95);
	gStyle->SetStatH(0.05);
	gStyle->SetStatW(0.18);

	h->Draw();

	//fit gain and retain fit values, stats
	statarr = FitGain(mac, h,chan,save_opt,gain_seed);

	//show fit results in stats box
	gStyle->SetOptStat(11);
	gStyle->SetOptFit(0);

	c->Update();

	if ( save_opt == 1 )
	{
		int ctr=1;
		TString suff = ".png";
		TString fname = calPlotDir+"/";
		fname+=FormatMacString(mac)+"_ch";
		if(chan<10) fname+="0";
		fname+=to_string(chan)+"_fit-spec";
		if(gain_seed>0) fname+="_gain-seeded";
		fname+="_1"+suff;

		while (!gSystem->AccessPathName(fname)&&ctr<10)
		{
			fname.Remove(fname.Length()-suff.Length(),fname.Length());
			fname.Remove(fname.Length()-1,fname.Length());
			ctr++;
			fname+=to_string(ctr)+suff;
		}

		//write image to file
		TImage *img = TImage::Create();
		img->FromPad(c);
		img->WriteImage(fname);

		//deallocate memory
		delete img;
		delete h;
		delete c;
	}

	//if(gain_seed<0)
	//	delete h;

	return statarr; 

}//end function definition

//*************************************************************************************
//** Fit function which takes a pointer to a max spectrrum to be fit as input param. **
//** Code copied from Igor's macro FEBDAQMULT.C					     **
//*************************************************************************************
Double_t* mppc::FitGain(Int_t mac, TH1F *hs, Int_t chan = 24, Int_t save_opt = 0, Double_t gain_seed=-2.0)
{
	TSpectrum *s = new TSpectrum();
	Int_t npeaks = s->Search(hs,1.75,"",0.18);//args=source histo, sigma of searched peaks, options, threshold
	//cout << "TSpectrum found " << npeaks << " peaks" << endl;
	Int_t ctr = 2;

	const size_t size = npeaks;
	Double_t peak_chi2[size];
	Double_t peak_ndf[size];
	Double_t peak_sigma[size];
	Double_t peak_sigma_err[size];
	Double_t * statarr = new Double_t[6];

	Double_t gain = gain_seed;
	Double_t *peds =GetPed(mac);
	Double_t *pwidths=GetPedWidth(mac);
	Double_t *gains;//=InitGain(mac);
	if(gain_seed<0){
		gain=gains[chan];
		gains=InitGain(mac);
	}
	Double_t *peaks = s->GetPositionX(); //candidate photopeaks ADC position

	//cout << "using ped mean(from fit): " << peds[chan] << " ADC" << endl;
	//cout << "using ped RMS(from fit): " << pwidths[chan] << " ADC" << endl;
	//cout << "using gain(guess): " << gains[chan] << " ADC/pe" << endl;

	//Sort peaks
	double temp;
	int nchanges;
	do
	{
		nchanges=0;
		for(int p=0;p<npeaks-1;p++)
			if(peaks[p]>peaks[p+1])
			{
				temp=peaks[p];
				peaks[p]=peaks[p+1];
				peaks[p+1]=temp;
				nchanges++;
			}
   	} while( nchanges != 0 );

	//ascending list of peak number(x) vs. ADC value(y) from TSpectrum
	Double_t x[10];
	Double_t y[10], ey[10];

	//same list with "bad" peaks excluded (what eventually goes into the gain fit)
	Double_t gx[11];
	Double_t gy[11], gey[11];
	int peak_offset = round((peaks[0]-peds[chan])/gain); //estimate peak number

	for (int j=0; j<npeaks&&j<10; j++) 
	{
		x[j] = j+peak_offset;
		y[j] = peaks[j];
	}

	int gg = 1; //index of passing peak array
	int nplow = 0; //no. peaks close to low hist edge
	int nphigh = 0; //no. peaks close to high hist edge

	//first peak in passing array is pedestal
	gx[0] = 0;
	gy[0] = peds[chan];
	gey[0] = pwidths[chan];

	//first and refit chi-squareds
	double chisqr, chisqr0;
	int rnum = 0; //refit number

	//fit peaks about TSpectrum values
	for (int g=0 ; g<npeaks; g++)
	{
		//initial gaus fit to peak from TSpectrum
		TF1 *gfit = new TF1("gfit","gaus",y[g]-20,y[g]+20);
		gfit->SetParameter(0,hs->GetBinContent(hs->FindBin(y[g])));//200);
		gfit->SetParameter(1,y[g]);
		gfit->SetParameter(2,12);
		gfit->SetParLimits(0,0,20000);
		gfit->SetParLimits(1,y[g]-15,y[g]+15);
		gfit->SetParLimits(2,8,40);
		//cout << "fitting gaussian to peak detected at " << y[g] << " ADC..." << endl;

		hs->Fit(gfit,"MQLR"); //fit peak
		chisqr0 = gfit->GetChisquare()/gfit->GetNDF(); //get reduced chi-square from fit

		//if chisqr is "large", try with smaller range
		/*if (chisqr0>1.9) 		
		{
			//cout << "X^2 too large (" << chisqr0 << "). refitting..." << endl;
			gfit->SetParameter(1,y[g]);            //reinitialize mean to TSpectrum peak
			gfit->SetParLimits(1,y[g]-10,y[g]+10); //reduce mean range by 10
			gfit->SetRange(y[g]-15,y[g]+15);       //reduce fit range by 10

			hs->Fit(gfit,"MQLR");                  //fit again
			chisqr = gfit->GetChisquare()/gfit->GetNDF();  //get chi-square

			//try another range
			gfit->SetParameter(1,y[g]);            //reinitialize mean to TSpectrum peak
			gfit->SetParLimits(1,y[g]-5,y[g]+5);   //reduce mean range by 10
			gfit->SetRange(y[g]-10,y[g]+10);       //reduce fit range by 10
			hs->Fit(gfit,"MQLR");                  //fit again

			//check if second refit has better chi-square than first refit
			if(gfit->GetChisquare()/gfit->GetNDF()<chisqr) {
				chisqr = gfit->GetChisquare()/gfit->GetNDF(); 
				rnum = 2; //flag fit result
			}
			//if refits produced better chi-square
			if(chisqr<chisqr0)
			{
				//cout << "new X^2 better (" << chisqr << "). keeping..." << endl;
				//if second refit best
				if (rnum==2){
					gfit->SetParameter(1,y[g]); 
					gfit->SetParLimits(1,y[g]-10,y[g]+10);
					gfit->SetRange(y[g]-15,y[g]+15);
				}
				//first refit best
				else {
					gfit->SetParameter(1,y[g]); 
					gfit->SetParLimits(1,y[g]-5,y[g]+5);
					gfit->SetRange(y[g]-10,y[g]+10);
				}
			}//end if lower chi-sqr
			else
			{
				gfit->SetParameter(1,y[g]);
				gfit->SetParLimits(1,y[g]-20,y[g]+20);
				gfit->SetParLimits(1,y[g]-15,y[g]+15);
			}//end else
			//fit and add results to fit list
			hs->Fit(gfit,"MQLR+");
			gfit->Draw("same");
		}//end if initial chi-square good

		//we are happy with this fit, add to fit list
		else {*/
			hs->Fit(gfit,"MQLR+");
			gfit->Draw("same");
		//}
		//check for peaks near the edges of the histogam
		if( y[g]<hs->GetBinLowEdge(1)+15) nplow++;
		if( y[g]>hs->GetBinLowEdge(hs->GetNbinsX())-15) nphigh++;

		//ignore edge peaks and peaks with fit mean > 15 from peak
		if( y[g]>hs->GetBinLowEdge(1)+15 && y[g]<hs->GetBinLowEdge(hs->GetNbinsX())-15&&abs(y[g] - gfit->GetParameter(1)) < 15) 
		{
			//skip false peaks (may need improvement - currently relies on init values)
			if(g!=0 && g!=npeaks-1 //check it's not first or last peaks
			 && (y[g+1]-y[g]<gain*0.7 || y[g]-y[g-1]<gain*0.7))
			{ //check peak within 30% of expected gain w.r.t adj.
				//cout << "determined peak " << g << " (" << y[g] << " ADC) is false peak. removing..." << endl;
				for (Int_t j=g; j<npeaks; j++) x[j] = x[j]-1; 
			}//overwrite current value, shifting higher values down by 1 index
			//if not missed peak, could there have been a skipped peak?
			else
			{
				if(g!=0&&g!=npeaks-1&&y[g+1]-y[g]<gain*1.2&&y[g]-y[g-1]>gain*1.5) //missed peak adjust
				{
					for (Int_t j=g; j<npeaks; j++) x[j] = x[j]+1;
					//cout << "missed peak detected. shifting peaks " << g << " and higher up by 1..." << endl;
				}
				//if last peak likely occuring after skipped peak
				if(g==npeaks-1 && (y[g]-y[g-1])/gain >1.5) {
					x[g]+=(int)((y[g]-y[g-1])/gain);
					//cout << "missed peak(s) detected before last peak...shifting last peak" << endl;
				}
						
				peak_chi2[gg] = gfit->GetChisquare();
				peak_ndf[gg] = gfit->GetNDF();
				peak_sigma[gg] =  gfit->GetParameter(2);
				peak_sigma_err[gg] = gfit->GetParError(2);

				gx[gg] = x[g];
				gy[gg] = gfit->GetParameter(1);
				gey[gg] = sqrt(gx[gg]+gey[0]*gey[0]);//gfit->GetParError(1);
				gg++;
			}//end else not missed but was there skip
		}//end if not edge peak, not far from TSpectrum value
	}//end for peaks

	//if(nplow>0) cout << "found and removed"  << nplow << " low edge peaks" << endl;
	//if(nphigh>0) cout << "found and removed"  << nplow << " high edge peaks" << endl;

	//if any low edge peaks, shift peaks to left in peak number
	if (nplow>0) for (Int_t i=1; i<gg; i++) gx[i] = gx[i]-(nplow-1);

	//graph of adc(y) vs. photo-peak number (x)
	TGraphErrors* gr_mean = new TGraphErrors(gg,gx,gy,0,gey);

	//linear fit function
	TF1 *fit = new TF1("fit","[0] + [1]*x",gx[0]-0.25,gx[gg-1]+0.25);

	//name, initialize gain fit parameters
	fit->SetParName(1,"Gain");
	fit->SetParName(0, "Pedestal");
	fit->SetParameter(1,gain);
	fit->SetParameter(0,peds[chan]);
	fit->SetParLimits(1,gain-20,gain+20);
	fit->SetParLimits(0,peds[chan]*0.8,peds[chan]*1.2);

	//perform gain fit
	gr_mean->Fit(fit, "QR");

	//check if gain fit is bad according to chi-square
	if (fit->GetChisquare()/fit->GetNDF()>5.0)
	{
		//cout << "gain fit X^2 too large...shifting all peaks by 1" << endl;
		chisqr=fit->GetChisquare();
		for(int i=1; i<gg; i++) gx[i]+=1;
		gr_mean = new TGraphErrors(gg,gx,gy,0,gey);
		fit->SetRange(gx[0]-0.25,gx[gg-1]+0.25);
		gr_mean->Fit(fit, "QR");
		if (fit->GetChisquare()<chisqr) chisqr=fit->GetChisquare();
		else
		{
			for(int i=1; i<gg; i++) gx[i]-=2;
			gr_mean = new TGraphErrors(gg,gx,gy,0,gey);
			fit->SetRange(gx[0]-0.25,gx[gg-1]+0.25);
			gr_mean->Fit(fit, "QR");
			if (fit->GetChisquare()<chisqr) chisqr=fit->GetChisquare();
			else
			{
				for(int i=1; i<gg; i++) gx[i]+=1;
				gr_mean = new TGraphErrors(gg,gx,gy,0,gey);
				fit->SetRange(gx[0]-0.25,gx[gg-1]+0.25);
				gr_mean->Fit(fit, "QR");
			}
		}
	}



	TString grtitle = "Mac5 "+FormatMacString(mac)+", Ch. "+to_string(chan)+" Gain Fit";

	gr_mean->SetTitle(grtitle);
	gr_mean->GetXaxis()->SetTitle("Peak #");
	gr_mean->GetYaxis()->SetTitle("ADC value");
	gr_mean->SetMarkerColor(4);
	gr_mean->SetMarkerStyle(20);
	gr_mean->SetFillColor(0);
	gr_mean->GetXaxis()->SetRangeUser(gx[0]-0.5,gx[gg-1]+0.5);

	//if (fit->GetChisquare()/fit->GetNDF()>30.0) cout << "Poor X^2! mac5 " << mac << "ch. " << chan << endl;

	gStyle->SetOptStat(0100);
	gStyle->SetOptFit(1111);

	TCanvas *c2 = new TCanvas();
	c2->cd();
	c2->SetGrid();

	gStyle->SetStatX(0.5);
	gStyle->SetStatY(0.9);
	gStyle->SetStatH(0.15);
	gStyle->SetStatW(0.2);

	gr_mean->Draw("ALP");
	fit->Draw("L,SAME");
	
	statarr[0] = fit->GetParameter(1); //gain
	statarr[1] = fit->GetParError(1);  //gain error
	statarr[2] = fit->GetParameter(0); //pedestal mean
	statarr[3] = fit->GetParError(0);  //pedestal mean error
	statarr[4] = fit->GetChisquare();  //X^2
	statarr[5] = fit->GetNDF();        //NDF

	Float_t chi2 = fit->GetChisquare();
	Int_t ndf = fit->GetNDF();
	Float_t chi2_reduced = chi2/ndf;

	if ( save_opt == 1 )
	{
		int ctr=1;
		TString suff = ".png";
		TString fname=calPlotDir+"/"; 
		fname+=FormatMacString(mac)+"_ch";
		if(chan<10) fname+="0";
		fname+=to_string(chan)+"_fit-gain";
		if(gain_seed>0) fname+="_gain-seeded";
		fname+="_1"+suff;
		while (!gSystem->AccessPathName(fname))
		{
			fname.Remove(fname.Length()-suff.Length(),fname.Length());
			fname.Remove(fname.Length()-1,fname.Length());
			ctr++;
			fname+=to_string(ctr)+suff;
		}

		//write image to file
		TImage *img = TImage::Create();
		img->FromPad(c2);
		img->WriteImage(fname);

		//deallocate memory
		delete img;
		delete c2;
		delete gr_mean;
		delete fit;
		delete s;
	}

	return statarr;

}//end function definition

//************************ generate plot of fit pedestal and return ped. mean value ******************
Double_t* mppc::PlotPedFit(Int_t mac=80, Int_t chan=30, Int_t save_opt=0, bool useInit=false)
{   
	//First, create histogram (a slight modification of GetHist, allows user to set min/max range)
	TTree* tree=0;
	if(mactolay[mac]=='i') tree = fTreeInner;
	else if(mactolay[mac]=='o') tree=fTreeOuter;
	const int nentries = tree->GetEntriesFast();

	Double_t* peds = new Double_t[32];
	Double_t* statarr = new Double_t[6];

	Int_t low=1, high=300, binning=2;
	if(useInit){ 
		peds = InitPed(mac);
		low = (Int_t)peds[chan]-20;
		high = (Int_t)peds[chan] + 20;
	}
	Int_t bins = (high-low)/binning;

	TString htitle = "mac5 "+to_string(mac)+" ch. "+to_string(chan)+" Pedestal Fit";
	TH1F *h = new TH1F("h", "", bins, low, high);
	h->GetXaxis()->SetTitle("signal amplitude [ADC]");
	TString ytitle = "counts";
	h->GetYaxis()->SetTitle(ytitle);
	h->GetYaxis()->SetTitleOffset(1.4);
	h->SetLineWidth(3);

	//loop over events in tree
	for (int ientry=0; ientry < nentries; ientry++ )
	{
		tree->GetEntry(ientry);
		if ( mac5==mac)  h->Fill(chg[chan]);
	}

	int max_bin = h->GetMaximumBin();
	int max_val = h->GetMaximum();
	int max_adc = h->GetBinLowEdge(max_bin);

	h->GetYaxis()->SetRangeUser(0,(max_val+h->GetBinError(max_bin))*1.05);
	h->SetTitle(htitle);

   	//Define fitting function:
	TF1 *gausfit = new TF1("gausfit","[0]*exp(-0.5*((x-[1])/[2])^2)",max_adc-12,max_adc+12);

	//Set parameter names:
	gausfit->SetParName(0,"Constant");
	gausfit->SetParName(1,"Peak value");
	gausfit->SetParName(2,"sigma");

	//Initial guesses for the parameters:
	gausfit->SetParameter(0,max_val);
	gausfit->SetParLimits(0,0.5*max_val,1000*max_val);
	gausfit->SetParameter(1,max_adc);
	gausfit->SetParLimits(1,max_adc-20,max_adc+20);
	gausfit->SetParameter(2,h->GetStdDev());
	gausfit->SetParLimits(2,1,50);
	
	TCanvas *c = new TCanvas();
	c->SetGrid();

	gStyle->Reset("Modern");
	if(h->GetMean()>(high-low)/2)
		gStyle->SetStatX(0.4);
	else
		gStyle->SetStatX(0.9);
	
	gStyle->SetStatY(0.89);
	gStyle->SetStatH(0.45);
	gStyle->SetStatW(0.28);
	h->Draw();
	h->Draw("E0same");//hist");
	h->Fit(gausfit,"RMQ");

	Float_t csq = gausfit->GetChisquare(); //chi-squared
	Float_t ndf = gausfit->GetNDF(); //NDF
	Float_t rcsq = csq/ndf; //reduced chi-squared

	if(rcsq>10.0)
	{
		//cout << "X^2/NDF > 10...recursive refit..." << endl;
		gausfit->SetRange(gausfit->GetParameter(1)-10,gausfit->GetParameter(1)+10);
		gausfit->SetParameter(0,gausfit->GetParameter(0));
		gausfit->SetParLimits(0,gausfit->GetParameter(0)-50,gausfit->GetParameter(0)+50);
		gausfit->SetParameter(1,gausfit->GetParameter(1));
		gausfit->SetParLimits(1,gausfit->GetParameter(1)-10,gausfit->GetParameter(1)+10);
		gausfit->SetParameter(2,gausfit->GetParameter(2));
		gausfit->SetParLimits(2,gausfit->GetParameter(2)-10,gausfit->GetParameter(2)+10);
		h->Fit(gausfit,"RMQ");
	}

	statarr[0]=gausfit->GetParameter(0); //const
	statarr[1]=gausfit->GetParError(0);  //const err
	statarr[2]=gausfit->GetParameter(1); //mean
	statarr[3]=gausfit->GetParError(1);  //mean err
	statarr[4]=gausfit->GetParameter(2); //sigma
	statarr[5]=gausfit->GetParError(2);  //sigma err

	//Adding chi^2/ndf and fit parameters to stat box:
	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(111);
	c->Update();

	csq = gausfit->GetChisquare();
	ndf = gausfit->GetNDF();
	rcsq = csq/ndf;
	if (rcsq>200.0) cout << "warning: possibly bad ped fit mac5 " << mac << ", ch. " 
			<< chan << " X^2/NDF=" << rcsq << " ADC" << endl;

	//save plot to file if desired
 	if ( save_opt == 1 )
	{
		int ctr = 1;
		TString suff=".png";
		TString fname = calPlotDir+"/";		 
		fname+=FormatMacString(mac)+ "_ch";
		if(chan<10) fname+="0";
		fname+=to_string(chan)+"_pedfit_"+to_string(ctr)+suff;

		while (!gSystem->AccessPathName(fname))
		{
			fname.Remove(fname.Length()-suff.Length(),fname.Length());
			fname.Remove(fname.Length()-1,fname.Length());
			ctr++;
			fname+=to_string(ctr)+suff;
		}

		//write image to file
		TImage *img = TImage::Create();
		img->FromPad(c);
		img->WriteImage(fname);

		delete img;
		delete[] peds;
		delete gausfit;
		delete h;
		delete c;
	}//end save opt

	return statarr;
}//end def

Double_t* mppc::InitGain(Int_t mac)
{
	if(!IsMacConfigured(mac)) {
		cout << "ERROR in InitGain: unconfigured mac addressed passed" << endl;
		return {};
	}

	Double_t *arr = new Double_t[32];
	Double_t val;
	Int_t nlines = 0;

	ifstream fgain;
	TString path = CAL+"conversion_init_"+FormatMacString(mac)+".txt";
	fgain.open(path);

	fgain >> val;

	if(fgain.good()) while (fgain.good())
	{
		if (nlines < 32 ) arr[nlines] = val;
		else cout << "Too many lines in gain cal file! Ignoring extra lines..." << endl;  
		nlines++;
		fgain >> val;
	}
	else {
		cout << "ERROR in InitGain: could not open file '" << path << "'" << endl;
		return {};
	}

	fgain.close();

	return arr;
}

Double_t* mppc::InitPed(Int_t mac)
{
	if(!IsMacConfigured(mac)) {
		cout << "ERROR in InitPed: unconfigured mac addressed passed" << endl;
		return {};
	}

	Double_t *arr = new Double_t[32];
	Double_t val;
	Int_t nlines = 0;

	TString fname = "pedestal_init_" + FormatMacString(mac) + ".txt";
	TString path = CAL+fname;
	ifstream fped;
	fped.open(path);

	fped >> val;

	if (fped.good()) {
		cout << "InitPed: reading from " << fname << endl;
		while (fped.good())
		{
			if (nlines < 32 ) arr[nlines] = val;
			else cout << "Too many lines in ped init file! Ignoring extra lines..." << endl;  
			nlines++;
			fped >> val;
		}
	}
	else
		cout << "ERROR in InitPed: could not open file " << path << endl;

	fped.close();

	return arr;
}

Double_t* mppc::GetGain(Int_t mac)
{
	Double_t *arr = new Double_t[32];
	Double_t val=0.0, chi=0.0, tmp=0.0;
	Int_t nlines = 0, ndf=0;

	TString fname = calDir+"/gains"+FormatMacString(mac)+".txt";
	ifstream fgain;

	if (!gSystem->AccessPathName(fname))
	{
		fgain.open(fname);

		fgain >> val; if(val!=0) {fgain >> tmp; fgain >> tmp; fgain >> tmp; fgain >> chi; fgain >> ndf;}

		while (fgain.good())
		{
			if (nlines < 32 )
			{
				arr[nlines] = val;
				//cout << "nline = " << nlines << " , val = " << val << endl;}
			}
			else cout << "Too many lines in gain cal file! Ignoring extra lines..." << endl;  
			nlines++;
		
			fgain >> val; if(val!=0) {fgain >> tmp; fgain >> tmp; fgain >> tmp; fgain >> chi; fgain >> ndf;}
		}

		fgain.close();
	}
	else {
		cout << "gain file - " << fname << " - not found!" << endl;
		return nullptr;
	}

	return arr;
}

Double_t* mppc::GetPed(Int_t mac)
{
	if(!IsMacConfigured(mac)) {
		cout << "ERROR in GetPed: unconfigured mac addressed passed" << endl;
		return nullptr;
	}

	Double_t *arr = new Double_t[32];
	Double_t val=0.0, tmp=0.0;
	Int_t nlines = 0;
	TString fname = calDir+"/peds"+FormatMacString(mac)+".txt";

	ifstream fped; 

	if (!gSystem->AccessPathName(fname))
	{
		fped.open(fname);
		fped >> val; if(val!=0) {fped >> tmp; fped >> tmp; fped>>tmp; fped>>tmp; fped>>tmp;}

		while (fped.good())
		{	
			if (nlines < 32 ) arr[nlines] = val;
			else cout << "Too many lines in ped cal file! Ignoring extra lines..." << endl;  
			nlines++;
			fped >> val; if(val!=0) {fped >> tmp; fped >> tmp; fped>>tmp; fped>>tmp; fped>>tmp;}
		}

		fped.close();
	}

	else {
		cout << "ped file - " << fname << " - not found!" << endl;
		return nullptr;
	}

	return arr;
}

Double_t* mppc::GetPedWidth(Int_t mac)
{
    if(!IsMacConfigured(mac)) {
        cout << "ERROR in GetPedWidth: unconfigured mac addressed passed" << endl;
        return {};
    }

    Double_t *arr = new Double_t[32];
    Double_t val=0.0, tmp=0.0;
    Int_t nlines = 0;
    TString fname = calDir+"/peds"+FormatMacString(mac)+".txt";

    ifstream fped;

    if (!gSystem->AccessPathName(fname))
    {
        fped.open(fname);
        fped >> tmp; if(tmp!=0) {fped >> tmp; fped >> val; fped>>tmp; fped>>tmp; fped>>tmp;}

        while (fped.good())
        {
            if (nlines < 32 ) arr[nlines] = val;
            else cout << "Too many lines in ped cal file! Ignoring extra lines..." << endl;
            nlines++;
            fped >> tmp; if(tmp!=0) {fped >> tmp; fped >> val; fped>>tmp; fped>>tmp; fped>>tmp;}
        }

        fped.close();
    }

    else cout << "file - " << fname << " - not found!" << endl;
    return arr;
}

//**************************************************************************************
//**************************************************************************************
//******		  END CALIBRATION FUNCTION DEFINITIONS		    ******
//**************************************************************************************
//**************************************************************************************

void mppc::FillLYArr(Int_t mac, Int_t (&n)[32], Double_t (&peak)[32], Double_t (&fwhm)[32], Double_t(&w)[32],
			Double_t (&werr)[32], Double_t (&mpv)[32], Double_t (&mpverr)[32], Double_t (&norm)[32],
			Double_t (&normerr)[32], Double_t (&gsig)[32], Double_t (&gsigerr)[32] )
{
	if(!IsMacConfigured(mac)) {
		cout << "ERROR in FillLYArr: unconfigured mac addressed passed" << endl;
		return;
	}

	Int_t nlines = 0;
	TString fname = anaDir+"/data_lyfit_thr2.5_mac"+FormatMacString(mac)+"_allchan.txt";

	if (!gSystem->AccessPathName(fname))
	{
		ifstream fdat;
		fdat.open(fname);

        while (fdat.good())
        {
            if(nlines==33){
                cout << "Too many lines in gain cal file! Ignoring extra lines..." << endl;
                break;
            }

            string tmp;
            fdat >> tmp;
            if(tmp!="0"&&tmp!=""&&nlines<32) {
                n[nlines]=(Int_t)atoi(tmp.c_str());
                fdat >> peak[nlines];
                fdat >> fwhm[nlines];
                fdat >> w[nlines];
                fdat >> werr[nlines];
				fdat >> mpv[nlines];
				fdat >> mpverr[nlines];
				fdat >> norm[nlines];
				fdat >> normerr[nlines];
				fdat >> gsig[nlines];
				fdat >> gsigerr[nlines];
            }
            if(tmp=="0"&&nlines<32){
                n[nlines]=0.0;
                peak[nlines]=0.0;
                fwhm[nlines]=0.0;
                w[nlines]=0.0;
                werr[nlines]=0.0;
                mpv[nlines]=0.0;
				mpverr[nlines]=0.0;
				norm[nlines]=0.0;
				normerr[nlines]=0.0;
				gsig[nlines]=0.0;
				gsigerr[nlines]=0.0;
            }
            if(tmp!="")
                nlines++;
        }

        fdat.close();
    }

    else cout << "ERROR: file '" << fname << "' not found! Arrays not filled!" << endl;
	return;
}

void mppc::FillGainArr( Int_t mac, Double_t (&g)[32], Double_t (&gerr)[32], Double_t (&gped)[32], 
					Double_t (&gpederr)[32], Double_t (&chi)[32], Double_t (&ndf)[32]) 
{
    if(!IsMacConfigured(mac)) {
        cout << "ERROR in FillGainArr: unconfigured mac addressed passed" << endl;
        return;
    }

    Int_t nlines = 0;
    TString fname = calDir+"/gains"+FormatMacString(mac)+".txt";

    if (!gSystem->AccessPathName(fname))
    {
        ifstream fgain;
        fgain.open(fname);

        while (fgain.good())
        {
            if(nlines==33){
                cout << "Too many lines in gain cal file! Ignoring extra lines..." << endl;
                break;
            }

            string tmp;
            fgain >> tmp;
            if(tmp!="0"&&tmp!=""&&nlines<32) {
                g[nlines]=(Double_t)atof(tmp.c_str());
                fgain >> gerr[nlines];
                fgain >> gped[nlines];
                fgain >> gpederr[nlines];
                fgain >> chi[nlines];
                fgain >> ndf[nlines];
            }
            if(tmp=="0"&&nlines<32){
                g[nlines]=0.0;
                gerr[nlines]=0.0;
                gped[nlines]=0.0;
                gpederr[nlines]=0.0;
                chi[nlines]=0.0;
                ndf[nlines]=0.0;
            }
            if(tmp!="")
                nlines++;
        }

        fgain.close();
    }

    else cout << "ERROR: file '" << fname << "' not found! Arrays not filled!" << endl;
    return;
}

void mppc::FillPedArr(Int_t mac, Double_t (&mean)[32], Double_t (&meanerr)[32],
					Double_t (&sig)[32], Double_t (&sigerr)[32], Double_t (&cons)[32], Double_t (&conserr)[32])
{
    if(!IsMacConfigured(mac)) {
        cout << "ERROR in FillPedArr: unconfigured mac addressed passed" << endl;
        return;
    }

    Int_t nlines = 0;
    TString fname = calDir+"/peds"+FormatMacString(mac)+".txt";

    if (!gSystem->AccessPathName(fname))
    {
		ifstream fped;
		fped.open(fname);

        while (fped.good())
        {
			if(nlines==33){
				cout << "Too many lines in ped cal file! Ignoring extra lines..." << endl;
				break;
			}

			string tmp;
			fped >> tmp;
			if(tmp!="0"&&tmp!=""&&nlines<32) {
				mean[nlines]=(Double_t)atof(tmp.c_str());
				fped >> meanerr[nlines];
				fped >> sig[nlines];
				fped >> sigerr[nlines];
				fped >> cons[nlines];
				fped >> conserr[nlines]; 
			}
			if(tmp=="0"&&nlines<32){
				mean[nlines]=0.0;
				meanerr[nlines]=0.0;
				sig[nlines]=0.0;
				sigerr[nlines]=0.0;
				cons[nlines]=0.0;
				conserr[nlines]=0.0;
			}
			if(tmp!="")
				nlines++;
		}

        fped.close();
    }

    else cout << "ERROR: file '" << fname << "' not found! Arrays not filled!" << endl;
	return;
}

Int_t mppc::GetNentries(char runtype='h', Int_t mac=85)
{
    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
	cout << "error in GetNentries: unknown runtype provided - " << runtype << endl;
	return -1; 
    } 

	Int_t nentries_tot = tree->GetEntriesFast();
	Int_t nentries = 0;

	for (Int_t i=0; i<nentries_tot; i++)
	{
		tree->GetEntry(i);
		if(mac5==mac) nentries++;
	}

	return nentries;
}

vector<Int_t>* mppc::GetRawArr(char runtype='h', Int_t mac=85)
{
    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
	cout << "error in GetRawArr: unknown runtype provided - " << runtype << endl;
	return {}; 
    } 

	Int_t nentries_tot = tree->GetEntriesFast();
	Int_t nentries = GetNentries(mac);
	vector<Int_t> *arr = new vector<Int_t>[20];		

	for (Int_t i=0; i<nentries_tot; i++)
	{
		tree->GetEntry(i);
		if(mac5==mac)
		{
			for(Int_t j=1; j<21; j++) arr[j-1].push_back(chg[j]);
		}
	}

	return arr;
}

vector<Double_t>* mppc::GetPEArr(char runtype = 'h', Int_t mac=85)
{
    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
	cout << "error in GetPEArr: unknown runtype provided - " << runtype << endl;
	return {}; 
    } 

	Int_t nentries_tot = tree->GetEntriesFast();
	Int_t nentries = GetNentries(mac);
	vector<Double_t> *arr = new vector<Double_t>[20];		

	Double_t *conv = new Double_t[32];
	Double_t *ped = new Double_t[32];

	conv=GetGain(mac); ped=GetPed(mac);

	for (Int_t i=0; i<nentries_tot; i++)
	{
		tree->GetEntry(i);
		if(mac5==mac)
		{
			for(Int_t j=1; j<21; j++) arr[j-1].push_back((chg[j]-ped[j])/conv[j]);
		}
	}

	return arr;
}

Float_t mppc::GetMaxInRange(TH1F *htemp, Float_t low, Float_t high)
{
	Int_t nbins = htemp->GetNbinsX(), val = 0;
	Float_t max = 0.0, edge = 0.0;

	for ( Int_t i = 1; i < nbins; i++ )
	{
		edge = htemp->GetBinLowEdge(i);
		val = htemp->GetBinContent(i);

		if ( edge > low && edge < high && val > max)
		{
			max = val;
		}
		else if (edge > high ) break;
	}

	return max;	
}

Int_t mppc::GetMaxBinInRange(TH1F *htemp, Float_t low, Float_t high)
{
	Int_t nbins = htemp->GetNbinsX(), val = 0, maxbin = 0;
	Float_t max = 0.0, edge = 0.0;

	for ( Int_t i = 1; i < nbins; i++ )
	{
		edge = htemp->GetBinLowEdge(i);
		val = htemp->GetBinContent(i);

		if ( edge > low && edge < high && val > max)
		{
			max = val;
			maxbin = i;
		}
		else if (edge > high ) break;
	}

	return maxbin;	
}


///////////////////////////////////////////////////////////////////////////////////
//	methods for dealing with time stamps
//////////////////////////////////////////////////////////////////////////////////
void mppc::CountEmptyEvents(char runtype='h', int mac=-1) 
{
/*    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
		cout << "error in ThresholdStats: unknown runtype provided - " << runtype << endl;
		return; 
    } 

	Int_t nentries = tree->GetEntriesFast();

	int nEmpty = 0;
	for (Int_t i=0; i<nentries; i++)
	{
		tree->GetEntry(i);
		if(mac>-1&&mac5!=mac) continue;

		bool isEmpty = true;
		for (Int_t chan=0; chan<32; chan++) {
		    if (chg[chan]!=0) {
			isEmpty = false;
			break;
		    }
		}

		if (isEmpty) nEmpty++;
	}

	cout << " found " << nEmpty << " empty events" << endl;
}

void mppc::DumpEvents(char runtype='h', Int_t mac=85) 
{
    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
        cout << "error in DumpEvents: unknown runtype provided - " << runtype << endl;
        return; 
    } 

	Int_t nentries_tot = tree->GetEntriesFast();

	cout << "|  entry   |  mac5  |  chg[1]   |  ts0   |   ts1  |  ts0_ref  |  ts1_ref  |" << endl;

	for (Int_t i=0; i<4000; i++)
	{
		tree->GetEntry(i);
		if (mac5!=mac&&mac!=-1) continue;

		cout << "| "  << i << " | " << (int)mac5 << " | " << chg[1]
		     << " | " << ts0 << " | " << ts1 << " | " << ts0_ref << " | " << ts1_ref << " |" << endl;

	}*/
}

void mppc::CheckEventPacket(char runtype='h') 
{/*
    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
        cout << "error in ThresholdStats: unknown runtype provided - " << runtype << endl;
        return; 
    } 

	Int_t nentries_tot = tree->GetEntriesFast();

	for (Int_t ientry=0; ientry<nentries_tot; ientry++) {
	    tree->GetEntry(ientry);
	    if(ts0==0) {ientry++; tree->GetEntry(ientry);}
	    else {
		cout << "broken event packet! entry " << ientry << endl;
		break;
	    }
	    if(ts0!=0&&ts1!=0) {ientry++; tree->GetEntry(ientry);}
	    else {
		cout << "broken event packet! entry " << ientry << endl;
		break;
	    }
	    if(ts1!=0) {
		cout << "broken event packet! entry " << ientry << endl;
		break;
	    }
	}*/
}

void mppc::PlotADCTimeSlice(char runtype='h', Int_t mac=213, Int_t adclow=0, Int_t adchigh=4080) 
{
/*
    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
        cout << "error in PlotADCTimeSlice: unknown runtype provided - " << runtype << endl;
        return; 
    } 

    Int_t nentries=tree->GetEntriesFast();

    UInt_t xlow=30120, xhigh=30180;
    UInt_t bins = xhigh-xlow;
    TH1F* h = new TH1F("h","ts1_ref-ts1",bins,xlow,xhigh);

    Bool_t foundT0 = kFALSE, foundT1 = kFALSE;
    UInt_t signalTime, dt;
    UShort_t sigmax;

    for(Int_t ientry=0; ientry<nentries; ientry++) {

	tree->GetEntry(ientry);
	if (mac5!=mac) continue;

	if(!foundT0&&ts0==0) foundT0 = kTRUE;
	if(!foundT1&&ts1==0) foundT1 = kTRUE;

	//not a time stamp event
	if (ts0!=0&&ts1!=0) {

	    if (!foundT0||!foundT1) continue;
	    signalTime = ts1;//*1e9/ts0_ref;
	    sigmax = chg[0];
	    for (Int_t i=0; i<10; i++) if(sigmax<chg[i]) sigmax=chg[i];
	    if(sigmax<adclow||sigmax>adchigh) continue;

	    for (Int_t jentry=ientry+1; jentry<ientry+5000; jentry++) {
		if (jentry>nentries-1) break;
		//if (jentry>entryHigh-1) break;
		tree->GetEntry(jentry);
		if (mac5!=mac) continue;
		//if time stamp event was missed
		if (ts1!=0&&ts1<signalTime) break;
		if (ts1==0) {
		    dt=ts1_ref-signalTime;
		    if(dt>xlow&&dt<xhigh) {
			h->Fill(dt);
		    }
		    break;
		}
		if(jentry==ientry+4999) cout << "LED time stamp not found for ientry " << ientry << "!" << endl;
	    }
	}

    }

    h->Draw();
    cout << "disribution RMS[ns]: " << h->GetRMS() << endl;
*/
}

void mppc::TimeResolution(char runtype='h', Int_t mac=85) 
{
/*    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
        cout << "error in TimeResolution: unknown runtype provided - " << runtype << endl;
        return; 
    } 

    Int_t nentries=tree->GetEntriesFast();

    UInt_t xlow=30120, xhigh=30180;
    UInt_t bins = xhigh-xlow;
    TString htitle = "signal delay: mac"; htitle+=mac;
    TH1F* h = new TH1F("h",htitle,bins,xlow,xhigh);
    TH2F* h2 = new TH2F("h2","",bins,xlow,xhigh,(4100-1000)/100,1000,4100);

    TH1F* hmaxsig_intime = new TH1F("hmaxsig_intime","maximum signal in event",4000/5,80,4080);
    TH1F* hmaxsig_outtime = new TH1F("hmaxsig_outtime","maximum signal in event",4000/5,80,4080);

    Bool_t foundT0 = kFALSE, foundT1 = kFALSE;
    UInt_t signalTime, dt;
    UShort_t sigmax;
    Int_t entryLow = 20e3, entryHigh=190e3; 
    //for(Int_t ientry=entryLow; ientry<entryHigh; ientry++) {
    for(Int_t ientry=0; ientry<nentries; ientry++) {

	tree->GetEntry(ientry);
	if (mac5!=mac) continue;

	if(!foundT0&&ts0==0) foundT0 = kTRUE;
	if(!foundT1&&ts1==0) foundT1 = kTRUE;

	//not a time stamp event
	if (ts0!=0&&ts1!=0) {
 
	    if (!foundT0||!foundT1) continue;
	    signalTime = ts1;//*1e9/ts0_ref;
	    sigmax = chg[0];
	    for (Int_t i=0; i<10; i++) if(sigmax<chg[i]) sigmax=chg[i]; 

	    //loop to find corresponding T1
	    for (Int_t jentry=ientry+1; jentry<ientry+5000; jentry++) {
		if (jentry>nentries-1) break;
		//if (jentry>entryHigh-1) break;
		tree->GetEntry(jentry);
		if (mac5!=mac) continue;
		//if time stamp event was missed
		if (ts1!=0&&ts1<signalTime) break;
		if (ts1==0) {
		    dt=ts1_ref-signalTime;
		    if(dt<xlow||dt>xhigh) {
			//cout << "dt out of bounds! dt=" << dt <<endl;
			hmaxsig_outtime->Fill(sigmax);   
		    }
		    else {
			h->Fill(dt);
			hmaxsig_intime->Fill(sigmax);
			h2->Fill(dt,sigmax);   
		    }
		    break;
		}
		if(jentry==ientry+4999) cout << "LED time stamp not found for ientry " << ientry << "!" << endl;
	    }
	}

    }

    const Int_t size = h2->GetNbinsY();
    Float_t meandelays[size], meanadc[size], errdelays[size], erradc[size];
    for (Int_t biny=1; biny<h2->GetNbinsY()+1; biny++) {
	Float_t mean=0.0, sig=0.0;
	Int_t counts=0;
	for (Int_t binx=1; binx<h2->GetNbinsX()+1; binx++) {
	    counts+=h2->GetBinContent(binx,biny);
	    mean = mean + h2->GetBinContent(binx,biny) * ((TAxis*)h2->GetXaxis())->GetBinCenter(binx);
	}
	if(counts>1) for (Int_t binx=1; binx<h2->GetNbinsX()+1; binx++) {
	    sig+=pow(mean/counts-((TAxis*)h2->GetXaxis())->GetBinCenter(binx),2);
	}
	meanadc[biny-1] = ((TAxis*)h2->GetYaxis())->GetBinCenter(biny);
	erradc[biny-1] = ((TAxis*)h2->GetYaxis())->GetBinWidth(biny)/sqrt(12);
	if(counts>1) {meandelays[biny-1] = mean/counts; errdelays[biny-1] = sqrt(sig/(counts-1));}
	else {meandelays[biny-1] = 0; errdelays[biny-1]=0;}
    }

    TF1 *fit = new TF1("f","[0]*log(x)+[1]",500,4050);
    fit->SetParameter(0,1);
    fit->SetParameter(1,30000);

    TF1 *fit2 = new TF1("f2","[0]+x^[1]",500,4050);
    fit2->SetParameter(0,31000);
    fit2->SetParameter(1,0.5);

    //TGraphErrors *g = new TGraphErrors(size,meanadc,meandelays,erradc,errdelays);
    TGraph *g = new TGraph(size,meanadc,meandelays);
    g->SetMarkerStyle(8);
    g->GetXaxis()->SetTitle("signal amplitude [ADC]");
    g->GetYaxis()->SetTitle("mean delay [ns]");
    TString gtitle = "time walk: mac"; gtitle+=mac;
    g->SetTitle(gtitle);

    h->GetXaxis()->SetTitle("delay [ns]");
    h->Draw();

    hmaxsig_intime->SetLineColor(kBlue);
    hmaxsig_outtime->SetLineColor(kRed);
    hmaxsig_intime->GetXaxis()->SetTitle("signal amplitude [ADC]");

    new TCanvas();
    hmaxsig_intime->Draw();
    hmaxsig_outtime->Draw("same");

    h2->GetXaxis()->SetTitle("signal delay [ns]");
    h2->GetYaxis()->SetTitle("signal amplitude [ADC]");
    TString h2title = "time walk with amplitude: mac"; h2title+=mac;
    h2->SetTitle(h2title);

    new TCanvas();
    h2->Draw("COLZ");

    new TCanvas();
    g->Draw("AP");
    //fit->Draw("same");
*/
}

void mppc::MissedTimeStamps(char runtype='h', Int_t mac=85) 
{
/*    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
        cout << "error in MissedTimeStamps: unknown runtype provided - " << runtype << endl;
        return;
    }

    Int_t nentries = tree->GetEntriesFast();
    cout << "reading " << nentries << " entries..." << endl;

    Int_t ientry=0;
    UInt_t ts0_current=0, nmiss0=0, nppssynch=0;
    UInt_t ts1_current=0, nmiss1=0;

    while(ientry<nentries) {

	//if (ientry>nentries-10) cout << ientry << endl;

	tree->GetEntry(ientry);
	if(ientry%100==0) cout << "ientry: " << ientry << endl;	

	if(mac5!=mac) {
	    ientry++;
	    continue;
	}

	if (ts0==0) nppssynch++;
	if (ts0==0) for(Int_t jentry = ientry+1; jentry<ientry+10000; jentry++) {

	     if (jentry==nentries-1) {
		ientry=jentry;
		break;
	    }

	    tree->GetEntry(jentry);
	    if(mac5!=mac) continue;

	    if (ts0==0) {
		ientry = jentry;
		break;
	    }

	    ts0_current = ts0;
	    tree->GetEntry(jentry+1);
	    if(ts0!=0&&ts0<ts0_current) {
		//cout << "missed time stamp during entry" << jentry << endl;
		nmiss0++;
	    }
	    if(jentry==ientry+9999) cout << "at end of jentry loop" << endl;
	}

	else if(ts1==0) for(Int_t jentry = ientry+1; jentry<ientry+10000; jentry++) {

	     if (jentry==nentries-1) {
		ientry=jentry;
		break;
	    }

	    tree->GetEntry(jentry);
	    if (mac5!=mac) continue;

	    if (ts1==0) {
		ientry = jentry;
		break;
	    }

	    ts1_current = ts1;
	    tree->GetEntry(jentry+1);
	    if(ts1!=0&&ts1<ts1_current) {
		cout << "missed time stamp during entry " << jentry << endl;
		nmiss1++;
	    }
	    if(jentry==ientry+9999) cout << "at end of jentry loop" << endl;
	}

	else ientry++;




    }

    cout << nppssynch << " (" << nppssynch/60.0 << ") sec(min) in synchronization period" << endl;
    cout << "PPS: missed " << nmiss0 << " time stamp events" << endl;
    cout << "T1 : missed " << nmiss1 << " time stamp events" << endl;*/
}

void mppc::RunTime(char runtype='h', int input=0) 
{
/*    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
        cout << "error in RunTime: unknown runtype provided - " << runtype << endl;
        return; 
    } 

    Int_t nentries = tree->GetEntriesFast();
    bool hasStart85 = false, hasStart213=false;
    UInt_t ts = 0, ts_ref = 0;
    uint64_t time85 = 0, time213 = 0; 
    size_t nstamp85=0, nstamp213=0;

    TH1F *h85 = new TH1F("h85","mac85",1000,67.054e6,67.072e6);

    for (Int_t i=0; i<nentries; i++)
    {
	tree->GetEntry(i);

	if (input==0) { ts = ts0; ts_ref = ts0_ref; }
	if (input==1) { ts = ts1; ts_ref = ts1_ref; }

	if (ts==0) {
	    if (mac5==85) {
		if (!hasStart85) hasStart85 = true;
		else {
		    time85+=ts_ref;
		    h85->Fill(ts_ref);
		}
		nstamp85++;
	    }

	    if (mac5==213) {
		if (!hasStart213) hasStart213 = true;
		else time213+=ts_ref;
		nstamp213++;
	    }
	}
    }

    h85->Draw();

    cout << "runtime mac85: "  << time85*1.0e-9  << " [s] , " << time85*1.0e-9/60.0  << " [min] with " << nstamp85  << " time stamps" << endl;
    cout << "runtime mac213: " << time213*1.0e-9 << " [s] , " << time213*1.0e-9/60.0 << " [min] with " << nstamp213 << " time stamps" << endl;
*/
}

void mppc::CorrelatedChannels(char runtype='h') 
{
    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
        cout << "error in CorrelatedChannels: unknown runtype provided - " << runtype << endl;
        return; 
    } 

    Int_t nentries = tree->GetEntriesFast();

    Double_t * conv213 = GetGain(213);
    Double_t * ped213 = GetPed(213);

    Float_t thresh[5] = {5.5,6.5,7.5,8.5,9.5};
    Float_t mult[5];
    Float_t adcall[10];
    Double_t coinch2trig[10];
    Double_t chans[10]; for (int i=0; i<10; i++) {chans[i]=(Double_t)i; coinch2trig[i]=0.0;}

    TH1F *hall = new TH1F("hall","signal with all channels > 9.5PE",4000/500,100,4100);
    hall->GetXaxis()->SetTitle("signal amplitude [ADC]");

    TH1F *histos[5];
    for (Int_t h=0; h<5; h++) {
	TString name = "h"; name+=h+5; name+="_5";
	histos[h] = new TH1F(name,"",11,0,11);
	if(h<4) histos[h]->SetLineColor(h+1);
	else histos[h]->SetLineColor(h+2);
    }

    for (Int_t ientry=0; ientry<nentries; ientry++) {
 
       tree->GetEntry(ientry);

       Int_t nabove = 0;
       Float_t ly = 0.0;
       Float_t lych2 = (chg[7]-ped213[7])/conv213[7];
       if (lych2>10.5) for (Int_t ch=0; ch<10; ch++) {
	   ly = (chg[ch]-ped213[ch])/conv213[ch];
	   if (ly>10.5) coinch2trig[ch]=coinch2trig[ch]+1.0;
	   
       }

       for (Int_t thr=0; thr<5; thr++) {
	   nabove=0;
	   for (Int_t ch=0; ch<10; ch++) {
	       ly = (chg[ch]-ped213[ch])/conv213[ch];
	       if (ly>thresh[thr]) {
		   nabove++;
		   if (thr==4) adcall[ch] = chg[ch];
	       }
	   } // loop over channels
	   mult[thr] = nabove;
       } // loop over thresholds

       if(mult[4]==10) for(Int_t i=0; i<10; i++) hall->Fill(adcall[i]);

       for (Int_t h=0; h<5; h++) {
	  histos[h]->Fill(mult[h]);
       }
    } //loop over entries

    histos[0]->GetXaxis()->SetTitle("no. channels above threshold");

    TLegend* leg = new TLegend(0.75,0.6,0.9,0.9);
    for (Int_t h=0; h<5; h++) {
       histos[h]->Sumw2();
       histos[h]->Scale(1.0/histos[h]->Integral());
       leg->AddEntry(histos[h],histos[h]->GetName(),"l");
    }

    Double_t nch2 = coinch2trig[7];
    for (int i=0; i<10; i++) {
	coinch2trig[i] = coinch2trig[i]/nch2;
    }

    TGraph* g = new TGraph(10,chans,coinch2trig);
    g->SetTitle("events with ch. 7 > 10.5PE");
    g->GetXaxis()->SetTitle("channel");
    g->GetYaxis()->SetTitle("fraction of events in coincidence with ch. 7");
    g->SetMarkerStyle(8);
    g->SetMarkerColor(kBlue);

    TCanvas *c = new TCanvas();
    c->SetLogy();
    histos[0]->Draw("e0hist");
    for (Int_t h=1; h<5; h++) histos[h]->Draw("e0histsame");
    leg->Draw();

    new TCanvas();
    hall->Draw("e0hist");

    new TCanvas();
    g->Draw("ap");
}

void mppc::CountClockResets(char runtype='h', Int_t mac=0) 
{
/*    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
        cout << "error in CountClockResets: unknown runtype provided - " << runtype << endl;
        return; 
    } 

    Int_t nentries = tree->GetEntriesFast();
    Int_t nreset=0,nreset_nonempty=0;
    Int_t nreset_hod=0,nmiss_hod=0;

    for(Int_t ientry=0; ientry<nentries-1; ientry++) {

	tree->GetEntry(ientry);

	if(mac5!=mac) continue;

	UInt_t ts0_current = ts0;
	UInt_t ts1_current = ts1;
	UInt_t ts0_ref_current = ts0_ref;
	Int_t jentry=ientry+1;
	//loop to find next entry with same mac5
	while(jentry<nentries) {
	    tree->GetEntry(jentry);
	    if (mac5!=mac) jentry++;
	    else {
		if(ts0<ts0_current) {
		    nreset++;
		    if (ts0!=0) {
			 nreset_nonempty++;
			 //cout << " missed PPS event!" << '\n' 
			 //     << "  ts0_curr-> " << ts0_current << " , ts0_next-> " << ts0 << '\n'
			 //     << "  ts0_ref_cur-> " << ts0_ref_current << " , ts0_ref_next-> " << ts0_ref << endl;
		    }
		    //cout << " clock reset at entry " << jentry <<": "<< '\n'
		    //	  << "  ts0_curr-> " << ts0_current << " , ts0_next-> " << ts0 << '\n'
		    //	  << "  ts0_ref_cur-> " << ts0_ref_current << " , ts0_ref_next-> " << ts0_ref << '\n' << endl;
		}
		if(ts1<ts1_current) nreset_hod++;
		if (ts1!=0&&ts1<ts1_current) nmiss_hod++;

		ientry=jentry-1; //will set ientry to jentry after incrementing at end of for loop
		break;
	    }//else
	}//while

    }//loop over entries

    cout << "found " << nreset << " clock resets amounting to " << nreset/60.0 << " (" << nreset/3600.0 << ") min(hr)" << endl;
    cout <<" found " << nreset_nonempty << " (" << 100.0*nreset_nonempty/nreset << "%) missing PPS events" << endl;
    cout <<" found " << nreset_hod << " T1 resets and " << nmiss_hod << " (" << 100.0*nmiss_hod/nreset_hod << "%) missing T1 events" << endl;
*/
}
//end file

void mppc::CountOverflows(char runtype='h') 
{
/*    TTree* tree;
    if(runtype=='h') tree = fTreeInner;
    else if(runtype=='s') tree = fTreeOuter;
    else {
        cout << "error in CountOverflows: unknown runtype provided - " << runtype << endl;
        return; 
    } 

    Int_t nentries = tree->GetEntriesFast();
    Int_t neve85_over_0=0, neve85_over_1=0;
    Int_t neve213_over_0=0, neve213_over_1=0;

    for(Int_t ientry=0; ientry<nentries; ientry++) {

	tree->GetEntry(ientry);

	if(mac5==85) {
	    if(ts0>1074e6) neve85_over_0++;
	    if(ts1>1074e6) neve85_over_1++;
	}

	if(mac5==213) {
	    if(ts0>1074e6) neve213_over_0++;
	    if(ts1>1074e6) neve213_over_1++;
	}

    }

    cout << "events with overflow of coarse counter:" << '\n'
	 << "  mac85:  ts0-> " << neve85_over_0 << " , ts1-> " << neve85_over_1 << '\n'
	 << "  mac213: ts0-> " << neve213_over_0 << " , ts1-> " << neve213_over_1 << endl;
*/
}

void mppc::PlotTS0(Int_t mac) {
        TTree* tree=0;
        if(mactolay[mac]=='i') tree = fTreeInner;
        else if(mactolay[mac]=='o') tree=fTreeOuter;
        const int nentries = tree->GetEntriesFast();
	int n=0, nflag=0, nmiss=0;

	for(int ientry=0; ientry<nentries; ientry++) {
		tree->GetEntry(ientry);
		if(mac5==mac) {
			n++; 
			if(flags==7) nflag++;
		}
	}

	const int size=n;
	const int sizeflag=nflag; 
	int* t = new int[size];
	int* ent = new int[size]; 
	int tflag[sizeflag], entflag[sizeflag]; 
	n=0, nflag=0;

        for(int ientry=0; ientry<nentries; ientry++) {
                tree->GetEntry(ientry);
                if(mac5==mac) {
			t[n] = ts0;
			ent[n] = n;
			n++;
			if(flags==7) {
				tflag[nflag]=ts0;
				entflag[nflag]=ent[n-1];
				nflag++;
			}
		}
        }

	vector<int> vtmiss, ventmiss;
	for(int i=0; i<size-1; i++) {
		bool isflag=false;
		if(1e9-t[i]<1e7&&t[i+1]<1e7) {
			for(int j=0; j<sizeflag; j++) {
				if(entflag[j]==ent[i]) {
					isflag=true;
					break;
				}
			}
			if(!isflag){
				nmiss++;
				vtmiss.push_back(t[i]);
				ventmiss.push_back(ent[i]);
			}
			else isflag=false;
		}
	}

	const int sizemiss=nmiss;
	int tmiss[sizemiss], entmiss[sizemiss];
	for(int i=0; i<nmiss; i++){
		tmiss[i]=vtmiss[i];
		entmiss[i]=ventmiss[i];
	}

	cout << "found " << nflag << " flags and " << nmiss << " missed flags" << endl;
	cout << " from this, estimated run time = " << 1.0*(nflag+nmiss)/60 << " min" << endl;

	TGraph *g = new TGraph(size,ent,t);
	TGraph* gflag = new TGraph(sizeflag,entflag,tflag);
	//TGraph* gmiss = new TGraph(sizemiss,entmiss,tmiss);
	gflag->SetMarkerStyle(8);
	gflag->SetMarkerColor(kBlue);
	//gmiss->SetMarkerStyle(8);
	//gmiss->SetMarkerColor(kRed);
	new TCanvas();
	g->Draw("ap");
	gflag->Draw("samep");
	for(int i=0; i<nmiss; i++){
		TLine* l = new TLine(entmiss[i],0,entmiss[i],999e6);
		l->SetLineColor(kRed);
		l->SetLineWidth(3);
		l->Draw("same");
	}
	
}

void mppc::EventRate(int mac){
        TTree* tree=0;
        if(mactolay[mac]=='i') tree = fTreeInner;
        else if(mactolay[mac]=='o') tree=fTreeOuter;
        const int nentries = tree->GetEntriesFast();

	const float tpoll=300e-3;
	int neve=0;
	bool ismac=false;
	std::vector<int> n;
	
	for(int ientry=0; ientry<nentries; ientry++) {
		tree->GetEntry(ientry);

		if(mac5==mac) {
			ismac=true;
			neve++;
		}
		else {
			if(ismac){
				n.push_back(neve);
				neve=0;
			}
			ismac=false;
		}

	}

	const size_t npoll = n.size();
	float rate[npoll/4], pollno[npoll/4];
	for(int i=0; i<npoll-4; i+=4){
		pollno[i/4] = i/4;
		rate[i/4] = (n[i]+n[i+1]+n[i+2]+n[i+3])/(4*tpoll);
	}

	TGraph* g = new TGraph(npoll/4,pollno,rate);
	string title = "Trigger Rate: FEB "+to_string(mac);
	g->SetTitle(title.c_str());
	g->GetXaxis()->SetTitle("poll number");
	g->GetYaxis()->SetTitle("TriggerRate [Hz]");


	new TCanvas();
	g->Draw("APL");
}
