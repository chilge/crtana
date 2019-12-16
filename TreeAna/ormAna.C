#define ormAna_cxx
#include "ormAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include "../src/mppc.h"

//   In a ROOT session, you can do:
//      root> .L ormAna.C
//      root> ormAna a("path/to/file")
//      root> a.UserFunction1();
//      root> a.UserFunction2();  
//

using namespace std;

//we know what FEBs we have, but his could change...be aware
//there's a better way to do this, but I'm lazy and this works for now
//vector<int> febs;
ormAna a("");

template <typename T>
double RMS(const int n, const T* arr){
	double mean=0., rms=0.;
	for(int i=0; i<n; i++) mean+=arr[i];
	mean=mean/n;
	for(int i=0; i<n; i++) rms+=pow(mean-arr[i],2);
	rms=sqrt(rms/(n-1));
	return rms;
}

void PlotFebOverlaySingle(int mac=1, int rebin=0){
	if(std::find(a.fMacs.begin(),a.fMacs.end(),mac)==a.fMacs.end()){
		cout << "provided mac address not found!" << endl;
		return;
	}
	string name = "spectra";
	if(mac<10) name+="00";
	else if(mac<100) name+="0";
	name+=to_string(mac);
	TList* l = (TList*)a.fList->FindObject(name.c_str());

        Int_t colors[32];
        for(int i=0; i<32; i++) {
                if(i<23) colors[i] = 49-i;
                else colors[i] = i-22;
        }

	gStyle->SetOptStat(0);
	gStyle->SetLegendTextSize(0.03);
	TCanvas* c = new TCanvas();
	TLegend* leg = new TLegend(0.75,0.65,0.98,0.98);
	leg->SetHeader("FEB Channels","C");
	leg->SetNColumns(3);

	for(int i=0; i<l->GetSize(); i++){
		TH1F* h = (TH1F*)l->At(i);
		h->SetLineColor(colors[i]);
		h->SetLineWidth(3);
		h->GetYaxis()->SetTitle("");
		if(rebin!=0) h->Rebin(rebin);

		string label = (string)h->GetName();
		label = label.substr(label.size()-2,label.size());
		if(mac==4&&label=="25") continue;
		if(mac==5&&label=="18") continue;
		if(i==0) {
			string title="Charge Spectra: Mac5 "+to_string(mac);
			h->SetTitle(title.c_str());
			h->GetXaxis()->SetTitle("charge amplitude [PE]");
			h->Draw();
		}
		else h->Draw("same");
                leg->AddEntry(h->GetName(),label.c_str(),"lC");

	}
	leg->Draw();
}

void PlotFebOverlayAll(char opt='a', bool norm=false, int rebin=0){

        gStyle->SetOptStat(0);
        gStyle->SetLegendTextSize(0.04);
        TCanvas* c = new TCanvas();
        TLegend* leg = new TLegend(0.75,0.65,0.98,0.98);
        leg->SetHeader("FEB mac5","C");
        leg->SetNColumns(2);

	for(int m=0; m<a.fMacs.size(); m++) {

		int mac = a.fMacs[m];
		if(opt=='m'&&(mac==3||mac==5||mac==7||mac==9)) continue; //mezzanine FEBs only
		if(opt=='p'&&(mac==1||mac==4||mac==6||mac==8)) continue; //pit FEBs only
	        string name = "spectra";
	        if(mac<10) name+="00";
	        else if(mac<100) name+="0";
	        name+=to_string(mac);
	        TList* l = (TList*)a.fList->FindObject(name.c_str());
	
	        for(int i=0; i<l->GetSize(); i++){
	                TH1F* h = (TH1F*)l->At(i);
	                h->SetLineColor(9-m);
	                h->SetLineWidth(3);
	                h->GetYaxis()->SetTitle("");
			//if(m!=0) h->SetLineStyle(11-m);
	                if(rebin!=0) h->Rebin(rebin);
			if(norm) h->Scale(1.0/h->Integral());
	
	                string label = to_string(mac);
	                if(mac==4&&label=="25") continue;
	                if(mac==5&&label=="18") continue;
	                if(i==0&&m==0) {
	                        string title="Charge Spectra: all FEBs";
	                        h->SetTitle(title.c_str());
	                        h->GetXaxis()->SetTitle("charge amplitude [PE]");
	                        h->Draw("hist");
	                }
	                else h->Draw("samehist");
	                if(i==0) leg->AddEntry(h->GetName(),label.c_str(),"lC");
	
	     	}
		
	}

        leg->Draw();
}

void PlotPedsByFeb() {
	const size_t size = a.fMacs.size()*32;
	Int_t feb[size], ped[size], n=0;
	TH1F* h = new TH1F("hpedsigma","FEB Pedestal Values RMS",15,0,15);

	for(Int_t ientry=0; ientry<a.fMacs.size(); ientry++){
		a.GetEntry(ientry);
		Int_t febpeds[32];
		for(int i=0; i<32; i++){
			feb[n] = a.mac5;
			ped[n] = a.pedMean[i];
			febpeds[i]=a.pedMean[i];
			n++;
		}
		h->Fill(RMS(32,febpeds));
	}

	h->GetXaxis()->SetTitle("RMS [ADC]");
	h->SetLineWidth(3);
	new TCanvas();
	h->Draw();

	TGraph *g = new TGraph(n,ped,feb);
	g->SetTitle("Pedestal Central Value by FEB");
	g->GetXaxis()->SetTitle("pedestal [ADC]");
	g->GetYaxis()->SetTitle("FEB mac5");
	g->SetMarkerStyle(8);
	g->SetMarkerSize(1);
	g->GetYaxis()->SetNdivisions(10,1,0);
	new TCanvas();
	g->Draw("AP");

}

void PlotGainsByFeb() {
        const size_t size = a.fMacs.size()*32;
	const size_t sizeorm = a.fMacs.size()*10;
        Int_t feb[size], g[size], n=0;
	Int_t feborm[sizeorm], gorm1[size], gorm2[size], gorm3[size];
	Int_t n1=0, n2=0, n3=0;
        TH1F* h = new TH1F("hgainrms","FEB Pedestal Values RMS",15,0,15);
	TH1F* hgain = new TH1F("hgain","FEB Gains",20,40,80);
	TH2F* hgain2 = new TH2F("hgain2","FEB Gains",20,40,80,9,0.5,9.5);

        for(Int_t ientry=0; ientry<a.fMacs.size(); ientry++){
                a.GetEntry(ientry);
                Int_t febgains[32], nchan=0;
                for(int i=0; i<32; i++){
			if(!a.chConfig[i]) continue;
                        feb[n] = a.mac5;
                        g[n] = a.gain[i];
                        febgains[nchan]=a.gain[i];
                        n++;
			nchan++;
			hgain->Fill(a.gain[i]);
			hgain2->Fill(a.gain[i],a.mac5);

			if(i>=2&&i<=11){
				feborm[n1] = a.mac5;
				gorm1[n1] = a.gain[i];
				n1++;
			}
                        if(i>=12&&i<=21){
                                gorm2[n2] = a.gain[i];
                                n2++;
                        }
                        if(i>=22&&i<=31){
                                gorm3[n3] = a.gain[i];
                                n3++;
                        }
                }
                h->Fill(RMS(nchan,febgains));
        }

        h->GetXaxis()->SetTitle("RMS [ADC]");
        h->SetLineWidth(3);
        new TCanvas();
        h->Draw();

        TGraph *gr = new TGraph(n,g,feb);
        gr->SetTitle("Gain by FEB");
        gr->GetXaxis()->SetTitle("gain [ADC/PE]");
        gr->GetYaxis()->SetTitle("FEB mac5");
        gr->SetMarkerStyle(8);
        gr->SetMarkerSize(1);
        gr->GetYaxis()->SetNdivisions(10,1,0);
        new TCanvas();
        gr->Draw("AP");
	cout << "graph points: " << gr->GetN() << endl;

	hgain->SetLineWidth(3);
	hgain->GetXaxis()->SetTitle("gain [ADC/PE]");
	new TCanvas();
	hgain->Draw();

	TGraph* gr1 = new TGraph(n1,gorm1,feborm);
	TGraph* gr2 = new TGraph(n2,gorm2,feborm);
	TGraph* gr3 = new TGraph(n3,gorm3,feborm);
	gr1->SetTitle("Gain by FEB, Grouped by ORM");
	gr1->GetXaxis()->SetTitle("gain [ADC/PE]");
	gr1->GetYaxis()->SetTitle("FEB mac5");
	gr1->SetMarkerStyle(8);
	gr1->SetMarkerSize(5);
	gr1->SetMarkerColor(kRed);
	gr1->GetYaxis()->SetNdivisions(10,1,0);
	gr2->SetMarkerStyle(8);
        gr2->SetMarkerSize(3);
        gr2->SetMarkerColor(kBlue);
        gr3->SetMarkerStyle(8);
        gr3->SetMarkerSize(1);
        gr3->SetMarkerColor(kGreen);

	new TCanvas();
	gr1->Draw("AP");
	gr2->Draw("sameP");
	gr3->Draw("sameP");

	gStyle->SetOptStat(0);
	hgain2->GetXaxis()->SetTitle("gain [ADC/PE]");
	hgain2->GetYaxis()->SetTitle("FEB mac5");
	hgain2->GetYaxis()->SetNdivisions(10,1,0);
	hgain2->SetMarkerSize(2);
	new TCanvas();
	hgain2->Draw("colztext");
}

void PlotFEBRate(){

	float ntot_in=0., ntot_out=0.;
	map<int,int> febcounts;

	TH2F* hratein = new TH2F("hFebRateIn","Trigger Fraction: Inner Layer",2,0,2,2,0,2);
	TH2F* hrateout = new TH2F("hFebRateOut","Trigger Fraction: Outer Layer",2,0,2,2,0,2);

	TH2F* hfebsin = new TH2F("hFebsIn","FEB Positions Facing South",2,0,2,2,0,2);
	TH2F* hfebsout = new TH2F("hFebsOut","FEB Positions Facing South",2,0,2,2,0,2);
	hfebsin->SetBinContent(1,1,7);
	hfebsin->SetBinContent(1,2,6);
	hfebsin->SetBinContent(2,1,3);
	hfebsin->SetBinContent(2,2,1);
	hfebsout->SetBinContent(1,1,9);
	hfebsout->SetBinContent(1,2,8);
	hfebsout->SetBinContent(2,1,5);
	hfebsout->SetBinContent(2,2,4);
	hfebsin->SetMarkerSize(5);
	hfebsout->SetMarkerSize(5);

	for(int i=0; i<a.fMacs.size(); i++) {
		int mac = a.fMacs[i];
	        string name = "spectra";
	        if(mac<10) name+="00";
        	else if(mac<100) name+="0";
	        name+=to_string(mac);	
	        TList* l = (TList*)a.fList->FindObject(name.c_str());
		TH1F* h = (TH1F*)l->At(0);
		febcounts[mac] = h->Integral();
	}

	ntot_in = febcounts[1]+febcounts[3]+febcounts[6]+febcounts[7];
	ntot_out = febcounts[4]+febcounts[5]+febcounts[8]+febcounts[9];

	cout << "mac -> tot (frac)" << endl;
	cout << " in (" << ntot_in << " tot): " << endl;
	cout << "   - 1 -> " << febcounts[1] << "(" << febcounts[1]/ntot_in << ")" << endl;
	cout << "   - 3 -> " << febcounts[3] << "(" << febcounts[3]/ntot_in << ")" << endl;
	cout << "   - 6 -> " << febcounts[6] << "(" << febcounts[6]/ntot_in << ")" << endl;
	cout << "   - 7 -> " << febcounts[7] << "(" << febcounts[7]/ntot_in << ")" << endl;
        cout << " out (" << ntot_out << " tot): " << endl;
        cout << "   - 4 -> " << febcounts[4] << "(" << febcounts[4]/ntot_out << ")" << endl;
        cout << "   - 5 -> " << febcounts[5] << "(" << febcounts[5]/ntot_out << ")" << endl;
        cout << "   - 8 -> " << febcounts[8] << "(" << febcounts[8]/ntot_out << ")" << endl;
        cout << "   - 9 -> " << febcounts[9] << "(" << febcounts[9]/ntot_out << ")" << endl;	

	hratein->SetBinContent(1,1,febcounts[7]/ntot_in);
	hratein->SetBinContent(1,2,febcounts[6]/ntot_in);
	hratein->SetBinContent(2,1,febcounts[3]/ntot_in);
	hratein->SetBinContent(2,2,febcounts[1]/ntot_in);
	hrateout->SetBinContent(1,1,febcounts[9]/ntot_out);
	hrateout->SetBinContent(1,2,febcounts[8]/ntot_out);
	hrateout->SetBinContent(2,1,febcounts[5]/ntot_out);
	hrateout->SetBinContent(2,2,febcounts[4]/ntot_out);

	hratein->GetXaxis()->SetBinLabel(1,"east");
	hratein->GetXaxis()->SetBinLabel(2,"west");
	hratein->GetYaxis()->SetBinLabel(1,"pit");
	hratein->GetYaxis()->SetBinLabel(2,"mezz");
	hratein->GetXaxis()->SetLabelSize(0.08);
	hratein->GetYaxis()->SetLabelSize(0.08);

        hrateout->GetXaxis()->SetBinLabel(1,"east");
        hrateout->GetXaxis()->SetBinLabel(2,"west");
        hrateout->GetYaxis()->SetBinLabel(1,"pit");
        hrateout->GetYaxis()->SetBinLabel(2,"mezz");
        hrateout->GetXaxis()->SetLabelSize(0.08);
        hrateout->GetYaxis()->SetLabelSize(0.08);

	gStyle->SetOptStat(0);

	new TCanvas();
	hratein->Draw("colz");
	hfebsin->Draw("sameTEXT");

	new TCanvas();
	hrateout->Draw("colz");
	hfebsout->Draw("sameTEXT");

}

void PlotThreshMult(int mac){


}	
