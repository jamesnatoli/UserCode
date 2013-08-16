//////////////////////////////////////////////////////////////////////////////
// Assumptions:
//
// - MC files are in ../root directory, and their names have following pattern:
//     "sample_name" + user-provided string
//   For example, tt_staco, zmm_staco -> user will provide string "_staco"
//
// - Users need a cutflow histogram which contains the total number of
//   MC events which were used to fill the MC histograms, before any
//   cut, to correctly normalize with respect to luminosity (and to
//   correctly normalize each MC sample w.r.t. each other).
//   The cutflow histogram name is "Events", in current macro
//
// - Canvases will be saved in ../figs directory
// - Possible interventions:
//   - change stack order: search for "stackingOrder" vector
//   - data vs lumi normalization: "Renormalization of data vs MC"
//   - stack Y axis max: "Special requests for Y axis range!"
//   - definition of Poisson errors (!): GetPoissonizedGraph
//
// Usage:
//
// - compile the macro in root:
//  [] .L doDataMCPlots.C++
// - launch functions for 1D plotting:
//  [] doDataMCPlots1D("_8TeV","myDataFile.root",19600,1,0,0)
//
// Support:
//
//  alberto.belloni at cern.ch
//
//////////////////////////////////////////////////////////////////////////////
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TFile.h>
#include <TMath.h>
#include <TLine.h>
#include <TList.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>

#include "globals.h"

using namespace std;

double error_poisson_up(double data);
double error_poisson_down(double data);

double compatibilityTest(TH1F* hist, THStack* hstack);

TGraphAsymmErrors* GetPoissonizedGraph(TH1F* histo);

void myText(Double_t x,Double_t y,Color_t color,const char *text);
void myCmsLabel(Double_t x,Double_t y, int type);
void myOtherText(Double_t x,Double_t y,double tsize,const char *text);

void doDataMCPlots1D(const char*, const char*,double, 
		     /*double,*/ int, int, int);

multimap<int,int> mapSamplesToCategories(multimap<int,int>&, int);


//////////////////////////////////////////////////////////////////////////////
// Definitions for Collective samples
//////////////////////////////////////////////////////////////////////////////
// Wgg, for example, combines W->e+gg and W->tau+gg
// WWg combines W(e)W(e/mu/tau)g and W(tau)W(e/mu/tau)g, 6 samples!
//////////////////////////////////////////////////////////////////////////////

const int MCGRP = 25;

enum {WAAEF = 0,  // WAA signal MC, FSR/ISR-enhanced
      WAAEI,
      WAAMF,
      WAAMI,
      ZAAE,      // Z + AA (no Drell-Yan)
      ZAAM,
      TOP,        // all top contributions
      AAJE,       // AA + jets
      AAJM,
      VVJE,       // diboson + jets
      VVJM,
      WWAE,       // WWA
      WWAM,
      WAJE,       // WA + jets
      WAJM,
      ZAJE,       // Z + A + jets (including Drell-Yan)
      ZAJM,
      ZJJE,       // Z + jets (including Drell-Yan)
      ZJJM,
      WJJE,       // W + jets (well, at least 2!)
      WJJM,
      AJJE,       // A + jets (one fakes other photon, the other the lepton)
      AJJM,
      VVVE,       // triboson
      VVVM
};

const char* legend[MCGRP] = {
  "W(e#nu)#gamma#gamma FSR",
  "W(e#nu)#gamma#gamma ISR",
  "W(#mu#nu)#gamma#gamma FSR",
  "W(#mu#nu)#gamma#gamma ISR",
  "Z#gamma#gamma",
  "Z#gamma#gamma",
  "top",
  "#gamma#gamma + jets",
  "#gamma#gamma + jets",
  "VV + jets",
  "VV + jets",
  "WW#gamma",
  "WW#gamma",
  "W#gamma + jets",
  "W#gamma + jets",
  "Drell-Yan + #gamma + jets",
  "Drell-Yan + #gamma + jets",
  "Drell-Yan + jets",
  "Drell-Yan + jets",
  "W + jets",
  "W + jets",
  "#gamma + jets",
  "#gamma + jets",
  "VVV",
  "VVV"
};
  
const int mccolor[MCGRP] = {
  10, // transparent?
  10,
  10,
  10,
  kAzure-9,
  kAzure-9,
  kGreen+1,
  kOrange+1,
  kOrange+1,
  kYellow-9,
  kYellow-9,
  kViolet,
  kViolet,
  kViolet-2,
  kViolet-2,
  kAzure-11,
  kAzure-11,
  kAzure-13,
  kAzure-13,
  kViolet-4,
  kViolet-4,
  kOrange-1,
  kOrange-1,
  kYellow-7,
  kYellow-7
};
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
void doDataMCPlots1D(const char* mcTag,
		     const char* dtFile,
		     double lumi,
		     /*double QCDSCALE,*/
		     int CmsLabel,
		     int decayChannel,
		     int doDataNorm) {

  // Group-sample vector: it defines how the individual MC files
  // are to be added up
  multimap<int, int> allowedSample;
  mapSamplesToCategories(allowedSample,decayChannel);

  // Open files /////////////////////////////////////////////////////////////
  TFile *dt = TFile::Open(dtFile);
  TFile *mc[MCID];
  for (multimap<int,int>::iterator it=allowedSample.begin();
       it!=allowedSample.end();
       ++it)
    mc[it->second] = TFile::Open(Form("root/%s%s.root",
				     mcsample[it->second],mcTag));
    

  // Open and book histograms ///////////////////////////////////////////////
  const int N1DHISTS = 7;

  TH1F *h_dt1d[N1DHISTS];
  TH1F *h_mc1d[N1DHISTS][MCID];

  // Cutflow histogram for normalizations
  TH1F* h_cutflow[MCID];


  // Output file names
  const char* file1dname[N1DHISTS]= {
    "PtLead",  //  0
    "PtSubL",  //  1
    "MEt",     //  2
    "MLAA",    //  3
    "PtLAA",   //  4
    "MtLAAMEt",//  5
    "PtLAAMet" //  6
  };

  // Histogram names
  const char* hist1dname[N1DHISTS] ={
    "PtLead",  //  0
    "PtSubL",  //  1
    "MEt",     //  2
    "MLAA",    //  3
    "PtLAA",   //  4
    "MtLAAMEt",//  5
    "PtLAAMEt" //  6
  };

  // List of Log scale plots:
  vector<int> logScale;
  logScale.push_back(0);
  logScale.push_back(1);
  logScale.push_back(2);
  logScale.push_back(3);
  logScale.push_back(4);
  logScale.push_back(5);
  logScale.push_back(6);

  // A couple of strings that will be used later to identify the plots
  const char* leptons[2] = {"_e","_m"};
  const char* titlept[2] = {"e","#mu"};

  // The following code may crash during run; it compiled, but
  // I am not sure that the Form(..) will work inside an array
  // initialization...
  const char* xtitle1d[N1DHISTS]= {
    "p_{T}(leading lepton) [GeV]",
    "p_{T}(second lepton) [GeV]",
    "E_{T}^{miss} [GeV]",
    Form("m(%s#gamma#gamma) [GeV]",titlept[decayChannel]),
    Form("p_{T}(%s#gamma#gamma) [GeV]",titlept[decayChannel]),
    Form("m_{T}(%s,E_{T}^{miss}) [GeV]",titlept[decayChannel]),
    Form("p_{T}(%s,E_{T}^{miss}) [GeV]",titlept[decayChannel])
  };
  
  // Actually, just the units!
  const char* ytitle1d[N1DHISTS]= {
    "GeV",
    "GeV",
    "GeV",
    "GeV",
    "GeV",
    "GeV",
    "GeV"
  };
  
  // Book histograms
  for (int i=0;i<N1DHISTS;++i) {
    h_dt1d[i] = (TH1F*)dt->Get(Form("%s%s",
				    hist1dname[i],
				    leptons[decayChannel]));
    cout << "DT: " << h_dt1d[i] << " " 
	 << Form("%s%s",hist1dname[i],leptons[decayChannel]) 
	 << " " << (h_dt1d[i]?h_dt1d[i]->GetNbinsX():0) << endl;
    
    for (multimap<int,int>::iterator it=allowedSample.begin();
	 it!=allowedSample.end();
	 ++it) {
      
      h_mc1d[i][it->second] =
	(TH1F*)mc[it->second]->Get(Form("%s%s",
					hist1dname[i],
					leptons[decayChannel]));
      cout << "MC" << it->second << ": " <<h_mc1d[i][it->second] << " " 
	   << Form("%s%s",hist1dname[i],leptons[decayChannel]) 
	   << " " 
	   << (h_mc1d[i][it->second]?h_mc1d[i][it->second]->GetNbinsX():0)
	   << endl;
      
      if(!h_mc1d[i][it->second])
	h_mc1d[i][it->second] = new TH1F();


      h_mc1d[i][it->second]->
	SetName(Form("%s%s_%s",
		     mcsample[it->second],
		     hist1dname[i],
		     leptons[decayChannel]));
      
      // Get cutflow histogram on first loop
      if (i==0)
	h_cutflow[it->second] =
	  (TH1F*)mc[it->second]->Get("Events");

    }
  }

  // End of open and book histograms ////////////////////////////////////////

  // Scale+Add MC histograms ////////////////////////////////////////////////
  for (int i=0;i<N1DHISTS;++i)
    for (multimap<int,int>::iterator it=allowedSample.begin();
	 it!=allowedSample.end();
	 ++it) {
      
      // Normalize to 1nb^-1
      if (h_cutflow[it->second]->GetBinContent(1)>0)
	h_mc1d[i][it->second]->
	  Scale(scaleFactor[it->second]/
		h_cutflow[it->second]->GetBinContent(1));

      if (i==0)
	cout << endl << "Expected in " << lumi
	     << " pb^-1, sample " << mcsample[it->second] 
	     << ": " << lumi*h_mc1d[0][it->second]->
	  Integral(1,h_mc1d[0][it->second]->GetNbinsX())
	     << endl;
      
      // Add, for each key, all histograms to first one in category
      if (it->second!=allowedSample.find(it->first)->second &&
	  h_mc1d[i][it->second]->GetEntries()!=0)
	h_mc1d[i][allowedSample.find(it->first)->second]
	  ->Add(h_mc1d[i][it->second]);
    }
  // End of scale+add MC histograms /////////////////////////////////////////
  

  // Renormalization of data vs MC /////////////////////////////////////////
  cout << "Renormalization of data vs MC" << endl;
  for (int i=0;i<N1DHISTS;++i) {

    // Calculate scale w.r.t. data entries
    double scale = 0.;
    
    // Loop only on first histogram in each category
    for (multimap<int,int>::iterator it=allowedSample.begin();
	 it!=allowedSample.end();
	 ++it) {
      if (it->second!=allowedSample.find(it->first)->second)
	continue;
      
      scale+=h_mc1d[i][it->second]->Integral();
    }
    
    scale  = 1./scale; // scale was sum of MC
    scale *=  h_dt1d[i]->Integral();
    
    cout << hist1dname[i] << " " << scale << " " << lumi << endl;
    
    if (doDataNorm!=0) {
      
      // Normalize MC samples to data
      for (multimap<int,int>::iterator it=allowedSample.begin();
	   it!=allowedSample.end();
	   ++it) {
	if (it->second!=allowedSample.find(it->first)->second)
	  continue;
	
	h_mc1d[i][it->second]->Scale(scale);
      }
    }
    else {

      // Normalize MC samples to luminosity
      for (multimap<int,int>::iterator it=allowedSample.begin();
	   it!=allowedSample.end();
	   ++it) {
	if (it->second!=allowedSample.find(it->first)->second)
	  continue;
	
	h_mc1d[i][it->second]->Scale(lumi);
      }
    }
  }
  
  // Print expected number of MC events, after all normalizations
  for (multimap<int,int>::iterator it=allowedSample.begin();
       it!=allowedSample.end();
       ++it) {
    if (it->second!=allowedSample.find(it->first)->second)
      continue;
    
    cout << endl << "Expected " << legend[it->first] 
	 << ": " << h_mc1d[0][it->second]->
      Integral(1,h_mc1d[0][it->second]->GetNbinsX())
	 << endl;
  }
  cout << endl;
  // End of renormalization of data vs MC ////////////////////////////////

  // Style up histograms ////////////////////////////////////////////////////
  cout << "Styling 1D" << endl;

  for (int i=0;i<N1DHISTS;++i) {
    h_dt1d[i]->SetLineColor(kBlack);
    h_dt1d[i]->SetMarkerColor(kBlack);
    h_dt1d[i]->SetFillColor(0);
    h_dt1d[i]->SetFillStyle(0);

    // Loop only on first histogram in each category
    for (multimap<int,int>::iterator it=allowedSample.begin();
	 it!=allowedSample.end();
	 ++it) {
      if (it->second!=allowedSample.find(it->first)->second)
	continue;
     
      h_mc1d[i][it->second]->SetLineColor(kBlack);
      h_mc1d[i][it->second]->SetMarkerColor(mccolor[it->first]);
      h_mc1d[i][it->second]->SetFillColor(mccolor[it->first]);
      h_mc1d[i][it->second]->SetFillStyle(1001);
      
    }
  }
  // End of style up histograms /////////////////////////////////////////////

  // Add titles to axes /////////////////////////////////////////////////////
  cout << "Titles and axes..." << endl;
  for (int i=0;i<N1DHISTS;++i) {
    h_dt1d[i]->GetXaxis()->SetTitle(Form("%s",
					 xtitle1d[i]));
    h_dt1d[i]->GetXaxis()->SetTitleOffset(1.3);
    if (strstr(hist1dname[i],"Eta")== NULL)
      h_dt1d[i]->GetYaxis()->SetTitle(Form("Entries / %3.1f %s",
					   h_dt1d[i]->GetBinWidth(1),
					   ytitle1d[i]));
    else
      h_dt1d[i]->GetYaxis()->SetTitle("Entries per bin");

    if (find(logScale.begin(),logScale.end(),i)!=logScale.end()) {
      h_dt1d[i]->GetYaxis()->SetTitleOffset(1.3);
    }
    else {
      h_dt1d[i]->GetYaxis()->SetTitleOffset(1.7);
      if (i==18||i==5||i==6||i==7)
	h_dt1d[i]->GetYaxis()->SetTitleOffset(1.9);
      if (i==10||i==8||i==23)
	h_dt1d[i]->GetYaxis()->SetTitleOffset(2.0);
    }

    // Loop only on first histogram in each category
    for (multimap<int,int>::iterator it=allowedSample.begin();
	 it!=allowedSample.end();
	 ++it) {
      if (it->second!=allowedSample.find(it->first)->second)
	continue;

      h_mc1d[i][it->second]->GetXaxis()->SetTitle(Form("%s",
						       xtitle1d[i]));
      h_mc1d[i][it->second]->GetXaxis()->SetTitleOffset(1.3);
      if (strstr(hist1dname[i],"Eta")== NULL)
	h_mc1d[i][it->second]->GetYaxis()->
	  SetTitle(Form("Entries / %3.1f %s",
			h_dt1d[i]->GetBinWidth(1),
			ytitle1d[i]));
      else
	h_mc1d[i][it->second]->GetYaxis()->SetTitle("Entries per bin");
      if (find(logScale.begin(),logScale.end(),i)!=logScale.end()) {
	h_mc1d[i][it->second]->GetYaxis()->SetTitleOffset(1.3);
      }
      else {
	h_mc1d[i][it->second]->GetYaxis()->SetTitleOffset(1.7);
	if (i==18||i==5||i==6||i==7)
	  h_mc1d[i][it->second]->GetYaxis()->SetTitleOffset(1.9);
	if (i==10||i==8||i==23)
	  h_mc1d[i][it->second]->GetYaxis()->SetTitleOffset(2.0);
      }
    }
  }
  // End of add titles to axes /////////////////////////////////////////////

  // Stack histograms //////////////////////////////////////////////////////
  THStack *s_mc1d[N1DHISTS];
  
  //////////////////////////////////////////////////////////////////////////
  // Here we define the stacking order
  //////////////////////////////////////////////////////////////////////////
  vector<int> stackingOrder;
  if      (decayChannel==0) { // e
    stackingOrder.push_back(AJJE);
    stackingOrder.push_back(WJJE);
    stackingOrder.push_back(ZJJE);
    stackingOrder.push_back(ZAJE); 
    stackingOrder.push_back(WAJE);
    stackingOrder.push_back(WWAE);
    stackingOrder.push_back(VVJE);
    stackingOrder.push_back(AAJE);
    stackingOrder.push_back(TOP);
    stackingOrder.push_back(ZAAE);
    stackingOrder.push_back(WAAEF);
    stackingOrder.push_back(WAAEI);
  }
  else if (decayChannel==1) { // mu
    stackingOrder.push_back(AJJM);
    stackingOrder.push_back(WJJM);
    stackingOrder.push_back(ZJJM);
    stackingOrder.push_back(ZAJM); 
    stackingOrder.push_back(WAJM);
    stackingOrder.push_back(WWAM);
    stackingOrder.push_back(VVJM);
    stackingOrder.push_back(AAJM);
    stackingOrder.push_back(TOP);
    stackingOrder.push_back(ZAAM);
    stackingOrder.push_back(WAAMF);
    stackingOrder.push_back(WAAMI);
  }
  
  for (int i=0;i<N1DHISTS;++i) {
    s_mc1d[i] = new THStack(hist1dname[i],hist1dname[i]);
    s_mc1d[i]->SetName(Form("stack_%s",hist1dname[i]));
    for (unsigned int j=0;j<stackingOrder.size();++j) {
      
      s_mc1d[i]->
	Add(h_mc1d[i][allowedSample.find(stackingOrder.at(j))->second]);
    }
  }
  // End of stack histograms ///////////////////////////////////////////////

  
  // Legends!
  cout << "Doing Legends" << endl;
  TLegend *leg1d[N1DHISTS];

  // Will make them in same position, put "if(i==..)" to change them
  for (int i=0;i<N1DHISTS;++i) {
    if   (i==6||i==7||i==25)
      leg1d[i] = new TLegend(0.20,0.64,0.56,0.915);
    else
      leg1d[i] = new TLegend(0.56,0.64,0.92,0.915);

    leg1d[i]->SetName(Form("leg_%s",h_dt1d[i]->GetName()));
    leg1d[i]->SetTextSize(0.03);
    leg1d[i]->SetFillColor(0);
    leg1d[i]->AddEntry(h_dt1d[i],"Data 2010 (#sqrt{s} = 7 TeV)","lp");

    
    // Always MC signal first!
    ////////////////////////////////////////////////////////////////////////
    // Bold assumption: signal is always going to be first sample
    // in definition 
    ////////////////////////////////////////////////////////////////////////
    leg1d[i]->AddEntry(h_mc1d[i][allowedSample.begin()->second],
		       legend[allowedSample.begin()->first],"f");

    for (unsigned int j=stackingOrder.size();j>0;--j) {
      if (stackingOrder.at(j-1)!=allowedSample.begin()->first)
	leg1d[i]->AddEntry
	  (h_mc1d[i][allowedSample.find(stackingOrder.at(j-1))->second],
	   legend[stackingOrder.at(j-1)],
	   "f");

    }
  }

  // Canvases, at last!
  TCanvas *canv1d[N1DHISTS];

  for (int i=0;i<N1DHISTS;++i) {

    canv1d[i] = new TCanvas(Form("c_%s%s",file1dname[i],
				 leptons[decayChannel]),
			    Form("c_%s%s",file1dname[i],
				 leptons[decayChannel]),
			    500,500);
    canv1d[i]->SetLeftMargin(0.18);
    canv1d[i]->SetRightMargin(0.05);

    // Set log scale for selected histograms
    if (find(logScale.begin(),logScale.end(),i)!=logScale.end()) {
      canv1d[i]->SetLogy();
      s_mc1d[i]->SetMaximum(h_dt1d[i]->GetMaximum()*5.e1);
      s_mc1d[i]->SetMinimum(5e-1);
    }      
    else {
      s_mc1d[i]->SetMaximum(h_dt1d[i]->GetMaximum()*1.5);
      if (i>=5&&i<=10) // eta and phi plots
	s_mc1d[i]->SetMaximum(h_dt1d[i]->GetMaximum()*1.8);
      s_mc1d[i]->SetMinimum(1.e-2);
    }

    s_mc1d[i]->Draw();

    // Set ranges of axes
    double x_max = 
      find(logScale.begin(),logScale.end(),i)!=logScale.end()?
      200.:120.;

    if (i==0||i==1)
      s_mc1d[i]->GetXaxis()->SetRangeUser(25.,x_max); 

    s_mc1d[i]->Draw("hist");

    TGraphAsymmErrors* g_data = GetPoissonizedGraph(h_dt1d[i]);
    g_data->SetMarkerStyle(20);
    g_data->SetMarkerSize(0.9);
    g_data->SetLineWidth(2);
    g_data->Draw("z,same");
    h_dt1d[i]->SetMarkerStyle(20);
    h_dt1d[i]->SetMarkerSize(0.9);
    h_dt1d[i]->Draw("p,same");
    h_dt1d[i]->Draw("axis,same");

    leg1d[i]->Draw("same");

    // Labels! 
    if      (i==16) {
      myCmsLabel(0.57,0.56,CmsLabel);
    }
    else if (i==13) {
      myCmsLabel(0.2,0.88,CmsLabel);
    }
    else if (i==3) {
      myCmsLabel(0.2,0.88,CmsLabel);
    }
    else if (i==20) {
      myCmsLabel(0.57,0.56,CmsLabel);
    }

    if      (i<5)
      myText(0.61, 0.54,1,Form("#int L dt = %.0f pb^{-1}",lumi));
    else if (i==6||i==7||i==25)
      myText(0.61, 0.78,1,Form("#int L dt = %.0f pb^{-1}",lumi));
    else if (i==24)
      myText(0.20, 0.78,1,Form("#int L dt = %.0f pb^{-1}",lumi));
    else if (i>11)
      myText(0.61, 0.54,1,Form("#int L dt = %.0f pb^{-1}",lumi));
    else
      myText(0.20, 0.78,1,Form("#int L dt = %.0f pb^{-1}",lumi));

    // Vertical dashed line in pT presel plot, marking the cut:
    //TLine* line = new TLine();
    //line->SetLineStyle(kDashed);
    //line->SetLineWidth(2);
    //if (i==19) {
    //  line->SetX1(20);
    //  line->SetY1(0);
    //  line->SetX2(20);
    //  line->SetY2(1.e7); // or some other high number!
    //  line->DrawClone();
    //  canv1d[i]->Update();
    //}

    if(strstr(hist1dname[i],"Pos")!=NULL) {
      if (i==6||i==7||i==25)
	myOtherText(0.40,0.73,0.12,"#mu^{+}");
      else
	myOtherText(0.77,0.73,0.12,"#mu^{+}");
    }
    if(strstr(hist1dname[i],"Neg")!=NULL) {
      if (i==6||i==7)
	myOtherText(0.40,0.73,0.12,"#mu^{-}");
      else
	myOtherText(0.77,0.73,0.12,"#mu^{-}");
    }

    canv1d[i]->Print(Form("figs/%s%s.C",
			  file1dname[i],
			  leptons[decayChannel]));
    canv1d[i]->Print(Form("figs/%s%s.png",
			  file1dname[i],
			  leptons[decayChannel]));
    canv1d[i]->Print(Form("figs/%s%s.pdf",
			  file1dname[i],
			  leptons[decayChannel]));
    canv1d[i]->Print(Form("figs/%s%s.eps",
			  file1dname[i],
			  leptons[decayChannel]));


    cout << " Probability " << file1dname[i] << leptons[decayChannel]
	 << " : "
	 << compatibilityTest(h_dt1d[i],s_mc1d[i]) << endl;

    //delete canv1d[i];
  }

}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* GetPoissonizedGraph(TH1F* histo) {

  TGraphAsymmErrors* graph = new TGraphAsymmErrors();
  graph->SetName(Form("g_%s",histo->GetName()));

  int j=1;
  for (int i=1;i<=histo->GetNbinsX();++i) {

    if (histo->GetBinContent(i)!=0) {
  
      graph->SetPoint(j,histo->GetBinCenter(i),histo->GetBinContent(i));

      graph->SetPointError(j,
			   histo->GetBinWidth(i)/2.,
			   histo->GetBinWidth(i)/2.,
			   error_poisson_down(histo->GetBinContent(i)),
			   error_poisson_up(histo->GetBinContent(i)));

      ++j;
    }
  }
  return graph;

}

double error_poisson_up(double data) {
  double y1 = data + 1.0;
  double d = 1.0 - 1.0/(9.0*y1) + 1.0/(3*TMath::Sqrt(y1));
  return y1*d*d*d-data;
}

double error_poisson_down(double data) {
  double y = data;
  if (y == 0.0) return 0.0;
  double d = 1.0 - 1.0/(9.0*y) - 1.0/(3.0*TMath::Sqrt(y));
  return data-y*d*d*d;
}
//////////////////////////////////////////////////////////////////////////////


// New function for comparing errors

double compatibilityTest(TH1F* hist, THStack* hstack) {

  // A bit of a disaster: need to loop on list of histograms in stack

  TList* stack_hists = hstack->GetHists();

  int ndf = 0;
  double chisq = 0.;
  double h_entries = 0.; // entries in histo
  double s_entries = 0.; // entries in stack
  double s_sqerror = 0.;
  for (int i = 1; i<=hist->GetNbinsX(); ++i) {
    
    h_entries += hist->GetBinContent(i);
    
    TIter iter(stack_hists->MakeIterator());
    while (TH1F* stack_hist = dynamic_cast<TH1F*>(iter())) {
      s_entries += stack_hist->GetBinContent(i);
      s_sqerror += (stack_hist->GetBinError(i)*
		    stack_hist->GetBinError(i));
    }
    
    if (h_entries<20 && i<hist->GetNbinsX()) 
      continue;
    else {
      ndf++;
      
      chisq += (h_entries-s_entries)*(h_entries-s_entries)/
	(s_sqerror + h_entries);
      

      h_entries = s_entries = s_sqerror = 0.;
    }
    
  } 
  cout << "Chisq, ndf: " << chisq << ", " << ndf << endl; 
  return TMath::Prob(chisq,ndf);
}


void myText(Double_t x,Double_t y,Color_t color,const char *text)
{

  TLatex l;
  l.SetNDC();
  l.SetTextColor(color);
  l.SetTextSize(0.05);
  l.DrawLatex(x,y,text);
}

void myOtherText(Double_t x,Double_t y,double tsize,const char *text)
{
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.DrawLatex(x,y,text);
}

void myCmsLabel(Double_t x,Double_t y, int type) {

  // CMS Label types:
  // 0: no label
  // 1: CMS

  if (type==0)
    return;
  
  TLatex l;
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextSize(0.04);
  l.DrawLatex(x,y,"CMS");
  l.SetTextFont(42);
  l.SetTextSize(0.04);
  l.DrawLatex(x+0.13,y,"Preliminary");

}


multimap<int,int> mapSamplesToCategories(multimap<int,int>& chMap, int decay)
{

  if      (decay==0) { // ee
    // Signal
    chMap.insert(pair<int,int>(WAAEI,qqWAAEI));
    chMap.insert(pair<int,int>(WAAEI,qqWAATI));
    chMap.insert(pair<int,int>(WAAEF,qqWAAEF));
    chMap.insert(pair<int,int>(WAAEF,qqWAATF));
    // Z + AA
    chMap.insert(pair<int,int>(ZAAE,qqZAAE));
    chMap.insert(pair<int,int>(ZAAE,qqZAAT));
    // Top
    chMap.insert(pair<int,int>(TOP,TTBAR));
    chMap.insert(pair<int,int>(TOP,TTAA));
    chMap.insert(pair<int,int>(TOP,TtE));
    chMap.insert(pair<int,int>(TOP,TtT));
    chMap.insert(pair<int,int>(TOP,TsE));
    chMap.insert(pair<int,int>(TOP,TsT));
    chMap.insert(pair<int,int>(TOP,WT));
    // AA + jets
    chMap.insert(pair<int,int>(AAJE,qqAA_10TO25));
    chMap.insert(pair<int,int>(AAJE,qqAA_25TO250));
    chMap.insert(pair<int,int>(AAJE,qqAA_250TOINF));
    chMap.insert(pair<int,int>(AAJE,DIPHOTONJETS));
    // VV + jets
    chMap.insert(pair<int,int>(VVJE,qqWZ2L2Q));
    chMap.insert(pair<int,int>(VVJE,qqWZ3L1N));
    chMap.insert(pair<int,int>(VVJE,qqZZ4E));
    chMap.insert(pair<int,int>(VVJE,qqZZ2E2M));
    chMap.insert(pair<int,int>(VVJE,qqZZ2E2T));
    chMap.insert(pair<int,int>(VVJE,qqZZ2M2T));
    chMap.insert(pair<int,int>(VVJE,qqZZ4T));
    chMap.insert(pair<int,int>(VVJE,qqWW2L2N));
    // WWA
    chMap.insert(pair<int,int>(WWAE,qqWWA));
    // WA + jets
    chMap.insert(pair<int,int>(WAJE,qqWA));
    // ZA + jets
    chMap.insert(pair<int,int>(ZAJE,qqZA));
    // Z + jets
    chMap.insert(pair<int,int>(ZJJE,qqZ)); // inclusive Z sample? Alpgen NP1, NP2?
    // W + jets
    chMap.insert(pair<int,int>(WJJE,qqW)); // inclusive W sample? Alpgen NP1, NP2?
    // A + jets
    chMap.insert(pair<int,int>(AJJE,qqA)); // inclusive A sample? Alpgen NP1, NP2?
    // VVV
    chMap.insert(pair<int,int>(VVVE,qqZZZ));
    chMap.insert(pair<int,int>(VVVE,qqWZZ));
    chMap.insert(pair<int,int>(VVVE,qqWWZ));
    chMap.insert(pair<int,int>(VVVE,qqWWW));
  }
  else if (decay==1) { // mumu
    // Signal
    chMap.insert(pair<int,int>(WAAMI,qqWAAMI));
    chMap.insert(pair<int,int>(WAAMI,qqWAATI));
    chMap.insert(pair<int,int>(WAAMF,qqWAAMF));
    chMap.insert(pair<int,int>(WAAMF,qqWAATF));
    // Z + AA
    chMap.insert(pair<int,int>(ZAAM,qqZAAM));
    chMap.insert(pair<int,int>(ZAAM,qqZAAT));
    // Top
    chMap.insert(pair<int,int>(TOP,TTBAR));
    chMap.insert(pair<int,int>(TOP,TTAA));
    chMap.insert(pair<int,int>(TOP,TtM));
    chMap.insert(pair<int,int>(TOP,TtT));
    chMap.insert(pair<int,int>(TOP,TsM));
    chMap.insert(pair<int,int>(TOP,TsT));
    chMap.insert(pair<int,int>(TOP,WT));
    // AA + jets
    chMap.insert(pair<int,int>(AAJM,qqAA_10TO25));
    chMap.insert(pair<int,int>(AAJM,qqAA_25TO250));
    chMap.insert(pair<int,int>(AAJM,qqAA_250TOINF));
    chMap.insert(pair<int,int>(AAJM,DIPHOTONJETS));
    // VV + jets
    chMap.insert(pair<int,int>(VVJM,qqWZ2L2Q));
    chMap.insert(pair<int,int>(VVJM,qqWZ3L1N));
    chMap.insert(pair<int,int>(VVJM,qqZZ4M));
    chMap.insert(pair<int,int>(VVJM,qqZZ2E2M));
    chMap.insert(pair<int,int>(VVJM,qqZZ2E2T));
    chMap.insert(pair<int,int>(VVJM,qqZZ2M2T));
    chMap.insert(pair<int,int>(VVJM,qqZZ4T));
    chMap.insert(pair<int,int>(VVJM,qqWW2L2N));
    // WWA
    chMap.insert(pair<int,int>(WWAM,qqWWA));
    // WA + jets
    chMap.insert(pair<int,int>(WAJM,qqWA));
    // ZA + jets
    chMap.insert(pair<int,int>(ZAJM,qqZA));
    // Z + jets
    chMap.insert(pair<int,int>(ZJJM,qqZ)); // inclusive Z sample? Alpgen NP1, NP2?
    // W + jets
    chMap.insert(pair<int,int>(WJJM,qqW)); // inclusive W sample? Alpgen NP1, NP2?
    // A + jets
    chMap.insert(pair<int,int>(AJJM,qqA)); // inclusive A sample? Alpgen NP1, NP2?
    // VVV
    chMap.insert(pair<int,int>(VVVM,qqZZZ));
    chMap.insert(pair<int,int>(VVVM,qqWZZ));
    chMap.insert(pair<int,int>(VVVM,qqWWZ));
    chMap.insert(pair<int,int>(VVVM,qqWWW));
  }
  else
    cout << endl << "Returning empty map!" << endl << endl;

  return chMap;
}
