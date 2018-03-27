#include "slimanalysis.h"
#include "helper_funcs.C"

void doAlignmentPlots(bool debug, const char* dir) {

  TChain* chain = new TChain("slim");
  chain->Add(Form("%s/*_slim.root",dir));

  // Plots without correction
  const int NPLANES = 4; // uhm, technically it is something different...
  TH1F* hist[NPLANES]; // There will be four histograms
  const string var[NPLANES] = {
    "xC-xA",
    "yC-yA",
    Form("xC-xA%s%f",wc_xca>0?"+":"",wc_xca),
    Form("yC-yA%s%f",wc_yca>0?"+":"",wc_yca)
  };
  const int color[NPLANES] = {
    kBlue,
    kBlue,
    kBlue,
    kBlue
  };

  double mean[NPLANES], sigma[NPLANES];
  for (int i = 0; i < NPLANES; ++i) {
    hist[i] = new TH1F(var[i].c_str(),"",240,-30,30);
    chain->Project(var[i].c_str(),var[i].c_str(),
                   Form("abs(%s)<30",var[i].c_str()));
    mean[i] = hist[i]->Fit("gaus","SQN")->GetParams()[1];
    // really not best way: why re-doing fit? I just need to
    // check if resolution of wire chambers are similar or not
    // (if different, will need to include uncertainties in tracking fit
    sigma[i] = hist[i]->Fit("gaus","SQN")->GetParams()[2];

    hist[i]->SetLineWidth(2);
    hist[i]->SetLineColor(color[i]);
    if (var[i].find('x')!=std::string::npos)
      hist[i]->GetXaxis()->SetTitle("#deltax [mm]");
    else
      hist[i]->GetXaxis()->SetTitle("#deltay [mm]");
    hist[i]->GetXaxis()->SetTitleOffset(1.2);
    hist[i]->GetYaxis()->SetTitle("Events");
    hist[i]->GetYaxis()->SetTitleOffset(1.6);
  }

  for (int i = 0; i < NPLANES; ++i)
    cout << var[i].c_str() << " mean  -> " << mean[i] << endl;
  for (int i = 0; i < NPLANES; ++i)
    cout << var[i].c_str() << " sigma -> " << sigma[i] << endl;

  TCanvas* canv[NPLANES];
  TLegend* leg[NPLANES];
  const string entry[NPLANES] = {
    "x_{C}-x_{A}",
    "y_{C}-y_{A}",
    "x_{C}-x_{A}",
    "y_{C}-y_{A}"
  };

  for (int i = 0; i < NPLANES; ++i) {
    canv[i] = new TCanvas(Form("align_%d",i),"",500,500); // makes the new TCanvases                                                            
    leg[i] = new TLegend(0.7,0.7,0.9,0.9,"","brNDC");
    leg[i]->SetTextSize(0.05);
    leg[i]->AddEntry(hist[i],entry[i].c_str(),"l");

    hist[i]->Draw();
    leg[i]->Draw("same");
    hist[i]->Draw("axis,same");
    canv[i]->Print(Form("Alignment_Plots/align_%d.png",i));
  }
}

void doMaps(bool debug, const char* dir) {
  // Identify channels we need to use
  struct channel {
    int chan;
    int ieta;
    int idepth;
    string name;
  };

  bool eff, eff_rot, effX, effY, effX_rot, effY_rot;
  effX = effY = effX_rot = effY_rot = !debug;
  eff = eff_rot = !debug;

  // Fill the Rotated Arrays
  // This fills both the Theta array, and the fiducial arrays
  fill_Rot_Array();

  // Let us define the histogram and book them
  TH2F *hist_eff[NUMCHAN]; 
  TH2F *hist_den[NUMCHAN];
  TH1F *hist_effX[NUMCHAN]; 
  TH1F *hist_denX[NUMCHAN];
  TH1F *hist_effY[NUMCHAN]; 
  TH1F *hist_denY[NUMCHAN];

  // Rotated Maps
  TH2F *hist_eff_rot[NUMCHAN];
  TH2F *hist_den_rot[NUMCHAN];
  TH1F *hist_effX_rot[NUMCHAN];
  TH1F *hist_denX_rot[NUMCHAN];
  TH1F *hist_effY_rot[NUMCHAN];
  TH1F *hist_denY_rot[NUMCHAN];
  
  if (channels.size()!= NUMCHAN) {
    cout << "Argh, something wrong!" << endl;
    return;
  }

  for (unsigned int i = 0; i < channels.size(); ++i) {
    // Horrible way to replace - with _ in the channel's name, so that we can print canvas to .C file
    // and obtain a usable macro
    hist_eff[i] = new TH2F(TString(Form("%s_eff",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                           "",350,-75,75,350,-75,75);
    hist_den[i] = new TH2F(TString(Form("%s_den",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                           "",350,-75,75,350,-75,75);
    hist_eff[i]->Sumw2();
    hist_den[i]->Sumw2();

    hist_effX[i] = new TH1F(TString(Form("%s_effX",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",350,-100,100);
    hist_denX[i] = new TH1F(TString(Form("%s_denX",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",350,-100,100);
    hist_effX[i]->Sumw2();
    hist_denX[i]->Sumw2();

    hist_effY[i] = new TH1F(TString(Form("%s_effY",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",400,-100,100);
    hist_denY[i] = new TH1F(TString(Form("%s_denY",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",400,-100,100);
    hist_effY[i]->Sumw2();
    hist_denY[i]->Sumw2();

    // Rotated Plots
    // 2D 
    hist_eff_rot[i] = new TH2F(TString(Form("%s_eff_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                           "",350,-75,75,350,-75,75);
    hist_den_rot[i] = new TH2F(TString(Form("%s_den_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                           "",350,-75,75,350,-75,75);
    hist_eff_rot[i]->Sumw2();
    hist_den_rot[i]->Sumw2();
    // X Efficiency
    hist_effX_rot[i] = new TH1F(TString(Form("%s_effX_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",350,-100,100);
    hist_denX_rot[i] = new TH1F(TString(Form("%s_denX_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",350,-100,100);
    hist_effX_rot[i]->Sumw2();
    hist_denX_rot[i]->Sumw2();
    // Y Efficiency
    hist_effY_rot[i] = new TH1F(TString(Form("%s_effY_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",400,-100,100);
    hist_denY_rot[i] = new TH1F(TString(Form("%s_denY_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",400,-100,100);
    hist_effY_rot[i]->Sumw2();
    hist_denY_rot[i]->Sumw2();
  }

  // Now we get the data
  TChain* chain = new TChain("slim");
  chain->Add(Form("%s/*_slim.root",dir));

  // Get the branches I need:
  vector<double> *xa = 0, *xc = 0,*ya = 0, *yc = 0;
  chain->SetBranchAddress("xA", &xa);
  chain->SetBranchAddress("xC", &xc);
  chain->SetBranchAddress("yA", &ya);
  chain->SetBranchAddress("yC", &yc);

  double pulse[NCH][NTS], ped[NCH];
  chain->SetBranchAddress("pulse", &pulse);
  chain->SetBranchAddress("ped", &ped);

  double intercept_X, slope_X;
  double intercept_Y, slope_Y;
  chain->SetBranchAddress("interceptX", &intercept_X);
  chain->SetBranchAddress("slopeX", &slope_X);
  chain->SetBranchAddress("interceptY", &intercept_Y);
  chain->SetBranchAddress("slopeY", &slope_Y);

  int run;
  chain->SetBranchAddress("run", &run);

  // Event Loop
  int counter = 0;
  for (unsigned int j = 0; j < chain->GetEntries(); ++j) {
    if (debug && j > 1000)
      continue;
    
    // Do the alignment
    chain->GetEntry(j);
    
    double x_hit = intercept_X + z_ex*slope_X;
    double y_hit = intercept_Y + z_ex*slope_Y;
    
    if (j%100000==0)
      std::cout << "  Done with " << j << " / "
                << chain->GetEntriesFast() << " events\n";
    //cout << "(x,y) = (" << x_hit << "," << y_hit << ")\n";
    
    // loop on the various channels
    for (unsigned int i = 0; i < channels.size(); ++i) {
      // For Finger tiles
      if (run < 3410 && i >= 3 && i <= 6)
	continue;
      
      double x_hit_rot = rotate_Point(x_hit, y_hit, i, 'X');
      double y_hit_rot = rotate_Point(x_hit, y_hit, i, 'Y');
      
      if ( isFiducial( i, x_hit, y_hit) ^ isRotFiducial( i, x_hit_rot, y_hit_rot)) {
	cout << "ISSUE WITH" << "\n" 
	     << "i = " << i << "\n"
	     << "x_hit = " << x_hit << "\n"
	     << "y_hit = " << y_hit << "\n"
	     << endl;
	counter++;
      }
      
      double energy_ps = 0;
      for (auto ts : TIMESLICES)
        energy_ps+=pulse[channels[i].chan][ts];
      energy_ps-=TIMESLICES.size()*ped[channels[i].chan];
      
      if (energy_ps>25) {
        hist_eff[i]->Fill(x_hit,y_hit);
        hist_effX[i]->Fill(x_hit);
        hist_effY[i]->Fill(y_hit);
	// Rotated
	hist_eff_rot[i]->Fill(x_hit_rot, y_hit_rot);
        hist_effX_rot[i]->Fill(x_hit_rot);
        hist_effY_rot[i]->Fill(y_hit_rot);
      }
      hist_den[i]->Fill(x_hit,y_hit);
      hist_denX[i]->Fill(x_hit);
      hist_denY[i]->Fill(y_hit);
      // Rotated
      hist_den_rot[i]->Fill(x_hit_rot, y_hit_rot);
      hist_denX_rot[i]->Fill(x_hit_rot);
      hist_denY_rot[i]->Fill(y_hit_rot);
      
    } // loop on channels
    
  } // loop on events
  cout << "counter = " << counter << endl;

  // here I should plot the efficiency maps, after some beautification
  TCanvas *canv[NUMCHAN], *canvX[NUMCHAN], *canvY[NUMCHAN];
  // Plot the rotated hists on these canvases
  TCanvas *canv_rot[NUMCHAN], *canvX_rot[NUMCHAN], *canvY_rot[NUMCHAN];
  
  for (unsigned int i = 0; i < channels.size(); ++i) {
    if (eff) {
      // 2D Efficiency
      // Set up line
      TLine* fid_line = new TLine();
      fid_line->SetLineWidth(2);
      fid_line->SetLineStyle(kDashed);
      
      // Divide to get efficiency
      hist_eff[i]->Divide(hist_eff[i],hist_den[i],1,1,"b");
      // Set Axis
      hist_eff[i]->GetXaxis()->SetTitle("x [mm]");
      hist_eff[i]->GetYaxis()->SetTitle("y [mm]");
      
      canv[i] = new TCanvas(TString(channels[i].name.c_str()).ReplaceAll("-","_").Data(), "", 550, 500);
      canv[i]->SetRightMargin(canv[i]->GetLeftMargin());
      hist_eff[i]->Draw("colz");
      // Draw Lines
      fid_line->DrawLine(fiducialX[i][0],fiducialY[i][0],
			 fiducialX[i][1],fiducialY[i][1]);
      fid_line->DrawLine(fiducialX[i][0],fiducialY[i][0],
			 fiducialX[i][2],fiducialY[i][2]);
      fid_line->DrawLine(fiducialX[i][3],fiducialY[i][3],
			 fiducialX[i][1],fiducialY[i][1]);
      fid_line->DrawLine(fiducialX[i][3],fiducialY[i][3],
			 fiducialX[i][2],fiducialY[i][2]);
      
      TLatex label;
      label.SetNDC();
      label.SetTextSize(0.05);
      label.SetTextAlign(30);
      //label.DrawLatex(0.92,0.875,entry[i].c_str());
      label.DrawLatex(0.8,0.875,entry[i].c_str());
      
      canv[i]->Print(Form("Original_Images/Efficiency_Maps_2D/efficiency_map_%s.png",channels[i].name.c_str()));
      canv[i]->Print(Form("Original_Images/Efficiency_Maps_2D/efficiency_map_%s.C",channels[i].name.c_str()));
    }
    if (eff_rot) {
      // 2D Rotated Efficiency
      TLine* fid_line_rot = new TLine();
      fid_line_rot -> SetLineWidth(2);
      fid_line_rot -> SetLineStyle(kDashed);
      
      hist_eff_rot[i] -> Divide( hist_eff_rot[i], hist_den_rot[i],1,1,"b");
      hist_eff_rot[i] -> GetXaxis() -> SetTitle("x [mm]");
      hist_eff_rot[i] -> GetYaxis() -> SetTitle("y [mm]");
      
      canv_rot[i] = new TCanvas( TString((channels[i].name + "_rot").c_str()).ReplaceAll("-","_").Data(), "", 550, 500);
      canv_rot[i] -> SetRightMargin(canv[i] -> GetLeftMargin());
      hist_eff_rot[i] -> Draw("colz");
      fid_line_rot -> DrawLine( rot_fiducialX[i][0], rot_fiducialY[i][0],
			        rot_fiducialX[i][1], rot_fiducialY[i][1]);
      fid_line_rot -> DrawLine( rot_fiducialX[i][0], rot_fiducialY[i][0],
			        rot_fiducialX[i][2], rot_fiducialY[i][2]);
      fid_line_rot -> DrawLine( rot_fiducialX[i][3], rot_fiducialY[i][3],
			        rot_fiducialX[i][1], rot_fiducialY[i][1]);
      fid_line_rot -> DrawLine( rot_fiducialX[i][3], rot_fiducialY[i][3],
			        rot_fiducialX[i][2], rot_fiducialY[i][2]);
      
      TLatex label_rot;
      label_rot.SetNDC();
      label_rot.SetTextSize(0.05);
      label_rot.SetTextAlign(30);
      //label.DrawLatex(0.92,0.875,entry[i].c_str());
      label_rot.DrawLatex(0.8,0.875, entry[i].c_str());
      
      canv_rot[i] -> Print(Form("Rotated_Images/Efficiency_Maps_2D/efficiency_map_rot%s.png", channels[i].name.c_str()));
      canv_rot[i] -> Print(Form("Rotated_Images/Efficiency_Maps_2D/efficiency_map_rot%s.C", channels[i].name.c_str()));
    }
    
    if (effX) {
      hist_effX[i]->Divide(hist_effX[i],hist_denX[i],1,1,"b");
      hist_effX[i]->GetXaxis()->SetTitle("x [mm]");
      hist_effX[i]->GetYaxis()->SetTitle("Efficiency");
      
      canvX[i] = new TCanvas(TString(Form("effX_%s", channels[i].name.c_str())).ReplaceAll("-","_").Data(),
			     "", 500, 500);
      hist_effX[i]->Draw();
      
      TLatex labelX;
      labelX.SetNDC();
      labelX.SetTextSize(0.05);
      labelX.SetTextAlign(30);
      labelX.DrawLatex(0.92,0.875,entry[i].c_str());
      
      canvX[i]->Print(Form("Original_Images/Efficiency_Maps_X/efficiency_x_%s.png", channels[i].name.c_str()));
      canvX[i]->Print(Form("Original_Images/Efficiency_Maps_X/efficiency_x_%s.C", channels[i].name.c_str()));
    }
    
    if (effX_rot) {
      hist_effX_rot[i]->Divide(hist_effX_rot[i], hist_denX_rot[i], 1, 1, "b");
      hist_effX_rot[i]->GetXaxis()->SetTitle("x [mm]");
      hist_effX_rot[i]->GetYaxis()->SetTitle("Efficiency");
      
      canvX_rot[i] = new TCanvas(TString(Form("effX_%s", (channels[i].name + "_rot").c_str())).ReplaceAll("-","_").Data(),
			     "", 500, 500);
      hist_effX_rot[i]->Draw();
      
      TLatex labelX_rot;
      labelX_rot.SetNDC();
      labelX_rot.SetTextSize(0.05);
      labelX_rot.SetTextAlign(30);
      labelX_rot.DrawLatex(0.92, 0.875, entry[i].c_str());
      
      canvX[i]->Print(Form("Rotated_Images/Efficiency_Maps_X/efficiency_x_rot%s.png", channels[i].name.c_str()));
      canvX[i]->Print(Form("Rotated_Images/Efficiency_Maps_X/efficiency_x_rot%s.C", channels[i].name.c_str()));
    }
    
    if ( effY) {
      hist_effY[i]->Divide(hist_effY[i],hist_denY[i],1,1,"b");
      hist_effY[i]->GetXaxis()->SetTitle("y [mm]");
      hist_effY[i]->GetYaxis()->SetTitle("Efficiency");
      
      canvY[i] = new TCanvas(TString(Form("effY_%s",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
			     "",500,500);
      hist_effY[i]->Draw();
      
      TLatex labelY;
      labelY.SetNDC();
      labelY.SetTextSize(0.05);
      labelY.SetTextAlign(30);
      labelY.DrawLatex(0.92,0.875,entry[i].c_str());
      
      canvY[i]->Print(Form("Original_Images/Efficiency_Maps_Y/efficiency_y_%s.png", channels[i].name.c_str()));
      canvY[i]->Print(Form("Original_Images/Efficiency_Maps_Y/efficiency_y_%s.C", channels[i].name.c_str()));
    }

    if (effY_rot) {
      hist_effY_rot[i]->Divide(hist_effY_rot[i], hist_denY_rot[i], 1, 1, "b");
      hist_effY_rot[i]->GetXaxis()->SetTitle("y [mm]");
      hist_effY_rot[i]->GetYaxis()->SetTitle("Efficiency");
      
      canvY_rot[i] = new TCanvas(TString(Form("effY_%s", (channels[i].name + "_rot").c_str())).ReplaceAll("-","_").Data(),
			     "",500,500);
      hist_effY_rot[i]->Draw();
      
      TLatex labelY_rot;
      labelY_rot.SetNDC();
      labelY_rot.SetTextSize(0.05);
      labelY_rot.SetTextAlign(30);
      labelY_rot.DrawLatex(0.92,0.875, entry[i].c_str());
      
      canvY[i]->Print(Form("Rotated_Images/Efficiency_Maps_Y/efficiency_y_rot%s.png", channels[i].name.c_str()));
      canvY[i]->Print(Form("Rotated_Images/Efficiency_Maps_Y/efficiency_y_rot%s.C", channels[i].name.c_str()));
    }
  }
}

// Energy; time-slice
void doEnergyTS(bool debug, const char* dir) {
  // Make a TFile
  TFile *energy_hists = new TFile("energy_hists.root", "RECREATE");
  if (!energy_hists->IsOpen()) {
    cout << "Something went wrong..." << endl;
    return;
  }

  // Let us define the histogram and book them
  TH1F *hist_en[NUMCHAN]; // 8 is the number of tiles; channels.size() == 8
  TH1F *hist_ts[NUMCHAN];
  TH1F *hist_tsF[NUMCHAN];

  if (channels.size()!= NUMCHAN) {
    cout << "Argh, something wrong!" << endl;
    return;
  }

  for (unsigned int i = 0; i < channels.size(); ++i) {
    hist_en[i] = new TH1F(TString(Form("en_%s",channels[i].name.c_str())).ReplaceAll("-","_").Data(),"",247,edges);
    hist_ts[i] = new TH1F(TString(Form("ts_%s",channels[i].name.c_str())).ReplaceAll("-","_").Data(),"",
                          10,0.5,10.5);
    hist_tsF[i] = new TH1F(TString(Form("ts_%s",channels[i].name.c_str())).ReplaceAll("-","_").Data(),"",
                          10,0.5,10.5);
  }

  // Now we get the data
  TChain* chain = new TChain("slim");
  chain->Add(Form("%s/*_slim.root",dir));

  // Get the branches I need:
  vector<double> *xa = 0, *xc = 0, *ya = 0, *yc = 0;
  chain->SetBranchAddress("xA",&xa);
  chain->SetBranchAddress("xC",&xc);
  chain->SetBranchAddress("yA",&ya);
  chain->SetBranchAddress("yC",&yc);

  double pulse[NCH][NTS], ped[NCH];
  chain->SetBranchAddress("pulse",&pulse);
  chain->SetBranchAddress("ped",&ped);

  double intercept_X, slope_X;
  double intercept_Y, slope_Y;
  chain->SetBranchAddress("interceptX",&intercept_X);
  chain->SetBranchAddress("slopeX",&slope_X);
  chain->SetBranchAddress("interceptY",&intercept_Y);
  chain->SetBranchAddress("slopeY",&slope_Y);

  int run;
  chain->SetBranchAddress("run", &run);

  for (unsigned int i=0;i<chain->GetEntries();++i) {
    if (debug && i>1000)
      continue;


    // Do the alignment
    chain->GetEntry(i);

    if (i%100000==0)
      std::cout << "  Done with " << i << " / "
                << chain->GetEntriesFast() << " events\n";

    double x_hit = intercept_X + z_ex*slope_X;
    double y_hit = intercept_Y + z_ex*slope_Y;

    // loop on the various channels
    for (unsigned int i = 0; i < channels.size(); ++i) {
      // This ensures that the finger tiles are using the correct runs
      if (run < 3410 && i >= 3 && i <= 6)
        continue;

      // Use the time slices numbers contained in TIMESLICES
      // and remove pedestal times # of time slices used 
      double energy_ps = 0;
      for (auto ts : TIMESLICES)
        energy_ps+=pulse[channels[i].chan][ts];
      energy_ps-=TIMESLICES.size()*ped[channels[i].chan];

      if(!isFiducial(i,x_hit,y_hit)) // THIS IS THE MOST IMPORTANT STEP 
        continue;

      hist_en[i]->Fill(energy_ps);

      for (int t=0;t<NTS;++t) {
        if (i < 3 || i == 7) {
          hist_ts[i]->Fill(t+1,
                           pulse[channels[i].chan][t]-
                           ped[channels[i].chan]);
        }
        else {
          hist_tsF[i]->Fill(t+1,
                            pulse[channels[i].chan][t]-
                            ped[channels[i].chan]);
        }
      }
    } // loop on channels
  } // loop on events

  // Make energy plots
  TCanvas* canv[NUMCHAN];

  for (unsigned int i = 0; i < channels.size(); ++i) {

    hist_en[i]->GetXaxis()->SetTitle("Charge [fC]");
    hist_en[i]->GetYaxis()->SetTitle("Events");
    hist_en[i]->SetLineWidth(2);
    hist_en[i]->SetLineColor(color[i]);

    canv[i] = new TCanvas(TString(channels[i].name.c_str()).ReplaceAll("-","_").Data(), "", 500, 500);
    canv[i]->SetLogx();
    canv[i]->SetLogy();
    hist_en[i]->Draw("colz");

    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.05);
    label.SetTextAlign(30);
    label.DrawLatex(0.92,0.875,entry[i].c_str());

    label.SetTextAlign(11);
    
    float eff = hist_en[i]->Integral(hist_en[i]->FindBin(25),
                                     hist_en[i]->GetNbinsX()) / hist_en[i]->GetEntries();
    
    //float eff_err = TMath::Sqrt(eff*(1-eff)/hist_en[i]->GetEntries());
    eff*=100;
    //eff_err*=100;
    label.DrawLatex(0.20,0.225,Form("#splitline{#epsilon=%4.1f%%}"
                                    "{Mean=%3.1f#pm%3.1ffC}",
                                    eff,
                                    hist_en[i]->GetMean(),
                                    hist_en[i]->GetMeanError()));
    
    canv[i]->Print(Form("Energy_Plots/energy_PS_%s.png",channels[i].name.c_str()));
  }
  

  TCanvas* canv_ts = new TCanvas("ts","",500,500);
  
  TLegend* leg = new TLegend(0.2,0.7,0.4,0.9,"","brNDC");
  leg->SetTextSize(0.05);
  TLegend* leg_fing = new TLegend(0.2,0.7,0.4,0.9,"","brNDC");
  leg_fing->SetTextSize(0.05);

  // Let us also get the histogram with max value...
  int max_index = 0;
  int max_indexF = 0;
  float max_value = -99;
  float max_valueF = -99;
  
  for (unsigned int i = 0; i < channels.size(); ++i) {
    if ((i < 3) || (i == 7)) {
      hist_ts[i]->SetLineWidth(2);
      hist_ts[i]->SetLineColor(color[i]);
      hist_ts[i]->SetLineStyle(style[i]);

      hist_ts[i]->GetXaxis()->SetTitle("Time Slice [25ns]");
      hist_ts[i]->GetXaxis()->SetNdivisions(NTS);
      hist_ts[i]->GetYaxis()->SetTitle("Charge [fC]");
      // changed this from 1.6
      hist_ts[i]->GetYaxis()->SetTitleOffset(3.0);
      leg->AddEntry(hist_ts[i],entry[i].c_str(),"l");

      if (hist_ts[i]->GetMaximum()>max_value) {
        max_index = i;
        max_value = hist_ts[i]->GetMaximum();
      }
    }
    else { // For finger tiles                                                                                                                  
      hist_tsF[i]->SetLineWidth(2);
      hist_tsF[i]->SetLineColor(color[i]);
      hist_tsF[i]->SetLineStyle(style[i]);

      hist_tsF[i]->GetXaxis()->SetTitle("Time Slice [25ns]");
      hist_tsF[i]->GetXaxis()->SetNdivisions(NTS);
      hist_tsF[i]->GetYaxis()->SetTitle("Charge [fC]");
      // changed this from 1.6
      hist_tsF[i]->GetYaxis()->SetTitleOffset(3.0);
      leg_fing->AddEntry(hist_tsF[i],entry[i].c_str(),"l");

      if (hist_tsF[i]->GetMaximum()>max_valueF) {
        max_indexF = i;
        max_valueF = hist_tsF[i]->GetMaximum();
      }
    }
  }
  hist_ts[max_index]->Draw("hist");
  for (unsigned int i=0;i<channels.size(); i++) {
    if ((i < 3) || (i == 7))
      hist_ts[i]->Draw("hist,same");
  }
  leg->Draw("same");
  hist_ts[max_index]->Draw("axis,same");

  canv_ts->Print("Time_Slice_Plots/ts.png");
  canv_ts->Print("Time_Slice_Plots/ts.C");


  TCanvas* canv_tsF = new TCanvas("ts","",500,500);
  hist_tsF[max_indexF]->Draw("hist");
  for (unsigned int i=0;i<channels.size(); i++) {
    if (i > 2 || i < 7)
      hist_tsF[i]->Draw("hist,same");
  }

  leg_fing->Draw("same");
  hist_tsF[max_indexF]->Draw("axis,same");

  canv_tsF->Print("Time_Slice_Plots/tsF.png");
  canv_tsF->Print("Time_Slice_Plots/tsF.C");
  energy_hists->Close();
}

// THIS IS NOT NEEDED
/*
// delta-t w.r.t. SCSN-81
void doTime(const char*
            dir="/data/users/abelloni/CERN_TB_Aug15_slim_ntuples",
            bool debug=false) {

  // Let us define the histogram and book them
  TH1F *hist_dt[NUMCHAN]; // 8 is the number of tiles; channels.size() == 8
  if (channels.size()!= NUMCHAN) {
    cout << "Argh, something wrong!" << endl;
    return;
  }
  for (unsigned int i = 0; i < channels.size(); ++i) {
    hist_dt[i] = new TH1F(Form("dt_%s",channels[i].name.c_str()),"",
                          320,-160,160);
  }

  // Now we get the data
  TChain* chain = new TChain("slim");
*/
  //chain->Add(Form("%s/*_slim.root",dir));
/*
  // Get the branches I need:
  vector<double> *xa = 0,*xb = 0,*xc = 0,*ya = 0,*yb = 0,*yc = 0;
  chain->SetBranchAddress("xA",&xa);
  //chain->SetBranchAddress("xB",&xb);
  chain->SetBranchAddress("xC",&xc);
  chain->SetBranchAddress("yA",&ya);
  //chain->SetBranchAddress("yB",&yb);
  chain->SetBranchAddress("yC",&yc);

  double pulse_tdc[29][10];
  chain->SetBranchAddress("pulse_tdc",&pulse_tdc);

  double pulse[29][10],ped[29];
  chain->SetBranchAddress("pulse",&pulse);
  chain->SetBranchAddress("ped",&ped);

  double intercept_X, slope_X;
  double intercept_Y, slope_Y;
  chain->SetBranchAddress("interceptX",&intercept_X);
  chain->SetBranchAddress("slopeX",&slope_X);
  chain->SetBranchAddress("interceptY",&intercept_Y);
  chain->SetBranchAddress("slopeY",&slope_Y);

  int run;
  chain->SetBranchAddress("run", &run);

  for (unsigned int i=0;i<chain->GetEntries();++i) {
    if (debug && i>1000)
      continue;

    // Do the alignment
    chain->GetEntry(i);

    if (i%100000==0)
      std::cout << "  Done with " << i << " / "
                << chain->GetEntriesFast() << " events\n";

    double x_hit = intercept_X + z_ex*slope_X;
    double y_hit = intercept_Y + z_ex*slope_Y;

    // time_ref is the SCSN-81 time, used in all delta(t)
    double time_ref = -99;

    // loop on the various channels
    for (unsigned int i = 0; i < channels.size(); ++i) {

      if(!isFiducial(i,x_hit,y_hit))
        continue;

      // Use TS = 5, 6, 7, 8; remove 4 times the pedestal
      double energy_ps =
	pulse[channels[i].chan][5]+
	pulse[channels[i].chan][6]+
	pulse[channels[i].chan][7]+
	pulse[channels[i].chan][8]-
	4*ped[channels[i].chan];
      if (energy_ps<=25) // require >25fC hits in both scintillators
	continue;
      
      double time = -99;
      // Get best time for current channel
      for (int ts=0;ts<10;++ts)
        if (pulse_tdc[channels[i].chan][ts]>-99 &&
            fabs(pulse_tdc[channels[i].chan][ts])>0.01) {
          time = pulse_tdc[channels[i].chan][ts];
          if(i==0)
            time_ref = time;
        }

      // Now I have the good time, can plot difference
      // Plot will naturally be 0 in SCSN case)
      if (time>-99 && time_ref>-99)
        hist_dt[i]->Fill(time-time_ref);

    } // loop on channels

  } // loop on events

  // here I make the plot, after some beautification
  TCanvas* canv[NUMCHAN];

  for (unsigned int i = 0; i < channels.size(); ++i) {

    hist_dt[i]->GetXaxis()->
      SetTitle(Form("t_{%s}-t_{SCSN-81} [ns]",entry[i].c_str()));
    hist_dt[i]->GetYaxis()->SetTitle("Events");
    hist_dt[i]->SetLineWidth(2);
    hist_dt[i]->SetLineColor(color[i]);

    hist_dt[i]->SetMaximum(hist_dt[i]->GetMaximum()*
                           TMath::Power(10,
                                        TMath::Log(hist_dt[i]->GetMaximum())/
                                        TMath::Log(10)/10.*2));

    canv[i] = new TCanvas(channels[i].name.c_str(),"",500,500);
    canv[i]->SetLogy();
    hist_dt[i]->Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.05);
    label.DrawLatex(0.20,0.875,Form("<#deltat>=%3.1f#pm%3.1fns",
                                    hist_dt[i]->GetMean(),
                                    hist_dt[i]->GetMeanError()));

    canv[i]->Print(Form("Time_Slice_Plots/time_%s.png",channels[i].name.c_str()));

  }


}

bool isFiducial(int i, float x_hit, float y_hit) {

  // Find if it is fiducial!
  int polyCorners = 4;
  int k, j=polyCorners-1 ;
  bool oddNodes=false;
  for (k=0; k<polyCorners; k++) {
    if (((fiducialY[i][k]< y_hit && fiducialY[i][j]>=y_hit) ||
         (fiducialY[i][j]< y_hit && fiducialY[i][k]>=y_hit)) &&
        (fiducialX[i][k]<=x_hit || fiducialX[i][j]<=x_hit)) {
      oddNodes^=(fiducialX[i][k]+(y_hit-fiducialY[i][k])/
                 (fiducialY[i][j]-fiducialY[i][k])*
                 (fiducialX[i][j]-fiducialX[i][k])<x_hit);
    }
    j=k;
  }

  return oddNodes;

}

// Corrections to get alignment:
// # Wire chamber means and standard deviations (xA-xC, xA-xC, yA-yC, etc.)
// wc_res = {}
// wc_res["x", "BC", "mean"] = -5.83e-01
// wc_res["y", "BC", "mean"] = -1.75e+01
// wc_res["x", "AC", "mean"] = -1.24e+00
// wc_res["y", "AC", "mean"] = -8.78e+00
// wc_res["x", "BC", "rms" ] =  3.96e+00
// wc_res["y", "BC", "rms" ] =  3.88e+00
// wc_res["x", "AC", "rms" ] =  4.30e+00
// wc_res["y", "AC", "rms" ] =  5.08e+00

// Based on numbers above, I get:
//     "xA-xB+0.657",
//     "xB-xC+0.583",
//     "xC-xA-1.24",
//     "yA-yB-1.12",
//     "yB-yC-3.96",
//     "yC-yA+5.08"

*/
