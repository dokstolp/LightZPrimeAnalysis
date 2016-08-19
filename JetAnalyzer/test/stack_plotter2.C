#include <fstream>
#include <vector>
#include <iomanip>
#include "TFile.h"
#include "TH2.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TLegend.h"

void stack_plotter2()
{
  TCanvas *c = new TCanvas("c", "canvas",700,640);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  c->SetLeftMargin(0.15);
  c->cd();
  
  //opening the data file and adding "METValue_4" histogram

  TFile *f_datafile_0 = new TFile("postMETdata_0.root");
  TFile *f_datafile_1 = new TFile("postMETdata_1.root");
  TH1F *histo_j1EtaWidth_data_0 = (TH1F*)f_datafile_0->Get("METValue_4");
  TH1F *histo_j1EtaWidth_data_1 = (TH1F*)f_datafile_1->Get("METValue_4");
  histo_j1EtaWidth_data_0->Add(histo_j1EtaWidth_data_1);
  histo_j1EtaWidth_data_0->SetStats(0);
  histo_j1EtaWidth_data_0->SetLineWidth(2);
  histo_j1EtaWidth_data_0->SetLineColor(kWhite);
  histo_j1EtaWidth_data_0->SetTitle("");
  histo_j1EtaWidth_data_0->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_data_0->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_data_0->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_data_0->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_data_0->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_data_0->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_data_0->Draw("");

  //opening background WJets Sample file
  TFile *f_WJets_0 = new TFile("postWJets_MLM_0.root");
  TFile *f_WJets_1 = new TFile("postWJets_MLM_1.root");
  TFile *f_W1Jets = new TFile("postW100to200_0.root");
  TFile *f_W1Jets_1 = new TFile("postW100to200_1.root");
  TFile *f_W1Jets_2 = new TFile("postW100to200_2.root");
  TFile *f_W2Jets = new TFile("postW200to400_0.root");
  TFile *f_W2Jets_1 = new TFile("postW200to400_1.root");
  TFile *f_W3Jets = new TFile("postW400to600_0.root");
  TFile *f_W4Jets = new TFile("postW600to800_0.root");
  TFile *f_W4Jets_1 = new TFile("postW600to800_1.root");
  TFile *f_W5Jets = new TFile("postW800to1200_0.root");
  TFile *f_W6Jets = new TFile("postW1200to2500_backup.root");
  TFile *f_W7Jets = new TFile("postW2500toInf_0.root");

  TH1F *histo_j1EtaWidth_WJets_0 = (TH1F*)f_WJets_0->Get("METValue_4");
  TH1F *histo_j1EtaWidth_WJets_1 = (TH1F*)f_WJets_1->Get("METValue_4");
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_WJets_1);
  TH1F *histo_j1EtaWidth_W1Jets = (TH1F*)f_W1Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_W1Jets_1 = (TH1F*)f_W1Jets_1->Get("METValue_4");
  TH1F *histo_j1EtaWidth_W1Jets_2 = (TH1F*)f_W1Jets_2->Get("METValue_4");
  histo_j1EtaWidth_W1Jets->Add(histo_j1EtaWidth_W1Jets_1);
  histo_j1EtaWidth_W1Jets->Add(histo_j1EtaWidth_W1Jets_2);
  TH1F *histo_j1EtaWidth_W2Jets = (TH1F*)f_W2Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_W2Jets_1 = (TH1F*)f_W2Jets_1->Get("METValue_4");
  histo_j1EtaWidth_W2Jets->Add(histo_j1EtaWidth_W2Jets_1);
  TH1F *histo_j1EtaWidth_W3Jets = (TH1F*)f_W3Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_W4Jets = (TH1F*)f_W4Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_W4Jets_1 = (TH1F*)f_W4Jets_1->Get("METValue_4");
  histo_j1EtaWidth_W4Jets->Add(histo_j1EtaWidth_W4Jets_1);
  TH1F *histo_j1EtaWidth_W5Jets = (TH1F*)f_W5Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_W6Jets = (TH1F*)f_W6Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_W7Jets = (TH1F*)f_W7Jets->Get("METValue_4");
  
  histo_j1EtaWidth_WJets_0->SetStats(0);
  histo_j1EtaWidth_W1Jets->SetStats(0);
  histo_j1EtaWidth_W2Jets->SetStats(0);
  histo_j1EtaWidth_W3Jets->SetStats(0);
  histo_j1EtaWidth_W4Jets->SetStats(0); 
  histo_j1EtaWidth_W5Jets->SetStats(0);
  histo_j1EtaWidth_W6Jets->SetStats(0);
  histo_j1EtaWidth_W7Jets->SetStats(0);

   // Scaling = (1/Totalevents)*Luminosity*NNLO-cross-section
  histo_j1EtaWidth_WJets_0->Scale((1.0/27320430)*2579.525*50690*1.21);
  histo_j1EtaWidth_W1Jets->Scale((1.0/28999950)*2579.525*1345*1.21);
  histo_j1EtaWidth_W2Jets->Scale((1.0/15004280)*2579.525*359.7*1.21);
  histo_j1EtaWidth_W3Jets->Scale((1.0/5445331)*2579.525*48.91*1.21);
  histo_j1EtaWidth_W4Jets->Scale((1.0/14490730)*2579.525*12.05*1.21);
  histo_j1EtaWidth_W5Jets->Scale((1.0/6314257)*2579.525*5.501*1.21);
  histo_j1EtaWidth_W6Jets->Scale((1.0/6857102)*2579.525*1.329*1.21);
  histo_j1EtaWidth_W7Jets->Scale((1.0/2254248)*2579.525*0.03216*1.21);
  
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W1Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W2Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W3Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W4Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W5Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W6Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W7Jets);  

  histo_j1EtaWidth_WJets_0->SetTitle("");
  histo_j1EtaWidth_WJets_0->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_WJets_0->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_WJets_0->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WJets_0->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_WJets_0->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_WJets_0->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WJets_0->SetFillColor(kRed);

  TFile *f_Zvv_100to200 = new TFile("postZ100to200_final.root");
  TFile *f_Zvv_200to400 = new TFile("postZ200to400_final.root");
  TFile *f_Zvv_400to600 = new TFile("postZ400to600_0.root");
  TFile *f_Zvv_600to800 = new TFile("postZ600to800_0.root");
  TFile *f_Zvv_800to1200 = new TFile("postZ800to1200_0.root");
  TFile *f_Zvv_1200to2500 = new TFile("postZ1200to2500_0.root");
  TFile *f_Zvv_2500toInf = new TFile("postZ2500toInf_0.root");

  TH1F *histo_j1EtaWidth_100to200 = (TH1F*)f_Zvv_100to200->Get("METValue_4");
  TH1F *histo_j1EtaWidth_200to400 = (TH1F*)f_Zvv_200to400->Get("METValue_4");
  TH1F *histo_j1EtaWidth_400to600 = (TH1F*)f_Zvv_400to600->Get("METValue_4");
  TH1F *histo_j1EtaWidth_600to800 = (TH1F*)f_Zvv_600to800->Get("METValue_4");
  TH1F *histo_j1EtaWidth_800to1200 = (TH1F*)f_Zvv_800to1200->Get("METValue_4");
  TH1F *histo_j1EtaWidth_1200to2500 = (TH1F*)f_Zvv_1200to2500->Get("METValue_4");
  TH1F *histo_j1EtaWidth_2500toInf = (TH1F*)f_Zvv_2500toInf->Get("METValue_4");
  histo_j1EtaWidth_100to200->SetStats(0);
  histo_j1EtaWidth_200to400->SetStats(0);
  histo_j1EtaWidth_400to600->SetStats(0);
  histo_j1EtaWidth_600to800->SetStats(0);
  histo_j1EtaWidth_800to1200->SetStats(0);
  histo_j1EtaWidth_1200to2500->SetStats(0);
  histo_j1EtaWidth_2500toInf->SetStats(0);

  // Scaling = (1/Totalevents)*Luminosity*LO-cross-section
  histo_j1EtaWidth_100to200->Scale((1.0/18887290)*2579.525*280.35*1.23);
  histo_j1EtaWidth_200to400->Scale((1.0/19881040)*2579.525*77.67*1.23);
  histo_j1EtaWidth_400to600->Scale((1.0/8692385)*2579.525*10.73*1.23);
  histo_j1EtaWidth_600to800->Scale((1.0/5712221)*2579.525*2.5611*1.23);
  histo_j1EtaWidth_800to1200->Scale((1.0/2141386)*2579.525*1.1826*1.23);
  histo_j1EtaWidth_1200to2500->Scale((1.0/369514)*2579.525*0.2874*1.23);
  histo_j1EtaWidth_2500toInf->Scale((1.0/405752)*2579.525*0.006924*1.23);

  //Add the ZJetsToNuNu histograms to the first one
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_200to400);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_400to600);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_600to800);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_800to1200);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_1200to2500);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_2500toInf);
  histo_j1EtaWidth_100to200->SetTitle("");
  histo_j1EtaWidth_100to200->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_100to200->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_100to200->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_100to200->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_100to200->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_100to200->GetYaxis()->SetLabelOffset(999);  
  histo_j1EtaWidth_100to200->SetFillColor(kBlue);

  //opening background samples Gamma+jets
  TFile *f_G1Jets = new TFile("postGJets40to100.root");
  TFile *f_G2Jets = new TFile("postGJets100to200.root");
  TFile *f_G3Jets = new TFile("postGJets200to400.root");
  TFile *f_G4Jets = new TFile("postGJets400to600.root ");
  TFile *f_G5Jets = new TFile("postGJets600toInf.root");
  
  TH1F *histo_j1EtaWidth_G1Jets = (TH1F*)f_G1Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_G2Jets = (TH1F*)f_G2Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_G3Jets = (TH1F*)f_G3Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_G4Jets = (TH1F*)f_G4Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_G5Jets = (TH1F*)f_G5Jets->Get("METValue_4");

  histo_j1EtaWidth_G1Jets->SetStats(0);
  histo_j1EtaWidth_G2Jets->SetStats(0);
  histo_j1EtaWidth_G3Jets->SetStats(0);
  histo_j1EtaWidth_G4Jets->SetStats(0);
  histo_j1EtaWidth_G5Jets->SetStats(0);

  //Scaling
  histo_j1EtaWidth_G1Jets->Scale((1.0/4468724.0)*2579.525*20790);
  histo_j1EtaWidth_G2Jets->Scale((1.0/5142780.0)*2579.525*9238);
  histo_j1EtaWidth_G3Jets->Scale((1.0/10322700.0)*2579.525*2305);
  histo_j1EtaWidth_G4Jets->Scale((1.0/2512170.0)*2579.525*274.4);
  histo_j1EtaWidth_G5Jets->Scale((1.0/2459260.0)*2579.525*93.46);
  
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G2Jets);
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G3Jets);
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G4Jets);
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G5Jets);

  histo_j1EtaWidth_G1Jets->SetTitle("");
  histo_j1EtaWidth_G1Jets->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_G1Jets->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_G1Jets->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_G1Jets->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_G1Jets->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_G1Jets->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_G1Jets->SetFillColor(kGreen);

  //opening background TTJets
  TFile *f_TTJets = new TFile("postTTJets_MLM.root");
  TH1F *histo_j1EtaWidth_TTJets = (TH1F*)f_TTJets->Get("METValue_4");
  histo_j1EtaWidth_TTJets->SetStats(0);
  histo_j1EtaWidth_TTJets->Scale((1.0/10273350)*2579.525*502.2);
  histo_j1EtaWidth_TTJets->SetTitle("");
  histo_j1EtaWidth_TTJets->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_TTJets->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_TTJets->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_TTJets->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_TTJets->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_TTJets->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_TTJets->SetFillColor(kOrange);

 //addding some backgrounds like ZllG, WW, WZ, ZZ	
  TFile *f_ZllG = new TFile("postZllG.root");
  TH1F *histo_j1EtaWidth_ZllG = (TH1F*)f_ZllG->Get("METValue_4");
  histo_j1EtaWidth_ZllG->SetStats(0);
  histo_j1EtaWidth_ZllG->Scale((1.0/489448)*2579.525*0.143 );
  histo_j1EtaWidth_ZllG->SetTitle("");
  histo_j1EtaWidth_ZllG->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_ZllG->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_ZllG->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_ZllG->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_ZllG->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_ZllG->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_ZllG->SetFillColor(kSpring);

  TFile *f_WW = new TFile("postWW.root");
  TH1F *histo_j1EtaWidth_WW = (TH1F*)f_WW->Get("METValue_4");
  histo_j1EtaWidth_WW->SetStats(0);
  histo_j1EtaWidth_WW->Scale((1.0/993214)*2579.525*118.7);
  histo_j1EtaWidth_WW->SetTitle("");
  histo_j1EtaWidth_WW->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_WW->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_WW->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WW->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_WW->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_WW->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WW->SetFillColor(kAzure+10);
 
  TFile *f_WZ = new TFile("postWZ.root");
  TH1F *histo_j1EtaWidth_WZ = (TH1F*)f_WZ->Get("METValue_4");
  histo_j1EtaWidth_WZ->SetStats(0);
  histo_j1EtaWidth_WZ->Scale((1.0/1000000)*2579.525*47.2);
  histo_j1EtaWidth_WZ->SetTitle("");
  histo_j1EtaWidth_WZ->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_WZ->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_WZ->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WZ->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_WZ->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_WZ->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WZ->SetFillColor(kCyan);

  TFile *f_ZZ = new TFile("postZZ.root");
  TH1F *histo_j1EtaWidth_ZZ = (TH1F*)f_ZZ->Get("METValue_4");
  histo_j1EtaWidth_ZZ->SetStats(0);
  histo_j1EtaWidth_ZZ->Scale((1.0/989312)*2579.525*16.6);
  histo_j1EtaWidth_ZZ->SetTitle("");
  histo_j1EtaWidth_ZZ->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_ZZ->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_ZZ->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_ZZ->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_ZZ->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_ZZ->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_ZZ->SetFillColor(kCyan+2);

 //opening QCD background files (HT-binned samples)
  TFile *f_Q1Jets = new TFile("postQCD100to200_final.root");
  TFile *f_Q2Jets = new TFile("postQCD200to300_final.root");
  TFile *f_Q3Jets = new TFile("postQCD300to500_final.root");
  TFile *f_Q4Jets = new TFile("postQCD500to700_final.root");
  TFile *f_Q5Jets = new TFile("postQCD700to1000_final.root");
  TFile *f_Q6Jets = new TFile("postQCD1000to1500_final.root");
  TFile *f_Q7Jets = new TFile("postQCD1500to2000_0.root");
  TFile *f_Q8Jets = new TFile("postQCD2000toInf_0.root");

  TH1F *histo_j1EtaWidth_Q1Jets = (TH1F*)f_Q1Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_Q2Jets = (TH1F*)f_Q2Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_Q3Jets = (TH1F*)f_Q3Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_Q4Jets = (TH1F*)f_Q4Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_Q5Jets = (TH1F*)f_Q5Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_Q6Jets = (TH1F*)f_Q6Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_Q7Jets = (TH1F*)f_Q7Jets->Get("METValue_4");
  TH1F *histo_j1EtaWidth_Q8Jets = (TH1F*)f_Q8Jets->Get("METValue_4");

  histo_j1EtaWidth_Q1Jets->SetStats(0);
  histo_j1EtaWidth_Q2Jets->SetStats(0);
  histo_j1EtaWidth_Q3Jets->SetStats(0);
  histo_j1EtaWidth_Q4Jets->SetStats(0);
  histo_j1EtaWidth_Q5Jets->SetStats(0);
  histo_j1EtaWidth_Q6Jets->SetStats(0);
  histo_j1EtaWidth_Q7Jets->SetStats(0);
  histo_j1EtaWidth_Q8Jets->SetStats(0);

  histo_j1EtaWidth_Q1Jets->Scale((1.0/82185100)*2579.525*27850000);
  histo_j1EtaWidth_Q2Jets->Scale((1.0/38705000)*2579.525*1717000);
  histo_j1EtaWidth_Q3Jets->Scale((1.0/37696900)*2579.525*351300);
  histo_j1EtaWidth_Q4Jets->Scale((1.0/44138700)*2579.525*31630);
  histo_j1EtaWidth_Q5Jets->Scale((1.0/29832300)*2579.525*6802);
  histo_j1EtaWidth_Q6Jets->Scale((1.0/10336000)*2579.525*1206);
  histo_j1EtaWidth_Q7Jets->Scale((1.0/7859230)*2579.525*120.4);
  histo_j1EtaWidth_Q8Jets->Scale((1.0/4047530)*2579.525*25.25);  

  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q2Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q3Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q4Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q5Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q6Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q7Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q8Jets);
  histo_j1EtaWidth_Q1Jets->SetTitle("");
  histo_j1EtaWidth_Q1Jets->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_Q1Jets->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_Q1Jets->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_Q1Jets->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_Q1Jets->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_Q1Jets->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_Q1Jets->SetFillColor(kMagenta);

  //Stack histograms using THStack
  THStack *hs_datamc = new THStack("hs_datamc","Data/MC comparison"); 
  hs_datamc->Add(histo_j1EtaWidth_ZllG);
  hs_datamc->Add(histo_j1EtaWidth_WW);
  hs_datamc->Add(histo_j1EtaWidth_WZ);
  hs_datamc->Add(histo_j1EtaWidth_ZZ);
  hs_datamc->Add(histo_j1EtaWidth_G1Jets);
  hs_datamc->Add(histo_j1EtaWidth_TTJets);
  hs_datamc->Add(histo_j1EtaWidth_100to200);
  hs_datamc->Add(histo_j1EtaWidth_WJets_0);
  hs_datamc->Add(histo_j1EtaWidth_Q1Jets);
  hs_datamc->SetTitle("");
  hs_datamc->Draw("HIST SAMEP0E1");
  
  histo_j1EtaWidth_data_0->SetLineColor(kBlack);
  histo_j1EtaWidth_data_0->Draw("SAMEP0E1");
  
  //TLegend *leg = new TLegend(0.181948,0.663948,0.567335,0.836868,"");
  TLegend *leg = new TLegend(0.55,0.35,0.65,0.65,"");
  leg->AddEntry(histo_j1EtaWidth_data_0,"MET2016B Data");
  leg->AddEntry(histo_j1EtaWidth_WJets_0,"W#rightarrowl#nu MC","F");
  leg->AddEntry(histo_j1EtaWidth_100to200,"Z#rightarrow#nu#nu MC","F");
  leg->AddEntry(histo_j1EtaWidth_Q1Jets, "QCD MC","F");
  leg->AddEntry(histo_j1EtaWidth_TTJets, "Top Quark MC", "F");
  leg->AddEntry(histo_j1EtaWidth_G1Jets,"#gamma+jets MC", "F");
  leg->AddEntry(histo_j1EtaWidth_ZllG, "Z#rightarrowll#gamma","F");
  leg->AddEntry(histo_j1EtaWidth_WW,"WW","F");
  leg->AddEntry(histo_j1EtaWidth_WZ,"WZ","F");
  leg->AddEntry(histo_j1EtaWidth_ZZ,"ZZ","F");
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->Draw();  

  TLatex *texS = new TLatex(0.45,0.807173,"#sqrt{s} = 13 TeV, 2.58 fb^{-1}");
  texS->SetNDC();
  texS->SetTextFont(42);
  texS->SetTextSize(0.040);
  texS->Draw();
  TLatex *texS1 = new TLatex(0.16092,0.907173,"#bf{UsamaHussain:} #it{Analysis in Progress}");
  texS1->SetNDC();
  texS1->SetTextFont(42);
  texS1->SetTextSize(0.040);
  texS1->Draw();

  c->Update();
 // hs_datamc->GetXaxis()->SetTitle("pfMET [GeV]");
 // hs_datamc->GetYaxis()->SetTitle("Events");
 // hs_datamc->GetYaxis()->SetLabelOffset(999);
 
  double xmin = c->GetUxmin();
  double ymin = c->GetUymin();
  double xmax = c->GetUxmax();
  double ymax = c->GetUymax();

  TGaxis *xaxis = new TGaxis(xmin,ymin,xmax,ymin,xmin,xmax,510);
  xaxis->SetTitle("pfMET [GeV]");
  xaxis->SetLabelFont(42);
  xaxis->SetLabelSize(0.030);
  xaxis->SetTitleFont(42);
  xaxis->SetTitleSize(0.035);
  xaxis->Draw("SAME");

  TGaxis *xaxis_top = new TGaxis(xmin,ymax,xmax,ymax,xmin,xmax,510,"-");
  xaxis_top->SetTitle("");
  xaxis_top->SetLabelOffset(999);
  xaxis_top->Draw("SAME");

  TGaxis *yaxis = new TGaxis(xmin,ymin,xmin,ymax,ymin,ymax,510);
  yaxis->SetTitle("Events");
  yaxis->SetLabelFont(42);
  yaxis->SetLabelSize(0.030);
  yaxis->SetTitleFont(42);
  yaxis->SetTitleSize(0.035);
  yaxis->SetTitleOffset(1.8);
  yaxis->Draw("SAME");

  TGaxis *yaxis_right = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+");
  yaxis_right->SetTitle("");
  yaxis_right->SetLabelOffset(999);
  yaxis_right->Draw("SAME");  
 
 
  c->SaveAs("datamc_METValue_4.png");

}
