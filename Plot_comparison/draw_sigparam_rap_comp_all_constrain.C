#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <cstring>
#include <math.h>

#include <TROOT.h>
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TMath.h>
#include <math.h>
#include <TH1.h>
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TClonesArray.h"
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>

#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom.h>

#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooCategory.h>
#include <RooGenericPdf.h>
#include <RooFFTConvPdf.h>
#include <RooWorkspace.h>
#include <RooBinning.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooConstVar.h>
#include "../SONGKYO.h"
#include "../commonUtility.h"
#include "../tdrstyle.C"
#include "../CMS_lumi_internal.C"

using namespace std;
using namespace RooFit;

double *getAvg(double avg_alpha = 0, double avg_n =0, double avg_sigma=0, double avg_f=0, double avg_x=0, int states=1, TString szAA = "PP");
int draw_sigparam_rap_comp_all_constrain(TString szAA = "AA", int states =1, int DrawOpt = 0, int nFit=3) 
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod;
  if(szAA=="PP") iPeriod=1; 
  else if(szAA=="AA") iPeriod=2; 
  int iPos = 33;
  /////////////////////////////////////////////////////////
  //// set style
  /////////////////////////////////////////////////////////
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);

/*  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);
*/
  gStyle->SetTextSize(0.04);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.048,"xyz");
  gStyle->SetTitleOffset(1.6,"y"); // KYO for yield

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.12) ; 
  gStyle->SetPadLeftMargin(0.16) ; // KYO for yield

  gStyle->SetEndErrorSize(0);  
  /////////////////////////////////////////////////////////
  //// binning setting
  /////////////////////////////////////////////////////////
  double tmpArr1s[7] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  double tmpArr2s[4] = {0.0, 0.8, 1.6, 2.4};
  double tmpArr3s[3] = {0.0, 1.2, 2.4};

  int tmpBin;
  if ( states ==1) {
    cout << " ***** 1S *****" << endl; tmpBin = 6;
  }else if (states ==2){
    cout << " ***** 2S *****" << endl; tmpBin = 3;
  }else if (states ==3){
    cout << " ***** 3S *****" << endl; tmpBin = 2;
  }else {
    cout << " Error ::: Select among 1S, 2S, and 3S" << endl; return 0;
  }
 
  const int nStates = 3; 
  const int nBin = tmpBin; // number of bin 
  const int nArrNum = nBin+1; // number of array
  double binArr[nArrNum]; // array

  

  cout << "nBin = " << nBin << endl;
  for (int ib =0; ib < nArrNum; ib ++ ) {
    if (states ==1) { binArr[ib] = tmpArr1s[ib]; }
    else if (states ==2) { binArr[ib] = tmpArr2s[ib]; }
    else if (states ==3) { binArr[ib] = tmpArr3s[ib]; }
    cout << ib <<"th bin = " << binArr[ib] << endl;
  }

  /////////////////////////////////////////////////////////
  //// Open RooDataFile
  /////////////////////////////////////////////////////////


  //file and ws
  TFile *fileIn[nFit][nBin];
  RooWorkspace* ws[nFit][nBin];
  // parameters 
  double alpha[nFit][nBin];
  double alphaErr[nFit][nBin];
  double n1s[nFit][nBin];
  double n1sErr[nFit][nBin];
  double sigma[nFit][nBin];
  double sigmaErr[nFit][nBin];
  double f1s[nFit][nBin];
  double f1sErr[nFit][nBin];
  double x1s[nFit][nBin];
  double x1sErr[nFit][nBin];
 
  double avg_alpha=0; double avg_n=0; double avg_sigma=0; double avg_f=0; double avg_x=0; 
  double avgParm[5]= {0.,0,0,0,0};// = getAvg(avg_alpha, avg_n, avg_sigma, avg_f, avg_x, states,szAA);
  
  char *Fit_loc[3] = {"../../../upsilonRAA5TeV/fitResults/Final_NomResult_170124/", "/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA_postCWR/UpsilonRAA_postCWR_fixed_nalphax/fitResults/FixedFit_n_alpha_x_Avg/", "/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA_postCWR/UpsilonRAA_postCWR_fixed_nalphax/fitResults/Constrain/"};
  char *Name_Fit[3] = {"All Free", "n_{Avg}, #alpha_{Avg}, x_{Avg} fixed", "n, #alpha, x constrained"};
  TString fileLoc[nFit];
  TString fitName[nFit]; 
  for(int i=0; i<nFit; i++)
  {
      fileLoc[i] = Fit_loc[i];
      fitName[i] = Name_Fit[i];
  }
 
  if(szAA=="AA"){fileLoc[1] = "/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA_postCWR/UpsilonRAA_postCWR_fixed_nalphax/fitResults/FixedFit_n_alpha_x_Avg_frompp/"; fitName[1] = "n_{Avg}, #alpha_{Avg}, x_{Avg}, f fixed"; fitName[2] = "n, #alpha, x, f constrained";}

  Int_t fitColorArr[] = { kGray+3, kRed+2, kBlue+1, kOrange+7, kGreen+3, kAzure+9, kViolet-1, kGreen+1,kBlack };

  TFile* rfp = new TFile("../AvgSigPar_PP.root","read");
  TH1D* hparn = (TH1D*) rfp->Get("hAvgn");
  TH1D* hparax = (TH1D*) rfp->Get("hAvgalphax");


  for (int ib =0; ib < nBin; ib ++ ) {
    for(int ifit = 0; ifit < nFit; ifit++){
      //// read files
      if (szAA == "PP" ) { 
        if(ifit==0) {fileIn[ifit][ib]= new TFile(Form("%sPAS_fitresults_upsilon_DoubleCB_%s_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0.root",fileLoc[ifit].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));}
        else if(ifit!=0) {fileIn[ifit][ib]= new TFile(Form("%sfitresults_upsilon_fixParm1_seed2_DoubleCB_%s_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0.root",fileLoc[ifit].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));}
      }
      else if (szAA == "AA" ) { 
        if(ifit==0) {fileIn[ifit][ib]= new TFile(Form("%sPAS_fitresults_upsilon_DoubleCB_%s_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",fileLoc[ifit].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));}
        else if(ifit!=0) {fileIn[ifit][ib]= new TFile(Form("%sfitresults_upsilon_fixParm1_seed2_DoubleCB_%s_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",fileLoc[ifit].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));}
      }
      else { cout << " Error ::: Select among PP and AA" << endl; return 0; }
      cout << ib << "th file = " << fileIn[ifit][ib]->GetName() << endl;
      if (fileIn[ifit][ib]->IsZombie()) { cout << "CANNOT open data root file\n"; return 1; }
      fileIn[ifit][ib]->cd();
      ws[ifit][ib]= (RooWorkspace*)fileIn[ifit][ib]->Get("workspace");
      //ws[ifit][ib]->Print();

      //// get parameters
      alpha[ifit][ib]=ws[ifit][ib]->var("alpha1s_1")->getVal();
      alphaErr[ifit][ib]=ws[ifit][ib]->var("alpha1s_1")->getError();
      n1s[ifit][ib]=ws[ifit][ib]->var("n1s_1")->getVal();
      n1sErr[ifit][ib]=ws[ifit][ib]->var("n1s_1")->getError();
      sigma[ifit][ib]=ws[ifit][ib]->var("sigma1s_1")->getVal();
      sigmaErr[ifit][ib]=ws[ifit][ib]->var("sigma1s_1")->getError();
      f1s[ifit][ib]=ws[ifit][ib]->var("f1s")->getVal();
      f1sErr[ifit][ib]=ws[ifit][ib]->var("f1s")->getError();
      x1s[ifit][ib]=ws[ifit][ib]->var("x1s")->getVal();
      x1sErr[ifit][ib]=ws[ifit][ib]->var("x1s")->getError();
      
      if(nFit==3 && ifit==1){
        alphaErr[ifit][ib]=hparax->GetBinContent(6);
        n1sErr[ifit][ib]=hparn->GetBinContent(7);
        x1sErr[ifit][ib] = hparax->GetBinContent(10);
      }
      //cout << ib << "th nSig1s = " << nSig1s[ifit][ib] << endl;
      //cout << ib << "th nSig2s = " << nSig2s[ifit][ib] << endl;
      //cout << ib << "th nSig3s = " << nSig3s[ifit][ib] << endl;
      //cout << ib << "th nBkg = " << nBkg[ifit][ib] << endl;
    }
  }
  //// histogram
  TH1D* h1_alpha[nFit]; 
  TH1D* h1_n[nFit]; 
  TH1D* h1_sigma[nFit]; 
  TH1D* h1_f[nFit]; 
  TH1D* h1_x[nFit]; 
  for(int ifit=0; ifit<nFit; ifit++){
    h1_alpha[ifit] = new TH1D(Form("h1_alpha%ds_%d",states,ifit+1),Form("h1_alpha%ds;|y|;#alpha",states),nBin,binArr); 
    h1_n[ifit] = new TH1D(Form("h1_n%ds_%d",states,ifit+1),Form("h1_n%ds;|y|;n",states),nBin,binArr); 
    h1_sigma[ifit] = new TH1D(Form("h1_sigma%ds_%d",states,ifit+1),Form("h1_sigma%d;|y|;#sigma",states),nBin,binArr); 
    h1_f[ifit] = new TH1D(Form("h1_f%ds_%d",states,ifit+1),Form("h1_f%ds;|y|;f",states),nBin,binArr); 
    h1_x[ifit] = new TH1D(Form("h1_x%ds_%d",states,ifit+1),Form("h1_x%ds;|y|;x",states),nBin,binArr); 
    h1_n[ifit]->GetYaxis()->SetTitleFont(32);
    h1_f[ifit]->GetYaxis()->SetTitleFont(32);
    h1_x[ifit]->GetYaxis()->SetTitleFont(32);
    for (int ib =0; ib < nBin; ib ++ ) {
      h1_alpha[ifit]->SetBinContent(ib+1,alpha[ifit][ib]);   
      h1_alpha[ifit]->SetBinError(ib+1,alphaErr[ifit][ib]);   
      h1_n[ifit]->SetBinContent(ib+1,n1s[ifit][ib]);   
      h1_n[ifit]->SetBinError(ib+1,n1sErr[ifit][ib]);   
      h1_sigma[ifit]->SetBinContent(ib+1,sigma[ifit][ib]);   
      h1_sigma[ifit]->SetBinError(ib+1,sigmaErr[ifit][ib]);   
      h1_f[ifit]->SetBinContent(ib+1,f1s[ifit][ib]);   
      h1_f[ifit]->SetBinError(ib+1,f1sErr[ifit][ib]);   
      h1_x[ifit]->SetBinContent(ib+1,x1s[ifit][ib]);   
      h1_x[ifit]->SetBinError(ib+1,x1sErr[ifit][ib]);   
    }
  }

  //// normalization
  for(int ifit=0; ifit<nFit; ifit++){
    SetHistStyle(h1_alpha[ifit],ifit, ifit);
    SetHistStyle(h1_n[ifit],ifit, ifit);
    SetHistStyle(h1_sigma[ifit],ifit, ifit);
    SetHistStyle(h1_f[ifit],ifit, ifit);
    SetHistStyle(h1_x[ifit],ifit, ifit);
  }

  int binmax[5] = { h1_alpha[0]->GetMaximumBin(),h1_n[0]->GetMaximumBin(), h1_sigma[0]->GetMaximumBin(), h1_f[0]->GetMaximumBin(), h1_x[0]->GetMaximumBin()};
  double valmax[5] = { h1_alpha[0]->GetBinContent(binmax[0]), h1_n[0]->GetBinContent(binmax[1]), h1_sigma[0]->GetBinContent(binmax[2]), h1_f[0]->GetBinContent(binmax[3]), h1_x[0]->GetBinContent(binmax[4])};
  h1_alpha[0]->GetYaxis()->SetRangeUser(0,valmax[0]*3);
  h1_n[0]->GetYaxis()->SetRangeUser(0,valmax[1]*3);
  h1_sigma[0]->GetYaxis()->SetRangeUser(0,valmax[2]*3);
  h1_f[0]->GetYaxis()->SetRangeUser(0,valmax[3]*3);
  h1_x[0]->GetYaxis()->SetRangeUser(0,valmax[4]*1.5);


  TGraphErrors *g_alpha[nFit];
  TGraphErrors *g_n[nFit];
  TGraphErrors *g_sigma[nFit]; 
  TGraphErrors *g_f[nFit];
  TGraphErrors *g_x[nFit];

  double pxtmp, pytmp, extmp, eytmp;
  double shift_x=0.03;
  double shift_x_diff=shift_x*2/3;

  for(int i=0;i<nFit;i++){
    g_alpha[i] = new TGraphErrors(h1_alpha[i]);
    g_n[i] = new TGraphErrors(h1_n[i]);
    g_sigma[i] = new TGraphErrors(h1_sigma[i]);
    g_f[i] = new TGraphErrors(h1_f[i]);
    g_x[i] = new TGraphErrors(h1_x[i]);
    for(int j=0;j<g_alpha[i]->GetN();j++){
      pxtmp=0; pytmp=0; extmp=0; eytmp=0;
      g_alpha[i]->GetPoint(j,pxtmp,pytmp);
      extmp=g_alpha[i]->GetErrorX(j);
      eytmp=g_alpha[i]->GetErrorY(j);
      g_alpha[i]->SetPoint(j,pxtmp-shift_x+shift_x_diff*i,pytmp);
      g_alpha[i]->SetPointError(j,0,eytmp);
    }
    for(int j=0;j<g_x[i]->GetN();j++){
      pxtmp=0; pytmp=0; extmp=0; eytmp=0;
      g_x[i]->GetPoint(j,pxtmp,pytmp);
      extmp=g_x[i]->GetErrorX(j);
      eytmp=g_x[i]->GetErrorY(j);
      g_x[i]->SetPoint(j,pxtmp-shift_x+shift_x_diff*i,pytmp);
      g_x[i]->SetPointError(j,0,eytmp);
    }
    for(int j=0;j<g_f[i]->GetN();j++){
      pxtmp=0; pytmp=0; extmp=0; eytmp=0;
      g_f[i]->GetPoint(j,pxtmp,pytmp);
      extmp=g_f[i]->GetErrorX(j);
      eytmp=g_f[i]->GetErrorY(j);
      g_f[i]->SetPoint(j,pxtmp-shift_x+shift_x_diff*i,pytmp);
      g_f[i]->SetPointError(j,0,eytmp);
    }
    for(int j=0;j<g_sigma[i]->GetN();j++){
      pxtmp=0; pytmp=0; extmp=0; eytmp=0;
      g_sigma[i]->GetPoint(j,pxtmp,pytmp);
      extmp=g_sigma[i]->GetErrorX(j);
      eytmp=g_sigma[i]->GetErrorY(j);
      g_sigma[i]->SetPoint(j,pxtmp-shift_x+shift_x_diff*i,pytmp);
      g_sigma[i]->SetPointError(j,0,eytmp);
    }
    for(int j=0;j<g_n[i]->GetN();j++){
      pxtmp=0; pytmp=0; extmp=0; eytmp=0;
      g_n[i]->GetPoint(j,pxtmp,pytmp);
      extmp=g_n[i]->GetErrorX(j);
      eytmp=g_n[i]->GetErrorY(j);
      g_n[i]->SetPoint(j,pxtmp-shift_x+shift_x_diff*i,pytmp);
      g_n[i]->SetPointError(j,0,eytmp);
    }
  }
  g_alpha[0]->SetMinimum(0);
  g_alpha[0]->SetMaximum(valmax[0]*3);
  g_n[0]->SetMinimum(0);
  g_n[0]->SetMaximum(valmax[1]*3);
  g_sigma[0]->SetMinimum(0);
  g_sigma[0]->SetMaximum(valmax[2]*3);
  g_f[0]->SetMinimum(0);
  g_f[0]->SetMaximum(valmax[3]*3);
  g_x[0]->SetMinimum(0);
  g_x[0]->SetMaximum(valmax[4]*2);

  g_alpha[0]->GetXaxis()->SetTitle("|y|^{#varUpsilon}");
  g_alpha[0]->GetYaxis()->SetTitle("#alpha");
  g_n[0]->GetXaxis()->SetTitle("|y|^{#varUpsilon}");
  g_n[0]->GetYaxis()->SetTitle("n");
  g_sigma[0]->GetXaxis()->SetTitle("|y|^{#varUpsilon}");
  g_sigma[0]->GetYaxis()->SetTitle("#sigma");
  g_f[0]->GetXaxis()->SetTitle("|y|^{#varUpsilon}");
  g_f[0]->GetYaxis()->SetTitle("f");
  g_x[0]->GetXaxis()->SetTitle("|y|^{#varUpsilon}");
  g_x[0]->GetYaxis()->SetTitle("x");

  g_alpha[0]->GetXaxis()->CenterTitle();
  g_alpha[0]->GetYaxis()->CenterTitle();
  g_n[0]->GetXaxis()->CenterTitle();
  g_n[0]->GetYaxis()->CenterTitle();
  g_sigma[0]->GetXaxis()->CenterTitle();
  g_sigma[0]->GetYaxis()->CenterTitle();
  g_f[0]->GetXaxis()->CenterTitle();
  g_f[0]->GetYaxis()->CenterTitle();
  g_x[0]->GetXaxis()->CenterTitle();
  g_x[0]->GetYaxis()->CenterTitle();



 

  //// actual draw
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextFont(32);
  latex->SetTextAlign(12);
  latex->SetTextSize(0.035);

  double pos_y_diff = 0.05; 
  double pos_y = 0.837; 
  double pos_x = 0.51; 
  double line_x_diff = 0.05;
  
  double leg_posx1 = 0.28;
  double leg_posy1 = 0.52;
  double leg_posx2 = 0.53;
  double leg_posy2 = 0.70;

  //Get Integrated bin
  TFile *fInt;
  if(szAA=="PP" && nFit!=1) fInt = new TFile(Form("%sfitresults_upsilon_fixParm1_seed2_DoubleCB_%s_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0.root",fileLoc[nFit-1].Data(),szAA.Data()));
  else if(szAA=="PP" && nFit==1) fInt = new TFile(Form("%sPAS_fitresults_upsilon_DoubleCB_%s_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0.root",fileLoc[nFit-1].Data(),szAA.Data()));
  if(szAA=="AA" && nFit!=1) fInt = new TFile(Form("%sfitresults_upsilon_fixParm1_seed2_DoubleCB_%s_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",fileLoc[nFit-1].Data(),szAA.Data()));
  else if(szAA=="AA" && nFit==1) fInt = new TFile(Form("%sPAS_fitresults_upsilon_DoubleCB_%s_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",fileLoc[nFit-1].Data(),szAA.Data()));
  RooWorkspace *ws_int = (RooWorkspace*) fInt->Get("workspace");
  double alpha_int = ws_int->var("alpha1s_1")->getVal();
  double alpha_int_err = ws_int->var("alpha1s_1")->getError();
  double n_int = ws_int->var("n1s_1")->getVal();
  double n_int_err = ws_int->var("n1s_1")->getError();
  double sigma_int = ws_int->var("sigma1s_1")->getVal();
  double sigma_int_err = ws_int->var("sigma1s_1")->getError();
  double f_int = ws_int->var("f1s")->getVal();
  double f_int_err = ws_int->var("f1s")->getError();
  double x_int = ws_int->var("x1s")->getVal();
  double x_int_err = ws_int->var("x1s")->getError();

  double IntVal[5] = {alpha_int,n_int,sigma_int,f_int,x_int};
  double IntValErr[5] = {alpha_int_err,n_int_err,sigma_int_err,f_int_err,x_int_err};

  TLegend* fitleg = new TLegend(leg_posx1,leg_posy1,leg_posx2,leg_posy2); fitleg->SetTextSize(19);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);

  TCanvas* c_alpha = new TCanvas("c_alpha","c_alpha",600,600);
  c_alpha->cd();
  for(int ifit=0; ifit<nFit; ifit++){ if(ifit==0) g_alpha[ifit]->Draw("AP"); else if(ifit>0) g_alpha[ifit]->Draw("P");}
  latex->SetTextColor(kBlack);
  latex->DrawLatex(pos_x,pos_y,Form("%s #Upsilon(%dS)",szAA.Data(),states));
  latex->DrawLatex(pos_x,pos_y-pos_y_diff*0.8,"param : #alpha");
  
/*  latex->DrawLatex(leg_posx1+line_x_diff,leg_posy1-pos_y_diff*0.6,"Int. Bin");
  latex->DrawLatex(leg_posx1+line_x_diff,leg_posy1-pos_y_diff*1.5,"Avg. Val");
  dashedLine(1.64,valmax[0]*3/1.9,1.72,valmax[0]*3/1.9,1,1);
  dashedLine(1.64,valmax[0]*3/2.15,1.72,valmax[0]*3/2.15,2,2);
  jumSun(0,alpha_int,2.4,alpha_int,1,1);
  jumSun(0,avgParm[0],2.4,avgParm[0],2,1);
*/
  for(int ifit=0;ifit<nFit; ifit++){
    if(ifit==0) fitleg->AddEntry(h1_alpha[ifit],Form("%s fit",fitName[ifit].Data()),"pe");
    else if(ifit!=0) fitleg->AddEntry(h1_alpha[ifit],Form("%s fit",fitName[ifit].Data()),"pe");
  }
  fitleg->Draw("same");
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.1);
  c_alpha->SetTicks(1,1);
  c_alpha->Modified();
  c_alpha->Update();
  CMS_lumi_internal( c_alpha, iPeriod, iPos );
  c_alpha->SaveAs(Form("sigparam_all/Constrain_rap_alpha_Upsilon%ds_%s_nFit%d.pdf",states,szAA.Data(),nFit));
  
  TCanvas* c_n = new TCanvas("c_n","c_n",600,600);
  c_n->cd();
  for(int ifit=0; ifit<nFit; ifit++){ if(ifit==0) g_n[ifit]->Draw("AP"); else if(ifit>0) g_n[ifit]->Draw("P");}
  latex->SetTextColor(kBlack);
  latex->DrawLatex(pos_x,pos_y,Form("%s #Upsilon(%dS)",szAA.Data(),states));
  latex->DrawLatex(pos_x,pos_y-pos_y_diff*0.8,"param : n");/*
  latex->DrawLatex(leg_posx1+line_x_diff,leg_posy1-pos_y_diff*0.6,"Int. Bin");
  latex->DrawLatex(leg_posx1+line_x_diff,leg_posy1-pos_y_diff*1.5,"Avg. Val");
  dashedLine(1.64,valmax[1]*3/1.9,1.72,valmax[1]*3/1.9,1,1);
  dashedLine(1.64,valmax[1]*3/2.15,1.72,valmax[1]*3/2.15,2,2);
  jumSun(0,n_int,2.4,n_int,1,1);
  jumSun(0,avgParm[1],2.4,avgParm[1],2,1);
*/
  fitleg->Draw("same");
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.1);
  c_n->SetTicks(1,1);
  c_n->Modified();
  c_n->Update();
  CMS_lumi_internal( c_n, iPeriod, iPos ); 
  c_n->SaveAs(Form("sigparam_all/Constrain_rap_n_Upsilon%ds_%s_nFit%d.pdf",states,szAA.Data(),nFit));
  
  TCanvas* c_sigma = new TCanvas("c_sigma","c_sigma",600,600);
  c_sigma->cd();
  for(int ifit=0; ifit<nFit; ifit++){ if(ifit==0) g_sigma[ifit]->Draw("AP"); else if(ifit>0) g_sigma[ifit]->Draw("P");}
  latex->SetTextColor(kBlack);
  latex->DrawLatex(pos_x,pos_y,Form("%s #Upsilon(%dS)",szAA.Data(),states));
  latex->DrawLatex(pos_x,pos_y-pos_y_diff*0.8,"param : #sigma");
  /*
  latex->DrawLatex(leg_posx1+line_x_diff,leg_posy1-pos_y_diff*0.6,"Int. Bin");
  latex->DrawLatex(leg_posx1+line_x_diff,leg_posy1-pos_y_diff*1.5,"Avg. Val");
  dashedLine(1.64,valmax[2]*3/1.9,1.72,valmax[2]*3/1.9,1,1);
  dashedLine(1.64,valmax[2]*3/2.15,1.72,valmax[2]*3/2.15,2,2);
  jumSun(0,sigma_int,2.4,sigma_int,1,1);
  jumSun(0,avgParm[2],2.4,avgParm[2],2,1);
*/
  fitleg->Draw("same");
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.1);
  c_sigma->SetTicks(1,1);
  c_sigma->Modified();
  c_sigma->Update();
  CMS_lumi_internal( c_sigma, iPeriod, iPos ); 
  c_sigma->SaveAs(Form("sigparam_all/Constrain_rap_sigma_Upsilon%ds_%s_nFit%d.pdf",states,szAA.Data(),nFit));
  
  TCanvas* c_f = new TCanvas("c_f","c_f",600,600);
  c_f->cd();
  for(int ifit=0; ifit<nFit; ifit++){ if(ifit==0) g_f[ifit]->Draw("AP"); else if(ifit>0) g_f[ifit]->Draw("P");}
  latex->SetTextColor(kBlack);
  latex->DrawLatex(pos_x,pos_y,Form("%s #Upsilon(%dS)",szAA.Data(),states));
  latex->DrawLatex(pos_x,pos_y-pos_y_diff*0.8,"param : f");
  /*
  latex->DrawLatex(leg_posx1+line_x_diff,leg_posy1-pos_y_diff*0.6,"Int. Bin");
  latex->DrawLatex(leg_posx1+line_x_diff,leg_posy1-pos_y_diff*1.5,"Avg. Val");
  dashedLine(1.64,valmax[3]*3/1.9,1.72,valmax[3]*3/1.9,1,1);
  dashedLine(1.64,valmax[3]*3/2.15,1.72,valmax[3]*3/2.15,2,2);
  jumSun(0,f_int,2.4,f_int,1,1);
  jumSun(0,avgParm[3],2.4,avgParm[3],2,1);
*/

  fitleg->Draw("same");
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.1);
  c_f->SetTicks(1,1);
  c_f->Modified();
  c_f->Update();
  CMS_lumi_internal( c_f, iPeriod, iPos ); 
  c_f->SaveAs(Form("sigparam_all/Constrain_rap_f_Upsilon%ds_%s_nFit%d.pdf",states,szAA.Data(),nFit));
  
  TCanvas* c_x = new TCanvas("c_x","c_x",600,600);
  c_x->cd();
  for(int ifit=0; ifit<nFit; ifit++){ if(ifit==0) g_x[ifit]->Draw("AP"); else if(ifit>0) g_x[ifit]->Draw("P same");}
  latex->SetTextColor(kBlack);
  latex->DrawLatex(pos_x,pos_y,Form("%s #Upsilon(%dS)",szAA.Data(),states));
  latex->DrawLatex(pos_x,pos_y-pos_y_diff*0.8,"param : x");
  /*
  latex->DrawLatex(leg_posx1+line_x_diff,leg_posy1-pos_y_diff*0.6,"Int. Bin");
  latex->DrawLatex(leg_posx1+line_x_diff,leg_posy1-pos_y_diff*1.5,"Avg. Val");
  dashedLine(1.64,valmax[4]*3/1.9,1.72,valmax[4]*3/1.9,1,1);
  dashedLine(1.64,valmax[4]*3/2.15,1.72,valmax[4]*3/2.15,2,2);
  jumSun(0,x_int,2.4,x_int,1,1);
  jumSun(0,avgParm[4],2.4,avgParm[4],2,1);
  */

  fitleg->Draw("same");
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.1);
  c_x->SetTicks(1,1);
  c_x->Modified();
  c_x->Update();
  CMS_lumi_internal( c_x, iPeriod, iPos ); 
  c_x->SaveAs(Form("sigparam_all/Constrain_rap_x_Upsilon%ds_%s_nFit%d.pdf",states,szAA.Data(),nFit));
  
  TFile *wf = new TFile(Form("sigparam_all/f_rap_Ups%ds_%s_wf.root",states,szAA.Data()),"recreate");
  TH1D* hInt = new TH1D("hInt",";Val;;",5,0,5);
  TH1D* hAvg = new TH1D("hAvg",";Val;;",5,0,5);

  for(int i=0;i<5;i++){
    hInt->SetBinContent(i+1,IntVal[i]);  
    hAvg->SetBinContent(i+1,avgParm[i]); 
  } 
   
  hInt->GetXaxis()->SetBinLabel(1,"#alpha");
  hInt->GetXaxis()->SetBinLabel(2,"n");
  hInt->GetXaxis()->SetBinLabel(3,"#sigma");
  hInt->GetXaxis()->SetBinLabel(4,"f");
  hInt->GetXaxis()->SetBinLabel(5,"x");
  hAvg->GetXaxis()->SetBinLabel(1,"#alpha");
  hAvg->GetXaxis()->SetBinLabel(2,"n");
  hAvg->GetXaxis()->SetBinLabel(3,"#sigma");
  hAvg->GetXaxis()->SetBinLabel(4,"f");
  hAvg->GetXaxis()->SetBinLabel(5,"x");

  wf->cd();
  hInt->Write();
  hAvg->Write();
  return 0;
}

double *getAvg(double avg_alpha, double avg_n , double avg_sigma, double avg_f, double avg_x, int states, TString szAA)
{
  int nFit = 2; 
  char *Fit_loc[2] = {" ../fitResults/AllParmFreeFit_pasVer/", "../MCFixedFit/"};
  char *Name_Fit[2] = {"Data", "MC"};
  TString fileLoc[nFit];
  TString fitName[nFit]; 
  for(int i=0; i<nFit; i++)
  {
      fileLoc[i] = Fit_loc[i];
      fitName[i] = Name_Fit[i];
  }
  
  double tmpArr1s[7] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  double tmpArr2s[4] = {0.0, 0.8, 1.6, 2.4};
  double tmpArr3s[3] = {0.0, 1.2, 2.4};
  
  int tmpBin;
  if ( states ==1) {
    cout << " ***** 1S *****" << endl; tmpBin = 6;  
  }else if (states ==2){
    cout << " ***** 2S *****" << endl; tmpBin = 3;  
  }else if (states ==3){
    cout << " ***** 3S *****" << endl; tmpBin = 2;  
  }else {
    cout << " Error ::: Select among 1S, 2S, and 3S" << endl; 
  }
 
  double tmpArr1s_[7] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  double tmpArr2s_[4] = {0.0, 0.8, 1.6, 2.4};
  double tmpArr3s_[3] = {0.0, 1.2, 2.4};

  int tmpBin_;
  if ( states ==1) {
    cout << " ***** 1S *****" << endl; tmpBin_ = 6;
  }else if (states ==2){
    cout << " ***** 2S *****" << endl; tmpBin_ = 3;
  }else if (states ==3){
    cout << " ***** 3S *****" << endl; tmpBin_ = 2;
  }else {
    cout << " Error ::: Select among 1S, 2S, and 3S" << endl;
  }
  const int nStates = 3; 
  const int nBin = tmpBin; // number of bin 
  const int nArrNum = nBin+1; // number of array
  double binArr[nArrNum]; // array
  
  const int nStates_ = 3; 
  const int nBin_ = tmpBin_; // number of bin 
  const int nArrNum_ = nBin_+1; // number of array
  double binArr_[nArrNum_]; // array
  
  TFile *fileIn[nBin];
  RooWorkspace* ws[nBin];
  double alpha[nBin];
  double n1s[nBin];
  double sigma[nBin];
  double f1s[nBin];
  double x1s[nBin];

  double alphaErr=0;
  double nErr=0;
  double xErr=0;
  double fErr=0;
  double sigmaErr=0;

  cout << "nBin = " << nBin << endl;
  for (int ib =0; ib < nArrNum; ib ++ ) {
    if (states ==1) { binArr[ib] = tmpArr1s[ib]; }
    else if (states ==2) { binArr[ib] = tmpArr2s[ib]; }
    else if (states ==3) { binArr[ib] = tmpArr3s[ib]; }
    cout << ib <<"th bin = " << binArr[ib] << endl;
  }
  for (int ib =0; ib < nBin; ib ++ ) {
    //// read files
    if (szAA == "PP" ) fileIn[ib]= new TFile(Form("%sPAS_fitresults_upsilon_DoubleCB_%s_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0.root",fileLoc[0].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));
    else if (szAA == "AA" ) fileIn[ib]= new TFile(Form("%sPAS_fitresults_upsilon_DoubleCB_%s_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",fileLoc[0].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));
    else { cout << " Error ::: Select among PP and AA" << endl;  }
    cout << ib << "th file = " << fileIn[ib]->GetName() << endl;
    if (fileIn[ib]->IsZombie()) { cout << "CANNOT open data root file\n"; }
    fileIn[ib]->cd();
    ws[ib]= (RooWorkspace*)fileIn[ib]->Get("workspace");

    //// get parameters
    alpha[ib]=ws[ib]->var("alpha1s_1")->getVal();
    n1s[ib]=ws[ib]->var("n1s_1")->getVal();
    sigma[ib]=ws[ib]->var("sigma1s_1")->getVal();
    f1s[ib]=ws[ib]->var("f1s")->getVal();
    x1s[ib]=ws[ib]->var("x1s")->getVal();

    avg_alpha = avg_alpha + alpha[ib];
    avg_n = avg_n + n1s[ib];
    avg_sigma = avg_sigma + sigma[ib];
    avg_f = avg_f + f1s[ib];
    avg_x = avg_x + x1s[ib];
    cout << "avg_alpha : " << avg_alpha << endl;  
  }
 
  avg_alpha = avg_alpha/(nBin);
  avg_n = avg_n/(nBin);
  cout << "avg_n : " << avg_n << endl;  
  avg_sigma  = avg_sigma/(nBin);
  avg_f  = avg_f/(nBin);
  avg_x  = avg_x/(nBin);
  
  for (int ib =0; ib < nBin; ib ++ ) {
    alphaErr += (alpha[ib]-avg_alpha)*(alpha[ib]-avg_alpha);
  }
  alphaErr = TMath::Sqrt( (alphaErr) / (nBin*(nBin-1)) );


  double* avgParm = new double[10];
  avgParm[0] = avg_alpha; avgParm[1] =  avg_n; avgParm[2] = avg_sigma; avgParm[3] = avg_f; avgParm[4] = avg_x;
  return avgParm;
}
