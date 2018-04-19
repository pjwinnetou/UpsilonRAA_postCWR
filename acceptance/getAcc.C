#include <iostream>
#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TROOT.h"
#include "TFile.h"
#include "../rootFitHeaders.h"
#include "../cutsAndBin.h"
#include <fstream>
//#include "SimplePtFit.C"
using namespace std;

valErr getYield(int state= 0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0, int cLow=0, int cHigh=0,   float dphiEp2Low=0,  float dphiEp2High=0) ;
void getMaxTH1D_eight ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0,TH1D* h5=0,TH1D* h6=0, TH1D* h7=0,TH1D* h8=0) ;
void getMaxTH1D_four ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0) ;

void getAcc(int state= 2, int collId= kAADATA) {

  cout<<"dN/dp_{T} starting "<< collId <<" state: " << endl;

  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1);

  TH1::SetDefaultSumw2();

  TString fCollId;
  if(collId == kPPDATA) fCollId = "PP";
  else if(collId == kAADATA) fCollId = "AA";

  //// modify by hand according to the pt range of the sample
  int nPtBins=0;
  double* ptBin;
  int nPtBinsMC=0;
  double* ptBinMC;
  int nYBins=0;
  double* yBin;

  if ( state == 1 ) {
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S;    yBin = yBin1S;
  }
  else if ( state == 2 ) {
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nYBins = nYBins1S;    yBin = yBin1S;
  }
  else if ( state == 3 ) {
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nYBins = nYBins1S;    yBin = yBin1S;
  }
  //For 2D plot :
  float ptLo = 0; float ptHi = 30;
  float yLo  = 0; float  yHi = 2.4;

  // Get MC :
  float massLow = 8; float massHigh = 14;
  double ptMin = ptLo; double ptMax = ptHi;
  double yMin = yLo;     double yMax = yHi;

  TH1D* hintGen = new TH1D("hintGen",";p_{T}(GeV/c);",1,0,30);
  TH1D *hintGenAcc  = new TH1D("hintGenAcc","; p_{T} (GeV/c) ; ", 1,0,30);

  TH1D* hptGen=new TH1D("hptGen",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D *hptGenAcc  = new TH1D("hptGenAcc","; p_{T} (GeV/c) ; ", nPtBins, ptBin);

  TH1D* hrapGen=new TH1D("hrapGen",";|y|;",nYBins,yBin);
  TH1D *hrapGenAcc = new TH1D("hrapGenAcc","; |y| ; ", nYBins, yBin);

  TH1D* hptGenAccSys[8];
  TH1D* hrapGenAccSys[8];
  TH1D* hintGenAccSys[8];
  TH1D* hptGenAccSys_allGen[8];
  TH1D* hrapGenAccSys_allGen[8];
  TH1D* hintGenAccSys_allGen[8];
  for(int i=0;i<8;i++){
    hptGenAccSys[i] = (TH1D*) hptGen->Clone(Form("hptGenAccSys_fill%d",i+1));
    hrapGenAccSys[i] = (TH1D*) hrapGen->Clone(Form("hrapGenAccSys_fill%d",i+1));
    hintGenAccSys[i] = (TH1D*) hintGen->Clone(Form("hintGenAccSys_fill%d",i+1));
    hptGenAccSys_allGen[i] = (TH1D*) hptGen->Clone(Form("hptGenAccSys_allGen_fill%d",i+1));
    hrapGenAccSys_allGen[i] = (TH1D*) hrapGen->Clone(Form("hrapGenAccSys_allGen_fill%d",i+1));
    hintGenAccSys_allGen[i] = (TH1D*) hintGen->Clone(Form("hintGenAccSys_allGen_fill%d",i+1));
  }

  TChain *mmGen = new TChain("mmGen");

  if(state==1){
    if(collId==kPPDATA){
      mmGen->Add("../skimmedFiles/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161281226_.root");
    }
    else if(collId==kAADATA){
      mmGen->Add("../skimmedFiles/yskimAA_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161281233_.root");
    }
  }
  else if(state==2){
    if(collId==kPPDATA){
      mmGen->Add("../skimmedFiles/yskimPP_MC_Ups2S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161281228_.root");
    }
    else if(collId==kAADATA){
      mmGen->Add("../skimmedFiles/yskimAA_MC_Ups2S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161281234_.root");
    }
  }
  else if(state==3){
    if(collId==kPPDATA){
      mmGen->Add("../skimmedFiles/yskimPP_MC_Ups3S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161281230_.root");
    }
    else if(collId==kAADATA){
      mmGen->Add("../skimmedFiles/yskimAA_MC_Ups3S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161281235_.root");
    }
  }

  DiMuon  dmGen;
  TBranch *b_dmGen;

  mmGen->SetBranchAddress("mmGen",&dmGen,&b_dmGen);

  float ptWeight;
  float ptWeight_func;
  double ptWeight_func_sys[8];
  
  int nVar = 4;
  if(state==1) nVar=8;
  
  // test code
  double parm_AA1[4] = {255.074, -93.4016, 44.2256, -4.81048};
  double parm_AA2[2] = {0.778896, 0.0209981};
  double parm_PP1[4] = {200.759, 7.09569, 25.3727, -4.51979};
  double parm_PP2[2] = {0.569212, 0.0637386};

  double parm_AA1_up[4] = {320.738, -77.2557, 48.3587, -4.21744};
  double parm_AA2_up[2] = {1.08026, 0.0468831};
  double parm_PP1_up[4] = {254.905, 17.5645, 26.6928, -4.04741};
  double parm_PP2_up[2] = {0.592215, 0.066485};

  double parm_AA1_do[4] = {189.409, -109.547, 40.0926, -5.40352};
  double parm_AA2_do[2] = {0.477534, -0.00488679};
  double parm_PP1_do[4] = {146.613, -3.37308, 24.0526, -4.99218};
  double parm_PP2_do[2] = { 0.546209, 0.0609922};

  TFile *fweight;
  if(state==1) {fweight = new TFile(Form("../CompareDataMC/WeightedFcN_fit/ratioDataMC_%s_DATA_1s_20180418.root",fCollId.Data()),"read");}
  else if(state!=1) {fweight = new TFile(Form("../CompareDataMC/WeightedFcN_fit/ratioDataMC_%s_DATA_2s_20180418.root",fCollId.Data()),"read");}

  TF1 *wfnom = (TF1*) fweight -> Get("dataMC_Ratio_norm");
  TF1 *wfAup = (TF1*) fweight -> Get("dataMC_Ratio_Ap");
  TF1 *wfAdo = (TF1*) fweight -> Get("dataMC_Ratio_Am");
  TF1 *wfBup = (TF1*) fweight -> Get("dataMC_Ratio_Bp");
  TF1 *wfBdo = (TF1*) fweight -> Get("dataMC_Ratio_Bm");
  TF1 *wfCup; 
  TF1 *wfCdo;  
  TF1 *wfDup;  
  TF1 *wfDdo; 

  if(state==1){
    wfCup = (TF1*) fweight->Get("dataMC_Ratio_Cp");
    wfCdo = (TF1*) fweight->Get("dataMC_Ratio_Cm");
    wfDup = (TF1*) fweight->Get("dataMC_Ratio_Dp");
    wfDdo = (TF1*) fweight->Get("dataMC_Ratio_Dm");
 /* 
    wfnom->SetParameters(parm_AA1[0],parm_AA1[1],parm_AA1[2],parm_AA1[3]);
    wfAup->SetParameters(parm_AA1_up[0],parm_AA1[1],parm_AA1[2],parm_AA1[3]);
    wfAdo->SetParameters(parm_AA1_do[0],parm_AA1[1],parm_AA1[2],parm_AA1[3]);
    wfBup->SetParameters(parm_AA1[0],parm_AA1_up[1],parm_AA1[2],parm_AA1[3]);
    wfBdo->SetParameters(parm_AA1[0],parm_AA1_do[1],parm_AA1[2],parm_AA1[3]);
    wfCup->SetParameters(parm_AA1[0],parm_AA1[1],parm_AA1_up[2],parm_AA1[3]);
    wfCdo->SetParameters(parm_AA1[0],parm_AA1[1],parm_AA1_do[2],parm_AA1[3]);
    wfDup->SetParameters(parm_AA1[0],parm_AA1[1],parm_AA1[2],parm_AA1_up[3]);
    wfDdo->SetParameters(parm_AA1[0],parm_AA1[1],parm_AA1[2],parm_AA1_do[3]); 
    wfnom->SetParameters(parm_AA1[0],parm_AA1[1],parm_AA1[2],parm_AA1[3]);
   */ 
  }

  for(int iev=0; iev<mmGen->GetEntries() ; ++iev)
  {
    mmGen->GetEntry(iev);    
    ptWeight = dmGen.weight0;
    ptWeight_func = wfnom->Eval(dmGen.pt);
    ptWeight_func_sys[0] = wfAup->Eval(dmGen.pt);
    ptWeight_func_sys[1] = wfAdo->Eval(dmGen.pt);
    ptWeight_func_sys[2] = wfBup->Eval(dmGen.pt);
    ptWeight_func_sys[3] = wfBdo->Eval(dmGen.pt);

    if(state==1){
      ptWeight_func_sys[4] = wfCup->Eval(dmGen.pt);
      ptWeight_func_sys[5] = wfCdo->Eval(dmGen.pt);
      ptWeight_func_sys[6] = wfDup->Eval(dmGen.pt);
      ptWeight_func_sys[7] = wfDdo->Eval(dmGen.pt);
    }

    if( !( (dmGen.pt > ptMin) && (dmGen.pt < ptMax) && (fabs(dmGen.y) > yMin) && (fabs(dmGen.y) < yMax) ) ) continue;
    hptGen -> Fill(dmGen.pt, ptWeight * ptWeight_func);
    hrapGen -> Fill(dmGen.y, ptWeight * ptWeight_func);
    hintGen -> Fill(dmGen.pt, ptWeight * ptWeight_func);
    
    //Sys
    for(int i=0;i<nVar;i++){
      hptGenAccSys_allGen[i]->Fill(dmGen.pt, ptWeight * ptWeight_func_sys[i]);
      hrapGenAccSys_allGen[i]->Fill(dmGen.y, ptWeight * ptWeight_func_sys[i]);
      hintGenAccSys_allGen[i]->Fill(dmGen.pt, ptWeight * ptWeight_func_sys[i]);
    }

    if( !(dmGen.pt1>4 && dmGen.pt2>4 && (fabs(dmGen.eta1)<2.4) && (fabs(dmGen.eta2)<2.4)) ) continue;
    hptGenAcc -> Fill(dmGen.pt, ptWeight * ptWeight_func);
    hrapGenAcc -> Fill(dmGen.y, ptWeight * ptWeight_func);
    hintGenAcc -> Fill(dmGen.pt, ptWeight * ptWeight_func);
 
    //Sys
    for(int i=0;i<nVar;i++){
      hptGenAccSys[i] -> Fill(dmGen.pt, ptWeight * ptWeight_func_sys[i]);
      hrapGenAccSys[i] -> Fill(dmGen.y, ptWeight * ptWeight_func_sys[i]);
      hintGenAccSys[i] -> Fill(dmGen.pt, ptWeight * ptWeight_func_sys[i]);
    }
  }

  
  /////////////////////////////////// Style ///////////////////////////////////
 
  TH1D* hptAcc  = (TH1D*) hptGenAcc  -> Clone(Form("hptAcc%s%dS",fCollId.Data(),state));
  TH1D* hrapAcc = (TH1D*) hrapGenAcc -> Clone(Form("hrapAcc%s%dS",fCollId.Data(),state));
  TH1D* hintAcc = (TH1D*) hintGenAcc -> Clone(Form("hIntAcc%s%dS",fCollId.Data(),state));
  
  hptAcc->Divide(hptGen); 
  hrapAcc->Divide(hrapGen); 
  hintAcc->Divide(hintGen); 
  
  handsomeTH1(hptAcc,1);     
  handsomeTH1(hrapAcc,1);   
  handsomeTH1(hintAcc,1);
  
  TH1D* hptAccSys[10];
  TH1D* hrapAccSys[10];
  TH1D* hintAccSys[10];
  
  for(int i=1; i<=nVar; i++){
    hptAccSys[i]  = (TH1D*) hptGenAccSys[i-1]  -> Clone(Form("hptAccSys%s%dS_%dvar",fCollId.Data(),state,i));
    hrapAccSys[i] = (TH1D*) hrapGenAccSys[i-1] -> Clone(Form("hrapAccSys%s%dS_%dvar",fCollId.Data(),state,i));
    hintAccSys[i] = (TH1D*) hintGenAccSys[i-1] -> Clone(Form("hintAccSys%s%dS_%dvar",fCollId.Data(),state,i));
    
    hptAccSys[i]  -> Divide(hptGenAccSys_allGen[i-1]);
    hrapAccSys[i] -> Divide(hrapGenAccSys_allGen[i-1]);
    hintAccSys[i] -> Divide(hintGenAccSys_allGen[i-1]);
    
    hptAccSys[i]  -> Add(hptAcc,-1);   hptAccSys[i]   -> Divide(hptAcc);
    hrapAccSys[i] -> Add(hrapAcc,-1);  hrapAccSys[i]  -> Divide(hrapAcc);
    hintAccSys[i] -> Add(hintAcc,-1);  hintAccSys[i]  -> Divide(hintAcc);

    handsomeTH1(hptAccSys[i],1);
    handsomeTH1(hrapAccSys[i],1);
    handsomeTH1(hintAccSys[i],1);
  }  
   
  hptAccSys[0] = (TH1D*) hptAcc -> Clone(Form("hptSys%s",fCollId.Data()));  hptAccSys[0]->Reset();
  hrapAccSys[0] = (TH1D*) hrapAcc -> Clone(Form("hrapSys%s",fCollId.Data()));  hrapAccSys[0]->Reset();
  if(collId == kPPDATA){ hintAccSys[0] = (TH1D*) hintAcc -> Clone(Form("hcentSys%s",fCollId.Data()));  hintAccSys[0]->Reset();}
  else if(collId == kAADATA){ hintAccSys[0] = (TH1D*) hintAcc -> Clone(Form("hcentSys%s_int",fCollId.Data()));  hintAccSys[0]->Reset();}

  if(state==1){
    getMaxTH1D_eight(hptAccSys[0], hptAccSys[1], hptAccSys[2], hptAccSys[3], hptAccSys[4],hptAccSys[5], hptAccSys[6], hptAccSys[7], hptAccSys[8]);
    getMaxTH1D_eight(hrapAccSys[0], hrapAccSys[1], hrapAccSys[2], hrapAccSys[3], hrapAccSys[4], hrapAccSys[5], hrapAccSys[6], hrapAccSys[7], hrapAccSys[8]);
    getMaxTH1D_eight(hintAccSys[0], hintAccSys[1], hintAccSys[2], hintAccSys[3], hintAccSys[4], hintAccSys[5], hintAccSys[6], hintAccSys[7], hintAccSys[8]);
  }
  else if(state!=1){
    getMaxTH1D_four(hptAccSys[0], hptAccSys[1], hptAccSys[2], hptAccSys[3], hptAccSys[4]);
    getMaxTH1D_four(hrapAccSys[0], hrapAccSys[1], hrapAccSys[2], hrapAccSys[3], hrapAccSys[4]);
    getMaxTH1D_four(hintAccSys[0], hintAccSys[1], hintAccSys[2], hintAccSys[3], hintAccSys[4]);
  }

  handsomeTH1(hptGen,1);
  handsomeTH1(hrapGen,1);
  handsomeTH1(hintGen,1);
  handsomeTH1(hptGenAcc,2);
  handsomeTH1(hrapGenAcc,2);
  handsomeTH1(hintGenAcc,2);

  //////////////////////////////////////////////HISTOGRAM/////////////////////////////////////////////////////////
  
  hptAcc->GetYaxis()->SetRangeUser(0,1.05);
  hrapAcc->GetYaxis()->SetRangeUser(0,1.05);
  hintAcc->GetYaxis()->SetRangeUser(0,1.05);

  TCanvas* c1 =  new TCanvas("c1","", 800, 600);
  c1->Divide(1,2);
  c1->cd(1);
  hptGen->Draw();
  hptGenAcc->Draw("same"); 
  c1->cd(2);
  hptAcc->Draw();

  TCanvas* c2 =  new TCanvas("c2","", 800, 600);
  c2->Divide(1,2);
  c2->cd(1);
  hrapGen->Draw();
  hrapGenAcc->Draw("same"); 
  c2->cd(2);
  hrapAcc->Draw();

  TCanvas* c3 =  new TCanvas("c3","", 800, 600);
  c3->Divide(1,2);
  c3->cd(1);
  hintGen->Draw();
  hintGenAcc->Draw("same"); 
  c3->cd(2);
  hintAcc->Draw();

  TCanvas* c4 =  new TCanvas("c3","", 800, 600);
//  hptAccSys[1]->Draw();
  hptGenAccSys[0]->Draw();


  TFile *rf_nom = new TFile(Form("acc_single_%s_%dS.root",fCollId.Data(),state),"recreate"); 
  TFile *rf_sys = new TFile(Form("sys_single_%s_%dS.root",fCollId.Data(),state),"recreate"); 

  rf_nom->cd();
  hptGen->Write();
  hptGenAcc->Write();
  hptAcc->Write();
  hrapGen->Write();
  hrapGenAcc->Write();
  hrapAcc->Write();
  hintGen->Write();
  hintGenAcc->Write();
  hintAcc->Write();
  rf_nom->Close(); 

  rf_sys->cd();
  hptAccSys[0]->Write();
  hrapAccSys[0]->Write();
  hintAccSys[0]->Write();
  rf_sys->Close();
  
}

//Get Yield
valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh, float dphiEp2Low,  float dphiEp2High) {
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, glbMuPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TFile* inf = new TFile(Form("../fitResults/Constrain/fitresults_upsilon_fixParm1_seed2_DoubleCB_%s.root",kineLabel.Data())); //Free Parameter
  cout<<kineLabel.Data()<<endl;
  RooWorkspace* ws = (RooWorkspace*)inf->Get("workspace");
  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  valErr ret;
  ret.val = fitResults->GetBinContent(state);
  ret.err = fitResults->GetBinError(state);
  //cout << kineLabel << ": " << " & " << ret.val << " $\pm$ " << ret.err << " & " <<ws->var("nBkg")->getVal() << " $\pm$ "<< ws->var("nBkg")->getError() << "\\\\" << endl;
  return ret;
}

void getMaxTH1D_eight ( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7, TH1D* h8) {
  if ( h0->GetNbinsX() != h1->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h2->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h3->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h4->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h5->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h6->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h7->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h8->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;

  for ( int i=1 ;  i<=h0->GetNbinsX(); i++) {
    float x1 =  fabs(h1->GetBinContent(i));
    float x2 =  fabs(h2->GetBinContent(i));
    float x3 =  fabs(h3->GetBinContent(i));
    float x4 =  fabs(h4->GetBinContent(i));
    float x5 =  fabs(h5->GetBinContent(i));
    float x6 =  fabs(h6->GetBinContent(i));
    float x7 =  fabs(h7->GetBinContent(i));
    float x8 =  fabs(h8->GetBinContent(i));
    float x = 0;
    if ( x1 > x2 ) x = x1;
    else x = x2;
    if ( x <= x3) x = x3;
    if ( x <= x4) x = x4;
    if ( x <= x5) x = x5;
    if ( x <= x6) x = x6;
    if ( x <= x7) x = x7;
    if ( x <= x8) x = x8;
    h0->SetBinContent(i, x );
  }
}

void getMaxTH1D_four ( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4) {
  if ( h0->GetNbinsX() != h1->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h2->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h3->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h4->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;

  for ( int i=1 ;  i<=h0->GetNbinsX(); i++) {
    float x1 =  fabs(h1->GetBinContent(i));
    float x2 =  fabs(h2->GetBinContent(i));
    float x3 =  fabs(h3->GetBinContent(i));
    float x4 =  fabs(h4->GetBinContent(i));
    float x = 0;
    if ( x1 > x2 ) x = x1;
    else x = x2;
    if ( x <= x3) x = x3;
    if ( x <= x4) x = x4;
    h0->SetBinContent(i, x );
  }
}

