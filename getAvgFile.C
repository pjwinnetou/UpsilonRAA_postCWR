#include <iostream>
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TH1.h"
using namespace std;

double *getAvg_n(double avg_alpha = 0, double avg_n =0, double avg_sigma=0, double avg_f=0, double avg_x=0, int states=1, TString szAA = "PP");
double *getAvg_ax(double avg_alpha = 0, double avg_n =0, double avg_sigma=0, double avg_f=0, double avg_x=0, int states=1, TString szAA = "PP");

void getAvgFile(TString szAA = "PP", int states =1) 
{
  TFile *rf = new TFile(Form("AvgSigPar_%s.root",szAA.Data()),"recreate");
  TH1D* hAvg_n = new TH1D("hAvgn",";;;",10,0,10);
  TH1D* hAvg_alphax = new TH1D("hAvgalphax",";;;",10,0,10);

  double avg_alpha=0; double avg_n=0; double avg_sigma=0; double avg_f=0; double avg_x=0; 
  double avg_alpha_=0; double avg_n_=0; double avg_sigma_=0; double avg_f_=0; double avg_x_=0; 
  double* avgParm_n = getAvg_n(avg_alpha, avg_n, avg_sigma, avg_f, avg_x, states,szAA);
  double* avgParm_alphax = getAvg_ax(avg_alpha_, avg_n_, avg_sigma_, avg_f_, avg_x_, states,szAA);
  
  for(int i=0;i<10;i++){
    hAvg_n->SetBinContent(i+1,avgParm_n[i]);
    hAvg_alphax->SetBinContent(i+1,avgParm_alphax[i]);
  }
  rf->cd();
  hAvg_n->Write();
  hAvg_alphax->Write();

  

}

double *getAvg_n(double avg_alpha, double avg_n , double avg_sigma, double avg_f, double avg_x, int states, TString szAA)
{
  int nFit = 2; 
  char *Fit_loc[2] = {"../../upsilonRAA5TeV/fitResults/Final_NomResult_170124/", "MCFixedFit/"};
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

  TFile *fileIn_[nBin_];
  RooWorkspace* ws_[nBin_];
  double alpha_[nBin_];
  double n1s_[nBin_];
  double sigma_[nBin_];
  double f1s_[nBin_];
  double x1s_[nBin_];

  double alphaErr=0;
  double nErr=0;
  double fErr=0;
  double xErr=0;
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
    else if (szAA == "AA" ) fileIn[ib]= new TFile(Form("%sfitresults_upsilon_fixParm1_seed2_DoubleCB_%s_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",fileLoc[0].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));
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
  cout << "avg_x : " << avg_x << endl;

  for (int ib =0; ib < nBin; ib ++ ) {
    alphaErr += (alpha[ib]-avg_alpha)*(alpha[ib]-avg_alpha);
    nErr += (n1s[ib]-avg_n)*(n1s[ib]-avg_n);
    sigmaErr += (sigma[ib]-avg_sigma)*(sigma[ib]-avg_sigma);
    fErr += (f1s[ib]-avg_f)*(f1s[ib]-avg_f);
    xErr += (x1s[ib]-avg_x)*(x1s[ib]-avg_x);
  }
  alphaErr = TMath::Sqrt( (alphaErr) / (nBin*(nBin-1)) ) ;
  nErr = TMath::Sqrt( (nErr) / (nBin*(nBin-1)) ) ;
  sigmaErr = TMath::Sqrt( (sigmaErr) / (nBin*(nBin-1)) ) ;
  fErr = TMath::Sqrt( (fErr) / (nBin*(nBin-1)) ) ;
  xErr = TMath::Sqrt( (xErr) / (nBin*(nBin-1)) ) ;
    

  double* avgParm = new double[10];
  avgParm[0] = avg_alpha; avgParm[1] =  avg_n; avgParm[2] = avg_sigma; avgParm[3] = avg_f; avgParm[4] = avg_x;
  avgParm[5] = alphaErr; avgParm[6] = nErr; avgParm[7] = sigmaErr; avgParm[8] = fErr; avgParm[9] = xErr;
  return avgParm;
}

double *getAvg_ax(double avg_alpha, double avg_n , double avg_sigma, double avg_f, double avg_x, int states, TString szAA)
{
  int nFit = 2; 
  char *Fit_loc[2] = {"../UpsilonRAA_postCWR_fixed_n/fitResults/FixedFit_n_Avg/", "MCFixedFit/"};
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

  TFile *fileIn_[nBin_];
  RooWorkspace* ws_[nBin_];
  double alpha_[nBin_];
  double n1s_[nBin_];
  double sigma_[nBin_];
  double f1s_[nBin_];
  double x1s_[nBin_];

  double alphaErr=0;
  double nErr=0;
  double fErr=0;
  double xErr=0;
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
    if (szAA == "PP" ) fileIn[ib]= new TFile(Form("%sfitresults_upsilon_fixParm1_seed2_DoubleCB_%s_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0.root",fileLoc[0].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));
    else if (szAA == "AA" ) fileIn[ib]= new TFile(Form("%sfitresults_upsilon_fixParm1_seed2_DoubleCB_%s_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",fileLoc[0].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));
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
  cout << "avg_x : " << avg_x << endl;

  for (int ib =0; ib < nBin; ib ++ ) {
    alphaErr += (alpha[ib]-avg_alpha)*(alpha[ib]-avg_alpha);
    nErr += (n1s[ib]-avg_n)*(n1s[ib]-avg_n);
    sigmaErr += (sigma[ib]-avg_sigma)*(sigma[ib]-avg_sigma);
    fErr += (f1s[ib]-avg_f)*(f1s[ib]-avg_f);
    xErr += (x1s[ib]-avg_x)*(x1s[ib]-avg_x);
  }
  alphaErr = TMath::Sqrt( (alphaErr) / (nBin*(nBin-1)) ) ;
  nErr = TMath::Sqrt( (nErr) / (nBin*(nBin-1)) ) ;
  sigmaErr = TMath::Sqrt( (sigmaErr) / (nBin*(nBin-1)) ) ;
  fErr = TMath::Sqrt( (fErr) / (nBin*(nBin-1)) ) ;
  xErr = TMath::Sqrt( (xErr) / (nBin*(nBin-1)) ) ;
    

  double* avgParm = new double[10];
  avgParm[0] = avg_alpha; avgParm[1] =  avg_n; avgParm[2] = avg_sigma; avgParm[3] = avg_f; avgParm[4] = avg_x;
  avgParm[5] = alphaErr; avgParm[6] = nErr; avgParm[7] = sigmaErr; avgParm[8] = fErr; avgParm[9] = xErr;
  return avgParm;
}
