#include <iostream>
#include "TH1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"
#include "../commonUtility.h"
#include "../SONGKYO.h"

using namespace std;

void merge()
{
  TFile *fpp1s = new TFile("acc_single_PP_1S.root","read");
  TFile *fpp2s = new TFile("acc_single_PP_2S.root","read");
  TFile *fpp3s = new TFile("acc_single_PP_3S.root","read");
  TFile *faa1s = new TFile("acc_single_AA_1S.root","read");
  TFile *faa2s = new TFile("acc_single_AA_2S.root","read");
  TFile *faa3s = new TFile("acc_single_AA_3S.root","read");

  TH1D* hptAccAA1S = (TH1D*) faa1s->Get("hptAccAA1S");
  TH1D* hptAccPP1S = (TH1D*) fpp1s->Get("hptAccPP1S");
  TH1D* hrapAccAA1S = (TH1D*) faa1s->Get("hrapAccAA1S");
  TH1D* hrapAccPP1S = (TH1D*) fpp1s->Get("hrapAccPP1S");
  TH1D* hintAccAA1S = (TH1D*) faa1s->Get("hIntAccAA1S");
  TH1D* hintAccPP1S = (TH1D*) fpp1s->Get("hIntAccPP1S");

  TH1D* hptAccAA2S = (TH1D*) faa2s->Get("hptAccAA2S");
  TH1D* hptAccPP2S = (TH1D*) fpp2s->Get("hptAccPP2S");
  TH1D* hrapAccAA2S = (TH1D*) faa2s->Get("hrapAccAA2S");
  TH1D* hrapAccPP2S = (TH1D*) fpp2s->Get("hrapAccPP2S");
  TH1D* hintAccAA2S = (TH1D*) faa2s->Get("hIntAccAA2S");
  TH1D* hintAccPP2S = (TH1D*) fpp2s->Get("hIntAccPP2S");

  TH1D* hptAccAA3S = (TH1D*) faa3s->Get("hptAccAA3S");
  TH1D* hptAccPP3S = (TH1D*) fpp3s->Get("hptAccPP3S");
  TH1D* hrapAccAA3S = (TH1D*) faa3s->Get("hrapAccAA3S");
  TH1D* hrapAccPP3S = (TH1D*) fpp3s->Get("hrapAccPP3S");
  TH1D* hintAccAA3S = (TH1D*) faa3s->Get("hIntAccAA3S");
  TH1D* hintAccPP3S = (TH1D*) fpp3s->Get("hIntAccPP3S");

  TFile *rfacc1s = new TFile("acceptance_wgt_norm_1S.root","recreate");
  rfacc1s->cd();
  hptAccAA1S->Write();
  hptAccPP1S->Write();
  hrapAccAA1S->Write();
  hrapAccPP1S->Write();
  hintAccAA1S->Write();
  hintAccPP1S->Write();

  TFile *rfacc2s = new TFile("acceptance_wgt_norm_2S.root","recreate");
  rfacc2s->cd();
  hptAccAA2S->Write();
  hptAccPP2S->Write();
  hrapAccAA2S->Write();
  hrapAccPP2S->Write();
  hintAccAA2S->Write();
  hintAccPP2S->Write();

  TFile *rfacc3s = new TFile("acceptance_wgt_norm_3S.root","recreate");
  rfacc3s->cd();
  hptAccAA3S->Write();
  hptAccPP3S->Write();
  hrapAccAA3S->Write();
  hrapAccPP3S->Write();
  hintAccAA3S->Write();
  hintAccPP3S->Write();

  handsomeTH1(hptAccPP1S,1);
  handsomeTH1(hptAccAA1S,2);
  handsomeTH1(hptAccPP2S,1);
  handsomeTH1(hptAccAA2S,2);
  handsomeTH1(hptAccPP3S,1);
  handsomeTH1(hptAccAA3S,2);

  handsomeTH1(hrapAccPP1S,1);
  handsomeTH1(hrapAccAA1S,2);
  handsomeTH1(hrapAccPP2S,1);
  handsomeTH1(hrapAccAA2S,2);
  handsomeTH1(hrapAccPP3S,1);
  handsomeTH1(hrapAccAA3S,2);

  TLegend *leg1 = new TLegend(0.5,0.6,0.8,0.8);
  SetLegendStyle(leg1);
  leg1->AddEntry(hptAccPP1S,"PP #varUpsilon(1S)","l");
  leg1->AddEntry(hptAccAA1S,"AA #varUpsilon(1S)","l");
  TLegend *leg2 = new TLegend(0.5,0.6,0.8,0.8);
  SetLegendStyle(leg2);
  leg2->AddEntry(hptAccPP2S,"PP #varUpsilon(2S)","l");
  leg2->AddEntry(hptAccAA2S,"AA #varUpsilon(2S)","l");
  TLegend *leg3 = new TLegend(0.5,0.6,0.8,0.8);
  SetLegendStyle(leg3);
  leg3->AddEntry(hptAccPP3S,"PP #varUpsilon(3S)","l");
  leg3->AddEntry(hptAccAA3S,"AA #varUpsilon(3S)","l");

  TLegend *lg1 = new TLegend(0.5,0.6,0.8,0.8);
  SetLegendStyle(lg1);
  lg1->AddEntry(hrapAccPP1S,"PP #varUpsilon(1S)","l");
  lg1->AddEntry(hrapAccAA1S,"AA #varUpsilon(1S)","l");
  TLegend *lg2 = new TLegend(0.5,0.6,0.8,0.8);
  SetLegendStyle(lg2);
  lg2->AddEntry(hrapAccPP2S,"PP #varUpsilon(2S)","l");
  lg2->AddEntry(hrapAccAA2S,"AA #varUpsilon(2S)","l");
  TLegend *lg3 = new TLegend(0.5,0.6,0.8,0.8);
  SetLegendStyle(lg3);
  lg3->AddEntry(hrapAccPP3S,"PP #varUpsilon(3S)","l");
  lg3->AddEntry(hrapAccAA3S,"AA #varUpsilon(3S)","l");

  hptAccPP1S->GetYaxis()->SetTitle("Acceptance");
  hptAccPP2S->GetYaxis()->SetTitle("Acceptance");
  hptAccPP3S->GetYaxis()->SetTitle("Acceptance");
  hrapAccPP1S->GetYaxis()->SetTitle("Acceptance");
  hrapAccPP2S->GetYaxis()->SetTitle("Acceptance");
  hrapAccPP3S->GetYaxis()->SetTitle("Acceptance");

  TCanvas *c1 = new TCanvas("c1","",1200,500);
  c1->Divide(3,1);
  c1->cd(1);
  hptAccPP1S->Draw();
  hptAccAA1S->Draw("same");
  leg1->Draw("same");
  c1->cd(2);
  hptAccPP2S->Draw();
  hptAccAA2S->Draw("same");
  leg2->Draw("same");
  c1->cd(3);
  hptAccPP3S->Draw();
  hptAccAA3S->Draw("same");
  leg3->Draw("same");
  c1->SaveAs("plots/c_acc_pt_comp.pdf");  

  TCanvas *c2 = new TCanvas("c2","",1200,500);
  c2->Divide(3,1);
  c2->cd(1);
  hrapAccPP1S->Draw();
  hrapAccAA1S->Draw("same");
  lg1->Draw("same");
  c2->cd(2);
  hrapAccPP2S->Draw();
  hrapAccAA2S->Draw("same");
  lg2->Draw("same");
  c2->cd(3);
  hrapAccPP3S->Draw();
  hrapAccAA3S->Draw("same");
  lg3->Draw("same");
  c2->SaveAs("plots/c_acc_rap_comp.pdf");  



  TFile *f_syspp1s = new TFile("sys_single_PP_1S.root","read");
  TFile *f_syspp2s = new TFile("sys_single_PP_2S.root","read");
  TFile *f_syspp3s = new TFile("sys_single_PP_3S.root","read");
  TFile *f_sysaa1s = new TFile("sys_single_AA_1S.root","read");
  TFile *f_sysaa2s = new TFile("sys_single_AA_2S.root","read");
  TFile *f_sysaa3s = new TFile("sys_single_AA_3S.root","read");

  TH1D* hptSysAA1S  = (TH1D*) f_sysaa1s->Get("hptSysAA");
  TH1D* hptSysPP1S  = (TH1D*) f_syspp1s->Get("hptSysPP");
  TH1D* hrapSysAA1S = (TH1D*) f_sysaa1s->Get("hrapSysAA");
  TH1D* hrapSysPP1S = (TH1D*) f_syspp1s->Get("hrapSysPP");
  TH1D* hintSysAA1S = (TH1D*) f_sysaa1s->Get("hcentSysAA_int");
  TH1D* hintSysPP1S = (TH1D*) f_syspp1s->Get("hcentSysPP");

  TH1D* hptSysAA2S  = (TH1D*) f_sysaa2s->Get("hptSysAA");
  TH1D* hptSysPP2S  = (TH1D*) f_syspp2s->Get("hptSysPP");
  TH1D* hrapSysAA2S = (TH1D*) f_sysaa2s->Get("hrapSysAA");
  TH1D* hrapSysPP2S = (TH1D*) f_syspp2s->Get("hrapSysPP");
  TH1D* hintSysAA2S = (TH1D*) f_sysaa2s->Get("hcentSysAA_int");
  TH1D* hintSysPP2S = (TH1D*) f_syspp2s->Get("hcentSysPP");

  TH1D* hptSysAA3S  = (TH1D*) f_sysaa3s->Get("hptSysAA");
  TH1D* hptSysPP3S  = (TH1D*) f_syspp3s->Get("hptSysPP");
  TH1D* hrapSysAA3S = (TH1D*) f_sysaa3s->Get("hrapSysAA");
  TH1D* hrapSysPP3S = (TH1D*) f_syspp3s->Get("hrapSysPP");
  TH1D* hintSysAA3S = (TH1D*) f_sysaa3s->Get("hcentSysAA_int");
  TH1D* hintSysPP3S = (TH1D*) f_syspp3s->Get("hcentSysPP");


  TFile *rfsys1s = new TFile("sys_acceptance_ups1S_1804.root","recreate");
  rfsys1s->cd();
  hptSysAA1S->Write();
  hptSysPP1S->Write();
  hrapSysAA1S->Write();
  hrapSysPP1S->Write();
  hintSysAA1S->Write();
  hintSysPP1S->Write();

  TFile *rfsys2s = new TFile("sys_acceptance_ups2S_1804.root","recreate");
  rfsys2s->cd();
  hptSysAA2S->Write();
  hptSysPP2S->Write();
  hrapSysAA2S->Write();
  hrapSysPP2S->Write();
  hintSysAA2S->Write();
  hintSysPP2S->Write();

  TFile *rfsys3s = new TFile("sys_acceptance_ups3S_1804.root","recreate");
  rfsys3s->cd();
  hptSysAA3S->Write();
  hptSysPP3S->Write();
  hrapSysAA3S->Write();
  hrapSysPP3S->Write();
  hintSysAA3S->Write();
  hintSysPP3S->Write();
}
  



