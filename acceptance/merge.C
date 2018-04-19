#include <iostream>
#include "TH1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"

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
  



