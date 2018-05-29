#include "SONGKYO.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "../cutsAndBin.h"

void compare_15001_RAA_cent_newglobal_asym_postCWRConstrain(int istate=2) //1 or 2 (1S or 2S)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 100; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nfile = 2; // 0: 15001, 1: ours
  double xmax = 420;
  double xmin_int = 0.5;
  double xmax_int = 1.5;
  double boxw = 6.5; // for syst. box (vs cent only)
  double boxw_int = 0.09;
//  double relsys = 0.1;
  
  //// 15001 values
  const int cn_1s =  8;
  double cpx_1s[cn_1s] =  {8.8, 42.0, 86.3, 130, 187, 261, 329, 381};
  double cpy_1s[cn_1s] =  {1.271, 0.815, 0.605, 0.546, 0.476, 0.442, 0.408, 0.330};
  double cex_1s[cn_1s] =  {0., 0., 0., 0., 0., 0., 0., 0.};
  double cey_1s[cn_1s] =  {0.229, 0.083, 0.065, 0.047, 0.035, 0.028, 0.033, 0.029};
  double cexsys_1s[cn_1s] =  {boxw, boxw, boxw, boxw, boxw, boxw, boxw, boxw};
  double ceysys_1s_1[cn_1s] =  {0.32, 0.134, 0.119, 0.084, 0.056, 0.055, 0.057, 0.052};
  double ceysys_1s_2[cn_1s] =  {0.077, 0.051, 0.038, 0.034, 0.030, 0.028, 0.025, 0.021};
  double ceysys_1s[cn_1s] = {0.401, 0.137, 0.144, 0.120, 0.059, 0.093, 0.073, 0.070};
 // for (int it=0; it < cn_1s ; it++) {
 //   ceysys_1s[it] = TMath::Sqrt(ceysys_1s_1[it]*ceysys_1s_1[it]+ceysys_1s_2[it]*ceysys_1s_2[it]);
 // }

  const int cn_2s =  4;
  double cpx_2s[cn_2s] =  {22.1, 108, 224, 355};
  double cpy_2s[cn_2s] =  {0.235, 0.294, 0.092, 0.076};
  double cex_2s[cn_2s] =  {0., 0., 0., 0.};
  double cey_2s[cn_2s] =  {0.155, 0.079, 0.043, 0.045};
  double cexsys_2s[cn_2s] =  {boxw, boxw, boxw, boxw};
  double ceysys_2s_1[cn_2s] =  {0.106, 0.081, 0.052, 0.035};
  double ceysys_2s_2[cn_2s] =  {0.015, 0.019, 0.006, 0.005};
  double ceysys_2s[cn_2s] = {0.128, 0.112, 0.065, 0.038};
//  for (int it=0; it < cn_2s ; it++) {
//    ceysys_2s[it] = TMath::Sqrt(ceysys_2s_1[it]*ceysys_2s_1[it]+ceysys_2s_2[it]*ceysys_2s_2[it]);
//  }

  const int cn_int = 1;
  double cpx_1s_int[cn_int] = {1};
  double cpy_1s_int[cn_int] = {0.453};
  double cex_1s_int[cn_int] = {0.};
  double cey_1s_int[cn_int] = {0.014};
  double cexsys_1s_int[cn_int] = {boxw_int};
  double ceysys_1s_int[cn_int] = {0.046};//0.029
  
  double cpx_2s_int[cn_int] = {1};
  double cpy_2s_int[cn_int] = {0.119};
  double cex_2s_int[cn_int] = {0.};
  double cey_2s_int[cn_int] = {0.028};
  double cexsys_2s_int[cn_int] = {boxw_int};
  double ceysys_2s_int[cn_int] = {0.015};//0.008

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
	TGraphErrors* gRAA[nfile]; // vs centrality
  TGraphAsymmErrors* gRAA_sys[nfile];
	TGraphErrors* gRAA_int[nfile]; // centrality-integrated
  TGraphAsymmErrors* gRAA_int_sys[nfile];
  //// 1) 15001
  if (istate==1) {
    gRAA[0] = new TGraphErrors(cn_1s, cpx_1s, cpy_1s, cex_1s, cey_1s); 
    gRAA_sys[0] = new TGraphAsymmErrors(cn_1s, cpx_1s, cpy_1s, cexsys_1s, cexsys_1s, ceysys_1s, ceysys_1s); 
    gRAA_int[0] = new TGraphErrors(cn_int, cpx_1s_int, cpy_1s_int, cex_1s_int, cey_1s_int); 
    gRAA_int_sys[0] = new TGraphAsymmErrors(cn_int, cpx_1s_int, cpy_1s_int, cexsys_1s_int, cexsys_1s_int, ceysys_1s_int, ceysys_1s_int); 
  }   
  else {
    gRAA[0] = new TGraphErrors(cn_2s, cpx_2s, cpy_2s, cex_2s, cey_2s); 
    gRAA_sys[0] = new TGraphAsymmErrors(cn_2s, cpx_2s, cpy_2s, cexsys_2s, cexsys_2s, ceysys_2s, ceysys_2s); 
    gRAA_int[0] = new TGraphErrors(cn_int, cpx_2s_int, cpy_2s_int, cex_2s_int, cey_2s_int); 
    gRAA_int_sys[0] = new TGraphAsymmErrors(cn_int, cpx_2s_int, cpy_2s_int, cexsys_2s_int, cexsys_2s_int, ceysys_2s_int, ceysys_2s_int); 
  } 
  //// 2) ours
  TFile* fIn = new TFile(Form("Ups_%d_RAA.root",istate),"READ");
  gRAA[1]=(TGraphErrors*)fIn->Get("gRAA_cent");
  //gRAA_sys[1]=(TGraphErrors*)fIn->Get("gRAA_cent");
  gRAA_sys[1] = new TGraphAsymmErrors();
  gRAA_int[1]=(TGraphErrors*)fIn->Get("gRAA_int");
  gRAA_int_sys[1] = new TGraphAsymmErrors();
  //gRAA_int_sys[1]=(TGraphErrors*)fIn->Get("gRAA_int");

  //// read input file : syst.
  TFile* fInSys_Hi;
  TFile* fInSys_Lo;
  TH1D* hSys_Hi;
  TH1D* hSys_Lo;
  TH1D* hSys_int_Hi;
  TH1D* hSys_int_Lo;
  int npoint_Hi;
  int npoint_Lo;
  int npoint_int_Hi;
  int npoint_int_Lo;
    fInSys_Hi = new TFile(Form("../Systematic/postCWR_mergedSys_constrain_ups%ds_asymHi.root",istate),"READ");
    fInSys_Lo = new TFile(Form("../Systematic/postCWR_mergedSys_constrain_ups%ds_asymLo.root",istate),"READ");
    hSys_Hi=(TH1D*)fInSys_Hi->Get("hcentRAA_merged");
    hSys_Lo=(TH1D*)fInSys_Lo->Get("hcentRAA_merged");
    npoint_Hi = hSys_Hi->GetSize()-2;
    npoint_Lo = hSys_Lo->GetSize()-2;
    hSys_int_Hi=(TH1D*)fInSys_Hi->Get("hintRAA_merged");
    hSys_int_Lo=(TH1D*)fInSys_Lo->Get("hintRAA_merged");
    npoint_int_Hi = hSys_int_Hi->GetSize()-2;
    npoint_int_Lo = hSys_int_Lo->GetSize()-2;
     

  //// set bin width and calculate systematic uncertainties 
  double pxtmp, pytmp, extmp, eytmp;
  double relsys_Lo, relsys_Hi;
  //// --- vs centrality
  if (npoint_Lo != gRAA[1]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
  for (int ipt=0; ipt< npoint_Lo; ipt++) {
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys_Lo=0; relsys_Hi=0;
    gRAA[1]->GetPoint(ipt, pxtmp, pytmp);
    extmp=gRAA[1]->GetErrorX(ipt);
    eytmp=gRAA[1]->GetErrorY(ipt);
    relsys_Hi=hSys_Hi->GetBinContent(npoint_Hi-ipt);
    relsys_Lo=hSys_Lo->GetBinContent(npoint_Lo-ipt);
    // 1) remove ex from gRAA
    gRAA[1]->SetPointError(ipt, 0, eytmp);
    // 2) set ey for gRAA_sys
    //if (istate==1) gRAA_sys[1]->SetPointError(ipt, exsys_1s[ipt], pytmp*relsys);
    //else gRAA_sys[1]->SetPointError(ipt, exsys_2s[ipt], pytmp*relsys);
    gRAA_sys[1]->SetPoint(ipt,pxtmp,pytmp);
    gRAA_sys[1]->SetPointError(ipt, boxw, boxw, pytmp*relsys_Lo, pytmp*relsys_Hi);
  }
  //// --- centrality-integrated
  if (npoint_int_Lo != gRAA_int[1]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
  for (int ipt=0; ipt< npoint_int_Lo; ipt++) {
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys_Hi=0;relsys_Lo=0;
    gRAA_int[1]->GetPoint(ipt, pxtmp, pytmp);
    extmp=gRAA_int[1]->GetErrorX(ipt);
    eytmp=gRAA_int[1]->GetErrorY(ipt);
    relsys_Hi = hSys_int_Hi->GetBinContent(ipt+1);
    relsys_Hi = TMath::Sqrt(relsys_Hi*relsys_Hi+lumi_unc_pp*lumi_unc_pp + nMB_unc*nMB_unc);
    relsys_Lo = hSys_int_Lo->GetBinContent(ipt+1);
    relsys_Lo = TMath::Sqrt(relsys_Lo*relsys_Lo+lumi_unc_pp*lumi_unc_pp + nMB_unc*nMB_unc);
    // 1) remove ex from gRAA
    gRAA_int[1]->SetPointError(ipt, 0, eytmp);
    // 2) set ey for gRAA_int_sys
    //if (istate==1) gRAA_int_sys[1]->SetPointError(ipt, exsys_1s[ipt], pytmp*relsys);
    //else gRAA_int_sys[1]->SetPointError(ipt, exsys_2s[ipt], pytmp*relsys);
    gRAA_int_sys[1]->SetPoint(ipt,pxtmp,pytmp);
    gRAA_int_sys[1]->SetPointError(ipt, boxw_int, boxw_int, pytmp*relsys_Lo,  pytmp*relsys_Hi);
  }
  
  ////////////////////////////////////////////////////////////////
 
  //// graph style 
  if(istate==1){
    SetGraphStyleOpen(gRAA[0], 4, 4, 0);
    SetGraphStyleSys(gRAA_sys[0], 4); 
    SetGraphStyleOpen(gRAA_int[0], 4, 4, 0); 
    SetGraphStyleSys(gRAA_int_sys[0], 4); 
    SetGraphStyle(gRAA[1], 0, 0); 
    SetGraphStyleSys(gRAA_sys[1], 0); 
    SetGraphStyle(gRAA_int[1], 0, 0); 
    SetGraphStyleSys(gRAA_int_sys[1], 0); 
  }
  else if(istate==2){
    SetGraphStyleOpen(gRAA[0], 4, 4, 1);
    SetGraphStyleSys(gRAA_sys[0], 4); 
    SetGraphStyleOpen(gRAA_int[0], 4, 4, 1); 
    SetGraphStyleSys(gRAA_int_sys[0], 4); 
    SetGraphStyle(gRAA[1], 1, 1); 
    SetGraphStyleSys(gRAA_sys[1], 1); 
    SetGraphStyle(gRAA_int[1], 1, 1); 
    SetGraphStyleSys(gRAA_int_sys[1], 1); 
  }

  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.038);
  
  //// legend
  TLegend *leg= new TLegend(0.57, 0.5045, 0.92, 0.6665);
  SetLegendStyle(leg);
  leg -> SetHeader(Form("#varUpsilon(%dS)",istate));
  leg -> AddEntry(gRAA[1],"#sqrt{s_{NN}} = 5.02 TeV","lp");
  leg -> AddEntry(gRAA[0],"#sqrt{s_{NN}} = 2.76 TeV","lp");
  TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
  header->SetTextFont(62);

  //// axis et. al
  gRAA_sys[0]->GetXaxis()->SetTitle("N_{part}");
  gRAA_sys[0]->GetXaxis()->CenterTitle();
  gRAA_sys[0]->GetYaxis()->SetTitle("R_{AA}");
  gRAA_sys[0]->GetYaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetLimits(0.,xmax);
  gRAA_sys[0]->SetMinimum(0.0);
//  gRAA_sys[0]->SetMaximum(1.3);
  gRAA_sys[0]->SetMaximum(1.7);
  //// for cent
  gRAA_sys[0]->GetXaxis()->SetTitleSize(0.06*1.0);
  gRAA_sys[0]->GetYaxis()->SetTitleSize(0.06*1.0);
  gRAA_sys[0]->GetXaxis()->SetLabelSize(0.05*1.0);
  gRAA_sys[0]->GetYaxis()->SetLabelSize(0.05*1.0);

  //// draw  
  double xlonger = 120; 
  TCanvas* c1 = new TCanvas("c1","c1",600+xlonger,600);
  TPad* pad_diff = new TPad("pad_diff", "",0, 0, 600/(600.+xlonger), 1.0); // vs centrality
  pad_diff->SetRightMargin(0);
  TPad* pad_int = new TPad("pad_int", "",600/(600.+xlonger), 0, 1.0, 1.0); // centrality-integrated
  pad_int->SetLeftMargin(0);
  pad_int->SetRightMargin(0.032*600/xlonger);

  //// --- 1st pad!!!
  c1->cd();
  pad_diff->Draw(); 
  pad_diff->cd(); 
  for (int is=0; is<nfile; is++){
    if ( is==0) gRAA_sys[is]->Draw("A5");
    else gRAA_sys[is]->Draw("5");
    gRAA[is]->Draw("P");
	}
  dashedLine(0.,1.,xmax,1.,1,1);
  leg->Draw();

  //// drwa text
  double sz_init = 0.886; double sz_step = 0.0558;
//  globtex->DrawLatex(0.22+0.04, sz_init, "p_{T}^{#mu} > 4 GeV/c");
  globtex->DrawLatex(0.22+0.04, sz_init, "p_{T} < 30 GeV");
//  globtex->DrawLatex(0.46+0.04, sz_init+0.002, "|#eta|^{#mu} < 2.4");
  globtex->DrawLatex(0.22+0.04, sz_init-sz_step, "|y| < 2.4");
  
  //Global Unc.
  TFile* fppInt = new TFile("../fitResults/Constrain/fitresults_upsilon_fixParm1_seed2_DoubleCB_PP_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0.root");
  TH1D* hppInt = (TH1D*) fppInt -> Get("fitResults");
  double relsys_ppInt1S = hppInt->GetBinError(1)/hppInt->GetBinContent(1);
  double relsys_ppInt2S = hppInt->GetBinError(2)/hppInt->GetBinContent(2);
  TH1D* hSys_glb;
  double sys_global_pp;
  double sys_global_val;
  double accept_sys;
  hSys_glb = (TH1D*) fInSys_Lo->Get("hintPP_merged");
  TFile* f_acc; 
  TH1D* hSys_glb_acc;
  if(istate==1) {f_acc = new TFile("../acceptance/sys_acceptance_ups1S_1804.root"); hSys_glb_acc = (TH1D*) f_acc->Get("hcentSysPP"); accept_sys = hSys_glb_acc->GetBinContent(1);}
  else if(istate==2) {f_acc = new TFile("../acceptance/sys_acceptance_ups2S_1804.root"); hSys_glb_acc = (TH1D*) f_acc->Get("hcentSysPP"); accept_sys = hSys_glb_acc->GetBinContent(1);}
  sys_global_pp = TMath::Sqrt(hSys_glb->GetBinContent(1)*hSys_glb->GetBinContent(1)+accept_sys*accept_sys);
  
  sys_global_val = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp + nMB_unc*nMB_unc);
  double sys_global_y = TMath::Sqrt(sys_global_val*sys_global_val + sys_global_pp*sys_global_pp); 
  double sys_global_y_15001;
  if(istate==1) sys_global_y_15001 = TMath::Sqrt(0.032*0.032+0.063*0.063); 
  else if(istate==2) sys_global_y_15001 = TMath::Sqrt(0.032*0.032+0.069*0.069); 
  double sys_global_x = 15;

  if(istate==1) sys_global_y = TMath::Sqrt(sys_global_y*sys_global_y+relsys_ppInt1S*relsys_ppInt1S);
  else if(istate==2) sys_global_y = TMath::Sqrt(sys_global_y*sys_global_y+relsys_ppInt2S*relsys_ppInt2S);

  TBox *globalUncBox = new TBox(xmax-sys_global_x*2,1-sys_global_y,xmax-sys_global_x,1+sys_global_y);
  if(istate==1){
    globalUncBox -> SetLineColor(kRed-2);
    globalUncBox -> SetFillColorAlpha(kPink-6,0.6);
  }
  else if(istate==2){
    globalUncBox -> SetLineColor(kBlue-3);
    globalUncBox -> SetFillColorAlpha(kBlue-3,0.6);
  }
  globalUncBox -> SetLineWidth(1);
  globalUncBox -> Draw("l same");

  TBox *ppRefUncBox1S = new TBox(xmax-sys_global_x,1-sys_global_y_15001,xmax,1+sys_global_y_15001);
  ppRefUncBox1S -> SetLineColor(kBlack);
  ppRefUncBox1S -> SetFillColorAlpha(kGray+3,0.6);
  ppRefUncBox1S -> SetLineWidth(1);
  ppRefUncBox1S -> Draw("same");

  CMS_lumi( pad_diff, iPeriod, iPos );
  
  //// --- 2nd pad!!!
  c1->cd();
  pad_int->Draw();
  pad_int->cd();

  //// for int
  gRAA_int_sys[0]->GetXaxis()->SetLimits(xmin_int,xmax_int);
  gRAA_int_sys[0]->SetMinimum(0.0);
//  gRAA_int_sys[0]->SetMaximum(1.3);
  gRAA_int_sys[0]->SetMaximum(1.7);
  gRAA_int_sys[0]->GetXaxis()->SetNdivisions(101);
  gRAA_int_sys[0]->GetXaxis()->SetLabelSize(0);
  gRAA_int_sys[0]->GetYaxis()->SetTickLength(0.03*600/xlonger);
  gRAA_int_sys[0]->GetYaxis()->SetLabelSize(0);
  
  for (int is=0; is<nfile; is++){
    if ( is==0) gRAA_int_sys[is]->Draw("A5");
    else gRAA_int_sys[is]->Draw("5");
    gRAA_int[is]->Draw("P");
	}
  dashedLine(0.,1.,xmax,1.,1,1);

  //// draw text 
  double sz_allign = 0.1897;
  globtex->SetTextAlign(22); //center-center
  globtex->SetTextSize(0.038*600./xlonger);
  globtex->DrawLatex(0.5*(1-0.032*600/xlonger), sz_init-sz_step-sz_allign, "Cent.");
  globtex->DrawLatex(0.5*(1-0.032*600/xlonger), sz_init-sz_step*2-sz_allign, "0-100%"); 
  
  c1->SaveAs(Form("plots/%dS_comp15001_RAA_vs_cent_newglobal_asym_Constrain.pdf",istate));
  c1->SaveAs(Form("plots/%dS_comp15001_RAA_vs_cent_newglobal_asym_Constrain.png",istate));

	return;

} // end of main func.

