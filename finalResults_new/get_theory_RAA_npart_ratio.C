#include "SONGKYO.h"
#include "tdrstyle.C"
#include "TROOT.h"
#include "CMS_lumi_raaCent_rat.C"
#include "../cutsAndBin.h"
#include "../commonUtility.h"

void get_theory_RAA_npart_ratio(bool isArrow =true, int drawState=2)
{
  setTDRStyle();
  writeExtraText = false;       // if extra text
  int iPeriod = 100; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nState = 3; // Y(1S), Y(2S), and Y(3S)
  double xmax = 420.0;
  double xmin_int = 0.6;
  double xmax_int = 1.8;
  double boxw = 6.5; // for syst. box (vs cent)
  double boxw_int = 0.09;


  ///////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn[nState];
	TGraphErrors* gRAA[nState]; // vs centrality
  TGraphAsymmErrors* gRAA_sys[nState];
  for (int is=0; is<nState; is++){
  	fIn[is] = new TFile(Form("Ups_%d_RAA.root",is+1),"READ");
    gRAA[is]=(TGraphErrors*)fIn[is]->Get("gRAA_cent");
    gRAA_sys[is]= new TGraphAsymmErrors();
  }
  //// read input file : syst.
  TFile* fInSys_Hi[nState];
  TFile* fInSys_Lo[nState];
  TH1D* hSys_Hi[nState];
  TH1D* hSys_Lo[nState];
  int npoint_Hi[nState];
  int npoint_Lo[nState];
  for (int is=0; is<nState; is++){
    fInSys_Hi[is] = new TFile(Form("../Systematic_new/postCWR_mergedSys_constrain_ups%ds_asymHi.root",is+1),"READ");
    fInSys_Lo[is] = new TFile(Form("../Systematic_new/postCWR_mergedSys_constrain_ups%ds_asymLo.root",is+1),"READ");
    hSys_Hi[is]=(TH1D*)fInSys_Hi[is]->Get("hcentRAA_merged");
    hSys_Lo[is]=(TH1D*)fInSys_Lo[is]->Get("hcentRAA_merged");
    npoint_Hi[is] = hSys_Hi[is]->GetSize()-2;
    npoint_Lo[is] = hSys_Lo[is]->GetSize()-2;
    cout << "*** CENT Hi *** Y("<<is+1<<") : # of point = " << npoint_Hi[is] << endl;
    cout << "*** CENT Lo *** Y("<<is+1<<") : # of point = " << npoint_Lo[is] << endl;
  }   
  
  //// set bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys_Lo, relsys_Hi;
  double xshift = 0.2;
  //// --- vs centrality
  for (int is=0; is<nState; is++){
    cout << is+1 <<"th state***************" << endl;
    if (npoint_Lo[is] != gRAA[is]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
    for (int ipt=0; ipt< npoint_Lo[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys_Lo=0; relsys_Hi=0;
      gRAA[is]->GetPoint(ipt, pxtmp, pytmp);
      extmp=gRAA[is]->GetErrorX(ipt);
      eytmp=gRAA[is]->GetErrorY(ipt);
      relsys_Hi=hSys_Hi[is]->GetBinContent(npoint_Hi[is]-ipt);
      relsys_Lo=hSys_Lo[is]->GetBinContent(npoint_Lo[is]-ipt);
      cout << ipt <<"th bin RAA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      //cout << ipt <<"th bin rel. syst. = " << relsys << endl;
      cout << ipt <<"th bin Hi syst. = " << pytmp*relsys_Hi << endl; 
      cout << ipt <<"th bin Lo syst. = " << pytmp*relsys_Lo << endl; 
      //// 1) remove ex from gRAA
      gRAA[is]->SetPointError(ipt, 0, eytmp);
      //// 2) set ey for gRAA_sys 
      gRAA_sys[is]->SetPoint(ipt, pxtmp,pytmp);
      gRAA_sys[is]->SetPointError(ipt, boxw, boxw, pytmp*relsys_Lo,pytmp*relsys_Hi);
    }
  }

  TGraphErrors *gRAA_fR[nState];
  TGraphAsymmErrors *gRAA_sysfR[nState];
  gRAA_fR[0] = (TGraphErrors*) gRAA[0]->Clone("gRAA_fR1S");
  gRAA_fR[1] = (TGraphErrors*) gRAA[1]->Clone("gRAA_fR2S");
  gRAA_fR[2] = (TGraphErrors*) gRAA[2]->Clone("gRAA_fR3S");
  gRAA_sysfR[0] = (TGraphAsymmErrors*) gRAA_sys[0]->Clone("gRAA_fR1S");
  gRAA_sysfR[1] = (TGraphAsymmErrors*) gRAA_sys[1]->Clone("gRAA_fR2S");
  gRAA_sysfR[2] = (TGraphAsymmErrors*) gRAA_sys[2]->Clone("gRAA_fR3S");


  double lower68_2s = 0;
  double lower95_2s = 0;
  double upper68_2s = upper68_2s_c1;
  double upper95_2s = upper95_2s_c1;

  TBox *box68per_2s;
  TArrow *arr95per_2s;
  pxtmp=0; pytmp=0; extmp=0; eytmp=0; 
  //remove Y(3S) points
  if(drawState==2){
    gRAA[drawState-1]->GetPoint(8,pxtmp,pytmp);
    box68per_2s = new TBox(pxtmp-boxw,lower68_2s,pxtmp+boxw,upper68_2s);
    arr95per_2s = new TArrow(pxtmp,lower95_2s,pxtmp,upper95_2s,0.025,"<-|");
    gRAA[drawState-1]->SetPoint(8,-10,-10);
    gRAA[drawState-1]->SetPointError(8,0,0);
    gRAA_sys[drawState-1]->SetPoint(8,-10,-10);
    gRAA_sys[drawState-1]->SetPointError(8,0,0,0,0);
    box68per_2s->SetLineColor(kBlack);
    box68per_2s->SetFillColor(kBlack);
    box68per_2s->SetLineWidth(1);
    arr95per_2s->SetLineWidth(1);
    arr95per_2s->SetLineColor(kBlack);
    box68per_2s->SetFillStyle(3005);
  }
   

  //// graph style 
  for (int is=0; is<nState; is++){
    SetGraphStyle(gRAA[is], 4, 4);  // to keep it black/gray
    gRAA[is]->SetMarkerStyle(kFullCircle);
    SetGraphStyleSys(gRAA_sys[is], 4); 
	}
   
  double xlonger = 60;
  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.0551);
   
  //// axis et. al
  gRAA_sys[drawState-1]->GetXaxis()->SetTitle(" ");
  gRAA_sys[drawState-1]->GetXaxis()->CenterTitle();
  gRAA_sys[drawState-1]->GetYaxis()->SetTitle("R_{AA}");
  gRAA_sys[drawState-1]->GetYaxis()->CenterTitle();
  gRAA_sys[drawState-1]->GetXaxis()->SetLimits(0.,xmax);
  gRAA_sys[drawState-1]->SetMinimum(0.0);
  gRAA_sys[drawState-1]->SetMaximum(1.3);
  //// for cent
  gRAA_sys[drawState-1]->GetXaxis()->SetTitleSize(0);
  gRAA_sys[drawState-1]->GetYaxis()->SetTitleSize(0.076*1.0);
  gRAA_sys[drawState-1]->GetYaxis()->SetTitleOffset(0.8*1.0);
  gRAA_sys[drawState-1]->GetXaxis()->SetLabelSize(0);
  gRAA_sys[drawState-1]->GetYaxis()->SetLabelSize(0.07*1.0);
  
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600+xlonger,600);
  TPad* pad_diff = new TPad("pad_diff", "",0, 0.4, 1.0, 1.0); // vs centrality
  pad_diff->SetRightMargin(0.05);
  pad_diff->SetBottomMargin(0.024);
  pad_diff->SetTopMargin(0.080);

  TPad* pad_rat = new TPad("pad_rat","",0,0, 1.0,0.4);
  pad_rat->SetBottomMargin(0.3);
  pad_rat->SetTopMargin(0);
  pad_rat->SetRightMargin(0.05);
    
  //// --- 1st pad!!!   
  c1->cd();
  pad_diff->Draw(); 
  pad_diff->cd(); 

  gRAA_sys[drawState-1]->Draw("A5");
  
  dashedLine(0.,1.,xmax,1.,1,1);

  //// legend
  double data_leg_xpos1 = 0.76;
  double data_leg_xpos2 = 0.98;
  double data_leg_ypos1 = 0.251;
  double data_leg_ypos2 = 0.528;
  double legtextsize = 0.0571;


  if(drawState==1){
    data_leg_xpos1 = 0.20;
    data_leg_xpos2 = 0.42;
    data_leg_ypos1 = 0.03;
    data_leg_ypos2 = 0.265;
  }

  TLegend *leg= new TLegend(data_leg_xpos1,data_leg_ypos1,data_leg_xpos2,data_leg_ypos2);
  SetLegendStyle(leg);
  leg->SetTextSize(legtextsize);
  leg -> AddEntry(gRAA[drawState-1],"Data","lp");

  TArrow *arrLeg_2s = new TArrow(331.,0.358,331.,0.418,0.019,"<-|");
/*  if(drawState==1){
    TLegendEntry *ent2s=leg->AddEntry("","","");
    ent2s->SetLineColor(kGray+2);
    ent2s->SetFillColor(kGray+2);
    ent2s->SetFillStyle(3005);
    ent2s=leg->AddEntry("","","");

    arrLeg_2s->SetLineColor(kBlack);
    arrLeg_2s->SetLineWidth(2);
  }*/
  if(drawState==2){
    TLegendEntry *ent2s=leg->AddEntry("ent"," 68\% CL","f");
    ent2s->SetLineColor(kBlack);
    ent2s->SetLineWidth(2);
    ent2s->SetFillColor(kBlack);
    ent2s->SetFillStyle(3005);
    ent2s=leg->AddEntry("ent"," 95\% CL","");

    arrLeg_2s->SetLineColor(kBlack);
    arrLeg_2s->SetLineWidth(2);
  }
  leg->Draw("same");
  if(drawState==2) arrLeg_2s->Draw();

  TLatex* globtexT = new TLatex();
  globtexT->SetNDC();
  globtexT->SetTextAlign(12); //left-center
  globtexT->SetTextFont(42);
  globtexT->SetTextSize(0.0751);
  if(drawState==1) globtexT->DrawLatex(0.207, 0.24, Form("#varUpsilon(%dS)",drawState));
  if(drawState==2) globtexT->DrawLatex(0.772, 0.58, Form("#varUpsilon(%dS)",drawState));
  
  
  //// draw text
  double sz_init = 0.847; double sz_step = 0.0558;
  //globtex->DrawLatex(0.2435+0.1, sz_init, "p_{T} < 30 GeV");
  //globtex->DrawLatex(0.2435+0.1, sz_init-sz_step, "|y| < 2.4");
  globtex->DrawLatex(0.1135+0.1, sz_init+0.01, "p_{T} < 30 GeV");
  globtex->DrawLatex(0.1135+0.1, sz_init-sz_step+0.01, "|y| < 2.4");



  //Add Rapp
  TFile *frapp = new TFile("TheoryCurve/Rapp_RAA_5023.root","READ");
  
  TGraphErrors *gRAA_max[3]; 
  TGraphErrors *gRAA_min[3]; 
  TGraphErrors *gRAA_shade[3]; 
  
  for(int i=0;i<3;i++)
  {
    gRAA_max[i] = (TGraphErrors*) frapp-> Get(Form("RAA_rapp_nPart_%dS_max",i+1));
    gRAA_min[i] = (TGraphErrors*) frapp-> Get(Form("RAA_rapp_nPart_%dS_min",i+1));
    gRAA_shade[i] = (TGraphErrors*) frapp-> Get(Form("RAA_%ds_shade",i+1));
    gRAA_max[i] -> SetLineWidth(2.);
    gRAA_min[i] -> SetLineWidth(2.0);
    gRAA_shade[i] -> SetLineWidth(0);
  }
  gRAA_max[0]->SetLineColor(kGreen+3);
  gRAA_max[1]->SetLineColor(kGreen+3);
  gRAA_max[2]->SetLineColor(kGreen+3);
  
  gRAA_min[0]->SetLineColor(kGreen+3);
  gRAA_min[1]->SetLineColor(kGreen+3);
  gRAA_min[2]->SetLineColor(kGreen+3);
   
  /*gRAA_shade[0]->SetFillStyle(3004);
  gRAA_shade[1]->SetFillStyle(3004);
  gRAA_shade[2]->SetFillStyle(3004);
  */
  
  int ci = TColor::GetColor("#299617");
  gRAA_shade[0]->SetFillColor(ci);
  gRAA_shade[1]->SetFillColor(ci);
  gRAA_shade[2]->SetFillColor(ci);
  gRAA_shade[drawState-1]->Draw("f");
  
  gRAA_sys[drawState-1]->Draw("5");

  double line_y = 0.99;
  double line_y_diff = 0.07;
  double line_y_diff_in = 0.025;
  double line_x_start = 108;//97
  double line_x_end = line_x_start+30;//122
  double linematch_diff_X = 125.2;
  double linematch_diff_Y = 0.054;

  double legstrick_pos_x1 = 0.597;
  double legstrick_pos_x2 = 0.93;
  double legstrick_pos_y1 = 0.61;
  double legstrick_pos_y2 = 0.705;

  double state_diff_x = 108;
  if(drawState==2){
    line_x_start = line_x_start - state_diff_x;
    line_x_end = line_x_end - state_diff_x;
    legstrick_pos_x1 = legstrick_pos_x1 -0.203;
    legstrick_pos_x2 = legstrick_pos_x2 -0.203;
  }

  TLegend *leg_strick= new TLegend(legstrick_pos_x1,legstrick_pos_y1,legstrick_pos_x2,legstrick_pos_y2);
  SetLegendStyle(leg_strick);
  leg_strick->SetTextSize(0.050);
  leg_strick->AddEntry(gRAA_shade[drawState-1]," ","f");
  leg_strick->Draw("same");

  drawText2("Krouppa, Strickland",line_x_start+165.3, line_y+0.084-0.29,19);
  drawText2("Du, He, Rapp",line_x_start+164.91, line_y+0.084-0.18,19);
   
  drawText2("4#pi #eta/s=1", line_x_start+165.5, line_y-0.295, 19);
  drawText2("4#pi #eta/s=2", line_x_start+165.5, line_y-line_y_diff*1-0.295 - line_y_diff_in+0.00075, 19);
  drawText2("4#pi #eta/s=3", line_x_start+165.5, line_y-line_y_diff*2-0.295 - line_y_diff_in*2, 19);
  
  TLine* t1 = new TLine(line_x_start+5+linematch_diff_X,line_y-linematch_diff_Y-line_y_diff_in*8.8 ,line_x_end+linematch_diff_X,line_y-linematch_diff_Y-line_y_diff_in*8.8);
  TLine* t2 = new TLine(line_x_start+5+linematch_diff_X,line_y-linematch_diff_Y-line_y_diff_in*12.38,line_x_end+linematch_diff_X,line_y-linematch_diff_Y-line_y_diff_in*12.38);
  TLine* t3 = new TLine(line_x_start+5+linematch_diff_X,line_y-linematch_diff_Y-line_y_diff_in*16.23,line_x_end+linematch_diff_X,line_y-linematch_diff_Y-line_y_diff_in*16.23);
  
  t1->SetLineColor(kOrange+7);
  t2->SetLineColor(kOrange+7);
  t3->SetLineColor(kOrange+7);

  t1->SetLineStyle(3);
  //t1->SetLineStyle(3);
  t1->SetLineWidth(3);
  t1->Draw("same");
  
  t2->SetLineStyle(1);
  t2->SetLineWidth(3);
  t2->Draw("same");

  t3->SetLineStyle(8);
  t3->SetLineWidth(3);
  t3->Draw("same");
  
  //Add Strickland 
  TFile *fstrickland = new TFile("TheoryCurve/StrickLand_RAA_5023.root","READ");
  
  TGraphErrors *gRAA_1S_strickland[3]; 
  TGraphErrors *gRAA_2S_strickland[3]; 
  TGraphErrors *gRAA_3S_strickland[3]; 
  
  for(int i=0;i<3;i++)
  {
    gRAA_1S_strickland[i] = (TGraphErrors*) fstrickland-> Get(Form("RAA_strick_nPart_1S_%d",i));
    gRAA_2S_strickland[i] = (TGraphErrors*) fstrickland-> Get(Form("RAA_strick_nPart_2S_%d",i));
    gRAA_3S_strickland[i] = (TGraphErrors*) fstrickland-> Get(Form("RAA_strick_nPart_3S_%d",i));
    gRAA_1S_strickland[i] -> SetLineWidth(3.);
    gRAA_2S_strickland[i] -> SetLineWidth(3.0);
    gRAA_3S_strickland[i] -> SetLineWidth(3);
  }
  gRAA_1S_strickland[0]->SetLineColor(kOrange+7);
  gRAA_1S_strickland[1]->SetLineColor(kOrange+7);
  gRAA_1S_strickland[2]->SetLineColor(kOrange+7);
  gRAA_1S_strickland[0]->SetLineStyle(3);
  gRAA_1S_strickland[1]->SetLineStyle(1);
  gRAA_1S_strickland[2]->SetLineStyle(8);
  
  gRAA_2S_strickland[0]->SetLineColor(kOrange+7);
  gRAA_2S_strickland[1]->SetLineColor(kOrange+7);
  gRAA_2S_strickland[2]->SetLineColor(kOrange+7);
  gRAA_2S_strickland[0]->SetLineStyle(3);
  gRAA_2S_strickland[1]->SetLineStyle(1);
  gRAA_2S_strickland[2]->SetLineStyle(8);
  
  gRAA_3S_strickland[0]->SetLineColor(kGreen+2);
  gRAA_3S_strickland[1]->SetLineColor(kGreen+2);
  gRAA_3S_strickland[2]->SetLineColor(kGreen+2);
  gRAA_3S_strickland[0]->SetLineStyle(3);
  gRAA_3S_strickland[1]->SetLineStyle(1);
  gRAA_3S_strickland[2]->SetLineStyle(8);
 
 Int_t np = gRAA_1S_strickland[0]->GetN();
 TGraph *g_st_sh = new TGraph(2*np);
 double rp_x_st, rp_y_st_max, rp_y_st_min;
 for(int ipt=0; ipt<np; ipt++){
   rp_x_st=0; rp_y_st_max=0; rp_y_st_min=0;
   if(drawState==1) gRAA_1S_strickland[0]->GetPoint(ipt,rp_x_st,rp_y_st_max);
   else if(drawState==2) gRAA_2S_strickland[0]->GetPoint(ipt,rp_x_st,rp_y_st_max);
   g_st_sh->SetPoint(ipt, rp_x_st, rp_y_st_max);
   rp_x_st=0; rp_y_st_max=0; rp_y_st_min=0;
   if(drawState==1) gRAA_1S_strickland[2]->GetPoint(np-ipt-1,rp_x_st,rp_y_st_min);
   else if(drawState==2) gRAA_2S_strickland[2]->GetPoint(np-ipt-1,rp_x_st,rp_y_st_min);
   g_st_sh->SetPoint(np+ipt,rp_x_st,rp_y_st_min);
 }

  g_st_sh->SetFillStyle(3005);
  g_st_sh->SetFillColor(kOrange+7);
  g_st_sh->SetLineWidth(1);
  if(drawState==1) {
    gRAA_1S_strickland[2]->Draw("same"); 
    gRAA_1S_strickland[0]->Draw("same"); 
    gRAA_1S_strickland[1]->Draw("same"); 
  }
  if(drawState==2) {
    gRAA_2S_strickland[2]->Draw("same"); 
    gRAA_2S_strickland[0]->Draw("same"); 
    gRAA_2S_strickland[1]->Draw("same"); 
  }
  //g_st_sh->Draw("f");

  TLegend *leg_st= new TLegend(0.37, 0.57, 0.65, 0.725);
  SetLegendStyle(leg_st);
  leg_st->SetTextSize(0.050);
//  leg_st->AddEntry(g_st_sh," ","f");
  leg_st->Draw("same");


  gRAA_shade[drawState-1]->SetFillColor(ci);
  

  if(drawState==2){
    box68per_2s->Draw("lf same");
    arr95per_2s->Draw();
  }
  
  gRAA[drawState-1]->Draw("P same");

  //Global Unc.
  TFile* fppInt = new TFile("../fitResults/Constrain/fitresults_upsilon_fixParm1_seed2_DoubleCB_PP_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0.root");
  TH1D* hppInt = (TH1D*) fppInt -> Get("fitResults");
  double relsys_ppInt1S = hppInt->GetBinError(1)/hppInt->GetBinContent(1);
  double relsys_ppInt2S = hppInt->GetBinError(2)/hppInt->GetBinContent(2);

  TH1D* hSys_glb[nState];
  double sys_global_pp[nState];
  double sys_global_val;
  double accept_sys;
  TFile* f_acc[nState]; 
  TH1D* hSys_glb_acc[nState];
  for(int is=0; is<nState; is++){
    hSys_glb[is] = (TH1D*) fInSys_Lo[is]->Get("hintPP_merged");
    f_acc[is] = new TFile(Form("../acceptance_new/sys_acceptance_ups%dS_1804.root",is+1));
    hSys_glb_acc[is] = (TH1D*) f_acc[is]->Get("hcentSysPP");
    accept_sys = hSys_glb_acc[is]->GetBinContent(1);
    sys_global_pp[is] = TMath::Sqrt(hSys_glb[is]->GetBinContent(1)*hSys_glb[is]->GetBinContent(1)+accept_sys*accept_sys);
  } 
  
  sys_global_val = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp + nMB_unc*nMB_unc);
  double sys_global_y = sys_global_val; 
  double sys_global_x = 15;
  double sys_pp_1S = sys_global_pp[0];
  double sys_pp_2S = sys_global_pp[1];

  sys_pp_1S = TMath::Sqrt(sys_pp_1S*sys_pp_1S+relsys_ppInt1S*relsys_ppInt1S);
  sys_pp_2S = TMath::Sqrt(sys_pp_2S*sys_pp_2S+relsys_ppInt2S*relsys_ppInt2S);

  TBox *globalUncBox = new TBox(xmax-sys_global_x*2,1-sys_global_y,xmax-sys_global_x,1+sys_global_y);
  globalUncBox -> SetFillColorAlpha(kGray+2,0);
  globalUncBox -> SetLineColor(kBlack);
  globalUncBox -> SetLineWidth(2);
//  globalUncBox -> Draw("l same");

  sys_pp_1S = TMath::Sqrt(sys_pp_1S*sys_pp_1S+sys_global_y*sys_global_y);
  sys_pp_2S = TMath::Sqrt(sys_pp_2S*sys_pp_2S+sys_global_y*sys_global_y);


  if(drawState==1){
    TBox *ppRefUncBox1S = new TBox(xmax-sys_global_x,1-sys_pp_1S,xmax,1+sys_pp_1S);
    ppRefUncBox1S -> SetFillColorAlpha(kGray+1,0.6);
    ppRefUncBox1S -> SetLineColor(kBlack);
    ppRefUncBox1S -> SetLineWidth(1);
    ppRefUncBox1S -> Draw("l same");
  }
  else if(drawState==2){
    TBox *ppRefUncBox2S = new TBox(xmax-sys_global_x,1-sys_pp_2S,xmax,1+sys_pp_2S);
    ppRefUncBox2S -> SetFillColorAlpha(kGray+1,0.6);
    ppRefUncBox2S -> SetLineColor(kBlack);
    ppRefUncBox2S -> SetLineWidth(1);
    ppRefUncBox2S -> Draw("l same");
  }

  pad_diff->Update();

  //Rapp ratio
  TFile* frapp_histo = new TFile("TheoryCurve/Rapp_5023_histo.root","read");
  TGraphErrors* g_st_hist[nState];
  TGraphErrors* g_st_hist_max[nState];
  TGraphErrors* g_st_hist_min[nState];

	TGraphAsymmErrors* gRAA_rat[nState]; 
  TGraphAsymmErrors* gRAA_sys_rat[nState];
    
  for(int is=0; is<nState;is++){
    g_st_hist[is] = (TGraphErrors*) frapp_histo -> Get(Form("g_val_%ds",is+1));
    g_st_hist_max[is] = (TGraphErrors*) frapp_histo -> Get(Form("g_val_%ds_max",is+1));
    g_st_hist_min[is] = (TGraphErrors*) frapp_histo -> Get(Form("g_val_%ds_min",is+1));
    gRAA_rat[is] = new TGraphAsymmErrors();
    gRAA_sys_rat[is] = new TGraphAsymmErrors();
  }

  double xplot_shift = 3;
  double np_ratx, np_raty, np_raty_min, np_raty_max, nerr_x, nerr_y, nerr_sysx, nerr_sysy, nerr_sysxl, nerr_sysyl , nerr_sysxh, nerr_sysyh; 
  double p_ratx, p_raty, p_valx, p_valy, err_valx, err_valy, err_sys_valxl, err_sys_valyl,err_sys_valxh, err_sys_valyh;
  for(int is = 0; is<nState; is++){
    for(int ipt=0; ipt < g_st_hist[is]->GetN(); ipt++){
      np_ratx = 0; np_raty =0; np_raty_min = 0; np_raty_max = 0; nerr_x =0; nerr_y =0; nerr_sysx =0; nerr_sysy = 0;
      p_ratx = 0; p_raty = 0; p_valx =0; p_valy =0; err_valx =0; err_valy = 0; err_sys_valxl = 0; err_sys_valyl =0; err_sys_valxh = 0; err_sys_valyh =0;
      g_st_hist[is]->GetPoint(ipt, p_ratx, p_raty);
      gRAA_fR[is] -> GetPoint(ipt, p_valx, p_valy);
      err_valx = gRAA_fR[is]->GetErrorX(ipt);
      err_valy = gRAA_fR[is]->GetErrorY(ipt);
      err_sys_valxl = gRAA_sysfR[is]->GetErrorXlow(ipt);
      err_sys_valyl = gRAA_sysfR[is]->GetErrorYlow(ipt);
      err_sys_valxh = gRAA_sysfR[is]->GetErrorXhigh(ipt);
      err_sys_valyh = gRAA_sysfR[is]->GetErrorYhigh(ipt);
      np_ratx = p_valx;
      np_raty = p_valy/p_raty;
      nerr_x = err_valx;
      nerr_y = err_valy / p_valy;
      nerr_sysxl = err_sys_valxl;
      nerr_sysyl = err_sys_valyl / p_valy ;
      nerr_sysxh = err_sys_valxh;
      nerr_sysyh = err_sys_valyh / p_valy ;
      nerr_sysyl = TMath::Sqrt(nerr_sysyl*nerr_sysyl+nerr_y*nerr_y);
      nerr_sysyh = TMath::Sqrt(nerr_sysyh*nerr_sysyh+nerr_y*nerr_y);
     
      gRAA_sys_rat[is] -> SetPoint(ipt, np_ratx, 1);
      gRAA_sys_rat[is] -> SetPointError(ipt, nerr_sysxl, nerr_sysxh, nerr_sysyl, nerr_sysyh);

      p_ratx = 0; p_raty = 0; p_valx =0; p_valy =0; 
      g_st_hist_max[is]->GetPoint(ipt, p_ratx, p_raty);
      gRAA_fR[is] -> GetPoint(ipt, p_valx, p_valy);
      np_ratx = p_valx;
      np_raty_max = p_valy/p_raty;
      p_ratx = 0; p_raty = 0; p_valx =0; p_valy =0; 
      g_st_hist_min[is]->GetPoint(ipt, p_ratx, p_raty);
      gRAA_fR[is] -> GetPoint(ipt, p_valx, p_valy);
      np_ratx = p_valx;
      np_raty_min = p_valy/p_raty;
  
      if( np_raty_min > np_raty_max){
        gRAA_rat[is] -> SetPoint(ipt, np_ratx+xplot_shift, np_raty);
        gRAA_rat[is] -> SetPointError(ipt, nerr_x, nerr_x, fabs(np_raty-np_raty_max), fabs(np_raty-np_raty_min));
      }
      else if(np_raty_min <= np_raty_max){
        gRAA_rat[is] -> SetPoint(ipt, np_ratx+xplot_shift, np_raty);
        gRAA_rat[is] -> SetPointError(ipt, nerr_x, nerr_x, fabs(np_raty-np_raty_min), fabs(np_raty-np_raty_max));
      }

    }
  }
  
  //Strickland ratio
  TFile* fstrick_histo = new TFile("TheoryCurve/Strickland_5023_histo.root","read");
  TGraphErrors* g_st_histst[nState];
  TGraphErrors* g_st_histst_max[nState];
  TGraphErrors* g_st_histst_min[nState];

  TGraphAsymmErrors* gRAA_ratst[nState];
    
  for(int is=0; is<nState;is++){
    g_st_histst[is] = (TGraphErrors*) fstrick_histo -> Get(Form("g_val_%ds_2",is+1));
    g_st_histst_min[is] = (TGraphErrors*) fstrick_histo -> Get(Form("g_val_%ds_1",is+1));
    g_st_histst_max[is] = (TGraphErrors*) fstrick_histo -> Get(Form("g_val_%ds_3",is+1));
    gRAA_ratst[is] = new TGraphAsymmErrors();
  }

  for(int is = 0; is<nState; is++){
    for(int ipt=0; ipt < g_st_histst[is]->GetN(); ipt++){
      np_ratx = 0; np_raty =0; np_raty_min = 0; np_raty_max = 0; nerr_x =0; nerr_y =0; nerr_sysx =0; nerr_sysy = 0;
      p_ratx = 0; p_raty = 0; p_valx =0; p_valy =0; err_valx =0; err_valy = 0; err_sys_valxl = 0; err_sys_valyl =0; err_sys_valxh = 0; err_sys_valyh =0;
      g_st_histst[is]->GetPoint(ipt, p_ratx, p_raty);
      gRAA_fR[is] -> GetPoint(ipt, p_valx, p_valy);
      err_valx = gRAA_fR[is]->GetErrorX(ipt);
      err_valy = gRAA_fR[is]->GetErrorY(ipt);
      np_ratx = p_valx;
      np_raty = p_valy/p_raty;
      nerr_x = err_valx;
      nerr_y = err_valy / p_valy;

      p_ratx = 0; p_raty = 0; p_valx =0; p_valy =0; 
      g_st_histst_max[is]->GetPoint(ipt, p_ratx, p_raty);
      gRAA_fR[is] -> GetPoint(ipt, p_valx, p_valy);
      np_ratx = p_valx;
      np_raty_max = p_valy/p_raty;

      p_ratx = 0; p_raty = 0; p_valx =0; p_valy =0; 
      g_st_histst_min[is]->GetPoint(ipt, p_ratx, p_raty);
      gRAA_fR[is] -> GetPoint(ipt, p_valx, p_valy);
      np_ratx = p_valx;
      np_raty_min = p_valy/p_raty;
  
      if(np_raty_min > np_raty_max){
        gRAA_ratst[is] -> SetPoint(ipt, np_ratx - xplot_shift, np_raty);
        gRAA_ratst[is] -> SetPointError(ipt, nerr_x, nerr_x, fabs(np_raty-np_raty_max), fabs(np_raty-np_raty_min));
      }
      else if(np_raty_min <= np_raty_max){
        gRAA_ratst[is] -> SetPoint(ipt, np_ratx - xplot_shift, np_raty);
        gRAA_ratst[is] -> SetPointError(ipt, nerr_x, nerr_x, fabs(np_raty-np_raty_min), fabs(np_raty-np_raty_max));
      }
    }
  }

  g_st_histst[drawState-1]->SetLineColor(kGreen+3);
  g_st_histst[drawState-1]->SetLineWidth(2);
//  g_st_histst[drawState-1]->Draw("same");  
  
  pad_diff->Update();
  CMS_lumi_raaCent_rat( pad_diff, iPeriod, iPos );
  pad_diff->Update();
  
  c1->Update();
  c1->cd();
  pad_rat->Draw();
  pad_rat->cd();
  for (int is=0; is<nState; is++){
    SetGraphStyleSys(gRAA_sys_rat[is], is); 
    gRAA_rat[is]-> SetMarkerColor(ci);
    gRAA_rat[is]-> SetMarkerStyle(kFullCircle);
    gRAA_rat[is]-> SetMarkerSize(1);
    gRAA_rat[is]-> SetLineColor(kGreen+2);
    gRAA_rat[is]-> SetLineWidth(1);
    gRAA_ratst[is]-> SetMarkerColor(kOrange+7);
    gRAA_ratst[is]-> SetLineColor(kOrange+7);
    gRAA_ratst[is]-> SetMarkerStyle(kFullCircle);
    gRAA_ratst[is]-> SetMarkerSize(1);
    gRAA_ratst[is]-> SetLineWidth(1);
	}

  double gl_unc_rt =0;
  if(drawState==1) gl_unc_rt = sys_pp_1S;
  else if(drawState==2) gl_unc_rt = sys_pp_2S;
  TBox *globalUncBox_rt = new TBox(xmax-sys_global_x+3, 1-gl_unc_rt, xmax, 1+gl_unc_rt);
  globalUncBox_rt -> SetFillColorAlpha(kGray+1,0.6);
  globalUncBox_rt -> SetLineColor(kBlack);
  globalUncBox_rt -> SetLineWidth(1);

  gRAA_sys_rat[drawState-1]->SetLineColor(kBlack); 
  gRAA_sys_rat[drawState-1]->SetLineWidth(2); 
  gRAA_sys_rat[drawState-1]->SetFillColorAlpha(kBlack,0.9); 
  gRAA_sys_rat[drawState-1]->SetFillStyle(3013); 
  gRAA_sys_rat[drawState-1]->Draw("A5");
  globalUncBox_rt -> Draw("l same");
//  gRAA_sys_rat[drawState-1]->Draw("5");
  gRAA_rat[drawState-1]->Draw("P same");
  gRAA_ratst[drawState-1]->Draw("P same");

  //  gRAA_rat[2]->Draw("P");
  gRAA_sys_rat[drawState-1]->GetXaxis()->SetTitle("<N_{part}>");
  gRAA_sys_rat[drawState-1]->GetXaxis()->CenterTitle();
  gRAA_sys_rat[drawState-1]->GetYaxis()->SetTitle("Data / Theory");
  gRAA_sys_rat[drawState-1]->GetYaxis()->CenterTitle();
  gRAA_sys_rat[drawState-1]->GetXaxis()->SetLimits(0.,xmax);
  gRAA_sys_rat[drawState-1]->SetMinimum(0.19);
  if(drawState==1) gRAA_sys_rat[drawState-1]->SetMaximum(2.199);
  else if(drawState==2) gRAA_sys_rat[drawState-1]->SetMaximum(2.999);
  gRAA_sys_rat[drawState-1]->GetXaxis()->SetTitleSize(0.12*1.0);
  gRAA_sys_rat[drawState-1]->GetYaxis()->SetTitleSize(0.1*1.0);
  gRAA_sys_rat[drawState-1]->GetYaxis()->SetTitleOffset(.59);
  gRAA_sys_rat[drawState-1]->GetXaxis()->SetTitleOffset(.99);
  gRAA_sys_rat[drawState-1]->GetXaxis()->SetLabelSize(0.1*1.0);
  gRAA_sys_rat[drawState-1]->GetYaxis()->SetLabelSize(0.1007*1.0);
  gRAA_sys_rat[drawState-1]->GetYaxis()->SetNdivisions(205);
  gRAA_sys_rat[drawState-1]->GetYaxis()->SetTickSize(0.02);
  gRAA_sys_rat[drawState-1]->GetXaxis()->SetTickSize(0.06);

  TLine* t_1 = new TLine(0,1,xmax,1);
  t_1->SetLineWidth(1);
  t_1->SetLineStyle(2);
  t_1->SetLineColor(1);
  t_1->Draw("same");

  if(drawState==1){
  TLegend *leg_model= new TLegend(0.52, 0.706, 0.75, 0.946);
  SetLegendStyle(leg_model);
  leg_model->SetTextSize(0.070);
  leg_model->AddEntry(gRAA_ratst[drawState-1],"Krouppa, Strickland","pe");
  leg_model->AddEntry(gRAA_rat[drawState-1],"Du, He, Rapp","pe");
  leg_model->Draw("same");
  }
  else if(drawState==2){
  TLegend *leg_model= new TLegend(0.52, 0.746, 0.75, 0.967);
  SetLegendStyle(leg_model);
  leg_model->SetTextSize(0.070);
  leg_model->AddEntry(gRAA_ratst[drawState-1],"Krouppa, Strickland","pe");
  leg_model->AddEntry(gRAA_rat[drawState-1],"Du, Rapp","pe");
  leg_model->Draw("same");
  }

  if(drawState==1){
    TLegend *leg_dataunc= new TLegend(0.18, 0.834, 0.5, 0.934);
    SetLegendStyle(leg_dataunc);
    leg_dataunc->SetTextSize(0.0705);
    leg_dataunc->AddEntry(gRAA_sys_rat[drawState-1],"Data total uncorrel. unc.","f");
    leg_dataunc->Draw("same");
  }
  else if(drawState==2){
    TLegend *leg_dataunc= new TLegend(0.18, 0.854, 0.5, 0.954);
    SetLegendStyle(leg_dataunc);
    leg_dataunc->SetTextSize(0.0705);
    leg_dataunc->AddEntry(gRAA_sys_rat[drawState-1],"Data total uncorrel. unc.","f");
    leg_dataunc->Draw("same");
  }
/*
  if(drawState==1){
    TBox *glt = new TBox(16,1.51,39.5,1.7);
    glt -> SetFillColorAlpha(kGray+3,0.4);
    glt -> SetLineColor(kBlack);
    glt -> SetLineWidth(1);
    glt -> Draw("l same");
    drawText2("Data correl. unc.",45.6,1.53,16);
  }

  if(drawState==2){
    TBox *glt = new TBox(16,2.06,39.5,2.32);
    glt -> SetFillColorAlpha(kGray+3,0.4);
    glt -> SetLineColor(kBlack);
    glt -> SetLineWidth(1);
    glt -> Draw("l same");
    drawText2("Data correl. unc.",45.35,2.088,16);
  }
*/
  pad_rat->Update();
	c1->Update();
//	c1->Modified();
  c1->SaveAs(Form("plots/theory_raa_npart_%dS_ratio.png",(int)drawState));
  c1->SaveAs(Form("plots/theory_raa_npart_%dS_ratio.pdf",(int)drawState));
  c1->SaveAs(Form("plots/theory_raa_npart_%dS_ratio.C",(int)drawState));

  /*
	///////////////////////////////////////////////////////////////////
	//// save as a root file
	TFile *outFile = new TFile("RAA_vs_cent.root", "RECREATE");
	outFile->cd();
	for (int is=0; is<nState; is++){
		gRAA_sys[is]->Write();	
		gRAA[is]->Write();	
	}
	outFile->Close();
*/	
	return;

} // end of main func.

